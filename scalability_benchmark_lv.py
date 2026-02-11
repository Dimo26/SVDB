#!/usr/bin/env python3

import sys
import os
import time
import glob
import psutil
import numpy as np
import matplotlib.pyplot as plt


def hamming_distance(seq1, seq2):
    if not seq1 or not seq2:
        return 1.0
    seq1, seq2 = str(seq1).upper(), str(seq2).upper()
    max_len = max(len(seq1), len(seq2))
    min_len = min(len(seq1), len(seq2))
    mismatches = sum(1 for i in range(min_len) if seq1[i] != seq2[i])
    return (mismatches + (max_len - min_len)) / max_len if max_len > 0 else 0.0


def levenshtein_distance(seq1, seq2):
    if not seq1 or not seq2:
        return 1.0
    
    seq1 = str(seq1).upper()
    seq2 = str(seq2).upper()
    
    len1 = len(seq1)
    len2 = len(seq2)
    max_len = max(len1, len2)
    
    dp = [[0] * (len2 + 1) for _ in range(len1 + 1)]
    
    for i in range(len1 + 1):
        dp[i][0] = i
    for j in range(len2 + 1):
        dp[0][j] = j
    
    for i in range(1, len1 + 1):
        for j in range(1, len2 + 1):
            if seq1[i-1] == seq2[j-1]:
                dp[i][j] = dp[i-1][j-1]
            else:
                dp[i][j] = 1 + min(dp[i-1][j], dp[i][j-1], dp[i-1][j-1])
    
    return dp[len1][len2] / max_len if max_len > 0 else 0.0


def apply_sequence_reclustering(labels, variants, max_threshold=0.2, distance_func='hamming'):
    if labels is None:
        return labels
    
    labels = np.array(labels)
    new_labels = np.full_like(labels, -1)
    next_cluster_id = 0
    
    dist_func = levenshtein_distance if distance_func == 'levenshtein' else hamming_distance
    
    for spatial_label in sorted(set(labels.tolist())):
        if spatial_label == -1:
            continue
        
        indices = np.where(labels == spatial_label)[0]
        ins_with_seq, other_indices = [], []
        
        for idx in indices:
            var = variants[idx]
            if var['type'] == 'INS' and var['sequence']:
                ins_with_seq.append(idx)
            else:
                other_indices.append(idx)
        
        if other_indices:
            new_labels[other_indices] = next_cluster_id
            next_cluster_id += 1
        
        if ins_with_seq:
            assigned = set()
            for i in range(len(ins_with_seq)):
                if i in assigned:
                    continue
                
                group = [i]
                assigned.add(i)
                idx_i = ins_with_seq[i]
                seq_i = variants[idx_i]['sequence']
                
                for j in range(i + 1, len(ins_with_seq)):
                    if j in assigned:
                        continue
                    
                    idx_j = ins_with_seq[j]
                    seq_j = variants[idx_j]['sequence']
                    
                    if dist_func(seq_i, seq_j) <= max_threshold:
                        group.append(j)
                        assigned.add(j)
                
                for g in group:
                    new_labels[ins_with_seq[g]] = next_cluster_id
                next_cluster_id += 1
    
    new_labels[labels == -1] = -1 
    return new_labels


def load_database(db_file, chromosome='1'):
    """Load Chr 1 variants from database."""
    try:
        from svdb import database
        db = database.DB(db_file)
    except:
        sys.path.insert(0, os.getcwd())
        from svdb.database import DB
        db = DB(db_file)
    
    coordinates, variants = [], []
    
    query = f'SELECT * FROM SVDB WHERE chrA = "{chromosome}" AND chrB = "{chromosome}"'
    
    for row in db.query(query):
        var_type, chrA, chrB = row[0], row[1], row[2]
        posA, posB = int(row[3]), int(row[6])
        seq = row[11] if len(row) > 11 else ''
        
        coordinates.append([posA, posB])
        variants.append({
            'type': var_type,
            'sequence': seq,
            'posA': posA,
            'posB': posB
        })
    
    return np.array(coordinates), variants


def get_sv_statistics(variants):
    """Calculate SV type statistics."""
    sv_types = {}
    for var in variants:
        var_type = var['type']
        sv_types[var_type] = sv_types.get(var_type, 0) + 1
    return sv_types


def benchmark_algorithm(coordinates, variants, algorithm, use_sequence=False, distance_func='hamming'):
    """Run clustering algorithm and measure performance."""
    process = psutil.Process()
    mem_before = process.memory_info().rss
    
    start_time = time.time()
    
    try:
        if algorithm == 'DBSCAN':
            from svdb.export_module import DBSCAN
            labels = DBSCAN.cluster(coordinates, epsilon=500, m=2)
        elif algorithm == 'OPTICS':
            from svdb.optics_clustering import optics_cluster
            labels = optics_cluster(coordinates, min_samples=2, max_eps=500)
        elif algorithm == 'INTERVAL_TREE':
            from svdb.interval_tree_overlap import interval_tree_cluster
            labels = interval_tree_cluster(coordinates, max_distance=500)
        else:
            labels = None
    except Exception as e:
        print(f"Error in {algorithm}: {e}")
        labels = None
    
    if use_sequence and labels is not None:
        labels = apply_sequence_reclustering(labels, variants, max_threshold=0.2, distance_func=distance_func)
    
    elapsed_time = time.time() - start_time
    mem_after = process.memory_info().rss
    memory_mb = max(0.001, (mem_after - mem_before) / 1024 / 1024)
    
    n_clusters = 0
    if labels is not None:
        unique = set(labels.tolist())
        n_clusters = len(unique) - (1 if -1 in unique else 0)
    
    return {
        'time': elapsed_time,
        'memory': memory_mb,
        'clusters': n_clusters
    }


def get_sample_count(db_file):
    """Extract sample count from filename."""
    basename = os.path.basename(db_file)
    parts = basename.replace('.db', '').split('_')
    for part in parts:
        clean = part.replace('samples', '').replace('sample', '')
        if clean.isdigit():
            return int(clean)
    return 0


def run_scalability_benchmark(db_files, algorithms):
    """Benchmark all algorithms across database sizes."""
    db_files_sorted = sorted(db_files, key=get_sample_count)
    
    results = {}
    for algo in algorithms:
        for mode in ['NO_SEQ', 'HAMMING', 'LEVENSHTEIN']:
            key = f"{algo}_{mode}"
            results[key] = {'sizes': [], 'times': [], 'memory': [], 'variants': [], 'clusters': []}
    
    sv_stats_per_db = {}
    
    for db_file in db_files_sorted:
        sample_count = get_sample_count(db_file)
        print(f"\nProcessing: {os.path.basename(db_file)} ({sample_count} samples)")
        
        coordinates, variants = load_database(db_file, chromosome='1')
        
        if coordinates is None or len(coordinates) == 0:
            print("  Skipped (no variants)")
            continue
        
        n_variants = len(variants)
        sv_stats = get_sv_statistics(variants)
        sv_stats_per_db[sample_count] = {'total': n_variants, 'breakdown': sv_stats}
        
        print(f"  Variants: {n_variants}")
        print(f"  SV breakdown: {sv_stats}")
        
        for algo in algorithms:
            print(f"  {algo}...", end=' ', flush=True)
            
            result = benchmark_algorithm(coordinates, variants, algo, use_sequence=False)
            key = f"{algo}_NO_SEQ"
            results[key]['sizes'].append(sample_count)
            results[key]['times'].append(result['time'])
            results[key]['memory'].append(result['memory'])
            results[key]['variants'].append(n_variants)
            results[key]['clusters'].append(result['clusters'])
            
            result = benchmark_algorithm(coordinates, variants, algo, use_sequence=True, distance_func='hamming')
            key = f"{algo}_HAMMING"
            results[key]['sizes'].append(sample_count)
            results[key]['times'].append(result['time'])
            results[key]['memory'].append(result['memory'])
            results[key]['variants'].append(n_variants)
            results[key]['clusters'].append(result['clusters'])
            
            result = benchmark_algorithm(coordinates, variants, algo, use_sequence=True, distance_func='levenshtein')
            key = f"{algo}_LEVENSHTEIN"
            results[key]['sizes'].append(sample_count)
            results[key]['times'].append(result['time'])
            results[key]['memory'].append(result['memory'])
            results[key]['variants'].append(n_variants)
            results[key]['clusters'].append(result['clusters'])
            
            print("done")
    
    return results, sv_stats_per_db


def plot_scalability_comparison(results, algorithms):
    """Generate scalability comparison plots."""
    colors = {
        'NO_SEQ': '#E74C3C',
        'HAMMING': '#3498DB',
        'LEVENSHTEIN': '#2ECC71'
    }
    
    markers = {
        'NO_SEQ': 'o',
        'HAMMING': 's',
        'LEVENSHTEIN': '^'
    }
    
    linestyles = {
        'NO_SEQ': '-',
        'HAMMING': '--',
        'LEVENSHTEIN': ':'
    }
    
    output_files = []
    
    for algo in algorithms:
        fig, ax = plt.subplots(figsize=(14, 8))
        
        for mode in ['NO_SEQ', 'HAMMING', 'LEVENSHTEIN']:
            key = f"{algo}_{mode}"
            if key not in results or not results[key]['sizes']:
                continue
            
            x_data = np.array(results[key]['sizes'])
            y_data = np.array(results[key]['times'])
            
            label_suffix = {'NO_SEQ': 'Spatial only', 'HAMMING': '+ Hamming', 'LEVENSHTEIN': '+ Levenshtein'}[mode]
            label = label_suffix
            
            ax.plot(x_data, y_data,
                   color=colors[mode],
                   marker=markers[mode],
                   markersize=8,
                   linewidth=2.5,
                   linestyle=linestyles[mode],
                   label=label,
                   alpha=0.9)
        
        ax.set_xlabel('Database Size (Number of Samples)', fontsize=14, fontweight='bold')
        ax.set_ylabel('Time (seconds)', fontsize=14, fontweight='bold')
        ax.set_yscale('log')
        ax.legend(fontsize=12, frameon=True, loc='best')
        ax.grid(False)
        
        plt.tight_layout()
        filename = f"scalability_time_{algo.lower()}.png"
        plt.savefig(filename, dpi=300, bbox_inches='tight')
        plt.close()
        output_files.append(filename)
        print(f"  ✓ {filename}")
        
        fig, ax = plt.subplots(figsize=(14, 8))
        
        for mode in ['NO_SEQ', 'HAMMING', 'LEVENSHTEIN']:
            key = f"{algo}_{mode}"
            if key not in results or not results[key]['sizes']:
                continue
            
            x_data = np.array(results[key]['sizes'])
            y_data = np.array(results[key]['memory'])
            
            label_suffix = {'NO_SEQ': 'Spatial only', 'HAMMING': '+ Hamming', 'LEVENSHTEIN': '+ Levenshtein'}[mode]
            label = label_suffix
            
            ax.plot(x_data, y_data,
                   color=colors[mode],
                   marker=markers[mode],
                   markersize=8,
                   linewidth=2.5,
                   linestyle=linestyles[mode],
                   label=label,
                   alpha=0.9)
        
        ax.set_xlabel('Database Size (Number of Samples)', fontsize=14, fontweight='bold')
        ax.set_ylabel('Memory (MB)', fontsize=14, fontweight='bold')
        ax.set_yscale('log')
        ax.legend(fontsize=12, frameon=True, loc='best')
        ax.grid(False)
        
        plt.tight_layout()
        filename = f"scalability_memory_{algo.lower()}.png"
        plt.savefig(filename, dpi=300, bbox_inches='tight')
        plt.close()
        output_files.append(filename)
        print(f"  ✓ {filename}")
    
    return output_files


def print_summary_table(results, algorithms, sv_stats_per_db):
    """Print summary table of results."""

    print("SCALABILITY BENCHMARK SUMMARY")

    
    print("\nSV STATISTICS PER DATABASE:")
    print(f"{'Sample Count':<15} {'Total SVs':<15} {'SV Type Breakdown':<50}")
    print("-"*90)
    for sample_count in sorted(sv_stats_per_db.keys()):
        stats = sv_stats_per_db[sample_count]
        breakdown_str = ', '.join([f"{k}: {v}" for k, v in sorted(stats['breakdown'].items())])
        print(f"{sample_count:<15} {stats['total']:<15} {breakdown_str:<50}")
    
    for algo in algorithms:
        print(f"\n{algo}:")
        print(f"{'Size':<10} {'Mode':<15} {'Time (s)':<12} {'Memory (MB)':<15} {'Clusters':<12}")
        print("-"*70)
        
        sizes = results[f"{algo}_NO_SEQ"]['sizes']
        for i, size in enumerate(sizes):
            for mode in ['NO_SEQ', 'HAMMING', 'LEVENSHTEIN']:
                key = f"{algo}_{mode}"
                if i < len(results[key]['times']):
                    mode_label = {'NO_SEQ': 'Spatial', 'HAMMING': 'Hamming', 'LEVENSHTEIN': 'Levenshtein'}[mode]
                    print(f"{size if mode == 'NO_SEQ' else '':<10} {mode_label:<15} "
                          f"{results[key]['times'][i]:<12.4f} "
                          f"{results[key]['memory'][i]:<15.2f} "
                          f"{results[key]['clusters'][i]:<12}")
            print()


def main():

    print(f"{'SCALABILITY BENCHMARK: Hamming vs Levenshtein':^90}")

    
    db_files = glob.glob("*.db")
    sample_dbs = [f for f in db_files 
                  if 'sample' in f.lower() and 'all' not in f.lower()]
    
    if not sample_dbs:
        print("Error: No sample databases found (e.g., *_10samples.db, *_100samples.db)")
        sys.exit(1)
    
    print(f"\nFound {len(sample_dbs)} sample databases:")
    for db in sorted(sample_dbs, key=get_sample_count):
        print(f"  • {os.path.basename(db)} ({get_sample_count(db)} samples)")
    
    algorithms = ['DBSCAN', 'OPTICS', 'INTERVAL_TREE']
    
    print("\nRunning benchmarks...")
    
    results, sv_stats_per_db = run_scalability_benchmark(sample_dbs, algorithms)
    
    print("\nGenerating plots...")
    
    plot_files = plot_scalability_comparison(results, algorithms)
    
    print_summary_table(results, algorithms, sv_stats_per_db)
    

    print(f"{'✓ SCALABILITY BENCHMARK COMPLETE':^90}")

    print(f"\nGenerated {len(plot_files)} plots:")
    for f in plot_files:
        print(f"  • {f}")


if __name__ == "__main__":
    main()