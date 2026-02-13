#!/usr/bin/env python3
"""
Comparison of Spatial-only vs Levenshtein distance clustering on pairwise databases.
Analyzes clustering ratios and statistics for specific sample pairs.
"""

import sys
import os
import time
import psutil
import numpy as np
import matplotlib.pyplot as plt


def levenshtein_distance(seq1, seq2):
    """Calculate normalized Levenshtein distance."""
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


def apply_sequence_reclustering(labels, variants, max_threshold=0.2):
    """Re-cluster INS variants using Levenshtein sequence similarity."""
    if labels is None:
        return labels
    
    labels = np.array(labels)
    new_labels = np.full_like(labels, -1)
    next_cluster_id = 0
    
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
                    
                    if levenshtein_distance(seq_i, seq_j) <= max_threshold:
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


def get_clustering_statistics(labels):
    """Calculate clustering statistics."""
    if labels is None:
        return {'total': 0, 'clustered': 0, 'unclustered': 0, 'ratio': 0.0, 'n_clusters': 0}
    
    labels = np.array(labels)
    total = len(labels)
    unclustered = np.sum(labels == -1)
    clustered = total - unclustered
    ratio = clustered / total if total > 0 else 0.0
    
    unique = set(labels.tolist())
    n_clusters = len(unique) - (1 if -1 in unique else 0)
    
    return {
        'total': total,
        'clustered': clustered,
        'unclustered': unclustered,
        'ratio': ratio,
        'n_clusters': n_clusters
    }


def benchmark_algorithm(coordinates, variants, algorithm, use_levenshtein=False):
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
    
    if use_levenshtein and labels is not None:
        labels = apply_sequence_reclustering(labels, variants, max_threshold=0.2)
    
    elapsed_time = time.time() - start_time
    mem_after = process.memory_info().rss
    memory_mb = max(0.001, (mem_after - mem_before) / 1024 / 1024)
    
    clustering_stats = get_clustering_statistics(labels)
    
    return {
        'time': elapsed_time,
        'memory': memory_mb,
        'labels': labels,
        **clustering_stats
    }


def get_db_label(db_file):
    """Extract sample pair name from database filename."""
    basename = os.path.basename(db_file).replace('.db', '')
    return basename


def run_comparison_benchmark(db_files, algorithms):
    """Benchmark algorithms on specific databases (Spatial vs Levenshtein only)."""
    results = {}
    for algo in algorithms:
        for mode in ['SPATIAL', 'LEVENSHTEIN']:
            key = f"{algo}_{mode}"
            results[key] = {
                'db_names': [], 
                'times': [], 
                'memory': [], 
                'total_variants': [],
                'clustered': [],
                'unclustered': [],
                'cluster_ratio': [],
                'n_clusters': [],
                'sv_stats': []
            }
    
    for db_file in db_files:
        db_label = get_db_label(db_file)
        print(f"\n{'='*80}")
        print(f"Processing: {db_label}")
        print(f"{'='*80}")
        
        coordinates, variants = load_database(db_file, chromosome='1')
        
        if coordinates is None or len(coordinates) == 0:
            print("  Skipped (no variants)")
            continue
        
        n_variants = len(variants)
        sv_stats = get_sv_statistics(variants)
        
        print(f"Total variants: {n_variants}")
        print(f"SV breakdown: {sv_stats}")
        
        for algo in algorithms:
            print(f"\n{algo}:")
            
            # Spatial only (no sequence)
            print(f"  Running spatial-only clustering...", end=' ', flush=True)
            result = benchmark_algorithm(coordinates, variants, algo, use_levenshtein=False)
            key = f"{algo}_SPATIAL"
            results[key]['db_names'].append(db_label)
            results[key]['times'].append(result['time'])
            results[key]['memory'].append(result['memory'])
            results[key]['total_variants'].append(n_variants)
            results[key]['clustered'].append(result['clustered'])
            results[key]['unclustered'].append(result['unclustered'])
            results[key]['cluster_ratio'].append(result['ratio'])
            results[key]['n_clusters'].append(result['n_clusters'])
            results[key]['sv_stats'].append(sv_stats)
            print(f"✓ (Clustered: {result['clustered']}/{result['total']}, "
                  f"Ratio: {result['ratio']:.2%}, Clusters: {result['n_clusters']})")
            
            # With Levenshtein
            print(f"  Running with Levenshtein...", end=' ', flush=True)
            result = benchmark_algorithm(coordinates, variants, algo, use_levenshtein=True)
            key = f"{algo}_LEVENSHTEIN"
            results[key]['db_names'].append(db_label)
            results[key]['times'].append(result['time'])
            results[key]['memory'].append(result['memory'])
            results[key]['total_variants'].append(n_variants)
            results[key]['clustered'].append(result['clustered'])
            results[key]['unclustered'].append(result['unclustered'])
            results[key]['cluster_ratio'].append(result['ratio'])
            results[key]['n_clusters'].append(result['n_clusters'])
            results[key]['sv_stats'].append(sv_stats)
            print(f"✓ (Clustered: {result['clustered']}/{result['total']}, "
                  f"Ratio: {result['ratio']:.2%}, Clusters: {result['n_clusters']})")
    
    return results


def plot_comparison(results, algorithms):
    """Generate comparison plots for Spatial vs Levenshtein."""
    colors = {
        'SPATIAL': '#E74C3C',
        'LEVENSHTEIN': '#2ECC71'
    }
    
    markers = {
        'SPATIAL': 'o',
        'LEVENSHTEIN': '^'
    }
    
    linestyles = {
        'SPATIAL': '-',
        'LEVENSHTEIN': ':'
    }
    
    output_files = []
    
    # Clustering Ratio Comparison
    for algo in algorithms:
        fig, ax = plt.subplots(figsize=(14, 8))
        
        for mode in ['SPATIAL', 'LEVENSHTEIN']:
            key = f"{algo}_{mode}"
            if key not in results or not results[key]['db_names']:
                continue
            
            x_pos = np.arange(len(results[key]['db_names']))
            y_data = np.array(results[key]['cluster_ratio']) * 100  # Convert to percentage
            
            label = {'SPATIAL': 'Spatial only', 'LEVENSHTEIN': '+ Levenshtein'}[mode]
            
            ax.plot(x_pos, y_data,
                   color=colors[mode],
                   marker=markers[mode],
                   markersize=10,
                   linewidth=2.5,
                   linestyle=linestyles[mode],
                   label=label,
                   alpha=0.9)
        
        ax.set_xlabel('Sample Pair', fontsize=14, fontweight='bold')
        ax.set_ylabel('Clustering Ratio (%)', fontsize=14, fontweight='bold')
        ax.set_title(f'{algo} - Clustering Ratio Comparison', fontsize=16, fontweight='bold')
        
        # Set x-tick labels
        db_names = results[f"{algo}_SPATIAL"]['db_names']
        ax.set_xticks(range(len(db_names)))
        ax.set_xticklabels(db_names, rotation=45, ha='right')
        
        ax.legend(fontsize=12, frameon=True, loc='best')
        ax.grid(True, alpha=0.3)
        
        plt.tight_layout()
        filename = f"clustering_ratio_{algo.lower()}.png"
        plt.savefig(filename, dpi=300, bbox_inches='tight')
        plt.close()
        output_files.append(filename)
        print(f"  ✓ {filename}")
        
        # Time Comparison
        fig, ax = plt.subplots(figsize=(14, 8))
        
        for mode in ['SPATIAL', 'LEVENSHTEIN']:
            key = f"{algo}_{mode}"
            if key not in results or not results[key]['db_names']:
                continue
            
            x_pos = np.arange(len(results[key]['db_names']))
            y_data = np.array(results[key]['times'])
            
            label = {'SPATIAL': 'Spatial only', 'LEVENSHTEIN': '+ Levenshtein'}[mode]
            
            ax.plot(x_pos, y_data,
                   color=colors[mode],
                   marker=markers[mode],
                   markersize=10,
                   linewidth=2.5,
                   linestyle=linestyles[mode],
                   label=label,
                   alpha=0.9)
        
        ax.set_xlabel('Sample Pair', fontsize=14, fontweight='bold')
        ax.set_ylabel('Time (seconds)', fontsize=14, fontweight='bold')
        ax.set_title(f'{algo} - Processing Time Comparison', fontsize=16, fontweight='bold')
        ax.set_yscale('log')
        
        ax.set_xticks(range(len(db_names)))
        ax.set_xticklabels(db_names, rotation=45, ha='right')
        
        ax.legend(fontsize=12, frameon=True, loc='best')
        ax.grid(True, alpha=0.3)
        
        plt.tight_layout()
        filename = f"processing_time_{algo.lower()}.png"
        plt.savefig(filename, dpi=300, bbox_inches='tight')
        plt.close()
        output_files.append(filename)
        print(f"  ✓ {filename}")
    
    return output_files


def print_summary_table(results, algorithms):
    """Print detailed summary table of results."""

    print(f"{'SPATIAL vs LEVENSHTEIN CLUSTERING COMPARISON':^120}")
    
    for algo in algorithms:
        print(f"\n{algo} RESULTS:")

        print(f"{'Database':<25} {'Mode':<15} {'Total SVs':<12} {'Clustered':<12} "
              f"{'Unclustered':<12} {'Ratio':<10} {'Clusters':<10} {'Time (s)':<12}")

        
        db_names = results[f"{algo}_SPATIAL"]['db_names']
        for i, db_name in enumerate(db_names):
            for mode in ['SPATIAL', 'LEVENSHTEIN']:
                key = f"{algo}_{mode}"
                mode_label = {'SPATIAL': 'Spatial only', 'LEVENSHTEIN': '+ Levenshtein'}[mode]
                
                print(f"{db_name if mode == 'SPATIAL' else '':<25} "
                      f"{mode_label:<15} "
                      f"{results[key]['total_variants'][i]:<12} "
                      f"{results[key]['clustered'][i]:<12} "
                      f"{results[key]['unclustered'][i]:<12} "
                      f"{results[key]['cluster_ratio'][i]:.2%}{'':<5} "
                      f"{results[key]['n_clusters'][i]:<10} "
                      f"{results[key]['times'][i]:<12.4f}")
            print()
        
        # Print SV type breakdown
        print(f"\nSV TYPE BREAKDOWN:")
        print("-" * 120)
        for i, db_name in enumerate(db_names):
            key = f"{algo}_SPATIAL"
            sv_stats = results[key]['sv_stats'][i]
            breakdown_str = ', '.join([f"{k}: {v}" for k, v in sorted(sv_stats.items())])
            print(f"{db_name:<25} {breakdown_str}")
        print()


def main():
    print(f"\n{'='*120}")
    print(f"{'SPATIAL vs LEVENSHTEIN CLUSTERING COMPARISON':^120}")
    print(f"{'='*120}\n")
    

    base_dir = "/proj/sens2023005/nobackup/wharf/dimam/dimam-sens2023005/Degree_project/SVDB/Platinum_experiments"
    
    db_files = [
        os.path.join(base_dir, "NA12877_NA12878.db"),
        os.path.join(base_dir, "NA12877_NA12879.db"),
        os.path.join(base_dir, "NA12879_NA12881.db")
    ]
    
    # Verify files exist
    print("Checking for database files:")
    missing_files = []
    for db_file in db_files:
        if os.path.exists(db_file):
            print(f"  ✓ {os.path.basename(db_file)}")
        else:
            print(f"  ✗ {os.path.basename(db_file)} - NOT FOUND")
            missing_files.append(db_file)
    
    if missing_files:
        print(f"\nError: {len(missing_files)} database file(s) not found!")
        sys.exit(1)
    
    algorithms = ['DBSCAN', 'OPTICS', 'INTERVAL_TREE']
    
    print(f"\nAlgorithms to test: {', '.join(algorithms)}")
    print("Modes: Spatial only, + Levenshtein")
    

    print("Running benchmarks...")

    
    results = run_comparison_benchmark(db_files, algorithms)

    print("Generating plots...")
 
    
    plot_files = plot_comparison(results, algorithms)
    
    print_summary_table(results, algorithms)
    
    print(f"{'✓ COMPARISON COMPLETE':^120}")
    print(f"\nGenerated {len(plot_files)} plots:")
    for f in plot_files:
        print(f"  • {f}")
    print()


if __name__ == "__main__":
    main()