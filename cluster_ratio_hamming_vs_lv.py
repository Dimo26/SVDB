#!/usr/bin/env python3


import sys
import os
import time
import psutil
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict


def hamming_distance(seq1, seq2):
    if not seq1 or not seq2:
        return 1.0
    
    seq1 = str(seq1).upper()
    seq2 = str(seq2).upper()
    min_len = min(len(seq1), len(seq2))
    max_len = max(len(seq1), len(seq2))
    

    mismatches = sum(1 for i in range(min_len) if seq1[i] != seq2[i])
    

    length_diff = max_len - min_len
    total_dist = mismatches + length_diff
    
    return total_dist / max_len if max_len > 0 else 0.0


def levenshtein_distance(seq1, seq2):
    """
    Calculate normalized Levenshtein (edit) distance.
    Accounts for insertions, deletions, and substitutions.
    """
    if not seq1 or not seq2:
        return 1.0
    
    seq1 = str(seq1).upper()
    seq2 = str(seq2).upper()
    
    len1 = len(seq1)
    len2 = len(seq2)
    max_len = max(len1, len2)
    
    # Dynamic programming table
    dp = [[0] * (len2 + 1) for _ in range(len1 + 1)]
    
    # Initialize base cases
    for i in range(len1 + 1):
        dp[i][0] = i
    for j in range(len2 + 1):
        dp[0][j] = j
    
    # Fill table
    for i in range(1, len1 + 1):
        for j in range(1, len2 + 1):
            if seq1[i-1] == seq2[j-1]:
                dp[i][j] = dp[i-1][j-1]  # No operation
            else:
                dp[i][j] = 1 + min(
                    dp[i-1][j],      # Deletion
                    dp[i][j-1],      # Insertion
                    dp[i-1][j-1]     # Substitution
                )
    
    return dp[len1][len2] / max_len if max_len > 0 else 0.0


def apply_sequence_reclustering(labels, variants, distance_func, max_threshold=0.2, method_name="Unknown"):

    if labels is None:
        return labels, {}
    
    labels = np.array(labels)
    new_labels = np.full_like(labels, -1)
    next_cluster_id = 0
    
    # Track how sequences are split
    split_info = {
        'splits_performed': 0,
        'insertions_processed': 0,
        'sequences_compared': 0,
        'method': method_name
    }
    
    for spatial_label in sorted(set(labels.tolist())):
        if spatial_label == -1:
            continue
        
        indices = np.where(labels == spatial_label)[0]
        ins_with_seq = []
        other_indices = []
        
        for idx in indices:
            var = variants[idx]
            if var['type'] == 'INS' and var['sequence']:
                ins_with_seq.append(idx)
            else:
                other_indices.append(idx)
        
        # Keep non-insertions together
        if other_indices:
            new_labels[other_indices] = next_cluster_id
            next_cluster_id += 1
        
        # Re-cluster insertions by sequence
        if len(ins_with_seq) > 0:
            split_info['insertions_processed'] += len(ins_with_seq)
            
            if len(ins_with_seq) == 1:
                # Single insertion, keep it
                new_labels[ins_with_seq[0]] = next_cluster_id
                next_cluster_id += 1
            else:
                # Multiple insertions - check sequences
                original_cluster_count = 1
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
                        
                        split_info['sequences_compared'] += 1
                        
                        # THIS IS THE KEY DIFFERENCE!
                        dist = distance_func(seq_i, seq_j)
                        
                        if dist <= max_threshold:
                            group.append(j)
                            assigned.add(j)
                    
                    # Assign this group to a cluster
                    for g in group:
                        new_labels[ins_with_seq[g]] = next_cluster_id
                    next_cluster_id += 1
                
                # Track if we split the spatial cluster
                final_cluster_count = len(set(new_labels[ins_with_seq]))
                if final_cluster_count > original_cluster_count:
                    split_info['splits_performed'] += 1
    
    # Preserve noise
    new_labels[labels == -1] = -1
    
    return new_labels, split_info


def load_database(db_file, chromosome='1'):
    """Load Chr 1 variants from database."""
    try:
        from svdb import database
        db = database.DB(db_file)
    except:
        sys.path.insert(0, os.getcwd())
        from svdb.database import DB
        db = DB(db_file)
    
    coordinates = []
    variants = []
    
    query = f'SELECT * FROM SVDB WHERE chrA = "{chromosome}" AND chrB = "{chromosome}"'
    
    for row in db.query(query):
        var_type = row[0]
        chrA, chrB = row[1], row[2]
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
        return {
            'total': 0,
            'clustered': 0,
            'unclustered': 0,
            'ratio': 0.0,
            'n_clusters': 0
        }
    
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


def benchmark_algorithm(coordinates, variants, algorithm, mode='spatial'):
    process = psutil.Process()
    mem_before = process.memory_info().rss
    start_time = time.time()
    
    split_info = {}
    
    try:
        # Step 1: Spatial clustering
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
        
        # Step 2: Apply sequence-based re-clustering if requested
        if mode == 'hamming' and labels is not None:
            labels, split_info = apply_sequence_reclustering(
                labels, variants, hamming_distance, 
                max_threshold=0.2, method_name="Hamming"
            )
        elif mode == 'levenshtein' and labels is not None:
            labels, split_info = apply_sequence_reclustering(
                labels, variants, levenshtein_distance,
                max_threshold=0.2, method_name="Levenshtein"
            )
        
    except Exception as e:
        print(f"Error in {algorithm} ({mode}): {e}")
        import traceback
        traceback.print_exc()
        labels = None
    
    elapsed_time = time.time() - start_time
    mem_after = process.memory_info().rss
    memory_mb = max(0.001, (mem_after - mem_before) / 1024 / 1024)
    
    clustering_stats = get_clustering_statistics(labels)
    
    return {
        'time': elapsed_time,
        'memory': memory_mb,
        'labels': labels,
        'split_info': split_info,
        **clustering_stats
    }


def get_db_label(db_file):
    """Extract sample pair name from database filename."""
    basename = os.path.basename(db_file).replace('.db', '')
    return basename


def run_comparison_benchmark(db_files, algorithms):
    """Benchmark algorithms with three modes: spatial, hamming, levenshtein."""
    
    modes = ['spatial', 'hamming', 'levenshtein']
    results = {}
    
    for algo in algorithms:
        for mode in modes:
            key = f"{algo}_{mode.upper()}"
            results[key] = {
                'db_names': [],
                'times': [],
                'memory': [],
                'total_variants': [],
                'clustered': [],
                'unclustered': [],
                'cluster_ratio': [],
                'n_clusters': [],
                'sv_stats': [],
                'split_info': []
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
        
        n_variants = len(coordinates)
        sv_stats = get_sv_statistics(variants)
        
        print(f"  Loaded {n_variants} variants from Chr 1")
        print(f"  SV types: {sv_stats}")
        
        for algo in algorithms:
            print(f"\n{algo}:")
            
            for mode in modes:
                mode_display = {
                    'spatial': 'Spatial only',
                    'hamming': 'Spatial + Hamming',
                    'levenshtein': 'Spatial + Levenshtein'
                }[mode]
                
                print(f"  Running {mode_display}...", end=' ', flush=True)
                
                result = benchmark_algorithm(coordinates, variants, algo, mode=mode)
                
                key = f"{algo}_{mode.upper()}"
                results[key]['db_names'].append(db_label)
                results[key]['times'].append(result['time'])
                results[key]['memory'].append(result['memory'])
                results[key]['total_variants'].append(n_variants)
                results[key]['clustered'].append(result['clustered'])
                results[key]['unclustered'].append(result['unclustered'])
                results[key]['cluster_ratio'].append(result['ratio'])
                results[key]['n_clusters'].append(result['n_clusters'])
                results[key]['sv_stats'].append(sv_stats)
                results[key]['split_info'].append(result.get('split_info', {}))
                
                print(f"✓ (Clustered: {result['clustered']}/{result['total']}, "
                      f"Ratio: {result['ratio']:.2%}, Clusters: {result['n_clusters']})")
                
                if result.get('split_info'):
                    info = result['split_info']
                    if info.get('insertions_processed', 0) > 0:
                        print(f"    └─ {info['method']}: {info['insertions_processed']} insertions, "
                              f"{info['sequences_compared']} comparisons, "
                              f"{info['splits_performed']} splits")
    
    return results


def plot_comparison(results, algorithms, output_dir='/mnt/user-data/outputs'):
    """Generate comprehensive comparison plots."""
    
    os.makedirs(output_dir, exist_ok=True)
    
    # Color scheme for three modes
    colors = {
        'SPATIAL': '#E74C3C',      # Red
        'HAMMING': '#F39C12',      # Orange
        'LEVENSHTEIN': '#2ECC71'   # Green
    }
    
    markers = {
        'SPATIAL': 'o',
        'HAMMING': 's',
        'LEVENSHTEIN': '^'
    }
    
    linestyles = {
        'SPATIAL': '-',
        'HAMMING': '--',
        'LEVENSHTEIN': ':'
    }
    
    output_files = []
    
    for algo in algorithms:
        # Plot 1: Clustering Ratio Comparison
        fig, ax = plt.subplots(figsize=(14, 8))
        
        for mode in ['SPATIAL', 'HAMMING', 'LEVENSHTEIN']:
            key = f"{algo}_{mode}"
            if key not in results or not results[key]['db_names']:
                continue
            
            x_pos = np.arange(len(results[key]['db_names']))
            y_data = np.array(results[key]['cluster_ratio']) * 100
            
            label = {
                'SPATIAL': 'Spatial only',
                'HAMMING': 'Spatial + Hamming',
                'LEVENSHTEIN': 'Spatial + Levenshtein'
            }[mode]
            
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
        ax.set_title(f'{algo} - Clustering Ratio: Spatial vs Hamming vs Levenshtein',
                    fontsize=16, fontweight='bold')
        
        db_names = results[f"{algo}_SPATIAL"]['db_names']
        ax.set_xticks(range(len(db_names)))
        ax.set_xticklabels(db_names, rotation=45, ha='right')
        
        ax.legend(fontsize=12, frameon=True, loc='best')

        
        plt.tight_layout()
        filename = os.path.join(output_dir, f"clustering_ratio_{algo.lower()}.png")
        plt.savefig(filename, dpi=300, bbox_inches='tight')
        plt.close()
        output_files.append(filename)
        print(f"  ✓ {os.path.basename(filename)}")
        
        # Plot 2: Number of Clusters Comparison
        fig, ax = plt.subplots(figsize=(14, 8))
        
        for mode in ['SPATIAL', 'HAMMING', 'LEVENSHTEIN']:
            key = f"{algo}_{mode}"
            if key not in results or not results[key]['db_names']:
                continue
            
            x_pos = np.arange(len(results[key]['db_names']))
            y_data = np.array(results[key]['n_clusters'])
            
            label = {
                'SPATIAL': 'Spatial only',
                'HAMMING': 'Spatial + Hamming',
                'LEVENSHTEIN': 'Spatial + Levenshtein'
            }[mode]
            
            ax.plot(x_pos, y_data,
                   color=colors[mode],
                   marker=markers[mode],
                   markersize=10,
                   linewidth=2.5,
                   linestyle=linestyles[mode],
                   label=label,
                   alpha=0.9)
        
        ax.set_xlabel('Sample Pair', fontsize=14, fontweight='bold')
        ax.set_ylabel('Number of Clusters', fontsize=14, fontweight='bold')
        ax.set_title(f'{algo} - Cluster Count: Impact of Sequence Distance Metric',
                    fontsize=16, fontweight='bold')
        
        ax.set_xticks(range(len(db_names)))
        ax.set_xticklabels(db_names, rotation=45, ha='right')
        
        ax.legend(fontsize=12, frameon=True, loc='best')

        plt.tight_layout()
        filename = os.path.join(output_dir, f"cluster_count_{algo.lower()}.png")
        plt.savefig(filename, dpi=300, bbox_inches='tight')
        plt.close()
        output_files.append(filename)
        print(f"  ✓ {os.path.basename(filename)}")
        
        # Plot 3: Time Comparison (INCLUDING HAMMING!)
        fig, ax = plt.subplots(figsize=(14, 8))
        
        for mode in ['SPATIAL', 'HAMMING', 'LEVENSHTEIN']:
            key = f"{algo}_{mode}"
            if key not in results or not results[key]['db_names']:
                continue
            
            x_pos = np.arange(len(results[key]['db_names']))
            y_data = np.array(results[key]['times'])
            
            label = {
                'SPATIAL': 'Spatial only',
                'HAMMING': 'Spatial + Hamming',
                'LEVENSHTEIN': 'Spatial + Levenshtein'
            }[mode]
            
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
        ax.set_title(f'{algo} - Processing Time: All Three Methods',
                    fontsize=16, fontweight='bold')
        ax.set_yscale('log')
        
        ax.set_xticks(range(len(db_names)))
        ax.set_xticklabels(db_names, rotation=45, ha='right')
        
        ax.legend(fontsize=12, frameon=True, loc='best')
        
        plt.tight_layout()
        filename = os.path.join(output_dir, f"processing_time_{algo.lower()}.png")
        plt.savefig(filename, dpi=300, bbox_inches='tight')
        plt.close()
        output_files.append(filename)
        print(f"  ✓ {os.path.basename(filename)}")
    
    return output_files


def print_summary_table(results, algorithms):
    """Print detailed summary table."""
    print(f"\n{'='*140}")
    print(f"{'THREE-WAY CLUSTERING COMPARISON: SPATIAL vs HAMMING vs LEVENSHTEIN':^140}")
    print(f"{'='*140}\n")
    
    for algo in algorithms:
        print(f"\n{algo} RESULTS:")
        print("-" * 140)
        
        print(f"{'Database':<20} {'Mode':<20} {'Total':<8} {'Clustered':<10} "
              f"{'Unclustered':<12} {'Ratio':<10} {'Clusters':<10} {'Time (s)':<10} {'Insertions':<12}")
        print("-" * 140)
        
        db_names = results[f"{algo}_SPATIAL"]['db_names']
        
        for i, db_name in enumerate(db_names):
            for mode in ['SPATIAL', 'HAMMING', 'LEVENSHTEIN']:
                key = f"{algo}_{mode}"
                mode_label = {
                    'SPATIAL': 'Spatial only',
                    'HAMMING': '+ Hamming',
                    'LEVENSHTEIN': '+ Levenshtein'
                }[mode]
                
                split_info = results[key]['split_info'][i] if results[key]['split_info'] else {}
                ins_info = f"{split_info.get('insertions_processed', 0)}" if split_info else "N/A"
                
                print(f"{db_name if mode == 'SPATIAL' else '':<20} "
                      f"{mode_label:<20} "
                      f"{results[key]['total_variants'][i]:<8} "
                      f"{results[key]['clustered'][i]:<10} "
                      f"{results[key]['unclustered'][i]:<12} "
                      f"{results[key]['cluster_ratio'][i]:.2%}{'':<5} "
                      f"{results[key]['n_clusters'][i]:<10} "
                      f"{results[key]['times'][i]:<10.4f} "
                      f"{ins_info:<12}")
            
            print()
        
        # Print difference analysis
        print(f"\nDIFFERENCE ANALYSIS:")

        
        for i, db_name in enumerate(db_names):
            spatial_key = f"{algo}_SPATIAL"
            hamming_key = f"{algo}_HAMMING"
            lev_key = f"{algo}_LEVENSHTEIN"
            
            spatial_clusters = results[spatial_key]['n_clusters'][i]
            hamming_clusters = results[hamming_key]['n_clusters'][i]
            lev_clusters = results[lev_key]['n_clusters'][i]
            
            hamming_diff = hamming_clusters - spatial_clusters
            lev_diff = lev_clusters - spatial_clusters
            hamming_vs_lev = hamming_clusters - lev_clusters
            
            print(f"{db_name:<20}")
            print(f"  Spatial → Hamming:     {hamming_diff:+4d} clusters ({hamming_clusters} total)")
            print(f"  Spatial → Levenshtein: {lev_diff:+4d} clusters ({lev_clusters} total)")
            print(f"  Hamming vs Levenshtein: {hamming_vs_lev:+4d} clusters")
            
            if hamming_clusters == lev_clusters:
                print(f"      This suggests insertions may not have sequences, or threshold is too loose.")
            elif abs(hamming_vs_lev) > 0:
                print(f"  ✓ Hamming and Levenshtein produce DIFFERENT results (as expected)")
            
            print()


def analyze_sequence_differences(results, algorithms):
    print(f"{'SEQUENCE DISTANCE ANALYSIS':^100}")

    
    print("Expected behavior:")
    print("  • Hamming: Counts mismatches + length differences")
    print("  • Levenshtein: Counts minimum edits (insertions/deletions/substitutions)")
    print("\nFor insertions with variable lengths or internal indels:")
    print("  • Levenshtein should be MORE ACCURATE")
    print("  • May produce DIFFERENT cluster assignments\n")


def main():
    print(f"\n{'='*100}")
    print(f"{'COMPREHENSIVE CLUSTERING COMPARISON':^100}")
    print(f"{'Spatial vs Hamming vs Levenshtein':^100}")
    print(f"{'='*100}\n")
    

    base_dir = "/proj/sens2023005/nobackup/wharf/dimam/dimam-sens2023005/Degree_project/SVDB/Platinum_experiments"
    
    db_files = [
        os.path.join(base_dir, "NA12877_NA12878.db"),
        os.path.join(base_dir, "NA12877_NA12879.db"),
        os.path.join(base_dir, "NA12879_NA12881.db")
    ]
    
    # Verify files
    print("Checking database files:")
    missing = []
    for db_file in db_files:
        if os.path.exists(db_file):
            print(f"  ✓ {os.path.basename(db_file)}")
        else:
            print(f"  ✗ {os.path.basename(db_file)} - NOT FOUND")
            missing.append(db_file)
    
    if missing:
        print(f"\nError: {len(missing)} file(s) not found!")
        sys.exit(1)
    
    algorithms = ['DBSCAN', 'OPTICS', 'INTERVAL_TREE']
    
    print(f"\nAlgorithms: {', '.join(algorithms)}")
    print("Modes: Spatial only, Hamming sequence, Levenshtein sequence")
    print("Threshold: 0.2 (20% difference allowed)\n")
    

    print("\nRunning benchmarks...")
    results = run_comparison_benchmark(db_files, algorithms)
    

    print("\nGenerating plots...")
    plot_files = plot_comparison(results, algorithms)
    

    print_summary_table(results, algorithms)
    
    analyze_sequence_differences(results, algorithms)
    

    print(f"{' COMPARISON COMPLETE':^100}")
    print(f"Generated {len(plot_files)} plots:")
    for f in plot_files:
        print(f"  • {os.path.basename(f)}")
    print()


if __name__ == "__main__":
    main()