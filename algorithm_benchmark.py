
#!/usr/bin/env python3

import sys
import os
import time
import glob
import psutil
import numpy as np
import matplotlib.pyplot as plt

# Add current directory to Python path
sys.path.insert(0, os.getcwd())

# Import SVDB modules
from svdb.database import DB
from svdb.export_module import DBSCAN
from svdb.optics_clustering import optics_cluster
from svdb.interval_tree_overlap import interval_tree_cluster


def load_database(db_file, chromosome='1'):
    """Load Chr 1 variants from database."""
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


def hamming_distance(seq1, seq2):
    """Normalized Hamming distance between sequences."""
    if not seq1 or not seq2:
        return 1.0
    seq1, seq2 = str(seq1).upper(), str(seq2).upper()
    max_len = max(len(seq1), len(seq2))
    min_len = min(len(seq1), len(seq2))
    mismatches = sum(1 for i in range(min_len) if seq1[i] != seq2[i])
    return (mismatches + (max_len - min_len)) / max_len if max_len > 0 else 0.0


def apply_hamming_reclustering(labels, variants, max_hamming=0.2):
    """Re-cluster insertions by sequence similarity using Hamming distance."""
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
        
        # Non-insertion variants keep cluster
        if other_indices:
            new_labels[other_indices] = next_cluster_id
            next_cluster_id += 1
        
        # Re-cluster insertions by sequence
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
                    
                    if hamming_distance(seq_i, seq_j) <= max_hamming:
                        group.append(j)
                        assigned.add(j)
                
                for g in group:
                    new_labels[ins_with_seq[g]] = next_cluster_id
                next_cluster_id += 1
    
    return new_labels


def benchmark_algorithm(coordinates, variants, algorithm, use_hamming=False):
    """Run clustering algorithm and measure time + memory."""
    process = psutil.Process()
    mem_before = process.memory_info().rss
    
    start_time = time.time()
    
    # Run clustering
    if algorithm == 'DBSCAN':
        labels = DBSCAN.cluster(coordinates, epsilon=500, m=2)
    elif algorithm == 'OPTICS':
        labels = optics_cluster(coordinates, min_samples=2, max_eps=2000)
    elif algorithm == 'INTERVAL_TREE':
        labels = interval_tree_cluster(coordinates, max_distance=1000)
    
    # Apply Hamming if requested
    if use_hamming and labels is not None:
        labels = apply_hamming_reclustering(labels, variants, max_hamming=0.2)
    
    elapsed_time = time.time() - start_time
    mem_after = process.memory_info().rss
    memory_mb = (mem_after - mem_before) / 1024 / 1024
    
    n_clusters = len(set(labels)) - (1 if -1 in labels else 0) if labels is not None else 0
    
    return {
        'time': elapsed_time,
        'memory': memory_mb,
        'clusters': n_clusters
    }

def get_sample_count(db_file):
    """Extract sample count from filename."""
    basename = os.path.basename(db_file)
    if 'all' in basename.lower():
        return 99999
    
    parts = basename.replace('.db', '').split('_')
    for part in parts:
        clean = part.replace('samples', '').replace('sample', '')
        if clean.isdigit():
            return int(clean)
    return 0


def run_benchmark(db_files, algorithms, use_hamming=False):
    """Benchmark all algorithms across all database sizes."""
    mode = "WITH Hamming" if use_hamming else "WITHOUT Hamming"
    print(f"BENCHMARKING: {mode}")

    
    # Sort databases by size
    db_files_sorted = sorted(db_files, key=get_sample_count)
    results = {algo: {'sizes': [], 'times': [], 'memory': [], 'variants': []} 
               for algo in algorithms}
    
    for db_file in db_files_sorted:
        sample_count = get_sample_count(db_file)
        db_name = os.path.basename(db_file)
        
        print(f"\n{db_name} ({sample_count if sample_count < 99999 else 'ALL'} samples):")
        
        coordinates, variants = load_database(db_file, chromosome='1')
        
        if coordinates is None or len(coordinates) == 0:
            print("  ⚠ No data - skipping")
            continue
        
        n_variants = len(variants)
        print(f"  Chr 1: {n_variants:,} variants")
        
        for algo in algorithms:
            print(f"    {algo}...", end=' ', flush=True)
            
            result = benchmark_algorithm(coordinates, variants, algo, use_hamming)
            
            # Store datapoint for this algorithm
            results[algo]['sizes'].append(sample_count)
            results[algo]['times'].append(result['time'])
            results[algo]['memory'].append(result['memory'])
            results[algo]['variants'].append(n_variants)
            
            print(f"{result['time']:.3f}s | {result['memory']:.1f}MB | {result['clusters']} clusters")
    
    return results


def calculate_big_o(x, y):
    """Calculate Big O complexity from log-log slope."""
    if len(x) < 2:
        return None, None, "N/A"
    
    log_x, log_y = np.log10(x), np.log10(y)
    coeffs = np.polyfit(log_x, log_y, 1)
    slope = coeffs[0]
    
    y_pred = np.polyval(coeffs, log_x)
    ss_res = np.sum((log_y - y_pred) ** 2)
    ss_tot = np.sum((log_y - np.mean(log_y)) ** 2)
    r2 = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0
    
    if slope < 1.15:
        complexity = "O(n)"
    elif slope < 1.5:
        complexity = "O(n log n)"
    elif slope < 2.5:
        complexity = "O(n²)"
    else:
        complexity = "O(n³+)"
    
    return slope, r2, complexity


def plot_results(results_no_hamming, results_with_hamming, algorithms):
    
    # Algorithm colors
    colors = {'DBSCAN': '#E74C3C','OPTICS': '#3498DB','INTERVAL_TREE': '#2ECC71'}
    
    markers = {'DBSCAN': 'o', 'OPTICS': 's','INTERVAL_TREE': '^'
    }
    
    output_files = []
    
    # Create 4 plots: Time (no/with Hamming) + Memory (no/with Hamming)
    plot_configs = [
        ('times', 'Time (seconds)', results_no_hamming, 'WITHOUT Hamming'),
        ('times', 'Time (seconds)', results_with_hamming, 'WITH Hamming'),
        ('memory', 'Memory (MB)', results_no_hamming, 'WITHOUT Hamming'),
        ('memory', 'Memory (MB)', results_with_hamming, 'WITH Hamming'),
    ]
    
    for metric, ylabel, results, mode in plot_configs:
        fig, ax = plt.subplots(figsize=(14, 8))
        
        # Plot each algorithm as a curve
        for algo in algorithms:
            x_data = np.array(results[algo]['sizes'])
            y_data = np.array(results[algo][metric])
            
            if len(x_data) == 0:
                continue
            
            # Calculate Big O
            slope, r2, complexity = calculate_big_o(x_data, y_data)
            label = f"{algo} ({complexity}, slope={slope:.2f}, R²={r2:.3f})" if slope else algo
            
            # Plot the curve: one line connecting all database sizes
            ax.plot(x_data, y_data,
                   color=colors[algo],
                   marker=markers[algo],
                   markersize=12,
                   linewidth=3.5,
                   label=label,
                   alpha=0.9,
                   zorder=3)
        
        # Configure plot
        ax.set_xlabel('Database Size (Number of Samples)', fontsize=15, fontweight='bold')
        ax.set_ylabel(ylabel, fontsize=15, fontweight='bold')
        ax.set_title(f'Algorithm Scalability: {ylabel}\n{mode}',
                     fontsize=17, fontweight='bold', pad=20)
        
        # Log scale on Y for Big O analysis
        ax.set_yscale('log')
        
        # Configure X-axis
        all_sizes = sorted(set([s for algo in algorithms for s in results[algo]['sizes']]))
        ax.set_xticks(all_sizes)
        ax.set_xticklabels([str(s) if s < 99999 else 'ALL' for s in all_sizes], fontsize=12)
        
        ax.grid(True, alpha=0.3, which='both', linestyle='--')
        ax.legend(fontsize=12, frameon=True, shadow=True, loc='best')
        
        # Add explanation box
        explanation = (
            f'X-axis: Database size progression\n'
            f'Y-axis: {ylabel} (log scale)\n'
            f'Each colored curve = one algorithm\n'
            f'Curve shape shows scalability'
        )
        ax.text(0.02, 0.98, explanation,
               transform=ax.transAxes,
               fontsize=10,
               verticalalignment='top',
               bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.6))
        
        plt.tight_layout()
        
        # Save
        metric_name = metric.replace('times', 'time')
        mode_str = 'with_hamming' if 'WITH' in mode else 'no_hamming'
        filename = f"scalability_{metric_name}_{mode_str}.png"
        plt.savefig(filename, dpi=300, bbox_inches='tight')
        plt.close()
        
        output_files.append(filename)
        print(f"  ✓ {filename}")
    
    return output_files


def main():
    # Find databases
    db_files = glob.glob("*.db")
    sample_dbs = [f for f in db_files if 'sample' in f.lower()]
    
    if not sample_dbs:
        print("ERROR: No sample databases found")
        sys.exit(1)
    
    print(f"\nFound databases:")
    for f in sorted(sample_dbs, key=get_sample_count):
        count = get_sample_count(f)
        print(f"  • {f} ({count if count < 99999 else 'ALL'} samples)")
    
    algorithms = ['DBSCAN', 'OPTICS', 'INTERVAL_TREE']
    
    # Run benchmarks
    results_no_hamming = run_benchmark(sample_dbs, algorithms, use_hamming=False)
    results_with_hamming = run_benchmark(sample_dbs, algorithms, use_hamming=True)
    
    # Show what we'll plot
    sample_sizes = sorted(set([s for algo in algorithms 
                               for s in results_no_hamming[algo]['sizes']]))
    

    print("PLOT CONFIGURATION")

    print(f"\n X-axis (database sizes): {[s if s < 99999 else 'ALL' for s in sample_sizes]}")
    print(f" Y-axis: Time (seconds) or Memory (MB)")
    print(f"\n Each plot will show 3 curves (one per algorithm)")
    print(f"   connecting all {len(sample_sizes)} database size points")
    
    # Generate plots
    plot_files = plot_results(results_no_hamming, results_with_hamming, algorithms)

    print(f"COMPLETE - Generated {len(plot_files)} scalability plots")

    
if __name__ == "__main__":
    main()