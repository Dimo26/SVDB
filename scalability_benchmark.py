#!/usr/bin/env python3
import sys
import os
import time
import glob
import psutil
import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0, os.getcwd())

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
    """Calculate normalized Hamming distance."""
    if not seq1 or not seq2:
        return 1.0
    seq1, seq2 = str(seq1).upper(), str(seq2).upper()
    max_len = max(len(seq1), len(seq2))
    min_len = min(len(seq1), len(seq2))
    mismatches = sum(1 for i in range(min_len) if seq1[i] != seq2[i])
    return (mismatches + (max_len - min_len)) / max_len if max_len > 0 else 0.0


def apply_hamming_reclustering(labels, variants, max_hamming=0.2):
    """
    Re-cluster INS variants within spatial clusters using Hamming distance.
    
    SVDB Clustering Rules:
    1. First: Spatial clustering by position (DBSCAN/OPTICS/INTERVAL_TREE)
    2. Then: For INS variants in each spatial cluster, re-cluster by sequence similarity
    3. Non-INS variants keep their spatial cluster assignment
    """
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
                    
                    if hamming_distance(seq_i, seq_j) <= max_hamming:
                        group.append(j)
                        assigned.add(j)
                
                for g in group:
                    new_labels[ins_with_seq[g]] = next_cluster_id
                next_cluster_id += 1
    
    return new_labels


def benchmark_algorithm(coordinates, variants, algorithm, use_hamming=False):
    """Run clustering algorithm and measure performance."""
    process = psutil.Process()
    mem_before = process.memory_info().rss
    
    start_time = time.time()
    
    if algorithm == 'DBSCAN':
        labels = DBSCAN.cluster(coordinates, epsilon=500, m=2)
    elif algorithm == 'OPTICS':
        labels = optics_cluster(coordinates, min_samples=2, max_eps=500)
    elif algorithm == 'INTERVAL_TREE':
        labels = interval_tree_cluster(coordinates, max_distance=500)
    else:
        labels = None
    
    if use_hamming and labels is not None:
        labels = apply_hamming_reclustering(labels, variants, max_hamming=0.2)
    
    elapsed_time = time.time() - start_time
    mem_after = process.memory_info().rss
    memory_mb = (mem_after - mem_before) / 1024 / 1024
    
    n_clusters = 0
    if labels is not None:
        unique = set(labels.tolist())
        n_clusters = len(unique) - (1 if -1 in unique else 0)
    
    return {
        'time': elapsed_time,
        'memory': memory_mb,
        'clusters': n_clusters
    }
''

def get_sample_count(db_file):
    """Extract sample count from filename."""
    basename = os.path.basename(db_file)
    parts = basename.replace('.db', '').split('_')
    for part in parts:
        clean = part.replace('samples', '').replace('sample', '')
        if clean.isdigit():
            return int(clean)
    return 0


def run_benchmark(db_files, algorithms, use_hamming=False):
    """Benchmark all algorithms across database sizes."""
    db_files_sorted = sorted(db_files, key=get_sample_count)
    
    results = {algo: {'sizes': [], 'times': [], 'memory': [], 'variants': [], 'clusters': []} 
               for algo in algorithms}
    
    for db_file in db_files_sorted:
        sample_count = get_sample_count(db_file)
        
        coordinates, variants = load_database(db_file, chromosome='1')
        
        if coordinates is None or len(coordinates) == 0:
            continue
        
        n_variants = len(variants)
        
        for algo in algorithms:
            result = benchmark_algorithm(coordinates, variants, algo, use_hamming)
            
            results[algo]['sizes'].append(sample_count)
            results[algo]['times'].append(result['time'])
            results[algo]['memory'].append(result['memory'])
            results[algo]['variants'].append(n_variants)
            results[algo]['clusters'].append(result['clusters'])
    
    return results


def plot_scalability_curves(results_no_hamming, results_with_hamming, algorithms):
    """Generate scalability curve plots."""
    colors = {
        'DBSCAN': '#E74C3C',
        'OPTICS': '#3498DB',
        'INTERVAL_TREE': '#2ECC71'
    }
    
    markers = {
        'DBSCAN': 'o',
        'OPTICS': 's',
        'INTERVAL_TREE': '^'
    }
    
    descriptions = {
        'DBSCAN': 'DBSCAN',
        'OPTICS': 'OPTICS',
        'INTERVAL_TREE': 'Interval Tree'
    }
    
    output_files = []
    
    plot_configs = [
        ('times', 'Time (seconds)', results_no_hamming, 'no_hamming'),
        ('times', 'Time (seconds)', results_with_hamming, 'with_hamming'),
        ('memory', 'Memory (MB)', results_no_hamming, 'no_hamming'),
        ('memory', 'Memory (MB)', results_with_hamming, 'with_hamming'),
    ]
    
    for metric, ylabel, results, mode in plot_configs:
        fig, ax = plt.subplots(figsize=(12, 7))
        
        hamming_label = "WITH Hamming" if 'with' in mode else "WITHOUT Hamming"
        
        algo_list = list(algorithms) if isinstance(algorithms, set) else algorithms
        for algo in algo_list:
            x_data = np.array(results[algo]['sizes'])
            y_data = np.array(results[algo][metric])
            
            if len(x_data) == 0:
                continue
            
            label = f"{algo}: {descriptions.get(algo, algo)} ({hamming_label})"
            
            ax.plot(x_data, y_data,
                   color=colors[algo],
                   marker=markers[algo],
                   markersize=10,
                   linewidth=3,
                   label=label,
                   alpha=0.95)
        
        ax.set_xlabel('Database Size (Number of Samples)', fontsize=13, fontweight='bold')
        ax.set_ylabel(ylabel, fontsize=13, fontweight='bold')
        ax.set_yscale('log')
        
        all_sizes = sorted(set([s for algo in algorithms for s in results[algo]['sizes']]))
        if all_sizes:
            ax.set_xticks(all_sizes)
            ax.set_xticklabels([str(s) for s in all_sizes], fontsize=11)
        ax.legend(fontsize=11, frameon=True, loc='best')
        
        plt.tight_layout()
        
        metric_name = metric.replace('times', 'time')
        filename = f"scalability_{metric_name}_{mode}.png"
        plt.savefig(filename, dpi=300, bbox_inches='tight')
        plt.close()
        
        output_files.append(filename)
    
    return output_files


def plot_1000_sample_comparison(results_no_hamming, results_with_hamming, algorithms):
    """
    Generate bar plots comparing WITH/WITHOUT Hamming for 1000 sample database.
    Two plots: Time comparison and Memory comparison.
    """
    colors = {
        'DBSCAN': '#E74C3C',
        'OPTICS': '#3498DB',
        'INTERVAL_TREE': '#2ECC71'
    }
    
    target_size = 1000
    output_files = []
    
    time_data = {'no_hamming': [], 'with_hamming': []}
    memory_data = {'no_hamming': [], 'with_hamming': []}
    
    for algo in algorithms:
        for idx, size in enumerate(results_no_hamming[algo]['sizes']):
            if size == target_size:
                time_data['no_hamming'].append(results_no_hamming[algo]['times'][idx])
                memory_data['no_hamming'].append(results_no_hamming[algo]['memory'][idx])
                break
        else:
            time_data['no_hamming'].append(0)
            memory_data['no_hamming'].append(0)
        
        for idx, size in enumerate(results_with_hamming[algo]['sizes']):
            if size == target_size:
                time_data['with_hamming'].append(results_with_hamming[algo]['times'][idx])
                memory_data['with_hamming'].append(results_with_hamming[algo]['memory'][idx])
                break
        else:
            time_data['with_hamming'].append(0)
            memory_data['with_hamming'].append(0)
    
    x = np.arange(len(algorithms))
    width = 0.35
    
    fig, ax = plt.subplots(figsize=(12, 7))
    
    bars1 = ax.bar(x - width/2, time_data['no_hamming'], width, 
                   label='WITHOUT Hamming', alpha=0.9,
                   color=[colors[algo] for algo in algorithms],
                   edgecolor='black', linewidth=1.5)
    
    bars2 = ax.bar(x + width/2, time_data['with_hamming'], width,
                   label='WITH Hamming', alpha=0.6,
                   color=[colors[algo] for algo in algorithms],
                   edgecolor='black', linewidth=1.5, hatch='//')
    
    ax.set_xlabel('Algorithm', fontsize=14, fontweight='bold')
    ax.set_ylabel('Time (seconds)', fontsize=14, fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels(algorithms, fontsize=12)
    ax.legend(fontsize=12, frameon=True)

    plt.tight_layout()
    filename = 'scalability_time_comparison_1000samples.png'
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.close()
    output_files.append(filename)
    
    fig, ax = plt.subplots(figsize=(12, 7))
    
    bars1 = ax.bar(x - width/2, memory_data['no_hamming'], width,
                   label='WITHOUT Hamming', alpha=0.9,
                   color=[colors[algo] for algo in algorithms],
                   edgecolor='black', linewidth=1.5)
    
    bars2 = ax.bar(x + width/2, memory_data['with_hamming'], width,
                   label='WITH Hamming', alpha=0.6,
                   color=[colors[algo] for algo in algorithms],
                   edgecolor='black', linewidth=1.5, hatch='//')
    
    ax.set_xlabel('Algorithm', fontsize=14, fontweight='bold')
    ax.set_ylabel('Memory (MB)', fontsize=14, fontweight='bold')

    ax.set_xticks(x)
    ax.set_xticklabels(algorithms, fontsize=12)
    ax.legend(fontsize=12, frameon=True)
    
    plt.tight_layout()
    filename = 'scalability_memory_comparison_1000samples.png'
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.close()
    output_files.append(filename)
    
    return output_files


def print_summary(results_no_hamming, results_with_hamming, algorithms):
    """Print cluster & variant counts per algorithm and per DB (no plotting)."""
    print("\nBenchmark summary:")
    for mode_label, results in [('WITHOUT Hamming', results_no_hamming), ('WITH Hamming', results_with_hamming)]:
        print(f"\n{mode_label}:")
        for algo in algorithms:
            sizes = results[algo]['sizes']
            variants = results[algo]['variants']
            clusters = results[algo]['clusters']
            print(f"  {algo}:")
            if not sizes:
                print("    (no data)")
                continue
            for s, v, c in zip(sizes, variants, clusters):
                print(f"    Size {s}: Variants={v}, Clusters={c}")
            print(f"    Total variants (sum across DBs)={sum(variants)}, Total clusters={sum(clusters)}")


def main():
    db_files = glob.glob("*.db")
    
    sample_dbs = [f for f in db_files 
                  if 'sample' in f.lower() and 'all' not in f.lower()]
    
    if not sample_dbs:
        sys.exit(1)
    
    algorithms = ['DBSCAN', 'OPTICS', 'INTERVAL_TREE']
    
    results_no_hamming = run_benchmark(sample_dbs, algorithms, use_hamming=False)
    results_with_hamming = run_benchmark(sample_dbs, algorithms, use_hamming=True)
    
    # print cluster/variant counts (simple textual stats)
    print_summary(results_no_hamming, results_with_hamming, algorithms)
    
    plot_files = plot_scalability_curves(results_no_hamming, results_with_hamming, algorithms)
    
    comparison_files = plot_1000_sample_comparison(results_no_hamming, results_with_hamming, algorithms)
    
    all_files = plot_files + comparison_files


if __name__ == "__main__":
    main()
    