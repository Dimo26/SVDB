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



def benchmark_algorithm(coordinates, variants, algorithm):
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


def run_benchmark(db_files, algorithms):
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
            result = benchmark_algorithm(coordinates, variants, algo)
            
            results[algo]['sizes'].append(sample_count)
            results[algo]['times'].append(result['time'])
            results[algo]['memory'].append(result['memory'])
            results[algo]['variants'].append(n_variants)
            results[algo]['clusters'].append(result['clusters'])
    
    return results


def plot_scalability_curves(results, algorithms):
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
    
    output_files = []
    
    plot_configs = [
        ('times', 'Time (seconds)'),
        ('memory', 'Memory (MB)'),
    ]
    
    for metric, ylabel in plot_configs:
        fig, ax = plt.subplots(figsize=(12, 7))
        
        algo_list = list(algorithms) if isinstance(algorithms, set) else algorithms
        for algo in algo_list:
            x_data = np.array(results[algo]['sizes'])
            y_data = np.array(results[algo][metric])
            
            if len(x_data) == 0:
                continue
            
            ax.plot(x_data, y_data,
                   color=colors[algo],
                   marker=markers[algo],
                   markersize=10,
                   linewidth=3,
                   label=algo,
                   alpha=0.95)
        
        ax.set_xlabel('Database Size (Number of Samples)', fontsize=13, fontweight='bold')
        ax.set_ylabel(ylabel, fontsize=13, fontweight='bold')
        ax.set_yscale('log')
        
        all_sizes = sorted(set([s for algo in algorithms for s in results[algo]['sizes']]))
        if all_sizes:
            ax.set_xticks(all_sizes)
            ax.set_xticklabels([str(s) for s in all_sizes], fontsize=11)
        ax.legend(fontsize=10, frameon=True, loc='best')
        
        plt.tight_layout()
        
        metric_name = metric.replace('times', 'time')
        filename = f"scalability_{metric_name}.png"
        plt.savefig(filename, dpi=300, bbox_inches='tight')
        plt.close()
        
        output_files.append(filename)
    
    return output_files



def print_summary(results, algorithms):
    """Print cluster & variant counts per algorithm and per DB (no plotting)."""
    print("\nBenchmark summary:")
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
    
    results = run_benchmark(sample_dbs, algorithms)
    
    # print cluster/variant counts (simple textual stats)
    print_summary(results, algorithms)
    
    plot_files = plot_scalability_curves(results, algorithms)


if __name__ == "__main__":
    main()
    