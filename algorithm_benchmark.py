
#!/usr/bin/env python3
"""
Algorithm Benchmark for SVDB
Compares DBSCAN, OPTICS, and Interval Tree clustering algorithms
"""

import sys
import os
import time
import psutil
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict

# Ensure current directory is in path
sys.path.insert(0, os.getcwd())

# Import SVDB modules
try:
    from svdb.database import DB
    from svdb.export_module import DBSCAN
    from svdb.optics_clustering import optics_cluster
    from svdb.interval_tree_overlap import interval_tree_cluster
except ImportError as e:
    print(f"Error importing SVDB modules: {e}")
    print("Make sure you're running from the SVDB directory")
    sys.exit(1)


class BenchmarkRunner:
    def __init__(self, db_path):
        self.db_path = db_path
        self.db = DB(db_path, memory=False)
        self.results = defaultdict(dict)
        
    def get_variant_data(self, variant_type, chrA, chrB):
        """Fetch variant coordinates from database"""
        query = f"""
            SELECT posA, posB FROM SVDB 
            WHERE var = '{variant_type}' 
            AND chrA = '{chrA}' 
            AND chrB = '{chrB}'
        """
        hits = self.db.query(query)
        if not hits:
            return None
        
        coordinates = np.array([[int(h[0]), int(h[1])] for h in hits])
        return coordinates
    
    def benchmark_algorithm(self, algorithm_name, coordinates, **params):
        """Benchmark a single algorithm"""
        print(f"  Testing {algorithm_name}...", end='', flush=True)
        
        # Memory before
        process = psutil.Process()
        mem_before = process.memory_info().rss / 1024 / 1024  # MB
        
        # Time execution
        start_time = time.time()
        
        try:
            if algorithm_name == 'DBSCAN':
                labels = DBSCAN.cluster(
                    coordinates,
                    epsilon=params.get('epsilon', 500),
                    m=params.get('min_pts', 2)
                )
            elif algorithm_name == 'OPTICS':
                labels = optics_cluster(
                    coordinates,
                    min_samples=params.get('min_samples', 2),
                    max_eps=params.get('max_eps', 2000)
                )
            elif algorithm_name == 'INTERVAL_TREE':
                labels = interval_tree_cluster(
                    coordinates,
                    max_distance=params.get('max_distance', 1000)
                )
            else:
                raise ValueError(f"Unknown algorithm: {algorithm_name}")
                
        except Exception as e:
            print(f" FAILED ({e})")
            return None
        
        end_time = time.time()
        
        # Memory after
        mem_after = process.memory_info().rss / 1024 / 1024  # MB
        
        # Calculate metrics
        runtime = end_time - start_time
        memory_used = mem_after - mem_before
        
        n_clusters = len(set(labels)) - (1 if -1 in labels else 0)
        n_noise = list(labels).count(-1)
        
        print(f" {runtime:.3f}s, {n_clusters} clusters, {n_noise} noise")
        
        return {
            'runtime': runtime,
            'memory': memory_used,
            'n_clusters': n_clusters,
            'n_noise': n_noise,
            'labels': labels
        }
    
    def run_benchmarks(self):
        """Run benchmarks on all variant types"""
        print("\nFetching variant types...")
        
        # Get unique variant types and chromosomes
        var_query = "SELECT DISTINCT var, chrA, chrB FROM SVDB"
        variants = self.db.query(var_query)
        
        if not variants:
            print("No variants found in database")
            return
        
        print(f"Found {len(variants)} variant type/chromosome combinations\n")
        
        # Test each variant type
        for var_type, chrA, chrB in variants:
            print(f"\n{var_type} on {chrA}-{chrB}:")
            
            coordinates = self.get_variant_data(var_type, chrA, chrB)
            if coordinates is None or len(coordinates) < 2:
                print(f"  Skipped (insufficient data: {len(coordinates) if coordinates is not None else 0} variants)")
                continue
            
            print(f"  Dataset size: {len(coordinates)} variants")
            
            # Test all algorithms
            algorithms = {
                'DBSCAN': {'epsilon': 500, 'min_pts': 2},
                'OPTICS': {'min_samples': 2, 'max_eps': 2000},
                'INTERVAL_TREE': {'max_distance': 1000}
            }
            
            for alg_name, params in algorithms.items():
                result = self.benchmark_algorithm(alg_name, coordinates, **params)
                
                if result:
                    key = f"{var_type}_{chrA}_{chrB}"
                    self.results[key][alg_name] = result
    
    def plot_results(self):
        """Generate comparison plots"""
        if not self.results:
            print("\nNo results to plot")
            return
        
        print("\nGenerating plots...")
        
        # Prepare data for plotting
        algorithms = ['DBSCAN', 'OPTICS', 'INTERVAL_TREE']
        variant_keys = list(self.results.keys())
        
        runtimes = {alg: [] for alg in algorithms}
        memories = {alg: [] for alg in algorithms}
        clusters = {alg: [] for alg in algorithms}
        
        for key in variant_keys:
            for alg in algorithms:
                if alg in self.results[key]:
                    runtimes[alg].append(self.results[key][alg]['runtime'])
                    memories[alg].append(self.results[key][alg]['memory'])
                    clusters[alg].append(self.results[key][alg]['n_clusters'])
                else:
                    runtimes[alg].append(np.nan)
                    memories[alg].append(np.nan)
                    clusters[alg].append(np.nan)
        
        # Create figure with subplots
        fig, axes = plt.subplots(1, 3, figsize=(15, 5))
        fig.suptitle('SVDB Clustering Algorithm Comparison', fontsize=14, fontweight='bold')
        
        x_pos = np.arange(len(variant_keys))
        width = 0.25
        
        colors = {'DBSCAN': '#1f77b4', 'OPTICS': '#ff7f0e', 'INTERVAL_TREE': '#2ca02c'}
        
        # Runtime plot
        for i, alg in enumerate(algorithms):
            axes[0].bar(x_pos + i*width, runtimes[alg], width, 
                       label=alg, color=colors[alg], alpha=0.8)
        axes[0].set_xlabel('Variant Type')
        axes[0].set_ylabel('Runtime (seconds)')
        axes[0].set_title('Execution Time')
        axes[0].set_xticks(x_pos + width)
        axes[0].set_xticklabels([k.replace('_', '\n') for k in variant_keys], rotation=45, ha='right')
        axes[0].legend()
        axes[0].grid(axis='y', alpha=0.3)
        
        # Memory plot
        for i, alg in enumerate(algorithms):
            axes[1].bar(x_pos + i*width, memories[alg], width,
                       label=alg, color=colors[alg], alpha=0.8)
        axes[1].set_xlabel('Variant Type')
        axes[1].set_ylabel('Memory Usage (MB)')
        axes[1].set_title('Memory Consumption')
        axes[1].set_xticks(x_pos + width)
        axes[1].set_xticklabels([k.replace('_', '\n') for k in variant_keys], rotation=45, ha='right')
        axes[1].legend()
        axes[1].grid(axis='y', alpha=0.3)
        
        # Clusters plot
        for i, alg in enumerate(algorithms):
            axes[2].bar(x_pos + i*width, clusters[alg], width,
                       label=alg, color=colors[alg], alpha=0.8)
        axes[2].set_xlabel('Variant Type')
        axes[2].set_ylabel('Number of Clusters')
        axes[2].set_title('Clustering Results')
        axes[2].set_xticks(x_pos + width)
        axes[2].set_xticklabels([k.replace('_', '\n') for k in variant_keys], rotation=45, ha='right')
        axes[2].legend()
        axes[2].grid(axis='y', alpha=0.3)
        
        plt.tight_layout()
        
        # Save plot
        output_file = 'algorithm_comparison.png'
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"✓ Saved plot: {output_file}")
        plt.close()
        
        # Print summary statistics
        self.print_summary()
    
    def print_summary(self):
        """Print summary statistics"""
        print("\n" + "="*60)
        print("BENCHMARK SUMMARY")
        print("="*60)
        
        algorithms = ['DBSCAN', 'OPTICS', 'INTERVAL_TREE']
        
        for alg in algorithms:
            runtimes = []
            memories = []
            
            for key in self.results:
                if alg in self.results[key]:
                    runtimes.append(self.results[key][alg]['runtime'])
                    memories.append(self.results[key][alg]['memory'])
            
            if runtimes:
                print(f"\n{alg}:")
                print(f"  Avg Runtime: {np.mean(runtimes):.3f}s (±{np.std(runtimes):.3f}s)")
                print(f"  Avg Memory:  {np.mean(memories):.2f}MB (±{np.std(memories):.2f}MB)")
                print(f"  Tests run:   {len(runtimes)}")


def main():
    # Find database file
    db_files = [f for f in os.listdir('.') if f.endswith('.db')]
    
    if not db_files:
        print("Error: No database files found")
        sys.exit(1)
    
    db_path = db_files[0]
    print(f"Using database: {db_path}")
    
    # Run benchmark
    runner = BenchmarkRunner(db_path)
    runner.run_benchmarks()
    runner.plot_results()
    
    print("\n✓ Benchmark complete")


if __name__ == '__main__':
    main()