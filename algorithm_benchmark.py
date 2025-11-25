import time
import numpy as np
import matplotlib.pyplot as plt
from memory_profiler import profile  # You'll need: pip install memory-profiler --break-system-packages
from sklearn.metrics import silhouette_score, adjusted_rand_score
import pandas as pd
import gzip 
from svdb.readVCF import readVCFLine

# TODO: Import algos once youve written them
# from dbscan_improved import DBSCAN_Cluster
# from optics_clustering import OPTICS_Cluster
# from interval_tree_overlap import IntervalTree_Overlap
# from rtree_spatial import RTree_Spatial
# from graph_clustering import ConnectedComponents_Cluster



def generate_test_sv_data(n_variants=1000, n_clusters=5, noise_ratio=0.1):

    coordinates = []
    true_labels = []

    cluster_centers = np.random.randint(1000000, 10000000, size=(n_clusters, 2))
    # 
    for cluster_id in range(n_clusters):
        center = cluster_centers[cluster_id]
        n_points = (n_variants * (1 - noise_ratio)) // n_clusters
        cluster_points = np.random.normal(loc=center, scale=5000, size=(n_points, 2))
        coordinates.extend(cluster_points)
        true_labels.extend([cluster_id] * n_points)
    n_noise = int(n_variants * noise_ratio)
    noise_points = np.random.randint(1000000, 10000000, size=(n_noise, 2))
    coordinates.extend(noise_points)
    true_labels.extend([-1] * n_noise)
    return np.array(coordinates), np.array(true_labels)

def load_real_vcf_data(vcf_path, max_variants=None):
    coordinates = []

    chr_types = {'chrA': [], 'chrB': [], 'variant_type': []}
    opener = gzip.open if vcf_path.endswith('.gz') else open
    with opener(vcf_path, 'rt') as f:
         for line in f:
              if line.startswith('#'):
                  continue
              record = readVCFLine(line)
    posA = record['POS']
    posB = record['INFO'].get('END', posA)
    chrA = record['CHROM']
    chrB = record['INFO'].get('CHR2', chrA)
    sv_type = record['INFO'].get('SVTYPE', 'UNK')
    if posB >= posA:
        coordinates.append([posA, posB])
        chr_types['chrA'].append(chrA)
        chr_types['chrB'].append(chrB)
        chr_types['variant_type'].append(sv_type)
    if max_variants and len(coordinates) >= max_variants:
        break

    return np.array(coordinates), chr_types
    pass  # Remove this when you implement



class AlgorithmBenchmark:
    def __init__(self, coordinates, true_labels=None):
        self.coordinates = coordinates
        self.true_labels = true_labels
        self.results = {}
    
    def benchmark_algorithm(self, algorithm_name, algorithm_func, **algorithm_params):
        start_time = time.time()
        try:
              labels = algorithm_func(self.coordinates, **algorithm_params)
              runtime = time.time() - start_time
       
              # Calculate metrics
              n_clusters = len(set(labels)) - (1 if -1 in labels else 0)

              if n_clusters >= 2:
                  silhouette = silhouette_score(self.coordinates, labels)
              else:
                  silhouette = -1
             
              if self.true_labels is not None:
                  ari = adjusted_rand_score(self.true_labels, labels)
              else:
                 ari = None
             
              self.results[algorithm_name] = {
                  'runtime': runtime,
                 'n_clusters': n_clusters,
                 'silhouette': silhouette,
                  'ari': ari,
                 'labels': labels,
                  'params': algorithm_params
             }
            
        except Exception as e:
              print(f"Algorithm {algorithm_name} failed: {e}")
              self.results[algorithm_name] = {'error': str(e)}
        
        pass  # Remove this when you implement
    
    def run_all_benchmarks(self, algorithms_dict):
        for name, (func, params) in algorithms_dict.items():
             print(f"Benchmarking {name}...")
             self.benchmark_algorithm(name, func, **params)
        
        pass  # Remove this when you implement
    
    def get_results_dataframe(self):

        data = []
        for alg_name, metrics in self.results.items():
              if 'error' not in metrics:
                  data.append({
                      'Algorithm': alg_name,
                      'Runtime (s)': metrics['runtime'],
                      'N_Clusters': metrics['n_clusters'],
                      'Silhouette': metrics['silhouette'],
                      'ARI': metrics['ari']
                  })        
        df = pd.DataFrame(data)
        return df.sort_values('Runtime (s)')
        
        pass  # Remove this when you implement

def plot_runtime_comparison(results_df, output_path='runtime_comparison.png'):
    output_path = 'Assets/' + output_path
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.bar(results_df['Algorithm'], results_df['Runtime (s)'])
    ax.set_xlabel('Algorithm')
    ax.set_ylabel('Runtime (seconds)')
    ax.set_title('Algorithm Runtime Comparison')
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    print(f"Saved runtime comparison to {output_path}")
    
    pass  # Remove this when you implement

def plot_quality_metrics(results_df, output_path='quality_metrics.png'):
    output_path = 'Assets/' + output_path
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
    ax1.bar(results_df['Algorithm'], results_df['Silhouette'], color='skyblue')
    ax1.set_xlabel('Algorithm')
    ax1.set_ylabel('Silhouette Score')
    ax1.set_title('Silhouette Score Comparison')
    ax1.tick_params(axis='x', rotation=45)

    ax2.bar(results_df['Algorithm'], results_df['ARI'], color='lightgreen')
    ax2.set_xlabel('Algorithm')
    ax2.set_ylabel('Adjusted Rand Index (ARI)')
    ax2.set_title('Adjusted Rand Index Comparison')
    ax2.tick_params(axis='x', rotation=45)

    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    print(f"Saved quality metrics comparison to {output_path}")



def plot_scalability_test(n_variants_list, results_dict, output_path='scalability.png'):
    """
    Test how algorithms scale with increasing data size
    
    WHAT YOU NEED TO DO:
    1. Run each algorithm with increasing dataset sizes
    2. Plot runtime vs. dataset size for each algorithm
    3. This shows which algorithms scale better (O(n), O(n log n), O(n²))
    
    Parameters:
    -----------
    n_variants_list : list
        List of dataset sizes to test (e.g., [100, 500, 1000, 5000, 10000])
    results_dict : dict
        Dictionary mapping algorithm name to list of runtimes
    output_path : str
        Where to save the plot
    
    IMPLEMENTATION HINTS:
    - Use plt.plot() for line plots
    - Different color for each algorithm
    - Log scale might be useful: plt.yscale('log')
    - Add legend with plt.legend()
    """
    
    # TODO: YOUR CODE HERE
    fig, ax = plt.subplots(figsize=(10, 6))
    for alg_name, runtimes in results_dict.items():
         ax.plot(n_variants_list, runtimes, marker='o', label=alg_name)
         ax.set_yscale('log')
         ax.set_xlabel('Number of Variants')
         ax.set_ylabel('Runtime (seconds)')
         ax.set_title('Algorithm Scalability')
         ax.legend()
         plt.savefig(output_path, dpi=300)
    
    pass  # Remove this when you implement


def plot_clustering_results(coordinates, labels, algorithm_name, output_path=None):
    """
    Visualize clustering results as scatter plot
    
    WHAT YOU NEED TO DO:
    1. Create scatter plot with posA on x-axis, posB on y-axis
    2. Color points by cluster label
    3. Mark noise points (label -1) differently
    
    IMPLEMENTATION HINTS:
    - Use plt.scatter() with c=labels for coloring
    - Use different marker for noise: labels == -1
    - Add colorbar if many clusters
    """
    
    # TODO: YOUR CODE HERE
    plt.figure(figsize=(8, 8))
    plt.scatter(coordinates[:, 0], coordinates[:, 1], c=labels, cmap='tab20', s=10)
    plt.xlabel('Position A')
    plt.ylabel('Position B')
    plt.title(f'Clustering Results: {algorithm_name}')
    if output_path:
         plt.savefig(output_path, dpi=300)
         print(f"Saved clustering plot to {output_path}")
    pass  # Remove this when you implement

def main():

    print("\n1. Generating test data...")
    coordinates, true_labels = generate_test_sv_data(n_variants=1000, n_clusters=5)
    print(f"   Generated {len(coordinates)} variants with {len(set(true_labels))} clusters")

    print("\n2. Setting up algorithms...")
    algorithms = {
         #'DBSCAN': (dbscan_function, {'epsilon': 1000, 'min_pts': 2}),
         #'OPTICS': (optics_function, {'min_samples': 2, 'max_eps': 2000}),
         #'IntervalTree': (interval_tree_function, {'max_distance': 1000}),
         'RTree': (rtree_function, {'max_distance': 1000})}
    
    print("\n3. Running benchmarks...")
    benchmark = AlgorithmBenchmark(coordinates, true_labels)
    benchmark.run_all_benchmarks(algorithms)
    print("\n4. Analyzing results...")
    results_df = benchmark.get_results_dataframe()
    print("\nResults Summary:")
    print(results_df)
    print("\n5. Generating visualizations...")
    plot_runtime_comparison(results_df)
    plot_quality_metrics(results_df)
    print("\n6. Testing scalability...")

n_variants_list = [100, 500, 1000, 2000, 5000]
scalability_results = {}
for alg_name, (alg_func, params) in algorithms.items():
         runtimes = []
         for n in n_variants_list:
            test_data, _ = generate_test_sv_data(n_variants=n)

         scalability_results[alg_name] = runtimes
plot_scalability_test(n_variants_list, scalability_results)
print("\nBenchmarking complete.")

if __name__ == "__main__":
    main()