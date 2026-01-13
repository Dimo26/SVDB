#!/usr/bin/env python3
"""
Algorithm Benchmark for SVDB Clustering
Compares: DBSCAN, OPTICS, INTERVAL_TREE, and overlap-based methods for report 
With and without Hamming distance post-processing for insertions which helps long reads insertions clustering. 
"""

import time
import sys
import os
import glob
import numpy as np
import psutil

def load_db_samples(db_file):
    """Load variants from SVDB database."""
    try:
        from svdb import database
        db = database.DB(db_file)
        
        coordinates = []
        variants = []
        
        for row in db.query('SELECT * FROM SVDB WHERE chrA = "1" AND chrB = "1"'):
            var_type = row[0]
            chrA = row[1]  
            chrB = row[2]
            posA = int(row[3])
            posB = int(row[6])
            seq = row[11] if len(row) > 11 else ''
            
            coordinates.append([posA, posB])
            variants.append({
                'posA': posA,
                'posB': posB,
                'type': var_type,
                'sequence': seq
            })
        
        return (np.array(coordinates) if coordinates else None, variants)
    except Exception as e:
        print(f"Error loading database {db_file}: {e}")
        return None, None


def _hamming_distance(seq1, seq2):
    """Calculate normalized Hamming distance."""
    if seq1 is None or seq2 is None or len(seq1) == 0 or len(seq2) == 0:
        return 1.0
    seq1 = str(seq1).upper()
    seq2 = str(seq2).upper()
    min_len = min(len(seq1), len(seq2))
    max_len = max(len(seq1), len(seq2))
    mismatches = sum(1 for i in range(min_len) if seq1[i] != seq2[i])
    length_diff = max_len - min_len
    total_dist = mismatches + length_diff
    return total_dist / max_len if max_len > 0 else 0.0


def apply_hamming_post_processing(labels, variants, max_hamming=0.2):
    """Apply Hamming distance to insertions AFTER spatial."""
    if labels is None or variants is None:
        return labels
    
    labels = np.array(labels)
    new_labels = np.full_like(labels, -1)
    next_cluster_id = 0
    
    unique_spatial_labels = sorted(set(labels.tolist()))
    
    for spatial_label in unique_spatial_labels:
        if spatial_label == -1:  # Noise
            continue
        
        indices = np.where(labels == spatial_label)[0]
        if len(indices) == 0:
            continue
        
        # Separate insertions with sequences from others
        ins_with_seq = []
        other_indices = []
        
        for idx in indices:
            if variants[idx]['type'] == 'INS' and variants[idx]['sequence']:
                ins_with_seq.append(idx)
            else:
                other_indices.append(idx)
        
        # Non-insertions keep their spatial cluster
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
                    dist = _hamming_distance(seq_i, variants[idx_j]['sequence'])
                    if dist <= max_hamming:
                        group.append(j)
                        assigned.add(j)
                
                for g in group:
                    new_labels[ins_with_seq[g]] = next_cluster_id
                next_cluster_id += 1
    
    return new_labels


def benchmark_clustering_algorithm(coordinates, variants, algorithm_name, apply_hamming=False, max_hamming=0.2):
    """
    Benchmark a clustering algorithm with optional Hamming post-processing.
    Returns the time, mem used and amount of clusters as well as labels 
    """
    from svdb.export_module import DBSCAN
    from svdb.optics_clustering import OPTICS
    from svdb.interval_tree_overlap import interval_tree_cluster
    from svdb import overlap_module
    
    if coordinates is None or len(coordinates) == 0:
        return None, None, 0, None
    
    distance_threshold = 10000
    proc = psutil.Process()
    
    # Measure spatial clustering
    mem_before = proc.memory_info().rss
    start_time = time.time()
    
    try:
        if algorithm_name == 'DBSCAN':
            labels = DBSCAN.cluster(coordinates, distance_threshold, 3)
        elif algorithm_name == 'OPTICS':
            optics = OPTICS(min_samples=3, max_eps=distance_threshold)
            labels = optics.fit_predict(coordinates)
        elif algorithm_name == 'INTERVAL_TREE':
            labels = interval_tree_cluster(coordinates, distance_threshold)
        elif algorithm_name == 'OVERLAP':
            # Overlap-based clustering (baseline)
            labels = np.full(len(coordinates), -1)
            # Simple overlap clustering logic
            cluster_id = 0
            for i in range(len(coordinates)):
                if labels[i] != -1:
                    continue
                labels[i] = cluster_id
                for j in range(i+1, len(coordinates)):
                    if labels[j] == -1:
                        dist = max(abs(coordinates[i][0] - coordinates[j][0]), 
                                  abs(coordinates[i][1] - coordinates[j][1]))
                        if dist <= distance_threshold:
                            labels[j] = cluster_id
                cluster_id += 1
        else:
            return None, None, 0, None
    except Exception as e:
        print(f"Error in {algorithm_name}: {e}")
        return None, None, 0, None
    
    spatial_time = time.time() - start_time
    mem_after_spatial = proc.memory_info().rss
    spatial_memory = mem_after_spatial - mem_before
    
    # Apply Hamming if requested
    hamming_time = 0
    hamming_memory = 0
    if apply_hamming:
        mem_before_hamming = proc.memory_info().rss
        hamming_start = time.time()
        labels = apply_hamming_post_processing(labels, variants, max_hamming)
        hamming_time = time.time() - hamming_start
        mem_after_hamming = proc.memory_info().rss
        hamming_memory = mem_after_hamming - mem_before_hamming
    
    total_time = spatial_time + hamming_time
    total_memory = spatial_memory + hamming_memory
    n_clusters = len(set(labels)) - (1 if -1 in labels else 0)
    
    return total_time, total_memory, n_clusters, labels


def main():

    print("SVDB Clustering Algorithm Benchmark")
    print("Benchmarking on CHROMOSOME 1 only")


    
    # Find database files
    db_files = glob.glob('*.db') + glob.glob('../*.db')
    
    if not db_files:
        print("No .db files found. Provide SVDB database files.")
        sys.exit(1)
    
    print(f"\nFound {len(db_files)} database file(s)")
    
    algorithms = ['DBSCAN', 'OPTICS', 'INTERVAL_TREE', 'OVERLAP']
    results = {}
    
    for db_file in db_files[:3]:  # Limit to first 3 files
        print(f"\n{'='*60}")
        print(f"Benchmarking: {os.path.basename(db_file)}")
        print(f"{'='*60}")

        coordinates, variants = load_db_samples(db_file)
        
        if coordinates is None or len(coordinates) == 0:
            print(f"  Skipped (no valid variants on chromosome 1)")
            continue
        
        # Count insertions
        insertions = [v for v in variants if v['type'] == 'INS']
        ins_with_seq = [v for v in insertions if v['sequence']]
        
        print(f"\nDataset Statistics (Chromosome 1):")
        print(f"  Total variants: {len(coordinates)}")
        print(f"  Insertions: {len(insertions)} ({len(ins_with_seq)} with sequences)")
        print(f"  Other SVs: {len(coordinates) - len(insertions)}")
        
        results[db_file] = {}
        


        print("WITHOUT Hamming (spatial clustering only):")


        for algo in algorithms:
            elapsed, memory, n_clusters, labels = benchmark_clustering_algorithm(
                coordinates, variants, algo, apply_hamming=False
            )
            
            if elapsed is not None:
                results[db_file][algo] = {
                    'time': elapsed,
                    'memory': memory,
                    'clusters': n_clusters
                }
                print(f"  {algo:15s}: {elapsed:7.4f}s | {n_clusters:4d} clusters | mem: {memory/1024/1024:6.2f} MB")
        

        if ins_with_seq:
            print("WITH Hamming (post-clustering for insertions):")

      
            for algo in algorithms:
                elapsed, memory, n_clusters, labels = benchmark_clustering_algorithm(
                    coordinates, variants, algo, apply_hamming=True, max_hamming=0.2
                )
                
                if elapsed is not None:
                    algo_key = f"{algo}_HAMMING"
                    results[db_file][algo_key] = {
                        'time': elapsed,
                        'memory': memory,
                        'clusters': n_clusters
                    }
                    print(f"  {algo:15s}: {elapsed:7.4f}s | {n_clusters:4d} clusters | mem: {memory/1024/1024:6.2f} MB")
    

    print("\n" + "="*60)
    print("BENCHMARK SUMMARY")
    print("="*60)

    
    for db_file in results:
        print(f"\n{os.path.basename(db_file)}:")

        for algo in algorithms:
            without = results[db_file].get(algo, None)
            with_h = results[db_file].get(f"{algo}_HAMMING", None)
            
            if without:
                print(f"\n{algo}:")
                print(f"  Without Hamming: {without['time']:.4f}s | {without['clusters']} clusters | {without['memory']/1024/1024:.2f} MB")
                
                if with_h:
                    overhead = (with_h['time'] - without['time']) / without['time'] * 100
                    cluster_diff = with_h['clusters'] - without['clusters']
                    print(f"  With Hamming:    {with_h['time']:.4f}s | {with_h['clusters']} clusters | {with_h['memory']/1024/1024:.2f} MB")
                    print(f"  Hamming Overhead: +{overhead:.1f}% time | {cluster_diff:+d} clusters")
    
    # Plotting
    try:
        import matplotlib.pyplot as plt
        
        fig, axes = plt.subplots(2, 2, figsize=(14, 10))
        fig.suptitle('SVDB Clustering Algorithm Benchmark (Chromosome 1)', fontsize=16, fontweight='bold')
        
        for db_idx, db_file in enumerate(list(results.keys())[:1]):  # First file only for clarity
            db_results = results[db_file]
            
            # Plot 1: Time comparison of the algorithms 
            ax1 = axes[0, 0]
            algos = []
            times_no_ham = []
            times_ham = []
            
            for algo in algorithms:
                if algo in db_results:
                    algos.append(algo)
                    times_no_ham.append(db_results[algo]['time'])
                    ham_key = f"{algo}_HAMMING"
                    times_ham.append(db_results.get(ham_key, {}).get('time', 0))
            
            x = np.arange(len(algos))
            width = 0.35
            ax1.bar(x - width/2, times_no_ham, width, label='Without Hamming', color='#FF6B6B')
            ax1.bar(x + width/2, times_ham, width, label='With Hamming', color='#4ECDC4')
            ax1.set_ylabel('Time (seconds)')
            ax1.set_title('Clustering Time Comparison')
            ax1.set_xticks(x)
            ax1.set_xticklabels(algos, rotation=45, ha='right')
            ax1.legend()

            
            # Plot 2: Memory comparison
            ax2 = axes[0, 1]
            mem_no_ham = [db_results[algo]['memory']/1024/1024 for algo in algos if algo in db_results]
            mem_ham = [db_results.get(f"{algo}_HAMMING", {}).get('memory', 0)/1024/1024 for algo in algos]
            
            ax2.bar(x - width/2, mem_no_ham, width, label='Without Hamming', color="#6DDAC6")
            ax2.bar(x + width/2, mem_ham, width, label='With Hamming', color="#A35E92")
            ax2.set_ylabel('Memory (MB)')
            ax2.set_title('Memory Usage Comparison')
            ax2.set_xticks(x)
            ax2.set_xticklabels(algos, rotation=45, ha='right')
            ax2.legend()

            
            # Plot 3: Clusters comparison
            ax3 = axes[1, 0]
            clusters_no_ham = [db_results[algo]['clusters'] for algo in algos if algo in db_results]
            clusters_ham = [db_results.get(f"{algo}_HAMMING", {}).get('clusters', 0) for algo in algos]
            
            ax3.bar(x - width/2, clusters_no_ham, width, label='Without Hamming', color='#AA96DA')
            ax3.bar(x + width/2, clusters_ham, width, label='With Hamming', color='#FCBAD3')
            ax3.set_ylabel('Number of Clusters')
            ax3.set_title('Cluster Count Comparison')
            ax3.set_xticks(x)
            ax3.set_xticklabels(algos, rotation=45, ha='right')
            ax3.legend()

            
            # Plot 4: Hamming overhead percentage
            ax4 = axes[1, 1]
            overhead = []
            for algo in algos:
                if algo in db_results and f"{algo}_HAMMING" in db_results:
                    base_time = db_results[algo]['time']
                    ham_time = db_results[f"{algo}_HAMMING"]['time']
                    overhead.append((ham_time - base_time) / base_time * 100 if base_time > 0 else 0)
                else:
                    overhead.append(0)
            
            ax4.bar(x, overhead, color='#A8D8EA')
            ax4.set_ylabel('Overhead (%)')
            ax4.set_title('Hamming Distance Overhead')
            ax4.set_xticks(x)
            ax4.set_xticklabels(algos, rotation=45, ha='right')
            ax4.axhline(y=0, color='black', linestyle='-', linewidth=0.5)
        
        plt.tight_layout()
        plt.savefig('algorithm_benchmark_results.png', dpi=150, bbox_inches='tight')
        print(f"\n✓ Benchmark plot saved to: algorithm_benchmark_results.png")

        #plot 5: cluster plots for each of the clustering algorithms with and without hamming
        print(f"\nGenerating cluster visualization plots...")
        for db_file in results:
            coordinates, variants = load_db_samples(db_file)
            if coordinates is None or len(coordinates) == 0:
                continue
            for algo in algorithms:
                for apply_hamming in [False, True]:
                    elapsed, memory, n_clusters, labels = benchmark_clustering_algorithm(
                        coordinates, variants, algo, apply_hamming=apply_hamming, max_hamming=0.2
                    )
                    if labels is None:
                        continue
                    
                    plt.figure(figsize=(8, 6))
                    unique_labels = set(labels)
                    colors = plt.cm.get_cmap('tab20', len(unique_labels))
                    
                    for k in unique_labels:
                        class_member_mask = (labels == k)
                        xy = coordinates[class_member_mask]
                        plt.scatter(xy[:, 0], xy[:, 1], s=10, color=colors(k), label=f'Cluster {k}' if k != -1 else 'Noise')
                    
                   #plt.title(f'{algo} Clustering {"with Hamming" if apply_hamming else "without Hamming"}\n{os.path.basename(db_file)}')
                    plt.xlabel('Position A')
                    plt.ylabel('Position B')
                    plot_filename = f"{os.path.basename(db_file).replace('.db','')}_chr1_{algo}_{'HAMMING' if apply_hamming else 'NO_HAMMING'}.png"
                    plt.savefig(plot_filename, dpi=150, bbox_inches='tight')
                    plt.close()

                    print(f"  ✓ Cluster plot saved to: {plot_filename}")
    except ImportError as e:
        print(f"\nMatplotlib not available; skipping visualization ({e})")
    except Exception as e:
        print(f"\nError creating plots: {e}")


if __name__ == "__main__":
    main()