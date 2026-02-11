#!/usr/bin/env python3


import time
import sys
import os
import glob
import numpy as np
import psutil
import csv
from datetime import datetime

def calculate_clustering_ratio(labels):
    """Calculate ratio of clustered/(clustered+unclustered)"""
    if labels is None or len(labels) == 0:
        return 0.0, 0, 0
    
    labels = np.array(labels)
    n_clustered = np.sum(labels != -1)
    n_unclustered = np.sum(labels == -1)
    total = len(labels)
    
    ratio = n_clustered / total if total > 0 else 0.0
    
    return ratio, n_clustered, n_unclustered


def hamming_distance(seq1, seq2):
    """Calculate normalized Hamming distance between sequences"""
    if not seq1 or not seq2:
        return 1.0
    seq1, seq2 = str(seq1).upper(), str(seq2).upper()
    max_len = max(len(seq1), len(seq2))
    min_len = min(len(seq1), len(seq2))
    mismatches = sum(1 for i in range(min_len) if seq1[i] != seq2[i])
    return (mismatches + (max_len - min_len)) / max_len if max_len > 0 else 0.0


def apply_hamming_reclustering(labels, variants, max_hamming=0.2):
    if labels is None:
        return labels
    
    labels = np.array(labels)
    new_labels = np.full_like(labels, -1)
    next_id = 0
    
    ins_with_seq = 0
    ins_no_seq = 0
    seq_clusters = 0
    
    # Process each spatial cluster
    for spatial_label in sorted(set(labels.tolist())):
        if spatial_label == -1:  # Skip noise
            continue
        
        indices = np.where(labels == spatial_label)[0]
        ins_with, others = [], []
        
        # Separate insertions with sequences from other variants
        for idx in indices:
            if variants[idx]['type'] == 'INS' and variants[idx]['sequence']:
                ins_with.append(idx)
            else:
                others.append(idx)
                if variants[idx]['type'] == 'INS':
                    ins_no_seq += 1
        
        # Non-insertion variants keep their spatial cluster
        if others:
            new_labels[others] = next_id
            next_id += 1
        
        # Insertions WITH sequences get re-clustered by sequence similarity
        if ins_with:
            ins_with_seq += len(ins_with)
            assigned = set()
            before_id = next_id
            
            for i in range(len(ins_with)):
                if i in assigned:
                    continue
                
                group = [i]
                assigned.add(i)
                seq_i = variants[ins_with[i]]['sequence']
                
                # Find all insertions with similar sequences
                for j in range(i + 1, len(ins_with)):
                    if j not in assigned:
                        seq_j = variants[ins_with[j]]['sequence']
                        if hamming_distance(seq_i, seq_j) <= max_hamming:
                            group.append(j)
                            assigned.add(j)
                
                # Assign new cluster ID to this sequence group
                for g in group:
                    new_labels[ins_with[g]] = next_id
                next_id += 1
            
            seq_clusters += (next_id - before_id)
    
    # *** CRITICAL FIX: Preserve noise labels ***
    new_labels[labels == -1] = -1
    
    # Print diagnostic information
    print(f"    Hamming re-clustering applied:")
    print(f"      • {ins_with_seq} insertions WITH sequences → {seq_clusters} sequence-based clusters")
    print(f"      • {ins_no_seq} insertions WITHOUT sequences → kept spatial clusters")
    
    return new_labels


def overlap_cluster(coordinates, variants, distance=500, overlap=0.6):
    """
    Overlap-based clustering using connected components
    """
    from svdb.overlap_module import isSameVariation, precise_overlap
    
    n = len(coordinates)
    adjacency = {i: set() for i in range(n)}
    
    # Build adjacency graph
    for i in range(n):
        chrA_i, posA_i = "1", int(coordinates[i][0])
        posB_i = int(coordinates[i][1])
        
        for j in range(i + 1, n):
            chrA_j, posA_j = "1", int(coordinates[j][0])
            posB_j = int(coordinates[j][1])
            
            # Quick distance filter
            if abs(posA_i - posA_j) > distance or abs(posB_i - posB_j) > distance:
                continue
            
            # Check for overlap
            if variants[i]['type'] == 'INS':
                sim, match = precise_overlap(posA_i, posB_i, posA_j, posB_j, distance)
            else:
                sim, match = isSameVariation(posA_i, posB_i, posA_j, posB_j, overlap, distance)
            
            if match:
                adjacency[i].add(j)
                adjacency[j].add(i)
    
    # Connected components (DFS)
    labels = np.full(n, -1)
    current_cluster = 0
    visited = set()
    
    def dfs(node):
        visited.add(node)
        labels[node] = current_cluster
        for neighbor in adjacency[node]:
            if neighbor not in visited:
                dfs(neighbor)
    
    for i in range(n):
        if i not in visited:
            dfs(i)
            current_cluster += 1
    
    return labels


def load_database(db_file):
    """Load variants from SVDB database (chromosome 1 only)"""
    try:
        from svdb import database
        db = database.DB(db_file)
        
        coordinates, variants = [], []
        
        print(f"Loading: {os.path.basename(db_file)}")
        
        for row in db.query('SELECT * FROM SVDB WHERE chrA = "1" AND chrB = "1"'):
            var_type, posA, posB = row[0], int(row[3]), int(row[6])
            seq = row[11] if len(row) > 11 else ''
            
            sv_size = len(seq) if var_type == 'INS' and seq else abs(posB - posA)
            
            coordinates.append([posA, posB])
            variants.append({
                'posA': posA, 'posB': posB, 'type': var_type,
                'sequence': seq, 'size': sv_size, 'sample': row[9]
            })
        
        print(f"  Loaded {len(variants)} variants from chromosome 1")
        
        return (np.array(coordinates), variants)
    
    except Exception as e:
        print(f"  ERROR loading database: {e}")
        return None, None


def benchmark_algorithm(coordinates, variants, algorithm, use_hamming, max_hamming=0.2):
    """
    Benchmark a single algorithm
    """
    try:
        process = psutil.Process()
        mem_before = process.memory_info().rss
        
        start_time = time.time()
        
        distance = 500  # Default distance threshold
        
        # Run clustering algorithm
        if algorithm == 'DBSCAN':
            from svdb.export_module import DBSCAN
            labels = DBSCAN.cluster(coordinates, distance, 2)
        
        elif algorithm == 'OPTICS':
            from svdb.optics_clustering import OPTICS
            labels = OPTICS(min_samples=2, max_eps=distance).fit_predict(coordinates)
        
        elif algorithm == 'INTERVAL_TREE':
            from svdb.interval_tree_overlap import interval_tree_cluster
            labels = interval_tree_cluster(coordinates, distance)
        
        elif algorithm == 'OVERLAP':
            labels = overlap_cluster(coordinates, variants, distance, overlap=0.6)
        
        else:
            print(f"  Unknown algorithm: {algorithm}")
            return None
        
        # Apply Hamming re-clustering if requested
        if use_hamming:
            labels = apply_hamming_reclustering(labels, variants, max_hamming)
        
        elapsed_time = time.time() - start_time
        mem_after = process.memory_info().rss
        memory_used = mem_after - mem_before
        
        # Calculate statistics
        labels = np.array(labels)
        n_clustered = np.sum(labels != -1)
        n_unclustered = np.sum(labels == -1)
        n_clusters = len(set(labels)) - (1 if -1 in labels else 0)
        
        return {
            'algorithm': algorithm,
            'use_hamming': use_hamming,
            'labels': labels,
            'time': elapsed_time,
            'memory': memory_used,
            'clustered': n_clustered,
            'unclustered': n_unclustered,
            'clusters': n_clusters
        }
    
    except Exception as e:
        print(f"  ERROR in {algorithm}: {e}")
        import traceback
        traceback.print_exc()
        return None


def analyze_single_database(db_file, algorithms, output_csv):
    """Analyze clustering ratios for a single database"""
    print(f"\n{'='*90}")
    print(f"PROCESSING: {os.path.basename(db_file)}")
    print(f"{'='*90}")
    
    coordinates, variants = load_database(db_file)
    if coordinates is None:
        return
    
    db_basename = os.path.basename(db_file).replace('.db', '')
    
    # Extract sample names from filename (e.g., NA12877_NA12878)
    sample_pair = db_basename
    
    results = []
    
    for algo in algorithms:
        print(f"\n{algo} Analysis:")
        
        # WITHOUT Hamming
        print(f"  WITHOUT Hamming...")
        result_no_ham = benchmark_algorithm(coordinates, variants, algo, False)
        if result_no_ham:
            ratio, clustered, unclustered = calculate_clustering_ratio(result_no_ham['labels'])
            print(f"    Clustered: {clustered:>5} | Unclustered: {unclustered:>5} | " 
                  f"Clusters: {result_no_ham['clusters']:>4} | Ratio: {ratio:.4f}")
            
            results.append({
                'Sample_Pair': sample_pair,
                'Algorithm': algo,
                'Hamming': 'No',
                'Clustered': clustered,
                'Unclustered': unclustered,
                'Total': clustered + unclustered,
                'Ratio': ratio,
                'Time_sec': result_no_ham['time'],
                'Memory_MB': result_no_ham['memory'] / 1024 / 1024,
                'N_Clusters': result_no_ham['clusters']
            })
        
        # WITH Hamming
        print(f"  WITH Hamming...")
        result_with_ham = benchmark_algorithm(coordinates, variants, algo, True)
        if result_with_ham:
            ratio, clustered, unclustered = calculate_clustering_ratio(result_with_ham['labels'])
            print(f"    Clustered: {clustered:>5} | Unclustered: {unclustered:>5} | "
                  f"Clusters: {result_with_ham['clusters']:>4} | Ratio: {ratio:.4f}")
            
            results.append({
                'Sample_Pair': sample_pair,
                'Algorithm': algo,
                'Hamming': 'Yes',
                'Clustered': clustered,
                'Unclustered': unclustered,
                'Total': clustered + unclustered,
                'Ratio': ratio,
                'Time_sec': result_with_ham['time'],
                'Memory_MB': result_with_ham['memory'] / 1024 / 1024,
                'N_Clusters': result_with_ham['clusters']
            })
        
        # Calculate improvement
        if result_no_ham and result_with_ham:
            ratio_no, _, _ = calculate_clustering_ratio(result_no_ham['labels'])
            ratio_with, _, _ = calculate_clustering_ratio(result_with_ham['labels'])
            improvement = ratio_with - ratio_no
            cluster_improvement = result_with_ham['clusters'] - result_no_ham['clusters']
            
            print(f"    Hamming Improvement:")
            print(f"      Ratio:    {improvement:+.4f} ({improvement*100:+.2f}%)")
            print(f"      Clusters: {cluster_improvement:+d} (from {result_no_ham['clusters']} to {result_with_ham['clusters']})")
    
    # Write results to CSV
    if results:
        write_results_to_csv(results, output_csv)
    
    return results


def write_results_to_csv(results, output_csv):
    """Write results to CSV file"""
    file_exists = os.path.exists(output_csv)
    
    with open(output_csv, 'a', newline='') as f:
        fieldnames = ['Sample_Pair', 'Algorithm', 'Hamming', 'Clustered', 'Unclustered', 
                      'Total', 'Ratio', 'Time_sec', 'Memory_MB', 'N_Clusters']
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        
        if not file_exists:
            writer.writeheader()
        
        for result in results:
            writer.writerow(result)


def create_summary_table(csv_file):
    """Create a formatted summary table from CSV results"""
    if not os.path.exists(csv_file):
        print(f"No results file found: {csv_file}")
        return
    
    print(f"{'CLUSTERING RATIO SUMMARY':^120}")

    
    # Read CSV
    with open(csv_file, 'r') as f:
        reader = csv.DictReader(f)
        data = list(reader)
    
    # Group by sample pair
    sample_pairs = sorted(set(row['Sample_Pair'] for row in data))
    algorithms = ['DBSCAN', 'OPTICS', 'INTERVAL_TREE', 'OVERLAP']
    
    for sample_pair in sample_pairs:
        print(f"\n{sample_pair}:")
        print(f"  {'Algorithm':<15} {'Hamming':<10} {'Clustered':>10} {'Unclustered':>12} {'Clusters':>10} {'Ratio':>10}")
        print(f"  {'-'*75}")
        
        for algo in algorithms:
            for hamming in ['No', 'Yes']:
                matches = [r for r in data 
                          if r['Sample_Pair'] == sample_pair 
                          and r['Algorithm'] == algo 
                          and r['Hamming'] == hamming]
                
                if matches:
                    r = matches[0]
                    print(f"  {algo:<15} {hamming:<10} {r['Clustered']:>10} {r['Unclustered']:>12} "
                          f"{r['N_Clusters']:>10} {float(r['Ratio']):>10.4f}")


def main():
    if len(sys.argv) < 2:
        print("Usage: python ratio_stats_plat.py <database1.db> <database2.db> ...")
        print("\nExample:")
        print("  python ratio_stats_plat.py NA12877_NA12878.db NA12877_NA12879.db")
        print("  python ratio_stats_plat.py *.db")
        sys.exit(1)
    
    # Collect database files
    db_files = []
    for arg in sys.argv[1:]:
        if '*' in arg:
            db_files.extend(glob.glob(arg))
        elif os.path.exists(arg):
            db_files.append(arg)
    
    if not db_files:
        print("ERROR: No valid database files found!")
        sys.exit(1)
    
    # Setup output CSV
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_csv = f"clustering_ratios_{timestamp}.csv"
    

    print(f"{'SVDB CLUSTERING RATIO ANALYSIS - FIXED VERSION':^90}")

    print(f"Databases to analyze: {len(db_files)}")
    for f in db_files:
        print(f"  • {os.path.basename(f)}")
    print(f"\nOutput CSV: {output_csv}")
    print(f"\nFixes applied:")
    print(f"  ✓ Hamming re-clustering now preserves noise labels")
    print(f"  ✓ OPTICS uses class-based implementation")
    print(f"  ✓ OVERLAP uses proper overlap_cluster function")
    print(f"  ✓ Diagnostic output shows sequence-based clusters")
    
    algorithms = ['DBSCAN', 'OPTICS', 'INTERVAL_TREE', 'OVERLAP']
    
    # Analyze each database
    all_results = []
    for db_file in db_files:
        results = analyze_single_database(db_file, algorithms, output_csv)
        if results:
            all_results.extend(results)
    
    # Create summary table
    create_summary_table(output_csv)
    
    print(f"{'✓ ANALYSIS COMPLETE':^90}")
    print(f"Results saved to: {output_csv}")
    print(f"Total comparisons: {len(all_results)}")


if __name__ == "__main__":
    main()

