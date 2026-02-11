#!/usr/bin/env python3
"""
Hamming vs Levenshtein Distance Comparison for SVDB Clustering
Implements the full SVDB workflow: spatial clustering + sequence-based re-clustering
"""

import sys
import os
import time
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Add SVDB to path
sys.path.insert(0, os.getcwd())
from svdb.database import DB
from svdb.export_module import DBSCAN
from svdb.optics_clustering import optics_cluster
from svdb.interval_tree_overlap import interval_tree_cluster


def hamming_distance(seq1, seq2):
    """Calculate normalized Hamming distance."""
    if not seq1 or not seq2:
        return 1.0
    seq1, seq2 = str(seq1).upper(), str(seq2).upper()
    max_len = max(len(seq1), len(seq2))
    min_len = min(len(seq1), len(seq2))
    mismatches = sum(1 for i in range(min_len) if seq1[i] != seq2[i])
    return (mismatches + (max_len - min_len)) / max_len if max_len > 0 else 0.0


def levenshtein_distance(seq1, seq2):
    """Calculate normalized Levenshtein distance."""
    if not seq1 or not seq2:
        return 1.0
    
    seq1 = str(seq1).upper()
    seq2 = str(seq2).upper()
    
    len1, len2 = len(seq1), len(seq2)
    max_len = max(len1, len2)
    
    # DP table
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


def load_variants(db_file, chromosome='1'):
    """Load ALL variants (not just INS) from database."""
    db = DB(db_file)
    variants = []
    
    query = f'SELECT var, posA, posB, sequence FROM SVDB WHERE chrA = "{chromosome}" AND chrB = "{chromosome}"'
    
    for row in db.query(query):
        var_type = row[0]
        posA, posB = int(row[1]), int(row[2])
        seq = row[3] if len(row) > 3 and row[3] else ''
        
        variants.append({
            'type': var_type,
            'posA': posA,
            'posB': posB,
            'sequence': seq
        })
    
    return variants


def apply_sequence_reclustering(labels, variants, max_threshold, distance_func):
    """
    Re-cluster INS variants within spatial clusters using sequence similarity.
    This is SVDB's actual workflow.
    """
    if labels is None:
        return labels
    
    labels = np.array(labels)
    new_labels = np.full_like(labels, -1)
    next_cluster_id = 0
    
    # Get unique spatial cluster IDs
    unique_labels = sorted(set(labels.tolist()))
    
    for spatial_label in unique_labels:
        if spatial_label == -1:  # Noise
            continue
        
        # Get all variants in this spatial cluster
        indices = np.where(labels == spatial_label)[0]
        if len(indices) == 0:
            continue
        
        # Separate INS with sequences from others
        ins_with_seq = []
        other_indices = []
        
        for idx in indices:
            var = variants[idx]
            if var['type'] == 'INS' and var['sequence'] and len(var['sequence']) >= 10:
                ins_with_seq.append(idx)
            else:
                other_indices.append(idx)
        
        # Non-INS variants keep their spatial cluster assignment
        if other_indices:
            new_labels[other_indices] = next_cluster_id
            next_cluster_id += 1
        
        # Re-cluster INS by sequence similarity
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
                    
                    # Use specified distance function
                    if distance_func == 'hamming':
                        dist = hamming_distance(seq_i, seq_j)
                    else:  # levenshtein
                        dist = levenshtein_distance(seq_i, seq_j)
                    
                    if dist <= max_threshold:
                        group.append(j)
                        assigned.add(j)
                
                # Assign new cluster ID for this sequence group
                for g in group:
                    new_labels[ins_with_seq[g]] = next_cluster_id
                next_cluster_id += 1
    
    return new_labels


def cluster_with_metric(variants, algorithm='DBSCAN', distance_func='hamming', max_threshold=0.2):
    """
    Full SVDB clustering pipeline:
    1. Spatial clustering
    2. Sequence-based re-clustering for INS
    """
    # Extract coordinates
    coordinates = np.array([[v['posA'], v['posB']] for v in variants])
    
    print(f"\n  Stage 1: Spatial clustering with {algorithm}...")
    start = time.time()
    
    # Stage 1: Spatial clustering
    if algorithm == 'DBSCAN':
        spatial_labels = DBSCAN.cluster(coordinates, epsilon=500, m=2)
    elif algorithm == 'OPTICS':
        spatial_labels = optics_cluster(coordinates, min_samples=2, max_eps=500)
    elif algorithm == 'INTERVAL_TREE':
        spatial_labels = interval_tree_cluster(coordinates, max_distance=500)
    else:
        raise ValueError(f"Unknown algorithm: {algorithm}")
    
    spatial_time = time.time() - start
    
    # Count spatial clusters
    unique_spatial = set(spatial_labels.tolist())
    n_spatial_clusters = len(unique_spatial) - (1 if -1 in unique_spatial else 0)
    n_spatial_noise = (spatial_labels == -1).sum()
    
    print(f"    Spatial clusters: {n_spatial_clusters}")
    print(f"    Spatial noise: {n_spatial_noise}")
    print(f"    Time: {spatial_time:.2f}s")
    
    # Stage 2: Sequence-based re-clustering
    print(f"  Stage 2: Re-clustering INS with {distance_func} (threshold={max_threshold})...")
    start = time.time()
    
    final_labels = apply_sequence_reclustering(
        spatial_labels, variants, max_threshold, distance_func
    )
    
    sequence_time = time.time() - start
    
    # Count final clusters
    unique_final = set(final_labels.tolist())
    n_final_clusters = len(unique_final) - (1 if -1 in unique_final else 0)
    n_final_noise = (final_labels == -1).sum()
    
    print(f"    Final clusters: {n_final_clusters}")
    print(f"    Final noise: {n_final_noise}")
    print(f"    Time: {sequence_time:.2f}s")
    
    return {
        'spatial_labels': spatial_labels,
        'final_labels': final_labels,
        'n_spatial_clusters': n_spatial_clusters,
        'n_final_clusters': n_final_clusters,
        'n_spatial_noise': n_spatial_noise,
        'n_final_noise': n_final_noise,
        'spatial_time': spatial_time,
        'sequence_time': sequence_time,
        'total_time': spatial_time + sequence_time
    }


def compare_clustering_results(hamming_results, lev_results, variants):
    """Compare clustering outcomes between Hamming and Levenshtein."""
    h_labels = hamming_results['final_labels']
    l_labels = lev_results['final_labels']
    
    print(f"\n{'-'*80}")
    print("CLUSTERING COMPARISON")
    print(f"{'-'*80}")
    
    # Basic metrics
    print(f"{'Metric':<40} {'Hamming':>15} {'Levenshtein':>15}")
    print(f"{'-'*80}")
    print(f"{'Spatial clusters (same for both)':<40} {hamming_results['n_spatial_clusters']:>15} "
          f"{lev_results['n_spatial_clusters']:>15}")
    print(f"{'Final clusters':<40} {hamming_results['n_final_clusters']:>15} "
          f"{lev_results['n_final_clusters']:>15}")
    print(f"{'Noise points':<40} {hamming_results['n_final_noise']:>15} "
          f"{lev_results['n_final_noise']:>15}")
    print(f"{'-'*80}")
    
    # Agreement analysis
    agreement = (h_labels == l_labels).sum()
    total = len(h_labels)
    agreement_pct = 100 * agreement / total
    
    print(f"\nCluster assignment agreement: {agreement}/{total} ({agreement_pct:.2f}%)")
    
    # Disagreement analysis
    disagreement_mask = (h_labels != l_labels)
    disagreements = disagreement_mask.sum()
    
    if disagreements > 0:
        print(f"Disagreements: {disagreements} ({100*disagreements/total:.2f}%)")
        
        # Analyze disagreement types
        ins_disagreements = sum(1 for i in range(len(variants)) 
                               if disagreement_mask[i] and variants[i]['type'] == 'INS')
        other_disagreements = disagreements - ins_disagreements
        
        print(f"  INS disagreements: {ins_disagreements}")
        print(f"  Other disagreements: {other_disagreements}")
    
    # INS-specific analysis
    ins_indices = [i for i, v in enumerate(variants) if v['type'] == 'INS' and v['sequence']]
    if ins_indices:
        ins_agreement = sum(1 for i in ins_indices if h_labels[i] == l_labels[i])
        ins_agreement_pct = 100 * ins_agreement / len(ins_indices)
        
        print(f"\nINS with sequences: {len(ins_indices)}")
        print(f"INS agreement: {ins_agreement}/{len(ins_indices)} ({ins_agreement_pct:.2f}%)")
        
        # Cluster size differences for INS
        h_ins_clusters = len(set(h_labels[ins_indices])) - (1 if -1 in h_labels[ins_indices] else 0)
        l_ins_clusters = len(set(l_labels[ins_indices])) - (1 if -1 in l_labels[ins_indices] else 0)
        
        print(f"INS clusters (Hamming): {h_ins_clusters}")
        print(f"INS clusters (Levenshtein): {l_ins_clusters}")
    
    # Timing
    print(f"\n{'-'*80}")
    print("PERFORMANCE")
    print(f"{'-'*80}")
    print(f"{'Stage':<40} {'Hamming':>15} {'Levenshtein':>15}")
    print(f"{'-'*80}")
    print(f"{'Spatial clustering (s)':<40} {hamming_results['spatial_time']:>15.2f} "
          f"{lev_results['spatial_time']:>15.2f}")
    print(f"{'Sequence re-clustering (s)':<40} {hamming_results['sequence_time']:>15.2f} "
          f"{lev_results['sequence_time']:>15.2f}")
    print(f"{'Total time (s)':<40} {hamming_results['total_time']:>15.2f} "
          f"{lev_results['total_time']:>15.2f}")
    print(f"{'-'*80}")
    
    return {
        'agreement': agreement,
        'agreement_pct': agreement_pct,
        'disagreements': disagreements,
        'ins_agreement_pct': ins_agreement_pct if ins_indices else 0
    }


def plot_comparison(hamming_results, lev_results, variants, output_prefix):
    """Generate comparison plots."""
    h_labels = hamming_results['final_labels']
    l_labels = lev_results['final_labels']
    
    # Plot 1: Cluster count comparison
    fig, ax = plt.subplots(figsize=(10, 6))
    
    categories = ['Spatial\nClusters', 'Final\nClusters', 'Noise\nPoints']
    hamming_vals = [
        hamming_results['n_spatial_clusters'],
        hamming_results['n_final_clusters'],
        hamming_results['n_final_noise']
    ]
    lev_vals = [
        lev_results['n_spatial_clusters'],
        lev_results['n_final_clusters'],
        lev_results['n_final_noise']
    ]
    
    x = np.arange(len(categories))
    width = 0.35
    
    ax.bar(x - width/2, hamming_vals, width, label='Hamming', 
           color='#3498DB', alpha=0.9, edgecolor='black', linewidth=1.5)
    ax.bar(x + width/2, lev_vals, width, label='Levenshtein',
           color='#E74C3C', alpha=0.9, edgecolor='black', linewidth=1.5)
    
    ax.set_xlabel('Metric', fontsize=12, fontweight='bold')
    ax.set_ylabel('Count', fontsize=12, fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels(categories, fontsize=11)
    ax.legend(fontsize=11, frameon=True)
    ax.grid(True, alpha=0.3, axis='y')
    
    plt.tight_layout()
    filename = f"{output_prefix}_cluster_comparison.png"
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  ✓ Saved: {filename}")
    
    # Plot 2: Agreement visualization
    fig, ax = plt.subplots(figsize=(10, 6))
    
    agreement = (h_labels == l_labels).sum()
    disagreement = (h_labels != l_labels).sum()
    
    colors = ['#2ECC71', '#E74C3C']
    labels_pie = ['Agreement', 'Disagreement']
    values = [agreement, disagreement]
    
    ax.pie(values, labels=labels_pie, autopct='%1.1f%%', startangle=90,
           colors=colors, textprops={'fontsize': 12, 'fontweight': 'bold'})
    
    plt.tight_layout()
    filename = f"{output_prefix}_agreement.png"
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  ✓ Saved: {filename}")
    
    # Plot 3: Timing comparison
    fig, ax = plt.subplots(figsize=(10, 6))
    
    stages = ['Spatial\nClustering', 'Sequence\nRe-clustering', 'Total\nTime']
    hamming_times = [
        hamming_results['spatial_time'],
        hamming_results['sequence_time'],
        hamming_results['total_time']
    ]
    lev_times = [
        lev_results['spatial_time'],
        lev_results['sequence_time'],
        lev_results['total_time']
    ]
    
    x = np.arange(len(stages))
    width = 0.35
    
    ax.bar(x - width/2, hamming_times, width, label='Hamming',
           color='#3498DB', alpha=0.9, edgecolor='black', linewidth=1.5)
    ax.bar(x + width/2, lev_times, width, label='Levenshtein',
           color='#E74C3C', alpha=0.9, edgecolor='black', linewidth=1.5)
    
    ax.set_xlabel('Stage', fontsize=12, fontweight='bold')
    ax.set_ylabel('Time (seconds)', fontsize=12, fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels(stages, fontsize=11)
    ax.legend(fontsize=11, frameon=True)
    ax.grid(True, alpha=0.3, axis='y')
    
    plt.tight_layout()
    filename = f"{output_prefix}_timing.png"
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  ✓ Saved: {filename}")


def main():
    import glob
    
    # Find databases
    db_pattern = "./Platinum_experiments/NA12*_NA12*.db"
    db_files = sorted(glob.glob(db_pattern))
    
    if not db_files:
        print(f"ERROR: No databases found matching {db_pattern}")
        sys.exit(1)
    
    print(f"\n{'='*80}")
    print(f"HAMMING vs LEVENSHTEIN CLUSTERING COMPARISON")
    print(f"Full SVDB workflow: Spatial clustering + Sequence re-clustering")
    print(f"{'='*80}")
    print(f"Found {len(db_files)} pairwise databases")
    
    # Configuration
    algorithm = 'DBSCAN'  # or 'OPTICS' or 'INTERVAL_TREE'
    threshold = 0.2
    
    print(f"\nConfiguration:")
    print(f"  Spatial clustering: {algorithm}")
    print(f"  Sequence threshold: {threshold}")
    
    # Process each database
    all_results = []
    
    for db_idx, db_file in enumerate(db_files, 1):
        db_name = os.path.basename(db_file).replace('.db', '')
        print(f"\n{'='*80}")
        print(f"DATABASE {db_idx}/{len(db_files)}: {db_name}")
        print(f"{'='*80}")
        
        # Load variants
        print(f"Loading variants from chromosome 1...")
        variants = load_variants(db_file, chromosome='1')
        print(f"  Total variants: {len(variants)}")
        
        # Count by type
        type_counts = {}
        for v in variants:
            type_counts[v['type']] = type_counts.get(v['type'], 0) + 1
        
        print(f"  Variant types:")
        for vtype, count in sorted(type_counts.items()):
            print(f"    {vtype}: {count}")
        
        ins_with_seq = sum(1 for v in variants if v['type'] == 'INS' and v['sequence'])
        print(f"  INS with sequences: {ins_with_seq}")
        
        if len(variants) < 2:
            print(f"Skipping {db_name} - insufficient variants")
            continue
        
        # Cluster with Hamming
        print(f"\n{'─'*80}")
        print("HAMMING DISTANCE CLUSTERING")
        print(f"{'─'*80}")
        hamming_results = cluster_with_metric(
            variants, algorithm=algorithm, 
            distance_func='hamming', max_threshold=threshold
        )
        
        # Cluster with Levenshtein
        print(f"\n{'─'*80}")
        print("LEVENSHTEIN DISTANCE CLUSTERING")
        print(f"{'─'*80}")
        lev_results = cluster_with_metric(
            variants, algorithm=algorithm,
            distance_func='levenshtein', max_threshold=threshold
        )
        
        # Compare results
        comparison = compare_clustering_results(hamming_results, lev_results, variants)
        
        # Generate plots
        print(f"\nGenerating plots...")
        output_prefix = f"clustering_comparison_{db_name}"
        plot_comparison(hamming_results, lev_results, variants, output_prefix)
        
        all_results.append({
            'db_name': db_name,
            'n_variants': len(variants),
            'n_ins': ins_with_seq,
            'hamming': hamming_results,
            'levenshtein': lev_results,
            'comparison': comparison
        })
    
    # Summary across all databases
    if all_results:
        print(f"\n{'='*80}")
        print(f"SUMMARY ACROSS ALL DATABASES")
        print(f"{'='*80}")
        
        total_variants = sum(r['n_variants'] for r in all_results)
        total_ins = sum(r['n_ins'] for r in all_results)
        avg_agreement = np.mean([r['comparison']['agreement_pct'] for r in all_results])
        avg_ins_agreement = np.mean([r['comparison']['ins_agreement_pct'] for r in all_results if r['n_ins'] > 0])
        
        print(f"Total databases analyzed: {len(all_results)}")
        print(f"Total variants: {total_variants}")
        print(f"Total INS with sequences: {total_ins}")
        print(f"Average overall agreement: {avg_agreement:.2f}%")
        print(f"Average INS agreement: {avg_ins_agreement:.2f}%")
        
        print(f"\n{'-'*80}")
        print(f"{'Database':<30} {'Variants':>12} {'INS':>12} {'Agreement %':>15}")
        print(f"{'-'*80}")
        for r in all_results:
            print(f"{r['db_name']:<30} {r['n_variants']:>12} {r['n_ins']:>12} "
                  f"{r['comparison']['agreement_pct']:>15.2f}")
    
    print(f"\n{'='*80}")
    print(f"ANALYSIS COMPLETE")
    print(f"{'='*80}\n")


if __name__ == "__main__":
    main()
