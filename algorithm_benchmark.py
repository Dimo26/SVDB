#!/usr/bin/env python3
"""
SVDB Clustering Algorithm Benchmark 
"""

import time
import sys
import os
import glob
import numpy as np
import psutil

def print_sv_statistics(variants):
    if not variants:
        print("No variants found")
        return None
    
    stats = {}
    for v in variants:
        sv_type = v['type']
        if sv_type not in stats:
            stats[sv_type] = {'count': 0, 'sizes': [], 'with_seq': 0, 'without_seq': 0}
        stats[sv_type]['count'] += 1
        stats[sv_type]['sizes'].append(v['size'])
        if sv_type == 'INS':
            if v['sequence'] and len(v['sequence']) > 0:
                stats[sv_type]['with_seq'] += 1
            else:
                stats[sv_type]['without_seq'] += 1
    

    print(f"{'STRUCTURAL VARIANT STATISTICS - Chromosome 1':^90}")
    print(f"{'SV TYPE':<10} {'COUNT':>8} {'%':>7} {'AVG_SIZE':>12} {'MIN':>10} {'MAX':>12} {'WITH_SEQ':>10} {'NO_SEQ':>8}")

    
    total = len(variants)
    for sv_type in sorted(stats.keys()):
        data = stats[sv_type]
        pct = data['count'] / total * 100 if total > 0 else 0
        avg_size = np.mean(data['sizes']) if data['sizes'] else 0
        min_size = min(data['sizes']) if data['sizes'] else 0
        max_size = max(data['sizes']) if data['sizes'] else 0
        
        print(f"{sv_type:<10} {data['count']:>8} {pct:>6.1f}% {avg_size:>12.1f} "
              f"{min_size:>10} {max_size:>12} {data['with_seq']:>10} {data['without_seq']:>8}")
    
    print(f"{'TOTAL':<10} {total:>8} {'100.0%':>7}")
    
    if 'INS' in stats:
        ins = stats['INS']
        print(f"{'INSERTION DETAILS':^90}")
        print(f"  Total insertions: {ins['count']}")
        if ins['count'] > 0:
            print(f"  With sequence:    {ins['with_seq']:>5} ({ins['with_seq']/ins['count']*100:>5.1f}%)  ← Can be Hamming-clustered")
            print(f"  Without sequence: {ins['without_seq']:>5} ({ins['without_seq']/ins['count']*100:>5.1f}%)  ← Spatial only")
            
            seq_sizes = [s for s in ins['sizes'] if s > 0]
            if seq_sizes:
                print(f"\n  Sequence lengths:")
                print(f"    Mean:   {np.mean(seq_sizes):>8.1f} bp")
                print(f"    Median: {np.median(seq_sizes):>8.1f} bp")
                print(f"    Std:    {np.std(seq_sizes):>8.1f} bp")
                print(f"    Range:  {min(seq_sizes)}-{max(seq_sizes)} bp")
    
    return stats

def analyze_clustering(coordinates, variants, labels, algorithm_name):
    if labels is None:
        return
    
    n_clustered = np.sum(labels != -1)
    n_noise = np.sum(labels == -1)
    n_clusters = len(set(labels)) - (1 if -1 in labels else 0)

    print(f"{algorithm_name} - CLUSTERING ANALYSIS")
    print(f"  Total variants:    {len(variants)}")
    print(f"  Clustered:         {n_clustered:>5} ({n_clustered/len(variants)*100:>5.1f}%)")
    print(f"  Unclustered/Noise: {n_noise:>5} ({n_noise/len(variants)*100:>5.1f}%) ← Why so many?")
    print(f"  Number of clusters: {n_clusters}")
    
    # Breakdown by SV type
    print(f"\n  Noise breakdown by SV type:")
    sv_stats = {}
    for v, label in zip(variants, labels):
        sv_type = v['type']
        if sv_type not in sv_stats:
            sv_stats[sv_type] = {'total': 0, 'noise': 0}
        sv_stats[sv_type]['total'] += 1
        if label == -1:
            sv_stats[sv_type]['noise'] += 1
    
    for sv_type in sorted(sv_stats.keys()):
        data = sv_stats[sv_type]
        pct = data['noise']/data['total']*100 if data['total'] > 0 else 0
        print(f"    {sv_type:<6} {data['noise']:>5}/{data['total']:<5} ({pct:>5.1f}%) are noise")
    
    # Distance analysis - WHY are they noise?
    noise_indices = np.where(labels == -1)[0]
    if len(noise_indices) > 0:
        print(f"  NOISE DIAGNOSIS - Why aren't these variants clustered?")
 
        sample_size = min(100, len(noise_indices))
        sample_indices = np.random.choice(noise_indices, sample_size, replace=False) if len(noise_indices) > 100 else noise_indices
        
        try:
            from scipy.spatial.distance import cdist
            noise_coords = coordinates[sample_indices]
            all_coords = coordinates

            min_distances = []
            for nc in noise_coords:
                distances = np.sqrt(np.sum((all_coords - nc)**2, axis=1))
                distances_nonzero = distances[distances > 0]
                if len(distances_nonzero) > 0:
                    min_distances.append(np.min(distances_nonzero))
            
            if min_distances:
                print(f"\n  Distance to nearest neighbor (sampled {len(min_distances)} noise points):")
                print(f"    Mean:   {np.mean(min_distances):>12,.1f} bp")
                print(f"    Median: {np.median(min_distances):>12,.1f} bp")
                print(f"    Min:    {np.min(min_distances):>12,.1f} bp")
                print(f"    Max:    {np.max(min_distances):>12,.1f} bp")
                print(f"\n  Clustering threshold: 500 bp (current setting)")
                
                within_500 = sum(1 for d in min_distances if d <= 500)
                within_1000 = sum(1 for d in min_distances if d <= 1000)
                within_5000 = sum(1 for d in min_distances if d <= 5000)
                
                print(f"\n  How many noise points have neighbors within:")
                print(f"    ≤500 bp:   {within_500:>5}/{len(min_distances):<5} ({within_500/len(min_distances)*100:>5.1f}%)")
                print(f"    ≤1000 bp:  {within_1000:>5}/{len(min_distances):<5} ({within_1000/len(min_distances)*100:>5.1f}%)")
                print(f"    ≤5000 bp:  {within_5000:>5}/{len(min_distances):<5} ({within_5000/len(min_distances)*100:>5.1f}%)")
                print(f"  DIAGNOSIS:")
                
                if within_500 < len(min_distances) * 0.1:
                    print(f"    ✓ Most noise points are >500bp from neighbors")
                    print(f"    → This is EXPECTED - they are genuinely isolated")
                    print(f"    → To reduce noise: increase distance threshold")
                    print(f"       Example: Try 1000bp or 5000bp threshold")
                elif within_500 > len(min_distances) * 0.3:
                    print(f"    ⚠ WARNING: {within_500/len(min_distances)*100:.0f}% of noise points ARE close (<500bp)!")
                    print(f"    → Algorithm may have issues (min_samples too high?)")
                    print(f"    → DBSCAN requires min_samples=2 neighbors")
                    print(f"    → These points might be singletons or have only 1 neighbor")
                else:
                    print(f"    • Moderate isolation - mix of close and distant points")
                    print(f"    → Some genuinely isolated, some could cluster with larger threshold")
        except Exception as e:
            print(f"  (Distance analysis failed: {e})")
    
def load_database(db_file):
    """Load variants from SVDB database."""
    try:
        from svdb import database
        db = database.DB(db_file)
        
        coordinates, variants = [], []
        
        print(f"\nLoading: {os.path.basename(db_file)}")
        
        for row in db.query('SELECT * FROM SVDB WHERE chrA = "1" AND chrB = "1"'):
            var_type, posA, posB = row[0], int(row[3]), int(row[6])
            seq = row[11] if len(row) > 11 else ''
            
            sv_size = len(seq) if var_type == 'INS' and seq else abs(posB - posA)
            
            coordinates.append([posA, posB])
            variants.append({
                'posA': posA, 'posB': posB, 'type': var_type,
                'sequence': seq, 'size': sv_size, 'sample': row[9]
            })
        
        print(f"✓ Loaded {len(variants)} variants")
        
        stats = print_sv_statistics(variants)
        
        return (np.array(coordinates), variants, stats)
    
    except Exception as e:
        print(f"Error loading database: {e}")
        import traceback
        traceback.print_exc()
        return None, None, None


def hamming_distance(seq1, seq2):
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
    
    ins_with_seq, ins_no_seq = 0, 0
    seq_clusters = 0
    
    for spatial_label in sorted(set(labels.tolist())):
        if spatial_label == -1:
            continue
        
        indices = np.where(labels == spatial_label)[0]
        ins_with, others = [], []
        
        for idx in indices:
            if variants[idx]['type'] == 'INS' and variants[idx]['sequence']:
                ins_with.append(idx)
            else:
                others.append(idx)
                if variants[idx]['type'] == 'INS':
                    ins_no_seq += 1
        
        if others:
            new_labels[others] = next_id
            next_id += 1
        
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
                
                for j in range(i+1, len(ins_with)):
                    if j not in assigned:
                        if hamming_distance(seq_i, variants[ins_with[j]]['sequence']) <= max_hamming:
                            group.append(j)
                            assigned.add(j)
                
                for g in group:
                    new_labels[ins_with[g]] = next_id
                next_id += 1
            
            seq_clusters += (next_id - before_id)
    
    new_labels[labels == -1] = -1
    
    print(f"  Hamming re-clustering:")
    print(f"    • {ins_with_seq} insertions WITH sequences → {seq_clusters} clusters")
    print(f"    • {ins_no_seq} insertions WITHOUT sequences → kept spatial clusters")
    
    return new_labels


def benchmark_algorithm(coordinates, variants, algorithm_name, apply_hamming=False):
    """Run clustering algorithm and measure performance."""
    from svdb.export_module import DBSCAN
    from svdb.optics_clustering import OPTICS
    from svdb.interval_tree_overlap import interval_tree_cluster
    
    if coordinates is None or len(coordinates) == 0:
        return None
    
    proc = psutil.Process()
    mem_before = proc.memory_info().rss
    start = time.time()
    
    try:
        distance = 500
        if algorithm_name == 'DBSCAN':
            labels = DBSCAN.cluster(coordinates, distance, 2)
        elif algorithm_name == 'OPTICS':
            labels = OPTICS(min_samples=2, max_eps=distance).fit_predict(coordinates)
        elif algorithm_name == 'INTERVAL_TREE':
            labels = interval_tree_cluster(coordinates, distance)
        else:
            return None
        
        if apply_hamming:
            labels = apply_hamming_reclustering(labels, variants)
        
        elapsed = time.time() - start
        memory = proc.memory_info().rss - mem_before
        
        return {
            'time': elapsed,
            'memory': memory,
            'labels': labels,
            'clusters': len(set(labels)) - (1 if -1 in labels else 0),
            'clustered': np.sum(labels != -1),
            'noise': np.sum(labels == -1)
        }
    except Exception as e:
        print(f"Error in {algorithm_name}: {e}")
        return None


def create_single_plot(coordinates, variants, labels, algorithm_name, hamming_state, db_basename, runtime):
    """Create single plot showing ALL SV types together with legend at bottom."""
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch
    
    colors = {'INS': '#FFD700', 'DEL': '#32CD32', 'DUP': '#4169E1',
              'INV': '#FF4500', 'BND': '#8B008B', 'NOISE': '#808080'}
    
    fig, ax = plt.subplots(figsize=(16, 12))
    
    # Plot all points
    for coord, variant, label in zip(coordinates, variants, labels):
        if label == -1:
            color, alpha, size, marker = colors['NOISE'], 0.4, 25, 'x'
        else:
            color = colors.get(variant['type'], '#000000')
            alpha, size, marker = 0.75, 55, 'o'
        
        ax.scatter(coord[0], coord[1], c=[color], s=size, alpha=alpha,
                  edgecolors='black' if label != -1 else 'none',
                  linewidths=0.7, marker=marker, zorder=2)
    
    # Legend at bottom
    legend_elements = []
    for sv_type in sorted(set(v['type'] for v in variants)):
        count = sum(1 for v in variants if v['type'] == sv_type)
        legend_elements.append(
            Patch(facecolor=colors.get(sv_type, '#000000'),
                  edgecolor='black', label=f'{sv_type} (n={count})')
        )
    
    noise_count = sum(1 for l in labels if l == -1)
    legend_elements.append(
        Patch(facecolor=colors['NOISE'],
              edgecolor='black', label=f'Noise (n={noise_count})')
    )
    
    ax.legend(handles=legend_elements, loc='upper center',
             bbox_to_anchor=(0.5, -0.05), ncol=len(legend_elements),
             fontsize=11, frameon=True)
    
    ax.set_xlabel('Position A (bp)', fontsize=14, fontweight='bold')
    ax.set_ylabel('Position B (bp)', fontsize=14, fontweight='bold')

    
    n_clusters = len(set(labels)) - (1 if -1 in labels else 0)
    n_clustered = sum(1 for l in labels if l != -1)
    
    ax.set_title(f'{algorithm_name} | {hamming_state}\n'
                 f'{n_clusters} clusters | {n_clustered} clustered | {noise_count} noise | {runtime:.3f}s',
                 fontsize=14, fontweight='bold', pad=15)
    
    plt.tight_layout()
    filename = f"{db_basename}_chr1_{algorithm_name}_{hamming_state}.png"
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.close()
    
    return filename


def create_bar_plots(results, algorithms, db_basename):
    """Create bar plots for time and memory comparison."""
    import matplotlib.pyplot as plt
    
    times_no = [results[a]['time'] for a in algorithms if a in results]
    times_yes = [results[f"{a}_HAMMING"]['time'] for a in algorithms if f"{a}_HAMMING" in results]
    mem_no = [results[a]['memory']/1024/1024 for a in algorithms if a in results]
    mem_yes = [results[f"{a}_HAMMING"]['memory']/1024/1024 for a in algorithms if f"{a}_HAMMING" in results]
    
    x = np.arange(len(algorithms))
    width = 0.35
    
    # Time plot
    fig, ax = plt.subplots(figsize=(12, 7))
    bars1 = ax.bar(x - width/2, times_no, width, label='WITHOUT Hamming',
                   color='#FF6B6B', edgecolor='black', linewidth=1.5)
    bars2 = ax.bar(x + width/2, times_yes, width, label='WITH Hamming',
                   color='#4ECDC4', edgecolor='black', linewidth=1.5)
    
    for bars in [bars1, bars2]:
        for bar in bars:
            h = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2., h, f'{h:.3f}s',
                   ha='center', va='bottom', fontsize=10, fontweight='bold')
    
    ax.set_ylabel('Time (seconds)', fontsize=13, fontweight='bold')
    ax.set_title(f'Runtime Comparison - {db_basename}', fontsize=14, fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels(algorithms, fontsize=12)
    ax.legend(loc='upper left', fontsize=11)

    plt.tight_layout()
    time_file = f"{db_basename}_time_comparison.png"
    plt.savefig(time_file, dpi=300, bbox_inches='tight')
    plt.close()

    fig, ax = plt.subplots(figsize=(12, 7))
    bars1 = ax.bar(x - width/2, mem_no, width, color='#AA96DA',
                   edgecolor='black', linewidth=1.5)
    bars2 = ax.bar(x + width/2, mem_yes, width, color='#FCBAD3',
                   edgecolor='black', linewidth=1.5)
    
    for bars in [bars1, bars2]:
        for bar in bars:
            h = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2., h, f'{h:.1f}MB',
                   ha='center', va='bottom', fontsize=10, fontweight='bold')
    
    ax.set_ylabel('Memory (MB)', fontsize=13, fontweight='bold')
    ax.set_title(f'Memory Comparison - {db_basename}', fontsize=14, fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels(algorithms, fontsize=12)
    
    plt.tight_layout()
    mem_file = f"{db_basename}_memory_comparison.png"
    plt.savefig(mem_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    return time_file, mem_file


def main():
    if len(sys.argv) < 2:
        print("Usage: python algorithm_benchmark_final.py <database.db> [...]")
        sys.exit(1)
    
    # Get database files
    db_files = []
    for arg in sys.argv[1:]:
        if '*' in arg:
            db_files.extend(glob.glob(arg))
        elif os.path.exists(arg):
            db_files.append(arg)
    
    if not db_files:
        print("Error: No valid database files!")
        sys.exit(1)
    
    print(f"{'SVDB CLUSTERING BENCHMARK - FINAL VERSION':^90}")
    print(f"Databases: {len(db_files)}")
    for f in db_files:
        print(f"  • {os.path.basename(f)}")

    
    algorithms = ['DBSCAN', 'OPTICS', 'INTERVAL_TREE']
    
    for db_file in db_files:
        print(f"{'PROCESSING: ' + os.path.basename(db_file):^90}")

        coordinates, variants, stats = load_database(db_file)
        if coordinates is None:
            continue
        
        db_results = {}
        db_basename = os.path.basename(db_file).replace('.db', '')
        
        # Run benchmarks
        for algo in algorithms:
            print(f"\n{'─'*90}")
            print(f"{algo} CLUSTERING")
            print(f"{'─'*90}")
            
            # Without Hamming
            print(f"\n▶ {algo} WITHOUT Hamming...")
            result = benchmark_algorithm(coordinates, variants, algo, False)
            if result:
                db_results[algo] = result
                print(f"  Time: {result['time']:.4f}s | Clusters: {result['clusters']} | "
                      f"Clustered: {result['clustered']} | Noise: {result['noise']}")
                analyze_clustering(coordinates, variants, result['labels'], f"{algo} WITHOUT Hamming")
            
            # With Hamming
            print(f"\n▶ {algo} WITH Hamming...")
            result = benchmark_algorithm(coordinates, variants, algo, True)
            if result:
                db_results[f"{algo}_HAMMING"] = result
                print(f"  Time: {result['time']:.4f}s | Clusters: {result['clusters']} | "
                      f"Clustered: {result['clustered']} | Noise: {result['noise']}")
                analyze_clustering(coordinates, variants, result['labels'], f"{algo} WITH Hamming")
        
        # Create plots
        print("GENERATING VISUALIZATIONS")

        try:
            # Main plots
            for algo in algorithms:
                for suffix, state in [('', 'NO_HAMMING'), ('_HAMMING', 'WITH_HAMMING')]:
                    key = f"{algo}{suffix}"
                    if key in db_results:
                        filename = create_single_plot(
                            coordinates, variants, db_results[key]['labels'],
                            algo, state, db_basename, db_results[key]['time']
                        )
                        print(f"  ✓ {filename}")

            time_file, mem_file = create_bar_plots(db_results, algorithms, db_basename)
            print(f"\n  ✓ {time_file}")
            print(f"  ✓ {mem_file}")
            
        except Exception as e:
            print(f"Error creating plots: {e}")
            import traceback
            traceback.print_exc()
    print(f"{'✓ BENCHMARK COMPLETE':^90}")

if __name__ == "__main__":
    main()