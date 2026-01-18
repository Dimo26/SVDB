#!/usr/bin/env python3
"""
Algorithm Benchmark for SVDB Clustering
Compares: DBSCAN, OPTICS, INTERVAL_TREE with/without Hamming distance
"""

import time
import sys
import os
import glob
import numpy as np
import psutil


def calculate_sv_statistics(variants):
    """Calculate comprehensive SV statistics from variant list."""
    if not variants:
        return None
    
    stats = {}
    for v in variants:
        sv_type = v['type']
        if sv_type not in stats:
            stats[sv_type] = {
                'count': 0,
                'sizes': [],
                'with_seq': 0,
                'without_seq': 0
            }
        
        stats[sv_type]['count'] += 1
        stats[sv_type]['sizes'].append(v['size'])
        
        if sv_type == 'INS':
            if v['sequence'] and len(v['sequence']) > 0:
                stats[sv_type]['with_seq'] += 1
            else:
                stats[sv_type]['without_seq'] += 1
    
    return stats


def print_statistics(stats, title="SV Statistics"):
    """Print formatted statistics table."""
    if not stats:
        print("No statistics available")
        return
    
    print(f"\n{'='*80}")
    print(f"{title}")
    print(f"{'='*80}")
    print(f"{'SV_TYPE':<8} {'COUNT':>8} {'AVG_SIZE':>12} {'MIN_SIZE':>12} "
          f"{'MAX_SIZE':>12} {'WITH_SEQ':>10} {'NO_SEQ':>10}")
    print(f"{'-'*8} {'-'*8} {'-'*12} {'-'*12} {'-'*12} {'-'*10} {'-'*10}")
    
    total_count = 0
    total_with_seq = 0
    total_no_seq = 0
    
    for sv_type in sorted(stats.keys()):
        data = stats[sv_type]
        sizes = data['sizes']
        
        avg_size = np.mean(sizes) if sizes else 0
        min_size = min(sizes) if sizes else 0
        max_size = max(sizes) if sizes else 0
        
        print(f"{sv_type:<8} {data['count']:>8} {avg_size:>12.1f} {min_size:>12} "
              f"{max_size:>12} {data['with_seq']:>10} {data['without_seq']:>10}")
        
        total_count += data['count']
        total_with_seq += data['with_seq']
        total_no_seq += data['without_seq']
    
    print(f"{'-'*8} {'-'*8} {'-'*12} {'-'*12} {'-'*12} {'-'*10} {'-'*10}")
    print(f"{'TOTAL':<8} {total_count:>8} {'':>12} {'':>12} {'':>12} "
          f"{total_with_seq:>10} {total_no_seq:>10}")
    
    # Insertion details
    if 'INS' in stats:
        ins_data = stats['INS']
        print(f"\n--- Insertion Details ---")
        print(f"  Total: {ins_data['count']}")
        print(f"  With sequence: {ins_data['with_seq']} "
              f"({ins_data['with_seq']/ins_data['count']*100:.1f}%)")
        print(f"  Without sequence: {ins_data['without_seq']} "
              f"({ins_data['without_seq']/ins_data['count']*100:.1f}%)")
        
        seq_sizes = [s for s in ins_data['sizes'] if s > 0]
        if seq_sizes:
            print(f"\n  Sequence length stats (insertions with sequence):")
            print(f"    Mean: {np.mean(seq_sizes):.1f} bp")
            print(f"    Median: {np.median(seq_sizes):.1f} bp")
            print(f"    Min: {min(seq_sizes)} bp")
            print(f"    Max: {max(seq_sizes)} bp")
            print(f"    Std: {np.std(seq_sizes):.1f} bp")
    
    print(f"{'='*80}\n")


def load_db_samples(db_file):
    """Load variants from SVDB database with proper size calculation."""
    try:
        from svdb import database
        db = database.DB(db_file)
        
        coordinates = []
        variants = []
        
        # Query chromosome 1
        for row in db.query('SELECT * FROM SVDB WHERE chrA = "1" AND chrB = "1"'):
            var_type = row[0]
            chrA = row[1]  
            chrB = row[2]
            posA = int(row[3])
            ci_A_lower = int(row[4])
            ci_A_upper = int(row[5])
            posB = int(row[6])
            ci_B_lower = int(row[7])
            ci_B_upper = int(row[8])
            sample = row[9]
            idx = int(row[10])
            seq = row[11] if len(row) > 11 else ''
            
            # Calculate proper SV size
            if var_type == 'INS':
                sv_size = len(seq) if seq else 0
            elif var_type == 'INV':
                sv_size = abs(posB - posA)
            else:
                sv_size = abs(posB - posA)
            
            coordinates.append([posA, posB])
            variants.append({
                'posA': posA,
                'posB': posB,
                'type': var_type,
                'sequence': seq,
                'size': sv_size,
                'sample': sample,
                'ci_A_lower': ci_A_lower,
                'ci_A_upper': ci_A_upper,
                'ci_B_lower': ci_B_lower,
                'ci_B_upper': ci_B_upper
            })
        
        print(f"✓ Loaded {len(variants)} variants from chromosome 1")
        
        # Calculate and print statistics
        stats = calculate_sv_statistics(variants)
        print_statistics(stats, "Database Statistics (Chromosome 1)")
        
        return (np.array(coordinates) if coordinates else None, variants, stats)
    
    except Exception as e:
        print(f"Error loading database {db_file}: {e}")
        import traceback
        traceback.print_exc()
        return None, None, None


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
    """
    Apply Hamming distance re-clustering for insertions WITHIN spatial clusters.
    - Only affects insertions WITH sequences
    - Insertions WITHOUT sequences keep their spatial cluster
    - Noise points remain noise regardless of sequence
    """
    if labels is None or variants is None:
        return labels
    
    labels = np.array(labels)
    new_labels = np.full_like(labels, -1)
    next_cluster_id = 0
    
    unique_spatial_labels = sorted(set(labels.tolist()))
    
    ins_with_seq_processed = 0
    ins_without_seq = 0
    seq_clusters_created = 0

    for spatial_label in unique_spatial_labels:
        if spatial_label == -1:
            continue
        
        indices = np.where(labels == spatial_label)[0]
        if len(indices) == 0:
            continue
        
        ins_with_seq = []
        other_indices = []
        
        for idx in indices:
            if variants[idx]['type'] == 'INS':
                if variants[idx]['sequence']:
                    ins_with_seq.append(idx)
                else:
                    # Insertions without sequence stay in spatial cluster
                    other_indices.append(idx)
                    ins_without_seq += 1
            else:
                other_indices.append(idx)
        
        # Non-insertions and insertions-without-sequence keep spatial cluster
        if other_indices:
            new_labels[other_indices] = next_cluster_id
            next_cluster_id += 1
        
        # Re-cluster insertions WITH sequences by sequence similarity
        if ins_with_seq:
            ins_with_seq_processed += len(ins_with_seq)
            assigned = set()
            cluster_count_before = next_cluster_id
            
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
            
            seq_clusters_created += (next_cluster_id - cluster_count_before)
    
    # All noise points stay as noise
    noise_indices = np.where(labels == -1)[0]
    for idx in noise_indices:
        new_labels[idx] = -1
    
    print(f"  → Hamming re-clustering:")
    print(f"     • {ins_with_seq_processed} insertions WITH sequences → {seq_clusters_created} clusters")
    print(f"     • {ins_without_seq} insertions WITHOUT sequences → kept spatial clusters")
    
    return new_labels


def benchmark_clustering_algorithm(coordinates, variants, algorithm_name, 
                                   apply_hamming=False, max_hamming=0.2):
    """Benchmark a single clustering algorithm."""
    from svdb.export_module import DBSCAN
    from svdb.optics_clustering import OPTICS
    from svdb.interval_tree_overlap import interval_tree_cluster
    from svdb import overlap_module
    
    if coordinates is None or len(coordinates) == 0:
        return None, None, 0, None, 0, 0
    
    distance_threshold = 500
    proc = psutil.Process()
    
    mem_before = proc.memory_info().rss
    start_time = time.time()
    
    try:
        if algorithm_name == 'DBSCAN':
            labels = DBSCAN.cluster(coordinates, distance_threshold, 2)
        elif algorithm_name == 'OPTICS':
            optics = OPTICS(min_samples=2, max_eps=distance_threshold)
            labels = optics.fit_predict(coordinates)
        elif algorithm_name == 'INTERVAL_TREE':
            labels = interval_tree_cluster(coordinates, distance_threshold)
        elif algorithm_name == 'OVERLAP':
            labels = np.full(len(coordinates), -1)
            cluster_id = 0
            for i in range(len(coordinates)):
                if labels[i] != -1:
                    continue
                labels[i] = cluster_id
                for j in range(i+1, len(coordinates)):
                    if labels[j] == -1:
                        similar, match = overlap_module.isSameVariation(
                            coordinates[i][0], coordinates[i][1],
                            coordinates[j][0], coordinates[j][1],
                            0.6, distance_threshold
                        )
                        if match:
                            labels[j] = cluster_id
                cluster_id += 1
        else:
            return None, None, 0, None, 0, 0
        
        if apply_hamming:
            labels = apply_hamming_post_processing(labels, variants, max_hamming)
        
        elapsed = time.time() - start_time
        mem_after = proc.memory_info().rss
        memory_used = mem_after - mem_before
        
        n_clusters = len(set(labels)) - (1 if -1 in labels else 0)
        n_clustered = np.sum(labels != -1)
        n_unclustered = np.sum(labels == -1)
        
        return elapsed, memory_used, n_clusters, labels, n_clustered, n_unclustered
    
    except Exception as e:
        print(f"Error in {algorithm_name}: {e}")
        import traceback
        traceback.print_exc()
        return None, None, 0, None, 0, 0


def determine_smart_zoom(coordinates, variants):
    """
    Determine smart zoom ranges to show all SV types clearly.
    Returns (zoom_x, zoom_y, zoom_regions) where zoom_regions is a list of interesting areas.
    """
    if coordinates is None or len(coordinates) == 0:
        return None, None, []
    
    # Separate insertions (posA ≈ posB) from other SVs
    ins_coords = []
    other_coords = []
    
    for coord, var in zip(coordinates, variants):
        if var['type'] == 'INS':
            ins_coords.append(coord)
        else:
            other_coords.append(coord)
    
    # Define zoom regions
    zoom_regions = []
    
    # Full view
    all_x = coordinates[:, 0]
    all_y = coordinates[:, 1]
    zoom_regions.append({
        'name': 'full',
        'x': (np.min(all_x), np.max(all_x)),
        'y': (np.min(all_y), np.max(all_y)),
        'description': 'Full chromosome view'
    })
    
    # If we have insertions, create a zoom for them
    if ins_coords:
        ins_array = np.array(ins_coords)
        ins_x = ins_array[:, 0]
        ins_y = ins_array[:, 1]
        
        # Insertions appear on diagonal (posA ≈ posB)
        margin = max((np.max(ins_x) - np.min(ins_x)) * 0.1, 1000)
        zoom_regions.append({
            'name': 'insertions',
            'x': (np.min(ins_x) - margin, np.max(ins_x) + margin),
            'y': (np.min(ins_y) - margin, np.max(ins_y) + margin),
            'description': 'Insertion-rich region'
        })
    
    # If we have other SVs, create zoom for them
    if other_coords:
        other_array = np.array(other_coords)
        other_x = other_array[:, 0]
        other_y = other_array[:, 1]
        
        margin = max((np.max(other_x) - np.min(other_x)) * 0.1, 1000)
        zoom_regions.append({
            'name': 'other_svs',
            'x': (np.min(other_x) - margin, np.max(other_x) + margin),
            'y': (np.min(other_y) - margin, np.max(other_y) + margin),
            'description': 'DEL/DUP/INV/BND region'
        })
    
    return zoom_regions


def main():
    """Main benchmark function."""
    if len(sys.argv) < 2:
        print("Usage: python algorithm_benchmark.py <database.db> [--no-zoom] [--min-title]")
        print("\nOptions:")
        print("  --no-zoom    : Skip zoom views (only full chromosome view)")
        print("  --min-title  : Minimal plot titles")
        print("\nExample:")
        print("  python algorithm_benchmark.py SVDB.db")
        print("  python algorithm_benchmark.py SVDB.db --no-zoom --min-title")
        sys.exit(1)
    
    # Parse arguments
    db_files = []
    no_zoom = False
    min_title = False
    
    for arg in sys.argv[1:]:
        if arg == '--no-zoom':
            no_zoom = True
        elif arg == '--min-title':
            min_title = True
        elif '*' in arg or '?' in arg:
            db_files.extend(glob.glob(arg))
        elif not arg.startswith('--'):
            db_files.append(arg)
    
    if not db_files:
        print("Error: No database files provided!")
        sys.exit(1)
    
    # Validate files
    valid_db_files = [f for f in db_files if os.path.exists(f)]
    if not valid_db_files:
        print("Error: No valid database files!")
        sys.exit(1)
    
    print(f"\n{'='*80}")
    print(f"SVDB CLUSTERING ALGORITHM BENCHMARK")
    print(f"{'='*80}")
    print(f"Databases: {len(valid_db_files)}")
    for f in valid_db_files:
        print(f"  • {os.path.basename(f)}")
    print(f"Options: {'No zoom' if no_zoom else 'With zoom'}, "
          f"{'Minimal titles' if min_title else 'Full titles'}")
    print(f"{'='*80}\n")
    
    algorithms = ['DBSCAN', 'OPTICS', 'INTERVAL_TREE']
    results = {}
    
    # Run benchmarks
    for db_file in valid_db_files:
        print(f"\n{'='*80}")
        print(f"Processing: {os.path.basename(db_file)}")
        print(f"{'='*80}")
        
        coordinates, variants, stats = load_db_samples(db_file)
        if coordinates is None or len(coordinates) == 0:
            print(f"  ✗ No data found, skipping...")
            continue
        
        results[db_file] = {'stats': stats}
        
        # Test each algorithm
        for algo in algorithms:
            print(f"\n{algo}:")
            
            # Without Hamming
            elapsed, memory, n_clusters, labels, n_clustered, n_unclustered = benchmark_clustering_algorithm(
                coordinates, variants, algo, apply_hamming=False
            )
            
            if elapsed is not None:
                results[db_file][algo] = {
                    'time': elapsed,
                    'memory': memory,
                    'clusters': n_clusters,
                    'clustered': n_clustered,
                    'unclustered': n_unclustered
                }
                print(f"  Without Hamming: {elapsed:7.4f}s | {n_clusters:4d} clusters | "
                      f"{n_clustered:4d} clustered | {n_unclustered:4d} noise | "
                      f"{memory/1024/1024:6.2f} MB")
            
            # With Hamming
            elapsed, memory, n_clusters, labels, n_clustered, n_unclustered = benchmark_clustering_algorithm(
                coordinates, variants, algo, apply_hamming=True, max_hamming=0.2
            )
            
            if elapsed is not None:
                results[db_file][f"{algo}_HAMMING"] = {
                    'time': elapsed,
                    'memory': memory,
                    'clusters': n_clusters,
                    'clustered': n_clustered,
                    'unclustered': n_unclustered
                }
                print(f"  With Hamming:    {elapsed:7.4f}s | {n_clusters:4d} clusters | "
                      f"{n_clustered:4d} clustered | {n_unclustered:4d} noise | "
                      f"{memory/1024/1024:6.2f} MB")
    
    # Summary
    print(f"\n{'='*80}")
    print("BENCHMARK SUMMARY")
    print(f"{'='*80}")
    
    for db_file in results:
        if 'stats' not in results[db_file]:
            continue
        print(f"\n{os.path.basename(db_file)}:")
        for algo in algorithms:
            without = results[db_file].get(algo, None)
            with_h = results[db_file].get(f"{algo}_HAMMING", None)
            
            if without:
                print(f"\n  {algo}:")
                print(f"    Without Hamming: {without['time']:.4f}s | "
                      f"{without['clusters']} clusters")
                
                if with_h:
                    overhead = (with_h['time'] - without['time']) / without['time'] * 100 if without['time'] > 0 else 0
                    cluster_diff = with_h['clusters'] - without['clusters']
                    print(f"    With Hamming:    {with_h['time']:.4f}s | "
                          f"{with_h['clusters']} clusters")
                    print(f"    Overhead: +{overhead:.1f}% time | {cluster_diff:+d} clusters")
    
    # Plotting
    print(f"\n{'='*80}")
    print("GENERATING VISUALIZATIONS")
    print(f"{'='*80}\n")
    
    try:
        import matplotlib.pyplot as plt
        from matplotlib.patches import Patch
        
        # SV type colors (colorblind-friendly)
        sv_colors = {
            'INS': '#FFD700',  # Gold (Yellow)
            'DEL': '#32CD32',  # Lime (Green)
            'DUP': '#4169E1',  # Royal Blue
            'INV': '#FF4500',  # Orange Red
            'BND': '#8B008B',  # Purple
            'NOISE': '#808080'  # Gray
        }
        
        print("Color scheme:")
        for sv_type, color in sv_colors.items():
            print(f"  {sv_type:<6} = {color}")
        print()
        
        for db_file in results:
            if 'stats' not in results[db_file]:
                continue
            
            coordinates, variants, stats = load_db_samples(db_file)
            if coordinates is None:
                continue
            
            db_basename = os.path.basename(db_file).replace('.db', '')
            
            # Determine zoom regions
            zoom_regions = [] if no_zoom else determine_smart_zoom(coordinates, variants)
            if not zoom_regions:
                zoom_regions = [{'name': 'full', 'x': None, 'y': None, 'description': 'Full view'}]
            
            print(f"Generating plots for {db_basename}:")
            print(f"  Zoom regions: {len(zoom_regions)}")
            for zr in zoom_regions:
                print(f"    • {zr['description']}")
            print()
            
            for algo in algorithms:
                for apply_hamming in [False, True]:
                    elapsed, memory, n_clusters, labels, n_clustered, n_unclustered = benchmark_clustering_algorithm(
                        coordinates, variants, algo, apply_hamming=apply_hamming, max_hamming=0.2
                    )
                    if labels is None:
                        continue
                    
                    hamming_suffix = "WITH_HAMMING" if apply_hamming else "NO_HAMMING"
                    
                    # Generate plot for each zoom region
                    for zoom_idx, zoom_region in enumerate(zoom_regions):
                        fig, ax = plt.subplots(figsize=(14, 10))
                        
                        # Plot each SV
                        for idx, (coord, variant, label) in enumerate(zip(coordinates, variants, labels)):
                            sv_type = variant['type']
                            
                            if label == -1:
                                color = sv_colors['NOISE']
                                alpha = 0.4
                                size = 25
                                marker = 'x'
                                linewidth = 1.0
                            else:
                                color = sv_colors.get(sv_type, '#000000')
                                alpha = 0.75
                                size = 50
                                marker = 'o'
                                linewidth = 0.8
                            
                            ax.scatter(coord[0], coord[1], c=[color], s=size, alpha=alpha, 
                                      edgecolors='black', linewidth=linewidth, marker=marker, zorder=2)
                        
                        # Legend
                        legend_elements = []
                        sv_types_present = set(v['type'] for v in variants)
                        
                        for sv_type in sorted(sv_types_present):
                            count = sum(1 for v in variants if v['type'] == sv_type)
                            legend_elements.append(
                                Patch(facecolor=sv_colors.get(sv_type, '#000000'), 
                                      edgecolor='black', linewidth=1.5,
                                      label=f'{sv_type} (n={count})')
                            )
                        
                        noise_count = sum(1 for l in labels if l == -1)
                        legend_elements.append(
                            Patch(facecolor=sv_colors['NOISE'], 
                                  edgecolor='black', linewidth=1.5,
                                  label=f'Unclustered (n={noise_count})')
                        )
                        
                        ax.legend(handles=legend_elements, loc='upper right', 
                                 title='SV Types', fontsize=11, title_fontsize=12, 
                                 framealpha=0.95, edgecolor='black', fancybox=False)
                        
                        # Apply zoom
                        if zoom_region['x'] is not None:
                            ax.set_xlim(zoom_region['x'])
                        if zoom_region['y'] is not None:
                            ax.set_ylim(zoom_region['y'])
                        
                        # Labels
                        ax.set_xlabel('Position A (bp)', fontsize=13, fontweight='bold')
                        ax.set_ylabel('Position B (bp)', fontsize=13, fontweight='bold')
                        
                        # Title
                        if not min_title:
                            hamming_str = "WITH Hamming" if apply_hamming else "WITHOUT Hamming"
                            title_str = (f'{algo} Clustering {hamming_str}\n'
                                       f'{zoom_region["description"]} | '
                                       f'{n_clusters} clusters | {n_clustered} clustered | '
                                       f'{n_unclustered} noise')
                            ax.set_title(title_str, fontsize=14, fontweight='bold', pad=15)
                        
                        ax.grid(alpha=0.25, linestyle='--', linewidth=0.6, zorder=1)
                        
                        # Save
                        zoom_suffix = f"_{zoom_region['name']}" if len(zoom_regions) > 1 else ""
                        plot_filename = f"{db_basename}_chr1_{algo}_{hamming_suffix}{zoom_suffix}.png"
                        plt.savefig(plot_filename, dpi=300, bbox_inches='tight')
                        plt.close()
                        
                        print(f"  ✓ {plot_filename}")
        
        print(f"\n{'='*80}")
        print("✓ All visualizations generated!")
        print(f"{'='*80}\n")
    
    except ImportError as e:
        print(f"\n✗ Matplotlib not available: {e}")
    except Exception as e:
        print(f"\n✗ Error creating plots: {e}")
        import traceback
        traceback.print_exc()


if __name__ == "__main__":
    main()