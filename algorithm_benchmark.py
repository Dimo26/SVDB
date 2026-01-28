
#!/usr/bin/env python3

import time
import sys
import os
import glob
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import psutil
from collections import defaultdict


def load_database(db_file, chromosome=None):
    """
    Load variants from SVDB database.
    
    Args:
        db_file: Path to .db file
        chromosome: Specific chromosome to query (None = all chromosomes)
    
    Returns:
        coordinates, variants, stats
    """
    try:
        sys.path.insert(0, '')
        from database import DB
        
        db = DB(db_file)
        coordinates, variants = [], []
        
        db_name = os.path.basename(db_file)
        print(f"\n{'='*90}")
        print(f"Loading: {db_name}")
        
        # Build query based on chromosome filter
        if chromosome:
            query = f'SELECT * FROM SVDB WHERE chrA = "{chromosome}" AND chrB = "{chromosome}"'
            print(f"Filtering: Chromosome {chromosome} only")
        else:
            query = 'SELECT * FROM SVDB'
            print(f"Loading: ALL chromosomes")
        
        for row in db.query(query):
            var_type = row[0]
            chrA, chrB = row[1], row[2]
            posA, posB = int(row[3]), int(row[6])
            seq = row[11] if len(row) > 11 else ''
            sample = row[9]
            
            # Calculate SV size
            if var_type == 'INS' and seq:
                sv_size = len(seq)
            else:
                sv_size = abs(posB - posA)
            
            coordinates.append([posA, posB])
            variants.append({
                'chrA': chrA, 'chrB': chrB,
                'posA': posA, 'posB': posB,
                'type': var_type,
                'sequence': seq,
                'size': sv_size,
                'sample': sample
            })
        
        print(f"✓ Loaded {len(variants):,} variants")
        
        # Print statistics
        stats = calculate_sv_statistics(variants)
        print_sv_statistics(stats, variants)
        
        return np.array(coordinates), variants, stats
    
    except Exception as e:
        print(f"ERROR loading database: {e}")
        import traceback
        traceback.print_exc()
        return None, None, None


def calculate_sv_statistics(variants):
    """Calculate comprehensive statistics for variants."""
    stats = defaultdict(lambda: {
        'count': 0,
        'sizes': [],
        'with_seq': 0,
        'without_seq': 0,
        'by_chr': defaultdict(int)
    })
    
    for v in variants:
        sv_type = v['type']
        stats[sv_type]['count'] += 1
        stats[sv_type]['sizes'].append(v['size'])
        stats[sv_type]['by_chr'][v['chrA']] += 1
        
        if sv_type == 'INS':
            if v['sequence'] and len(v['sequence']) > 0:
                stats[sv_type]['with_seq'] += 1
            else:
                stats[sv_type]['without_seq'] += 1
    
    return dict(stats)


def print_sv_statistics(stats, variants):
    """Print formatted SV statistics table."""
    if not stats:
        print("No variants found")
        return
    
    print(f"\n{'SV TYPE':<10} {'COUNT':>10} {'%':>8} {'AVG SIZE':>12} {'MIN':>12} {'MAX':>15}")
    print("="*70)
    
    total = len(variants)
    for sv_type in sorted(stats.keys()):
        data = stats[sv_type]
        pct = (data['count'] / total * 100) if total > 0 else 0
        
        if data['sizes']:
            avg_size = np.mean(data['sizes'])
            min_size = min(data['sizes'])
            max_size = max(data['sizes'])
        else:
            avg_size = min_size = max_size = 0
        
        print(f"{sv_type:<10} {data['count']:>10,} {pct:>7.1f}% "
              f"{avg_size:>11.1f} {min_size:>12,} {max_size:>15,}")
    
    print(f"{'TOTAL':<10} {total:>10,} {'100.0%':>8}")
    
    # Insertion details
    if 'INS' in stats:
        ins = stats['INS']
        print(f"\n{'INSERTION DETAILS':^70}")
        print(f"  Total insertions:     {ins['count']:>8,}")
        if ins['count'] > 0:
            pct_with = (ins['with_seq'] / ins['count'] * 100)
            pct_without = (ins['without_seq'] / ins['count'] * 100)
            print(f"  With sequence:        {ins['with_seq']:>8,} ({pct_with:>5.1f}%)")
            print(f"  Without sequence:     {ins['without_seq']:>8,} ({pct_without:>5.1f}%)")


# ============================================================================
# CLUSTERING ALGORITHMS
# ============================================================================

def dbscan_cluster(coordinates, epsilon=500, min_pts=2):
    """DBSCAN clustering implementation."""
    sys.path.insert(0, '/mnt/project')
    from export_module import DBSCAN
    return DBSCAN.cluster(coordinates, epsilon, min_pts)


def optics_cluster(coordinates, min_samples=2, max_eps=2000):
    """OPTICS clustering implementation."""
    sys.path.insert(0, '/mnt/project')
    from optics_clustering import optics_cluster
    return optics_cluster(coordinates, min_samples=min_samples, max_eps=max_eps)


def interval_tree_cluster(coordinates, max_distance=1000):
    """Interval tree clustering implementation."""
    sys.path.insert(0, '/mnt/project')
    from interval_tree_overlap import interval_tree_cluster
    return interval_tree_cluster(coordinates, max_distance=max_distance)


def overlap_cluster(coordinates, variants, distance=500, overlap=0.6):
    """Overlap-based clustering (simpler baseline)."""
    sys.path.insert(0, '/mnt/project')
    from overlap_module import isSameVariation
    
    n = len(coordinates)
    labels = np.full(n, -1)
    cluster_id = 0
    
    for i in range(n):
        if labels[i] != -1:
            continue
        
        labels[i] = cluster_id
        for j in range(i + 1, n):
            if labels[j] != -1:
                continue
            
            posA_i, posB_i = coordinates[i]
            posA_j, posB_j = coordinates[j]
            
            dist_A = abs(posA_i - posA_j)
            dist_B = abs(posB_i - posB_j)
            
            if dist_A <= distance and dist_B <= distance:
                _, match = isSameVariation(posA_i, posB_i, posA_j, posB_j, overlap, distance)
                if match:
                    labels[j] = cluster_id
        
        cluster_id += 1
    
    return labels


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
    """Apply Hamming distance re-clustering for insertions."""
    if labels is None:
        return labels
    
    labels = np.array(labels)
    new_labels = np.full_like(labels, -1)
    next_cluster_id = 0
    
    unique_spatial_labels = sorted(set(labels.tolist()))
    
    for spatial_label in unique_spatial_labels:
        if spatial_label == -1:
            continue
        
        indices = np.where(labels == spatial_label)[0]
        if len(indices) == 0:
            continue
        
        # Separate insertions with sequences from others
        ins_with_seq = []
        other_indices = []
        
        for idx in indices:
            var = variants[idx]
            if var['type'] == 'INS' and var['sequence']:
                ins_with_seq.append(idx)
            else:
                other_indices.append(idx)
        
        # Non-insertion variants keep their cluster
        if other_indices:
            new_labels[other_indices] = next_cluster_id
            next_cluster_id += 1
        
        # Re-cluster insertions by sequence similarity
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
                    
                    dist = hamming_distance(seq_i, seq_j)
                    if dist <= max_hamming:
                        group.append(j)
                        assigned.add(j)
                
                for g in group:
                    new_labels[ins_with_seq[g]] = next_cluster_id
                next_cluster_id += 1
    
    return new_labels


# ============================================================================
# BENCHMARKING
# ============================================================================

def benchmark_algorithm(coordinates, variants, algorithm, use_hamming=False, params=None):
    """
    Benchmark a single clustering algorithm.
    
    Returns dict with: time, memory, labels, clusters, clustered, noise
    """
    if params is None:
        params = {}
    
    process = psutil.Process()
    mem_before = process.memory_info().rss
    
    start_time = time.time()
    
    try:
        # Run spatial clustering
        if algorithm == 'DBSCAN':
            labels = dbscan_cluster(coordinates, **params)
        elif algorithm == 'OPTICS':
            labels = optics_cluster(coordinates, **params)
        elif algorithm == 'INTERVAL_TREE':
            labels = interval_tree_cluster(coordinates, **params)
        elif algorithm == 'OVERLAP':
            labels = overlap_cluster(coordinates, variants, **params)
        else:
            raise ValueError(f"Unknown algorithm: {algorithm}")
        
        # Apply Hamming re-clustering if requested
        if use_hamming and labels is not None:
            labels = apply_hamming_reclustering(labels, variants, max_hamming=0.2)
        
        elapsed_time = time.time() - start_time
        mem_after = process.memory_info().rss
        memory_used = mem_after - mem_before
        
        # Calculate statistics
        if labels is not None:
            n_clustered = np.sum(labels != -1)
            n_noise = np.sum(labels == -1)
            n_clusters = len(set(labels)) - (1 if -1 in labels else 0)
        else:
            n_clustered = n_noise = n_clusters = 0
            labels = np.full(len(coordinates), -1)
        
        return {
            'time': elapsed_time,
            'memory': memory_used,
            'labels': labels,
            'clusters': n_clusters,
            'clustered': n_clustered,
            'noise': n_noise
        }
    
    except Exception as e:
        print(f"  ERROR: {e}")
        import traceback
        traceback.print_exc()
        return None


# ============================================================================
# SCALABILITY ANALYSIS
# ============================================================================

def analyze_scalability(db_files, algorithms, use_hamming=False, chromosome='1'):
    """
    Analyze algorithm performance across different database sizes.
    
    Args:
        db_files: List of database files with different sample counts
        algorithms: List of algorithm names to test
        use_hamming: Whether to use Hamming re-clustering
        chromosome: Which chromosome to analyze (for consistency)
    
    Returns:
        results dict with db_size -> algorithm -> metrics
    """
    print(f"\n{'='*90}")
    print(f"SCALABILITY ANALYSIS - {'WITH' if use_hamming else 'WITHOUT'} Hamming")
    print(f"{'='*90}")
    
    results = {}
    
    # Extract sample count from database filename
    def get_sample_count(db_file):
        basename = os.path.basename(db_file)
        if 'all_samples' in basename:
            return 10000  # Placeholder for "all"
        for part in basename.split('_'):
            if part.replace('samples', '').isdigit():
                return int(part.replace('samples', ''))
        return 0
    
    # Sort databases by sample count
    db_files_sorted = sorted(db_files, key=get_sample_count)
    
    for db_file in db_files_sorted:
        sample_count = get_sample_count(db_file)
        db_name = os.path.basename(db_file)
        
        print(f"\n{'─'*90}")
        print(f"Database: {db_name} (n={sample_count})")
        print(f"{'─'*90}")
        
        # Load database
        coordinates, variants, stats = load_database(db_file, chromosome=chromosome)
        if coordinates is None or len(coordinates) == 0:
            print(f"  ⚠ Skipping {db_name} - no data")
            continue
        
        n_variants = len(variants)
        results[sample_count] = {
            'n_variants': n_variants,
            'algorithms': {}
        }
        
        # Benchmark each algorithm
        for algo in algorithms:
            print(f"\n  {algo}...", end=' ', flush=True)
            
            result = benchmark_algorithm(coordinates, variants, algo, use_hamming)
            
            if result:
                results[sample_count]['algorithms'][algo] = result
                print(f"✓ {result['time']:.3f}s | "
                      f"{result['clusters']} clusters | "
                      f"{result['noise']} noise")
            else:
                print(f"✗ FAILED")
    
    return results


def plot_scalability_curves(results, algorithms, output_prefix, use_hamming=False):
    """
    Create time complexity curves showing algorithm performance vs database size.
    
    Creates plots:
    1. Time vs Number of Samples
    2. Time vs Number of Variants
    3. Memory vs Number of Samples
    """
    if not results:
        print("No results to plot")
        return []
    
    sample_counts = sorted(results.keys())
    
    # Prepare data for plotting
    plot_data = {algo: {
        'samples': [],
        'variants': [],
        'times': [],
        'memory': []
    } for algo in algorithms}
    
    for sample_count in sample_counts:
        if sample_count not in results:
            continue
        
        n_variants = results[sample_count]['n_variants']
        
        for algo in algorithms:
            if algo in results[sample_count]['algorithms']:
                result = results[sample_count]['algorithms'][algo]
                plot_data[algo]['samples'].append(sample_count)
                plot_data[algo]['variants'].append(n_variants)
                plot_data[algo]['times'].append(result['time'])
                plot_data[algo]['memory'].append(result['memory'] / 1024 / 1024)  # MB
    
    # Color scheme
    colors = {
        'DBSCAN': '#FF6B6B',
        'OPTICS': '#4ECDC4',
        'INTERVAL_TREE': '#45B7D1',
        'OVERLAP': '#96CEB4'
    }
    
    hamming_label = " (WITH Hamming)" if use_hamming else " (NO Hamming)"
    
    output_files = []
    
    # Plot 1: Time vs Number of Samples
    fig, ax = plt.subplots(figsize=(12, 7))
    
    for algo in algorithms:
        if plot_data[algo]['samples']:
            ax.plot(plot_data[algo]['samples'], plot_data[algo]['times'],
                   marker='o', linewidth=2.5, markersize=8,
                   color=colors.get(algo, '#333333'),
                   label=algo)
    
    ax.set_xlabel('Number of Samples', fontsize=14, fontweight='bold')
    ax.set_ylabel('Time (seconds)', fontsize=14, fontweight='bold')
    ax.set_title(f'Algorithm Performance vs Database Size{hamming_label}',
                 fontsize=16, fontweight='bold', pad=20)
    ax.legend(fontsize=12, frameon=True, shadow=True)
    ax.set_xscale('log')
    ax.set_yscale('log')
    
    plt.tight_layout()
    filename = f"{output_prefix}_scalability_time_vs_samples{'_hamming' if use_hamming else ''}.png"
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.close()
    output_files.append(filename)
    print(f"  ✓ {filename}")
    
    # Plot 2: Time vs Number of Variants
    fig, ax = plt.subplots(figsize=(12, 7))
    
    for algo in algorithms:
        if plot_data[algo]['variants']:
            ax.plot(plot_data[algo]['variants'], plot_data[algo]['times'],
                   marker='s', linewidth=2.5, markersize=8,
                   color=colors.get(algo, '#333333'),
                   label=algo)
    
    ax.set_xlabel('Number of Variants in chromosome 1', fontsize=14, fontweight='bold')
    ax.set_ylabel('Time (seconds)', fontsize=14, fontweight='bold')
    ax.set_title(f'Algorithm Performance vs Variant Count{hamming_label}',
                 fontsize=16, fontweight='bold', pad=20)
    ax.legend(fontsize=12, frameon=True, shadow=True)
    ax.grid(True, alpha=0.3, linestyle='--')
    ax.set_xscale('log')
    ax.set_yscale('log')
    
    plt.tight_layout()
    filename = f"{output_prefix}_scalability_time_vs_variants{'_hamming' if use_hamming else ''}.png"
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.close()
    output_files.append(filename)
    print(f"  ✓ {filename}")
    
    # Plot 3: Memory vs Number of Samples
    fig, ax = plt.subplots(figsize=(12, 7))
    
    for algo in algorithms:
        if plot_data[algo]['samples']:
            ax.plot(plot_data[algo]['samples'], plot_data[algo]['memory'],
                   marker='D', linewidth=2.5, markersize=8,
                   color=colors.get(algo, '#333333'),
                   label=algo)
    
    ax.set_xlabel('Number of Samples', fontsize=14, fontweight='bold')
    ax.set_ylabel('Memory Usage (MB)', fontsize=14, fontweight='bold')
    ax.set_title(f'Memory Usage vs Database Size{hamming_label}',
                 fontsize=16, fontweight='bold', pad=20)
    ax.legend(fontsize=12, frameon=True, shadow=True)
    ax.set_xscale('log')
    
    plt.tight_layout()
    filename = f"{output_prefix}_scalability_memory_vs_samples{'_hamming' if use_hamming else ''}.png"
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.close()
    output_files.append(filename)
    print(f"  ✓ {filename}")
    
    return output_files


# ============================================================================
# SV SIZE DISTRIBUTION PLOTS
# ============================================================================

def create_sv_distribution_plot(variants, output_prefix):
    """
    Create comprehensive SV size distribution plot (like the example image).
    Separate subplot for each SV type showing distribution.
    """
    # Group variants by type
    sv_by_type = defaultdict(list)
    for v in variants:
        sv_by_type[v['type']].append(v['size'])
    
    # Create subplots (2x2 or 2x3 depending on number of types)
    n_types = len(sv_by_type)
    if n_types <= 4:
        nrows, ncols = 2, 2
    else:
        nrows, ncols = 2, 3
    
    fig, axes = plt.subplots(nrows, ncols, figsize=(16, 10))
    axes = axes.flatten()
    
    colors = {
        'BND': '#8B008B',
        'DEL': '#32CD32', 
        'DUP': '#FFD700',
        'INV': '#FF4500',
        'INS': '#4169E1'
    }
    
    for idx, (sv_type, sizes) in enumerate(sorted(sv_by_type.items())):
        if idx >= len(axes):
            break
        
        ax = axes[idx]
        color = colors.get(sv_type, '#808080')
        
        # Create histogram
        counts, bins, patches = ax.hist(sizes, bins=50, color=color, 
                                       edgecolor='black', alpha=0.7, linewidth=0.5)
        
        # Calculate statistics
        mean_size = np.mean(sizes)
        median_size = np.median(sizes)
        min_size = np.min(sizes)
        max_size = np.max(sizes)
        
        # Format axis
        ax.set_xlabel('SV Size (bp)', fontsize=12, fontweight='bold')
        ax.set_ylabel('Count', fontsize=12, fontweight='bold')
        ax.set_title(f'{sv_type} (n={len(sizes):,})', 
                    fontsize=14, fontweight='bold', pad=10)
        
        # Add statistics box
        stats_text = f'Mean: {mean_size:,.0f} bp\n'
        stats_text += f'Median: {median_size:,.0f} bp\n'
        stats_text += f'Min: {min_size:,.0f} bp\n'
        stats_text += f'Max: {max_size:,.0f} bp'
        
        ax.text(0.97, 0.97, stats_text,
               transform=ax.transAxes,
               fontsize=10,
               verticalalignment='top',
               horizontalalignment='right',
               bbox=dict(boxstyle='round', facecolor='white', alpha=0.8, edgecolor='black'))
        
        
        # Use log scale for y-axis if range is large
        if max(counts) > 100:
            ax.set_yscale('log')
    
    # Turn off unused subplots
    for idx in range(len(sv_by_type), len(axes)):
        axes[idx].axis('off')
    
    plt.suptitle('SV Size Distribution by Type', 
                fontsize=18, fontweight='bold', y=0.995)
    plt.tight_layout()
    
    filename = f"{output_prefix}_sv_size_distribution.png"
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.close()
    
    return filename


# ============================================================================
# MAIN FUNCTION
# ============================================================================

def main():
    """Main benchmark orchestration."""
    
    if len(sys.argv) < 2:
        print("Usage:")
        print("  Single database:   python svdb_comprehensive_benchmark.py <database.db>")
        print("  Multiple databases: python svdb_comprehensive_benchmark.py <db1.db> <db2.db> ...")
        print("  Scalability test:  python svdb_comprehensive_benchmark.py --scalability <db_dir>")
        sys.exit(1)
    
    # Check if scalability analysis mode
    if '--scalability' in sys.argv:
        scalability_mode = True
        # Find all database files in specified directory
        db_dir = sys.argv[sys.argv.index('--scalability') + 1]
        
        if not os.path.isdir(db_dir):
            print(f"ERROR: {db_dir} is not a directory")
            sys.exit(1)
        
        db_files = glob.glob(os.path.join(db_dir, '*samples*.db'))
        db_files = sorted(db_files)
        
        if not db_files:
            print(f"ERROR: No *samples*.db files found in {db_dir}")
            sys.exit(1)
        
        print(f"{'='*90}")
        print(f"SVDB SCALABILITY ANALYSIS")
        print(f"{'='*90}")
        print(f"Database directory: {db_dir}")
        print(f"Found {len(db_files)} databases:")
        for f in db_files:
            print(f"  • {os.path.basename(f)}")
        
    else:
        scalability_mode = False
        # Single/multiple database mode
        db_files = []
        for arg in sys.argv[1:]:
            if '*' in arg:
                db_files.extend(glob.glob(arg))
            elif os.path.exists(arg):
                db_files.append(arg)
        
        if not db_files:
            print("ERROR: No valid database files found")
            sys.exit(1)
        
        print(f"{'='*90}")
        print(f"SVDB COMPREHENSIVE BENCHMARK")
        print(f"{'='*90}")
        print(f"Databases: {len(db_files)}")
        for f in db_files:
            print(f"  • {os.path.basename(f)}")
    
    # Define algorithms to test
    algorithms = ['DBSCAN', 'OPTICS', 'INTERVAL_TREE', 'OVERLAP']
    
    if scalability_mode:
        # Scalability analysis mode
        output_prefix = os.path.join(os.getcwd(), 'scalability_analysis')
        
        # Run analysis WITHOUT Hamming
        print(f"\n{'#'*90}")
        print(f"PHASE 1: Scalability WITHOUT Hamming")
        print(f"{'#'*90}")
        results_no_hamming = analyze_scalability(
            db_files, algorithms, use_hamming=False, chromosome='1'
        )
        
        # Run analysis WITH Hamming
        print(f"\n{'#'*90}")
        print(f"PHASE 2: Scalability WITH Hamming")
        print(f"{'#'*90}")
        results_with_hamming = analyze_scalability(
            db_files, algorithms, use_hamming=True, chromosome='1'
        )
        
        # Generate plots
        print(f"\n{'='*90}")
        print(f"GENERATING SCALABILITY PLOTS")
        print(f"{'='*90}")
        
        plot_files = []
        plot_files.extend(plot_scalability_curves(
            results_no_hamming, algorithms, output_prefix, use_hamming=False
        ))
        plot_files.extend(plot_scalability_curves(
            results_with_hamming, algorithms, output_prefix, use_hamming=True
        ))
        
        print(f"\n{'='*90}")
        print(f"✓ SCALABILITY ANALYSIS COMPLETE")
        print(f"{'='*90}")
        print(f"Generated {len(plot_files)} plots")
        
    else:
        # Single/multiple database mode - detailed analysis
        for db_file in db_files:
            db_basename = os.path.basename(db_file).replace('.db', '')
            
            print(f"\n{'#'*90}")
            print(f"PROCESSING: {db_basename}")
            print(f"{'#'*90}")
            
            # Load database (ALL chromosomes for accurate stats)
            coordinates, variants, stats = load_database(db_file, chromosome=None)
            if coordinates is None or len(coordinates) == 0:
                print(f"  ⚠ Skipping {db_basename} - no data")
                continue
            
            # Create SV distribution plot
            print(f"\n{'─'*90}")
            print(f"GENERATING SV DISTRIBUTION PLOT")
            print(f"{'─'*90}")
            dist_plot = create_sv_distribution_plot(variants, db_basename)
            print(f"  ✓ {dist_plot}")
            
            # For algorithm benchmarking, use only Chr1 for speed
            print(f"\n{'─'*90}")
            print(f"ALGORITHM BENCHMARKING (Chromosome 1 only)")
            print(f"{'─'*90}")
            coords_chr1, vars_chr1, _ = load_database(db_file, chromosome='1')
            
            if coords_chr1 is None or len(coords_chr1) == 0:
                print(f"  ⚠ No Chr1 data for benchmarking")
                continue
            
            # Benchmark algorithms
            results = {}
            for algo in algorithms:
                print(f"\n  {algo}...")
                
                # WITHOUT Hamming
                result = benchmark_algorithm(coords_chr1, vars_chr1, algo, use_hamming=False)
                if result:
                    results[algo] = result
                    print(f"    NO Hamming:  {result['time']:.4f}s | "
                          f"{result['clusters']:>4} clusters | {result['noise']:>5} noise")
                
                # WITH Hamming
                result = benchmark_algorithm(coords_chr1, vars_chr1, algo, use_hamming=True)
                if result:
                    results[f"{algo}_HAMMING"] = result
                    print(f"    WITH Hamming: {result['time']:.4f}s | "
                          f"{result['clusters']:>4} clusters | {result['noise']:>5} noise")
            
            print(f"\n{'='*90}")
            print(f"✓ PROCESSING COMPLETE: {db_basename}")
            print(f"{'='*90}")
    
    print(f"\n{'█'*90}")
    print(f"{'✓ BENCHMARK COMPLETE':^90}")
    print(f"{'█'*90}")


if __name__ == "__main__":
    main()