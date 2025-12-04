import pysam
from collections import defaultdict
import numpy as np

from .interval_tree_overlap import interval_tree_cluster 
from .export_module import DBSCAN

class SVEvidence:
    def __init__(self, chrom_a, pos_a, chrom_b, pos_b, sv_type, support_reads):
        self.chrA = chrom_a
        self.posA = pos_a
        self.chrB = chrom_b
        self.posB = pos_b
        self.sv_type = sv_type
        self.support_reads = support_reads
        self.ci_a_lower = 0
        self.ci_a_upper = 0
        self.ci_b_lower = 0
        self.ci_b_upper = 0
        self.confidence = None  # Added for storing mapping quality
        self.inserted_sequence = None

def extract_sv_from_split_reads(bam_file, min_mapq=20, min_sv_size=50):
    sv_evidence = []
    try:
        samfile = pysam.AlignmentFile(bam_file, "rb")
    except Exception as e:
        print(f"Error: Could not open BAM file: {e}")
        return sv_evidence
    
    read_alignments = defaultdict(list)
    
    for read in samfile.fetch():
        if read.is_unmapped or read.mapping_quality < min_mapq or read.is_secondary:
            continue
        read_alignments[read.query_name].append(read)
 
        if read.has_tag("SA"):
            sa_tags = read.get_tag("SA").split(";")
            for sa in sa_tags:
                if sa:
                    fields = sa.split(",")
                    if len(fields) >= 6:
                        sa_chrom = fields[0]
                        sa_pos = int(fields[1])
                        sa_strand = fields[2]
                        sa_mapq = int(fields[4])
                        
                        if sa_mapq >= min_mapq:
                            read_alignments[read.query_name].append({
                                'chrom': sa_chrom,
                                'pos': sa_pos,
                                'is_reverse': (sa_strand == '-'),
                                'mapq': sa_mapq,
                                'query_name': read.query_name,
                                'is_sa': True
                            })
    
    for read_name, alignments in read_alignments.items():
        if len(alignments) < 2:
            continue

        def get_sort_key(aln):
            if isinstance(aln, dict):
                return (aln['chrom'], aln['pos'])
            else:
                return (aln.reference_name, aln.reference_start)
        
        alignments.sort(key=get_sort_key)
        
        for i in range(len(alignments) - 1):
            primary = alignments[i]
            supplementary = alignments[i + 1]
            sv_type, evidence = _analyze_alignment_pair(primary, supplementary, min_sv_size)
            if sv_type and evidence:
                sv_evidence.append(evidence)
    
    samfile.close()
    return sv_evidence

def _analyze_alignment_pair(aln1, aln2, min_sv_size):
    # Handle both AlignedSegment and dict types
    if isinstance(aln1, dict):
        chr1 = aln1['chrom']
        pos1 = aln1['pos']
        is_reverse1 = aln1['is_reverse']
        query_name = aln1['query_name']
    else:
        chr1 = aln1.reference_name
        pos1 = aln1.reference_start
        is_reverse1 = aln1.is_reverse
        query_name = aln1.query_name
    
    if isinstance(aln2, dict):
        chr2 = aln2['chrom']
        pos2 = aln2['pos']
        is_reverse2 = aln2['is_reverse']
    else:
        chr2 = aln2.reference_name
        pos2 = aln2.reference_start
        is_reverse2 = aln2.is_reverse

    if chr1 != chr2:
        sv_type = "BND"
        if chr1 > chr2:
            chr1, chr2 = chr2, chr1
            pos1, pos2 = pos2, pos1
        evidence = SVEvidence(chr1, pos1, chr2, pos2, sv_type, [query_name])
        return sv_type, evidence

    distance = abs(pos2 - pos1)
    if distance < min_sv_size:
        return None, None
    
    same_strand = is_reverse1 == is_reverse2
    if same_strand:
        sv_type = "DEL"
        evidence = SVEvidence(chr1, min(pos1, pos2), chr2, max(pos1, pos2), sv_type, [query_name])
        return sv_type, evidence

    return None, None

def extract_sv_from_cigar(bam_file, min_indel_size=50):
    sv_evidence = []
    try:
        samfile = pysam.AlignmentFile(bam_file, "rb")
    except Exception as e:
        print(f"Error: Could not open BAM file: {e}")
        return sv_evidence

    for read in samfile.fetch():
        if read.is_unmapped or read.is_secondary:
            continue
        ref_pos = read.reference_start
        
        soft_clip_total = 0
        for op, length in read.cigartuples:
            if op == 4:  # Soft clip
                soft_clip_total += length
        if soft_clip_total > 0.5 * read.query_length:
            continue
            
        for op, length in read.cigartuples:
            if op == 2 and length >= min_indel_size:  # Deletion
                evidence = SVEvidence(read.reference_name, ref_pos, read.reference_name, ref_pos + length, "DEL", [read.query_name])
                evidence.confidence = read.mapping_quality
                sv_evidence.append(evidence)
                ref_pos += length
            elif op == 1 and length >= min_indel_size:  # Insertion
                read_pos = 0
                for operation, op_length in read.cigartuples:
                    if operation in [0,1,4,7,8]:
                        read_pos += op_length
                ins_start = read_pos - length
                ins_end = read_pos
                inserted_seq = read.query_sequence[ins_start: ins_end] if read.query_sequence else ""
                evidence = SVEvidence(read.reference_name, ref_pos, read.reference_name, ref_pos, "INS", [read.query_name])
                evidence.inserted_sequence = inserted_seq
                sv_evidence.append(evidence)

            elif op in [0, 2, 3, 7, 8]:  # Consumes reference
                ref_pos += length
    samfile.close()
    return sv_evidence

def bam_to_vcf_format(sv_evidence_list, sample_name):
    vcf_like_variants = []
    for evidence in sv_evidence_list:
        INFO = {}
        FORMAT = {'GT': ['0/1']}
        if evidence.chrA == evidence.chrB:
            INFO['SVLEN'] = str(abs(evidence.posB - evidence.posA))
            INFO['END'] = str(evidence.posB)

        INFO['SVTYPE'] = evidence.sv_type
        INFO['support'] = str(len(evidence.support_reads))

        # Include insertion sequence in INFO field for tracking
        if evidence.sv_type == 'INS' and hasattr(evidence, 'inserted_sequence'):
            if evidence.inserted_sequence is not None and len(str(evidence.inserted_sequence)) > 0:
                INFO['INSSEQ'] = str(evidence.inserted_sequence)

        vcf_like_variants.append((evidence.chrA, evidence.posA, evidence.chrB, evidence.posB, evidence.sv_type, INFO, FORMAT, sample_name))
    return vcf_like_variants


def cluster_sv_evidence(sv_evidence_list, distance_threshold=500, algorithm='interval_tree'):
    if not sv_evidence_list:
        return []

    # Separate insertions from other variants for sequence-based pre-filtering
    insertions = [e for e in sv_evidence_list if e.sv_type == 'INS']
    non_insertions = [e for e in sv_evidence_list if e.sv_type != 'INS']
    
    # Pre-filter insertions by sequence similarity (Hamming distance)
    if insertions:
        print(f"\n=== Insertion Sequence Pre-filtering (Hamming distance) ===")
        print(f"Starting with {len(insertions)} insertions")
        insertions = _prefilter_insertions_by_sequence(insertions, max_hamming=0.2)
        print(f"After sequence pre-filtering: {len(insertions)} insertions remain (grouped by similarity)\n")
    
    # Recombine for spatial clustering
    filtered_evidence = insertions + non_insertions
    
    if not filtered_evidence:
        return []
    
    coordinates = np.array([[e.posA, e.posB] for e in filtered_evidence])
    
    if algorithm == 'interval_tree':
        labels = interval_tree_cluster(coordinates, distance_threshold)
    elif algorithm == 'optics':
        from .optics_clustering import optics_cluster
        labels = optics_cluster(coordinates, min_samples=2, max_eps=distance_threshold)
    elif algorithm == 'dbscan':
        labels = DBSCAN.cluster(coordinates, distance_threshold, 2)
    
    clusters = {}
    for i, label in enumerate(labels):
          if label not in clusters:
              clusters[label] = []
          clusters[label].append(filtered_evidence[i])
    
    clustered = []
    for cluster_evidence in clusters.values():
          merged = _merge_evidence_cluster(cluster_evidence)
          clustered.append(merged)
    
    return clustered


def _hamming_distance(seq1, seq2):
    """Calculate normalized Hamming distance between two sequences."""
    if seq1 is None or seq2 is None or len(seq1) == 0 or len(seq2) == 0:
        return 1.0
    
    seq1 = str(seq1).upper()
    seq2 = str(seq2).upper()
    min_len = min(len(seq1), len(seq2))
    max_len = max(len(seq1), len(seq2))
    
    mismatches = sum(1 for i in range(min_len) if seq1[i] != seq2[i])
    length_diff = max_len - min_len
    total_dist = mismatches + length_diff
    
    normalized = total_dist / max_len if max_len > 0 else 0.0
    return normalized


def _prefilter_insertions_by_sequence(insertions, max_hamming=0.2):
    """
    Group insertions by sequence similarity using Hamming distance.
    Prints matching groups for testing and returns grouped insertions.
    """
    if not insertions:
        return insertions
    
    # Filter out insertions with no sequence
    insertions_with_seq = [e for e in insertions if e.inserted_sequence and len(str(e.inserted_sequence)) > 0]
    insertions_without_seq = [e for e in insertions if not e.inserted_sequence or len(str(e.inserted_sequence)) == 0]
    
    if not insertions_with_seq:
        print(f"  No insertions with sequences to filter. Using {len(insertions_without_seq)} sequences-less insertions.")
        return insertions
    
    # Build sequence groups using Hamming distance
    groups = []
    assigned = set()
    
    for i, ins in enumerate(insertions_with_seq):
        if i in assigned:
            continue
        
        group = [i]
        assigned.add(i)
        
        for j in range(i + 1, len(insertions_with_seq)):
            if j in assigned:
                continue
            
            dist = _hamming_distance(insertions_with_seq[i].inserted_sequence, 
                                     insertions_with_seq[j].inserted_sequence)
            if dist <= max_hamming:
                group.append(j)
                assigned.add(j)
        
        groups.append(group)
    
    # Print groups for testing
    for group_id, group_indices in enumerate(groups):
        if len(group_indices) > 1:
            seqs = [str(insertions_with_seq[idx].inserted_sequence)[:20] for idx in group_indices]
            print(f"  Group {group_id}: {len(group_indices)} insertions with similar sequences:")
            for idx, seq_preview in zip(group_indices, seqs):
                print(f"    - {seq_preview}... (support: {len(insertions_with_seq[idx].support_reads)})")
        else:
            seq_preview = str(insertions_with_seq[group_indices[0]].inserted_sequence)[:20]
            print(f"  Group {group_id}: {seq_preview}... (singleton, {len(insertions_with_seq[group_indices[0]].support_reads)} support)")
    
    # Return all insertions (with and without seq) — grouping is done, clustering will use spatial coords
    return insertions_with_seq + insertions_without_seq


def _merge_evidence_cluster(evidence_list):
    pos_a_list = [e.posA for e in evidence_list]
    pos_b_list = [e.posB for e in evidence_list]
    
    pos_a_median = sorted(pos_a_list)[len(pos_a_list) // 2]
    pos_b_median = sorted(pos_b_list)[len(pos_b_list) // 2]
    
    all_reads = []
    for evidence in evidence_list:
        all_reads.extend(evidence.support_reads)
    
    merged = SVEvidence(evidence_list[0].chrA, pos_a_median, evidence_list[0].chrB, pos_b_median, evidence_list[0].sv_type, all_reads)
    
    merged.ci_a_lower = pos_a_median - min(pos_a_list)
    merged.ci_a_upper = max(pos_a_list) - pos_a_median
    merged.ci_b_lower = pos_b_median - min(pos_b_list)
    merged.ci_b_upper = max(pos_b_list) - pos_b_median
    
    # For insertions, preserve a representative sequence
    if evidence_list[0].sv_type == 'INS' and evidence_list[0].inserted_sequence:
        merged.inserted_sequence = evidence_list[0].inserted_sequence
    
    return merged

def read_bam_file(bam_file, sample_name, min_sv_size=50, min_mapq=20, algorithm='interval_tree'):

    split_read_evidence = extract_sv_from_split_reads(bam_file, min_mapq, min_sv_size)
    cigar_evidence = extract_sv_from_cigar(bam_file, min_sv_size)
    all_evidence = split_read_evidence + cigar_evidence

    print(f"Found {len(all_evidence)} raw SV evidence")

    clustered_evidence = cluster_sv_evidence(all_evidence, algorithm=algorithm)

    print(f"Clustered into {len(clustered_evidence)} consensus SVs using {algorithm}")
    
    # Print insertion sequence frequency for testing
    _print_insertion_sequence_frequency(clustered_evidence)
    
    vcf_format_variants = bam_to_vcf_format(clustered_evidence, sample_name)
    
    return vcf_format_variants


def _print_insertion_sequence_frequency(clustered_evidence):
    """
    Print frequency of insertion sequences for testing/debugging.
    Shows how many times each unique insertion sequence was found.
    """
    insertions = [e for e in clustered_evidence if e.sv_type == 'INS']
    if not insertions:
        return
    
    seq_counts = {}
    for ins in insertions:
        if ins.inserted_sequence and len(str(ins.inserted_sequence)) > 0:
            seq = str(ins.inserted_sequence)
            seq_counts[seq] = seq_counts.get(seq, 0) + 1
    
    if seq_counts:
        print(f"\n=== Insertion Sequence Frequency Summary ===")
        sorted_seqs = sorted(seq_counts.items(), key=lambda x: x[1], reverse=True)
        for seq, count in sorted_seqs:
            seq_display = seq if len(seq) <= 30 else seq[:27] + "..."
            print(f"  Sequence: {seq_display}")
            print(f"    Frequency: {count} (support reads: {count})")
        print(f"Total unique insertion sequences: {len(seq_counts)}\n")
    else:
        print(f"No insertion sequences captured (insertions may lack sequence data).\n")