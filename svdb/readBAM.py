import pysam
from collections import defaultdict

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

    for read_name, alignments in read_alignments.items():
        if len(alignments) < 2:
            continue
        alignments.sort(key=lambda x: (x.reference_name, x.reference_start))
        for i in range(len(alignments) - 1):
            primary = alignments[i]
            supplementary = alignments[i + 1]
            sv_type, evidence = _analyze_alignment_pair(primary, supplementary, min_sv_size)
            if sv_type and evidence:
                sv_evidence.append(evidence)
    samfile.close()
    return sv_evidence


def _analyze_alignment_pair(aln1, aln2, min_sv_size):
    chr1 = aln1.reference_name
    pos1 = aln1.reference_start
    chr2 = aln2.reference_name
    pos2 = aln2.reference_start

    if chr1 == chr2:
        distance = abs(pos2 - pos1)
        if distance < min_sv_size:
            return None, None
    same_strand = aln1.is_reverse == aln2.is_reverse
    if same_strand:
        if distance > min_sv_size:
            sv_type = "DEL" 
            evidence = SVEvidence(chr1, min(pos1, pos2), chr2, max(pos1, pos2), sv_type, [aln1.query_name])
            return sv_type, evidence
    else:
        sv_type = "BND"  
        if chr1 > chr2:
            chr1, chr2 = chr2, chr1
            pos1, pos2 = pos2, pos1
        evidence = SVEvidence(chr1, pos1, chr2, pos2, sv_type, [aln1.query_name])
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
        for op, length in read.cigartuples:
            if op == 2 and length >= min_indel_size:  
                evidence = SVEvidence(read.reference_name, ref_pos, read.reference_name, ref_pos + length, "DEL", [read.query_name])
                sv_evidence.append(evidence)
                ref_pos += length
            elif op == 1 and length >= min_indel_size:  
                evidence = SVEvidence(read.reference_name, ref_pos, read.reference_name, ref_pos, "INS", [read.query_name])
                sv_evidence.append(evidence) 
            elif op in [0,2,3,7,8]:
                ref_pos += length
    samfile.close()
    return sv_evidence

def bam_to_vcf_format(sv_evidence_list, sample_name):
    vcf_like_variants = []
    for evidence in sv_evidence_list:
        INFO = {}
        FORMAT = {'GT': ['0/1']} #Assuming heterozygous otherwise change to '1/1'
        if evidence.chrA == evidence.chrB:
            INFO['SVLEN'] = str(abs(evidence.posB - evidence.posA))
            INFO['END'] = str(evidence.posB)

        INFO['SVTYPE'] = evidence.sv_type
        INFO['support'] = str(len(evidence.support_reads))

        vcf_like_variants.append((evidence.chrA, evidence.posA, evidence.chrB, evidence.posB, evidence.sv_type, INFO, FORMAT, sample_name))
    return vcf_like_variants


def cluster_sv_evidence(sv_evidence_list, distance_threshold=500):
    if not sv_evidence_list:
        return []
    
    groups = defaultdict(list)
    for evidence in sv_evidence_list:
        key = (evidence.chrA, evidence.chrB, evidence.sv_type)
        groups[key].append(evidence)
    
    clustered = []
    
    for key, evidence_group in groups.items():
        evidence_group.sort(key=lambda x: (x.posA, x.posB))
        
        current_cluster = [evidence_group[0]]
        
        for evidence in evidence_group[1:]:
            last = current_cluster[-1]
        
            if (abs(evidence.posA - last.posA) <= distance_threshold and
                abs(evidence.posB - last.posB) <= distance_threshold):
                current_cluster.append(evidence)
            else:
                merged = _merge_evidence_cluster(current_cluster)
                clustered.append(merged)
                current_cluster = [evidence]
        
        if current_cluster:
            merged = _merge_evidence_cluster(current_cluster)
            clustered.append(merged)
    
    return clustered


def _merge_evidence_cluster(evidence_list):

    pos_a_list = [e.posA for e in evidence_list]
    pos_b_list = [e.posB for e in evidence_list]
    
    pos_a_median = sorted(pos_a_list)[len(pos_a_list) // 2]
    pos_b_median = sorted(pos_b_list)[len(pos_b_list) // 2]
    
    all_reads = []
    for evidence in evidence_list:
        all_reads.extend(evidence.support_reads)
    
    merged = SVEvidence( evidence_list[0].chrA, pos_a_median, evidence_list[0].chrB, pos_b_median, evidence_list[0].sv_type, all_reads)
    
    merged.ci_a_lower = pos_a_median - min(pos_a_list)
    merged.ci_a_upper = max(pos_a_list) - pos_a_median
    merged.ci_b_lower = pos_b_median - min(pos_b_list)
    merged.ci_b_upper = max(pos_b_list) - pos_b_median
    
    return merged

def read_bam_file(bam_file, sample_name, min_sv_size=50, min_mapq=20):
    print(f"Extracting SV evidence from BAM file: {bam_file}")

    split_read_evidence = extract_sv_from_split_reads(
        bam_file, min_mapq, min_sv_size
    )

    cigar_evidence = extract_sv_from_cigar(bam_file, min_sv_size)

    all_evidence = split_read_evidence + cigar_evidence

    print(f"Found {len(all_evidence)} raw SV evidence")

    clustered_evidence = cluster_sv_evidence(all_evidence)

    print(f"Clustered into {len(clustered_evidence)} consensus SVs")
    
    vcf_format_variants = bam_to_vcf_format(clustered_evidence, sample_name)
    
    return vcf_format_variants
