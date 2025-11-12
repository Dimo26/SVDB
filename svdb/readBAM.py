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
            if op == 2 and length >= min_indel_size:  # DELETION
                evidence = SVEvidence(read.reference_name, ref_pos, read.reference_name, ref_pos + length, "DEL", [read.query_name])
                sv_evidence.append(evidence)
                ref_pos += length
            elif op == 1 and length >= min_indel_size:  # INSERTION
                evidence = SVEvidence(read.reference_name, ref_pos, read.reference_name, ref_pos, "INS", [read.query_name])
                sv_evidence.append(evidence) 
            elif op in [0,2,3,7,8]:
                ref_pos += length
    samfile.close()
    return sv_evidence
