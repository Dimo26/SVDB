import re
import numpy as np

class sveEvidence:
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


def cluster_sv_evidence(sv_evidence_list, distance_threshold=500, algorithm='interval_tree'):
    if not sv_evidence_list:
        return []

    coordinates = np.array([[e.posA, e.posB] for e in sv_evidence_list])
    
    if algorithm == 'interval_tree':
        from interval_tree_overlap import interval_tree_cluster
        labels = interval_tree_cluster(coordinates, distance_threshold)
    elif algorithm == 'optics':
        from optics_clustering import optics_cluster
        labels = optics_cluster(coordinates, min_samples=2, max_eps=distance_threshold)
    elif algorithm == 'dbscan':
        from export_module import DBSCAN
        labels = DBSCAN.cluster(coordinates, distance_threshold, 2)
    else:
        from interval_tree_overlap import interval_tree_cluster
        labels = interval_tree_cluster(coordinates, distance_threshold)

    clusters = {}
    for i, label in enumerate(labels):
          if label not in clusters:
              clusters[label] = []
          clusters[label].append(sv_evidence_list[i])
    
    clustered = []
    for cluster_evidence in clusters.values():
          merged = _merge_evidence_cluster(cluster_evidence)
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
    
    merged = sveEvidence(evidence_list[0].chrA, pos_a_median, evidence_list[0].chrB, pos_b_median, evidence_list[0].sv_type, all_reads)
    
    merged.ci_a_lower = pos_a_median - min(pos_a_list)
    merged.ci_a_upper = max(pos_a_list) - pos_a_median
    merged.ci_b_lower = pos_b_median - min(pos_b_list)
    merged.ci_b_upper = max(pos_b_list) - pos_b_median
    
    return merged


def readVCFLine(line):
    if line[0].startswith("#"):
        return None

    variation = line.strip().split("\t")
    event_type = ""
    chrA = variation[0].replace("chr", "").replace("Chr", "").replace("CHR", "")
    posA = int(variation[1])
    posB = 0

    description = {}
    INFO = variation[7].split(";")
    for tag in INFO:
        tag = tag.split("=")
        if(len(tag) > 1):
            description[tag[0]] = tag[1]

    format = {}
    format_keys = {}
    if len(variation) > 8:
        format_string = variation[8].split(":")

        i = 0
        for key in format_string:
            format_keys[i] = key
            format[key] = []
            i += 1
        format_fields = variation[9:]
        for sample in format_fields:
            i = 0
            format_string = sample.split(":")
            for i in range(0, len(format_keys)):
                format[format_keys[i]].append(format_string[i])
                i += 1

    # Delly translocations
    if "TRA" in variation[4]:
        event_type = "BND"
        chrB = description["CHR2"]
        posB = int(description["END"])
        if chrA > chrB:
            chrT = chrA
            chrA = chrB
            chrB = chrT
            posB = posA

    # intrachromosomal variant
    elif "]" not in variation[4] and "[" not in variation[4]:

        chrB = chrA
        posB = posA

        if "END" in description:
            posB = int(description["END"])
        elif "SVLEN" in description:
            posB = posA + abs(int(description["SVLEN"]))

        if posB < posA:
            tmp = posB
            posB = posA
            posA = tmp

        nucleotides=set(["A","T","C","G","N"])
        #SVtype is given in alt column
        if "<" in variation[4] and ">" in variation[4]:
            event_type = variation[4].strip("<").rstrip(">")
            if "DUP" in event_type:
                event_type = "DUP"

        #SVtype present in INFO
        elif "SVTYPE" in description:
            event_type = description["SVTYPE"]

        #no SVTYPE, actual sequence is written in alt
        elif set(list(variation[4])).union(nucleotides) == nucleotides:
            if len(variation[4]) > len(variation[3]):
                event_type="INS"
            else:
                event_type="DEL"
                posB=posA+len(variation[3])-1

        #treat the insertion as single points
        if "INS" in event_type:
            posA=int(variation[1])
            posB=int(variation[1])
            insertion_seq = None
            
            if "<INS>" not in variation[4]:

                  insertion_seq = variation[4]
            else:
                  if "SEQ" in description:
                      insertion_seq = description["SEQ"]
            
            if insertion_seq is not None:
                  description["INSSEQ"] = insertion_seq
    else:
        B = variation[4]
        combinations={"[]":"]", "[[":"[", "]]":"]", "][":"["}

        for c in combinations:
            if B.startswith(c):
                B=B.replace(c,combinations[c])

        B = re.split("[],[]", B)
        chr_and_pos = B[1]
        chrB = ":".join(chr_and_pos.split(":")[:-1]).replace("chr", "").replace("Chr", "").replace("CHR", "")
        posB = int(chr_and_pos.split(":")[-1])
        if chrA > chrB:
            chrT = chrA
            chrA = chrB
            chrB = chrT
            posA, posB = posB, posA

        if chrA == chrB: #intrachromosomal
            if posB < posA:
                posA, posB = posB, posA

        event_type = "BND"

    return chrA, posA, chrB, posB, event_type, description, format