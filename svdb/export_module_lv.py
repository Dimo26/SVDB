from __future__ import absolute_import
import sys
import numpy as np
from . import database, overlap_module_lv


class DBSCAN:
    def x_coordinate_clustering(data, epsilon, m):
        clusters = np.zeros(len(data)) - 1
        cluster_id = -1
        cluster = False

        for i in range(len(data) - m + 1):
            current = data[i, :]
            points = data[i + 1:i + m, :]
            distances = [abs(point[0] - current[0]) for point in points]
            if max(distances) < epsilon:
                if cluster:  
                    clusters[i + m - 1] = cluster_id
                else:  
                    cluster_id += 1
                    cluster = True
                    for j in range(i, i + m):
                        clusters[j] = cluster_id
            else:
                cluster = False
        return clusters, cluster_id

    def y_coordinate_clustering(data, epsilon, m, cluster_id, clusters):

        cluster_id_list = set(clusters)
        for cluster in cluster_id_list:
            if cluster == -1:
                continue
            class_member_mask = (clusters == cluster)
            indexes = np.where(class_member_mask)[0]
            signals = data[class_member_mask]

            y_coordinates = [[signal[1], indexes[i]] for i, signal in enumerate(signals)]
            y_coordinates.sort(key=lambda x: x[0])

            sub_clusters = np.zeros(len(indexes)) - 1

            active_cluster = False
            sub_cluster_id = 0
            y_coordinates = np.array(y_coordinates)
            for i in range(len(y_coordinates) - m + 1):
                current = y_coordinates[i, :]
                distances = [abs(pos[0] - current[0]) for pos in y_coordinates[i + 1:i + m, :]]

                if max(distances) < epsilon:
                    if active_cluster:
                        sub_clusters[i + m - 1] = sub_cluster_id
                    else:
                        sub_cluster_id += 1
                        active_cluster = True
                        for j in range(i, i + m):
                            sub_clusters[j] = sub_cluster_id
                else:
                    active_cluster = False

            for i in range(len(sub_clusters)):
                if sub_clusters[i] == 1:
                    clusters[y_coordinates[i][1]] = cluster
                elif sub_clusters[i] > -1:
                    clusters[y_coordinates[i][1]] = sub_clusters[i] + cluster_id - 1
                elif sub_clusters[i] == -1:
                    clusters[y_coordinates[i][1]] = -1
            if sub_cluster_id > 1:
                cluster_id += sub_cluster_id - 1
        return clusters, cluster_id

    def cluster(data, epsilon, m):
        clusters, cluster_id = DBSCAN.x_coordinate_clustering(data, epsilon, m)
        clusters, cluster_id = DBSCAN.y_coordinate_clustering(data, epsilon, m, cluster_id, clusters)
        return clusters


def _levenshtein_distance(seq1, seq2):
    """Calculate normalized Levenshtein distance between two sequences."""
    if seq1 is None or seq2 is None or len(seq1) == 0 or len(seq2) == 0:
        return 1.0
    
    seq1 = str(seq1).upper()
    seq2 = str(seq2).upper()
    
    len1 = len(seq1)
    len2 = len(seq2)
    max_len = max(len1, len2)
    
    # Create DP table
    dp = [[0] * (len2 + 1) for _ in range(len1 + 1)]
    
    # Initialize
    for i in range(len1 + 1):
        dp[i][0] = i
    for j in range(len2 + 1):
        dp[0][j] = j
    
    # Fill DP table
    for i in range(1, len1 + 1):
        for j in range(1, len2 + 1):
            if seq1[i-1] == seq2[j-1]:
                dp[i][j] = dp[i-1][j-1]
            else:
                dp[i][j] = 1 + min(dp[i-1][j], dp[i][j-1], dp[i-1][j-1])
    
    normalized = dp[len1][len2] / max_len if max_len > 0 else 0.0
    return normalized


def apply_levenshtein_to_insertions(labels, variant_dict, max_edit=0.2):
    """
    Re-cluster insertions within spatial clusters using Levenshtein distance.
    
    This is the Levenshtein equivalent of apply_hamming_to_insertions.
    Maintains same logic: spatial clusters first, then sequence-based re-clustering.
    """
    if labels is None or variant_dict is None:
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
        
        # Separate insertions with sequences from others
        ins_with_seq = []
        other_indices = []
        
        for idx in indices:
            var_info = variant_dict.get(idx, {})
            if var_info.get('type') == 'INS' and var_info.get('sequence'):
                ins_with_seq.append(idx)
            else:
                other_indices.append(idx)
        
        # Non-insertion variants keep their cluster
        if other_indices:
            new_labels[other_indices] = next_cluster_id
            next_cluster_id += 1
        
        # Re-cluster insertions by sequence similarity using Levenshtein
        if ins_with_seq:
            assigned = set()
            for i in range(len(ins_with_seq)):
                if i in assigned:
                    continue
                
                group = [i]
                assigned.add(i)
                idx_i = ins_with_seq[i]
                seq_i = variant_dict[idx_i]['sequence']
                
                for j in range(i + 1, len(ins_with_seq)):
                    if j in assigned:
                        continue
                    
                    idx_j = ins_with_seq[j]
                    seq_j = variant_dict[idx_j]['sequence']
                    
                    dist = _levenshtein_distance(seq_i, seq_j)
                    if dist <= max_edit:
                        group.append(j)
                        assigned.add(j)
                
                # Assign new cluster ID for this sequence group
                for g in group:
                    new_labels[ins_with_seq[g]] = next_cluster_id
                next_cluster_id += 1
    
    return new_labels


def fetch_index_variant(db, index):
    A = 'SELECT posA, ci_A_lower, ci_A_upper, posB, ci_B_lower, ci_B_upper, sample, sequence FROM SVDB WHERE idx IN ({}) '.format(
        ", ".join([str(idx) for idx in index]))
    hits = db.query(A)
    variant = {}
    coordinates = []
    for i, hit in enumerate(hits):
        variant[i] = {}
        variant[i]["posA"] = int(hit[0])
        variant[i]["ci_A_start"] = int(hit[1])
        variant[i]["ci_A_end"] = int(hit[2])
        variant[i]["posB"] = int(hit[3])
        variant[i]["ci_B_start"] = int(hit[4])
        variant[i]["ci_B_end"] = int(hit[5])
        variant[i]["sample_id"] = hit[6]
        variant[i]["sequence"] = hit[7] if len(hit) > 7 else ""
        coordinates.append([i, int(hit[0]), int(hit[3])])
    return variant, np.array(coordinates)


def fetch_cluster_variant(db, index):
    query = 'SELECT posA, posB, sample, idx, sequence FROM SVDB WHERE idx IN ({}) '.format(
            ", ".join([str(idx) for idx in index]))
    hits = db.query(query)

    variant_dict = {}
    for hit in hits:
        variant_dict[int(hit[3])] = {}
        variant_dict[int(hit[3])]["posA"] = int(hit[0])
        variant_dict[int(hit[3])]["posB"] = int(hit[1])
        variant_dict[int(hit[3])]["sample_id"] = hit[2]
        variant_dict[int(hit[3])]["sequence"] = hit[4] if len(hit) > 4 else ""
    return variant_dict

def db_header(args):
    headerString = '##fileformat=VCFv4.1\n'
    headerString += '##source=SVDB\n'
    headerString += '##ALT=<ID=DEL,Description="Deletion">\n'
    headerString += '##ALT=<ID=DUP,Description="Duplication">\n'
    headerString += '##ALT=<ID=INV,Description="Inversion">\n'
    headerString += '##ALT=<ID=INS,Description="Insertion">\n'
    headerString += '##ALT=<ID=BND,Description="Break end">\n'
    headerString += '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">\n'
    headerString += '##INFO=<ID=END,Number=1,Type=String,Description="End of an intra-chromosomal variant">\n'
    headerString += '##INFO=<ID=OCC,Number=1,Type=Integer,Description="The number of occurences of the event in the database">\n'
    headerString += '##INFO=<ID=NSAMPLES,Number=1,Type=Integer,Description="the number of samples within the database">\n'
    headerString += '##INFO=<ID=VARIANTS,Number=1,Type=Integer,Description="a| separated list of the positions of the clustered variants">\n'
    headerString += '##INFO=<ID=FRQ,Number=1,Type=Float,Description="the frequency of the variant">\n'
    headerString += '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">\n'
    headerString += '##INFO=<ID=CIPOS,Number=2,Type=Integer,Description="Confidence interval around POS">\n'
    headerString += '##INFO=<ID=CIEND,Number=2,Type=Integer,Description="Confidence interval around END">\n'
    headerString += '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
    headerString += '##SVDB_version={} cmd=\"{}\"'.format(args.version, " ".join(sys.argv))
    return headerString

def vcf_line(cluster, id_tag, sample_IDs):
    info_field = "SVTYPE={};".format(cluster[0]["type"])
    vcf_line = []
    vcf_line.append(cluster[0]["chrA"])
    vcf_line.append(str(cluster[0]["posA"]))
    vcf_line.append(id_tag)
    vcf_line.append("N")
    if cluster[0]["chrA"] == cluster[0]["chrB"] and cluster[0]["type"] != "BND":
        vcf_line.append("<" + cluster[0]["type"] + ">")
        info_field += "END={};SVLEN={};".format(cluster[0]["posB"], abs(cluster[0]["posA"] - cluster[0]["posB"]))
    else:
        vcf_line.append("N[{}:{}[".format(cluster[0]["chrB"], cluster[0]["posB"]))

    sample_set = set([])
    CIPOS = []
    CIEND = []
    for variant in cluster[1]:
        CIPOS.append(cluster[1][variant]["posA"])
        CIEND.append(cluster[1][variant]["posB"])
        sample_set.add(cluster[1][variant]["sample_id"])

    CIPOS_start = -abs(cluster[0]["posA"] - min(CIPOS))
    CIPOS_end = abs(cluster[0]["posA"] - max(CIPOS))

    CIEND_start = -abs(cluster[0]["posB"] - min(CIEND))
    CIEND_end = abs(cluster[0]["posB"] - max(CIEND))

    info_field += "NSAMPLES={};OCC={};FRQ={};CIPOS={},{};CIEND={},{};".format(len(sample_IDs), len(
        sample_set), round(len(sample_set) / float(len(sample_IDs)), 4), CIPOS_start, CIPOS_end, CIEND_start, CIEND_end)
    variant_field = "VARIANTS="
    for variant in cluster[1]:
        variant_field += "|{}:{}:{}".format(cluster[1][variant]["sample_id"], cluster[1][variant]["posA"], cluster[1][variant]["posB"])
    info_field += variant_field
    vcf_line.append(".")
    vcf_line.append("PASS")
    vcf_line.append(info_field)
    zygosity_list = {}
    for sample in sample_IDs:
        zygosity_list[sample] = "0/0"

    for variant in cluster[1]:
        zygosity_list[cluster[1][variant]["sample_id"]] = "./1"
    format_cols = []
    for sample in sample_IDs:
        format_cols.append(zygosity_list[sample])
    vcf_line.append("GT")
    vcf_line.append("\t".join(format_cols))
    return "\t".join(vcf_line)

def expand_chain(chain, coordinates, chrA, chrB, distance, overlap):
    from .interval_tree_overlap import IntervalTree
    tree = IntervalTree()
    for i, coord in enumerate(coordinates):
         tree.add(coord[1], coord[2], index=i)
    tree.build()
    
    chain_data = {}

    for i, idx in enumerate(chain):
          chain_data[i] = []
          variant = chain[idx]

          overlaps = tree.query(variant["posA"] - distance, variant["posB"] + distance)

          for interval in overlaps:
              candidate_idx = interval.index
              var = chain[candidate_idx]

              similar = False
              match = False
              if chrA != chrB:
                 similar = True
                 match = True
              else:
                  similar, match = overlap_module_lv.isSameVariation(variant["posA"], variant["posB"], var["posA"], var["posB"], overlap, distance)
              if match:
                  chain_data[i].append(candidate_idx)
         
          chain_data[i] = np.array(chain_data[i])
    
    return chain_data

def cluster_variants(variant_dictionary, similarity_matrix):
    cluster_sizes = [[i, len(similarity_matrix[i])] for i in range(len(variant_dictionary))]

    clusters = []
    for i, _ in sorted(cluster_sizes, key=lambda x: (x[1]), reverse=True):
        if similarity_matrix[i][0] == -1:
            continue

        cluster_dictionary = {}
        for var in similarity_matrix[i]:
            similarity_matrix[var][0] = -1
            cluster_dictionary[var] = variant_dictionary[var]
        variant = variant_dictionary[i]

        clusters.append([variant, cluster_dictionary])
    return clusters

def fetch_variants(variant, chrA, chrB, db):
    chr_db = {}
    chr_db[variant] = {}

    hits = db.query('SELECT posA,posB,sample,idx,var,sequence FROM SVDB WHERE var == \'{}\'AND chrA == \'{}\' AND chrB == \'{}\''.format(
        variant, chrA, chrB))
    if not hits:
        return False

    x = [v[0] for v in hits]
    y = [v[1] for v in hits]

    chr_db[variant]["coordinates"] = np.column_stack((x, y))
    chr_db[variant]["var_info"] = np.array([v[2] for v in hits])
    chr_db[variant]["index"] = np.array([v[3] for v in hits])
    chr_db[variant]["sequences"] = [v[5] if len(v) > 5 else "" for v in hits]
    return chr_db


def overlap_cluster(db, indexes, variant, chrA, chrB, sample_IDs, args, f, i):
    variant_dictionary, coordinates = fetch_index_variant(db, indexes)
    if "INS" in variant:
        similarity_matrix = expand_chain(
           variant_dictionary, coordinates, chrA, chrB, args.ins_distance, -1)
    else:
        similarity_matrix = expand_chain(
           variant_dictionary, coordinates, chrA, chrB, args.bnd_distance, args.overlap)

    clusters = cluster_variants(variant_dictionary, similarity_matrix)
    for clustered_variants in clusters:
        clustered_variants[0]["type"] = variant
        clustered_variants[0]["chrA"] = chrA
        clustered_variants[0]["chrB"] = chrB
        f.write(vcf_line(clustered_variants, "cluster_{}".format(i), sample_IDs) + "\n")
    return i + len(clusters)


def svdb_cluster_main(chrA, chrB, variant, sample_IDs, args, db, i, algorithm='DBSCAN'):
    f = open(args.prefix + ".vcf", 'a')
    chr_db = fetch_variants(variant, chrA, chrB, db)
    if not chr_db:
        f.close()
        return i
    
    # Step 1: Spatial clustering (DBSCAN, OPTICS, or INTERVAL_TREE)
    if algorithm == 'OPTICS':
        from .optics_clustering import optics_cluster
        dbscan = optics_cluster(
            chr_db[variant]["coordinates"],
            min_samples=2,
            max_eps=args.bnd_distance
        )
    elif algorithm == 'INTERVAL_TREE':
        from .interval_tree_overlap import interval_tree_cluster
        dbscan = interval_tree_cluster(
            chr_db[variant]["coordinates"],
            max_distance=args.bnd_distance
        )
    elif args.DBSCAN or algorithm == 'DBSCAN':
        if args.DBSCAN:
            dbscan = DBSCAN.cluster(chr_db[variant]["coordinates"], args.epsilon, args.min_pts)
        elif "INS" in variant:
            dbscan = DBSCAN.cluster(chr_db[variant]["coordinates"], args.ins_distance, 2)        
        else:
            dbscan = DBSCAN.cluster(chr_db[variant]["coordinates"], args.bnd_distance, 2)
    else:
        # Default to overlap-based clustering
        dbscan = None
    
    # Step 2: Apply Levenshtein distance for insertions (if requested)
    if getattr(args, 'use_levenshtein', False) and "INS" in variant:
        # Build variant dictionary for Levenshtein distance
        variant_dict = {}
        for idx, seq in zip(chr_db[variant]["index"], chr_db[variant]["sequences"]):
            variant_dict[idx] = {
                'type': variant,
                'sequence': seq
            }
        
        # Apply Levenshtein-based re-clustering
        max_edit = getattr(args, 'max_edit', 0.2)
        dbscan = apply_levenshtein_to_insertions(dbscan, variant_dict, max_edit)
        print(f"  Applied Levenshtein distance clustering for {variant} (max_edit={max_edit})")
    
    if dbscan is None:
        # Fall back to overlap-based clustering
        unique_xy = chr_db[variant]["coordinates"]
        unique_index = chr_db[variant]["index"]
        i = overlap_cluster(db, unique_index, variant, chrA, chrB, sample_IDs, args, f, i)
        f.close()
        return i
    
    unique_labels = set(dbscan)
    unique_xy = chr_db[variant]["coordinates"][dbscan == -1]
    unique_index = chr_db[variant]["index"][dbscan == -1]
    
    # Handle noise points
    for xy, indexes in zip(unique_xy, unique_index):
        variant_dictionary = fetch_cluster_variant(db, [indexes])
        representing_var = {}
        representing_var["type"] = variant
        representing_var["chrA"] = chrA
        representing_var["chrB"] = chrB
        representing_var["posA"] = xy[0]
        representing_var["ci_A_start"] = xy[0]
        representing_var["ci_A_end"] = xy[0]
        representing_var["posB"] = xy[1]
        representing_var["ci_B_start"] = xy[1]
        representing_var["ci_B_end"] = xy[1]

        cluster = [representing_var, variant_dictionary]
        f.write(vcf_line(cluster, "cluster_{}".format(i), sample_IDs) + "\n")
        i += 1
    del unique_xy
    del unique_index

    # Handle clustered points
    for unique_label in unique_labels:
        if unique_label == -1:
            continue
        class_member_mask = (dbscan == unique_label)
        xy = chr_db[variant]["coordinates"][class_member_mask]
        indexes = chr_db[variant]["index"][class_member_mask]

        if args.DBSCAN or algorithm in ['DBSCAN', 'OPTICS', 'INTERVAL_TREE']:
            avg_point = np.array([np.mean(xy[:, 0]), np.mean(xy[:, 1])])

            variant_dictionary = fetch_cluster_variant(db, indexes)

            representing_var = {}
            representing_var["type"] = variant
            representing_var["chrA"] = chrA
            representing_var["chrB"] = chrB
            representing_var["posA"] = int(avg_point[0])
            representing_var["ci_A_start"] = np.amin(xy[:, 0])
            representing_var["ci_A_end"] = np.amax(xy[:, 0])
            representing_var["posB"] = int(avg_point[1])
            representing_var["ci_B_start"] = np.amin(xy[:, 1])
            representing_var["ci_B_end"] = np.amax(xy[:, 1])

            cluster = [representing_var, variant_dictionary]
            f.write(vcf_line(cluster, "cluster_{}".format(i), sample_IDs) + "\n")
            i += 1

        else:
            i = overlap_cluster(db, indexes, variant, chrA,
                                chrB, sample_IDs, args, f, i)

    f.close()
    return i

def export(args, sample_IDs):
    db = database.DB(args.db, memory=args.memory)

    chrA_list = []
    for chrA in db.query('SELECT DISTINCT chrA FROM SVDB'):
        chrA_list.append(chrA[0])

    chrB_list = []
    for chrB in db.query('SELECT DISTINCT chrB FROM SVDB'):
        chrB_list.append(chrB[0])

    var_list = []
    for variant in db.query('SELECT DISTINCT var FROM SVDB'):
        var_list.append(variant[0])

    i = 0
    algorithm = getattr(args, 'algorithm', 'DBSCAN')
    
    for chrA in chrA_list:
        for chrB in chrB_list:
            for variant in var_list:
                i = svdb_cluster_main(chrA, chrB, variant, sample_IDs, args, db, i, algorithm=algorithm)

def main(args):
    sample_IDs = []
    if not args.prefix:
        args.prefix = args.db.replace(".db", "")

    db = database.DB(args.db)

    for sample in db.query('SELECT DISTINCT sample FROM SVDB'):
        sample_IDs.append(sample[0])

    with open(args.prefix + ".vcf", 'w') as f:
        f.write(db_header(args) + "\n")
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{}\n".format("\t".join(sample_IDs)))
    export(args, sample_IDs)
