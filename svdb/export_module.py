from __future__ import absolute_import
import sys
import numpy as np
from . import database
from . import overlap_module
from . import DBSCAN

try:
    from interval_tree_overlap import IntervalTree, interval_tree_cluster
    HAS_INTERVAL_TREE = True
except ImportError:
    HAS_INTERVAL_TREE = False

try:
    from optics_clustering import optics_cluster
    HAS_OPTICS = True
except ImportError:
    HAS_OPTICS = False


def fetch_index_variant(db, index):
    A = 'SELECT posA, ci_A_lower, ci_A_upper, posB, ci_B_lower, ci_B_upper, sample FROM SVDB WHERE idx IN ({}) '.format(
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
        coordinates.append([i, int(hit[0]), int(hit[3])])
    return variant, np.array(coordinates)


def fetch_cluster_variant(db, index):
    query = 'SELECT posA, posB, sample, idx FROM SVDB WHERE idx IN ({}) '.format(
            ", ".join([str(idx) for idx in index]))
    hits = db.query(query)

    variant_dict = {}
    for hit in hits:
        variant_dict[int(hit[3])] = {}
        variant_dict[int(hit[3])]["posA"] = int(hit[0])
        variant_dict[int(hit[3])]["posB"] = int(hit[1])
        variant_dict[int(hit[3])]["sample_id"] = hit[2]
    return variant_dict


def fetch_variant_sequences(db, index_list):
    """
    Fetch sequence information for insertion variants.
    Returns: dict mapping idx -> sequence string
    """
    query = 'SELECT idx, sequence FROM SVDB WHERE idx IN ({})'.format(
        ", ".join([str(idx) for idx in index_list]))
    hits = db.query(query)
    
    sequences = {}
    for hit in hits:
        sequences[int(hit[0])] = hit[1] if hit[1] else ""
    
    return sequences


def cluster_insertions_by_sequence(variant_dictionary, coordinates, db, max_hamming_distance=0.2):
    variant_indices = list(variant_dictionary.keys())
    sequences = fetch_variant_sequences(db, variant_indices)
    
    similarity_matrix = {}
    

    for i in variant_indices:
        similar_variants = [] 
        
        seq_i = sequences.get(i, "")

        if not seq_i:
            similarity_matrix[i] = np.array([i])  
            continue
        
        for j in variant_indices:
            seq_j = sequences.get(j, "")

            if not seq_j:
                continue
            

            raw_dist, normalized_dist = overlap_module.hamming_distance(seq_i, seq_j)
            

            if normalized_dist <= max_hamming_distance:
                similar_variants.append(j)

        if not similar_variants:
            similar_variants = [i]
        
        similarity_matrix[i] = np.array(similar_variants)
    
    return similarity_matrix


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

    if HAS_INTERVAL_TREE:
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
                    similar, match = overlap_module.isSameVariation(
                        variant["posA"], variant["posB"], 
                        var["posA"], var["posB"], 
                        overlap, distance)
                if match:
                    chain_data[i].append(candidate_idx)
            
            chain_data[i] = np.array(chain_data[i])
        
        return chain_data
    else:
        # Original O(n²) method (fallback)
        chain_data = {}
        for i, idx in enumerate(chain):
            chain_data[i] = []
            variant = chain[idx]

            rows = coordinates[(distance >= abs(coordinates[:, 1] - variant["posA"]))
                               & (distance >= abs(coordinates[:, 2] - variant["posB"]))]

            candidates = rows[:, 0]
            for candidate in candidates:
                var = chain[candidate]
                similar = False
                match = False
                if chrA != chrB:
                    similar = True
                    match = True
                else:
                    similar, match = overlap_module.isSameVariation(
                        variant["posA"], variant["posB"], 
                        var["posA"], var["posB"], 
                        overlap, distance)
                if match:
                    chain_data[i].append(candidate)

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

    hits = db.query('SELECT posA,posB,sample,idx,var FROM SVDB WHERE var == \'{}\'AND chrA == \'{}\' AND chrB == \'{}\''.format(
        variant, chrA, chrB))
    if not hits:
        return False

    x = [v[0] for v in hits]
    y = [v[1] for v in hits]

    chr_db[variant]["coordinates"] = np.column_stack((x, y))
    chr_db[variant]["var_info"] = np.array([v[2] for v in hits])
    chr_db[variant]["index"] = np.array([v[3] for v in hits])
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


def svdb_cluster_main(chrA, chrB, variant, sample_IDs, args, db, i):
  
    f = open(args.prefix + ".vcf", 'a')
    chr_db = fetch_variants(variant, chrA, chrB, db)
    if not chr_db:
        f.close()
        return i

    algorithm = getattr(args, 'algorithm', 'DBSCAN')
    
    cluster_labels = None
    
    if "INS" in variant:
        variant_dictionary, coordinates = fetch_index_variant(db, chr_db[variant]["index"])

        similarity_matrix = cluster_insertions_by_sequence(
            variant_dictionary, 
            coordinates, 
            db,
            max_hamming_distance=0.2  # 80% similarity threshold
        )
        
        dbscan = np.zeros(len(chr_db[variant]["coordinates"])) - 1
        cluster_id = 0
        processed = set()
        
        for idx in variant_dictionary.keys():
            if idx in processed:
                continue
            
            similar_indices = similarity_matrix[idx]
            for sim_idx in similar_indices:
                if sim_idx not in processed:
                    pos = list(variant_dictionary.keys()).index(sim_idx)
                    dbscan[pos] = cluster_id
                    processed.add(sim_idx)
            
            cluster_id += 1
    
    elif algorithm == 'OPTICS' and HAS_OPTICS:
        cluster_labels = optics_cluster(
            chr_db[variant]["coordinates"], 
            min_samples=args.min_pts if args.DBSCAN else 2, 
            max_eps=args.epsilon if args.DBSCAN else args.bnd_distance
        )
        dbscan = cluster_labels
    
    elif algorithm == 'INTERVAL_TREE' and HAS_INTERVAL_TREE:
        cluster_labels = interval_tree_cluster(
            chr_db[variant]["coordinates"], 
            max_distance=args.bnd_distance
        )
        dbscan = cluster_labels
    
    elif args.DBSCAN:
        cluster_labels = DBSCAN.main(
            chr_db[variant]["coordinates"], 
            args.epsilon, 
            args.min_pts
        )
        dbscan = cluster_labels
    
    else:
        cluster_labels = DBSCAN.main(
            chr_db[variant]["coordinates"], 
            args.bnd_distance, 
            2
        )
        dbscan = cluster_labels
    
    unique_labels = set(dbscan)
    
    unique_xy = chr_db[variant]["coordinates"][dbscan == -1]
    unique_index = chr_db[variant]["index"][dbscan == -1]
    
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

    for unique_label in unique_labels:
        if unique_label == -1:
            continue
        
        class_member_mask = (dbscan == unique_label)
        xy = chr_db[variant]["coordinates"][class_member_mask]
        indexes = chr_db[variant]["index"][class_member_mask]

        if args.DBSCAN or algorithm in ['OPTICS', 'INTERVAL_TREE'] or "INS" in variant:
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
            i = overlap_cluster(db, indexes, variant, chrA, chrB, sample_IDs, args, f, i)

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
    for chrA in chrA_list:
        for chrB in chrB_list:
            for variant in var_list:
                i = svdb_cluster_main(chrA, chrB, variant, sample_IDs, args, db, i)

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