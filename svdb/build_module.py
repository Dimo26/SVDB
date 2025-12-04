from __future__ import absolute_import

import glob
import gzip
import os

from . import database
from . import readVCF
from . import readBAM


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


def _prefilter_vcf_insertions_by_sequence(variants_list, max_hamming=0.2):
    """Pre-filter VCF insertion variants by sequence similarity."""
    insertions = [(i, v) for i, v in enumerate(variants_list) if len(v) > 11 and v[4] == 'INS' and v[11]]
    
    if not insertions:
        return
    
    print(f"\n=== VCF Insertion Sequence Pre-filtering (Hamming distance) ===")
    print(f"Found {len(insertions)} insertions with sequences")
    
    # Build groups
    groups = []
    assigned = set()
    
    for i, (idx_i, var_i) in enumerate(insertions):
        if i in assigned:
            continue
        group = [i]
        assigned.add(i)
        seq_i = var_i[11]
        
        for j in range(i + 1, len(insertions)):
            if j in assigned:
                continue
            idx_j, var_j = insertions[j]
            dist = _hamming_distance(seq_i, var_j[11])
            if dist <= max_hamming:
                group.append(j)
                assigned.add(j)
        groups.append(group)
    
    # Print groups
    for group_id, group_indices in enumerate(groups):
        if len(group_indices) > 1:
            print(f"  Group {group_id}: {len(group_indices)} insertions with similar sequences")
            for idx in group_indices[:3]:  # Show first 3
                seq = str(insertions[idx][1][11])[:20]
                print(f"    - {seq}...")
            if len(group_indices) > 3:
                print(f"    ... and {len(group_indices) - 3} more")
        else:
            seq = str(insertions[group_indices[0]][1][11])[:20]
            print(f"  Group {group_id}: {seq}... (singleton)")
    
    print(f"VCF insertions grouped into {len(groups)} sequence clusters\n")


def populate_db(args):
    db = database.DB(args.db)
    tables = db.tables

    idx = 0

    if "SVDB" in tables:
        db.drop("DROP TABLE IF EXISTS SVDB")
        db.drop("DROP INDEX IF EXISTS SV")
        db.drop("DROP INDEX IF EXISTS IDX")
        db.drop("DROP INDEX IF EXISTS CHR")


    query = "CREATE TABLE SVDB (var TEXT, chrA TEXT, chrB TEXT, posA INT, ci_A_lower INT, ci_A_upper INT, posB INT, ci_B_lower INT, ci_B_upper INT, sample TEXT, idx INT, sequence TEXT)"
    db.create(query)
    sample_IDs = []

    # Populate the tables
    for input_file in args.files:
        sample_name = input_file.split("/")[-1].split(".vcf")[0].split(".bam")[0]
        sample_name = sample_name.replace(".", "_")
        sample_IDs.append(sample_name)
        
        A = 'SELECT sample FROM SVDB WHERE sample == \'{}\' '.format(sample_name)
        hits = [hit for hit in db.query(A)]
        if hits:
            continue
            
        if not os.path.exists(input_file):
            print("error: unable to open {}".format(input_file))
            continue

        var = []
        
        if input_file.endswith('.bam'):
            print(f'Processing BAM file: {input_file}')
            variants = readBAM.read_bam_file(input_file, sample_name, 
                                             min_sv_size=getattr(args, "min_sv_size", 50),
                                             min_mapq=getattr(args, "min_mapq", 20))
            
            for chrA, posA, chrB, posB, event_type, INFO, FORMAT, _ in variants:
                if hasattr(args, 'passonly') and args.passonly:
                    pass  # BAM doesn't have FILTER field
                
                ci_A_lower = 0
                ci_A_upper = 0
                ci_B_lower = 0
                ci_B_upper = 0
                sequence = ""
                
                # Get sequence for insertions
                if "INS" in event_type and "INSSEQ" in INFO:
                    sequence = INFO.get("INSSEQ", "")
                
                if "GT" not in FORMAT:
                    var.append((event_type, chrA, chrB, posA, ci_A_lower, ci_A_upper, 
                               posB, ci_B_lower, ci_B_upper, sample_name, idx, sequence))
                    idx += 1
                else:
                    for genotype in FORMAT["GT"]:
                        if genotype not in ["0/0", "./."]:
                            var.append((event_type, chrA, chrB, posA, ci_A_lower, ci_A_upper,
                                       posB, ci_B_lower, ci_B_upper, sample_name, idx, sequence))
                            idx += 1
        
        elif input_file.endswith('.vcf') or input_file.endswith('.vcf.gz'):
            print(f'Processing VCF file: {input_file}')
            vcf = input_file
            sample_names = []

            opener = gzip.open if vcf.endswith('.vcf.gz') else open
            with opener(vcf, 'rt') as lines:
                for line in lines:
                    if line.startswith("#"):
                        if "CHROM" in line:
                            content = line.strip().split()
                            if len(content) > 9:
                                sample_names = content[9:]
                        continue

                    if not len(line.strip()):
                        continue

                    chrA, posA, chrB, posB, event_type, INFO, FORMAT = readVCF.readVCFLine(line)
                    
                    if hasattr(args, 'passonly') and args.passonly:
                        FILTER = line.split("\t")[6]
                        if not (FILTER in ["PASS", "."]):
                            continue

                    ci_A_lower = 0
                    ci_A_upper = 0
                    ci_B_lower = 0
                    ci_B_upper = 0
                    sequence = ""
                    
                    if "CIPOS" in INFO:
                        ci = INFO["CIPOS"].replace('(','').replace(')','').split(",")
                        if len(ci) > 1:
                            ci_A_lower = abs(int(ci[0]))
                            ci_A_upper = abs(int(ci[1]))
                            ci_B_lower = abs(int(ci[0]))
                            ci_B_upper = abs(int(ci[1]))
                        else:
                            ci_A_lower = abs(int(ci[0]))
                            ci_A_upper = abs(int(ci[0]))
                            ci_B_lower = abs(int(ci[0]))
                            ci_B_upper = abs(int(ci[0]))

                    if "CIEND" in INFO:
                        ci = INFO["CIEND"].replace('(','').replace(')','').split(",")
                        if len(ci) > 1:
                            ci_B_lower = abs(int(ci[0]))
                            ci_B_upper = abs(int(ci[1]))
                        else:
                            ci_B_lower = abs(int(ci[0]))
                            ci_B_upper = abs(int(ci[0]))
                    
                    # Extract sequence for insertions
                    if "INS" in event_type:
                        vcf_columns = line.strip().split("\t")
                        alt_field = vcf_columns[4] if len(vcf_columns) > 4 else ""

                        if "<INS>" not in alt_field and alt_field not in ["", ".", "N"]:
                            sequence = alt_field
                        elif "INSSEQ" in INFO:
                            sequence = INFO["INSSEQ"]
                        elif "SEQ" in INFO:
                            sequence = INFO["SEQ"]
                    
                    if "GT" not in FORMAT or not len(sample_names):
                        var.append((event_type, chrA, chrB, posA, ci_A_lower,
                                    ci_A_upper, posB, ci_B_lower, ci_B_upper, sample_name, idx, sequence))
                        idx += 1
                    else:
                        sample_index = 0
                        for genotype in FORMAT["GT"]:
                            if genotype not in ["0/0", "./."]:
                                var.append((event_type, chrA, chrB, posA, ci_A_lower, ci_A_upper,
                                            posB, ci_B_lower, ci_B_upper, sample_names[sample_index], idx, sequence))
                                idx += 1
                            sample_index += 1
            
            # Pre-filter VCF insertions by sequence similarity (after all variants collected)
            if var:
                _prefilter_vcf_insertions_by_sequence(var, max_hamming=0.2)
        
        else:
            print(f"Error, the file format of {input_file} is not supported. Only .vcf, .vcf.gz and .bam are supported.")
            continue

        if var:
            db.insert_many(var)

    # Create indexes
    db.create_index(name='SV', columns='(var, chrA, chrB, posA, posB)')
    db.create_index(name='IDX', columns='(idx)')
    db.create_index(name='CHR', columns='(chrA, chrB)')
    
    return sample_IDs

def main(args):
    args.db = args.prefix
    if not args.files and hasattr(args, 'folder') and args.folder:
        args.files = glob.glob(os.path.join(args.folder, "*.vcf")) + \
                     glob.glob(os.path.join(args.folder, "*.vcf.gz")) + \
                     glob.glob(os.path.join(args.folder, "*.bam"))
    populate_db(args)