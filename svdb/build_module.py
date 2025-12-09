from __future__ import absolute_import

import glob
import gzip
import os

from . import database
from . import readVCF
# from . import readBAM  # BAM integration commented out - BCF support via pysam instead


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
        sample_name = input_file.split("/")[-1].split(".vcf")[0].split(".bcf")[0]  # .split(".bam")[0]
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
        
        # BAM processing commented out - BCF support implemented instead
        # if input_file.endswith('.bam'):
        #     print(f'Processing BAM file: {input_file}')
        #     variants = readBAM.read_bam_file(input_file, sample_name, 
        #                                      min_sv_size=getattr(args, "min_sv_size", 50),
        #                                      min_mapq=getattr(args, "min_mapq", 20))
        #     
        #     for chrA, posA, chrB, posB, event_type, INFO, FORMAT, _ in variants:
        #         # ... BAM processing code ...
        
        if input_file.endswith('.bcf'):
            # Use pysam to read BCF files
            import pysam
            print(f'Processing BCF file: {input_file}')
            
            try:
                vcf_reader = pysam.VariantFile(input_file, 'rb')
                sample_names = list(vcf_reader.header.samples) if vcf_reader.header.samples else []
                
                for record in vcf_reader.fetch():
                    if hasattr(args, 'passonly') and args.passonly:
                        if record.filter.keys() and 'PASS' not in record.filter.keys() and '.' not in record.filter.keys():
                            continue
                    
                    # Extract basic variant info
                    chrA = str(record.chrom).replace("chr", "").replace("Chr", "").replace("CHR", "")
                    posA = record.pos
                    
                    # Determine variant type and positions
                    event_type = record.info.get('SVTYPE', 'UNK') if 'SVTYPE' in record.info else 'UNK'
                    
                    # Handle different ALT formats
                    if record.alts:
                        alt = str(record.alts[0])
                        if '<' in alt and '>' in alt:
                            event_type = alt.strip('<').rstrip('>')
                            if 'DUP' in event_type:
                                event_type = 'DUP'
                    
                    chrB = chrA
                    posB = record.info.get('END', posA) if 'END' in record.info else posA
                    
                    # Handle BND/translocation
                    if 'CHR2' in record.info:
                        chrB = str(record.info['CHR2']).replace("chr", "").replace("Chr", "").replace("CHR", "")
                        event_type = 'BND'
                    
                    # Confidence intervals
                    ci_A_lower = ci_A_upper = ci_B_lower = ci_B_upper = 0
                    if 'CIPOS' in record.info:
                        cipos = record.info['CIPOS']
                        if isinstance(cipos, (list, tuple)) and len(cipos) >= 2:
                            ci_A_lower = abs(int(cipos[0]))
                            ci_A_upper = abs(int(cipos[1]))
                            ci_B_lower = ci_A_lower
                            ci_B_upper = ci_A_upper
                    
                    if 'CIEND' in record.info:
                        ciend = record.info['CIEND']
                        if isinstance(ciend, (list, tuple)) and len(ciend) >= 2:
                            ci_B_lower = abs(int(ciend[0]))
                            ci_B_upper = abs(int(ciend[1]))
                    
                    # Extract sequence for insertions
                    sequence = ""
                    if "INS" in event_type:
                        if record.alts and '<INS>' not in str(record.alts[0]):
                            sequence = str(record.alts[0])
                        elif 'INSSEQ' in record.info:
                            sequence = str(record.info['INSSEQ'])
                        elif 'SEQ' in record.info:
                            sequence = str(record.info['SEQ'])
                    
                    # Add variants based on genotypes
                    if sample_names and record.samples:
                        for sample_idx, sample in enumerate(sample_names):
                            gt = record.samples[sample].get('GT', None)
                            if gt and gt != (0, 0) and gt != (None, None):
                                var.append((event_type, chrA, chrB, posA, ci_A_lower, ci_A_upper,
                                           posB, ci_B_lower, ci_B_upper, sample, idx, sequence))
                                idx += 1
                    else:
                        var.append((event_type, chrA, chrB, posA, ci_A_lower, ci_A_upper,
                                   posB, ci_B_lower, ci_B_upper, sample_name, idx, sequence))
                        idx += 1
                
                vcf_reader.close()
            except Exception as e:
                print(f"Error processing BCF file {input_file}: {e}")
                continue
        
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
        
        else:
            print(f"Error, the file format of {input_file} is not supported. Only .vcf, .vcf.gz and .bcf are supported.")  # .bam removed
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
                     glob.glob(os.path.join(args.folder, "*.bcf"))  # + glob.glob(os.path.join(args.folder, "*.bam"))
    populate_db(args)