#!/usr/bin/env python3

import pysam
import random
import os

def create_synthetic_bam(output_bam="test_sample.bam", n_reads=1000):

    print(f"Creating synthetic BAM file: {output_bam}")
    chromosomes = [("chr1", 10000000),("chr2", 10000000), ("chr3", 5000000)]
    
    header = {
        'HD': {'VN': '1.0', 'SO': 'coordinate'},
        'SQ': [{'LN': length, 'SN': name} for name, length in chromosomes]
    }
    
    temp_bam = output_bam.replace('.bam', '_temp.bam')
    
    with pysam.AlignmentFile(temp_bam, "wb", header=header) as outf:
        
        read_id = 0

        print("  Generating normal reads...")
        for i in range(n_reads // 2):
            read = create_normal_read(f"read_{read_id}", "chr1", random.randint(1000, 9000000), header)
            outf.write(read)
            read_id += 1
        
        print("  Generating deletion evidence (split reads)...")
        deletion_positions = [ ("chr1", 1000000, 1005000), ("chr1", 2000000, 2010000), ("chr2", 3000000, 3002000) ]
        
        for chrom, start, end in deletion_positions:
            for i in range(5):
                read_name = f"split_{read_id}"
                
                read1 = create_read(read_name, chrom, start - 100 + random.randint(-10, 10), header,is_reverse=False,is_supplementary=False)
                read1.set_tag('SA', f'{chrom},{end + random.randint(-10, 10)},+,100M,60,0')
                outf.write(read1)
                
                read2 = create_read(read_name, chrom,end + random.randint(-10, 10),header,is_reverse=False,is_supplementary=True)
                outf.write(read2)
                
                read_id += 1
    
        print("  Generating deletion evidence (CIGAR)...")
        cigar_del_positions = [ ("chr1", 4000000, 500), ("chr1", 5000000, 1000), ("chr2", 1000000, 200) ]
        
        for chrom, pos, del_size in cigar_del_positions:
            for i in range(3):
                read = create_read_with_deletion(f"cigar_del_{read_id}", chrom, pos + random.randint(-10, 10), del_size, header)
                outf.write(read)
                read_id += 1
        
        print("  Generating insertion evidence (CIGAR)...")
        cigar_ins_positions = [("chr1", 6000000, 150),("chr2", 4000000, 100)]
        
        for chrom, pos, ins_size in cigar_ins_positions:
            for i in range(3):
                read = create_read_with_insertion(f"cigar_ins_{read_id}",chrom, pos + random.randint(-10, 10), ins_size, header)
                outf.write(read)
                read_id += 1
        
        print("  Generating translocation evidence...")
        translocations = [("chr1", 7000000, "chr2", 6000000), ("chr1", 8000000, "chr3", 2000000)]
        
        for chr1, pos1, chr2, pos2 in translocations:
            for i in range(4):
                read_name = f"trans_{read_id}"
                read1 = create_read(read_name, chr1, pos1 + random.randint(-10, 10), header, is_reverse=False, is_supplementary=False)
                read1.set_tag('SA', f'{chr2},{pos2 + random.randint(-10, 10)},+,100M,60,0')
                outf.write(read1)
                read2 = create_read(read_name,chr2, pos2 + random.randint(-10, 10), header, is_reverse=False, is_supplementary=True)
                outf.write(read2)
                
                read_id += 1
    
    print(f"  Generated {read_id} reads")
    
    print("  Sorting BAM file...")
    pysam.sort("-o", output_bam, temp_bam)
    
    print("  Indexing BAM file...")
    pysam.index(output_bam)
    
    os.remove(temp_bam)
    
    size_mb = os.path.getsize(output_bam) / (1024 * 1024)
    
def create_normal_read(read_name, chrom, pos, header_dict):
    read = pysam.AlignedSegment()
    read.query_name = read_name
    
    # Get chromosome index from header
    chrom_id = None
    for i, sq in enumerate(header_dict['SQ']):
        if sq['SN'] == chrom:
            chrom_id = i
            break
    
    read.reference_id = chrom_id
    read.reference_start = pos
    read.query_sequence = "ACGT" * 25  # 100bp read
    read.flag = 0
    read.mapping_quality = 60
    read.cigar = [(0, 100)]  # 100M
    read.query_qualities = [30] * 100
    return read

def create_read(read_name, chrom, pos, header_dict, is_reverse=False, is_supplementary=False):
    read = pysam.AlignedSegment()
    read.query_name = read_name
    
    # Get chromosome index from header
    chrom_id = None
    for i, sq in enumerate(header_dict['SQ']):
        if sq['SN'] == chrom:
            chrom_id = i
            break
    
    read.reference_id = chrom_id
    read.reference_start = pos
    read.query_sequence = "ACGT" * 25  # 100bp read
    read.mapping_quality = 60
    read.cigar = [(0, 100)]  # 100M
    read.query_qualities = [30] * 100
    
    read.flag = 0
    if is_reverse:
        read.flag |= 0x10  # Read reverse strand
    if is_supplementary:
        read.flag |= 0x800  # Supplementary alignment
    
    return read

def create_read_with_deletion(read_name, chrom, pos, deletion_size, header_dict):
    read = pysam.AlignedSegment()
    read.query_name = read_name
    
    # Get chromosome index from header
    chrom_id = None
    for i, sq in enumerate(header_dict['SQ']):
        if sq['SN'] == chrom:
            chrom_id = i
            break
    
    read.reference_id = chrom_id
    read.reference_start = pos
    

    read.query_sequence = "ACGT" * 25  
    read.mapping_quality = 60
    
    read.cigar = [(0, 50), (2, deletion_size), (0, 50)] 
    
    read.query_qualities = [30] * 100
    read.flag = 0
    
    return read

def create_read_with_insertion(read_name, chrom, pos, insertion_size, header_dict):
    read = pysam.AlignedSegment()
    read.query_name = read_name
    
    chrom_id = None
    for i, sq in enumerate(header_dict['SQ']):
        if sq['SN'] == chrom:
            chrom_id = i
            break
    
    read.reference_id = chrom_id
    read.reference_start = pos
    
    seq_length = 100 + insertion_size
    base_seq = "ACGT" * (seq_length // 4 + 1)  
    read.query_sequence = base_seq[:seq_length]  
    read.mapping_quality = 60
    

    read.cigar = [(0, 50), (1, insertion_size), (0, 50)]  
    read.query_qualities = [30] * len(read.query_sequence)
    read.flag = 0
    
    return read

def validate_bam(bam_file):
    print(f"\nValidating {bam_file}...")
    
    try:
        samfile = pysam.AlignmentFile(bam_file, "rb")

        n_reads = 0
        n_split = 0
        n_cigar_indels = 0
        
        for read in samfile.fetch():
            n_reads += 1
            
            # Check for split reads
            if read.has_tag('SA'):
                n_split += 1

            if read.cigar:
                for op, length in read.cigar:
                    if op in [1, 2] and length >= 50:  # Insertion or Deletion
                        n_cigar_indels += 1
                        break
        
        samfile.close()
        
        print(f"✓ Total reads: {n_reads}")
        print(f"✓ Split reads: {n_split}")
        print(f"✓ Reads with large CIGAR indels: {n_cigar_indels}")
        print("\n✓ BAM file is valid!")
        
        return True
        
    except Exception as e:
        print(f" Error validating BAM: {e}")
        return False

if __name__ == "__main__":
    create_synthetic_bam("test_sample.bam", n_reads=1000)
    validate_bam("test_sample.bam")
    