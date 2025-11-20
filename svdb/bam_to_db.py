#!/usr/bin/env python3


import sys
import os

# Import the BAM reading module
sys.path.insert(0, os.path.dirname(__file__))
from readBAM import read_bam_file

def write_vcf_header(f, sample_name):
    f.write("##fileformat=VCFv4.2\n")
    f.write("##source=SVDB_BAM_extraction\n")
    f.write("##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of SV\">\n")
    f.write("##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position\">\n")
    f.write("##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"SV length\">\n")
    f.write("##INFO=<ID=SUPPORT,Number=1,Type=Integer,Description=\"Number of supporting reads\">\n")
    f.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
    f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{}\n".format(sample_name))


def variants_to_vcf(variants, output_vcf, sample_name):
    
    print(f"Writing {len(variants)} variants to {output_vcf}")
    
    with open(output_vcf, 'w') as f:
        write_vcf_header(f, sample_name)
        for i, (chrA, posA, chrB, posB, sv_type, INFO, FORMAT) in enumerate(variants, 1):
            info_parts = [f"SVTYPE={sv_type}"]
            
            if "END" in INFO:
                info_parts.append(f"END={INFO['END']}")
            elif chrA == chrB and posB != posA:
                info_parts.append(f"END={posB}")
            
            if "SVLEN" in INFO:
                info_parts.append(f"SVLEN={INFO['SVLEN']}")
            elif chrA == chrB:
                svlen = abs(posB - posA)
                if sv_type != "INS":
                    info_parts.append(f"SVLEN={svlen}")
            
            if "SUPPORT" in INFO:
                info_parts.append(f"SUPPORT={INFO['SUPPORT']}")
            
            info_str = ";".join(info_parts)
            if sv_type == "DEL":
                alt = "<DEL>"
            elif sv_type == "DUP":
                alt = "<DUP>"
            elif sv_type == "INS":
                alt = "<INS>"
            elif sv_type == "INV":
                alt = "<INV>"
            elif sv_type == "BND":

                alt = f"N[{chrB}:{posB}["
            else:
                alt = f"<{sv_type}>"
            
            # Get genotype
            gt = FORMAT.get("GT", ["0/1"])[0] if FORMAT else "0/1"
            
            # Write VCF line
            vcf_line = [chrA, str(posA), f"sv_{i}", "N", alt,".", "PASS", info_str, "GT", gt]
            
            f.write("\t".join(vcf_line) + "\n")



def main():
    import argparse
    
    parser = argparse.ArgumentParser(description="Extract SVs from BAM and build SVDB database")
    parser.add_argument("bam_file", help="Input BAM file")
    parser.add_argument("--prefix", default="output", help="Output prefix (default: output)")
    parser.add_argument("--sample", default="sample", help="Sample name (default: sample)")
    parser.add_argument("--min-sv-size", type=int, default=50, help="Minimum SV size (default: 50)")
    parser.add_argument("--min-mapq", type=int, default=20, help="Minimum mapping quality (default: 20)")
    
    args = parser.parse_args()
    
    
    print(f"Extracting SVs from BAM: {args.bam_file}")
    print(f" Sample: {args.sample}")
    print(f" Min SV size: {args.min_sv_size} bp")
    print(f" Min MAPQ: {args.min_mapq}")
    print()
    
    try:
        variants = read_bam_file(
            args.bam_file,
            args.sample,
            min_sv_size=args.min_sv_size,
            min_mapq=args.min_mapq
        )
    except Exception as e:
        print(f"Error extracting SVs: {e}")
        return 1
    
    if not variants:
        print("No SVs extracted from BAM file")
        return 1
    
    print(f"Extracted {len(variants)} variants")
    print()
    
    print(f"Converting to VCF format")
    vcf_file = f"{args.prefix}_from_bam.vcf"
    
    try:
        variants_to_vcf(variants, vcf_file, args.sample)
    except Exception as e:
        print(f"✗ Error writing VCF: {e}")
        return 1
    
    print()
    
    print(f"Building SVDB database")
    print(f" Running: svdb --build --files {vcf_file} --prefix {args.prefix}")
    print()
    
    import subprocess
    try:
        result = subprocess.run(
            ["svdb", "--build", "--files", vcf_file, "--prefix", args.prefix],
            check=True,
            capture_output=True,
            text=True
        )
        print(result.stdout)
        if result.stderr:
            print(result.stderr)
    except subprocess.CalledProcessError as e:
        print(f" Error building database: {e}")
        print(e.stdout)
        print(e.stderr)
        return 1
    except FileNotFoundError:
        print("svdb command not found")
        print("   Try: python -m svdb --build --files {} --prefix {}".format(vcf_file, args.prefix))
        return 1

    
    return 0


if __name__ == "__main__":
    sys.exit(main())
