#!/usr/bin/env python3

import sys
import gzip
import numpy as np


def parse_info(info_str):
    """Parse INFO field into dict — same logic as readVCF.py lines 82-87."""
    description = {}
    for tag in info_str.split(";"):
        parts = tag.split("=")
        if len(parts) > 1:
            description[parts[0]] = parts[1]
    return description


def parse_vcf_line(line):
    cols = line.strip().split("\t")
    if len(cols) < 8:
        return None

    chrA = cols[0].replace("chr", "").replace("Chr", "").replace("CHR", "")
    posA = int(cols[1])
    posB = 0
    event_type = ""
    sequence = ""
    alt = cols[4]
    ref = cols[3]

    INFO = parse_info(cols[7])

    # BND with brackets (interchromosomal, won't pass chr1 filter)
    if "]" in alt or "[" in alt:
        event_type = "BND"
        chrB = "other"
        posB = 0
        genotypes = _extract_genotypes(cols)
        return {'chrA': chrA, 'posA': posA, 'chrB': chrB, 'posB': posB,
                'type': event_type, 'sequence': sequence, 'genotypes': genotypes}


    if "TRA" in alt:
        event_type = "BND"
        chrB = INFO.get("CHR2", chrA)
        posB = int(INFO.get("END", posA))
        if chrA > chrB:
            chrA, chrB = chrB, chrA
            posB = posA
        genotypes = _extract_genotypes(cols)
        return {'chrA': chrA, 'posA': posA, 'chrB': chrB, 'posB': posB,
                'type': event_type, 'sequence': sequence, 'genotypes': genotypes}


    chrB = chrA
    posB = posA

    if "END" in INFO:
        posB = int(INFO["END"])
    elif "SVLEN" in INFO:
        posB = posA + abs(int(INFO["SVLEN"]))

    if posB < posA:
        posA, posB = posB, posA

    nucleotides = set(["A", "T", "C", "G", "N"])

    if "<" in alt and ">" in alt:
        event_type = alt.strip("<").rstrip(">")
        if "DUP" in event_type:
            event_type = "DUP"
    elif "SVTYPE" in INFO:
        event_type = INFO["SVTYPE"]
    elif set(list(alt)).union(nucleotides) == nucleotides:
        if len(alt) > len(ref):
            event_type = "INS"
        else:
            event_type = "DEL"
            posB = posA + len(ref) - 1

    if "INS" in event_type:
        posA = int(cols[1])
        posB = int(cols[1])
        if "<INS>" not in alt and alt not in ["", ".", "N"]:
            sequence = alt
        elif "INSSEQ" in INFO:
            sequence = INFO["INSSEQ"]
        elif "SEQ" in INFO:
            sequence = INFO["SEQ"]

    genotypes = _extract_genotypes(cols)

    return {'chrA': chrA, 'posA': posA, 'chrB': chrB, 'posB': posB,
            'type': event_type, 'sequence': sequence, 'genotypes': genotypes}


def _extract_genotypes(cols):
    """Extract GT values from FORMAT/SAMPLE columns."""
    genotypes = []
    if len(cols) > 9:
        format_keys = cols[8].split(":")
        if "GT" in format_keys:
            gt_index = format_keys.index("GT")
            for sample_col in cols[9:]:
                fields = sample_col.split(":")
                if gt_index < len(fields):
                    genotypes.append(fields[gt_index])
                else:
                    genotypes.append("./.")
    return genotypes


def verify_vcf(vcf_file):
    opener = gzip.open if vcf_file.endswith('.vcf.gz') else open

    sample_names = []
    stats = {}
    total_lines = 0
    skipped_chr = 0
    skipped_gt = 0

    with opener(vcf_file, 'rt') as f:
        for line in f:
            if line.startswith("#"):
                if "#CHROM" in line:
                    cols = line.strip().split("\t")
                    if len(cols) > 9:
                        sample_names = cols[9:]
                continue

            if not line.strip():
                continue

            total_lines += 1
            parsed = parse_vcf_line(line)
            if parsed is None:
                continue

            if parsed['chrA'] != '1' or parsed['chrB'] != '1':
                skipped_chr += 1
                continue

            sv_type = parsed['type']
            if not sv_type:
                continue

            if sv_type not in stats:
                stats[sv_type] = {'count': 0, 'sizes': [], 'with_seq': 0, 'without_seq': 0}

            if parsed['genotypes']:
                any_passed = False
                for gt in parsed['genotypes']:
                    if gt not in ["0/0", "./."]:
                        any_passed = True
                        break
                if not any_passed:
                    skipped_gt += 1
                    continue


            if sv_type == 'INS' and parsed['sequence'] and len(parsed['sequence']) > 0:
                size = len(parsed['sequence'])
                stats[sv_type]['with_seq'] += 1
            else:
                size = abs(parsed['posB'] - parsed['posA'])
                if sv_type == 'INS':
                    stats[sv_type]['without_seq'] += 1

            stats[sv_type]['count'] += 1
            stats[sv_type]['sizes'].append(size)


    print(f"  VCF INPUT VERIFICATION: {vcf_file}")
    print(f"  Samples found: {sample_names if sample_names else '(no sample columns)'}")
    print(f"  Total variant lines in file: {total_lines}")
    print(f"  Skipped (not chr1 intrachromosomal): {skipped_chr}")
    print(f"  Skipped (all genotypes are 0/0 or ./.): {skipped_gt}")
    print(f"  PREDICTED .db CONTENT — run verify_benchmark_stats.py on your .db to compare")


    total = sum(d['count'] for d in stats.values())

    print(f"\n{'SV TYPE':<10} {'COUNT':>8} {'%':>7} {'AVG_SIZE':>12} {'MIN':>10} {'MAX':>12} {'WITH_SEQ':>10} {'NO_SEQ':>8}")


    for sv_type in sorted(stats.keys()):
        d = stats[sv_type]
        pct = d['count'] / total * 100 if total > 0 else 0
        avg = np.mean(d['sizes']) if d['sizes'] else 0
        mn  = min(d['sizes']) if d['sizes'] else 0
        mx  = max(d['sizes']) if d['sizes'] else 0
        print(f"{sv_type:<10} {d['count']:>8} {pct:>6.1f}% {avg:>12.1f} "
              f"{mn:>10} {mx:>12} {d['with_seq']:>10} {d['without_seq']:>8}")

    print(f"{'TOTAL':<10} {total:>8} {'100.0%':>7}")

    if 'INS' in stats:
        ins = stats['INS']
        print(f"\n{'INSERTION DETAILS':^90}")
        print(f"  Total insertions: {ins['count']}")
        if ins['count'] > 0:
            print(f"  With sequence:    {ins['with_seq']:>5} ({ins['with_seq']/ins['count']*100:>5.1f}%)")
            print(f"  Without sequence: {ins['without_seq']:>5} ({ins['without_seq']/ins['count']*100:>5.1f}%)")
            seq_sizes = [s for s in ins['sizes'] if s > 0]
            if seq_sizes:
                print(f"  Sequence lengths:")
                print(f"    Mean:   {np.mean(seq_sizes):>8.1f} bp")
                print(f"    Median: {np.median(seq_sizes):>8.1f} bp")
                print(f"    Range:  {min(seq_sizes)}-{max(seq_sizes)} bp")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python3 verify_vcf_input.py sample.vcf [sample2.vcf ...]")
        sys.exit(1)

    for vcf in sys.argv[1:]:
        verify_vcf(vcf)
