#!/usr/bin/env python3

import sys
import os
import sqlite3
import numpy as np


def load_and_stats(db_file):
    conn = sqlite3.connect(db_file)
    cur = conn.cursor()

    # Exact same filter as benchmark
    cur.execute('SELECT var, posA, posB, sequence FROM SVDB WHERE chrA = "1" AND chrB = "1"')
    rows = cur.fetchall()
    conn.close()

    if not rows:
        print(f"  No variants found for chrA=1, chrB=1")
        return

    stats = {}
    for var_type, posA, posB, seq in rows:
        if var_type not in stats:
            stats[var_type] = {'count': 0, 'sizes': [], 'with_seq': 0, 'without_seq': 0}

        stats[var_type]['count'] += 1

        # EXACT same size logic as benchmark
        if var_type == 'INS' and seq and len(seq) > 0:
            sv_size = len(seq)
            stats[var_type]['with_seq'] += 1
        else:
            sv_size = abs(posB - posA)
            if var_type == 'INS':
                stats[var_type]['without_seq'] += 1

        stats[var_type]['sizes'].append(sv_size)

    # Print in same format as benchmark
    total = len(rows)
    print(f"{'STRUCTURAL VARIANT STATISTICS - Chromosome 1':^90}")
    print(f"{'SV TYPE':<10} {'COUNT':>8} {'%':>7} {'AVG_SIZE':>12} {'MIN':>10} {'MAX':>12} {'WITH_SEQ':>10} {'NO_SEQ':>8}")
    print("-" * 90)

    for sv_type in sorted(stats.keys()):
        d = stats[sv_type]
        pct = d['count'] / total * 100
        avg  = np.mean(d['sizes'])
        mn   = min(d['sizes'])
        mx   = max(d['sizes'])
        print(f"{sv_type:<10} {d['count']:>8} {pct:>6.1f}% {avg:>12.1f} "
              f"{mn:>10} {mx:>12} {d['with_seq']:>10} {d['without_seq']:>8}")

    print(f"{'TOTAL':<10} {total:>8} {'100.0%':>7}")

    # Insertion detail block (same as benchmark)
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


def main():
    if len(sys.argv) < 2:
        print("Usage: python3 verify_benchmark_stats.py <file.db> [file2.db ...]")
        sys.exit(1)

    for db_file in sys.argv[1:]:
        if not os.path.exists(db_file):
            print(f"ERROR: {db_file} not found")
            continue

        print(f"\n{'=' * 90}")
        print(f"  DATABASE: {os.path.basename(db_file)}")
        print(f"{'=' * 90}")
        load_and_stats(db_file)


if __name__ == "__main__":
    main()