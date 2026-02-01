#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Quick verification script for Platinum sample databases
"""

import sys
import sqlite3
from pathlib import Path

def verify_database(db_path):
    """Run comprehensive checks on a database"""
    db_name = Path(db_path).stem
    print(f"VERIFYING: {db_name}")

    
    try:
        conn = sqlite3.connect(db_path)
        cursor = conn.cursor()
        
        # Check 1: Table exists
        cursor.execute("SELECT name FROM sqlite_master WHERE type='table'")
        tables = [row[0] for row in cursor.fetchall()]
        print(f"✓ Tables found: {', '.join(tables)}")
        
        if 'SVDB' not in tables:
            print("✗ ERROR: SVDB table not found!")
            return False
        
        # Check 2: Table structure
        cursor.execute("PRAGMA table_info(SVDB)")
        columns = [row[1] for row in cursor.fetchall()]
        expected_columns = ['var', 'chrA', 'chrB', 'posA', 'ci_A_lower', 'ci_A_upper', 
                          'posB', 'ci_B_lower', 'ci_B_upper', 'sample', 'idx', 'sequence']
        
        if set(columns) == set(expected_columns):
            print(f"✓ Table structure correct ({len(columns)} columns)")
        else:
            print(f"⚠ Warning: Column mismatch")
            print(f"  Expected: {expected_columns}")
            print(f"  Found: {columns}")
        
        # Check 3: Total variant count
        cursor.execute("SELECT COUNT(*) FROM SVDB")
        total = cursor.fetchone()[0]
        print(f"✓ Total rows in SVDB: {total:,}")
        
        # Check 4: Unique variants
        cursor.execute("SELECT COUNT(DISTINCT idx) FROM SVDB")
        unique = cursor.fetchone()[0]
        print(f"✓ Unique variants (distinct idx): {unique:,}")
        
        # Check 5: SV type distribution
        print(f"\n--- SV Type Distribution ---")
        cursor.execute("""
            SELECT var, COUNT(*) as count 
            FROM SVDB 
            GROUP BY var 
            ORDER BY count DESC
        """)
        
        sv_counts = {}
        for row in cursor.fetchall():
            svtype, count = row
            sv_counts[svtype] = count
            print(f"  {svtype:6} : {count:>6,} variants")
        
        # Check 6: Sequence field for insertions
        if 'INS' in sv_counts:
            cursor.execute("""
                SELECT 
                    COUNT(*) as total,
                    SUM(CASE WHEN sequence IS NOT NULL AND sequence != '' THEN 1 ELSE 0 END) as with_seq,
                    SUM(CASE WHEN sequence IS NULL OR sequence = '' THEN 1 ELSE 0 END) as without_seq
                FROM SVDB 
                WHERE var = 'INS'
            """)
            total, with_seq, without_seq = cursor.fetchone()
            
            print(f"\n--- Insertion Sequence Storage ---")
            print(f"  Total INS: {total}")
            print(f"  With sequence: {with_seq} ({100*with_seq/total:.1f}%)")
            print(f"  Without sequence: {without_seq} ({100*without_seq/total:.1f}%)")
            
            if without_seq > total * 0.5:
                print(f"  ⚠ Warning: >50% of insertions missing sequence data")
        
        # Check 7: Coordinate validity
        print(f"\n--- Coordinate Validity ---")
        cursor.execute("""
            SELECT COUNT(*) 
            FROM SVDB 
            WHERE posA IS NULL OR posB IS NULL OR posA < 0 OR posB < 0
        """)
        invalid_coords = cursor.fetchone()[0]
        
        if invalid_coords == 0:
            print(f"  ✓ All coordinates valid (no NULL or negative values)")
        else:
            print(f"  ✗ ERROR: {invalid_coords} variants with invalid coordinates")
        
        # Check 8: Chromosome distribution
        print(f"\n--- Chromosome Distribution (Top 5) ---")
        cursor.execute("""
            SELECT chrA, COUNT(*) as count 
            FROM SVDB 
            GROUP BY chrA 
            ORDER BY count DESC 
            LIMIT 5
        """)
        
        for row in cursor.fetchall():
            chrom, count = row
            print(f"  chr{chrom:2} : {count:>6,} variants")
        
        # Check 9: Sample information
        cursor.execute("SELECT COUNT(DISTINCT sample) FROM SVDB")
        num_samples = cursor.fetchone()[0]
        print(f"\n--- Sample Information ---")
        print(f"  Unique samples: {num_samples}")
        
        if num_samples <= 10:
            cursor.execute("""
                SELECT sample, COUNT(*) as count 
                FROM SVDB 
                GROUP BY sample 
                ORDER BY count DESC
            """)
            print(f"  Sample breakdown:")
            for row in cursor.fetchall():
                sample, count = row
                print(f"    {sample}: {count:,} variants")
        
        # Check 10: Index efficiency
        print(f"\n--- Database Indices ---")
        cursor.execute("""
            SELECT name, sql 
            FROM sqlite_master 
            WHERE type = 'index' AND tbl_name = 'SVDB'
        """)
        
        indices = cursor.fetchall()
        if indices:
            print(f"  ✓ Found {len(indices)} indices:")
            for name, sql in indices:
                print(f"    - {name}")
        else:
            print(f"  ⚠ Warning: No indices found (queries will be slow)")
        
        conn.close()
        
        print(f"✓ VERIFICATION COMPLETE: {db_name}")
        
        return True
        
    except sqlite3.Error as e:
        print(f"\n✗ Database error: {e}\n")
        return False
    except Exception as e:
        print(f"\n✗ Unexpected error: {e}\n")
        import traceback
        traceback.print_exc()
        return False

def compare_databases(db_paths):
    """Compare multiple databases to check for merging issues"""
    print("COMPARING DATABASES")

    stats = {}
    
    for db_path in db_paths:
        db_name = Path(db_path).stem
        try:
            conn = sqlite3.connect(db_path)
            cursor = conn.cursor()
            
            # Get SV counts by type
            cursor.execute("""
                SELECT var, COUNT(*) 
                FROM SVDB 
                GROUP BY var
            """)
            
            stats[db_name] = {row[0]: row[1] for row in cursor.fetchall()}
            conn.close()
            
        except Exception as e:
            print(f"Error reading {db_name}: {e}")
            continue
    
    # Print comparison table
    if len(stats) > 1:
        sv_types = set()
        for db_stats in stats.values():
            sv_types.update(db_stats.keys())
        
        sv_types = sorted(sv_types)
        
        print(f"{'Database':<15}", end='')
        for svtype in sv_types:
            print(f"{svtype:>10}", end='')
        print()
        print('-' * 70)
        
        for db_name, db_stats in stats.items():
            print(f"{db_name:<15}", end='')
            for svtype in sv_types:
                count = db_stats.get(svtype, 0)
                print(f"{count:>10,}", end='')
            print()
        
        # Check for merged databases
        print(f"\n--- Checking Merged Databases ---")
        for db_name in stats.keys():
            if '_' in db_name:
                parts = db_name.split('_')
                if len(parts) == 2:
                    db1_name = parts[0]
                    db2_name = parts[1]
                    
                    if db1_name in stats and db2_name in stats:
                        print(f"\nChecking {db_name} = {db1_name} + {db2_name}:")
                        
                        for svtype in sv_types:
                            count1 = stats[db1_name].get(svtype, 0)
                            count2 = stats[db2_name].get(svtype, 0)
                            expected = count1 + count2
                            actual = stats[db_name].get(svtype, 0)
                            
                            if actual == expected:
                                status = "✓"
                            else:
                                status = "✗"
                            
                            print(f"  {status} {svtype}: {actual:,} (expected {expected:,} = {count1:,} + {count2:,})")

def main():
    if len(sys.argv) < 2:
        print("Usage: python3 verify_databases.py <database1.db> [database2.db ...]")
        print("\nThis script will:")
        print("  1. Verify each database structure and content")
        print("  2. Check for common issues")
        print("  3. Compare databases if multiple provided")
        sys.exit(1)
    
    print("SVDB DATABASE VERIFICATION TOOL")

    # Verify each database
    valid_dbs = []
    for db_path in sys.argv[1:]:
        if not Path(db_path).exists():
            print(f"\n✗ Database not found: {db_path}")
            continue
        
        if verify_database(db_path):
            valid_dbs.append(db_path)
    
    # Compare databases if multiple valid ones
    if len(valid_dbs) > 1:
        compare_databases(valid_dbs)
    

    print("VERIFICATION COMPLETE")

    print(f"\nVerified {len(valid_dbs)}/{len(sys.argv)-1} databases successfully")

if __name__ == "__main__":
    main()
