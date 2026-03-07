#!/usr/bin/env python3
"""
Levenshtein test suite for insertion sequences: edit distance, thresholding and position+sequence overlap.
"""

import sys
import os

# Keep path setup minimal
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))


def levenshtein_distance(seq1, seq2):
    """Return (raw_distance, normalized_distance) between two strings."""
    if seq1 is None or seq2 is None:
        return (float('inf'), 1.0)
    if len(seq1) == 0 or len(seq2) == 0:
        return (float('inf'), 1.0)

    seq1 = str(seq1).upper()
    seq2 = str(seq2).upper()
    len1, len2 = len(seq1), len(seq2)

    dp = [[0] * (len2 + 1) for _ in range(len1 + 1)]
    for i in range(len1 + 1):
        dp[i][0] = i
    for j in range(len2 + 1):
        dp[0][j] = j

    for i in range(1, len1 + 1):
        for j in range(1, len2 + 1):
            if seq1[i - 1] == seq2[j - 1]:
                dp[i][j] = dp[i - 1][j - 1]
            else:
                dp[i][j] = 1 + min(
                    dp[i - 1][j],      # deletion
                    dp[i][j - 1],      # insertion
                    dp[i - 1][j - 1]   # substitution
                )

    raw_distance = dp[len1][len2]
    max_len = max(len1, len2)
    normalized_distance = raw_distance / max_len if max_len > 0 else 0.0
    return (raw_distance, normalized_distance)


def compare_insertion_sequences_lv(seq_query, seq_db, max_levenshtein_distance=0.2):
    """Compare two insertion sequences: return (similarity_score, is_match)."""
    if seq_query is None or seq_db is None:
        return (None, False)
    if len(seq_query) == 0 or len(seq_db) == 0:
        return (None, False)

    raw_dist, norm_dist = levenshtein_distance(seq_query, seq_db)
    if raw_dist == float('inf'):
        return (None, False)

    is_match = norm_dist <= max_levenshtein_distance
    similarity = 1.0 - norm_dist
    return (similarity, is_match)


def insertion_overlap_with_sequence_lv(chrApos_query, chrBpos_query, seq_query,
                                       chrApos_db, chrBpos_db, seq_db,
                                       distance, max_levenshtein=0.2):
    """Position + sequence overlap using Levenshtein; returns (combined_similarity, is_match) or (None, False)."""
    Adist = abs(chrApos_query - chrApos_db)
    Bdist = abs(chrBpos_query - chrBpos_db)
    pos_dist = max(Adist, Bdist)
    pos_match = pos_dist <= distance
    if not pos_match:
        return (None, False)

    seq_similarity, seq_match = compare_insertion_sequences_lv(seq_query, seq_db, max_levenshtein)
    overall_match = pos_match and seq_match

    if overall_match:
        pos_similarity = 1.0 - (pos_dist / distance) if distance > 0 else 0.0
        combined_similarity = (pos_similarity + seq_similarity) / 2.0
        return (combined_similarity, True)
    else:
        return (None, False)



def run_tests():
    # 1: Identical sequences
    print("\n[Test 1] Identical sequences")
    raw, norm = levenshtein_distance("ACGTACGT", "ACGTACGT")
    print(f"raw={raw}, norm={norm}")
    assert raw == 0 and norm == 0.0

    # 2: Single substitution
    print("\n[Test 2] Single substitution")
    raw, norm = levenshtein_distance("ACGTACGT", "ACTTACGT")
    print(f"raw={raw}, norm={norm}")
    assert raw == 1 and abs(norm - 0.125) < 1e-6

    # 3: Single deletion
    print("\n[Test 3] Single deletion")
    raw, norm = levenshtein_distance("ACGTACGT", "ACGACGT")
    print(f"raw={raw}, norm={norm}")
    assert raw == 1

    # 4: Single insertion
    print("\n[Test 4] Single insertion")
    raw, norm = levenshtein_distance("ACGTACGT", "ACGTTACGT")
    print(f"raw={raw}, norm={norm}")
    assert raw == 1

    # 5: Multiple edits
    print("\n[Test 5] Multiple edits")
    raw, norm = levenshtein_distance("ACGTACGT", "ACTTACTT")
    print(f"raw={raw}, norm={norm}")
    assert raw == 2 and abs(norm - 0.25) < 1e-6

    # 6: Completely different sequences
    print("\n[Test 6] Completely different sequences")
    raw, norm = levenshtein_distance("AAAA", "TTTT")
    print(f"raw={raw}, norm={norm}")
    assert raw == 4 and norm == 1.0

    # 7: Empty sequence handling
    print("\n[Test 7] Empty sequence handling")
    raw, norm = levenshtein_distance("ACGT", "")
    print(f"raw={raw}, norm={norm}")
    assert raw == float('inf') and norm == 1.0

    # 8: Similar insertions match at threshold
    print("\n[Test 8] Similar insertions threshold")
    similarity, match = compare_insertion_sequences_lv("ACGTACGT", "ACTTACGT", 0.2)
    print(f"similarity={similarity:.3f}, match={match}")
    assert match is True and similarity > 0.8

    # 9: Dissimilar insertions do not match
    print("\n[Test 9] Dissimilar insertions threshold")
    similarity, match = compare_insertion_sequences_lv("ACGTACGT", "TGCATGCA", 0.2)
    print(f"similarity={similarity:.3f}, match={match}")
    assert match is False

    # 10: Position + sequence matching (close + similar)
    print("\n[Test 10] Position + sequence overlap")
    similarity, match = insertion_overlap_with_sequence_lv(
        1000, 1000, "ACGTACGT",
        1005, 1005, "ACTTACGT",
        distance=10, max_levenshtein=0.2
    )
    print(f"similarity={similarity:.3f}, match={match}")
    assert match is True and similarity is not None

    # 11: Insertion vs deletion difference
    print("\n[Test 11] Insertion recognized as single edit")
    raw_lv, norm_lv = levenshtein_distance("ATCGATCG", "ATCGAATCG")
    print(f"raw={raw_lv}, norm={norm_lv}")
    assert raw_lv == 1

    # 12: Long sequences with indels
    print("\n[Test 12] Long-read scenario with indels")
    original = "ATCGATCGATCGATCGATCG"
    with_errors = "ATCGATCGTTCGATCGACG"
    raw, norm = levenshtein_distance(original, with_errors)
    similarity, match = compare_insertion_sequences_lv(original, with_errors, 0.15)
    print(f"raw={raw}, norm={norm:.3f}, similarity={similarity:.3f}, match={match}")
    assert raw == 2 and abs(norm - 0.1) < 0.02

    print("\nAll tests completed.")
    return True



if __name__ == "__main__":
    ok = False
    try:
        ok = run_tests()
    finally:
        sys.exit(0 if ok else 1)