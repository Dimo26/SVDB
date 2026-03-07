#!/usr/bin/env python3

def hamming_distance(seq1, seq2):
    if not seq1 or not seq2:
        return (float('inf'), 1.0)

    seq1 = seq1.upper()
    seq2 = seq2.upper()

    min_len = min(len(seq1), len(seq2))
    max_len = max(len(seq1), len(seq2))

    mismatches = 0
    for i in range(min_len):
        if seq1[i] != seq2[i]:
            mismatches += 1

    length_diff = max_len - min_len
    total_distance = mismatches + length_diff
    normalized = total_distance / max_len if max_len > 0 else 0.0

    return (total_distance, normalized)


def compare_insertion_sequences(seq_query, seq_db, max_hamming_distance=0.2):
    if seq_query is None or seq_db is None:
        return (None, False)

    _, norm_dist = hamming_distance(seq_query, seq_db)
    is_match = norm_dist <= max_hamming_distance
    similarity = 1.0 - norm_dist
    return (similarity, is_match)


def insertion_overlap_with_sequence(chrApos_query, chrBpos_query, seq_query,
                                    chrApos_db, chrBpos_db, seq_db,
                                    distance, max_hamming=0.2):
    Adist = abs(chrApos_query - chrApos_db)
    Bdist = abs(chrBpos_query - chrBpos_db)
    pos_dist = max(Adist, Bdist)
    pos_match = pos_dist <= distance

    if not pos_match:
        return (None, False)

    seq_similarity, seq_match = compare_insertion_sequences(
        seq_query, seq_db, max_hamming
    )

    overall_match = pos_match and seq_match
    if overall_match:
        pos_similarity = 1.0 - (pos_dist / distance)
        combined_similarity = (pos_similarity + seq_similarity) / 2.0
        return (combined_similarity, True)
    else:
        return (None, False)


def run_tests():
    # Test 1: Identical sequences
    raw, norm = hamming_distance("ACGTACGT", "ACGTACGT")
    assert raw == 0
    assert norm == 0.0

    # Test 2: Single mismatch
    raw, norm = hamming_distance("ACGTACGT", "ACTTACGT")
    assert raw == 1
    assert norm == 0.125

    # Test 3: Multiple mismatches (no strict assert on values)
    hamming_distance("ACGTACGT", "TGCATGCA")

    # Test 4: Different lengths
    raw, norm = hamming_distance("ACGTACGT", "ACGT")
    assert raw == 4
    assert norm == 0.5

    # Test 5: Empty sequence handling
    raw, norm = hamming_distance("ACGT", "")
    assert raw == float('inf')
    assert norm == 1.0

    # Test 6: Compare similar insertions (should match)
    similarity, match = compare_insertion_sequences("ACGTACGT", "ACTTACGT", max_hamming_distance=0.2)
    assert match is True
    assert similarity > 0.8

    # Test 7: Compare dissimilar insertions (should NOT match)
    similarity, match = compare_insertion_sequences("ACGTACGT", "TGCATGCA", max_hamming_distance=0.2)
    assert match is False

    # Test 8: Position + sequence matching
    similarity, match = insertion_overlap_with_sequence(
        1000, 1000, "ACGTACGT",
        1005, 1005, "ACTTACGT",
        distance=10, max_hamming=0.2
    )
    assert match is True
    assert similarity is not None

    # Test 9: Close position but different sequence
    similarity, match = insertion_overlap_with_sequence(
        1000, 1000, "ACGTACGT",
        1005, 1005, "TGCATGCA",
        distance=10, max_hamming=0.2
    )
    assert match is False
    assert similarity is None

    # Test 10: Far positions (should fail regardless of sequence)
    similarity, match = insertion_overlap_with_sequence(
        1000, 1000, "ACGTACGT",
        2000, 2000, "ACGTACGT",
        distance=10, max_hamming=0.2
    )
    assert match is False
    assert similarity is None


def usage_examples():
    seq1 = "ACGTACGTACGT"
    seq2 = "ACGTACTTACGT"
    _ = hamming_distance(seq1, seq2)

    ins1 = "ACGTACGT"
    ins2 = "ACTTACGT"
    threshold = 0.2
    _ = compare_insertion_sequences(ins1, ins2, threshold)

    pos1 = 1_000_000
    seqA = "ACGTACGTACGT"
    pos2 = 1_000_005
    seqB = "ACGTACTTACGT"
    _ = insertion_overlap_with_sequence(
        pos1, pos1, seqA,
        pos2, pos2, seqB,
        distance=50,
        max_hamming=0.2
    )


def benchmark():
    import time
    import random

    def random_dna(length):
        return ''.join(random.choice('ACGT') for _ in range(length))

    lengths = [10, 50, 100, 500, 1000]
    for length in lengths:
        seq1 = random_dna(length)
        seq2 = random_dna(length)
        start = time.time()
        for _ in range(1000):
            hamming_distance(seq1, seq2)
        elapsed = time.time() - start
        print(f"hamming len={length}: {elapsed:.4f}s for 1000 comps")

    length = 100
    n = 100
    sequences = [random_dna(length) for _ in range(n)]
    start = time.time()
    for i in range(n):
        for j in range(i + 1, n):
            insertion_overlap_with_sequence(
                1000, 1000, sequences[i],
                1010, 1010, sequences[j],
                distance=50, max_hamming=0.2
            )
    elapsed = time.time() - start
    total = n * (n - 1) // 2
    print(f"overlap n={n}: {elapsed:.4f}s for {total} pairs")


if __name__ == "__main__":
    import sys

    if len(sys.argv) > 1 and sys.argv[1] == "--benchmark":
        benchmark()
    elif len(sys.argv) > 1 and sys.argv[1] == "--examples":
        usage_examples()
    else:
        run_tests()