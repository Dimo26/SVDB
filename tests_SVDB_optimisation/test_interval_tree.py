#!/usr/bin/env python3

import numpy as np
import sys
sys.path.append('/Users/dee/Desktop/SVDB')

from svdb.interval_tree_overlap import interval_tree_cluster


def test_interval_tree_overlapping_deletions():
    """Test 1: Deletions that overlap"""
    coordinates = np.array([
        [1000, 2000],   # DEL 1
        [1500, 2500],   # DEL 2 (overlaps DEL 1)
        [2000, 3000],   # DEL 3 (overlaps DEL 2)
        [10000, 11000], # DEL 4 (no overlap)
    ])

    max_distance = 1000
    labels = interval_tree_cluster(coordinates, max_distance=max_distance)

    expected_clusters = 2
    actual_clusters = len(set(labels))
    return actual_clusters == expected_clusters


def test_interval_tree_non_overlapping():
    """Test 2: Deletions that don't overlap but are close"""
    coordinates = np.array([
        [1000, 2000],   # DEL 1
        [2100, 3100],   # DEL 2 (100bp gap after DEL 1)
        [3200, 4200],   # DEL 3 (100bp gap after DEL 2)
    ])

    # Small max_distance: should not connect (expect 3 clusters)
    max_distance_small = 50
    labels_small = interval_tree_cluster(coordinates, max_distance=max_distance_small)
    n_clusters_small = len(set(labels_small))
    ok_small = (n_clusters_small == 3)

    # Large max_distance: should connect all (expect 1 cluster)
    max_distance_large = 300
    labels_large = interval_tree_cluster(coordinates, max_distance=max_distance_large)
    n_clusters_large = len(set(labels_large))
    ok_large = (n_clusters_large == 1)

    return ok_small and ok_large


def test_interval_tree_insertions():
    """Test 3: Insertions treated as single points"""
    coordinates = np.array([
        [5000, 5000],
        [5050, 5050],
        [5100, 5100],
        [10000, 10000],
        [10050, 10050],
    ])

    max_distance = 150
    labels = interval_tree_cluster(coordinates, max_distance=max_distance)

    expected_clusters = 2
    actual_clusters = len(set(labels))
    return actual_clusters == expected_clusters


def test_interval_tree_different_sizes():
    """Test 4: SVs of very different sizes"""
    coordinates = np.array([
        [1000, 1100],   # 100bp
        [900, 1200],    # 300bp (overlaps small)
        [1150, 1250],   # 100bp (overlaps large)
    ])

    max_distance = 100
    labels = interval_tree_cluster(coordinates, max_distance=max_distance)

    expected_clusters = 1
    actual_clusters = len(set(labels))
    return actual_clusters == expected_clusters


def test_interval_tree_chain_connection():
    """Test 5: Chain of connections"""
    coordinates = np.array([
        [1000, 1500],   # SV1
        [1400, 1900],   # SV2 (connects to SV1)
        [1800, 2300],   # SV3 (connects to SV2)
        [2200, 2700],   # SV4 (connects to SV3)
    ])

    max_distance = 200
    labels = interval_tree_cluster(coordinates, max_distance=max_distance)

    expected_clusters = 1
    actual_clusters = len(set(labels))
    return actual_clusters == expected_clusters


def test_interval_tree_inversions():
    """Test 6: Inversions cluster like deletions"""
    coordinates = np.array([
        [1000, 3000],
        [1200, 3200],
        [1400, 3400],
    ])

    max_distance = 500
    labels = interval_tree_cluster(coordinates, max_distance=max_distance)

    expected_clusters = 1
    actual_clusters = len(set(labels))
    return actual_clusters == expected_clusters


def run_all_interval_tree_tests():
    """Run all Interval Tree tests and return (name, passed) tuples"""
    tests = [
        ("Overlapping Deletions", test_interval_tree_overlapping_deletions),
        ("Non-overlapping Close", test_interval_tree_non_overlapping),
        ("Insertions", test_interval_tree_insertions),
        ("Different Sizes", test_interval_tree_different_sizes),
        ("Chain Connection", test_interval_tree_chain_connection),
        ("Inversions", test_interval_tree_inversions),
    ]

    results = []
    for test_name, test_func in tests:
        try:
            passed = test_func()
        except Exception:
            passed = False
        results.append((test_name, passed))
    return results


if __name__ == "__main__":
    results = run_all_interval_tree_tests()
    total = len(results)
    passed = sum(1 for _, ok in results if ok)
    print(f"Interval Tree tests: {passed}/{total} passed")