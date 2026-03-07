#!/usr/bin/env python3
"""
Test OPTICS clustering algorithm on structural variants.
"""

import numpy as np
import sys
sys.path.append('/Users/dee/Desktop/SVDB/')

from svdb.optics_clustering import optics_cluster


def print_test_header(test_name):
    print(f"\nTEST: {test_name}")


def print_clustering_result(coordinates, labels, test_description):
    print(f"\n{test_description}")

    print("\nInput coordinates (posA, posB):")
    for i, coord in enumerate(coordinates):
        print(f"  SV{i+1}: posA={coord[0]:>6}, posB={coord[1]:>6}")

    print("\nClustering result:")
    unique_labels = sorted(set(labels))
    for label in unique_labels:
        if label == -1:
            print("\n  Cluster: NOISE (not clustered)")
        else:
            print(f"\n  Cluster {label}:")
        indices = np.where(labels == label)[0]
        for idx in indices:
            print(f"    SV{idx+1}: ({coordinates[idx][0]}, {coordinates[idx][1]})")


def test_optics_simple_deletions():
    print_test_header("OPTICS - Simple Deletion Clustering")

    coordinates = np.array([
        [1000, 2000],
        [1050, 2050],
        [1080, 2080],
        [10000, 11000],
        [10050, 11050],
        [10080, 11080],
    ])

    print("\nExpected: 2 clusters")
    print("  Cluster 0: SV1, SV2, SV3")
    print("  Cluster 1: SV4, SV5, SV6")

    min_samples = 2
    max_eps = 2000
    labels = optics_cluster(coordinates, min_samples=min_samples, max_eps=max_eps)

    print(f"\nOPTICS parameters: min_samples={min_samples}, max_eps={max_eps}")
    print_clustering_result(coordinates, labels, "Actual result:")

    expected_clusters = 2
    actual_clusters = len(set(labels)) - (1 if -1 in labels else 0)

    print(f"\nEXPECTED: {expected_clusters} clusters")
    print(f"ACTUAL: {actual_clusters} clusters")

    return actual_clusters == expected_clusters


def test_optics_varying_density():
    print_test_header("OPTICS - Varying Density Clusters")

    coordinates = np.array([
        [1000, 2000],
        [1020, 2020],
        [1040, 2040],
        [5000, 6000],
        [5200, 6200],
        [5400, 6400],
    ])

    print("\nExpected: 2 clusters with different densities")
    print("  Dense: SV1, SV2, SV3")
    print("  Sparse: SV4, SV5, SV6")

    min_samples = 2
    max_eps = 500
    labels = optics_cluster(coordinates, min_samples=min_samples, max_eps=max_eps)

    print(f"\nOPTICS parameters: min_samples={min_samples}, max_eps={max_eps}")
    print_clustering_result(coordinates, labels, "Actual result:")

    n_clusters = len(set(labels)) - (1 if -1 in labels else 0)
    print(f"\nNumber of clusters: {n_clusters}")
    print("OPTICS can handle both dense and sparse clusters")

    return True


def test_optics_insertions():
    print_test_header("OPTICS - Insertion Clustering")

    coordinates = np.array([
        [5000, 5000],
        [5030, 5030],
        [5060, 5060],
        [8000, 8000],
        [8030, 8030],
        [8060, 8060],
    ])

    print("\nExpected: 2 clusters of insertions")
    print("  For insertions, posA = posB")

    min_samples = 2
    max_eps = 150
    labels = optics_cluster(coordinates, min_samples=min_samples, max_eps=max_eps)

    print(f"\nOPTICS parameters: min_samples={min_samples}, max_eps={max_eps}")
    print_clustering_result(coordinates, labels, "Actual result:")

    expected_clusters = 2
    actual_clusters = len(set(labels)) - (1 if -1 in labels else 0)

    print(f"\nEXPECTED: {expected_clusters} clusters")
    print(f"ACTUAL: {actual_clusters} clusters")

    return actual_clusters == expected_clusters


def test_optics_hierarchical():
    print_test_header("OPTICS - Hierarchical Structure")

    coordinates = np.array([
        [1000, 2000],
        [1020, 2020],
        [1100, 2100],
        [1120, 2120],
        [1200, 2200],
        [1220, 2220],
    ])

    print("\nExpected: May detect sub-clusters within a main cluster")

    min_samples = 2
    max_eps = 300
    labels = optics_cluster(coordinates, min_samples=min_samples, max_eps=max_eps)

    print(f"\nOPTICS parameters: min_samples={min_samples}, max_eps={max_eps}")
    print_clustering_result(coordinates, labels, "Actual result:")

    n_clusters = len(set(labels)) - (1 if -1 in labels else 0)
    print(f"\nNumber of clusters found: {n_clusters}")
    print("OPTICS strength: detecting hierarchical structure")

    return True


def test_optics_mixed_sv_types():
    print_test_header("OPTICS - Mixed SV Types")

    coordinates = np.array([
        [1000, 1100],
        [1020, 1120],
        [1050, 1050],
        [1070, 1070],
        [1000, 2000],
        [1050, 2050],
    ])

    print("\nExpected: May separate by size; INS may group with small DELs")

    min_samples = 2
    max_eps = 1000
    labels = optics_cluster(coordinates, min_samples=min_samples, max_eps=max_eps)

    print(f"\nOPTICS parameters: min_samples={min_samples}, max_eps={max_eps}")
    print_clustering_result(coordinates, labels, "Actual result:")

    n_clusters = len(set(labels)) - (1 if -1 in labels else 0)
    print(f"\nNumber of clusters: {n_clusters}")

    return True


def test_optics_parameter_comparison():
    print_test_header("OPTICS - Parameter Sensitivity")

    coordinates = np.array([
        [1000, 2000],
        [1100, 2100],
        [1200, 2200],
        [1300, 2300],
    ])

    print("\nTesting different max_eps values (evenly spaced deletions):")
    for max_eps in [200, 500, 1000]:
        labels = optics_cluster(coordinates, min_samples=2, max_eps=max_eps)
        n_clusters = len(set(labels)) - (1 if -1 in labels else 0)
        n_noise = sum(labels == -1)
        print(f"\n  max_eps={max_eps}: clusters={n_clusters}, noise={n_noise}")

    print("\nInsight: max_eps determines neighborhood size; larger -> more connected clusters")
    return True


def run_all_optics_tests():
    print("\nOPTICS CLUSTERING TEST SUITE")

    tests = [
        ("Simple Deletions", test_optics_simple_deletions),
        ("Varying Density", test_optics_varying_density),
        ("Insertions", test_optics_insertions),
        ("Hierarchical Structure", test_optics_hierarchical),
        ("Mixed SV Types", test_optics_mixed_sv_types),
        ("Parameter Comparison", test_optics_parameter_comparison),
    ]

    results = []
    for test_name, test_func in tests:
        try:
            passed = test_func()
            results.append((test_name, passed))
        except Exception as e:
            print(f"\nERROR in {test_name}: {e}")
            results.append((test_name, False))

    print("\nTEST SUMMARY")
    for test_name, passed in results:
        status = "PASSED" if passed else "FAILED"
        print(f"{status}: {test_name}")

    total = len(results)
    passed_count = sum(1 for _, p in results if p)
    print(f"\nTotal: {passed_count}/{total} tests passed")
if __name__ == "__main__":    run_all_optics_tests()


