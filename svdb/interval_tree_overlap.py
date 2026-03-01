import numpy as np
from collections import namedtuple

# Named tuple for storing interval with metadata
Interval = namedtuple('Interval', ['start', 'end', 'index', 'data'])


class IntervalNode:

    def __init__(self, intervals, depth=0, max_depth=20):

        self.intervals_by_start = []  # Center intervals sorted by start
        self.intervals_by_end = []    # Center intervals sorted by end
        self.center = None
        self.left = None
        self.right = None
        
        if not intervals or depth >= max_depth:
            # Leaf node: store intervals sorted by start and end
            self.intervals_by_start = sorted(intervals, key=lambda x: x.start)
            self.intervals_by_end = sorted(intervals, key=lambda x: x.end)
            return

        # compute center as median of all interval endpoints
        all_points = []
        for interval in intervals:
            all_points.extend([interval.start, interval.end])
        self.center = float(np.median(all_points))

        left_intervals = []
        center_intervals = []
        right_intervals = []

        for interval in intervals:
            if interval.end < self.center:
                left_intervals.append(interval)
            elif interval.start > self.center:
                right_intervals.append(interval)
            else:
                center_intervals.append(interval)

        self.intervals_by_start = sorted(center_intervals, key=lambda x: x.start)
        self.intervals_by_end = sorted(center_intervals, key=lambda x: x.end)

        if left_intervals:
            self.left = IntervalNode(left_intervals, depth + 1, max_depth)
        if right_intervals:
            self.right = IntervalNode(right_intervals, depth + 1, max_depth)
    
    def query(self, start, end):
        results = []
        # leaf node: center is None  scan stored intervals_by_start with early stop
        if self.center is None:
            for interval in self.intervals_by_start:
                if interval.start > end:
                    break
                if interval.end >= start:
                    results.append(interval)
            return results

        # check center intervals (sorted by start) with early stopping
        for interval in self.intervals_by_start:
            if interval.start > end:
                break
            if self.intervals_overlap(start, end, interval.start, interval.end):
                results.append(interval)

        # recurse left/right as needed
        if self.left and start <= self.center:
            results.extend(self.left.query(start, end))
        if self.right and end >= self.center:
            results.extend(self.right.query(start, end))

        return results


    def intervals_overlap(self, start1, end1, start2, end2):
        return start1 <= end2 and start2 <= end1
    

class IntervalTree:
    
    def __init__(self):
        self.intervals = []
        self.root = None
    
    def add(self, start, end, index=None, data=None):
        interval = Interval(start=start, end=end, index = index, data=data)
        self.intervals.append(interval)
    
    def build(self):
        if self.intervals:
             self.root = IntervalNode(self.intervals)

    def query(self, start, end):
        if self.root is None:
             return []  
        return self.root.query(start,end)

    
    def query_indices(self, start, end):

        results = self.query(start, end)
        return [interval.index for interval in results]

class SVIntervalTree:

    def __init__(self):
        self.trees = {}
        self.variants = []

 
    def add_variant(self, chrA, posA, chrB, posB, variant_type, index=None, **metadata):
        key = (chrA, chrB)
        if key not in self.trees:
            self.trees[key] = IntervalTree()

        data = {'chrA': chrA, 'posA': posA, 'chrB': chrB, 'posB': posB,
            'variant_type': variant_type, **metadata}
        self.trees[key].add(posA, posB, index=index, data=data)
        self.variants.append(data)

    
    def build(self):
        for tree in self.trees.values():
             tree.build()

    
    def query_overlaps(self, chrA, posA, chrB, posB, distance=500):
        key = (chrA, chrB)
        if key not in self.trees:
            return []
        results = self.trees[key].query(posA-distance, posB+distance)
        overlapping = []
        for interval in results:
             var_data = interval.data
             if (abs(posA - var_data['posA']) <= distance and abs(posB - var_data['posB']) <= distance):
                  overlapping.append(var_data)
        return overlapping


def interval_tree_cluster(coordinates, max_distance=500):

    n = len(coordinates)
    tree = IntervalTree()
    for i, (posA, _) in enumerate(coordinates):
        tree.add(posA, posA, index=i)  # point interval on posA only
    tree.build()

    adjacency = {i: set() for i in range(n)}

    for i, (posA, posB) in enumerate(coordinates):
        candidates = tree.query(posA - max_distance, posA + max_distance)
        for interval in candidates:
            j = interval.index
            if i == j:
                continue
            if abs(coordinates[j][1] - posB) <= max_distance:
                adjacency[i].add(j)
                adjacency[j].add(i)

    labels = np.full(n, -1)
    current_cluster = 0
    visited = set()

    # Iterative DFS using stack to avoid recursion depth issues
    for i in range(n):
        if i not in visited:
            stack = [i]
            visited.add(i)
            labels[i] = current_cluster
            
            while stack:
                node = stack.pop()
                for neighbor in adjacency[node]:
                    if neighbor not in visited:
                        visited.add(neighbor)
                        labels[neighbor] = current_cluster
                        stack.append(neighbor)
            
            current_cluster += 1

    return labels