import numpy as np

class OPTICS:
    def __init__(self, min_samples=2, max_eps=2000, metric='euclidean'):
        self.min_samples = min_samples
        self.max_eps = max_eps
        self.metric = metric
        
        self.labels_ = None  # Cluster labels for each point
        self.reachability_ = None  # Reachability distance for each point
        self.ordering_ = None  # Order of points processed
        self.core_distances_ = None  # Core distance for each point
    
    def _calculate_core_distance(self, point_idx, data, neighbors):
        
        if len(neighbors) < self.min_samples:
            return np.inf

        distances = []
        for neighbor_idx in neighbors: #distance calc for all neigborinos
            dist = np.linalg.norm(data[point_idx] - data[neighbor_idx])
            distances.append(dist)

        distances = np.sort(distances)
        return distances[self.min_samples - 1]

    def _get_neighbors(self, point_idx, data):

        distances = np.linalg.norm(data - data[point_idx], axis=1)
        neighbors = np.where((distances <= self.max_eps) & (distances > 0))[0]
        return neighbors
        
    def _calculate_reachability_distance(self, point1_idx, point2_idx, data, core_dist1):
        
        actual_distance = np.linalg.norm(data[point1_idx] - data[point2_idx])
        return max(core_dist1, actual_distance)

    def _optics_algorithm(self, data):

        n_samples = len(data)
        processed = np.zeros(n_samples, dtype=bool)
        reachability = np.full(n_samples, np.inf)
        core_distances = np.full(n_samples, np.inf)
        ordering = [] # list for order of processed points
        seeds = [] #list for reachability and point indices
        for point_idx in range(n_samples):
              if processed[point_idx]:
                  continue

              neighbors = self._get_neighbors(point_idx, data)
              processed[point_idx] = True
              ordering.append(point_idx)

              core_dist = self._calculate_core_distance(point_idx, data, neighbors)
              core_distances[point_idx] = core_dist

              if core_dist != np.inf:
                  for neighbor_idx in neighbors:
                      if not processed[neighbor_idx]:
                          reach_dist = self._calculate_reachability_distance(
                              point_idx, neighbor_idx, data, core_dist
                          )

                          if reach_dist < reachability[neighbor_idx]:
                              reachability[neighbor_idx] = reach_dist

                              seeds.append((reach_dist, neighbor_idx))

              seeds.sort()  # Sort by reachability distance
              while seeds:
                  _, seed_idx = seeds.pop(0)
                  if processed[seed_idx]:
                      continue

        return ordering, reachability, core_distances
    
    def _extract_clusters(self, reachability, ordering, threshold=None, min_cluster_size=None):
        n_samples = len(reachability)
        labels = np.full(n_samples, -1)  # noise
        
        if min_cluster_size is None:
            min_cluster_size = self.min_samples
         
        if threshold is None:
            valid_reach = reachability[reachability < np.inf]
            if len(valid_reach) > 0:
                # Use adaptive threshold based on data characteristics
                q25 = np.percentile(valid_reach, 25)
                q75 = np.percentile(valid_reach, 75)
                q90 = np.percentile(valid_reach, 90)
                median = np.median(valid_reach)
                max_reach = np.max(valid_reach)
                
                iqr = q75 - q25
                if iqr < median * 0.6:  # dense data
                    # Use 90th percentile or max_eps, whichever is smaller
                    threshold = min(q90 * 1.2, self.max_eps)
                else:  # Mixed density
                    threshold = q75
            else:
                return labels
         
        current_cluster = 0
        in_cluster = False
        cluster_start = 0
        
        for i, point_idx in enumerate(ordering):
            reach_dist = reachability[point_idx]
             
            if reach_dist <= threshold:
                if not in_cluster:
                    # Starting a new cluster
                    cluster_start = i
                    in_cluster = True
                labels[point_idx] = current_cluster
            else:
                # Reachability exceeds threshold
                if in_cluster:
                    # Check if cluster is large enough
                    cluster_size = i - cluster_start
                    if cluster_size < min_cluster_size:
                        # Mark small cluster as noise
                        for j in range(cluster_start, i):
                            labels[ordering[j]] = -1
                    # Move to next cluster
                    current_cluster += 1
                    in_cluster = False
        
        if in_cluster:
            cluster_size = len(ordering) - cluster_start
            if cluster_size < min_cluster_size:
                for j in range(cluster_start, len(ordering)):
                    labels[ordering[j]] = -1
        
        return labels
    
    def fit(self, data):
      
        self.ordering, self.reachability_, self.core_distances_ = self._optics_algorithm(data)
        self.labels_ = self._extract_clusters(self.reachability_, self.ordering)
        return self

    def fit_predict(self, data):

        self.fit(data)
        return self.labels_
        
def optics_cluster(coordinates, min_samples=2, max_eps=1000):

    clusterer = OPTICS(min_samples=min_samples, max_eps=max_eps)
    labels = clusterer.fit_predict(coordinates)
    return labels
