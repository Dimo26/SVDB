
import numpy as np
import matplotlib.pyplot as plt

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
        for neighbor_idx in neighbors: #distance calc for all neighbors
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
    
    def _extract_clusters(self, reachability, ordering, threshold=None):
        n_samples = len(reachability)
        labels = np.full(n_samples, -1)  # noise
         
        if threshold is None:
            valid_reach = reachability[reachability < np.inf]
            if len(valid_reach) > 0:
                threshold = np.percentile(valid_reach, 75)
            else:
                return labels
         
        current_cluster = 0
        for i, point_idx in enumerate(ordering):
            reach_dist = reachability[point_idx]
             
            if reach_dist <= threshold:
                labels[point_idx] = current_cluster
            else:
                if i > 0 and reachability[ordering[i-1]] <= threshold:
                    current_cluster += 1
        
        return labels
    
    def fit(self, data):
      
        self.ordering, self.reachability_, self.core_distances_ = self._optics_algorithm(data)
        self.labels_ = self._extract_clusters(self.reachability_, self.ordering)
        return self

    def fit_predict(self, data):

        self.fit(data)
        return self.labels_
        
def optics_cluster(coordinates, min_samples=2, max_eps=2000):

    clusterer = OPTICS(min_samples=min_samples, max_eps=max_eps)
    labels = clusterer.fit_predict(coordinates)
    return labels

def test_optics():
    cluster1 = np.random.normal(loc=(1000,1000), scale=100, size=(30,2))
    cluster2 = np.random.normal(loc=(5000,5000), scale=100, size=(30,2))
    cluster3 = np.random.normal(loc=(9000,1000), scale=100, size=(30,2))
    noise = np.random.uniform([0,0], [10000,10000], size=(10,2))
    data = np.vstack([cluster1, cluster2, cluster3, noise])
    labels = optics_cluster(data)
    plt.scatter(data[:,0], data[:,1], c=labels)
    plt.title('OPTICS Clustering Result')
    plt.xlabel('Position A') 
    plt.ylabel('Position B')
    plt.colorbar(label='Clusters')
    plt.show()

if __name__ == "__main__":
    test_optics()
