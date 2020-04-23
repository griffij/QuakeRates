"""Calculate number of clusters by use of the 
Gap statistic. Uses: https://github.com/milesgranger/gap_statistic
and based on their Example.ipynb.
"""

import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from gap_statistic import OptimalK
from sklearn.datasets.samples_generator import make_blobs
from sklearn.cluster import KMeans

# Initialise OptimalK class
#optimalK = OptimalK(parallel_backend='rust')
optimalK = OptimalK() 
optimalK

# Make some test data 
#X, y = make_blobs(n_samples=int(1e5), n_features=2, centers=3, random_state=25)
#print('Data shape: ', X.shape)
#print(X, type(X))
#X = np.array([[100., 1.], [200.,1.],[220.,1.],[230.,1.], [500.,1.], [600.,1.]])
X = np.array([[100.],[200.],[220.],[230.], [580.], [600.]]) 
#X = np.array([[100.],[200.],[300.],[400.], [500.], [600.]])
#X = np.array([[100.],[180.],[300.],[410.], [500.], [610.]])
print(X, type(X)) 
# Call OptimalK to determine best number of clusters
print('Calculating optimal number of clusters')
n_clusters = optimalK(X, cluster_array=np.arange(1, 6), n_refs=100)
print('Optimal clusters: ', n_clusters)
print('Diff', optimalK.gap_df["diff"])
#sys.exit()
optimalK.plot_results()
# Plot some results
plt.plot(optimalK.gap_df.n_clusters, optimalK.gap_df.gap_value, linewidth=3)
plt.scatter(optimalK.gap_df[optimalK.gap_df.n_clusters == n_clusters].n_clusters,
            optimalK.gap_df[optimalK.gap_df.n_clusters == n_clusters].gap_value, s=250, c='r')
plt.grid(True)
plt.xlabel('Cluster Count')
plt.ylabel('Gap Value')
plt.title('Gap Values by Cluster Count')
plt.show()

# Now that we have the optimal clusters, n, we build our own KMeans model...
km = KMeans(n_clusters)
km.fit(X)

print('Cluster centres', km.cluster_centers_)
print('Cluster labels', km.labels_)
all_ie_times = []
for i, cc in enumerate(km.cluster_centers_):
    indices = np.argwhere(km.labels_ == i)
    cluster_events = X[indices].flatten()
    print(cluster_events)
    if len(cluster_events) > 1:
        print(np.sort(cluster_events))
        interevent_times = np.diff(np.sort(cluster_events), axis=0)
        for ie_t in interevent_times:
            all_ie_times.append(ie_t)
all_ie_times = np.array(all_ie_times)
print('all_ie_times', all_ie_times)
mean_within_cluster_ie_time = np.mean(all_ie_times)
print('mean_within_cluster_ie_time', mean_within_cluster_ie_time)

# Now calculate mean time between cluster centres
cc = np.sort(km.cluster_centers_.flatten())
print(cc)
cluster_ie_times = np.diff(cc)
mean_cluster_ie_time = np.mean(cluster_ie_times)
print('cluster_ie_times', cluster_ie_times)
print('mean_cluster_ie_time', mean_cluster_ie_time)

#df = pd.DataFrame(X, columns=['x','y'])
df = pd.DataFrame(X, columns=['x'])
df['label'] = km.labels_
data_ones = np.ones(len(X))
colors = plt.cm.Spectral(np.linspace(0, 1, len(df.label.unique())))

#for color, label in zip(colors, df.label.unique()):
    
#    tempdf = df[df.label == label]
#    tempdf.y=data_ones # Add some dummy data for plotting
#    plt.scatter(tempdf.x, tempdf.y, c=color)
plt.scatter(df.x, data_ones)
#plt.show()
cluster_ones = np.ones(len(km.cluster_centers_))
#plt.scatter(km.cluster_centers_[:,0], km.cluster_centers_[:, 1], c='r', s=500, alpha=0.7, )
plt.scatter(km.cluster_centers_, cluster_ones, c='r', s=500, alpha=0.7, )
plt.grid(True)
plt.show()

