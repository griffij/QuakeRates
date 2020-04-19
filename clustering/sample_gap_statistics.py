"""Calculate number of clusters by use of the 
Gap statistic.
Generate n chronology samples of each paleo-earthquake record and
then estimate the number of clusters in the record by maximising the 
Gap statistic. 

Gap statistic is from Tibshirani, R., Walther, G., and Hastie, T. (2001). 
Estimating the numbers of clusters in a data set via the gap statistic. 
J. R. Statist. Soc. B, 63(2): 411-423.

This code uses Miles Granger's implementation at: 
https://github.com/milesgranger/gap_statistic
and is based on their Example.ipynb.    
"""

import os, sys
from glob import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from gap_statistic import OptimalK
from sklearn.datasets.samples_generator import make_blobs
from sklearn.cluster import KMeans
from QuakeRates.dataman.parse_params import parse_param_file, \
    get_event_sets

# Define parameter files
filepath = '../params'
param_file_list = glob(os.path.join(filepath, 'A*.txt'))
n_samples = 500  # Number of Monte Carlo samples of the eq chronologies
half_n = int(n_samples/2)

plot_dir = './plots'
if not os.path.exists(plot_dir):
    os.makedirs(plot_dir)
# Define subset to take
#faulting_styles = ['Reverse']
#faulting_styles = ['Normal']
#faulting_styles = ['Strike_slip']
faulting_styles = ['all']
tectonic_regions = ['all']
#tectonic_regions = ['Plate_boundary_master', 'Plate_boundary_network']
min_number_events = 1

names, event_sets, event_certainties, num_events = \
    get_event_sets(param_file_list, tectonic_regions,
                   faulting_styles, min_number_events)  

# Now loop over paleo-earthquake records
for i, event_set in enumerate(event_sets):
    # Generate some chronologies
    event_set.gen_chronologies(n_samples, observation_end=2019, min_separation=1)
    print(num_events[i])
    optimal_ks = []
    for chron in event_set.chronology.T:
     #   print(chron)
        # Initialise OptimalK
        #optimalK = OptimalK(parallel_backend='rust')
        optimalK = OptimalK() 
        optimalK
        # Call OptimalK to determine best number of clusters
#        print('Calculating optimal number of clusters')
        cluster_array = np.arange(1, num_events[i])
        try:
            n_clusters = optimalK(chron, cluster_array=cluster_array, n_refs=50)
        except UserWarning:
            print('Re-running as did not get sensible results')
            # Try re-running once (will have different random initialisation
            # if get odd results
            try:
                n_clusters = optimalK(chron, cluster_array=cluster_array, n_refs=50)
            except UserWarning:
                msg = 'Could not run fault %s, continuing to next fault' % names[i]
                print(msg)
                continue
            
#        print('Optimal clusters: ', n_clusters)
        optimal_ks.append(n_clusters)
    print(names[i])
    print(optimal_ks)
    # Plot histogram of optimal k values
    plt.clf()
    bins = np.arange(1, num_events[i]+1) - 0.5
    plt.hist(optimal_ks, bins=bins,
             facecolor='0.6', edgecolor='0.2', density=True)
    plt.xlim([0, num_events[i]])
    figname = 'OptimalK_histogram_%s.png' % names[i]
    fig_filename = os.path.join(plot_dir, figname)
    plt.savefig(fig_filename)
sys.exit()
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

print(km.cluster_centers_)
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

