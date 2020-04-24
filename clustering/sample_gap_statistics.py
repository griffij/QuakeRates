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
param_file_list = glob(os.path.join(filepath, '*.txt'))
n_samples = 20  # Number of Monte Carlo samples of the eq chronologies
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
plot_colours = []

names, event_sets, event_certainties, num_events = \
    get_event_sets(param_file_list, tectonic_regions,
                   faulting_styles, min_number_events)  

# Now loop over paleo-earthquake records
within_clusters = []
between_clusters = []
long_term_rates = []
for i, event_set in enumerate(event_sets):
    print('Running %s fault' % names[i])
    # Generate some chronologies
    event_set.gen_chronologies(n_samples, observation_end=2019, min_separation=1)
    event_set.calculate_cov()
    event_set.cov_density()  
    event_set.basic_chronology_stats()
    long_term_rates.append(event_set.long_term_rates)
    print(num_events[i])
    optimal_ks = []
    mean_within_cluster_ie_times = []
    mean_between_cluster_times = []
    for chron in event_set.chronology.T:
        chron = np.expand_dims(chron, axis=1)
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
        # Now run k-means clustering algorithm with optimal k
        try:
            km = KMeans(n_clusters)
            km.fit(chron)
        except UserWarning:
            print('Re-running k-means algorithm')
            try:
                km.fit(chron)
            except UserWarning:
                msg = 'Did not fit clusters'
                print(msg)
                continue
        # Now get statistics of within cluster event times and time
        # between cluster centres
        if n_clusters > 1:
            all_ie_times = []
            cluster_bounds = []
            for j, cc in enumerate(km.cluster_centers_):
                indices = np.argwhere(km.labels_ == j)
                cluster_events = chron[indices].flatten()
                if len(cluster_events) > 1:
                    cluster_events = np.sort(cluster_events)
                    interevent_times = np.diff(cluster_events, axis=0)
                    for ie_t in interevent_times:
                        all_ie_times.append(ie_t)
                # Get start and end event of the cluster (which may be the same)
                cluster_bounds.append([cluster_events[0], cluster_events[-1]])
            all_ie_times = np.array(all_ie_times)
            mean_within_cluster_ie_time = np.mean(all_ie_times)
            cluster_bounds = np.array(cluster_bounds)
            
            # Update: calculate mean time bewteen cluster edge events
            #cluster_ind = np.argsort(km.cluster_centers_.flatten())
            cluster_bounds = np.sort(cluster_bounds, axis=0)
            cluster_ie_times = []
            for k, cb in enumerate(cluster_bounds):
                if k == 0:
                    pass
                else:
                    cluster_ie_time = cluster_bounds[k][0] - cluster_bounds[k-1][1]
                    cluster_ie_times.append(cluster_ie_time)
            cluster_ie_times = np.array(cluster_ie_times)
            mean_cluster_ie_time = np.mean(cluster_ie_times)
            mean_within_cluster_ie_times.append(mean_within_cluster_ie_time)
            mean_between_cluster_times.append(mean_cluster_ie_time)
        # Case where we only have one cluster (i.e. unclustered data)
        elif n_clusters == 1:
            ie_times = np.diff(chron.flatten())
            mean_ie_time = np.mean(ie_times)
            # In this case both the between and within cluster time are treated as equal
            mean_within_cluster_ie_times.append(mean_ie_time)
            mean_between_cluster_times.append(mean_ie_time)
            
        optimal_ks.append(n_clusters)
    within_clusters.append(mean_within_cluster_ie_times)
    between_clusters.append(mean_between_cluster_times)
    # Plot histogram of optimal k values
    plt.clf()
    bins = np.arange(1, num_events[i]+1) - 0.5
    plt.hist(optimal_ks, bins=bins,
             facecolor='0.6', edgecolor='0.2', density=True)
    plt.xlim([0, num_events[i]])
    figname = 'OptimalK_histogram_%s.png' % names[i]
    fig_filename = os.path.join(plot_dir, figname)
    plt.savefig(fig_filename)

    # Get colours for later plotting
    if event_set.faulting_style == 'Normal':
        plot_colours.append('r')
    elif event_set.faulting_style == 'Reverse':
        plot_colours.append('b')
    elif event_set.faulting_style == 'Strike_slip':
        plot_colours.append('g')
    else:
        plot_colours.append('k')
print(between_clusters)
print(within_clusters)

# Calculate means across samples
long_term_rates_T = np.array(long_term_rates).T 
mean_ltr = np.mean(long_term_rates_T, axis = 0)  
std_ltr = np.std(long_term_rates_T, axis = 0)
ltr_bounds = np.array([abs(mean_ltr - (np.percentile(long_term_rates_T, 2.5, axis=0))),
                       abs(mean_ltr - (np.percentile(long_term_rates_T, 97.5, axis=0)))])
within_clusters = np.array(within_clusters)
between_clusters = np.array(between_clusters)
ratio_within_between_clusters = np.mean(within_clusters.T, axis=0)/ \
    np.mean(between_clusters.T, axis=0)
print(ratio_within_between_clusters)

plt.clf()
ax = plt.subplot(111)

for i, wc in enumerate(within_clusters):
#    print(between_clusters[i])
#    print(wc)
    plt.scatter(np.mean(between_clusters[i]), np.mean(wc),
                marker = 's', c=plot_colours[i], s=25)
    ax.annotate(names[i][:4], (np.mean(between_clusters[i]), np.mean(wc)),
                fontsize=8)
# Add y=x lines
xvals = [1e2, 1e6]
yvals = xvals
plt.plot(xvals, yvals)
# Add y=0.1x line
xvals = [1e2, 1e7]
yvals = [1e1, 1e6]
plt.plot(xvals, yvals)
    
ax.set_xlabel('Mean time between clusters')
ax.set_ylabel('Mean within-cluster inter-event time')
ax.set_xscale('log')
ax.set_yscale('log')
plt.savefig('within_vs_between_cluster_time.png')

# Now plot as ratio against long-term rate 
plt.clf()
ax = plt.subplot(111)
plt.errorbar(mean_ltr, ratio_within_between_clusters,
             xerr = ltr_bounds,
             ecolor='0.4',
             linestyle="None")
plt.scatter(mean_ltr, ratio_within_between_clusters, marker = 's',
            c=plot_colours, s=25)
for i, name in enumerate(names):
    ax.annotate(name[:4], (mean_ltr[i], ratio_within_between_clusters[i]),
                fontsize=8)
ax.set_xlabel('Long-term rate')
ax.set_ylabel('Ratio of within cluster inter-event time to between cluster time')
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlim([1e-6, 1e-1])
plt.savefig('Ratio_within_btn_ie_time_ltr.png')

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

