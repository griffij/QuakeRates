"""Script forA sampling COV uncertainties on many faults and plotting them
"""

import os
import ast
from glob import glob
import numpy as np
from matplotlib import pyplot
from matplotlib.patches import PathPatch
from scipy.stats import kde
from QuakeRates.dataman.event_dates import EventSet
from QuakeRates.dataman.parse_oxcal import parse_oxcal

filepath = '../params'
param_file_list = glob(os.path.join(filepath, '*.txt'))
n_samples = 10000  # Number of Monte Carlo samples to take of the earthquake chronologies


params = {}


def parse_param_file(param_file_name):
    with open(param_file_name) as f_in:
        for line in f_in:
            var, value = line.strip().split('=')
            params[var.strip()] = ast.literal_eval(value.strip())
    return params

covs = []
long_term_rates = []
for param_file in param_file_list:
    params = parse_param_file(param_file)
    print(params)
    events = parse_oxcal(params['filename'], params['events'],
                         params['event_order'])
    event_set = EventSet(events)  
    event_set.gen_chronologies(n_samples)
    event_set.calculate_cov() 
    event_set.cov_density()
    covs.append(event_set.covs)
    long_term_rates.append(event_set.long_term_rates)

# Now do some plotting
pyplot.clf()
ax = pyplot.subplot(111)
nbins = 100
for i, cov_set in enumerate(covs):
    cov_samples = np.array([cov_set, long_term_rates[i]])                                                                               
    x, y = cov_samples
    # Evaluate a gaussian kde on a regular grid of nbins x nbins over data extents                                                          
    k = kde.gaussian_kde(cov_samples)
    xi, yi = np.mgrid[0.98*x.min():1.02*x.max():nbins*1j, 0.98*y.min():1.02*y.max():nbins*1j]
    zi = k(np.vstack([xi.flatten(), yi.flatten()]))
    # Calculate percentiles
    cumsum = np.cumsum(zi)
    # Normalise data and calculate bottom 5th percentile for drawing contour
    # that contains 95% of the distribution.
    zi_norm = zi/max(cumsum)
    perc_5th = 0.05*max(zi_norm)
    # Plot the data
    # Slightly extend the bounds of the data for smoother plotting
    
 #   m = ax.pcolormesh(xi, yi, zi_norm.reshape(xi.shape), shading='gouraud', cmap=pyplot.cm.BuGn)
    cs = ax.contour(xi, yi, zi_norm.reshape(xi.shape), [perc_5th], colors='k',
                    linewidths=0.4)
    for c in cs.collections:
        c_path = c.get_paths()[0]
        patch = PathPatch(c_path, transform=ax.transData, facecolor='none',
                          linewidth=0.4, linestyle='--')
##        print('Adding patch')
##        ax.add_patch(patch)
    m = ax.pcolormesh(xi, yi, zi_norm.reshape(xi.shape), shading='gouraud',
                      cmap=pyplot.cm.Greys, clip_path=patch)

#FIXME - make bounds parametric/automatic
ax.set_xlim([0, 2.5])
ax.set_ylim([1./100000, 1./100])
ax.set_yscale('log')
ax.set_xlabel('COV')
ax.set_ylabel('Long-term rate (events per year)')
pyplot.savefig('cov_vs_lt_rate.png')
