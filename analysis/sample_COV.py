"""Script forA sampling COV uncertainties on many faults and plotting them
"""

import os
import ast
from glob import glob
from operator import itemgetter 
import numpy as np
from matplotlib import pyplot
from matplotlib.patches import PathPatch
from scipy.stats import kde
from QuakeRates.dataman.event_dates import EventSet
from QuakeRates.dataman.parse_oxcal import parse_oxcal
from QuakeRates.dataman.parse_age_sigma import parse_age_sigma

filepath = '../params'
param_file_list = glob(os.path.join(filepath, 'Yam*.txt'))
n_samples = 1000  # Number of Monte Carlo samples of the eq chronologies


params = {}


def parse_param_file(param_file_name):
    params={}
    with open(param_file_name) as f_in:
        for line in f_in:
            var, value = line.strip().split('=')
            params[var.strip()] = ast.literal_eval(value.strip())
    return params

covs = []
long_term_rates = []
names = []
for param_file in param_file_list:
    name = param_file.split('/')[-1].split('_')[0]
    print(name)
    names.append(name)
    params = parse_param_file(param_file)
    print(params)
    # Deal with OxCal output and lists of dates with uncertainties
    # separately
    # Check that data file exists
    if not os.path.isfile(params['filename']):
        msg = 'Data file ' + params['filename'] + ' does not exist.' + \
            ' Continuing to next file.'
        print(msg)
        continue
    try:
        params['chron_type']
    except KeyError:
        msg = 'chron_type not defined in parameter file ' + param_file
        print(msg)
        raise
    if params['chron_type'] == 'OxCal':
        events = parse_oxcal(params['filename'], params['events'],
                             params['event_order'])
        event_set = EventSet(events)  
    elif params['chron_type'] == 'Age2Sigma':
        # Write method to sample these
        events, event_certainty = parse_age_sigma(params['filename'],
                                                  params['sigma_level'],
                                                  params['event_order'])
        event_set = EventSet(events)
    else:
        msg = 'Unknown form of chron_type defined in ' + param_file
        raise Exception(msg)
    # Handle cases with uncertain number of events. Where events identification is
    # unsure, event_certainty is given a value of 0, compared with 1 for certain
    # events
    # First generate chronologies assuming all events are certain
    event_set.gen_chronologies(n_samples, observation_end=2019, min_separation=1)
    event_set.calculate_cov()
    event_set.cov_density() 
    # Now generate chronologies assuming uncertain events did not occur
    if sum(event_certainty) < len(events):
        indices = np.where(event_certainty == 1)
        indices = list(indices[0])
        print(indices[0], type(indices))
        events_subset = list(itemgetter(*indices)(events)) 
        event_set_certain = EventSet(events_subset)
        event_set_certain.gen_chronologies(n_samples, observation_end=2019, min_separation=1)
        event_set_certain.calculate_cov()
        event_set_certain.cov_density()
        combined_covs = np.concatenate([event_set.covs, event_set_certain.covs])
        covs.append(combined_covs)
        combined_ltrs = np.concatenate([event_set.long_term_rates,
                                        event_set_certain.long_term_rates])
        long_term_rates.append(combined_ltrs)
    else:
        covs.append(event_set.covs)
        long_term_rates.append(event_set.long_term_rates)

# Now do some plotting
pyplot.clf()
ax = pyplot.subplot(111)
nbins = 100
for i, cov_set in enumerate(covs):
    print('Plotting ' + names[i])
    if names[i]=='NankaiTrough': # Deal with special case later
        print('Skipping as single data point')
        continue
    cov_samples = np.array([cov_set, long_term_rates[i]])                                                                               
    x, y = cov_samples
    # Evaluate a gaussian kde on a regular grid of nbins x nbins over
    # data extents
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
        try:
            c_path = c.get_paths()[0]
            patch = PathPatch(c_path, transform=ax.transData, facecolor='none',
                              linewidth=0.4, linestyle='--')
        except:
            print('Cannot plot this contour, does not exist')
##        print('Adding patch')
##        ax.add_patch(patch)
    m = ax.pcolormesh(xi, yi, zi_norm.reshape(xi.shape), shading='gouraud',
                      cmap=pyplot.cm.Greys, clip_path=patch)

#FIXME - make bounds parametric/automatic
ax.set_xlim([0, 2.5])
ax.set_ylim([1./1000000, 1./50])
ax.set_yscale('log')
ax.set_xlabel('COV')
ax.set_ylabel('Long-term rate (events per year)')
pyplot.savefig('cov_vs_lt_rate.png')
