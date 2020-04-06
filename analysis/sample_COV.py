"""Script for sampling COV uncertainties on many faults and plotting them
"""

import os
import ast
from glob import glob
from operator import itemgetter 
import numpy as np
from scipy.optimize import curve_fit
from matplotlib import pyplot
from matplotlib.patches import PathPatch
from scipy.stats import kde
from QuakeRates.dataman.event_dates import EventSet
from QuakeRates.dataman.parse_oxcal import parse_oxcal
from QuakeRates.dataman.parse_age_sigma import parse_age_sigma

filepath = '../params'
param_file_list = glob(os.path.join(filepath, '*.txt'))
n_samples = 10000  # Number of Monte Carlo samples of the eq chronologies
half_n = int(n_samples/2)
print(half_n)

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
max_interevent_times = []
min_interevent_times = []
min_paired_interevent_times = []
std_min_paired_interevent_times = []
std_min_interevent_times = []
std_max_interevent_times = []
max_interevent_times_bounds = []
min_interevent_times_bounds = []
min_paired_interevent_times_bounds = []
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
    # Now calculate some statistics on the sampled chronologies
    event_set.basic_chronology_stats()
    min_paired_interevent_times.append(event_set.mean_minimum_pair_interevent_time)
    max_interevent_times.append(event_set.mean_maximum_interevent_time)
    min_interevent_times.append(event_set.mean_minimum_interevent_time)  
    std_min_paired_interevent_times.append(event_set.std_minimum_pair_interevent_time)
    std_min_interevent_times.append(event_set.std_minimum_interevent_time)
    std_max_interevent_times.append(event_set.std_maximum_interevent_time)
    max_interevent_times_bounds.append([abs(event_set.mean_maximum_interevent_time -
                                            event_set.maximum_interevent_time_lb),
                                        abs(event_set.mean_maximum_interevent_time -
                                            event_set.maximum_interevent_time_ub)])
    min_interevent_times_bounds.append([abs(event_set.mean_minimum_interevent_time -
                                            event_set.minimum_interevent_time_lb),
                                        abs(event_set.mean_minimum_interevent_time -
                                            event_set.minimum_interevent_time_ub)])
    min_paired_interevent_times_bounds.append([abs(event_set.mean_minimum_pair_interevent_time -
                                                   event_set.minimum_pair_interevent_time_lb),
                                        abs(event_set.mean_minimum_pair_interevent_time -
                                            event_set.minimum_pair_interevent_time_ub)])
    # Now generate chronologies assuming uncertain events did not occur
    if sum(event_certainty) < len(events):
        indices = np.where(event_certainty == 1)
        indices = list(indices[0])
#        print(indices[0], type(indices))
        events_subset = list(itemgetter(*indices)(events)) 
        event_set_certain = EventSet(events_subset)
        event_set_certain.gen_chronologies(n_samples, observation_end=2019, min_separation=1)
        event_set_certain.calculate_cov()
        event_set_certain.cov_density()
        event_set.basic_chronology_stats()
        combined_covs = np.concatenate([event_set.covs[:half_n],
                                        event_set_certain.covs[:half_n]])
        covs.append(combined_covs)
        # Combine, taking n/2 samples from each set
        combined_ltrs = np.concatenate([event_set.long_term_rates[:half_n],
                                        event_set_certain.long_term_rates[:half_n]])
        print(len(combined_ltrs))
        long_term_rates.append(combined_ltrs)
    else:
        covs.append(event_set.covs)
        long_term_rates.append(event_set.long_term_rates)

max_interevent_times = np.array(max_interevent_times)
min_interevent_times = np.array(min_interevent_times)
min_paired_interevent_times = np.array(min_paired_interevent_times)
print('min_paired_interevent_times', min_paired_interevent_times)
max_interevent_times_bounds = np.array(max_interevent_times_bounds).T
min_interevent_times_bounds = np.array(min_interevent_times_bounds).T
min_paired_interevent_times_bounds = np.array(min_paired_interevent_times_bounds).T
long_term_rates_T = np.array(long_term_rates).T
mean_ltr = np.mean(long_term_rates_T, axis = 0)
print('mean_ltr', mean_ltr)
ltr_bounds = np.array([abs(mean_ltr - (np.percentile(long_term_rates_T, 2.5, axis=0))),
                       abs(mean_ltr - (np.percentile(long_term_rates_T, 97.5, axis=0)))])
print('ltr_bounds', ltr_bounds)
print(max_interevent_times_bounds, np.shape(max_interevent_times_bounds))
# Ratios 
#ratio_min_pair_max = min_paired_interevent_times/max_interevent_times
#print('ratio_min_pair_max', ratio_min_pair_max)
##ratio_min_pair_max_T = ratio_min_pair_max.T
#mean_ratio_min_pair_max = np.mean(ratio_min_pair_max, axis=0)
#print('mean_ratio_min_pair_max ', mean_ratio_min_pair_max)
#ratio_min_max = min_interevent_times/max_interevent_times
#ratio_min_max_T = ratio_min_max.T
#mean_ratio_min_max = np.mean(ratio_min_max_T, axis=0) 
#ratio_min_pair_max_bounds = \
#    np.array([abs(mean_ratio_min_pair_max - (np.percentile(ratio_min_pair_max_T, 2.5, axis=0))),
 #             abs(mean_ratio_min_pair_max - (np.percentile(ratio_min_pair_max_T, 97.5, axis=0)))])
#print('ratio_min_pair_max_bound', ratio_min_pair_max_bounds)
#ratio_min_max_bounds = \
#    np.array([abs(mean_ratio_min_max - (np.percentile(ratio_min_max_T, 2.5, axis=0))),
#              abs(mean_ratio_min_max - (np.percentile(ratio_min_max_T, 97.5, axis=0)))])
#sys.exit()
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
ax.set_ylim([1./1000000, 1./40])
ax.set_yscale('log')
ax.set_xlabel('COV')
ax.set_ylabel('Long-term rate (events per year)')
pyplot.savefig('cov_vs_lt_rate.png')


# Now just plot the means
pyplot.clf()
ax = pyplot.subplot(111)
mean_covs = []
#mean_ltrs = []
for i, cov_set in enumerate(covs):
    mean_cov = np.mean(cov_set)
    mean_covs.append(mean_cov)
#    mean_ltr = np.mean(long_term_rates[i])
#    mean_ltrs.append(mean_ltr)
#print(mean_covs, type(mean_covs))
#print(long_term_rates, type(long_term_rates))
pyplot.scatter(mean_covs, mean_ltr, marker = 's', c='0.1', s=25)
ax.set_xlim([0, 2.5])
ax.set_ylim([1./1000000, 1./40])
ax.set_yscale('log')
ax.set_xlabel('COV')
ax.set_ylabel('Long-term rate (events per year)')
pyplot.savefig('mean_cov_vs_lt_rate.png')

# Now plot basic statistics
pyplot.clf()
ax = pyplot.subplot(111)
pyplot.errorbar(max_interevent_times, min_interevent_times,
                yerr = min_interevent_times_bounds,
                ecolor = '0.4',
                linestyle="None")
pyplot.errorbar(max_interevent_times, min_interevent_times,
                   xerr = max_interevent_times_bounds,
                   ecolor = '0.6',
                   linestyle="None")
pyplot.scatter(max_interevent_times, min_interevent_times,
               marker = 's', c='0.1', s=25)
ax.set_xlabel('Maximum interevent time')
ax.set_ylabel('Minimum interevent time') 
ax.set_xscale('log')
ax.set_yscale('log')
# Label low-slip rate faults
for i, txt in enumerate(names):
    if max_interevent_times[i] > 10000:
        ax.annotate(txt, (max_interevent_times[i],
                          min_interevent_times[i]))  

# Linear fit only bottom end of data
indices = np.argwhere(max_interevent_times < 10000).flatten()
print('indices', indices)
print(max_interevent_times[indices])
print(min_paired_interevent_times[indices])
lf = np.polyfit(np.log10(max_interevent_times[indices]),
                   np.log10(min_interevent_times[indices]), 1)
xvals_short = np.arange(100, 1e4, 100)
log_yvals = lf[0]*np.log10(xvals_short) + lf[1]
yvals = np.power(10, log_yvals)
pyplot.plot(xvals_short, yvals)
# Add formula for linear fit to low-end of data
txt = 'Log(Y) = %.2fLog(x) + %.2f' % (lf[0], lf[1])
print(txt)
ax.annotate(txt, (800, 10000))

pyplot.savefig('min_vs_max_interevent_time.png')

# Plot minimum pairs
pyplot.clf()
ax = pyplot.subplot(111)
pyplot.errorbar(max_interevent_times, min_paired_interevent_times,
                yerr = min_paired_interevent_times_bounds,
                ecolor = '0.4',
                linestyle="None")
pyplot.errorbar(max_interevent_times, min_paired_interevent_times,
                xerr = max_interevent_times_bounds,
                ecolor = '0.6',
                linestyle="None")
pyplot.scatter(max_interevent_times, min_paired_interevent_times,
               marker = 's', c='0.1', s=25)
ax.set_xlabel('Maximum interevent time')
ax.set_ylabel('Minimum interevent time \n(mean of two shortest consecutive interevent times)')
ax.set_xscale('log')
ax.set_yscale('log') 
# Label low-slip rate faults
for i, txt in enumerate(names):
    if max_interevent_times[i] > 10000:
        ax.annotate(txt, (max_interevent_times[i],
                          min_paired_interevent_times[i]))
# Now fit with a regression in log-log space
xvals = np.arange(100, 2e6, 100) # For plotting
# Linear fit
lf = np.polyfit(np.log10(max_interevent_times),
                np.log10(min_paired_interevent_times), 1)
log_yvals = lf[0]*np.log10(xvals) + lf[1]
yvals = np.power(10, log_yvals)
pyplot.plot(xvals, yvals)

# Linear fit only bottom end of data
indices = np.argwhere(max_interevent_times < 10000).flatten()
print('indices', indices)
print(max_interevent_times[indices])
print(min_paired_interevent_times[indices])
lf = np.polyfit(np.log10(max_interevent_times[indices]),
                np.log10(min_paired_interevent_times[indices]), 1)
xvals_short = np.arange(100, 1e4, 100)
log_yvals = lf[0]*np.log10(xvals_short) + lf[1]
yvals = np.power(10, log_yvals)
pyplot.plot(xvals_short, yvals)
# Add formula for linear fit to low-end of data
txt = 'Log(Y) = %.2fLog(x) + %.2f' % (lf[0], lf[1])
print(txt)
ax.annotate(txt, (100, 10000))

# Quadratic fit
qf = np.polyfit(np.log10(max_interevent_times),
                      np.log10(min_paired_interevent_times), 2)
print(qf)
log_yvals = qf[0]*np.log10(xvals)**2 + qf[1]*np.log10(xvals) + qf[2]
yvals = np.power(10, log_yvals)
pyplot.plot(xvals, yvals)

# Power law fit
# Define function to fit
#def func_powerlaw(x, m, c, c0):
#    return c0 + x**m * c
#target_func = func_powerlaw
#popt, pcov = curve_fit(target_func, max_interevent_times, min_paired_interevent_times)
#print('popt', popt)
#pyplot.plot(xvals, target_func(xvals, *popt), '--')
pyplot.savefig('min_pair_vs_max_interevent_time.png')

# Similar plots, against long term rates
pyplot.clf()
ax = pyplot.subplot(111)  
print('mean_ltr', mean_ltr)
print('min_interevent_times', min_interevent_times)
print('ltr-bounds', ltr_bounds)
pyplot.errorbar(mean_ltr, min_interevent_times,
                yerr = min_interevent_times_bounds,
                ecolor = '0.4',
                linestyle="None")
pyplot.errorbar(mean_ltr, min_interevent_times,
                xerr = ltr_bounds,
                ecolor = '0.6',
                linestyle="None")
pyplot.scatter(mean_ltr, min_interevent_times,
               marker = 's', c='0.1', s=25)
ax.set_xlabel('Long-term rate')
ax.set_ylabel('Minimum interevent time')
ax.set_xscale('log')
ax.set_yscale('log') 
# Label low-slip rate faults
for i, txt in enumerate(names):
    if max_interevent_times[i] > 10000:
        ax.annotate(txt, (mean_ltr[i],
                          min_interevent_times[i]))
# Now fit with a regression in log-log space
#xvals = np.arange(10e-6, 10e-1, 1e-7) # For plotting
# Linear fit
#lf = np.polyfit(np.log10(mean_ltr),
#                np.log10(min_interevent_times), 1)
#log_yvals = lf[0]*np.log10(xvals) + lf[1]
#yvals = np.power(10, log_yvals)
#pyplot.plot(xvals, yvals)

# Linear fit only bottom end of data
indices = np.argwhere(mean_ltr > 5e-4).flatten()
print('indices', indices)
lf = np.polyfit(np.log10(mean_ltr[indices]),
                np.log10(min_interevent_times[indices]), 1)
xvals_short = np.arange(5e-4, 1e-2, 1e-4)
log_yvals = lf[0]*np.log10(xvals_short) + lf[1]
yvals = np.power(10, log_yvals)
pyplot.plot(xvals_short, yvals)
# Add formula for linear fit to low-end of data
txt = 'Log(Y) = %.2fLog(x) + %.2f' % (lf[0], lf[1])
ax.annotate(txt, (1e-4, 10000))      
pyplot.savefig('min_interevent_time_vs_ltr.png')

# Plot long term rate against minimum pair
pyplot.clf()
ax = pyplot.subplot(111)  
pyplot.errorbar(mean_ltr, min_paired_interevent_times,
                yerr = min_paired_interevent_times_bounds,
                ecolor = '0.4',
                linestyle="None")
pyplot.errorbar(mean_ltr, min_paired_interevent_times,
                xerr = ltr_bounds,
                ecolor = '0.6',
                linestyle="None")
pyplot.scatter(mean_ltr, min_paired_interevent_times,
               marker = 's', c='0.1', s=25)
ax.set_xlabel('Long-term rate')
ax.set_ylabel('Minimum interevent time \n(mean of two shortest consecutive interevent times)')
ax.set_xscale('log')
ax.set_yscale('log') 
# Label low-slip rate faults
for i, txt in enumerate(names):
    if max_interevent_times[i] > 10000:
        ax.annotate(txt, (mean_ltr[i],
                          min_interevent_times[i]))
# Now fit with a regression in log-log space
#xvals = np.arange(10e-6, 10e-1, 1e-7) # For plotting
# Linear fit
#lf = np.polyfit(np.log10(mean_ltr),
#                np.log10(min_paired_interevent_times), 1)
#log_yvals = lf[0]*np.log10(xvals) + lf[1]
#yvals = np.power(10, log_yvals)
#pyplot.plot(xvals, yvals)

# Linear fit only bottom end of data
indices = np.argwhere(mean_ltr > 5e-4).flatten()
print('indices', indices)
lf = np.polyfit(np.log10(mean_ltr[indices]),
                np.log10(min_paired_interevent_times[indices]), 1)
xvals_short = np.arange(5e-4, 1e-2, 1e-4)
log_yvals = lf[0]*np.log10(xvals_short) + lf[1]
yvals = np.power(10, log_yvals)
pyplot.plot(xvals_short, yvals)
# Add formula for linear fit to low-end of data
txt = 'Log(Y) = %.2fLog(x) + %.2f' % (lf[0], lf[1])
print(txt)
ax.annotate(txt, (1e-4, 10000))

pyplot.savefig('min_pair_vs_ltr.png')

# Now plot ratios against long term rates
pyplot.clf()
ax = pyplot.subplot(111)  
#ratio_min_pair_max = min_paired_interevent_times/max_interevent_times
print('ratio_min_pair_max_bounds', ratio_min_pair_max_bounds)
print(np.shape(ratio_min_pair_max_bounds))
pyplot.errorbar(mean_ltr, mean_ratio_min_pair_max,
                yerr = ratio_min_pair_max_bounds,
                ecolor = '0.4',
                linestyle="None")
pyplot.errorbar(mean_ltr, mean_ratio_min_pair_max,
                xerr = ltr_bounds,
                ecolor = '0.6',
                linestyle="None")
pyplot.scatter(mean_ltr, mean_ratio_min_pair_max,
               marker = 's', c='0.1', s=25)
ax.set_xlabel('Long-term rate')
ax.set_ylabel('Minimum pair interevent time: maximum interevent time')
ax.set_xscale('log')
ax.set_yscale('log') 
# Label low-slip rate faults
for i, txt in enumerate(names):
    if max_interevent_times[i] > 10000:
        ax.annotate(txt, (mean_ltr[i],
                          min_interevent_times[i]))

#log_yvals = lf[0]*np.log10(xvals) + lf[1]
#yvals = np.power(10, log_yvals)
#pyplot.plot(xvals, yvals)

# Linear fit only bottom end of data
indices = np.argwhere(mean_ltr > 5e-4).flatten()
print('indices', indices)
lf = np.polyfit(np.log10(mean_ltr[indices]),
                                np.log10(mean_ratio_min_pair_max[indices]), 1)
xvals_short = np.arange(5e-4, 1e-2, 1e-4)
log_yvals = lf[0]*np.log10(xvals_short) + lf[1]
yvals = np.power(10, log_yvals)
pyplot.plot(xvals_short, yvals)
# Add formula for linear fit to low-end of data
txt = 'Log(Y) = %.2fLog(x) + %.2f' % (lf[0], lf[1])
print(txt)
ax.annotate(txt, (1e-4, 10000))

pyplot.savefig('min_pair_max_ratio_vs_ltr.png')

