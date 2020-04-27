"""Script for sampling COV uncertainties on many faults and plotting them
"""

import os, sys
import ast
from glob import glob
from operator import itemgetter 
import numpy as np
from scipy.optimize import curve_fit
from scipy.odr import Model, RealData, ODR
from matplotlib import pyplot
from matplotlib.patches import PathPatch
from scipy.stats import kde
from QuakeRates.dataman.event_dates import EventSet
from QuakeRates.dataman.parse_oxcal import parse_oxcal
from QuakeRates.dataman.parse_age_sigma import parse_age_sigma
from QuakeRates.dataman.parse_params import parse_param_file, \
    get_event_sets, file_len

filepath = '../params'
param_file_list = glob(os.path.join(filepath, '*.txt'))
n_samples = 500  # Number of Monte Carlo samples of the eq chronologies
half_n = int(n_samples/2)
print(half_n)

# Define subset to take
faulting_styles = ['Reverse']
#faulting_styles = ['Normal']
#faulting_styles = ['Strike_slip'] 
#faulting_styles = ['all']
tectonic_regions = ['all']
#tectonic_regions = ['Plate_boundary_master', 'Plate_boundary_network']
min_number_events = 6

#Summarise for comment to add to figure filename
fig_comment = ''
for f in faulting_styles:
    fig_comment += f
    fig_comment += '_'
for t in tectonic_regions:
    fig_comment += t
    fig_comment += '_'
fig_comment += str(min_number_events)

covs = []
cov_bounds = []
burstinesses = []
burstiness_bounds = []
memory_coefficients = []
memory_bounds = []
long_term_rates = []
max_interevent_times = []
min_interevent_times = []
min_paired_interevent_times = []
std_min_paired_interevent_times = []
std_min_interevent_times = []
std_max_interevent_times = []
max_interevent_times_bounds = []
min_interevent_times_bounds = []
min_paired_interevent_times_bounds = []
ratio_min_pair_max = []
ratio_min_max = []
std_ratio_min_pair_max = []
std_ratio_min_max = []
ratio_min_pair_max_bounds =[]
ratio_min_max_bounds = []

names, event_sets, event_certainties, num_events = \
    get_event_sets(param_file_list, tectonic_regions,
                   faulting_styles, min_number_events)

for i, event_set in enumerate(event_sets):
    # Handle cases with uncertain number of events. Where events identification is
    # unsure, event_certainty is given a value of 0, compared with 1 for certain
    # events
    # First generate chronologies assuming all events are certain
    event_set.gen_chronologies(n_samples, observation_end=2019, min_separation=1)
    event_set.calculate_cov()
    event_set.cov_density()
    event_set.memory_coefficient()
    # Now calculate some statistics on the sampled chronologies
    event_set.basic_chronology_stats()
    min_paired_interevent_times.append(event_set.mean_minimum_pair_interevent_time)
    max_interevent_times.append(event_set.mean_maximum_interevent_time)
    min_interevent_times.append(event_set.mean_minimum_interevent_time)  
    std_min_paired_interevent_times.append(event_set.std_minimum_pair_interevent_time)
    std_min_interevent_times.append(event_set.std_minimum_interevent_time)
    std_max_interevent_times.append(event_set.std_maximum_interevent_time)
    if event_set.std_maximum_interevent_time == 0:
        print(name)
        sys.exit()
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
    ratio_min_pair_max.append(event_set.mean_ratio_min_pair_max)
    ratio_min_max.append(event_set.mean_ratio_min_max)
    std_ratio_min_pair_max.append(event_set.std_ratio_min_pair_max)
    std_ratio_min_max.append(event_set.std_ratio_min_max)
    ratio_min_pair_max_bounds.append([abs(event_set.mean_ratio_min_pair_max -
                                          event_set.ratio_min_pair_max_lb),
                                      abs(event_set.mean_ratio_min_pair_max -
                                          event_set.ratio_min_pair_max_ub)])
    ratio_min_max_bounds.append([abs(event_set.mean_ratio_min_max -
                                     event_set.ratio_min_max_lb),
                                 abs(event_set.mean_ratio_min_max -
                                     event_set.ratio_min_max_ub)])
    # Now generate chronologies assuming uncertain events did not occur
    if sum(event_certainties[i]) < event_set.num_events:
        indices = np.where(event_certainties[i] == 1)
        indices = list(indices[0])
#        print(indices[0], type(indices))
        events_subset = list(itemgetter(*indices)(event_set.event_list)) 
        event_set_certain = EventSet(events_subset)
        event_set_certain.gen_chronologies(n_samples, observation_end=2019, min_separation=1)
        event_set_certain.calculate_cov()
        event_set_certain.cov_density()
        event_set_certain.basic_chronology_stats()
        event_set_certain.memory_coefficient()
        combined_covs = np.concatenate([event_set.covs[:half_n],
                                        event_set_certain.covs[:half_n]])
        combined_burstiness = np.concatenate([event_set.burstiness[:half_n],
                                        event_set_certain.burstiness[:half_n]])
        combined_memory = np.concatenate([event_set.mem_coef[:half_n],
                                          event_set_certain.mem_coef[:half_n]])
        covs.append(combined_covs)
        burstinesses.append(combined_burstiness)
        memory_coefficients.append(combined_memory)
        cov_bounds.append([abs(np.mean(combined_covs) - \
                               min(event_set.cov_lb, event_set_certain.cov_lb)),
                           abs(np.mean(combined_covs) - \
                               max(event_set.cov_ub, event_set_certain.cov_ub))])
        burstiness_bounds.append([abs(np.mean(combined_burstiness) - \
                                      min(event_set.burstiness_lb,
                                          event_set_certain.burstiness_lb)),
                                  abs(np.mean(combined_burstiness) - \
                                      max(event_set.burstiness_ub,
                                          event_set_certain.burstiness_ub))])
        memory_bounds.append([abs(np.mean(combined_memory) - \
                                  min(event_set.memory_lb,
                                      event_set_certain.memory_lb)),
                              abs(np.mean(combined_memory) - \
                                  max(event_set.memory_ub,
                                      event_set_certain.memory_ub))])
        # Combine, taking n/2 samples from each set
        combined_ltrs = np.concatenate([event_set.long_term_rates[:half_n],
                                        event_set_certain.long_term_rates[:half_n]])
        print(len(combined_ltrs))
        long_term_rates.append(combined_ltrs)
    else:
        covs.append(event_set.covs)
        burstinesses.append(event_set.burstiness)
        memory_coefficients.append(event_set.mem_coef)
        long_term_rates.append(event_set.long_term_rates)
        cov_bounds.append([abs(event_set.mean_cov - event_set.cov_lb),
                          abs(event_set.mean_cov - event_set.cov_ub)])
        burstiness_bounds.append([abs(event_set.mean_burstiness - event_set.burstiness_lb),
                                  abs(event_set.mean_burstiness - event_set.burstiness_ub)])
        memory_bounds.append([abs(event_set.mean_mem_coef - event_set.memory_lb),
                              abs(event_set.mean_mem_coef - event_set.memory_ub)])
# Convert to numpy arrays and transpose where necessary
num_events = np.array(num_events)
max_interevent_times = np.array(max_interevent_times)
min_interevent_times = np.array(min_interevent_times)
min_paired_interevent_times = np.array(min_paired_interevent_times)
std_max_interevent_times = np.array(std_max_interevent_times)
std_min_interevent_times = np.array(std_min_interevent_times)
std_min_paired_interevent_times = np.array(std_min_paired_interevent_times)
print('std_max_interevent_times', std_max_interevent_times)
print('std_min_paired_interevent_times', std_min_paired_interevent_times)
max_interevent_times_bounds = np.array(max_interevent_times_bounds).T
min_interevent_times_bounds = np.array(min_interevent_times_bounds).T
min_paired_interevent_times_bounds = np.array(min_paired_interevent_times_bounds).T
long_term_rates_T = np.array(long_term_rates).T
mean_ltr = np.mean(long_term_rates_T, axis = 0)
std_ltr = np.std(long_term_rates_T, axis = 0)
ltr_bounds = np.array([abs(mean_ltr - (np.percentile(long_term_rates_T, 2.5, axis=0))),
                       abs(mean_ltr - (np.percentile(long_term_rates_T, 97.5, axis=0)))])
ratio_min_pair_max = np.array(ratio_min_pair_max)
ratio_min_max = np.array(ratio_min_max)
std_ratio_min_pair_max = np.array(std_ratio_min_pair_max)
std_ratio_min_max = np.array(std_ratio_min_max)
ratio_min_pair_max_bounds = np.array(ratio_min_pair_max_bounds).T
ratio_min_max_bounds = np.array(ratio_min_max_bounds).T
cov_bounds = np.array(cov_bounds).T
burstiness_bounds = np.array(burstiness_bounds).T
memory_bounds = np.array(memory_bounds).T 

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
figname = 'cov_vs_lt_rate_%s.png' % fig_comment
pyplot.savefig(figname)


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
colours = []
for mean_cov in mean_covs:
    if mean_cov <= 0.9:
        colours.append('b')
    elif mean_cov > 0.9 and mean_cov <= 1.1:
        colours.append('g')
    else:
        colours.append('r')
pyplot.errorbar(mean_covs, mean_ltr,
                yerr = ltr_bounds,
                ecolor = '0.4',
                linestyle="None")
pyplot.errorbar(mean_covs, mean_ltr,
                   xerr = cov_bounds,
                   ecolor = '0.6',
                   linestyle="None")
pyplot.scatter(mean_covs, mean_ltr, marker = 's', c=colours, s=25)
for i, txt in enumerate(names):
    if max_interevent_times[i] > 10:
        ax.annotate(txt[:4],
                    (mean_covs[i], mean_ltr[i]),
                    fontsize=8)
ax.set_xlim([0, 2.5])
ax.set_ylim([1./1000000, 1./40])
ax.set_yscale('log')
ax.set_xlabel('COV')
ax.set_ylabel('Long-term rate (events per year)')
figname = 'mean_cov_vs_lt_rate_%s.png' % fig_comment 
pyplot.savefig(figname)

# Plot burstiness against mean ltr
pyplot.clf()
ax = pyplot.subplot(111)
mean_bs = []
#mean_ltrs = []
for i, b_set in enumerate(burstinesses):
    mean_b = np.mean(b_set)
    mean_bs.append(mean_b)
colours = []
for mean_b in mean_bs:
    if mean_b <= -0.05:
        colours.append('b')
    elif mean_b > -0.05 and mean_b <= 0.05:
        colours.append('g')
    else:
        colours.append('r')

pyplot.errorbar(mean_bs, mean_ltr,
                yerr = ltr_bounds,
                ecolor = '0.4',
                linestyle="None")
pyplot.errorbar(mean_bs, mean_ltr,
                   xerr = burstiness_bounds,
                   ecolor = '0.6',
                   linestyle="None")
pyplot.scatter(mean_bs, mean_ltr, marker = 's', c=colours, s=25)
for i, txt in enumerate(names):
    if max_interevent_times[i] > 10:
        ax.annotate(txt[:4],
                    (mean_bs[i], mean_ltr[i]),
                    fontsize=8)
ax.set_xlim([-1, 1])
ax.set_ylim([1./1000000, 1./40])
ax.set_yscale('log')
ax.set_xlabel('B')
ax.set_ylabel('Long-term rate (events per year)')
figname = 'burstiness_vs_lt_rate_%s.png' % fig_comment 
pyplot.savefig(figname)

# Plot memory coefficients against long term rates
pyplot.clf()
ax = pyplot.subplot(111)
mean_mems = []
#mean_ltrs = []
for i, mem_set in enumerate(memory_coefficients):
    mean_mem = np.mean(mem_set)
    mean_mems.append(mean_mem)
colours = []
for mean_mem in mean_mems:
    if mean_mem <= -0.05:
        colours.append('b')
    elif mean_mem > -0.05 and mean_mem <= 0.05:
        colours.append('g')
    else:
        colours.append('r')

pyplot.errorbar(mean_mems, mean_ltr,
                yerr = ltr_bounds,
                ecolor = '0.4',
                linestyle="None")
pyplot.errorbar(mean_mems, mean_ltr,
                   xerr = memory_bounds,
                   ecolor = '0.6',
                   linestyle="None")
pyplot.scatter(mean_mems, mean_ltr, marker = 's', c=colours, s=25)
for i, txt in enumerate(names):
    if max_interevent_times[i] > 10:
        ax.annotate(txt[:4],
                    (mean_mems[i], mean_ltr[i]),
                    fontsize=8)
#ax.set_xlim([-1, 1])
ax.set_ylim([1./1000000, 1./40])
ax.set_yscale('log')
ax.set_xlabel('M')
ax.set_ylabel('Long-term rate (events per year)')
figname = 'memory_coefficient_vs_lt_rate_%s.png' % fig_comment
pyplot.savefig(figname)

# Plot burstiness against memory coefficient
# Plot burstiness against mean ltr
pyplot.clf()
ax = pyplot.subplot(111)
colours = []
for mean_b in mean_bs:
    if mean_b <= -0.05:
        colours.append('b')
    elif mean_b > -0.05 and mean_b <= 0.05:
        colours.append('g')
    else:
        colours.append('r')
pyplot.errorbar(mean_bs, mean_mems,
                yerr = memory_bounds,
                ecolor = '0.4',
                linestyle="None")
pyplot.errorbar(mean_bs, mean_mems,
                   xerr = burstiness_bounds,
                   ecolor = '0.6',
                   linestyle="None")
pyplot.scatter(mean_bs, mean_mems, marker = 's', c=colours, s=25)
for i, txt in enumerate(names):
    if max_interevent_times[i] > 10:
        ax.annotate(txt[:4],
                    (mean_bs[i], mean_mems[i]),
                    fontsize=8)
#ax.set_xlim([-1, 1])
#ax.set_ylim([1./1000000, 1./40])
#ax.set_yscale('log')
ax.set_xlabel('B')
ax.set_ylabel('M')
figname = 'burstiness_vs_memory_coefficient_%s.png' % fig_comment 
pyplot.savefig(figname)


# Plot COV against number of events to look at sampling biases
pyplot.clf()
ax = pyplot.subplot(111)
mean_covs = []
#mean_ltrs = []
for i, cov_set in enumerate(covs):
    mean_cov = np.mean(cov_set)
    mean_covs.append(mean_cov)
colours = []
for mean_cov in mean_covs:
    if mean_cov <= 0.9:
        colours.append('b')
    elif mean_cov > 0.9 and mean_cov <= 1.1:
        colours.append('g')
    else:
        colours.append('r')

pyplot.errorbar(mean_covs, num_events,
                   xerr = cov_bounds,
                   ecolor = '0.6',
                   linestyle="None")
pyplot.scatter(mean_covs, num_events, marker = 's', c=colours, s=25)
for i, txt in enumerate(names):
    if max_interevent_times[i] > 10:
        ax.annotate(txt[:4],
                    (mean_covs[i], num_events[i]),
                    fontsize=8)
#ax.set_xlim([0, 2.5])
#ax.set_ylim([1./1000000, 1./40])
#ax.set_yscale('log')
ax.set_xlabel('COV')
ax.set_ylabel('Number of events in earthquake record')
figname = 'mean_cov_vs_number_events_%s.png' % fig_comment
pyplot.savefig(figname)

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
               marker = 's', c=colours, s=25)
ax.set_xlabel('Maximum interevent time')
ax.set_ylabel('Minimum interevent time') 
ax.set_xscale('log')
ax.set_yscale('log')
# Label low-slip rate faults
for i, txt in enumerate(names):
    if max_interevent_times[i] > 10:
        ax.annotate(txt[:4],
                    (max_interevent_times[i], min_interevent_times[i]),
                    fontsize=8)

# Linear fit only bottom end of data
indices = np.argwhere(max_interevent_times < 10000).flatten()
indices_slow_faults = np.argwhere(max_interevent_times >= 10000).flatten()
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
figname = 'min_vs_max_interevent_time_%s.png' % fig_comment
pyplot.savefig(figname)

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
               marker = 's', c=colours, s=25)
ax.set_xlabel('Maximum interevent time')
ax.set_ylabel('Minimum interevent time \n(mean of two shortest consecutive interevent times)')
ax.set_xscale('log')
ax.set_yscale('log') 
# Label low-slip rate faults
for i, txt in enumerate(names):
    if max_interevent_times[i] > 10:
        ax.annotate(txt[:4],
                    (max_interevent_times[i], min_paired_interevent_times[i]),
                    fontsize=8)

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
# Linear fit with data uncertainties
def linear_func(B, x):
    return B[0]*x + B[1]
linear_model = Model(linear_func)
data = RealData(np.log10(max_interevent_times[indices]),
                np.log10(min_paired_interevent_times[indices]),
                sx = np.log10(std_max_interevent_times[indices]),
                sy = np.log10(std_min_paired_interevent_times[indices]))
# Set up ODR with the model and data
odr = ODR(data, linear_model, beta0=[1., 1.])
# Run the regression.
out = odr.run()
print('Regression output')
out.pprint()

log_y_fit = linear_func(out.beta, np.log10(xvals_short))
y_fit = np.power(10, log_y_fit)
#popt, pcov = curve_fit(target_func, max_interevent_times, min_paired_interevent_times)
#print('popt', popt)
pyplot.plot(xvals_short, y_fit, '--')
txt = 'Log(Y) =  %.2fLog(x) + %.2f' % (out.beta[0], out.beta[1])
print(txt)
ax.annotate(txt, (100, 40000))
figname = 'min_pair_vs_max_interevent_time_%s.png' % fig_comment
pyplot.savefig(figname)

# Similar plots, against long term rates
pyplot.clf()
ax = pyplot.subplot(111)  
pyplot.errorbar(mean_ltr, min_interevent_times,
                yerr = min_interevent_times_bounds,
                ecolor = '0.4',
                linestyle="None")
pyplot.errorbar(mean_ltr, min_interevent_times,
                xerr = ltr_bounds,
                ecolor = '0.6',
                linestyle="None")
pyplot.scatter(mean_ltr, min_interevent_times,
               marker='s', c=colours, s=25)
ax.set_xlabel('Long-term rate')
ax.set_ylabel('Minimum interevent time')
ax.set_xscale('log')
ax.set_yscale('log') 
# Label low-slip rate faults
for i, txt in enumerate(names):
    if max_interevent_times[i] > 10:
        ax.annotate(txt[:4],
                    (mean_ltr[i], min_interevent_times[i]),
                    fontsize=8)

# Now fit with a regression in log-log space
#xvals = np.arange(10e-6, 10e-1, 1e-7) # For plotting
# Linear fit
#lf = np.polyfit(np.log10(mean_ltr),
#                np.log10(min_interevent_times), 1)
#log_yvals = lf[0]*np.log10(xvals) + lf[1]
#yvals = np.power(10, log_yvals)
#pyplot.plot(xvals, yvals)

# Linear fit only bottom end of data
indices = np.argwhere(mean_ltr > 2e-4).flatten()
lf = np.polyfit(np.log10(mean_ltr[indices]),
                np.log10(min_interevent_times[indices]), 1)
xvals_short = np.arange(5e-4, 1e-2, 1e-4)
log_yvals = lf[0]*np.log10(xvals_short) + lf[1]
yvals = np.power(10, log_yvals)
pyplot.plot(xvals_short, yvals)
# Add formula for linear fit to low-end of data
txt = 'Log(Y) = %.2fLog(x) + %.2f' % (lf[0], lf[1])
ax.annotate(txt, (1e-4, 10000))
figname = 'min_interevent_time_vs_ltr_%s.png' % fig_comment
pyplot.savefig(figname)

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
               marker='s', c=colours, s=25)
#c='0.1', s=25)
ax.set_xlabel('Long-term rate')
ax.set_ylabel('Minimum interevent time \n(mean of two shortest consecutive interevent times)')
ax.set_xscale('log')
ax.set_yscale('log') 
# Label low-slip rate faults
for i, txt in enumerate(names):
    if max_interevent_times[i] > 10:
        ax.annotate(txt[:4],
                    (mean_ltr[i], min_paired_interevent_times[i]),
                    fontsize=8)
# Now fit with a regression in log-log space
#xvals = np.arange(10e-6, 10e-1, 1e-7) # For plotting
# Linear fit
#lf = np.polyfit(np.log10(mean_ltr),
#                np.log10(min_paired_interevent_times), 1)
#log_yvals = lf[0]*np.log10(xvals) + lf[1]
#yvals = np.power(10, log_yvals)
#pyplot.plot(xvals, yvals)

# Linear fit only bottom end of data
indices = np.argwhere(mean_ltr > 2e-4).flatten()
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
figname = 'min_pair_vs_ltr_%s.png' % fig_comment 
pyplot.savefig(figname)

# Plot long term rate against maximum interevent time
pyplot.clf()
ax = pyplot.subplot(111)  
pyplot.errorbar(mean_ltr, max_interevent_times,
                yerr = max_interevent_times_bounds,
                ecolor = '0.4',
                linestyle="None")
pyplot.errorbar(mean_ltr, max_interevent_times,
                xerr = ltr_bounds,
                ecolor = '0.6',
                linestyle="None")
pyplot.scatter(mean_ltr, max_interevent_times,
               marker='s', c=colours, s=25)
#c='0.1', s=25)
ax.set_xlabel('Long-term rate')
ax.set_ylabel('Maximum interevent time')
ax.set_xscale('log')
ax.set_yscale('log') 
# Label low-slip rate faults
for i, txt in enumerate(names):
    if max_interevent_times[i] > 10:
        ax.annotate(txt[:4],
                    (mean_ltr[i], max_interevent_times[i]),
                    fontsize=8)
# Now fit with a regression in log-log space
#xvals = np.arange(10e-6, 10e-1, 1e-7) # For plotting
# Linear fit
#lf = np.polyfit(np.log10(mean_ltr),
#                np.log10(min_paired_interevent_times), 1)
#log_yvals = lf[0]*np.log10(xvals) + lf[1]
#yvals = np.power(10, log_yvals)
#pyplot.plot(xvals, yvals)

# Linear fit only bottom end of data
indices = np.argwhere(mean_ltr > 2e-10).flatten() # All data for now
lf = np.polyfit(np.log10(mean_ltr[indices]),
                np.log10(max_interevent_times[indices]), 1)
xvals_short = np.arange(2e-6, 1e-2, 1e-6)
log_yvals = lf[0]*np.log10(xvals_short) + lf[1]
yvals = np.power(10, log_yvals)
pyplot.plot(xvals_short, yvals)
# Add formula for linear fit to low-end of data
txt = 'Log(Y) = %.2fLog(x) + %.2f' % (lf[0], lf[1])
print(txt)
ax.annotate(txt, (1e-4, 10000))
figname = 'max_interevent_time_vs_ltr_%s.png' % fig_comment
pyplot.savefig(figname)


# Now plot ratios against long term rates
pyplot.clf()
ax = pyplot.subplot(111)  
pyplot.errorbar(mean_ltr, ratio_min_pair_max,
                yerr = ratio_min_pair_max_bounds,
                ecolor = '0.4',
                linestyle="None")
pyplot.errorbar(mean_ltr, ratio_min_pair_max,
                xerr = ltr_bounds,
                ecolor = '0.6',
                linestyle="None")
pyplot.scatter(mean_ltr, ratio_min_pair_max,
               marker='s', c=colours, s=25)
ax.set_xlabel('Long-term rate')
ax.set_ylabel('Minimum pair interevent time: maximum interevent time')
ax.set_xscale('log')
ax.set_yscale('log') 
# Label low-slip rate faults
for i, txt in enumerate(names):
    if max_interevent_times[i] > 10:
        ax.annotate(txt[:4],
                    (mean_ltr[i], ratio_min_pair_max[i]),
                    fontsize=8)

#log_yvals = lf[0]*np.log10(xvals) + lf[1]
#yvals = np.power(10, log_yvals)
#pyplot.plot(xvals, yvals)

# Linear fit high and low long term rate data separately
indices = np.argwhere(mean_ltr > 2e-4).flatten()
indices_slow_faults = np.argwhere(mean_ltr <= 2e-4).flatten()
lf = np.polyfit(np.log10(mean_ltr[indices]),
                                np.log10(ratio_min_pair_max[indices]), 1)
xvals_short = np.arange(2e-4, 5e-2, 1e-4)
log_yvals = lf[0]*np.log10(xvals_short) + lf[1]
yvals = np.power(10, log_yvals)
pyplot.plot(xvals_short, yvals)
# Add formula for linear fit to low-end of data
txt = 'Log(Y) = %.2fLog(x) + %.2f' % (lf[0], lf[1])
print(txt)
ax.annotate(txt, (5e-4, 1e-2))

# Slow long-term rates
if len(indices_slow_faults) > 0:
    lf = np.polyfit(np.log10(mean_ltr[indices_slow_faults]),
                    np.log10(ratio_min_pair_max[indices_slow_faults]), 1)
    xvals_short = np.arange(2e-6, 2e-4, 1e-6)
    log_yvals = lf[0]*np.log10(xvals_short) + lf[1]
    yvals = np.power(10, log_yvals)
    pyplot.plot(xvals_short, yvals)
    # Add formula for linear fit to low-end of data
    txt = 'Log(Y) = %.2fLog(x) + %.2f' % (lf[0], lf[1])
    print(txt)
    ax.annotate(txt, (1e-5, 5e-3))

"""
#Orthogonal linear fit high long-term rate data
linear_model = Model(linear_func)
data = RealData(np.log10(mean_ltr[indices]),
                np.log10(ratio_min_pair_max[indices]),
                sx = np.log10(std_ltr[indices]),
                sy = np.log10(std_ratio_min_pair_max[indices]))
# Set up ODR with the model and data
odr = ODR(data, linear_model, beta0=[1., 1.])
# Run the regression.
out = odr.run()
print('Regression output')
out.pprint()
log_y_fit = linear_func(out.beta, np.log10(xvals_short))
y_fit = np.power(10, log_y_fit)
#popt, pcov = curve_fit(target_func, max_interevent_times, min_paired_interevent_times)
#print('popt', popt)
pyplot.plot(xvals_short, y_fit, '--')
txt = 'Log(Y) =  %.2fLog(x) + %.2f' % (out.beta[0], out.beta[1])
print(txt)
ax.annotate(txt, (1e-3, 1e-2))

#Orthogonal linear fit low long-term rate data
linear_model = Model(linear_func)
data = RealData(np.log10(mean_ltr[indices_slow_faults]),
                np.log10(ratio_min_pair_max[indices_slow_faults]),
                sx = np.log10(std_ltr[indices_slow_faults]),
                sy = np.log10(std_min_paired_interevent_times[indices_slow_faults]))
# Set up ODR with the model and data
odr = ODR(data, linear_model, beta0=[1., 1.])
# Run the regression.
out = odr.run()
print('Regression output')
out.pprint()
xvals_short = np.arange(2e-6, 2e-4, 1e-6)
log_y_fit = linear_func(out.beta, np.log10(xvals_short))
y_fit = np.power(10, log_y_fit)
#popt, pcov = curve_fit(target_func, max_interevent_times, min_paired_interevent_times)
#print('popt', popt)
pyplot.plot(xvals_short, y_fit, '--')
txt = 'Log(Y) =  %.2fLog(x) + %.2f' % (out.beta[0], out.beta[1])
print(txt)
ax.annotate(txt, (1e-5, 5e-3))
"""

figname = 'min_pair_max_ratio_vs_ltr_%s.png' % fig_comment 
pyplot.savefig(figname)

# Now plot ratios against long term rates
pyplot.clf()
ax = pyplot.subplot(111)  
pyplot.errorbar(mean_ltr, ratio_min_max,
                yerr = ratio_min_max_bounds,
                ecolor = '0.4',
                linestyle="None")
pyplot.errorbar(mean_ltr, ratio_min_max,
                xerr = ltr_bounds,
                ecolor = '0.6',
                linestyle="None")
pyplot.scatter(mean_ltr, ratio_min_max,
               marker = 's', c=colours, s=25)
ax.set_xlabel('Long-term rate')
ax.set_ylabel('Minimum interevent time: maximum interevent time')
ax.set_xscale('log')
ax.set_yscale('log') 
# Label low-slip rate faults
for i, txt in enumerate(names):
    if max_interevent_times[i] > 10:
        ax.annotate(txt[:4],
                    (mean_ltr[i], ratio_min_max[i]),
                    fontsize=8)

#log_yvals = lf[0]*np.log10(xvals) + lf[1]
#yvals = np.power(10, log_yvals)
#pyplot.plot(xvals, yvals)

# Linear fit only bottom end of data
indices = np.argwhere(mean_ltr > 5e-4).flatten()
lf = np.polyfit(np.log10(mean_ltr[indices]),
                                np.log10(ratio_min_max[indices]), 1)
xvals_short = np.arange(5e-4, 1e-2, 1e-4)
log_yvals = lf[0]*np.log10(xvals_short) + lf[1]
yvals = np.power(10, log_yvals)
pyplot.plot(xvals_short, yvals)
# Add formula for linear fit to low-end of data
txt = 'Log(Y) = %.2fLog(x) + %.2f' % (lf[0], lf[1])
print(txt)
ax.annotate(txt, (1e-4, 1e-3))
figname = 'min_max_ratio_vs_ltr_%s.png' % fig_comment
pyplot.savefig(figname)
