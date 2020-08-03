"""Script for sampling COV uncertaintiesA on many faults and plotting them
"""

import os, sys
import ast
from glob import glob
from operator import itemgetter 
from re import finditer
import numpy as np
from scipy.optimize import curve_fit
from scipy.odr import Model, RealData, ODR
import scipy.odr.odrpack as odrpack
from scipy.stats import expon, ks_2samp, kstest
from matplotlib import pyplot
from matplotlib.patches import PathPatch
import matplotlib.gridspec as gridspec
from scipy.stats import binom, kde
from adjustText import adjust_text
from QuakeRates.dataman.event_dates import EventSet
from QuakeRates.dataman.parse_oxcal import parse_oxcal
from QuakeRates.dataman.parse_age_sigma import parse_age_sigma
from QuakeRates.dataman.parse_params import parse_param_file, \
    get_event_sets, file_len
from QuakeRates.utilities.bilinear import bilinear_reg_zero_slope, \
    bilinear_reg_fix, bilinear_reg_fix_zero_slope
from QuakeRates.utilities.memory_coefficient import burstiness

filepath = '../params'
param_file_list = glob(os.path.join(filepath, '*.txt'))

param_file_list_NZ = ['Akatore4eventBdy_output.txt',
                      'AlpineHokuriCk_Berryman_2012_simple.txt',
                      'AlpineSouthWestland_Cochran_2017_simple.txt',
                      'AwatereEast_Nicol_2016_simple.txt',
                      'ClarenceEast_Nicol_2016_simple.txt',
                      'CloudyFault_Nicol_2016_simple.txt',
                      'Dunstan6_GNS_unpub_simple.txt',
                      'HopeConway_Hatem_2019_simple.txt',
                      'Hope_Khajavi_2016_simple.txt',
                      'Ihaia_Nicol_2016_simple.txt',
                      'Oaonui_Nicol_2016_simple.txt',
                      'Ohariu_Nicol_2016_simple.txt',
                      'Paeroa_Nicol_2016_simple.txt',
                      'Pihama_Nicol_2016_simple.txt',
                      'PortersPassEast_Nicol_2016_simple.txt',
                      'Ngakuru_Nicol_2016_simple.txt',
                      'Mangatete_Nicol_2016_simple.txt',
                      'Rangipo_Nicol_2016_simple.txt',
                      'Rotoitipakau_Nicol_2016_simple.txt',
                      'Rotohauhau_Nicol_2016_simple.txt',
                      'Snowden_Nicol_2016_simple.txt',
                      'Vernon_Nicol_2016_simple.txt',
                      'WairarapaSouth_Nicol_2016_simple.txt',
                      'Wairau_Nicol_2018_simple.txt',
                      'Waimana_Nicol_2016_simple.txt',
                      'Wellington_Langridge_2011_simple.txt',
                      'Waitangi_simple.txt',
                      'Whakatane_Nicol_2016_simple.txt',
                      'Whirinaki_Nicol_2016_simple.txt']
#param_file_list = []
#for f in param_file_list_NZ:
#    param_file_list.append(os.path.join(filepath, f))
n_samples = 1000  # Number of Monte Carlo samples of the eq chronologies
half_n = int(n_samples/2)
print(half_n)
annotate_plots = False # If True, lable each fault on the plot
plot_folder = './plots'
if not os.path.exists(plot_folder):
    os.makedirs(plot_folder)

# Define subset to take
#faulting_styles = ['Reverse']
#faulting_styles = ['Normal']
#faulting_styles = ['Strike_slip'] 
faulting_styles = ['all']
tectonic_regions = ['all']
#tectonic_regions = ['Intraplate_noncratonic', 'Intraplate_cratonic', 'Near_plate_boundary']
#tectonic_regions = ['Plate_boundary_master', 'Plate_boundary_network']
#tectonic_regions = ['Plate_boundary_network', 'Near_plate_boundary'] 
#tectonic_regions = ['Plate_boundary_master']
#tectonic_regions = ['Subduction']
#tectonic_regions = ['Near_plate_boundary']
min_number_events = 6

#Summarise for comment to add to figure filename
fig_comment = ''
#fig_comment = 'NZ_examples_'
for f in faulting_styles:
    fig_comment += f
    fig_comment += '_'
for t in tectonic_regions:
    fig_comment += t
    fig_comment += '_'
fig_comment += str(min_number_events)

def piecewise_linear(x, x0, y0, k1, k2):
    return np.piecewise(x, [x < x0], [lambda x:k1*x + y0-k1*x0, lambda x:k2*x + y0-k2*x0])
def camel_case_split(identifier):
    matches = finditer('.+?(?:(?<=[a-z])(?=[A-Z])|(?<=[A-Z])(?=[A-Z][a-z])|$)', identifier)
    return [m.group(0) for m in matches]

plot_colours = []
covs = []
cov_bounds = []
burstinesses = []
burstiness_bounds = []
burstiness_stds = []
burstinesses_expon = []
memory_coefficients = []
memory_bounds = []
memory_stds = []
memory_spearman_coefficients = []
memory_spearman_bounds = []
memory_spearman_lag2_coef = []
memory_spearman_lag2_bounds = []
long_term_rates = []
long_term_rate_stds = []
slip_rates = []
slip_rate_stds = []
slip_rate_bounds = []
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
references = []
# Get citations for each dataset from filename
for s in param_file_list:
    sp = s.split('_')
    if sp[0].split('/')[2] in names:
        references.append(sp[1] + ' ' + sp[2])
n_faults = len(names)

for i, event_set in enumerate(event_sets):
    # Handle cases with uncertain number of events. Where events identification is
    # unsure, event_certainty is given a value of 0, compared with 1 for certain
    # events
    # First generate chronologies assuming all events are certain
#    event_set.name = names[i]
    event_set.gen_chronologies(n_samples, observation_end=2019, min_separation=1)
    event_set.calculate_cov()
    event_set.cov_density()
    event_set.memory_coefficient()
    event_set.memory_spearman_rank_correlation()
    # Now calculate some statistics on the sampled chronologies
    event_set.basic_chronology_stats()
    # Plot histogram of interevent times
    figfile = os.path.join(plot_folder, ('interevent_times_%s.png' % names[i]))
    event_set.plot_interevent_time_hist(fig_filename=figfile)
    min_paired_interevent_times.append(event_set.mean_minimum_pair_interevent_time)
    max_interevent_times.append(event_set.mean_maximum_interevent_time)
    min_interevent_times.append(event_set.mean_minimum_interevent_time)  
    std_min_paired_interevent_times.append(event_set.std_minimum_pair_interevent_time)
    std_min_interevent_times.append(event_set.std_minimum_interevent_time)
    std_max_interevent_times.append(event_set.std_maximum_interevent_time)
    if event_set.std_maximum_interevent_time == 0:
        print('Zero std_maximum_interevent_time for ', names[i])

    slip_rates.append(event_set.slip_rates[0])
    slip_rate_bounds.append([event_set.slip_rates[1], event_set.slip_rates[2]])
    slip_rate_stds.append(abs(np.log10(event_set.slip_rates[2]) - \
                              np.log10(event_set.slip_rates[1]))/4) # Approx from 95% intervals
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
    # Generate random exponentially distributed samples of length num_events - 1
    # i.e. the number of inter-event times in the chronology. These will be used
    # later for testing
    scale = 100 # Fix scale, as burstiness is independent of scale for exponentiall distribution
    ie_times_expon = expon(scale=scale).rvs(size=(n_samples*(event_set.num_events-1)))
    ie_times_expon = np.reshape(np.array(ie_times_expon), (n_samples, (event_set.num_events-1)))
    ie_times_expon_T = ie_times_expon.T
    burst_expon = burstiness(ie_times_expon_T)
    # Now generate chronologies assuming uncertain events did not occur
    if sum(event_certainties[i]) < event_set.num_events:
        indices = np.where(event_certainties[i] == 1)
        indices = list(indices[0])
#        print(indices[0], type(indices))
        events_subset = list(itemgetter(*indices)(event_set.event_list)) 
        event_set_certain = EventSet(events_subset)
        event_set_certain.name = names[i] 
        event_set_certain.gen_chronologies(n_samples, observation_end=2019, min_separation=1)
        event_set_certain.calculate_cov()
        event_set_certain.cov_density()
        event_set_certain.basic_chronology_stats()
        event_set_certain.memory_coefficient()
        event_set_certain.memory_spearman_rank_correlation()
        # Generate random exponentially distributed samples of length num_events - 1
        # i.e. the number of inter-event times in the chronology. These will be used
        # later for testing
        ie_times_expon_certain = expon(scale=scale).rvs(size=(n_samples*(len(indices)-1)))
        ie_times_expon_certain = np.reshape(np.array(ie_times_expon_certain), (n_samples, (len(indices)-1)))
        ie_times_expon_certain_T = ie_times_expon_certain.T
        burst_expon_certain = burstiness(ie_times_expon_certain_T)
        # Now combine results from certain chronolgies with uncertain ones
        combined_covs = np.concatenate([event_set.covs[:half_n],
                                        event_set_certain.covs[:half_n]])
        combined_burstiness = np.concatenate([event_set.burstiness[:half_n],
                                        event_set_certain.burstiness[:half_n]])
        combined_memory = np.concatenate([event_set.mem_coef[:half_n],
                                          event_set_certain.mem_coef[:half_n]])
        combined_memory_spearman = np.concatenate([event_set.rhos[:half_n],
                                                   event_set_certain.rhos[:half_n]])
        combined_memory_spearman_lag2 = np.concatenate([event_set.rhos2[:half_n],
                                                        event_set_certain.rhos2[:half_n]])
        combined_burst_expon = np.concatenate([burst_expon[:half_n],
                                               burst_expon_certain[:half_n]])
        covs.append(combined_covs)
        burstinesses.append(combined_burstiness)
        memory_coefficients.append(combined_memory)
        memory_stds.append(np.std(np.array(combined_memory)))
        memory_spearman_coefficients.append(combined_memory_spearman)
        memory_spearman_lag2_coef.append(combined_memory_spearman_lag2)
        burstinesses_expon.append(combined_burst_expon)
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
        memory_spearman_bounds.append([abs(np.mean(combined_memory_spearman) - \
                                           min(event_set.rho_lb,
                                               event_set_certain.rho_lb)),
                                       abs(np.mean(combined_memory_spearman) - \
                                           max(event_set.rho_ub,
                                               event_set_certain.rho_ub))])
        memory_spearman_lag2_bounds.append([abs(np.mean(combined_memory_spearman_lag2) - \
                                                min(event_set.rho2_lb,
                                                    event_set_certain.rho2_lb)),
                                            abs(np.mean(combined_memory_spearman_lag2) - \
                                                max(event_set.rho2_ub,
                                                    event_set_certain.rho2_ub))])
        # Combine, taking n/2 samples from each set
        combined_ltrs = np.concatenate([event_set.long_term_rates[:half_n],
                                        event_set_certain.long_term_rates[:half_n]])
        burstiness_stds.append(np.std(combined_burstiness))
        print(len(combined_ltrs))
        long_term_rates.append(combined_ltrs)
        long_term_rate_stds.append(np.std(combined_ltrs))
    else:
        covs.append(event_set.covs)
        burstinesses.append(event_set.burstiness)
        memory_coefficients.append(event_set.mem_coef)
        memory_stds.append(np.std(np.array(event_set.mem_coef)))
        memory_spearman_coefficients.append(event_set.rhos)
        memory_spearman_lag2_coef.append(event_set.rhos2)
        long_term_rates.append(event_set.long_term_rates)
        burstinesses_expon.append(burst_expon)
        cov_bounds.append([abs(event_set.mean_cov - event_set.cov_lb),
                          abs(event_set.mean_cov - event_set.cov_ub)])
        burstiness_bounds.append([abs(event_set.mean_burstiness - event_set.burstiness_lb),
                                  abs(event_set.mean_burstiness - event_set.burstiness_ub)])
        memory_bounds.append([abs(event_set.mean_mem_coef - event_set.memory_lb),
                              abs(event_set.mean_mem_coef - event_set.memory_ub)])
        memory_spearman_bounds.append([abs(event_set.mean_rho - event_set.rho_lb),
                                        abs(event_set.mean_rho - event_set.rho_ub)])
        memory_spearman_lag2_bounds.append([abs(event_set.mean_rho2 - event_set.rho2_lb),
                                            abs(event_set.mean_rho2 - event_set.rho2_ub)])
        burstiness_stds.append(event_set.std_burstiness)
        long_term_rate_stds.append(np.mean(long_term_rates)) 
    # Get colours for plotting later
    if event_set.faulting_style == 'Normal':
        plot_colours.append('r')
    elif event_set.faulting_style == 'Reverse':
        plot_colours.append('b')
    elif event_set.faulting_style == 'Strike_slip':
        plot_colours.append('g')
    else:
        plot_colours.append('k')

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
long_term_rate_stds = np.array(long_term_rate_stds)
slip_rates = np.array(slip_rates).T
slip_rate_bounds = np.array(slip_rate_bounds).T
slip_rate_stds = np.array(slip_rate_stds).T
print('Mean_ltr', mean_ltr)
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
burstiness_stds = np.array(burstiness_stds)
burstiness_expon = np.array(burstinesses_expon)
memory_stds = np.array(memory_stds)
memory_bounds = np.array(memory_bounds).T 
memory_spearman_bounds = np.array(memory_spearman_bounds).T
memory_spearman_lag2_bounds = np.array(memory_spearman_lag2_bounds).T
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
    try:
        k = kde.gaussian_kde(cov_samples)
    except:
        msg = 'Skipping %s' % names[i]
        print(msg)
        continue
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


# Now just plot the means and 95% error bars
pyplot.clf()
ax = pyplot.subplot(111)
mean_covs = []
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
pyplot.errorbar(mean_ltr, mean_covs,
                xerr = ltr_bounds,
                ecolor = '0.3',
                elinewidth=0.7,
                linestyle="None",
                zorder=1)
pyplot.errorbar(mean_ltr, mean_covs,
                yerr = cov_bounds,
                ecolor = '0.3',
                elinewidth=0.7,
                linestyle="None",
                zorder=1)
pyplot.scatter(mean_ltr, mean_covs, marker = 's', c=plot_colours,
               s=25, zorder=2)
for i, txt in enumerate(names):
    if max_interevent_times[i] > 10 and annotate_plots:
        ax.annotate(txt[:4],
                    (mean_ltr[i], mean_covs[i]),
                    fontsize=8)
ax.set_ylim([0, 2.5])
ax.set_xlim([1./1000000, 1./40])
ax.set_xscale('log')
ax.set_xlabel('Long-term rate (events per year)')
ax.set_ylabel('COV')
figname = 'mean_cov_vs_lt_rate_%s.png' % fig_comment 
pyplot.savefig(figname)

# Plot burstiness against mean ltr
pyplot.clf()
ax = pyplot.subplot(111)
mean_bs = []
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

pyplot.errorbar(mean_ltr, mean_bs,
                xerr = ltr_bounds,
                ecolor = '0.3',
                elinewidth=0.7,
                linestyle="None",
                zorder=1)
pyplot.errorbar(mean_ltr, mean_bs,
                yerr = burstiness_bounds,
                ecolor = '0.3',
                elinewidth=0.7,
                linestyle="None",
                zorder=1)
pyplot.scatter(mean_ltr, mean_bs, marker = 's', c=plot_colours,
               s=25, zorder=2)
for i, txt in enumerate(names):
    if max_interevent_times[i] > 10 and annotate_plots:
        ax.annotate(txt[:4],
                    (mean_ltr[i], mean_bs[i]),
                    fontsize=8)
ax.set_ylim([-1, 1])
ax.set_xlim([1./1000000, 1./40])
# Add B=0 linear
pyplot.plot([1./1000000, 1./40], [0, 0], linestyle='dashed', linewidth=1, c='0.5')
ax.set_xscale('log')
ax.set_xlabel('Long-term rate (events per year)')
ax.set_ylabel('B')

# Now do a bi-linear fit to the data
#indices = np.argwhere(mean_ltr > 2e-4)#.flatten()
mean_bs = np.array(mean_bs)
indices = np.flatnonzero(mean_ltr > 3e-4)
indices = indices.flatten()
#print(indices, type(indices))

indices_slow_faults = np.flatnonzero(mean_ltr <= 3e-4)
indices_slow_faults = indices_slow_faults.flatten()
# Fit fast rate faults
lf = np.polyfit(np.log10(mean_ltr[indices]),
                   mean_bs[indices], 1)
# Now force to be a flat line1
lf[0] = 0.
lf[1] = np.mean(mean_bs[indices])
std_lf = np.std(mean_bs[indices])
xvals_short = np.arange(1.5e-4, 2e-2, 1e-4)
yvals = lf[0]*np.log10(xvals_short) + lf[1]
#yvals = np.power(10, log_yvals)
pyplot.plot(xvals_short, yvals, c='0.2')
# Fit slow faults
if len(indices_slow_faults > 1):
    lf_slow = np.polyfit(np.log10(mean_ltr[indices_slow_faults]),
                         mean_bs[indices_slow_faults], 1)
    xvals_short = np.arange(1e-6, 1.5e-4, 1e-6)
    yvals = lf_slow[0]*np.log10(xvals_short) + lf_slow[1]
    #print(yvals)
    #print(xvals)
    #yvals = np.power(10, log_yvals)
    pyplot.plot(xvals_short, yvals, c='0.2')
    # Add formula for linear fits of data
print('Fits for B vs LTR')
#txt = 'Y = %.2fLog(x) + %.2f +/- %.2f' % (lf[0], lf[1], std_lf)
txt = 'Y = {:=+6.2f} +/- {:4.2f}'.format(lf[1], std_lf)
print(txt)
ax.annotate(txt, (2e-4, 0.2), fontsize=8)
#txt = 'Y = %.2fLog(x) + %.2f' % (lf_slow[0], lf_slow[1])
try:
    txt = 'Y = {:4.2f}Log(x) {:=+6.2f}'.format(lf_slow[0], lf_slow[1]) 
    print(txt)
    ax.annotate(txt, (1.5e-6, 0.75), fontsize=8)
except:
    pass
    

# Now try bilinear ODR linear fit
data = odrpack.RealData(np.log10(mean_ltr), mean_bs,
                        sx=np.log10(long_term_rate_stds), sy=burstiness_stds)
bilin = odrpack.Model(bilinear_reg_zero_slope)
odr = odrpack.ODR(data, bilin, beta0=[-3, -1.0, -4]) # array are starting values
odr.set_job(fit_type=0)
out = odr.run()
out.pprint()
a = out.beta[0]
b = out.beta[1]
hx = out.beta[2]
xvals = np.arange(1.e-6, 2e-2, 1e-6)
yrng = a*np.log10(xvals) + b #10**(b + a * xvals)
ylevel = a*hx + b #10**(b + a * hx)
print('ylevel', ylevel)
print(10**ylevel)
idx = xvals > 10**hx
yrng[idx] = (ylevel)
print('yrng', yrng)
print('hx', hx)
pyplot.plot(xvals, yrng, c='g')

# Bilinear fixed hinge
hxfix = np.log10(2e-4)
#bilin_hxfix = odrpack.Model(bilinear_reg_fix)
#odr = odrpack.ODR(data, bilin_hxfix, beta0=[-3, -1.0, -4])
#odr.set_job(fit_type=0)
#out = odr.run()
#out.pprint()
#a = out.beta[0]
#b = out.beta[1]
#a2 = out.beta[2]
#yrng = a*np.log10(xvals) + b
#ylevel = a*hxfix + b
#print('ylevel', ylevel)
#idx = xvals > 10**hxfix
#x_trans = np.log10(xvals[idx]) - hxfix
#print(x_trans)
#yrng[idx] = a2*x_trans + ylevel
#pyplot.plot(xvals, yrng, c='b')

# Bilinear fixed hinge and constant slope
#hxfix = -4
bilin_hxfix_cons_slope = odrpack.Model(bilinear_reg_fix_zero_slope)
odr = odrpack.ODR(data, bilin_hxfix_cons_slope, beta0=[-3, -1.0])
odr.set_job(fit_type=0)
out = odr.run()
out.pprint()
a = out.beta[0]
b = out.beta[1]
yrng = a*np.log10(xvals) + b
ylevel = a*hxfix + b
print('ylevel hxfix zero slope', ylevel)
print(10**ylevel)
idx = xvals > 10**hxfix
yrng[idx] = (ylevel)
print('yrng', yrng)
print('hx', hxfix)
pyplot.plot(xvals, yrng, c='r')

# Get flat mean and std values
#idx = np.argwhere(long_term_rates_T > 10**(hxfix))
#mean_b = np.mean(np.array(burstinesses))
#print('mean b', mean_b)
#mean_b_fast = np.mean(np.array(burstinesses)[idx])
#print('mean_b_fast', mean_b_fast)
#std_b_fast = np.std(np.array(burstinesses)[idx])
#txt = 'Y = {:=+6.2f} +/- {:4.2f}'.format(mean_b_fast, std_b_fast)
#print(txt)
#ax.annotate(txt, (2e-4, 0.5), fontsize=8)                                                                                                 
figname = 'burstiness_vs_lt_rate_%s.png' % fig_comment 
pyplot.savefig(figname)

#

# Plot burstiness against slip rate
pyplot.clf()
ax = pyplot.subplot(111)
#mean_bs = []
#for i, b_set in enumerate(burstinesses):
#    mean_b = np.mean(b_set)
#    mean_bs.append(mean_b)
#colours = []
#for sr in slip_rates:
#    if m <= -0.05:
#        colours.append('b')
#    elif mean_b > -0.05 and mean_b <= 0.05:
#        colours.append('g')
#    else:
#        colours.append('r')

pyplot.errorbar(slip_rates, mean_bs,
                xerr = slip_rate_bounds,
                ecolor = '0.3',
                elinewidth=0.7,
                linestyle="None",
                zorder=1)
pyplot.errorbar(slip_rates, mean_bs,
                yerr = burstiness_bounds,
                ecolor = '0.3',
                elinewidth=0.7,
                linestyle="None",
                zorder=1)
pyplot.scatter(slip_rates, mean_bs, marker = 's', c=plot_colours,
               s=25, zorder=2)
#for i, txt in enumerate(names):
#    if max_interevent_times[i] > 10 and annotate_plots:
#        ax.annotate(txt[:4],
#                    (mean_ltr[i], mean_bs[i]),
#                    fontsize=8)
ax.set_ylim([-1, 1])
ax.set_xlim([1./1000, 100])
# Add B=0 linear
pyplot.plot([1./1000, 100], [0, 0], linestyle='dashed', linewidth=1, c='0.5')
ax.set_xscale('log')
ax.set_xlabel('Slip rate (mm/yr)')
ax.set_ylabel('B')

# Now do a bi-linear fit to the data
#indices = np.argwhere(mean_ltr > 2e-4)#.flatten()
#mean_bs = np.array(mean_bs)
#indices = np.flatnonzero(slip_rates > 5)
#indices = indices.flatten()

#indices_slow_faults = np.flatnonzero(slip_rates <= 5)
#indices_slow_faults = indices_slow_faults.flatten()
# Fit fast rate faults
#lf = np.polyfit(np.log10(slip_rates[indices]),
#                   mean_bs[indices], 1)
# Now force to be a flat line1
#lf[0] = 0.
#lf[1] = np.mean(mean_bs[indices])
#std_lf = np.std(mean_bs[indices])
#xvals_short = np.arange(5, 100, 1)
#yvals = lf[0]*np.log10(xvals_short) + lf[1]
#pyplot.plot(xvals_short, yvals, c='0.2')
# Fit slow faults
#if len(indices_slow_faults > 1):
#    lf_slow = np.polyfit(np.log10(slip_rates[indices_slow_faults]),
#                         mean_bs[indices_slow_faults], 1)
#    xvals_short = np.arange(1/1000., 5, 1e-3)
#    yvals = lf_slow[0]*np.log10(xvals_short) + lf_slow[1]
#    pyplot.plot(xvals_short, yvals, c='0.2')
# Add formula for linear fits of data
#print('Fits for B vs LTR')
#txt = 'Y = {:=+6.2f} +/- {:4.2f}'.format(lf[1], std_lf)
#print(txt)
#ax.annotate(txt, (10, 0.2), fontsize=8)
#try:
#    txt = 'Y = {:4.2f}Log(x) {:=+6.2f}'.format(lf_slow[0], lf_slow[1]) 
#    print(txt)
#    ax.annotate(txt, (1.e-2, 0.75), fontsize=8)
#except:
#    pass

# Now try linear ODR linear fit
def f(B, x):
    return B[0]*x + B[1]
print(slip_rates)
print(np.log10(slip_rates))
print(slip_rate_stds)
print(np.log10(slip_rate_stds))
print(burstiness_stds)
wd = 1./np.power(burstiness_stds, 2)
print(wd)
we = 1./np.power(slip_rate_stds, 2)
print(we)
# Std dev already in log-space
data = odrpack.RealData(np.log10(slip_rates), mean_bs,
                        sx=np.sqrt(slip_rate_stds), sy=np.sqrt(burstiness_stds))

#bilin = odrpack.Model(bilinear_reg_zero_slope)
#odr = odrpack.ODR(data, bilin, beta0=[-1, -1.0, -1]) # array are starting values
linear  = odrpack.Model(f)
odr = odrpack.ODR(data, linear, beta0=[-1, -1.0,]) 
odr.set_job(fit_type=0)
out = odr.run()
out.pprint()
a = out.beta[0]
b = out.beta[1]
#hx = out.beta[2]
xvals = np.arange(1.e-4, 1e2, 1e-2)
yrng = a*np.log10(xvals) + b #10**(b + a * xvals)
#yrng = 10**(b + a * xvals)
#ylevel = a*hx + b #10**(b + a * hx)
#print('ylevel', ylevel)
#print(10**ylevel)
#idx = xvals > 10**hx
#yrng[idx] = (ylevel)
#print('yrng', yrng)
#print('hx', hx)
pyplot.plot(xvals, yrng, c='0.6')
txt = 'Y = {:4.2f}Log(x) {:=+6.2f}'.format(a, b)
print(txt)                                                                                                                                 
ax.annotate(txt, (1e0, 0.9), color='0.6') 

# Now try bilinear fixed hinge
bilin = odrpack.Model(bilinear_reg_fix_zero_slope)
#hxfix = np.log10(1e0)
odr = odrpack.ODR(data, bilin, beta0=[-1, -1.0, -1])
odr.set_job(fit_type=0)
out = odr.run()
out.pprint()
a = out.beta[0]
b = out.beta[1]   
yrng = a*np.log10(xvals) + b
ylevel = a*hxfix + b
print('ylevel hxfix zero slope', ylevel)
print(10**ylevel)
idx = xvals > 10**hxfix
yrng[idx] = (ylevel)
print('yrng', yrng)
print('hx', hxfix)
#idx = np.argwhere(slip_rates > np.power(hxfix))
#std_dev = np.std(slip_rates
pyplot.plot(xvals, yrng, c='0.2')
txt = 'Y = {:4.2f}Log(x) {:=+6.2f}, x < {:4.2f}'.format(a, b, np.power(10,hxfix))
print(txt)
ax.annotate(txt, (2e-3, 0.9), color='0.2')
txt = 'Y = {:4.2f}, x >= {:4.2f}'.format(ylevel, np.power(10,hxfix)) 
print(txt)
ax.annotate(txt, (1.2e-2, 0.8), color='0.2')

# Simple least squares fit
#lf = np.polyfit(np.log10(slip_rates),                                                                                                 
#                mean_bs, 1)
#print(lf)
#yvals = lf[0]*np.log10(xvals) + lf[1]
#pyplot.plot(xvals, yvals, c='0.2')

figname = 'burstiness_vs_slip_rate_%s.png' % fig_comment   
pyplot.savefig(figname)


# Plot memory coefficients against long term rates
pyplot.clf()
ax = pyplot.subplot(111)
mean_mems = []
for i, mem_set in enumerate(memory_coefficients):
    mean_mem = np.mean(mem_set)
#    print('Mean memory coefficient combined', mean_mem)
    mean_mems.append(mean_mem)
colours = []
for mean_mem in mean_mems:
    if mean_mem <= -0.05:
        colours.append('b')
    elif mean_mem > -0.05 and mean_mem <= 0.05:
        colours.append('g')
    else:
        colours.append('r')
pyplot.errorbar(mean_ltr, mean_mems,
                xerr = ltr_bounds,
                ecolor = '0.3',
                elinewidth=0.7,
                linestyle="None",
                zorder=1)
pyplot.errorbar(mean_ltr, mean_mems,
                yerr = memory_bounds,
                ecolor = '0.3',
                elinewidth=0.7,
                linestyle="None",
                zorder=1)
pyplot.scatter(mean_ltr, mean_mems, marker = 's', c=plot_colours,
               s=25, zorder=2)
for i, txt in enumerate(names):
    if max_interevent_times[i] > 10 and annotate_plots:
        ax.annotate(txt[:4],
                    (mean_ltr[i], mean_mems[i]),
                    fontsize=8)
#ax.set_xlim([-1, 1])
ax.set_xlim([1./1000000, 1./40])
ax.set_xscale('log')
ax.set_xlabel('Long-term rate (events per year)')
ax.set_ylabel('M')
figname = 'memory_coefficient_vs_lt_rate_%s.png' % fig_comment
pyplot.savefig(figname)

# Plot Spearman Rank coefficients against long term rates
pyplot.clf()
ax = pyplot.subplot(111)
mean_mems_L1 = []
for i, mem_set in enumerate(memory_spearman_coefficients):
    mean_mem = np.mean(mem_set)
    mean_mems_L1.append(mean_mem)
colours = []
for mean_mem in mean_mems_L1:
    if mean_mem <= -0.05:
        colours.append('b')
    elif mean_mem > -0.05 and mean_mem <= 0.05:
        colours.append('g')
    else:
        colours.append('r')
pyplot.errorbar(mean_ltr, mean_mems_L1,
                xerr = ltr_bounds,
                ecolor = '0.3',
                elinewidth=0.7,
                linestyle="None",
                zorder=1)
pyplot.errorbar(mean_ltr, mean_mems_L1,
                yerr = memory_spearman_bounds,
                elinewidth=0.7,
                ecolor = '0.3',
                linestyle="None",
                zorder=1)
pyplot.scatter(mean_ltr, mean_mems_L1, marker = 's', c=plot_colours,
               s=25, zorder=2)
for i, txt in enumerate(names):
    if max_interevent_times[i] > 10 and annotate_plots:
        ax.annotate(txt[:4],
                    (mean_ltr[i], mean_mems_L1[i]),
                    fontsize=8)
ax.set_xlim([1./1000000, 1./40])
ax.set_xscale('log')
ax.set_xlabel('Long-term rate (events per year)')
ax.set_ylabel('M (Spearman Rank)')
figname = 'memory_coefficient_Spearman_vs_lt_rate_%s.png' % fig_comment
pyplot.savefig(figname)


# Plot Spearman Rank (Lag-2) coefficients against long term rates
pyplot.clf()
ax = pyplot.subplot(111)
mean_mems_L2 = []
for i, mem_set in enumerate(memory_spearman_lag2_coef):
    mean_mem = np.mean(mem_set)
    mean_mems_L2.append(mean_mem)
colours = []
for mean_mem in mean_mems_L2:
    if mean_mem <= -0.05:
        colours.append('b')
    elif mean_mem > -0.05 and mean_mem <= 0.05:
        colours.append('g')
    else:
        colours.append('r')
pyplot.errorbar(mean_ltr, mean_mems_L2,
                xerr = ltr_bounds,
                ecolor = '0.3',
                elinewidth=0.7,
                linestyle="None",
                zorder=1)
pyplot.errorbar(mean_ltr, mean_mems_L2,
                yerr = memory_spearman_lag2_bounds,
                ecolor = '0.3',
                elinewidth=0.7,
                linestyle="None",
                zorder=1)
pyplot.scatter(mean_ltr, mean_mems_L2, marker = 's', c=plot_colours,
               s=25, zorder=2)
for i, txt in enumerate(names):
    if max_interevent_times[i] > 10 and annotate_plots:
        ax.annotate(txt[:4],
                    (mean_ltr[i], mean_mems_L2[i]),
                    fontsize=8)
ax.set_xlim([1./1000000, 1./40])
ax.set_xscale('log')
ax.set_xlabel('Long-term rate (events per year)')
ax.set_ylabel('M (Spearman Rank Lag-2)')
figname = 'memory_coefficient_Spearman_Lag2_vs_lt_rate_%s.png' % fig_comment
pyplot.savefig(figname)

# Plot Spearman rank Lag-1 against Lag-2
# Plot Spearman Rank coefficients against long term rates
pyplot.clf()
ax = pyplot.subplot(111)
colours = []
for mean_mem in mean_mems_L1:
    if mean_mem <= -0.05:
        colours.append('b')
    elif mean_mem > -0.05 and mean_mem <= 0.05:
        colours.append('g')
    else:
        colours.append('r')
pyplot.errorbar(mean_mems_L1, mean_mems_L2,
                xerr = memory_spearman_bounds,
                ecolor = '0.3',
                elinewidth=0.7,
                linestyle="None",
                zorder=1)
pyplot.errorbar(mean_mems_L1, mean_mems_L2,
                yerr = memory_spearman_lag2_bounds,
                ecolor = '0.3',
                elinewidth=0.7,
                linestyle="None",
                zorder=1)
pyplot.scatter(mean_mems_L1, mean_mems_L2, marker = 's', c=plot_colours,
               s=25, zorder=2)
for i, txt in enumerate(names):
    if max_interevent_times[i] > 10 and annotate_plots:
        ax.annotate(txt[:4],
                    (mean_mems_L1[i], mean_mems_L2[i]),
                    fontsize=8)
#ax.set_xlim([1./1000000, 1./40])
#ax.set_xscale('log')
ax.set_xlabel('M (Spearman Rank Lag-1)')
ax.set_ylabel('M (Spearman Rank Lag-2)')
figname = 'memory_coefficient_Spearman_L1_vs_L2_%s.png' % fig_comment
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
pyplot.errorbar(mean_mems, mean_bs,
                xerr = memory_bounds,
                ecolor = '0.3',
                elinewidth=0.7,
                linestyle="None",
                zorder=1)
pyplot.errorbar(mean_mems, mean_bs,
                yerr = burstiness_bounds,
                ecolor = '0.3',
                elinewidth=0.7,
                linestyle="None",
                zorder=1)
pyplot.scatter(mean_mems, mean_bs, marker = 's', c=plot_colours,
               s=25, zorder=2)
for i, txt in enumerate(names):
    if annotate_plots:
        ax.annotate(txt,
                    (mean_mems[i], mean_bs[i]),
                    fontsize=8)
ax.set_xlim([-1, 1])
ax.set_ylim([-1, 1])
# Add y = 0, x=0 lines
pyplot.plot([0,0],[-1, 1], linestyle='dashed', linewidth=1, c='0.5')
pyplot.plot([-1,1],[0, 0], linestyle='dashed', linewidth=1, c='0.5')
ax.set_ylabel('B')
ax.set_xlabel('M')
figname = 'burstiness_vs_memory_coefficient_%s.png' % fig_comment 
pyplot.savefig(figname)


# Plot COV against number of events to look at sampling biases
pyplot.clf()
ax = pyplot.subplot(111)
mean_covs = []
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
pyplot.scatter(mean_covs, num_events, marker = 's', c=plot_colours, s=25)
for i, txt in enumerate(names):
    if max_interevent_times[i] > 10 and annotate_plots:
        ax.annotate(txt[:4],
                    (mean_covs[i], num_events[i]),
                    fontsize=8)
ax.set_xlabel('COV')
ax.set_ylabel('Number of events in earthquake record')
figname = 'mean_cov_vs_number_events_%s.png' % fig_comment
pyplot.savefig(figname)

# Now plot basic statistics
pyplot.clf()
ax = pyplot.subplot(111)
pyplot.errorbar(max_interevent_times, min_interevent_times,
                yerr = min_interevent_times_bounds,
                ecolor = '0.3',
                elinewidth=0.7,
                linestyle="None",
                zorder=1)
pyplot.errorbar(max_interevent_times, min_interevent_times,
                xerr = max_interevent_times_bounds,
                ecolor = '0.3',
                elinewidth=0.7,
                linestyle="None",
                zorder=1)
pyplot.scatter(max_interevent_times, min_interevent_times,
               marker = 's', c=colours, s=25, zorder=2)
ax.set_xlabel('Maximum interevent time')
ax.set_ylabel('Minimum interevent time') 
ax.set_xscale('log')
ax.set_yscale('log')
# Label low-slip rate faults
for i, txt in enumerate(names):
    if max_interevent_times[i] > 10 and annotate_plots:
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
                ecolor = '0.3',
                elinewidth=0.7,
                linestyle="None",
                zorder=1)
pyplot.errorbar(max_interevent_times, min_paired_interevent_times,
                xerr = max_interevent_times_bounds,
                ecolor = '0.3',
                elinewidth=0.7,
                linestyle="None",
                zorder=1)
pyplot.scatter(max_interevent_times, min_paired_interevent_times,
               marker = 's', c=colours, s=25, zorder=2)
ax.set_xlabel('Maximum interevent time')
ax.set_ylabel('Minimum interevent time \n(mean of two shortest consecutive interevent times)')
ax.set_xscale('log')
ax.set_yscale('log') 
# Label low-slip rate faults
for i, txt in enumerate(names):
    if max_interevent_times[i] > 10 and annotate_plots:
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
#def linear_func(B, x):
#    return B[0]*x + B[1]
#linear_model = Model(linear_func)
#data = RealData(np.log10(max_interevent_times[indices]),
#                np.log10(min_paired_interevent_times[indices]),
#                sx = np.log10(std_max_interevent_times[indices]),
#                sy = np.log10(std_min_paired_interevent_times[indices]))
# Set up ODR with the model and data
#odr = ODR(data, linear_model, beta0=[1., 1.])
# Run the regression.
#out = odr.run()
#print('Regression output')
#out.pprint()

#log_y_fit = linear_func(out.beta, np.log10(xvals_short))
#y_fit = np.power(10, log_y_fit)
#popt, pcov = curve_fit(target_func, max_interevent_times, min_paired_interevent_times)
#print('popt', popt)
#pyplot.plot(xvals_short, y_fit, '--')
#txt = 'Log(Y) =  %.2fLog(x) + %.2f' % (out.beta[0], out.beta[1])
#print(txt)
#ax.annotate(txt, (100, 40000))
figname = 'min_pair_vs_max_interevent_time_%s.png' % fig_comment
pyplot.savefig(figname)

# Similar plots, against long term rates
pyplot.clf()
ax = pyplot.subplot(111)  
pyplot.errorbar(mean_ltr, min_interevent_times,
                yerr = min_interevent_times_bounds,
                ecolor = '0.3',
                elinewidth=0.7,
                linestyle="None",
                zorder=1)
pyplot.errorbar(mean_ltr, min_interevent_times,
                xerr = ltr_bounds,
                ecolor = '0.3',
                elinewidth=0.7,
                linestyle="None",
                zorder=1)
pyplot.scatter(mean_ltr, min_interevent_times,
               marker='s', c=colours, s=25, zorder=2)
ax.set_xlabel('Long-term rate')
ax.set_ylabel('Minimum interevent time')
ax.set_xscale('log')
ax.set_yscale('log') 
# Label low-slip rate faults
for i, txt in enumerate(names):
    if max_interevent_times[i] > 10 and annotate_plots:
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
                ecolor = '0.3',
                elinewidth=0.7,
                linestyle="None",
                zorder=1)
pyplot.errorbar(mean_ltr, min_paired_interevent_times,
                xerr = ltr_bounds,
                ecolor = '0.3',
                elinewidth=0.7,
                linestyle="None",
                zorder=1)
pyplot.scatter(mean_ltr, min_paired_interevent_times,
               marker='s', c=colours, s=25, zorder=2)
#c='0.1', s=25)
ax.set_xlabel('Long-term rate')
ax.set_ylabel('Minimum interevent time \n(mean of two shortest consecutive interevent times)')
ax.set_xscale('log')
ax.set_yscale('log') 
# Label low-slip rate faults
for i, txt in enumerate(names):
    if max_interevent_times[i] > 10 and annotate_plots:
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
                ecolor = '0.3',
                elinewidth=0.7,
                linestyle="None",
                zorder=1)
pyplot.errorbar(mean_ltr, max_interevent_times,
                xerr = ltr_bounds,
                ecolor = '0.3',
                elinewidth=0.7,
                linestyle="None",
                zorder=1)
pyplot.scatter(mean_ltr, max_interevent_times,
               marker='s', c=plot_colours, s=25, zorder=2)
#c='0.1', s=25)
ax.set_xlabel('Long-term rate')
ax.set_ylabel('Maximum interevent time')
ax.set_xscale('log')
ax.set_yscale('log') 
# Label low-slip rate faults
for i, txt in enumerate(names):
    if max_interevent_times[i] > 10 and annotate_plots:
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
ax.annotate(txt, (1e-4, 100000))
figname = 'max_interevent_time_vs_ltr_%s.png' % fig_comment
pyplot.savefig(figname)


# Now plot ratios against long term rates
pyplot.clf()
ax = pyplot.subplot(111)  
pyplot.errorbar(mean_ltr, ratio_min_pair_max,
                yerr = ratio_min_pair_max_bounds,
                ecolor = '0.3',
                elinewidth=0.7,
                linestyle="None",
                zorder=1)
pyplot.errorbar(mean_ltr, ratio_min_pair_max,
                xerr = ltr_bounds,
                ecolor = '0.3',
                elinewidth=0.7,
                linestyle="None",
                zorder=1)
pyplot.scatter(mean_ltr, ratio_min_pair_max,
               marker='s', c=plot_colours, s=25, zorder=2)
ax.set_xlabel('Long-term rate')
ax.set_ylabel('Minimum pair interevent time: maximum interevent time')
ax.set_xscale('log')
ax.set_yscale('log') 
# Label low-slip rate faults
for i, txt in enumerate(names):
    if max_interevent_times[i] > 10 and annotate_plots:
        ax.annotate(txt[:4],
                    (mean_ltr[i], ratio_min_pair_max[i]),
                    fontsize=8)

#log_yvals = lf[0]*np.log10(xvals) + lf[1]
#yvals = np.power(10, log_yvals)
#pyplot.plot(xvals, yvals)

# Linear fit high and low long term rate data separately
indices = np.argwhere(mean_ltr > 4e-4).flatten()
indices_slow_faults = np.argwhere(mean_ltr <= 4e-4).flatten()
lf = np.polyfit(np.log10(mean_ltr[indices]),
                                np.log10(ratio_min_pair_max[indices]), 1)
xvals_short = np.arange(2e-4, 5e-2, 1e-4)
log_yvals = lf[0]*np.log10(xvals_short) + lf[1]
yvals = np.power(10, log_yvals)
pyplot.plot(xvals_short, yvals, c='k')
# Add formula for linear fit to low-end of data
txt = 'Log(Y) = %.2fLog(x) + %.2f' % (lf[0], lf[1])
print(txt)
ax.annotate(txt, (5e-4, 1e-2))

# Slow long-term rates
print('At if statement')
if len(indices_slow_faults) > 0:
    print('Plotting slow faults')
    lf = np.polyfit(np.log10(mean_ltr[indices_slow_faults]),
                    np.log10(ratio_min_pair_max[indices_slow_faults]), 1)
    xvals_short = np.arange(2e-6, 4e-4, 1e-6)
    log_yvals = lf[0]*np.log10(xvals_short) + lf[1]
    yvals = np.power(10, log_yvals)
    pyplot.plot(xvals_short, yvals, c='k')
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
                ecolor = '0.3',
                elinewidth=0.7,
                linestyle="None",
                zorder=1)
pyplot.errorbar(mean_ltr, ratio_min_max,
                xerr = ltr_bounds,
                ecolor = '0.3',
                elinewidth=0.7,
                linestyle="None", zorder=1)
pyplot.scatter(mean_ltr, ratio_min_max,
               marker = 's', c=plot_colours, s=25, zorder=2)
ax.set_xlabel('Long-term rate')
ax.set_ylabel('Minimum interevent time: maximum interevent time')
ax.set_xscale('log')
ax.set_yscale('log') 
# Label low-slip rate faults
for i, txt in enumerate(names):
    if max_interevent_times[i] > 10 and annotate_plots:
        ax.annotate(txt[:4],
                    (mean_ltr[i], ratio_min_max[i]),
                    fontsize=8)

#log_yvals = lf[0]*np.log10(xvals) + lf[1]
#yvals = np.power(10, log_yvals)
#pyplot.plot(xvals, yvals)

# Linear fit only bottom end of data
indices = np.argwhere(mean_ltr > 9e-5).flatten()
indices_slow_faults = np.argwhere(mean_ltr <= 9e-5).flatten()
lf = np.polyfit(np.log10(mean_ltr[indices]),
                                np.log10(ratio_min_max[indices]), 1)
# Now just plot as constant mean value
lf[0] = 0
lf[1] = np.mean(np.log10(ratio_min_max[indices]))
xvals_short = np.arange(3.46e-5, 1e-2, 1e-4)
log_yvals = lf[0]*np.log10(xvals_short) + lf[1]
yvals = np.power(10, log_yvals)
pyplot.plot(xvals_short, yvals, c='k')
# Add formula for linear fit to low-end of data
#txt = 'Log(Y) = %.2fLog(x) %+.2f' % (lf[0], lf[1])
txt = 'Log(Y) = {:4.2f}Log(x) {:=+6.2f}'.format(lf[0], lf[1])
print(txt)
ax.annotate(txt, (1e-4, 1e-3))
# Slow long-term rates
if len(indices_slow_faults) > 0:
    lf = np.polyfit(np.log10(mean_ltr[indices_slow_faults]),
                    np.log10(ratio_min_max[indices_slow_faults]), 1)
    xvals_short = np.arange(2e-6, 3.47e-5, 1e-6)
    log_yvals = lf[0]*np.log10(xvals_short) + lf[1]
    yvals = np.power(10, log_yvals)
    pyplot.plot(xvals_short, yvals, c='k')
    # Add formula for linear fit to low-end of data
#    txt = 'Log(Y) = %.2fLog(x) %+.2f' % (lf[0], lf[1])
    txt = 'Log(Y) = {:4.2f} {:=+6.2f}'.format(lf[0], lf[1])
    print(txt)
    ax.annotate(txt, (3e-6, 8e-1))

figname = 'min_max_ratio_vs_ltr_%s.png' % fig_comment
pyplot.savefig(figname)

#############################################
# Make multipanel figure plot
pyplot.clf()
fig = pyplot.figure(1)
# set up subplot grid
gridspec.GridSpec(3, 2)
#First plot
pyplot.subplot2grid((3, 2), (0,0), colspan=1, rowspan=1)
ax = pyplot.gca()
# Plot burstiness against mean ltr
pyplot.errorbar(mean_ltr, mean_bs,
                xerr = ltr_bounds,
                ecolor = '0.3',
                elinewidth=0.5,
                linestyle="None",
                zorder=1)
pyplot.errorbar(mean_ltr, mean_bs,
                yerr = burstiness_bounds,
                ecolor = '0.3',
                elinewidth=0.5,
                linestyle="None",
                zorder=1)
pyplot.scatter(mean_ltr, mean_bs, marker = 's', c=plot_colours,
               s=18, zorder=2)
ax.set_ylim([-1, 1])
ax.set_xlim([1./1000000, 1./40])
pyplot.plot([1./1000000, 1./40], [0, 0], linestyle='dashed', linewidth=1, c='0.5')
ax.set_xscale('log')
ax.set_xlabel('Long-term rate (events per year)', fontsize=10)
ax.set_ylabel('B', fontsize=10)
# Add a legend using some dummy data
line1 = ax.scatter([1], [100], marker = 's', c = 'r', s=18)
line2 = ax.scatter([1], [100], marker = 's', c = 'g', s=18)
line3 = ax.scatter([1], [100], marker = 's', c = 'b', s=18)
pyplot.legend((line1, line2, line3), ('Normal', 'Strike slip', 'Reverse'))

# Now do a bi-linear fit to the data
#mean_bs = np.array(mean_bs)
#indices = np.flatnonzero(mean_ltr > 3e-4)
#indices = indices.flatten()
#indices_slow_faults = np.flatnonzero(mean_ltr <= 3e-4)
#indices_slow_faults = indices_slow_faults.flatten()
# Fit fast rate faults
#lf = np.polyfit(np.log10(mean_ltr[indices]),
#                   mean_bs[indices], 1)
# Now force to be a flat line1
#lf[0] = 0.
#lf[1] = np.mean(mean_bs[indices])
#std_lf = np.std(mean_bs[indices])
#xvals_short = np.arange(1.5e-4, 2e-2, 1e-4)
#yvals = lf[0]*np.log10(xvals_short) + lf[1]
#pyplot.plot(xvals_short, yvals, c='0.4')
# Fit slow faults
#lf_slow = np.polyfit(np.log10(mean_ltr[indices_slow_faults]),
#                   mean_bs[indices_slow_faults], 1)
#xvals_short = np.arange(1e-6, 1.5e-4, 1e-6)
#yvals = lf_slow[0]*np.log10(xvals_short) + lf_slow[1]
#pyplot.plot(xvals_short, yvals, c='0.4')
# Add formula for linear fits of data
#print('Fits for B vs LTR')
#txt = 'Y = {:=+6.2f} +/- {:4.2f}'.format(lf[1], std_lf)
#print(txt)
#ax.annotate(txt, (2e-4, 0.2), fontsize=8)
#txt = 'Y = {:4.2f}Log(x) {:=+6.2f}'.format(lf_slow[0], lf_slow[1]) 
#print(txt)
#ax.annotate(txt, (1.6e-6, -0.75), fontsize=8)

# Bilinear fixed hinge and constant slope ODR
hxfix = np.log10(2e-4)
bilin_hxfix_cons_slope = odrpack.Model(bilinear_reg_fix_zero_slope)
data = odrpack.RealData(np.log10(mean_ltr), mean_bs,
                        sx=np.log10(long_term_rate_stds), sy=burstiness_stds)
odr = odrpack.ODR(data, bilin_hxfix_cons_slope, beta0=[-3, -1.0])
odr.set_job(fit_type=0)
out = odr.run()
out.pprint()
a = out.beta[0]
b = out.beta[1]
xvals = np.arange(1e-6, 2e-2, 1e-5)
yrng = a*np.log10(xvals) + b
ylevel = a*hxfix + b
print('ylevel hxfix zero slope', ylevel)
print(10**ylevel)
idx = xvals > 10**hxfix
yrng[idx] = (ylevel)
print('yrng', yrng)
print('hx', hxfix)
pyplot.plot(xvals, yrng, c='0.4')
txt = 'y = {:4.2f}Log(x) {:=+6.2f}, x < {:3.1E}'.format(a, b, np.power(10, hxfix))
ax.annotate(txt, (1.5e-6, -0.85), fontsize=8) 
txt = 'y = {:=+4.2f}, x >= {:3.1E}'.format(ylevel, np.power(10, hxfix))
ax.annotate(txt, (1.5e-6, -0.95), fontsize=8)

ax.annotate('a)', (-0.23, 0.98), xycoords = 'axes fraction', fontsize=10)

# Add second plot (Memory vs LTR)
pyplot.subplot2grid((3, 2), (0,1), colspan=1, rowspan=1) 
ax = pyplot.gca()
mean_mems = []
for i, mem_set in enumerate(memory_coefficients):
    mean_mem = np.mean(mem_set)
#    print('Mean memory coefficient combined', mean_mem)
    mean_mems.append(mean_mem)
pyplot.errorbar(mean_ltr, mean_mems,
                xerr = ltr_bounds,
                ecolor = '0.3',
                elinewidth=0.5,
                linestyle="None",
                zorder=1)
pyplot.errorbar(mean_ltr, mean_mems,
                yerr = memory_bounds,
                ecolor = '0.3',
                elinewidth=0.5,
                linestyle="None",
                zorder=1)
pyplot.scatter(mean_ltr, mean_mems, marker = 's', c=plot_colours,
               s=18, zorder=2)
for i, txt in enumerate(names):
    if max_interevent_times[i] > 10 and annotate_plots:
        ax.annotate(txt[:4],
                    (mean_ltr[i], mean_mems[i]),
                    fontsize=8)
#ax.set_xlim([-1, 1])
ax.set_xlim([1./1000000, 1./40])
ax.set_ylim([-1, 1]) 
pyplot.plot([1./1000000, 1./40], [0, 0], linestyle='dashed', linewidth=1, c='0.5')
ax.set_xscale('log')
ax.set_xlabel('Long-term rate (events per year)', fontsize=10)
ax.set_ylabel('M', fontsize=10)

# Bilinear fixed hinge and constant slope ODR
hxfix = np.log10(2e-4)
bilin_hxfix_cons_slope = odrpack.Model(bilinear_reg_fix_zero_slope)
data = odrpack.RealData(np.log10(mean_ltr), mean_mems,
                        sx=np.log10(long_term_rate_stds), sy=memory_stds)
odr = odrpack.ODR(data, bilin_hxfix_cons_slope, beta0=[-3, -1.0])
odr.set_job(fit_type=0)
out = odr.run()
out.pprint()
a = out.beta[0]
b = out.beta[1]
yrng = a*np.log10(xvals) + b
ylevel = a*hxfix + b
print('ylevel hxfix zero slope', ylevel)
print(ylevel)
idx = xvals > 10**hxfix
yrng[idx] = (ylevel)
print('yrng', yrng)
print('hx', hxfix)
pyplot.plot(xvals, yrng, c='0.4')
txt = 'y = {:4.2f}Log(x) {:=+6.2f}, x < {:3.1E}'.format(a, b, np.power(10, hxfix))
ax.annotate(txt, (1.5e-6, -0.85), fontsize=8) 
txt = 'y = {:4.2f}, x >= {:3.1E}'.format(ylevel, np.power(10, hxfix))
ax.annotate(txt, (1.5e-6, -0.95), fontsize=8)

# Now try inverting fo hinge point
#bilin = odrpack.Model(bilinear_reg_zero_slope)
#odr = odrpack.ODR(data, bilin, beta0=[-3, -1.0, -4]) # array are starting values
#odr.set_job(fit_type=0)
#out = odr.run()
#out.pprint()
#a = out.beta[0]
#b = out.beta[1]
#hx = out.beta[2]
#xvals = np.arange(1.e-6, 2e-2, 1e-6)
#yrng = a*np.log10(xvals) + b
#ylevel = a*hx + b
#print('ylevel', ylevel)
#print(ylevel)
#idx = xvals > 10**hx
#yrng[idx] = (ylevel)
#print('yrng', yrng)
#print('hx', hx)
#pyplot.plot(xvals, yrng, c='0.4')
#txt = 'Log(y) = {:4.2f}Log(x) {:=+6.2f}, x < {:3.1E}'.format(a, b, np.power(10, hx))
#ax.annotate(txt, (1.5e-6, 0.9), fontsize=8)
#txt = 'y = {:4.2f}, x >= {:3.1E}'.format(ylevel, np.power(10, hx))
#ax.annotate(txt, (1.3e-6, 0.8), fontsize=8)

# Linear ODR fit
#linear  = odrpack.Model(f)
#odr = odrpack.ODR(data, linear, beta0=[-1, -1.0,]) 
#odr.set_job(fit_type=0)
#out = odr.run()
#out.pprint()
#a = out.beta[0]
#b = out.beta[1]
#xvals = np.arange(1.e-4, 1e2, 1e-2)
#yrng = a*np.log10(xvals) + b #10**(b + a * xvals)
#pyplot.plot(xvals, yrng, c='0.6')
#txt = 'Y = {:4.2f}Log(x) {:=+6.2f}'.format(a, b)
#print(txt)                                                                                                                                 
#ax.annotate(txt, (1.3e-6, 0.6), color='0.6')

ax.annotate('b)', (-0.23, 0.98), xycoords = 'axes fraction', fontsize=10)

# Now do a bi-linear fit to the data
#mean_mems = np.array(mean_mems)
#indices = np.flatnonzero(mean_ltr > 3e-4)
#indices = indices.flatten()
#indices_slow_faults = np.flatnonzero(mean_ltr <= 3e-4)
#indices_slow_faults = indices_slow_faults.flatten()
# Fit fast rate faults
#lf = np.polyfit(np.log10(mean_ltr[indices]),
#                   mean_mems[indices], 1)
# Now force to be a flat line1
#lf[0] = 0.
#lf[1] = np.mean(mean_mems[indices])
#std_lf = np.std(mean_mems[indices])
#xvals_short = np.arange(1.5e-4, 2e-2, 1e-4)
#yvals = lf[0]*np.log10(xvals_short) + lf[1]
#pyplot.plot(xvals_short, yvals, c='0.4')
# Fit slow faults
#lf_slow = np.polyfit(np.log10(mean_ltr[indices_slow_faults]),
#                   mean_mems[indices_slow_faults], 1)
#xvals_short = np.arange(1e-6, 1.5e-4, 1e-6)
#yvals = lf_slow[0]*np.log10(xvals_short) + lf_slow[1]
#pyplot.plot(xvals_short, yvals, c='0.4')
# Add formula for linear fits of data
#print('Fits for B vs LTR')
#txt = 'Y = {:=+6.2f} +/- {:4.2f}'.format(lf[1], std_lf)
#print(txt)
#ax.annotate(txt, (2e-4, 0.2), fontsize=8)
#txt = 'Y = {:4.2f}Log(x) {:=+6.2f}'.format(lf_slow[0], lf_slow[1]) 
#print(txt)
#ax.annotate(txt, (1.6e-6, -0.75), fontsize=8)


# Add third plot
pyplot.subplot2grid((3, 2), (1,0), colspan=1, rowspan=1)
ax = pyplot.gca()
pyplot.errorbar(mean_mems, mean_bs,
                xerr = memory_bounds,
                ecolor = '0.3',
                elinewidth=0.5,
                linestyle="None",
                zorder=1)
pyplot.errorbar(mean_mems, mean_bs,
                yerr = burstiness_bounds,
                ecolor = '0.3',
                elinewidth=0.5,
                linestyle="None",
                zorder=1)
pyplot.scatter(mean_mems, mean_bs, marker = 's', c=plot_colours,
               s=18, zorder=2)
for i, txt in enumerate(names):
    if max_interevent_times[i] > 10 and annotate_plots:
        ax.annotate(txt[:7],
                    (mean_mems[i], mean_bs[i]),
                    fontsize=8)
ax.set_xlim([-1, 1])
ax.set_ylim([-1, 1])
# Add y = 0, x=0 lines
pyplot.plot([0,0],[-1, 1], linestyle='dashed', linewidth=1, c='0.5')
pyplot.plot([-1,1],[0, 0], linestyle='dashed', linewidth=1, c='0.5')
#ax.set_yscale('log')
ax.set_ylabel('B', fontsize=10)
ax.set_xlabel('M', fontsize=10)
ax.annotate('c)', (-0.23, 0.98), xycoords = 'axes fraction', fontsize=10)

# Add fourth plot
pyplot.subplot2grid((3, 2), (1,1), colspan=1, rowspan=1)
ax = pyplot.gca()
pyplot.errorbar(mean_ltr, max_interevent_times,
                yerr = max_interevent_times_bounds,
                ecolor = '0.3',
                elinewidth=0.5,
                linestyle="None",
                zorder=1)
pyplot.errorbar(mean_ltr, max_interevent_times,
                xerr = ltr_bounds,
                ecolor = '0.3',
                elinewidth=0.5,
                linestyle="None",
                zorder=1)
pyplot.scatter(mean_ltr, max_interevent_times,
               marker='s', c=plot_colours, s=18, zorder=2)
#c='0.1', s=25)
ax.set_xlabel('Long-term rate (events per year)', fontsize=10)
ax.set_ylabel(r'$\tau_{max}$', fontsize=10)
ax.set_xscale('log')
ax.set_yscale('log')
# Label low-slip rate faults
for i, txt in enumerate(names):
    if max_interevent_times[i] > 10 and annotate_plots:
        ax.annotate(txt[:4],
                    (mean_ltr[i], max_interevent_times[i]),
                    fontsize=8)
indices = np.argwhere(mean_ltr > 2e-10).flatten() # All data for now
lf = np.polyfit(np.log10(mean_ltr[indices]),
                np.log10(max_interevent_times[indices]), 1)
xvals_short = np.arange(2e-6, 2e-2, 1e-6)
log_yvals = lf[0]*np.log10(xvals_short) + lf[1]
yvals = np.power(10, log_yvals)
pyplot.plot(xvals_short, yvals, c='0.4')
# Add formula for linear fit to low-end of data
txt = 'Log(y) = %.2fLog(x) + %.2f' % (lf[0], lf[1])
print(txt)
ax.annotate(txt, (1e-5, 2000000), fontsize=8)
ax.annotate('d)', (-0.23, 0.98), xycoords = 'axes fraction', fontsize=10)

# Add fifth plot
pyplot.subplot2grid((3, 2), (2,0), colspan=1, rowspan=1)
ax = pyplot.gca()
pyplot.errorbar(mean_ltr, ratio_min_max,
                yerr = ratio_min_max_bounds,
                ecolor = '0.3',
                elinewidth=0.5,
                linestyle="None",
                zorder=1)
pyplot.errorbar(mean_ltr, ratio_min_max,
                xerr = ltr_bounds,
                ecolor = '0.3',
                elinewidth=0.5,
                linestyle="None",
                zorder=1)
pyplot.scatter(mean_ltr, ratio_min_max,
               marker='s', c=plot_colours, s=18, zorder=2)
ax.set_xlabel('Long-term rate (events per year)', fontsize=10)
ax.set_ylabel(r'$\tau_{min}$ / $\tau_{max}$', fontsize=10)
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlim([1e-6, 2e-2])
ax.set_ylim([5e-4, 2])
# Linear fit only bottom end of data
#indices = np.argwhere(mean_ltr > 9e-5).flatten()
#indices_slow_faults = np.argwhere(mean_ltr <= 9e-5).flatten()
#lf = np.polyfit(np.log10(mean_ltr[indices]),
#                                np.log10(ratio_min_max[indices]), 1)
# Now just plot as constant mean value
#lf[0] = 0
#lf[1] = np.mean(np.log10(ratio_min_max[indices]))
#std_lf = np.std(np.log10(ratio_min_max[indices]))
#xvals_short = np.arange(3.46e-5, 1e-2, 1e-4)
#log_yvals = lf[0]*np.log10(xvals_short) + lf[1]
#yvals = np.power(10, log_yvals)
#pyplot.plot(xvals_short, yvals, c='0.4')
# Add formula for linear fit to low-end of data
#txt = 'Log(Y) = %.2fLog(x) %+.2f' % (lf[0], lf[1])
#txt = 'Log(y) = {:=+6.2f} +/- {:4.2f}'.format(lf[1], std_lf)
#print(txt)
#ax.annotate(txt, (6e-5, 1e-3), fontsize=8)
# Slow long-term rates
#if len(indices_slow_faults) > 0:
#    lf = np.polyfit(np.log10(mean_ltr[indices_slow_faults]),
#                    np.log10(ratio_min_max[indices_slow_faults]), 1)
#    xvals_short = np.arange(2e-6, 3.47e-5, 1e-6)
#    log_yvals = lf[0]*np.log10(xvals_short) + lf[1]
#    yvals = np.power(10, log_yvals)
#    pyplot.plot(xvals_short, yvals, c='0.4')
#    # Add formula for linear fit to low-end of data
#    txt = 'Log(Y) = %.2fLog(x) %+.2f' % (lf[0], lf[1])
#    txt = 'Log(y) = {:4.2f}Log(x) {:=+6.2f}'.format(lf[0], lf[1])
#    print(txt)
#    ax.annotate(txt, (2e-6, 1.e-0), fontsize=8)

# Bilinear fixed hinge and constant slope ODR
hxfix = np.log10(2e-4)
bilin_hxfix_cons_slope = odrpack.Model(bilinear_reg_fix_zero_slope)
data = odrpack.RealData(np.log10(mean_ltr), np.log10(ratio_min_max),
                        sx=np.log10(np.sqrt(long_term_rate_stds)), sy=np.log10(np.sqrt(std_ratio_min_max)))
odr = odrpack.ODR(data, bilin_hxfix_cons_slope, beta0=[-3, -1.0])
odr.set_job(fit_type=0)
out = odr.run()
out.pprint()
a = out.beta[0]
b = out.beta[1]
log_y = a*np.log10(xvals) + b
yrng = np.power(10, log_y)
ylevel = np.power(10, (a*hxfix + b)) 
print('ylevel hxfix zero slope', ylevel)
print(10**ylevel)
idx = xvals > 10**hxfix
yrng[idx] = (ylevel)
print('yrng', yrng)
print('hx', hxfix)
#pyplot.plot(xvals, yrng, c='0.4')
#txt = 'Log(y) = {:4.2f}Log(x) {:=+6.2f}, x < {:3.1E}'.format(a, b, np.power(10, hxfix))
#ax.annotate(txt, (1.5e-6, 1.05), fontsize=8) 
#txt = 'y = {:4.2f}, x >= {:3.1E}'.format(ylevel, np.power(10, hxfix))
#ax.annotate(txt, (1.3e-6, 0.6), fontsize=8)

# Now try inverting fo hinge point
bilin = odrpack.Model(bilinear_reg_zero_slope)
odr = odrpack.ODR(data, bilin, beta0=[-3, -1.0, -4]) # array are starting values
odr.set_job(fit_type=0)
out = odr.run()
out.pprint()
a = out.beta[0]
b = out.beta[1]
hx = out.beta[2]
xvals = np.arange(1.e-6, 2e-2, 1e-6)
log_y = a*np.log10(xvals) + b
yrng = np.power(10, log_y)
ylevel = np.power(10, a*hx + b) #10**(b + a * hx)
print('ylevel', ylevel)
print(10**ylevel)
idx = xvals > 10**hx
yrng[idx] = (ylevel)
print('yrng', yrng)
print('hx', hx)
pyplot.plot(xvals, yrng, c='0.4')
txt = 'Log(y) = {:4.2f}Log(x) {:=+6.2f}, x < {:3.1E}'.format(a, b, np.power(10, hx))
ax.annotate(txt, (1.5e-6, 1.08), fontsize=8)
txt = 'y = {:4.2f}, x >= {:3.1E}'.format(ylevel, np.power(10, hx))
ax.annotate(txt, (1.5e-6, 0.6), fontsize=8) 

ax.annotate('e)', (-0.23, 0.98), xycoords = 'axes fraction', fontsize=10)

# Add sixth plot
pyplot.subplot2grid((3, 2), (2,1), colspan=1, rowspan=1)
ax = pyplot.gca()
pyplot.errorbar(mean_ltr, ratio_min_pair_max,
                yerr = ratio_min_pair_max_bounds,
                ecolor = '0.3',
                elinewidth=0.5,
                linestyle="None",
                zorder=1)
pyplot.errorbar(mean_ltr, ratio_min_pair_max,
                xerr = ltr_bounds,
                ecolor = '0.3',
                elinewidth=0.5,
                linestyle="None",
                zorder=1)
pyplot.scatter(mean_ltr, ratio_min_pair_max,
               marker='s', c=plot_colours, s=18, zorder=2)
ax.set_xlabel('Long-term rate (events per year)', fontsize=10)
ax.set_ylabel(r'$\bar{\tau}_{min(p)}$ / $\tau_{max}$', fontsize=10)
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlim([1e-6, 2e-2])
ax.set_ylim([5e-4, 2])
# Now just plot as constant mean value
#lf[0] = 0
#lf[1] = np.mean(np.log10(ratio_min_pair_max[indices]))
#std_lf = np.std(np.log10(ratio_min_pair_max[indices]))
#xvals_short = np.arange(3.46e-5, 1e-2, 1e-4)
#log_yvals = lf[0]*np.log10(xvals_short) + lf[1]
#yvals = np.power(10, log_yvals)
#pyplot.plot(xvals_short, yvals, c='0.4')
# Add formula for linear fit to low-end of data
#txt = 'Log(Y) = %.2fLog(x) %+.2f' % (lf[0], lf[1])
#txt = 'Log(y) = {:=+6.2f} +/- {:4.2f}'.format(lf[1], std_lf)
#print(txt)
#ax.annotate(txt, (6e-5, 1e-3), fontsize=8)

# Slow long-term rates
#if len(indices_slow_faults) > 0:
#    lf = np.polyfit(np.log10(mean_ltr[indices_slow_faults]),
#                    np.log10(ratio_min_pair_max[indices_slow_faults]), 1)
#    xvals_short = np.arange(2e-6, 3.47e-5, 1e-6)
#    log_yvals = lf[0]*np.log10(xvals_short) + lf[1]
#    yvals = np.power(10, log_yvals)
#    pyplot.plot(xvals_short, yvals, c='0.4')
#    # Add formula for linear fit to low-end of data
#    txt = 'Log(Y) = %.2fLog(x) %+.2f' % (lf[0], lf[1])
#    txt = 'Log(y) = {:4.2f}Log(x) \n{:=+6.2f}'.format(lf[0], lf[1])
#    print(txt)
#    ax.annotate(txt, (2e-6, 8.e-1), fontsize=8)

# Bilinear fixed hinge and constant slope ODR
hxfix = np.log10(2e-4)
bilin_hxfix_cons_slope = odrpack.Model(bilinear_reg_fix_zero_slope)
data = odrpack.RealData(np.log10(mean_ltr), np.log10(ratio_min_pair_max),
                        sx=np.log10(long_term_rate_stds), sy=np.log10(std_ratio_min_pair_max))
odr = odrpack.ODR(data, bilin_hxfix_cons_slope, beta0=[-3, -1.0])
odr.set_job(fit_type=0)
out = odr.run()
out.pprint()
a = out.beta[0]
b = out.beta[1]
log_y = a*np.log10(xvals) + b
yrng = np.power(10, log_y) 
ylevel = np.power(10, (a*hxfix + b))
print('ylevel hxfix zero slope', ylevel)
print(10**ylevel)
idx = xvals > 10**hxfix
yrng[idx] = (ylevel)
print('yrng', yrng)
print('hx', hxfix)
#pyplot.plot(xvals, yrng, c='0.4')
#txt = 'Log(y) = {:4.2f}Log(x) {:=+6.2f}, x < {:3.1E}'.format(a, b, np.power(10, hxfix))
#ax.annotate(txt, (1e-5, 5.e-3), fontsize=8) 
#txt = 'y = {:4.2f}, x >= {:3.1E}'.format(ylevel, np.power(10, hxfix))
#ax.annotate(txt, (1e-5, 1.e-3), fontsize=8)

# Now try inverting fo hinge point
bilin = odrpack.Model(bilinear_reg_zero_slope)
odr = odrpack.ODR(data, bilin, beta0=[-3, -1.0, -4]) # array are starting values
odr.set_job(fit_type=0)
out = odr.run()
out.pprint()
a = out.beta[0]
b = out.beta[1]
hx = out.beta[2]
xvals = np.arange(1.e-6, 2e-2, 1e-6)
log_y = a*np.log10(xvals) + b
yrng = np.power(10, log_y)
ylevel = np.power(10, a*hx + b) #10**(b + a * hx)
print('ylevel', ylevel)
print(10**ylevel)
idx = xvals > 10**hx
yrng[idx] = (ylevel)
print('yrng', yrng)
print('hx', hx)
pyplot.plot(xvals, yrng, c='0.4')
txt = 'Log(y) = {:4.2f}Log(x) {:=+6.2f}, x < {:3.1E}'.format(a, b, np.power(10, hx))
ax.annotate(txt, (2e-6, 2.e-3), fontsize=8)
txt = 'y = {:4.2f}, x >= {:3.1E}'.format(ylevel, np.power(10, hx))
ax.annotate(txt, (2e-6, 1.e-3), fontsize=8) 

ax.annotate('f)', (-0.23, 0.98), xycoords = 'axes fraction', fontsize=10)


fig.tight_layout(pad=1.2, w_pad=1.0, h_pad=-1)
fig.set_size_inches(w=7.5,h=10.5)
figname = 'combined_plot_%s.png' % fig_comment
pyplot.savefig(figname)


# Plot M-B phase diagram with stuff over the top
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
pyplot.errorbar(mean_mems, mean_bs,
                xerr = memory_bounds,
                ecolor = '0.8',
                elinewidth=0.7,
                linestyle="None",
                zorder=1)
pyplot.errorbar(mean_mems, mean_bs,
                yerr = burstiness_bounds,
                ecolor = '0.8',
                elinewidth=0.7,
                linestyle="None",
                zorder=1)
pyplot.scatter(mean_mems, mean_bs, marker = 's', c=plot_colours,
               s=25, zorder=2)
texts = []
for i, txt in enumerate(names):
    sp = camel_case_split(txt)
    # Handle special cases of two word fault names
    if sp[0] == 'San' or sp[0] == 'Dead' or sp[0] == 'Loma':
        flt_txt = sp[0] + ' '  + sp [1]
    else:
        flt_txt = sp[0]
    text = ax.annotate(flt_txt,
                        (mean_mems[i], mean_bs[i]),
                        fontsize=8, zorder=3)
    texts.append(text)

ax.set_xlim([-1, 1])
ax.set_ylim([-1, 1])
# Add y = 0, x=0 lines
pyplot.plot([0,0],[-1, 1], linestyle='dashed', linewidth=1, c='0.5')
pyplot.plot([-1,1],[0, 0], linestyle='dashed', linewidth=1, c='0.5')
#ax.set_yscale('log')
ax.set_ylabel('B')
ax.set_xlabel('M')

# Now we add on some of the points from seismicity catalogs from Chen et al 2020
# First regions
mem_seis = [0.07, 0.25, -0.11, 0.35, -0.02, 0.31, 0.0, 0.21, -0.23]
b_seis = [0.10, 0.23, 0.05, 0.08, 0.09, 0.12, 0.31, 0.06, 0.03]
labels = ['Global', 'Japan', 'Taiwan','California', 'New Zealand',
          'North China', 'East Africa', 'Tibet', 'Xinjiang']
pyplot.scatter(mem_seis, b_seis, marker = '^', s=25, zorder=2, c='k')
for i, txt in enumerate(labels):
    text = ax.annotate(txt, (mem_seis[i], b_seis[i]), fontsize=8, zorder=3, style='italic')
    texts.append(text)
# Now individual faults from Chen et al 2020
mem_seis = [-0.15, -0.06, 0.23, -0.34]
b_seis = [-0.05, 0.07, 0.01, 0.02]
labels = ['Great Sumatran', 'North Anatolian', 'Sagaing', 'Xianshuihe']
pyplot.scatter(mem_seis, b_seis, marker = 'v', s=25, zorder=2, c='k')
for i, txt in enumerate(labels):
    text = ax.annotate(txt, (mem_seis[i], b_seis[i]), fontsize=8, zorder=3, style='italic',
                       fontweight='bold')
    texts.append(text)

print('Adjusting label locations')
adjust_text(texts, arrowprops=dict(arrowstyle='->', color='k', linewidth=0.5))

# Now we add 95% limits from synthetically generated datasets
for p in [68, 95]:
    if p == 95:
        ls = 'solid'
    elif p == 68:
        ls = 'dashed'
    data = np.genfromtxt(('../plotting/Exponential_B_M_%iper_contour_nsim_100000_nevents_%i.txt' % (p, min_number_events)),
                         delimiter=',')
    lb = 'Exponential %i%%' % p
    pyplot.plot(data[:,0], data[:,1], c='orangered', linestyle = ls, linewidth=1, zorder=1, label=lb)
    data = np.genfromtxt(('../plotting/Gamma_B_M_%iper_contour_nsim_100000_nevents_%i.txt' % (p, min_number_events)),
                         delimiter=',')
    lb = 'Gamma %i%%' % p 
    pyplot.plot(data[:,0], data[:,1], c='orange', linestyle = ls, linewidth=1, zorder=1, label=lb)
    data = np.genfromtxt(('../plotting/Weibull_B_M_%iper_contour_nsim_100000_nevents_%i.txt' % (p, min_number_events)),
                         delimiter=',')
    lb = 'Weibull %i%%' % p 
    pyplot.plot(data[:,0], data[:,1], c='yellow', linestyle = ls, linewidth=1, zorder=1, label=lb)
pyplot.legend()
# Add a legend using some dummy data
line1 = ax.scatter([1], [100], marker = 's', c = 'r', s=25)
line2 = ax.scatter([1], [100], marker = 's', c = 'g', s=25)
line3 = ax.scatter([1], [100], marker = 's', c = 'b', s=25)
line4 = ax.scatter([1], [100], marker = '^', c = 'k', s=25)
line5 = ax.scatter([1], [100], marker = 'v', c = 'k', s=25)
line6, = ax.plot([1, 2], [100, 101], c='orangered', linewidth=1)
line7, = ax.plot([1, 2], [100, 101], c='orange', linewidth=1)
line8, = ax.plot([1, 2], [100, 101], c='yellow', linewidth=1)
pyplot.legend((line1, line2, line3, line4, line5, line6, line7, line8),
              ('Normal', 'Strike slip', 'Reverse', 'Instrumental - regional',
               'Instrumental - single fault', 'Exponential', 'Gamma', 'Weibull'))
figname = 'B_M_phase_comparison_%s.png' % fig_comment 
fig.set_size_inches(w=8,h=8.)
pyplot.savefig(figname)

# Dump all results to a csv file
results_filename = 'Results_summary_%s.csv' % fig_comment
#print(mean_bs, len(mean_bs))
#print(burstiness_bounds[:,0])
#print(burstiness_bounds)
#print(burstiness_bounds[0,:])
#print(references, len(references))
#print(names, len(names))

all_results = np.vstack([names, references, mean_ltr, ltr_bounds[0,:], ltr_bounds[1,:],
                         mean_bs, burstiness_bounds[0,:], burstiness_bounds[1,:],
                         mean_mems, memory_bounds[0,:], memory_bounds[1,:]]).T
header = 'Name, reference, mean long-term annual rate, long-term rate 2.5p, long-term rate 97.5p, mean burstiness,'\
    'burstiness 2.5p, burstiness 97.5p, mean memory coefficient, memory coefficient 2.5p,' \
    'memory coefficient 97.5p'
np.savetxt(results_filename, all_results, header = header, delimiter=',', fmt="%s")

#############################
# Calculate percentage of burstiness values below zero
print(np.array(burstinesses))
b_neg = np.where(np.array(burstinesses) < 0, 1., 0)
num_neg_idv_flt = sum(b_neg.T) # Number of negative values for each individual fault
print('num_neg_idv_flt', num_neg_idv_flt)
percent_neg_idv_flt = num_neg_idv_flt/n_samples*100
print('percent_neg_idv_flt', percent_neg_idv_flt)
maj_neg = np.where(percent_neg_idv_flt > 50, 1.0, 0.)
print('maj_neg', maj_neg)
num_neg_total = np.sum(maj_neg)
fraction_maj_neg = np.sum(maj_neg)/n_faults
percent_maj_neg = np.sum(maj_neg)/n_faults*100 # Get percentage of faults with majority of values of B < 0
print('num_neg_total', num_neg_total)
print('percent_maj_neg', percent_maj_neg)
# Binomial distribution - if it was a Poisson process, we'd expect roughly
# half of our faults to dominantly B <= 0, and half dominantly B >= 0
binom_mod = binom(n_faults, (1-0.744))
prob_poisson = binom_mod.cdf((n_faults - num_neg_total))
print('prob_poisson', prob_poisson)
total_samples = len(b_neg.flatten())
num_neg = sum(b_neg.flatten())
percent_neg = num_neg/total_samples*100
print('n_faults', n_faults)
print('Total samples', total_samples)
print('Num neg', num_neg)
print('Percent neg', percent_neg)

# Plot histogram of all burstiness values against all random exponentially
# sampled burstiness values
pyplot.clf()
#pyplot.hist(np.array(burstinesses).flatten(), bins = 40,
#            alpha=0.5, density=False, label = 'Data')
#pyplot.savefig('burstiness_hist.png')
pyplot.hist(np.array(burstiness_expon.flatten()), bins = 60,
            alpha=0.5, density=True, label = 'Random sample')
pyplot.hist(np.array(burstinesses).flatten(), bins = 60,
            alpha=0.5, density=True, label = 'Data')
ax = pyplot.gca()
ax.set_xlabel('B')
ax.set_ylabel('Density')

pyplot.legend()

# Perform Kolmogorov-Smirnov test to see if our real
# fault data is less bursty than our exponentially distributed
# random data
ks_stat = ks_2samp(np.array(burstinesses).flatten(), np.array(burstiness_expon).flatten())
print('Komogorov-Smirnov statistic, p-value', ks_stat)
# This doesn't work on my current version of scipy
#ks_stat_less = kstest(np.array(burstinesses).flatten(), np.array(burstiness_expon).flatten(), alternative='less') 
#print('Komogorov-Smirnov statistic for real distirbution less than exponential, p-value', ks_stat_less)
lab = 'KS = %.2f\np value = %.2E' % (ks_stat[0], ks_stat[1])
ax.annotate(lab, (-0.8, 0.8), fontsize=10)
figname = 'burstiness_hist_random_%i.png' % min_number_events
pyplot.savefig(figname) 
