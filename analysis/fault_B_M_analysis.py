"""Script for sampling COV, burstiness and memory coeficient, and 
their  uncertainties, on many faults and plotting them

Jonathan Griffin
University of Otago
2020
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
from scipy.stats import expon, gamma, weibull_min, ks_2samp, kstest
# !!! Dangerous hack to swap Weibull for gamma
#from scipy.stats import weibull_min as gamma #
# !!!
from matplotlib import pyplot
from matplotlib.patches import PathPatch
import matplotlib.gridspec as gridspec
from matplotlib.ticker import FormatStrFormatter
from scipy.stats import binom, kde
from adjustText import adjust_text
from QuakeRates.dataman.event_dates import EventSet
from QuakeRates.dataman.parse_oxcal import parse_oxcal
from QuakeRates.dataman.parse_age_sigma import parse_age_sigma
from QuakeRates.dataman.parse_params import parse_param_file, \
    get_event_sets, file_len
from QuakeRates.utilities.bilinear import bilinear_reg_zero_slope, \
    bilinear_reg_fix, bilinear_reg_fix_zero_slope
from QuakeRates.utilities.memory_coefficient import burstiness, memory_coefficient

filepath = '../params'
param_file_list = glob(os.path.join(filepath, '*.txt'))

param_file_list_NZ = ['Akatore_TaylorSilva_2019.txt',
                      'AlpineHokuriCk_Berryman_2012_simple.txt',
                      'AlpineSouthWestland_Cochran_2017_simple.txt',
                      'AwatereEast_Nicol_2016_simple.txt',
                      'ClarenceEast_Nicol_2016_simple.txt',
                      'CloudyFault_Nicol_2016_simple.txt',
                      'Dunstan_GNS_unpub_simple.txt',
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
                      'Waitangi_GNS_unpub_simple.txt',
                      'Whakatane_Nicol_2016_simple.txt',
                      'Whirinaki_Nicol_2016_simple.txt']
# List of faults in study by Williams et al 2019
# Note this is not entirely the same, as there are some records from
# that study that are not included in ours.
param_file_list_W = ['AlpineHokuriCk_Berryman_2012_simple.txt',
                     'HaywardTysons_Lienkaemper_2007_simple.txt',
                     'SanJacintoMysticLake_Onderdonk_2018_simple.txt',
                     'NorthAnatolianElmacik_Fraser_2010_simple.txt',
                     'SanAndreasWrightwood_Weldon_2004_simple.txt',
                     'SanAndreasCarizzo_Akciz_2010_simple.txt',
                     'SanJacintoHogLake_Rockwell_2015_simple.txt',
                     'SanAndreasMissionCk_Fumal_2002_simple.txt',
                     'SanAndreasPalletCk_Scharer_2011_simple.txt',
                     'Xorkoli_Altyn_Tagh_Yuan_2018.txt',
                     'NorthAnatolianYaylabeli_Kozaci_2011_simple.txt',
                     'ElsinoreTemecula_Vaughan_1999_simple.txt',
                     'DeadSeaJordan_Ferry_2011_simple.txt',
                     'SanAndreasBigBend_Scharer_2017_simple.txt',
                     'WasatchBrigham_McCalpin_1996_simple.txt',
                     'Irpinia_Pantosti_1993_simple.txt',
                     'WasatchWeber_Duross_2011_simple.txt',
                     'WasatchNilphi_Duross_2017_simple.txt',
                     'LomaBlanca_Williams_2017_simple.txt',
                     'AlaskaPWSCopper_Plafker_1994_simple.txt',
                     'NankaiTrough_Hori_2004_simple.txt',
                     'CascadiaNth_Adams_1994_simple.txt',
                     'CascadiaSth_Goldfinger_2003_simple.txt',
                     'JavonCanyon_SarnaWojicki_1987_simple.txt',
                     'NewGuinea_Ota_1996_simple.txt',
                     'ChileMargin_Moernaut_2018_simple.txt']


#param_file_list = []
#for f in param_file_list_NZ:
#for f in param_file_list_W:
#    param_file_list.append(os.path.join(filepath, f))
n_samples = 10000  # Number of Monte Carlo samples of the eq chronologies
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
min_number_events = 5 # Use for all other calculations.
min_num_events_mem = 6 # Use for memory coefficient

#Summarise for comment to add to figure filename
fig_comment = ''
#fig_comment = 'NZ_examples_'
#fig_comment = 'Williams2019_'

for f in faulting_styles:
    fig_comment += f
    fig_comment += '_'
for t in tectonic_regions:
    fig_comment += t
    fig_comment += '_'
fig_comment += str(min_number_events)
#fig_comment += 'test_add_event_data'

def piecewise_linear(x, x0, y0, k1, k2):
    return np.piecewise(x, [x < x0], [lambda x:k1*x + y0-k1*x0, lambda x:k2*x + y0-k2*x0])
def camel_case_split(identifier):
    matches = finditer('.+?(?:(?<=[a-z])(?=[A-Z])|(?<=[A-Z])(?=[A-Z][a-z])|$)', identifier)
    return [m.group(0) for m in matches]

plot_colours = []
all_ie_times = []
added_events = [] # Store names of records where we've added an event due to
#                   exceptionally long current open interval
covs = []
cov_bounds = []
burstinesses = []
burstiness_bounds = []
burstiness_stds = []
burstinesses_expon = []
burstinesses_gamma = []
ie_gamma_alpha = []
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

names, event_sets, event_certainties, num_events, tect_regions, fault_styles = \
    get_event_sets(param_file_list, tectonic_regions,
                   faulting_styles, min_number_events)

references = []
# Get citations for each dataset from filename
for s in param_file_list:
    sp = s.split('_')
    if sp[0].split('/')[2] in names:
        references.append(sp[1] + ' ' + sp[2])
n_faults = len(names)
print('Number of faults', n_faults)

for i, event_set in enumerate(event_sets):
    # Handle cases with uncertain number of events. Where events identification is
    # unsure, event_certainty is given a value of 0, compared with 1 for certain
    # events
    # First generate chronologies assuming all events are certain
#    event_set.name = names[i]
    event_set.gen_chronologies(n_samples, observation_end=2020, min_separation=1)
    event_set.calculate_cov()
    event_set.cov_density()
    event_set.memory_coefficient()
    event_set.memory_spearman_rank_correlation()
    # Store all inter-event times for global statistics
    all_ie_times.append(event_set.interevent_times)
    # Now calculate some statistics on the sampled chronologies
    event_set.basic_chronology_stats()
    # Plot histogram of interevent times
    figfile = os.path.join(plot_folder, ('interevent_times_%s.png' % names[i]))
    event_set.plot_interevent_time_hist(fig_filename=figfile)
    # Fit gamma distirbution to event set data
    event_set.fit_gamma()
    ie_gamma_alpha.append(event_set.mean_gamma_alpha_all) # Get mean estimate of alpha
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
    # Generate random exponentially and gamma distributed samples of length num_events - 1
    # i.e. the number of inter-event times in the chronology. These will be used
    # later for testing
    scale = 100 # Fix scale, as burstiness is independent of scale for exponentiall distribution
    ie_times_expon = expon(scale=scale).rvs(size=(n_samples*(event_set.num_events-1)))
    ie_times_expon = np.reshape(np.array(ie_times_expon), (n_samples, (event_set.num_events-1)))
    ie_times_expon_T = ie_times_expon.T
    burst_expon = burstiness(ie_times_expon_T)
    # Gamma
    alpha_g = 2.3 #2.2 #1.6 ##2.35 #2.4 #2.0
    ie_times_g = gamma(alpha_g, scale=scale).rvs(size=(n_samples*(event_set.num_events-1)))
    ie_times_g = np.reshape(np.array(ie_times_g), (n_samples, (event_set.num_events-1)))
    ie_times_g_T = ie_times_g.T
    burst_g = burstiness(ie_times_g_T)
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
        ie_times_g_certain = gamma(alpha_g, scale=scale).rvs(size=(n_samples*(event_set.num_events-1)))
        ie_times_g_certain = np.reshape(np.array(ie_times_g_certain), (n_samples, (event_set.num_events-1)))
        ie_times_g_certain_T = ie_times_g_certain.T
        burst_g_certain = burstiness(ie_times_g_T)
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
        combined_burst_g = np.concatenate([burst_g[:half_n],
                                               burst_g_certain[:half_n]])
        covs.append(combined_covs)
        burstinesses.append(combined_burstiness)
        memory_coefficients.append(combined_memory)
        memory_stds.append(np.std(np.array(combined_memory)))
        memory_spearman_coefficients.append(combined_memory_spearman)
        memory_spearman_lag2_coef.append(combined_memory_spearman_lag2)
        burstinesses_expon.append(combined_burst_expon)
        burstinesses_gamma.append(combined_burst_g)
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
        burstinesses_gamma.append(burst_g)
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
    if event_set.add_events: # List of records where we model long open interval
        added_events.append(event_set.name)
        
# Convert to numpy arrays and transpose where necessary
num_events = np.array(num_events)
all_ie_times = np.array(all_ie_times)
max_interevent_times = np.array(max_interevent_times)
min_interevent_times = np.array(min_interevent_times)
min_paired_interevent_times = np.array(min_paired_interevent_times)
std_max_interevent_times = np.array(std_max_interevent_times)
std_min_interevent_times = np.array(std_min_interevent_times)
std_min_paired_interevent_times = np.array(std_min_paired_interevent_times)
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
burstiness_gamma = np.array(burstinesses_gamma)
inds = np.where(num_events >= min_num_events_mem) # Get memory coefficients for more than 6 events
memory_coefficients = np.array(memory_coefficients)
memory_coefficients_min = memory_coefficients[inds]
memory_stds = np.array(memory_stds)
memory_stds_min = memory_stds[inds]
memory_bounds_min = np.array(memory_bounds)[inds].T
memory_bounds = np.array(memory_bounds).T
memory_spearman_bounds = np.array(memory_spearman_bounds).T
memory_spearman_lag2_bounds = np.array(memory_spearman_lag2_bounds).T
ie_gamma_alpha = np.array(ie_gamma_alpha)

# Now plot the means and 95% error bars of COV
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

################################

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
# Add B=0 linear
pyplot.plot([1./1000000, 1./40], [0, 0], linestyle='dashed', linewidth=1, c='0.5')
ax.set_xscale('log')
ax.set_xlabel('Long-term rate (events per year)')
ax.set_ylabel('B')

# Now do a bi-linear fit to the data
mean_bs = np.array(mean_bs)
indices = np.flatnonzero(mean_ltr > 3e-4)
indices = indices.flatten()

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
pyplot.plot(xvals_short, yvals, c='0.2')

# Fit slow faults
if len(indices_slow_faults > 1):
    lf_slow = np.polyfit(np.log10(mean_ltr[indices_slow_faults]),
                         mean_bs[indices_slow_faults], 1)
    xvals_short = np.arange(1e-6, 1.5e-4, 1e-6)
    yvals = lf_slow[0]*np.log10(xvals_short) + lf_slow[1]
    pyplot.plot(xvals_short, yvals, c='0.2')
# Add formula for linear fits of data
print('Fits for B vs LTR')
txt = 'Y = {:=+6.2f} +/- {:4.2f}'.format(lf[1], std_lf)
print(txt)
ax.annotate(txt, (2e-4, 0.2), fontsize=8)
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
print(out.sum_square)
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
bilin_hxfix_cons_slope = odrpack.Model(bilinear_reg_fix_zero_slope)
odr = odrpack.ODR(data, bilin_hxfix_cons_slope, beta0=[-3, -1.0])
odr.set_job(fit_type=0)
out = odr.run()
print('bilinear hxfix_cons_slope')
print(out.sum_square)
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

figname = 'burstiness_vs_lt_rate_%s.png' % fig_comment 
pyplot.savefig(figname)

#########################

# Plot burstiness against slip rate
pyplot.clf()
ax = pyplot.subplot(111)
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
ax.set_ylim([-1, 1])
ax.set_xlim([1./1000, 100])
# Add B=0 linear
pyplot.plot([1./1000, 100], [0, 0], linestyle='dashed', linewidth=1, c='0.5')
ax.set_xscale('log')
ax.set_xlabel('Slip rate (mm/yr)')
ax.set_ylabel('B')

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
linear  = odrpack.Model(f)
odr = odrpack.ODR(data, linear, beta0=[-1, -1.0,]) 
odr.set_job(fit_type=0)
out = odr.run()
out.pprint()
a = out.beta[0]
b = out.beta[1]
xvals = np.arange(1.e-4, 1e2, 1e-2)
yrng = a*np.log10(xvals) + b #10**(b + a * xvals)
pyplot.plot(xvals, yrng, c='0.6')
txt = 'Y = {:4.2f}Log(x) {:=+6.2f}'.format(a, b)
print(txt)                                                                                                                                 
ax.annotate(txt, (1e0, 0.9), color='0.6') 

# Now try bilinear fixed hinge
bilin = odrpack.Model(bilinear_reg_fix_zero_slope)
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
pyplot.plot(xvals, yrng, c='0.2')
txt = 'Y = {:4.2f}Log(x) {:=+6.2f}, x < {:4.2f}'.format(a, b, np.power(10,hxfix))
print(txt)
ax.annotate(txt, (2e-3, 0.9), color='0.2')
txt = 'Y = {:4.2f}, x >= {:4.2f}'.format(ylevel, np.power(10,hxfix)) 
print(txt)
ax.annotate(txt, (1.2e-2, 0.8), color='0.2')

figname = 'burstiness_vs_slip_rate_%s.png' % fig_comment   
pyplot.savefig(figname)
figname = 'burstiness_vs_slip_rate_%s.pdf' % fig_comment
pyplot.savefig(figname)

# Plot memory coefficients against long term rates
pyplot.clf()
ax = pyplot.subplot(111)
mean_mems = []
mean_ltr_mem = mean_ltr[inds]
ltr_bounds_mem = ltr_bounds.T[inds].T
for i, mem_set in enumerate(memory_coefficients):
    mean_mem = np.mean(mem_set)
#    print('Mean memory coefficient combined', mean_mem)
    mean_mems.append(mean_mem)
mean_mems = np.array(mean_mems)
colours = []
plot_colours_mem = list(np.array(plot_colours)[inds])
for mean_mem in mean_mems:
    if mean_mem <= -0.05:
        colours.append('b')
    elif mean_mem > -0.05 and mean_mem <= 0.05:
        colours.append('g')
    else:
        colours.append('r')
pyplot.errorbar(mean_ltr_mem, mean_mems[inds],
                xerr = ltr_bounds_mem,
                ecolor = '0.3',
                elinewidth=0.7,
                linestyle="None",
                zorder=1)
pyplot.errorbar(mean_ltr_mem, mean_mems[inds],
                yerr = memory_bounds_min,
                ecolor = '0.3',
                elinewidth=0.7,
                linestyle="None",
                zorder=1)
pyplot.scatter(mean_ltr_mem, mean_mems[inds], marker = 's', c=plot_colours_mem,
               s=25, zorder=2)
for i, txt in enumerate(names):
    if max_interevent_times[i] > 10 and annotate_plots:
        ax.annotate(txt[:4],
                    (mean_ltr[i], mean_mems[i]),
                    fontsize=8)
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
ax.set_xlabel('M (Spearman Rank Lag-1)')
ax.set_ylabel('M (Spearman Rank Lag-2)')
figname = 'memory_coefficient_Spearman_L1_vs_L2_%s.png' % fig_comment
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

# Bilinear fixed hinge and constant slope ODR
hxfix = np.log10(2e-4)
bilin_hxfix_cons_slope = odrpack.Model(bilinear_reg_fix_zero_slope)
data = odrpack.RealData(np.log10(mean_ltr), mean_bs,
                        sx=np.log10(np.sqrt(long_term_rate_stds)), sy=np.sqrt(burstiness_stds))
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
pyplot.errorbar(mean_ltr_mem, mean_mems[inds],
                xerr = ltr_bounds_mem,
                ecolor = '0.3',
                elinewidth=0.5,
                linestyle="None",
                zorder=1)
pyplot.errorbar(mean_ltr_mem, mean_mems[inds],
                yerr = memory_bounds_min,
                ecolor = '0.3',
                elinewidth=0.5,
                linestyle="None",
                zorder=1)
pyplot.scatter(mean_ltr_mem, mean_mems[inds], marker = 's', c=plot_colours_mem,
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

def linear_func(B, x):
    return B[0]*x + B[1]   
# Bilinear fixed hinge and constant slope ODR
hxfix = np.log10(2e-4)
bilin_hxfix_cons_slope = odrpack.Model(bilinear_reg_fix_zero_slope)
long_term_rate_stds_mem = long_term_rate_stds[inds]
data = odrpack.RealData(np.log10(mean_ltr_mem), mean_mems[inds],
                        sx=np.log10(np.sqrt(long_term_rate_stds_mem)), sy=np.sqrt(memory_stds_min))
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
pyplot.plot(xvals, yrng, c='0.4', linestyle = '--')
txt = 'y = {:4.2f}Log(x) {:=+6.2f}, x < {:3.1E}'.format(a, b, np.power(10, hxfix))
ax.annotate(txt, (1.5e-6, -0.85), fontsize=8) 
txt = 'y = {:4.2f}, x >= {:3.1E}'.format(ylevel, np.power(10, hxfix))
ax.annotate(txt, (1.5e-6, -0.95), fontsize=8)

# Linear ODR fit
linear  = odrpack.Model(linear_func)
odr = odrpack.ODR(data, linear, beta0=[-1, -1.0,]) 
odr.set_job(fit_type=0)
out = odr.run()
out.pprint()
a = out.beta[0]
b = out.beta[1]
xvals = np.arange(1.e-4, 1e2, 1e-2)
yrng = a*np.log10(xvals) + b #10**(b + a * xvals)
ax.annotate('b)', (-0.23, 0.98), xycoords = 'axes fraction', fontsize=10)

# Add third plot
pyplot.subplot2grid((3, 2), (1,0), colspan=1, rowspan=1)
ax = pyplot.gca()
mean_bs_mem = mean_bs[inds]
burstiness_bounds_mem = burstiness_bounds.T[inds].T
pyplot.errorbar(mean_mems[inds], mean_bs_mem,
                xerr = memory_bounds_min,
                ecolor = '0.3',
                elinewidth=0.5,
                linestyle="None",
                zorder=1)
pyplot.errorbar(mean_mems[inds], mean_bs_mem,
                yerr = burstiness_bounds_mem,
                ecolor = '0.3',
                elinewidth=0.5,
                linestyle="None",
                zorder=1)
pyplot.scatter(mean_mems[inds], mean_bs_mem, marker = 's', c=plot_colours_mem,
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

#Orthogonal linear fit
def linear_func(B, x):
    return B[0]*x + B[1]
linear_model = Model(linear_func)
burstiness_stds_mem = burstiness_stds[inds]
data = RealData(np.array(mean_mems[inds]).flatten(),
                np.array(mean_bs_mem).flatten(),
                sx = np.sqrt(memory_stds.flatten()),
                sy = np.sqrt(burstiness_stds_mem.flatten()))
# Set up ODR with the model and data
odr = ODR(data, linear_model, beta0=[1., -1.])
out = odr.run()
out.pprint()
xvals = np.arange(-0.75, 0.75, 0.01)
yvals = linear_func(out.beta, xvals)
pyplot.plot(xvals, yvals, c='0.4')
ax.set_ylabel('B', fontsize=10)
ax.set_xlabel('M', fontsize=10)
# Add formula for linear fit to low-end of data
txt = 'y = {:4.2f}Log(x) {:=+6.2f}'.format(out.beta[0], out.beta[1])
print(txt)
ax.annotate(txt, (-0.95, 0.8), fontsize=8)
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

# Bilinear fixed hinge and constant slope ODR
hxfix = np.log10(2e-4)
bilin_hxfix_cons_slope = odrpack.Model(bilinear_reg_fix_zero_slope)
data = odrpack.RealData(np.log10(mean_ltr), np.log10(ratio_min_pair_max),
                        sx=np.log10(np.sqrt(long_term_rate_stds)),
                        sy=np.log10(np.sqrt(std_ratio_min_pair_max)))
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
figname = 'combined_plot_%s.pdf' % fig_comment 
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
pyplot.errorbar(mean_mems[inds], mean_bs_mem,
                xerr = memory_bounds_min,
                ecolor = '0.8',
                elinewidth=0.7,
                linestyle="None",
                zorder=1)
pyplot.errorbar(mean_mems[inds], mean_bs_mem,
                yerr = burstiness_bounds_mem,
                ecolor = '0.8',
                elinewidth=0.7,
                linestyle="None",
                zorder=1)
pyplot.scatter(mean_mems[inds], mean_bs_mem, marker = 's', c=plot_colours_mem,
               s=25, zorder=2)
texts = []

for i, txt in enumerate(list(np.array(names)[inds])):
    sp = camel_case_split(txt)
    # Handle special cases of two word fault names
    if sp[0] == 'San' or sp[0] == 'Dead':
        flt_txt = sp[0] + ' '  + sp [1] #+ ' (' + sp [2] + ')' # Uncomment to get segment names
    elif sp[0] == 'Loma':
        flt_txt = sp[0] + ' '  + sp [1]
    else:
        flt_txt = sp[0]
    text = ax.annotate(flt_txt,
                        (mean_mems[inds][i], mean_bs_mem[i]),
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
mem_seis = [0.07]#, 0.25, -0.11, 0.35, -0.02, 0.31, 0.0, 0.21, -0.23]
b_seis = [0.10]#, 0.23, 0.05, 0.08, 0.09, 0.12, 0.31, 0.06, 0.03]
labels = ['Global']#, 'Japan', 'Taiwan','California', 'New Zealand',
#          'North China', 'East Africa', 'Tibet', 'Xinjiang']
pyplot.scatter(mem_seis, b_seis, marker = '^', s=25, zorder=2, c='k')
for i, txt in enumerate(labels):
    text = ax.annotate(txt, (mem_seis[i], b_seis[i]), fontsize=8, zorder=3, style='italic')
    texts.append(text)
# Now individual faults from Chen et al 2020
#mem_seis = [-0.15, -0.06, 0.23, -0.34]
#b_seis = [-0.05, 0.07, 0.01, 0.02]
#labels = ['Great Sumatran', 'North Anatolian', 'Sagaing', 'Xianshuihe']
#pyplot.scatter(mem_seis, b_seis, marker = 'v', s=25, zorder=2, c='k')
#for i, txt in enumerate(labels):
#    text = ax.annotate(txt, (mem_seis[i], b_seis[i]), fontsize=8, zorder=3, style='italic',
#                       fontweight='bold')
#    texts.append(text)

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
    pyplot.plot(data[:,0], data[:,1], c='slategrey', linestyle = ls, linewidth=1, zorder=1, label=lb)
pyplot.legend()
# Add a legend using some dummy data
line1 = ax.scatter([1], [100], marker = 's', c = 'r', s=25)
line2 = ax.scatter([1], [100], marker = 's', c = 'g', s=25)
line3 = ax.scatter([1], [100], marker = 's', c = 'b', s=25)
line4 = ax.scatter([1], [100], marker = '^', c = 'k', s=25)
#line5 = ax.scatter([1], [100], marker = 'v', c = 'k', s=25)
line6, = ax.plot([1, 2], [100, 101], c='orangered', linewidth=1)
line7, = ax.plot([1, 2], [100, 101], c='orange', linewidth=1)
line8, = ax.plot([1, 2], [100, 101], c='slategrey', linewidth=1)
pyplot.legend((line1, line2, line3, line4, line6, line7, line8),
              ('Normal', 'Strike slip', 'Reverse', 'Instrumental - global',
               'Exponential', 'Gamma', 'Weibull'))
figname = 'B_M_phase_comparison_%s.png' % fig_comment 
fig.set_size_inches(w=8,h=8.)
pyplot.savefig(figname)
figname = 'B_M_phase_comparison_%s.pdf' % fig_comment 
pyplot.savefig(figname)

###############################################################

# Dump all results to a csv file
results_filename = 'Results_summary_%s.csv' % fig_comment
all_results = np.vstack([names, references, num_events, mean_ltr, (mean_ltr-ltr_bounds[0,:]), (mean_ltr+ltr_bounds[1,:]),
                         mean_bs, (mean_bs-burstiness_bounds[0,:]), (mean_bs+burstiness_bounds[1,:]),
                         mean_mems, (mean_mems-memory_bounds[0,:]), (mean_mems+memory_bounds[1,:])]).T
header = 'Name, reference, number of events, mean long-term annual rate, long-term rate 2.5p, long-term rate 97.5p, mean burstiness,'\
    'burstiness 2.5p, burstiness 97.5p, mean memory coefficient, memory coefficient 2.5p,' \
    'memory coefficient 97.5p'
np.savetxt(results_filename, all_results, header = header, delimiter=',', fmt="%s")

################################################################

# Plot histogram of all burstiness values against all random exponentially
# sampled burstiness values
pyplot.clf()
burstiness_expon = np.array(burstiness_expon)
burstinesses = np.array(burstinesses)
pyplot.hist(np.array(burstiness_expon.flatten()), bins = 60,
            alpha=0.5, density=True, label = 'Random sample')
pyplot.hist(burstinesses.flatten(), bins = 60,
            alpha=0.5, density=True, label = 'Data')
ax = pyplot.gca()
ax.set_xlabel('B')
ax.set_ylabel('Density')

pyplot.legend()

# Perform Kolmogorov-Smirnov test to see if our real
# fault data is less bursty than our exponentially distributed
# random data
# All data first
ks_stat = ks_2samp(np.array(burstinesses).flatten(), np.array(burstiness_expon).flatten())
print('Komogorov-Smirnov statistic, p-value', ks_stat)

# Get proportion of overlap
b_p = burstinesses.flatten() - burstiness_expon.flatten()
b_neg = np.count_nonzero(b_p < 0)   
sum_neg = np.sum(b_neg)
print('Sum_neg Expon', sum_neg)
print('Total size', len(burstinesses.flatten()))
print('As percent', sum_neg/len(burstinesses.flatten()))
lab = 'KS = %.2f\np value = %.2E' % (ks_stat[0], ks_stat[1])
ax.annotate(lab, (-0.8, 0.8), fontsize=10)
figname = 'burstiness_hist_random_%s.png' % fig_comment
pyplot.savefig(figname)
# Dump out as text file
f_name = 'burstiness_hist_random_%s.txt' % fig_comment
data = np.array([burstinesses.flatten(), burstiness_expon.flatten()]).T 
np.savetxt(f_name, data, header='Data,Exponential', delimiter=',')

######################
#Do KS test sample by sample
ks_stats = []
p_values = []
for b in burstinesses.T:
    ks_stat = ks_2samp(b, burstiness_expon.flatten()) 
    ks_stats.append(ks_stat[0])
    p_values.append(ks_stat[1])
pyplot.clf()
pyplot.hist(p_values, bins=40, density=True)
ax = pyplot.gca()
ax.set_xlabel('p value')
ax.set_ylabel('Density')
figname = 'burstiness_hist_KS_pvalue_random_%s.png' % fig_comment  
pyplot.savefig(figname)
########

###########

# Now do for only high activity rate faults
indices = np.argwhere(mean_ltr > 2e-4)
burstiness_fast = np.array(burstinesses)[indices]
burstiness_expon_fast = np.array(burstiness_expon)[indices]

# Plot histogram of all burstiness values against all random exponentially
# sampled burstiness values
pyplot.clf()
pyplot.hist(np.array(burstiness_expon_fast.flatten()), bins = 40,
            alpha=0.5, density=True, label = 'Random sample')
pyplot.hist(np.array(burstiness_fast).flatten(), bins = 40,
            alpha=0.5, density=True, label = 'Data')
ax = pyplot.gca()
ax.set_xlabel('B')
ax.set_ylabel('Density')
ax.set_ylim([0.0, 4])
ax.set_xlim([-1., 0.5]) 
pyplot.legend()

figname = 'burstiness_hist_random_high_activity_%s.png' % fig_comment
ks_stat = ks_2samp(burstiness_fast.flatten(), burstiness_expon_fast.flatten())
print('Komogorov-Smirnov statistic for high activity rate faults, p-value', ks_stat)    
lab = 'KS = %.2f\np value = %.2E' % (ks_stat[0], ks_stat[1])
ax.annotate(lab, (-0.8, 0.8), fontsize=10)
pyplot.savefig(figname)

####################

# Now do only for low activity rate faults
indices_slow_faults = np.flatnonzero(mean_ltr <= 2e-4)
indices_slow_faults = indices_slow_faults.flatten()
burstiness_slow = burstinesses[indices_slow_faults]
burstiness_expon_slow = np.array(burstiness_expon)[indices_slow_faults]

# Plot histogram of all burstiness values against all random exponentially
# sampled burstiness values
pyplot.clf()
pyplot.hist(burstiness_expon_slow.flatten(), bins = 40,
            alpha=0.5, density=True, label = 'Random sample')
pyplot.hist(burstiness_slow.flatten(), bins = 40,
            alpha=0.5, density=True, label = 'Data')
ax = pyplot.gca()
ax.set_xlabel('B')
ax.set_ylabel('Density')
ax.set_ylim([0.0, 4])
ax.set_xlim([-1., 0.5])
pyplot.legend() 
figname = 'burstiness_hist_random_low_activity_%s.png' % fig_comment
# Calculate Kolmogorov-Smirnov statistic
ks_stat = ks_2samp(burstiness_slow.flatten(), burstiness_expon_slow.flatten())
print('Komogorov-Smirnov statistic for low activity rate faults, p-value', ks_stat)
lab = 'KS = %.2f\np value = %.2E' % (ks_stat[0], ks_stat[1])
ax.annotate(lab, (-0.8, 0.8), fontsize=10)
pyplot.savefig(figname) 

###########################################
# Plot histogram of all burstiness values against all random gamma distributions
# sampled burstiness values
burstiness_gamma = np.array(burstiness_gamma)
pyplot.clf()
pyplot.hist(burstiness_gamma.flatten(), bins = 60,
            alpha=0.5, density=True, label = 'Random sample')
pyplot.hist(burstinesses.flatten(), bins = 60,
            alpha=0.5, density=True, label = 'Data')
ax = pyplot.gca()
ax.set_xlabel('B')
ax.set_ylabel('Density')

pyplot.legend()

# Perform Kolmogorov-Smirnov test to see if our real
# fault data is less bursty than our gamma distributed
# random data
# All data first
ks_stat = ks_2samp(burstinesses.flatten(), burstiness_gamma.flatten())
print('Komogorov-Smirnov statistic for gamma distribution, p-value', ks_stat)
lab = 'KS = %.2f\np value = %.2E' % (ks_stat[0], ks_stat[1])
ax.annotate(lab, (-0.8, 0.8), fontsize=10)

# Get proportion of overlap
b_p = burstinesses.flatten() - burstiness_gamma.flatten()
b_neg = np.count_nonzero(b_p < 0)   
sum_neg = np.sum(b_neg)
print('Sum_neg Gamma', sum_neg)
print('Total size', len(burstinesses.flatten()))
print('As percent', sum_neg/len(burstinesses.flatten()))
      
figname = 'burstiness_hist_gamma_%s.png' % fig_comment
pyplot.savefig(figname) 
# Dump out as text file
f_name = 'burstiness_hist_gamma_%s.txt' % fig_comment
data = np.array([burstinesses.flatten(), burstiness_gamma.flatten()]).T
np.savetxt(f_name, data, header='Data,Gamma', delimiter=',')

######################
#Do KS test sample by sample
burstiness_gamma_fast = np.array(burstiness_gamma)[indices]

ks_stats = []
p_values = []
for b in burstinesses.T:
    ks_stat = ks_2samp(b.flatten(), burstiness_gamma.flatten()) 
    ks_stats.append(ks_stat[0])
    p_values.append(ks_stat[1])
pyplot.clf()
#p_values_fast = np.array(p_values)[indices]
pyplot.hist(p_values, bins=40, density=True)
ax = pyplot.gca()
ax.set_xlabel('p value')
ax.set_ylabel('Density')
figname = 'burstiness_hist_KS_pvalue_gamma_%s.png' % fig_comment  
pyplot.savefig(figname)
########

# Now do only for high activity rate faults
# Plot histogram of all burstiness values against all random exponentially
# sampled burstiness values
pyplot.clf()
pyplot.hist(burstiness_gamma_fast.flatten(), bins = 40,
            alpha=0.5, density=True, label = 'Random sample')
pyplot.hist(burstiness_fast.flatten(), bins = 40,
            alpha=0.5, density=True, label = 'Data')
ax = pyplot.gca()
ax.set_xlabel('B')
ax.set_ylabel('Density')
ax.set_ylim([0.0, 4])
ax.set_xlim([-1., 0.5]) 
pyplot.legend()


figname = 'burstiness_hist_gamma_high_activity_%s.png' % fig_comment
ks_stat = ks_2samp(burstiness_fast.flatten(), burstiness_gamma_fast.flatten())
print('Komogorov-Smirnov statistic for high activity rate faults, p-value', ks_stat)    
lab = 'KS = %.2f\np value = %.2E' % (ks_stat[0], ks_stat[1])
ax.annotate(lab, (-0.8, 0.8), fontsize=10)
pyplot.savefig(figname)

####################

# Now do only for low activity rate faults

burstiness_gamma_slow = np.array(burstiness_gamma)[indices_slow_faults]

# Plot histogram of all burstiness values against all random gamma
# sampled burstiness values
pyplot.clf()
pyplot.hist(burstiness_gamma_slow.flatten(), bins = 40,
            alpha=0.5, density=True, label = 'Random sample')
pyplot.hist(burstiness_slow.flatten(), bins = 40,
            alpha=0.5, density=True, label = 'Data')
ax = pyplot.gca()
ax.set_xlabel('B')
ax.set_ylabel('Density')
ax.set_ylim([0.0, 4])
ax.set_xlim([-1., 0.5])
pyplot.legend() 
figname = 'burstiness_hist_gamma_low_activity_%s.png' % fig_comment
# Calculate Kolmogorov-Smirnov statistic
ks_stat = ks_2samp(burstiness_slow.flatten(), burstiness_gamma_slow.flatten())
print('Komogorov-Smirnov statistic for low activity rate faults, p-value', ks_stat)
lab = 'KS = %.2f\np value = %.2E' % (ks_stat[0], ks_stat[1])
ax.annotate(lab, (-0.8, 0.8), fontsize=10)
pyplot.savefig(figname) 


#######################################3
# In this analysis, we now calculate the KS statistic
# for each fault individually, and plot these.
all_pvalues = []
print(np.shape(burstinesses))
print(np.shape(burstiness_expon))
for i, b in enumerate(burstinesses):
    ks = ks_2samp(b, burstiness_expon[i])
    all_pvalues.append(ks[0])
pyplot.clf()
pyplot.hist(all_pvalues, bins=50, density=True)
ax = pyplot.gca()
ax.set_xlabel('p value')
ax.set_ylabel('Density')
figname = 'ks_p_value_hist_%s.png' % fig_comment
pyplot.savefig(figname)

#######################################################3
# Now make a nice combined figure showing all the results
# 4 rows by 3 columns
# Plot results against expected distirbutions for Poisson and Gamma distributions.
# Do this for: All data; high activity rate data; low activity rate data;
# Strike-slip faults; normal faults; reverse faults;
pyplot.clf()
fig = pyplot.figure(1)
# set up subplot grid
gridspec.GridSpec(4, 3)
#First plot - all data against Poisson
pyplot.subplot2grid((4, 3), (0,0), colspan=1, rowspan=1)
#Do KS test sample by sample
ks_stats = []
p_values = []
for b in burstinesses.T:
    ks_stat = ks_2samp(b.flatten(), burstiness_expon.flatten()) 
    ks_stats.append(ks_stat[0])
    p_values.append(ks_stat[1])
p_reject = (np.array(p_values) < 0.05).sum() / len(p_values)
if p_reject < 0.95:
    rej = 'Accept'
else:
    rej = 'Reject'
pyplot.hist(np.array(burstiness_expon.flatten()), bins = 60,
            alpha=0.5, density=True, label = 'Exponential', color='#1f77b4')
pyplot.hist(burstinesses.flatten(), bins = 60, color='#ff7f0e',
            alpha=0.5, density=True, label = 'Data')
ax = pyplot.gca()
ax.set_xlim([-1.0, 0.5])
ax.set_ylim([0.0, 3.2])
ax.set_xlabel('B')
ax.set_ylabel('Density')
pyplot.legend(loc=1, fontsize=9, handlelength=1.5, framealpha=0.2)
# Annotate figure
txt = 'p reject: %.2f\n%s\nAll' % (p_reject, rej)
ax.annotate(txt, (0.03, 0.77), xycoords = 'axes fraction', fontsize = 10)
ax.annotate('a)', (-0.23, 0.98), xycoords = 'axes fraction', fontsize=10) 
##############
# Second subplot - high activity rate faults
pyplot.subplot2grid((4, 3), (0,1), colspan=1, rowspan=1)
#Do KS test sample by sample
ks_stats = []
p_values = []
for b in burstiness_fast.T:
    ks_stat = ks_2samp(b.flatten(), burstiness_expon_fast.flatten()) 
    ks_stats.append(ks_stat[0])
    p_values.append(ks_stat[1])
p_reject = (np.array(p_values) < 0.05).sum() / len(p_values)
if p_reject < 0.95:
    rej = 'Accept'
else:
    rej = 'Reject'
pyplot.hist(np.array(burstiness_expon_fast.flatten()), bins = 60,
            alpha=0.5, density=True, label = 'Exponential', color='#1f77b4')
pyplot.hist(burstiness_fast.flatten(), bins = 60,
            alpha=0.5, density=True, label = 'Data', color='#ff7f0e')
ax = pyplot.gca()
ax.set_xlim([-1.0, 0.5])
ax.set_xlabel('B')
ax.set_ylabel('Density')
# Annotate figure
txt = 'p reject: %.2f\n%s\nHigh rate' % (p_reject, rej)
ax.annotate(txt, (0.03, 0.77), xycoords = 'axes fraction', fontsize = 10)
ax.annotate('b)', (-0.23, 0.98), xycoords = 'axes fraction', fontsize=10)

# Third subplot - low activity rate faults
pyplot.subplot2grid((4, 3), (0,2), colspan=1, rowspan=1)
#Do KS test sample by sample
ks_stats = []
p_values = []
for b in burstiness_slow.T:
    ks_stat = ks_2samp(b.flatten(), burstiness_expon_slow.flatten()) 
    ks_stats.append(ks_stat[0])
    p_values.append(ks_stat[1])
p_reject = (np.array(p_values) < 0.05).sum() / len(p_values)
if p_reject < 0.95:
    rej = 'Accept'
else:
    rej = 'Reject'
pyplot.hist(np.array(burstiness_expon_slow.flatten()), bins = 60,
            alpha=0.5, density=True, label = 'Exponential', color='#1f77b4')
pyplot.hist(burstiness_slow.flatten(), bins = 60,
            alpha=0.5, density=True, label = 'Data', color='#ff7f0e')
ax = pyplot.gca()
ax.set_xlim([-1.0, 0.5])
ax.set_xlabel('B')
ax.set_ylabel('Density')
ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
# Annotate figure
txt = 'p reject: %.2f\n%s\nLow Rate' % (p_reject, rej)
ax.annotate(txt, (0.03, 0.77), xycoords = 'axes fraction', fontsize = 10)
ax.annotate('c)', (-0.23, 0.98), xycoords = 'axes fraction', fontsize=10)

##############
# Add fourth subplot - strike-slip faults
pyplot.subplot2grid((4, 3), (1,0), colspan=1, rowspan=1)
fault_styles = np.array(fault_styles)
indices_ss = np.argwhere(fault_styles == 'Strike_slip')
indices_ss_hs = np.intersect1d(indices, indices_ss)   
burstiness_ss = burstinesses[indices_ss_hs]
burstiness_expon_ss = burstiness_expon[indices_ss_hs]
#Do KS test sample by sample
ks_stats = []
p_values = []
for b in burstiness_ss.T:
    ks_stat = ks_2samp(b.flatten(), burstiness_expon_ss.flatten()) 
    ks_stats.append(ks_stat[0])
    p_values.append(ks_stat[1])
p_reject = (np.array(p_values) < 0.05).sum() / len(p_values)
if p_reject < 0.95:
    rej = 'Accept'
else:
    rej = 'Reject'
pyplot.hist(np.array(burstiness_expon_ss.flatten()), bins = 60,
            alpha=0.5, density=True, label = 'Exponential', color='#1f77b4')
pyplot.hist(burstiness_ss.flatten(), bins = 60,
            alpha=0.5, density=True, label = 'Data', color='#ff7f0e')
ax = pyplot.gca()
ax.set_xlim([-1.0, 0.5])
ax.set_xlabel('B')
ax.set_ylabel('Density')
# Annotate figure
txt = 'p reject: %.2f\n%s\nStrike-slip\n(High rate)' % (p_reject, rej)
ax.annotate(txt, (0.03, 0.72), xycoords = 'axes fraction', fontsize = 10)
ax.annotate('d)', (-0.23, 0.98), xycoords = 'axes fraction', fontsize=10)

##############
# Add fifth subplot - normal faults
pyplot.subplot2grid((4, 3), (1,1), colspan=1, rowspan=1)
indices_n = np.argwhere(fault_styles == 'Normal')
# Get indices of normal faults with high slip rates
indices_n_hs = np.intersect1d(indices, indices_n)
burstiness_n = burstinesses[indices_n_hs]
burstiness_expon_n = burstiness_expon[indices_n_hs]
#Do KS test sample by sample
ks_stats = []
p_values = []
for b in burstiness_n.T:
    ks_stat = ks_2samp(b.flatten(), burstiness_expon_n.flatten()) 
    ks_stats.append(ks_stat[0])
    p_values.append(ks_stat[1])
p_reject = (np.array(p_values) < 0.05).sum() / len(p_values)
if p_reject < 0.95:
    rej = 'Accept'
else:
    rej = 'Reject'
pyplot.hist(np.array(burstiness_expon_n.flatten()), bins = 60,
            alpha=0.5, density=True, label = 'Exponential', color='#1f77b4')
pyplot.hist(burstiness_n.flatten(), bins = 60,
            alpha=0.5, density=True, label = 'Data', color='#ff7f0e')
ax = pyplot.gca()
ax.set_xlim([-1.0, 0.5])
ax.set_xlabel('B')
ax.set_ylabel('Density')
# Annotate figure
txt = 'p reject: %.2f\n%s\nNormal\n(High rate)' % (p_reject, rej)
ax.annotate(txt, (0.03, 0.72), xycoords = 'axes fraction', fontsize = 10)
ax.annotate('e)', (-0.23, 0.98), xycoords = 'axes fraction', fontsize=10)

##############
# Add sixth subplot - reverse faults
pyplot.subplot2grid((4, 3), (1,2), colspan=1, rowspan=1)
indices_r = np.argwhere(fault_styles == 'Reverse')
indices_r_hs = np.intersect1d(indices, indices_r)
burstiness_r = burstinesses[indices_r_hs]
burstiness_expon_r = burstiness_expon[indices_r_hs]
#Do KS test sample by sample
ks_stats = []
p_values = []
for b in burstiness_r.T:
    ks_stat = ks_2samp(b.flatten(), burstiness_expon_r.flatten()) 
    ks_stats.append(ks_stat[0])
    p_values.append(ks_stat[1])
p_reject = (np.array(p_values) < 0.05).sum() / len(p_values)
if p_reject < 0.95:
    rej = 'Accept'
else:
    rej = 'Reject'
pyplot.hist(np.array(burstiness_expon_r.flatten()), bins = 60,
            alpha=0.5, density=True, label = 'Exponential', color='#1f77b4')
pyplot.hist(burstiness_r.flatten(), bins = 60,
            alpha=0.5, density=True, label = 'Data', color='#ff7f0e')
ax = pyplot.gca()
ax.set_xlim([-1.0, 0.5])
ax.set_xlabel('B')
ax.set_ylabel('Density')
# Annotate figure
txt = 'p reject: %.2f\n%s\nReverse\n(High rate)' % (p_reject, rej)
ax.annotate(txt, (0.03, 0.72), xycoords = 'axes fraction', fontsize = 10)
ax.annotate('f)', (-0.23, 0.98), xycoords = 'axes fraction', fontsize=10)
ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

#########################
# Now we add plots against gamma distribution
# Seventh plot - all data against Gamma
pyplot.subplot2grid((4, 3), (2,0), colspan=1, rowspan=1)
#Do KS test sample by sample
ks_stats = []
p_values = []
for b in burstinesses.T:
    ks_stat = ks_2samp(b.flatten(), burstiness_gamma.flatten()) 
    ks_stats.append(ks_stat[0])
    p_values.append(ks_stat[1])
p_reject = (np.array(p_values) < 0.05).sum() / len(p_values)
if p_reject < 0.95:
    rej = 'Accept'
else:
    rej = 'Reject'
pyplot.hist(np.array(burstiness_gamma.flatten()), bins = 60,
            alpha=0.5, density=True, label = 'Gamma', color='slategrey')
pyplot.hist(burstinesses.flatten(), bins = 60,
            alpha=0.5, density=True, label = 'Data', color='#ff7f0e')
ax = pyplot.gca()
ax.set_xlim([-1.0, 0.5])
ax.set_xlabel('B')
ax.set_ylabel('Density')
pyplot.legend(loc=1, fontsize=9, handlelength=1.5, framealpha=0.2)
# Annotate figure
txt = 'p reject: %.2f\n%s\nAll' % (p_reject, rej)
ax.annotate(txt, (0.03, 0.77), xycoords = 'axes fraction', fontsize = 10)
ax.annotate('g)', (-0.23, 0.98), xycoords = 'axes fraction', fontsize=10)

#############
# Eighth plot - high activity rate faults against gamma
pyplot.subplot2grid((4, 3), (2,1), colspan=1, rowspan=1)
#Do KS test sample by sample
ks_stats = []
p_values = []
for b in burstiness_fast.T:
    ks_stat = ks_2samp(b.flatten(), burstiness_gamma_fast.flatten()) 
    ks_stats.append(ks_stat[0])
    p_values.append(ks_stat[1])
p_reject = (np.array(p_values) < 0.05).sum() / len(p_values)
if p_reject < 0.95:
    rej = 'Accept'
else:
    rej = 'Reject'
pyplot.hist(np.array(burstiness_gamma_fast.flatten()), bins = 60,
            alpha=0.5, density=True, label = 'Gamma', color='slategrey')
pyplot.hist(burstiness_fast.flatten(), bins = 60,
            alpha=0.5, density=True, label = 'Data', color='#ff7f0e')
ax = pyplot.gca()
ax.set_xlim([-1.0, 0.5])
ax.set_xlabel('B')
ax.set_ylabel('Density')
# Annotate figure
txt = 'p reject: %.2f\n%s\nHigh Rate' % (p_reject, rej)
ax.annotate(txt, (0.03, 0.77), xycoords = 'axes fraction', fontsize = 10)
ax.annotate('h)', (-0.23, 0.98), xycoords = 'axes fraction', fontsize=10)

############# 
# nineth plot - low activity rate faults against gamma
pyplot.subplot2grid((4, 3), (2,2), colspan=1, rowspan=1)
#Do KS test sample by sample
ks_stats = []
p_values = []
for b in burstiness_slow.T:
    ks_stat = ks_2samp(b.flatten(), burstiness_gamma_slow.flatten()) 
    ks_stats.append(ks_stat[0])
    p_values.append(ks_stat[1])
p_reject = (np.array(p_values) < 0.05).sum() / len(p_values)
if p_reject < 0.95:
    rej = 'Accept'
else:
    rej = 'Reject'
pyplot.hist(np.array(burstiness_gamma_slow.flatten()), bins = 60,
            alpha=0.5, density=True, label = 'Gamma', color='slategrey')
pyplot.hist(burstiness_slow.flatten(), bins = 60,
            alpha=0.5, density=True, label = 'Data', color='#ff7f0e')
ax = pyplot.gca()
ax.set_xlim([-1.0, 0.5])
ax.set_xlabel('B')
ax.set_ylabel('Density')
# Annotate figure
txt = 'p reject: %.2f\n%s\nLow Rate' % (p_reject, rej)
ax.annotate(txt, (0.03, 0.77), xycoords = 'axes fraction', fontsize = 10)
ax.annotate('i)', (-0.23, 0.98), xycoords = 'axes fraction', fontsize=10)
ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

#############
# Tenth plot - strike-slip faults against gamma
pyplot.subplot2grid((4, 3), (3,0), colspan=1, rowspan=1)
burstiness_gamma_ss = burstiness_gamma[indices_ss_hs]
#Do KS test sample by sample
ks_stats = []
p_values = []
for b in burstiness_ss.T:
    ks_stat = ks_2samp(b.flatten(), burstiness_gamma_ss.flatten()) 
    ks_stats.append(ks_stat[0])
    p_values.append(ks_stat[1])
p_reject = (np.array(p_values) < 0.05).sum() / len(p_values)
if p_reject < 0.95:
    rej = 'Accept'
else:
    rej = 'Reject'
pyplot.hist(np.array(burstiness_gamma_ss.flatten()), bins = 60,
            alpha=0.5, density=True, label = 'Gamma', color='slategrey')
pyplot.hist(burstiness_ss.flatten(), bins = 60,
            alpha=0.5, density=True, label = 'Data', color='#ff7f0e')
ax = pyplot.gca()
ax.set_xlim([-1.0, 0.5])
ax.set_xlabel('B')
ax.set_ylabel('Density')
# Annotate figure
txt = 'p reject: %.2f\n%s\nStrike-slip\n(High rate)' % (p_reject, rej)
ax.annotate(txt, (0.03, 0.72), xycoords = 'axes fraction', fontsize = 10)
ax.annotate('j)', (-0.23, 0.98), xycoords = 'axes fraction', fontsize=10)

#############
# Eleventh plot - normal faults against gamma
pyplot.subplot2grid((4, 3), (3,1), colspan=1, rowspan=1)
burstiness_gamma_n = burstiness_gamma[indices_n_hs]
#Do KS test sample by sample
ks_stats = []
p_values = []
for b in burstiness_n.T:
    ks_stat = ks_2samp(b.flatten(), burstiness_gamma_n.flatten()) 
    ks_stats.append(ks_stat[0])
    p_values.append(ks_stat[1])
p_reject = (np.array(p_values) < 0.05).sum() / len(p_values)
if p_reject < 0.95:
    rej = 'Accept'
else:
    rej = 'Reject'
pyplot.hist(np.array(burstiness_gamma_n.flatten()), bins = 60,
            alpha=0.5, density=True, label = 'Gamma', color='slategrey')
pyplot.hist(burstiness_n.flatten(), bins = 60,
            alpha=0.5, density=True, label = 'Data', color='#ff7f0e')
ax = pyplot.gca()
ax.set_xlim([-1.0, 0.5])
ax.set_xlabel('B')
ax.set_ylabel('Density')
# Annotate figure
txt = 'p reject: %.2f\n%s\nNormal\n(High rate)' % (p_reject, rej)
ax.annotate(txt, (0.03, 0.72), xycoords = 'axes fraction', fontsize = 10)
ax.annotate('k)', (-0.23, 0.98), xycoords = 'axes fraction', fontsize=10)

#############
# Twelfth plot - reverse faults against gamma
pyplot.subplot2grid((4, 3), (3,2), colspan=1, rowspan=1)
burstiness_gamma_r = burstiness_gamma[indices_r_hs]
#Do KS test sample by sample
ks_stats = []
p_values = []
for b in burstiness_r.T:
    ks_stat = ks_2samp(b.flatten(), burstiness_gamma_r.flatten()) 
    ks_stats.append(ks_stat[0])
    p_values.append(ks_stat[1])
p_reject = (np.array(p_values) < 0.05).sum() / len(p_values)
if p_reject < 0.95:
    rej = 'Accept'
else:
    rej = 'Reject'
pyplot.hist(np.array(burstiness_gamma_r.flatten()), bins = 60,
            alpha=0.5, density=True, label = 'Gamma', color='slategrey')
pyplot.hist(burstiness_r.flatten(), bins = 60,
            alpha=0.5, density=True, label = 'Data', color='#ff7f0e')
ax = pyplot.gca()
ax.set_xlim([-1.0, 0.5])
ax.set_xlabel('B')
ax.set_ylabel('Density')
# Annotate figure
txt = 'p reject: %.2f\n%s\nReverse\n(High rate)' % (p_reject, rej)
ax.annotate(txt, (0.03, 0.72), xycoords = 'axes fraction', fontsize = 10)
ax.annotate('l)', (-0.23, 0.98), xycoords = 'axes fraction', fontsize=10) 
ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f')) 
#############
fig.set_size_inches(w=9,h=12.) 
pyplot.tight_layout()
figname = 'combined_burstiness_hist_%s.png' % fig_comment
pyplot.savefig(figname)
figname = 'combined_burstiness_hist_%s.pdf' % fig_comment
pyplot.savefig(figname)

#########################################
# In this analysis, we implement the method of Williams et al 2019,
# except we do not sample with replacement for our chronologies,
# because we want to keep in chronological order to be consistent
# with our memory coefficient analysis.
p_values = []
d = burstinesses - burstiness_expon
for i, dd in enumerate(d):
    pos_ind = (dd > 0).astype(int)
    p_value = np.sum(pos_ind)/len(dd)
    p_values.append(p_value)
p_values = np.array(p_values)
# Get number at various p level
p_1 = np.count_nonzero(p_values < 0.01)
p_5 = np.count_nonzero(p_values < 0.05)
p_10 = np.count_nonzero(p_values < 0.1)
p_20 = np.count_nonzero(p_values < 0.2)
print('pvalues: 1, 5, 10, 20')
print(p_1, p_5, p_10, p_20)
print('Number of faults: %i' % n_faults)
print('Percentages')
print(p_1/n_faults*100, p_5/n_faults*100, p_10/n_faults*100, p_20/n_faults*100)
pyplot.clf()
pyplot.hist(p_values, bins=50, density=True)
ax = pyplot.gca()
ax.set_xlabel('p value')
ax.set_ylabel('Density')
figname = 'williams_p_value_hist_%s.png' % fig_comment
pyplot.savefig(figname)  
figname = 'williams_p_value_hist_%s.pdf' % fig_comment
pyplot.savefig(figname)

##########################3
# Print some basic stats
mean_all_covs = np.mean(np.array(covs))
print('Mean COV, all records', mean_all_covs)
mean_all_bs = np.mean(burstinesses)
print('Mean burstiness, all records', mean_all_bs)
print('Std dev burstiness, all records', np.std(burstinesses)) 
print('Mean memory coefficient, all records', np.mean(memory_coefficients))
print('Std dev memory coefficient, all records', np.std(memory_coefficients))
print('Mean burstiness, fast faults', np.mean(burstinesses[indices]))
print('Std dev burstiness, fast faults', np.std(burstinesses[indices]))
print('Mean memory coefficient, fast faults', np.mean(memory_coefficients[indices])) 
print('Std dev memory coefficient, fast faults', np.std(memory_coefficients[indices])) 

# Get alpha only for high activity rate faults
alpha_fast_faults = ie_gamma_alpha[indices]
print('Mean alpha paramater for gamma distribution', np.mean(ie_gamma_alpha))
print('Median alpha paramater for gamma distribution', np.median(ie_gamma_alpha)) 
print('Mean alpha paramater for gamma distribution, high activity rate faults',
      np.mean(alpha_fast_faults))
print('Median alpha paramater for gamma distribution, high activity rate faults', np.median(alpha_fast_faults))
# Try excluding outliers
alpha_fast_faults_exclude_outliers = alpha_fast_faults[alpha_fast_faults < 10]
alpha_all_faults_exclude_outliers = ie_gamma_alpha[ie_gamma_alpha < 10] 
print('Mean alpha paramater for gamma distribution fast faults, exclude outliers',
      np.mean(alpha_fast_faults_exclude_outliers))
print('Median alpha paramater for gamma distribution fast faults, exclude outliers',
      np.median(alpha_fast_faults_exclude_outliers)) 
print('Mean alpha paramater for gamma distribution, all faults excluding outliers',
      np.mean(alpha_all_faults_exclude_outliers))
print('Median alpha paramater for gamma distribution, all faults excluding outliers',
      np.median(alpha_all_faults_exclude_outliers)) 

#################################
# Look at events where we've modelled the open interval because it's exceptionally long
print('Open interval has been modelled for these records:', added_events)

st = set(added_events)
# Get indices of faults with added events
idx = [i for i, e in enumerate(names) if e in st]
pyplot.clf()
fig = pyplot.figure(1)
# set up subplot grid
gridspec.GridSpec(2, 2)
labels = ['Teton', 'Loma Blanca', 'Wasatch (Brigham)', 'San Andreas (Coachella)']
for j,i in enumerate(idx):
    if j < 2:
        pyplot.subplot2grid((2, 2), (0,j), colspan=1, rowspan=1)
    else:
        pyplot.subplot2grid((2, 2), (1,(j-2)), colspan=1, rowspan=1) 
    last_ie_time = all_ie_times[i][-1]
    ax = pyplot.gca()
    pyplot.hist(last_ie_time, bins=40, density=True, color='0.5', label = labels[j])
    pyplot.legend()
    ax.set_xlabel('Length of final interevent time (years)')
    ax.set_ylabel('Density')
pyplot.tight_layout()
figname = 'Added_interval_histograms.png'
pyplot.savefig(figname)    
figname = 'Added_interval_histograms.pdf'
pyplot.savefig(figname)   
