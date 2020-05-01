"""Calculate conditional probability of a short interevent 
time being followed by another short interevent time, compared
with the unconditional probability.
This is used to test whether fault records have memory
"""

import os, sys
from glob import glob
import numpy as np
import matplotlib.pyplot as plt
from QuakeRates.dataman.parse_params import parse_param_file, \
    get_event_sets

# Define parameter files
filepath = '../params'
param_file_list = glob(os.path.join(filepath, '*.txt'))
n_samples = 500  # Number of Monte Carlo samples of the eq chronologies
half_n = int(n_samples/2)

plot_dir = './plots_conditional_probs'
if not os.path.exists(plot_dir):
    os.makedirs(plot_dir)
# Define subset to take
#faulting_styles = ['Reverse']
#faulting_styles = ['Normal']
#faulting_styles = ['Strike_slip']
faulting_styles = ['all']
tectonic_regions = ['all']
#tectonic_regions = ['Plate_boundary_master', 'Plate_boundary_network']
min_number_events = 10

names, event_sets, event_certainties, num_events = \
    get_event_sets(param_file_list, tectonic_regions,
                   faulting_styles, min_number_events)

# Now loop over paleo-earthquake records
for i, event_set in enumerate(event_sets):
    # Generate some chronologies
    event_set.gen_chronologies(n_samples, observation_end=2019, min_separation=1)
    print(num_events[i])
    event_set.calculate_cov() # Calculate  interevent times and mean as part of this
    # Lists to store results
    uncond_probs = []
    cond_probs = []
    for j, sample in enumerate(event_set.interevent_times.T):
        num_less_mean = len(np.argwhere(sample < event_set.means[j]))
        uncond_prob_less_mean = num_less_mean/event_set.num_events
        count_short = 0
        for k, ie_time in enumerate(sample):
            if k==0:
                ie_time_0 = ie_time
            else:
                if ie_time < event_set.means[i] and \
                   ie_time_0 < event_set.means[i]:
                    count_short += 1
                ie_time_0 = ie_time
        cond_prob_less_mean = count_short/num_less_mean
        uncond_probs.append(uncond_prob_less_mean)
        cond_probs.append(cond_prob_less_mean)
    print(uncond_probs)
    print(cond_probs)
    uncond_probs = np.array(uncond_probs)
    cond_probs = np.array(cond_probs)
    probs_ratio = cond_probs/uncond_probs
    print(probs_ratio)
    plt.clf()
    plt.hist(probs_ratio, bins = 10, facecolor='0.6',
             edgecolor='0.2', density=True)
    figname = 'conditional_prob_ratio_histogram_%s.png' % names[i]
    fig_filename = os.path.join(plot_dir, figname)
    plt.savefig(fig_filename) 
