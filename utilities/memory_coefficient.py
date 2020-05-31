"""Function to calcualte memory coefficients
"""

import numpy as np
from scipy.stats import spearmanr

def burstiness(interevent_times):
    """Calculate burstiness - See:
    Goh, K.-I. and A.-L. Barabási (2008). Burstiness and memory in 
    complex systems, Europhysics Lett., 81, no. 4, 48002.
    """
    stds = np.std(interevent_times, axis=0)
    means = np.mean(interevent_times, axis=0)
    burstiness = (stds - means)/(stds + means)
    return(burstiness)
    
def memory_coefficient(interevent_times):
    """ Calculate the memory coeeficient of 
    Goh, K.-I. and A.-L. Barabási (2008). Burstiness and memory in 
    complex systems, Europhysics Lett., 81, no. 4, 48002.
    """
    #print(self.interevent_times)
    m1 = np.mean(interevent_times[0:-1], axis=0)
    m2 = np.mean(interevent_times[1:], axis = 0)
    s1 = np.std(interevent_times[0:-1], axis=0)
    s2 = np.std(interevent_times[1:], axis = 0)
    s1s2 = s1*s2
    # Loop over chronologies
    ie_min_m1 = interevent_times - m1
    # Remove last element
    ie_min_m1 = ie_min_m1[0:-1]
    ie_min_m2 = interevent_times - m2
    # Remove first element
    ie_min_m2 = ie_min_m2[1:]
    numerator = ie_min_m1 * ie_min_m2
    numerator = np.sum(numerator, axis=0)
    sum_term = numerator/s1s2
    num_events = len(interevent_times)
    mem_coef = sum_term * (1/(num_events - 1))
#    print('mem_coef', mem_coef)
#    mean_mem_coef = np.mean(mem_coef)
#    print('Mean memory coefficient', mean_mem_coef)
#    memory_lb = np.percentile(mem_coef, 2.5)
#    memory_ub = np.percentile(mem_coef, 97.5)
    return(mem_coef)

def memory_spearman_rank_correlation(interevent_times):
    """ Calculate alternative memory coefficient using Spearman rank
    correlation. This may avoid some biases in the calculation. 
    See: Schleiss, M. and J. A. Smith (2016). Two Simple Metrics for 
    Quantifying Rainfall Intermittency: The Burstiness and Memory of
    Interamount Times, J. Hydrometeorol., 17, no. 1, 421–436.
    """
    # Lag-1
    a = interevent_times[0:-1].T
    b = interevent_times[1:].T
    # Lag-2
    a2 = interevent_times[0:-2].T
    b2 = interevent_times[2:].T 
    #        print(a, b)
    rhos = []
    rhos2 = []
    #        .rhos, .pvalues = spearmanr(a, b)
    for i, ie_t in enumerate(a):
        #            print(ie_t.T)
        #            print(b[i].T)
        rho, pvalue = spearmanr(ie_t.T, b[i].T)
        #            print('rho', rho)
        rhos.append(rho)
    for i, ie_t in enumerate(a2):
        rho, pvalue = spearmanr(ie_t.T, b2[i].T)
        rhos2.append(rho)
    rhos = np.array(rhos)
    rhos2 = np.array(rhos2)
    #        print(.rhos, .pvalues)
    mean_rho = np.mean(rhos)
    mean_rho2 = np.mean(rhos2)
    print('Mean Spearman Rank coefficient', mean_rho)
    print('Mean Spearman Rank coefficient', mean_rho2)
    rho_lb = np.percentile(rhos, 2.5)
    rho_ub = np.percentile(rhos, 97.5)
    rho2_lb = np.percentile(rhos2, 2.5)
    rho2_ub = np.percentile(rhos2, 97.5) 
    #   za     print(.rho_lb, .rho_ub)
