"""Parse data in the format:
Age Uncertainty
Sample data assumning a normal distrubtion with mean defined by Age and sigma 
defined by Uncertainty
"""

import numpy as np
from scipy.stats import norm
from QuakeRates.dataman.event_dates import EventDate, EventSet

def parse_age_sigma(filename, sigma_level, event_order, truncation=3,
                    delimiter=None, header = 1):
    """Parse a text file containing a list of event ages and associated
    uncertainties
    :param filename: string of path to input file 
    :param sigma_level: Number of sigmas represented by the uncertainty 
        columm
    :param event_order: String, 'Forwards' or 'Backwards' in time. I.e. if
        'Forwards', oldest event is in the first row of the file.
    :param truncation: Number of sigma levels to sample from
    :param delimiter: Delimiter of input text file.
    :param header: Number of header lines to discard
    """

    event_list = []
    data = np.genfromtxt(filename, delimiter=delimiter, skip_header=header)
    # We want time to be running forwards
    if event_order == 'Backwards':
        data = np.flip(data, axis=0)
    sigmas = data[:,1]/sigma_level #Convert, e.g. 2 sigma to 1 sigma
    for i,mean_age in enumerate(data[:,0]):
        event_id = i
        # Special case of zero uncertainty
        if sigmas[i]==0:
            ages = np.array([mean_age])
            probs = np.array([1.])
        else:
            ages = np.arange(mean_age-truncation*sigmas[i],
                             mean_age+truncation*sigmas[i]+1, 1)
            probs = norm.pdf(ages, mean_age, sigmas[i])
            # Normalise probs due to truncation of distribution
            probs = probs/sum(probs)
        event = EventDate(event_id, 'manual', 'age_sigma')
        event.add_dates_and_probs(ages, probs)
        print(event.dates)
        print(event.probabilities)
        event_list.append(event)
    return event_list

if __name__ == "__main__":
    filename = '../data/SanJacinto_Rockwell_2015_simple.txt'
    event_list = parse_age_sigma(filename, sigma_level=2, event_order='Backwards')
    event_set = EventSet(event_list)
    print(event_list)
    print(event_set)
    n_samples = 10000                                                                                            
    event_set.gen_chronologies(n_samples)
    event_set.calculate_cov()
    event_set.cov_density()
