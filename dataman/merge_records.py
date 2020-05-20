""""Contains functions to merge event sets from multiple records, 
e.g. to look at statistics for fault systems rather than single faults

Jonathan Griffin
University of Otago
"""

import os
import numpy as np
from glob import glob
from QuakeRates.dataman.event_dates import EventSet 
from QuakeRates.dataman.parse_oxcal import parse_oxcal
from QuakeRates.dataman.parse_params import parse_param_file, \
    get_event_sets, file_len 


def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx, array[idx]

def merge_event_sets(event_sets, n_samples, merge_tolerance = 1, order = None):
    """Function to merge event sets.
    params event_sets: List of EventSet objects, one for each fault segement
    params n_samples: Number of chronology samples to create for each segment.
    params merge_tolerance: If the number of years between simulated event dates is 
                           less than this, and the fault segements are adjacent, then
                           the events are considered the same multi-sgement event.
    params order: List of names of each fault segment describing such that fault 
                  segments are listed from one end of the fault system to the other.
    """

    # Create empty EventSet object for storing merged chronology
    combined_event_set = EventSet([])
    # Loop over list of fault segment event sets to generate chronologies
    chron_dict = {}
    for event_set in event_sets:
        event_set.gen_chronologies(n_samples, observation_end=2019, min_separation=1)
        print(event_set.name)
        chron_dict[event_set.name] = event_set.chronology
    print(chron_dict)
    for i, name in enumerate(order):
        if i == 0:
            chrons = chron_dict[name].T
            print('chrons', chrons)
        else:
            merged_chrons = []
            for j, chron in enumerate(chrons):
                print('chron',chron)
                updated_chron = []
                near_indices = []
                for k, realisation in enumerate(chron):
                    print(realisation)
#                    print('chron_dict[name]',chron_dict[name].T[j])
#                    for l, new_chron in enumerate(chron_dict[name][j]):
                    new_chron = chron_dict[name].T[j]
                    print('new_chron', new_chron)
                    near_ind, nearest = find_nearest(new_chron, realisation)
                    print('nearest', nearest)
                    if abs(nearest - realisation) < merge_tolerance:
                        if near_ind not in near_indices: # Ensure we don't use the same event twice
                            updated_chron.append(np.mean([nearest, realisation]))
                            near_indices.append(near_ind) 
                        else:
                            updated_chron.append(realisation)
                    else:
                        updated_chron.append(realisation)
                        # Now add in all other events that haven't been merged, and then sort
                new_chron = np.delete(new_chron, near_indices)
                print('updated_chron', updated_chron)
                print(new_chron)
                updated_chron += list(new_chron)
                updated_chron.sort()
                print('Original chron', chron)
                print('Second chron', chron_dict[name].T[j])
                print ('Merged chron', updated_chron)
                merged_chrons.append(updated_chron)
            chrons = np.array(merged_chrons).T
            print('Merged_chrons combined', chrons)
#            sys.exit()                   
    combined_event_set = chrons
    print(combined_event_set)
    return combined_event_set

if __name__ == "__main__":
    filepath = '../params'
    n_samples = 100
    start_year = 0 # Year (positive AD, negative BC) for which record is complete
    # for all segements
    faulting_styles = ['all']
    tectonic_regions = ['all']
    min_number_events = 6
    param_file_list = glob(os.path.join(filepath, 'SanAndreas*.txt')) 
    names, event_sets, event_certainties, num_events = \
        get_event_sets(param_file_list, tectonic_regions,
                       faulting_styles, min_number_events)
    # Give name attribute to EventSet objects
    print(names) # Test just using order the faults are read in as
    for i, name in enumerate(names):
        event_sets[i].name = name
    merge_event_sets(event_sets, n_samples, merge_tolerance=5, order = names)
    
