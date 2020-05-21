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
from QuakeRates.utilities.memory_coefficient import burstiness, memory_coefficient

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx, array[idx]

def merge_event_sets(event_sets, n_samples, merge_tolerance = 1,
                     order = None, start_year=-1e9):
    """Function to merge event sets.
    params event_sets: List of EventSet objects, one for each fault segement
    params n_samples: Number of chronology samples to create for each segment.
    params merge_tolerance: If the number of years between simulated event dates is 
                           less than this, and the fault segements are adjacent, then
                           the events are considered the same multi-sgement event.
    params order: List of names of each fault segment describing such that fault 
                  segments are listed from one end of the fault system to the other.
    params start_year: Trim records before this year, to deal with completeness
    """

    # Create empty EventSet object for storing merged chronology
    # Maybe not necessary - won't handle chronologies of difference lengths well
#    combined_event_set = EventSet([])
    # Loop over list of fault segment event sets to generate chronologies
    burstinesses = []
    mems = []
    chron_dict = {}
    for event_set in event_sets:
        event_set.gen_chronologies(n_samples, observation_end=2019, min_separation=1)
        print(event_set.name)
        chron_dict[event_set.name] = event_set.chronology
    print(chron_dict)
    for i, name in enumerate(order):
        if i == 0:
            chrons = chron_dict[name].T
        else:
            merged_chrons = []
            for j, chron in enumerate(chrons):
                updated_chron = []
                near_indices = []
                for k, realisation in enumerate(chron):
                    new_chron = chron_dict[name].T[j]
                    near_ind, nearest = find_nearest(new_chron, realisation)
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
                updated_chron += list(new_chron)
                updated_chron.sort()
                # Now just get events from complete record
                updated_chron = np.array(updated_chron)
                updated_chron = updated_chron[updated_chron > start_year]
                merged_chrons.append(updated_chron)
            chrons = np.array(merged_chrons)
#            print('Merged_chrons combined', chrons)
    combined_event_set = chrons
    # Now loop over merged chronologies and calculate memory and burstiness
    for chron in combined_event_set:
        ie_times = np.diff(np.array(chron))
        b = burstiness(ie_times)
        burstinesses.append(b)
        mem = memory_coefficient(ie_times)
        mems.append(mem)
#        print(ie_times)
#        print(np.mean(ie_times))
#        print(np.std(ie_times))
    print('Combined event set', combined_event_set)
    print('Burstiness', burstinesses)
    print('Memory coefficients', mems)
    print('Mean B', np.mean(burstinesses))
    print('Mean M', np.mean(mems))
    b_ub = np.percentile(burstinesses, 97.5)
    b_lb = np.percentile(burstinesses, 2.5)
    m_ub = np.percentile(mems, 97.5)
    m_lb = np.percentile(mems, 2.5)
    print('b_ub', b_ub)
    print('b_lb', b_lb)
    print('m_ub', m_ub)
    print('m_lb', m_lb)
    return combined_event_set

if __name__ == "__main__":
    filepath = '../params'
    n_samples = 100
    start_year = 0 # Year (positive AD, negative BC) for which record is complete
    # for all segements
    faulting_styles = ['all']
    tectonic_regions = ['all']
    min_number_events = 5
    param_file_list = glob(os.path.join(filepath, 'NorthA*.txt'))
    print(param_file_list)
    try:
        # Remove northern segemnt studies from San Andreas
        param_file_list.remove('../params/SanAndreasVedanta_Zhang_2005unpub_simple.txt')
        param_file_list.remove('../params/SanAndreasMendocino_Merritts_1996_simple.txt')
    except ValueError:
        pass
    try:
        # Exclude Taybeh site for longer record 
        param_file_list.remove('../params/DeadSeaTaybeh_Lefevre_2018_simple.txt')
    except ValueError:
        pass
    try:
        # Remove as only events until 600 AD
        param_file_list.remove('../params/NorthAnatolianElmacik_Fraser_2010_simple.txt')
        # Remove as distant and relatively short record
        param_file_list.remove('../params/NorthAnatolianElmacik_Fraser_2010_simple.txt') 
    except:
        ValueError
        pass
    print(param_file_list)
    
#    sys.exit()
    names, event_sets, event_certainties, num_events = \
        get_event_sets(param_file_list, tectonic_regions,
                       faulting_styles, min_number_events)
    # Give name attribute to EventSet objects
    print(names) # Test just using order the faults are read in as
#    sys.exit()
    SanAndreas_order = ['SanAndreasCoachella', 'SanAndreasThousandPalms', 'SanAndreasBurro', 'SanAndreasPittman', 'SanAndreasWrightwood', 'SanAndreasPalletCk', 'SanAndreasBigBend', 'SanAndreasCarizzo'] #, 'SanAndreasVedanta'] # South to north
    start_year_sa = 600
#    DeadSea_order = ['DeadSeaQatar', 'DeadSeaTaybeh', 'DeadSeaJordan','DeadSeaBeteiha', 'DeadSeaYammouneh']
    # South to north
    DeadSea_order = ['DeadSeaQatar', 'DeadSeaJordan','DeadSeaBeteiha', 'DeadSeaYammouneh'] # Exclude Taybeh site for longer record 
    start_year_ds = -3000
    # East to west
    NorthAnatolian_order = ['NorthAnatolianCukurcimen', 'NorthAnatolianYaylabeli', 'NorthAnatolianGunalan', 'NorthAnatolianLakeLadik']
    start_year_na = -1400
    for i, name in enumerate(names):
        event_sets[i].name = name
    merge_event_sets(event_sets, n_samples, merge_tolerance=2, order = NorthAnatolian_order, start_year=start_year_na)
    
