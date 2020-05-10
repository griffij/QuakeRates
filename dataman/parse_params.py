"""Contains basic functions for reading parameter file and creating 
initial event_sets from raw data
"""

import os, sys
import ast
import numpy as np
from QuakeRates.dataman.event_dates import EventSet
from QuakeRates.dataman.parse_oxcal import parse_oxcal
from QuakeRates.dataman.parse_age_sigma import parse_age_sigma

def file_len(fname):
    """
    Get the number of lines in a file
    """
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1 

def parse_param_file(param_file_name):
    """
    Read parameters for fault from parameter file
    """
    params={}
    with open(param_file_name) as f_in:
        for line in f_in:
            var, value = line.strip().split('=')
            params[var.strip()] = ast.literal_eval(value.strip())
    return params

def get_event_sets(param_file_list, tectonic_regions,
                   faulting_styles, min_number_events):
    """Function to take parameter file list, read data
    and return list of event_sets and fautl names
    """
    names = [] # Store names
    event_sets = [] # Store EventSet objects
    event_certainties = [] # Store list of 0/1 for uncertain/certain events
    num_events = []
    for param_file in param_file_list:
        name = param_file.split('/')[-1].split('_')[0]
        print(name)
        params = parse_param_file(param_file)
        print(params)
        # Now we want to take subsets of the data based on parameters
        if params['tectonic_region'] not in tectonic_regions:
            if tectonic_regions[0] == 'all':
                pass
            else:
                continue
        if params['faulting_style'] not in faulting_styles:
            if faulting_styles[0] == 'all':
                pass
            else:
                continue
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
            if len(params['event_order']) >= min_number_events:
                num_events.append(len(params['event_order']))
                events = parse_oxcal(params['filename'], params['events'],
                                     params['event_order'])
                event_set = EventSet(events)
                try:
                    event_certainty = np.array(params['event_certainty'])
                except KeyError: # If not specified, assume all events are certain
                    event_certainty = np.ones(len(events))
                event_certainties.append(event_certainty)
            else:
                continue
        elif params['chron_type'] == 'Age2Sigma':
            # Assume single header line
            if (file_len(params['filename']) - 1) >= min_number_events:
                num_events.append((file_len(params['filename']) - 1))
                events, event_certainty = parse_age_sigma(params['filename'],
                                                          params['sigma_level'],
                                                          params['event_order'])
                event_set = EventSet(events)
                event_certainties.append(event_certainty)
            else:
                continue
        else:
            msg = 'Unknown form of chron_type defined in ' + param_file
            raise Exception(msg)
        names.append(name)
        # Add some parameters to the event set
        event_set.faulting_style = params['faulting_style']
        event_set.tectonic_region = params['tectonic_region']
        event_sets.append(event_set)
    return names, event_sets, event_certainties, num_events
