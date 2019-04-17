"""Parse in csv output from ocxcal and extract the relevant date information

Jonathan Griffin
January 2019
"""

import os, sys
import csv
from QuakeRates.dataman.event_dates import EventDate, EventSet

def parse_oxcal(filename, key_dict, event_order=None):
    """Parse in csv output file from OxCal and extract data as desired.
    :param filename: string of path to input file
    :param key_dict: dictionary specifying which part of the OxCal file is
    desired, e.g. we may only want the caclulated event dates, not the
    raw C14 dates
    :param event_order: Ordered list of key for key_dict specifying the
    chronological order of the events (forward in time, i.e. oldest event
    is first).
    """

    # Dicts for storing dates and probabilities
    date_dict = {}
    prob_dict = {}
    for key in key_dict:                                                                                
        date_dict[key] = []
        prob_dict[key] = []
    with open(filename) as csv_file:
        csv_reader = csv.DictReader(csv_file)
        line_count = 0
        for row in csv_reader:
            if line_count == 0:
                print(f'Column names are {", ".join(row)}')
                line_count += 1
            else:
                if row['name'] in key_dict:
                    #print(row['op'])
                    #print(row['type'])
                    #print(key_dict[row['name']][0])
                    #print(key_dict[row['name']][1])
                    if row['op'] ==  key_dict[row['name']][0] and \
                       row['type'] == key_dict[row['name']][1]:
                        date_dict[row['name']].append(float(row['value']))
                        prob_dict[row['name']].append(float(row['probability']))

    event_list = []
    if event_order is not None:
        keys = event_order
    else:
        keys = list(key_dict.keys())
    for key in keys:
        event = EventDate(key, key_dict[key][0], key_dict[key][1])
        event.add_dates_and_probs(date_dict[key], prob_dict[key])
        event_list.append(event)
        
 #   print(event_list)
 #   print(event_list[0].dates)
 #   print(event_list[0].probabilities)
            
    return event_list

if __name__ == "__main__":
#    filename = '../data/Xorkoli_Altyn_Tagh_Yuan_2018.csv'
    filename = '../../OxCal/akatore/Akatore4event_Bdy/Akatore4eventBdy_output.csv'
    # key_dict for parseing OxCal output. The key specifies the event
    # name, while the two parameters within the list specify that we
    # want the calculated dates from the posterior distribution.
#    key_dict = {'I': ['Calculate', 'posterior'],
#                'H': ['Calculate', 'posterior'],
#                'G': ['Calculate', 'posterior'],
#                'F': ['Calculate', 'posterior'],
#                'E': ['Calculate', 'posterior'],
#                'D': ['Calculate', 'posterior'],
    key_dict = {'D': ['Calculate', 'posterior'], 
                'C': ['Calculate', 'posterior'],
                'B': ['Boundary', 'posterior'],
                'A': ['Boundary', 'posterior']}
    
    # The event order is defined separately from key_dict as may want to
    # try different orderings etc
    #    event_order = ['I', 'H', 'G', 'F', 'E', 'D', 'C', 'B', 'A']
    event_order = ['D', 'C', 'B', 'A'] 
# Parse OxCal file
    events = parse_oxcal(filename, key_dict, event_order)
    event_set = EventSet(events)
#    event_set.cov()
    n_samples = 10000
    event_set.gen_chronologies(n_samples)

    event_set.calculate_cov()
    figname = '%s_%i_chronologies.png' % (filename[:-3], n_samples)
    chron_filename = '%s_%i_chronologies.csv' % (filename[:-3], n_samples) 
    event_set.plot_chronology(figname, normalise=False)
    event_set.write_chronology(chron_filename)
    # for event in events:
   #     fig_filename = 'event_' + event.id + '_pdf.png'
   #     event.plot_date_pdf(fig_filename)
   #     event.random_sample(100, plot=True)

        
