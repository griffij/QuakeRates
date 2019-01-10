"""Parse in csv output from ocxcal and extract the relevant date information

Jonathan Griffin
January 2019
"""

import os, sys
import csv
from QuakeRates.dataman.event_dates import EventDate

def parse_oxcal(filename, key_dict):
    """Parse in csv output file from OxCal and extract data as desired.
    :param filename: string of path to input file
    :param key_dict: dictionary specifying which part of the OxCal file is
    desired, e.g. we may only want the caclulated event dates, not the
    raw C14 dates
    """

    with open(filename) as csv_file:
        csv_reader = csv.DictReader(csv_file)
        line_count = 0
        for row in csv_reader:
            if line_count == 0:
                print(f'Column names are {", ".join(row)}')
                line_count += 1
            else:
                if row['name'] in key_dict:
                    print(row['op'])
                    print(row['type'])
                    print(key_dict[row['name']][0])
                    print(key_dict[row['name']][1])
                    if row['op'] ==  key_dict[row['name']][0] and \
                       row['type'] == key_dict[row['name']][1]:
                        date_dict[row['name']].append(float(row['value']))
                        prob_dict[row['name']].append(float(row['probability']))

    event_list = []
    for key in key_dict:
        event = EventDate(key, key_dict[key][0], key_dict[key][1])
        event.add_dates_and_probs(date_dict[key], prob_dict[key])
        event_list.append(event)
    print(event_list)
    print(event_list[0].dates)
    print(event_list[0].probabilities)

    return event_list

if __name__ == "__main__":
    filename = '../data/Xorkoli_Altyn_Tagh_Yuan_2018.csv'
    key_dict = {'I': ['Calculate', 'posterior'],
                'H': ['Calculate', 'posterior']} 
    date_dict = {}                                                                                 
    prob_dict = {}
    for key in key_dict:
        date_dict[key] = []
        prob_dict[key] = []
    events = parse_oxcal(filename, key_dict)
    for event in events:
        fig_filename = 'event_' + event.id + '_pdf.png'
        event.plot_date_pdf(fig_filename)
        event.random_sample(10, plot=True)
    
