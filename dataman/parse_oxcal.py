"""Parse in csv output from ocxcal and extract the relevant date information

Jonathan Griffin
January 2019
"""

import os, sys
import csv
from QuakeRates.dataman.event_dates import EventDate

def parse_oxcal(filename, keys):
    """Parse in csv output file from OxCal and extract data as desired.
    :param filename: string of path to input file
    :param keys: dictionary specifying which part of the OxCal file is
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
                        date_dict[row['name']].append(row['value'])
                        prob_dict[row['name']].append(row['probability'])
    print(date_dict['I'])
    print(prob_dict['I'])


if __name__ == "__main__":
    filename = '../data/Xorkoli_Altyn_Tagh_Yuan_2018.csv'
    key_dict = {'I': ['Calculate', 'posterior'],
                'H': ['Calculate', 'posterior']} 
    date_dict = {}                                                                                 
    prob_dict = {}
    for key in key_dict:
        date_dict[key] = []
        prob_dict[key] = []
    parse_oxcal(filename, key_dict)
    
    
    
