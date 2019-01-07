"""Parse in csv output from ocxcal and extract the relevant date information

Jonathan Griffin
January 2019
"""

import os, sys
import csv
from QuakeRates.dataman.event_dates import EventDate

key_dict = {'I': ['Calculate', 'posterior']

def parse_oxcal(filename, keys):
    """Parse in csv output file from OxCal and extract data as desired.
    :param filename: string of path to input file
    :param keys: dictionary specifying which part of the OxCal file is
    desired, e.g. we may only want the caclulated event dates, not the
    raw C14 dates
    """

    with open(filename) as csvfile:
        csv_reader = csv.DictReader(csv_file)
        line_count = 0
        for row in csv_reader:
            if line_count == 0:
                print(f'Column names are {", ".join(row)}')
                line_count += 1
            else:
                if csv_reader['name'] in key_dict:
                    
    
    
    
