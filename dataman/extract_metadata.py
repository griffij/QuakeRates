"""Script to summarise metadata from param files into a single file
"""

from QuakeRates.dataman.parse_params import parse_param_file, get_event_sets
from glob import glob
import csv

paramfiles = glob('../params/*.txt')
# Remove unwanted records      
paramfiles.remove('../params/EmersonNorth_Rockwell_2000_simple.txt')

f_out = open('metadata_summary.csv', 'w')
writer = csv.writer(f_out)
header = ['Name', 'Tectonic_region', 'Faulting_style']
names, event_sets, event_certainties, num_events, \
    tect_regions, fault_styles = get_event_sets(paramfiles, ['all'], ['all'], 1)
for i, event_set in enumerate(event_sets):
    if num_events[i] >= 5:
        metadata = [names[i], tect_regions[i], fault_styles[i]]
        writer.writerow(metadata)

f_out.close()
