"""Script to summarise metadata from param files into a single file
"""

from QuakeRates.dataman.parse_params import parse_param_file, get_event_sets
from glob import glob
import csv
import numpy as np

paramfiles = glob('../params/*.txt')
# Remove unwanted records      
paramfiles.remove('../params/EmersonNorth_Rockwell_2000_simple.txt')

f_out = open('metadata_summary.csv', 'w')
writer = csv.writer(f_out)
header = ['Name', 'Tectonic_region', 'Faulting_style', 'Method', 'Longitude',
          'Latitude', 'Location_accuracy']
names, event_sets, event_certainties, num_events, \
    tect_regions, fault_styles, methods, \
    longitudes, latitudes, location_accuracies = get_event_sets(paramfiles, ['all'], ['all'], 1,
                                                return_full=True)
sort_index = np.argsort(names)
writer.writerow(header)
for i in sort_index:
#for i, event_set in enumerate(event_sets):
    if num_events[i] >= 5:# and longitudes[i] != 'NULL':
        metadata = [names[i], tect_regions[i], fault_styles[i], methods[i],
                    longitudes[i], latitudes[i], location_accuracies[i]]
        writer.writerow(metadata)

f_out.close()
