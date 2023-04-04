"""Calculate basic uncertainties to do with event accuracy
"""

import numpy as np
from glob import glob
datafiles = ['../chronologies_all_final/AlpineHokuriCkSouthWestland_10000_chronologies.csv',
            '../chronologies_all_final/AlaskaPWSCopper_10000_chronologies.csv']

datafiles = glob('../chronologies_all_final/*.csv')

all_rel_uncerts = []
all_age_uncerts = []
for datafile in datafiles:
    data = np.genfromtxt(datafile, delimiter=',')
#    print(data)
    means = np.mean(data, axis=0)
    mean_ages = 2020 - means
#    print(means)
#    print(mean_ages)
    # Mean value for each interevent time
    interevent_times = np.mean(np.diff(data, axis=1), axis=0)
    interevent_times_std = np.std(np.diff(data, axis=1), axis=0) 
#    print('Interevent times', interevent_times)
    std_devs = np.std(data, axis=0)
#    print(std_devs)
#    rel_uncert = (std_devs/mean_ages) * 100
#    rel_uncert = (np.mean(std_devs)/np.mean(interevent_times))*100
    rel_uncert = np.mean(interevent_times_std/interevent_times)*100
    print(datafile)
    print('Relative uncertainty', rel_uncert)
    all_rel_uncerts.append(rel_uncert)

mean_rel_uncert = np.mean(all_rel_uncerts)
std_dev_rel_uncert = np.std(all_rel_uncerts)
print('mean_rel_uncert', mean_rel_uncert)
print('std_dev_rel_uncert', std_dev_rel_uncert)
print('min uncert', np.min(all_rel_uncerts))
print('max uncert', np.max(all_rel_uncerts))
