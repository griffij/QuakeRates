"""PlAot exOAample daOAtasets in memory-burstiness spOAace
"""

import os, sys
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
from scipy.stats import expon, gamma, weibull_min  
from QuakeRates.utilities.memory_coefficient import memory_coefficient, burstiness

np.random.seed(123)
#ie_times = [100, 110, 120, 150, 200, 210, 220, 250]
ie_times = expon(scale=100).rvs(size=(10))
#zeros = np.zeros(len(ie_times))
#ones = np.ones(len(ie_times))

mems = []
bursts = []
labels = []

# Set up plot
plt.clf() 
fig = plt.figure(1)#, figsize=(6., 24.))
# set up subplot grid
gridspec.GridSpec(6, 2)
#First plot
plt.subplot2grid((6, 2), (0,0), colspan=1, rowspan=1)
#ax = plt.subplot(111)
ax = plt.gca()
ax.set_xlim([0, sum(ie_times)+50])
ax.set_ylim([0, 1]) 
time_sum = 0
for i, ie_time in enumerate(ie_times):
    time_sum += ie_time
    plt.plot([time_sum, time_sum], [0, 1], c='k')
#plt.show()
ax.annotate('a)', (-0.13, 0.9), xycoords = 'axes fraction', fontsize=10)
ax.set_title(r'Poisson $(\lambda = 100)$', fontsize=10)
burst = burstiness(ie_times)
memory = memory_coefficient(ie_times)
print(burst)
print(memory)
mems.append(memory)
bursts.append(burst)
labels.append('a')

# Next subplot
plt.subplot2grid((6, 2), (1,0), colspan=1, rowspan=1)
sorted_ie = np.sort(ie_times)
sorted_ie = np.concatenate(([sorted_ie[-2]], [sorted_ie[-3]], sorted_ie[0:-3],[sorted_ie[-1]]))
time_sum = 0
for i, ie_time in enumerate(sorted_ie):
    time_sum += ie_time
    plt.plot([time_sum, time_sum], [0, 1], c='k')
#plt.show()
ax = plt.gca()
ax.set_xlim([0, sum(ie_times)+50])
ax.set_ylim([0, 1])
ax.set_title('Shuffled, positive $M$', fontsize=10)
ax.annotate('b)', (-0.13, 0.9), xycoords = 'axes fraction', fontsize=10) 
burst = burstiness(sorted_ie)
memory = memory_coefficient(sorted_ie)
print(burst)
print(memory)
mems.append(memory)
bursts.append(burst)
labels.append('b')


# Next subplot
plt.subplot2grid((6, 2), (2,0), colspan=1, rowspan=1)
ies = np.sort(ie_times)
sorted_ie = np.array([ies[0], ies[-1], ies[1], ies[-2], ies[2], ies[-3], ies[3], ies[-4], ies[4], ies[-5]])
#sorted_ie = np.concatenate(([sorted_ie[-2]], [sorted_ie[-3]], sorted_ie[0:-3],[sorted_ie[-1]]))
time_sum = 0
for i, ie_time in enumerate(sorted_ie):
    time_sum += ie_time
    plt.plot([time_sum, time_sum], [0, 1], c='k')
#plt.show()
ax = plt.gca()
ax.set_ylim([0, 1]) 
ax.annotate('c)', (-0.13, 0.9), xycoords = 'axes fraction', fontsize=10) 
ax.set_title('Shuffled, negative $M$', fontsize=10)
ax.set_xlabel('Time (years)')
burst = burstiness(sorted_ie)
memory = memory_coefficient(sorted_ie)
print(burst)
print(memory)
mems.append(memory)
bursts.append(burst)
labels.append('c')

# Next subplot
plt.subplot2grid((6, 2), (0,1), colspan=1, rowspan=1)
# Based on alpine fault
af = [300, 500, 705, 860, 1155, 1364, 1777]
af_ie_times = []
for i, date in enumerate(af):
    plt.plot([date, date], [0, 1], c='k')
    if i > 0:
        ie_times = date - af[i-1]
        af_ie_times.append(ie_times)
ax = plt.gca()
ax.set_ylim([0, 1])
ax.annotate('d)', (-0.13, 0.9), xycoords = 'axes fraction', fontsize=10)
ax.set_title('Alpine Fault (Hokuri Creek)', fontsize=10)
burst = burstiness(af_ie_times)
memory = memory_coefficient(af_ie_times)
print(burst)
print(memory)
mems.append(memory)
bursts.append(burst)
labels.append('AF')

# Next subplot
plt.subplot2grid((6, 2), (1,1), colspan=1, rowspan=1)
# Based on Mentawai
sum_m = [1314, 1350, 1388, 1569, 1597, 1613, 1631, 1658, 1703, 1797, 1833, 2007, 2010]
sum_m_ie_times = []
for i, date in enumerate(sum_m):
    plt.plot([date, date], [0, 1], c='k')
    if i > 0:
        ie_times = date - sum_m[i-1]
        sum_m_ie_times.append(ie_times)
ax = plt.gca() 
ax.set_ylim([0, 1])
ax.annotate('e)', (-0.13, 0.9), xycoords = 'axes fraction', fontsize=10) 
ax.set_title('Sunda Arc (Mentawai Segment)', fontsize=10) 
burst = burstiness(sum_m_ie_times)
memory = memory_coefficient(sum_m_ie_times)
print(burst)
print(memory)
mems.append(memory)
bursts.append(burst)
labels.append('SA')

# Next subplot
plt.subplot2grid((6, 2), (2,1), colspan=1, rowspan=1)
# Based on Cadell fault
cf = [-4500000, -2000000, -1000000, -70000, -62500, -55000, -45000, -38500, -32000]
cf_ie_times = []
for i, date in enumerate(cf):
    plt.plot([date, date], [0, 1], c='k')
    if i > 0:
        ie_times = date - cf[i-1]
        cf_ie_times.append(ie_times)
ax = plt.gca()
ax.set_ylim([0, 1])
ax.annotate('f)', (-0.13, 0.9), xycoords = 'axes fraction', fontsize=10)
ax.set_title('Cadell Fault', fontsize=10)
ax.set_xlabel('Year')
plt.xticks([-4000000, -2000000, 0])
burst = burstiness(cf_ie_times)
memory = memory_coefficient(cf_ie_times)
print(burst)
print(memory)
mems.append(memory)
bursts.append(burst)
labels.append('CF')

# Now plot M_B phase diagram
plt.subplot2grid((6, 2), (3,0), colspan=2, rowspan=3)
#markers = ['s', 's', 's', '^', '^', '^']
plt.scatter(mems[0:3], bursts[0:3], marker = 's', c='k', s=20)
plt.scatter(mems[3:6], bursts[3:6], marker = '^', c='k', s=20)
#plt.show()
ax = plt.gca()
ax.set_xlabel('M')
ax.set_ylabel('B')
ax.set_xlim([-1, 1])
ax.set_ylim([-1, 1])
ax.annotate('g)', (-0.2, 0.98), xycoords = 'axes fraction', fontsize=10) 
# Add y = 0, x=0 lines
plt.plot([0,0],[-1, 1], linestyle='dashed', linewidth=1, c='0.5')
plt.plot([-1,1],[0, 0], linestyle='dashed', linewidth=1, c='0.5')
for i, txt in enumerate(labels):
    ax.annotate(txt,
                (mems[i]-0.05, bursts[i]+0.05),
                fontsize=10)
ax.annotate('Quasi-periodic', (-0.25, -0.8), fontsize=10, fontstyle='italic')
ax.annotate('Elastic \nrebound \ndominated', (-0.85, -0.70), fontsize=10, fontstyle='italic')
ax.annotate('Supercycles', (-0.25, 0.12), fontsize=10, fontstyle='italic')
ax.annotate('Poisson', (-0.25, -0.12), fontsize=10, fontstyle='italic')
ax.annotate('Clusters', (0.25, 0.25), fontsize=10, fontstyle='italic') 
ax.set_aspect('equal')

fig.tight_layout(pad=1.5, w_pad=1.5, h_pad=-0.4)
fig.set_size_inches(w=6,h=8) 
fig.savefig('B_M_phase_examples.png')#, bbox_inches='tight',pad_inches = 0)
