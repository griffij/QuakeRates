"""Plot example datasets in memory-burstiness space
"""

import os, sys
from matplotlib import pyplot as plt
import numpy as np
from scipy.stats import expon, gamma, weibull_min  
from QuakeRates.utilities.memory_coefficient import memory_coefficient, burstiness

# Generate some random data from an exponential distirbution
n_sim = 1000
n_events = 8
# We run nsim simulations, each with n_events
# Do efficiently by doing all at once and then reshaping,
# as simulations are independently distributed
ie_times = expon(scale=100).rvs(size=(n_sim*n_events))
print(ie_times)
ie_times = np.reshape(ie_times, (n_sim, n_events))
print(ie_times)
plt.hist(ie_times, bins=20)
plt.savefig('Exponential_hist.png')
ie_times_T = ie_times.T
mem = memory_coefficient(ie_times_T)
burst = burstiness(ie_times_T)
print(mem)
print(burst)
print(min(mem))
print(max(mem))
print(min(burst))
print(max(burst))

# Now plot in M-B space
plt.clf()
ax = plt.subplot(111) 
plt.scatter(mem, burst)
ax.set_xlim([-1, 1])
ax.set_ylim([-1, 1])
ax.set_xlabel('M')
ax.set_ylabel('B')
# Add y = 0, x=0 lines
plt.plot([0,0],[-1, 1], linestyle='dashed', linewidth=1, c='0.5')
plt.plot([-1,1],[0, 0], linestyle='dashed', linewidth=1, c='0.5') 
plt.savefig('Exponential_B_M_diagram.png')

# Now try gamma distribution
ie_times = gamma(2.0, scale=100).rvs(size=(n_sim*n_events))
print(ie_times)
ie_times = np.reshape(ie_times, (n_sim, n_events))
plt.hist(ie_times, bins=20)
plt.savefig('Gamma_hist.png') 
ie_times_T = ie_times.T
mem = memory_coefficient(ie_times_T)
burst = burstiness(ie_times_T)
# Now plot in M-B space
plt.clf()
ax = plt.subplot(111)
plt.scatter(mem, burst)
ax.set_xlim([-1, 1])
ax.set_ylim([-1, 1])
ax.set_xlabel('M')
ax.set_ylabel('B')
plt.plot([0,0],[-1, 1], linestyle='dashed', linewidth=1, c='0.5')
plt.plot([-1,1],[0, 0], linestyle='dashed', linewidth=1, c='0.5')
plt.savefig('Gamma_B_M_diagram.png')

# Now try Weibull distribution
ie_times = weibull_min(2.0, scale=100).rvs(size=(n_sim*n_events))
print(ie_times)
ie_times = np.reshape(ie_times, (n_sim, n_events))
plt.hist(ie_times, bins=20)
plt.savefig('Weibull_hist.png') 
ie_times_T = ie_times.T
mem = memory_coefficient(ie_times_T)
burst = burstiness(ie_times_T)
# Now plot in M-B space
plt.clf()
ax = plt.subplot(111)
plt.scatter(mem, burst)
ax.set_xlim([-1, 1])
ax.set_ylim([-1, 1])
ax.set_xlabel('M')
ax.set_ylabel('B')
plt.plot([0,0],[-1, 1], linestyle='dashed', linewidth=1, c='0.5')
plt.plot([-1,1],[0, 0], linestyle='dashed', linewidth=1, c='0.5')
plt.savefig('Weibull_B_M_diagram.png')
