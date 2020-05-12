"""Plot example datasets in memory-burstiness space
"""

import os, sys
from matplotlib import pyplot as plt
from matplotlib.patches import PathPatch
from scipy.stats import kde
import numpy as np
from scipy.stats import expon, gamma, weibull_min  
from QuakeRates.utilities.memory_coefficient import memory_coefficient, burstiness

# Generate some random data from an exponential distirbution
n_sim = 100000
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


# Evaluate a gaussian kde on a regular grid of nbins x nbins over
# data extents
nbins=100
data_ar = np.array([mem, burst])
x, y = data_ar
k = kde.gaussian_kde(data_ar)
xi, yi = np.mgrid[0.98*x.min():1.02*x.max():nbins*1j, 0.98*y.min():1.02*y.max():nbins*1j]
zi = k(np.vstack([xi.flatten(), yi.flatten()]))
# Calculate percentiles
cumsum = np.cumsum(zi)
# Normalise data and calculate bottom 5th percentile for drawing contour
# that contains 95% of the distribution.
zi_norm = zi/max(cumsum)
perc_5th = 0.05*max(zi_norm)
# Plot the data
# Slightly extend the bounds of the data for smoother plotting
cs = ax.contour(xi, yi, zi_norm.reshape(xi.shape), [perc_5th], colors='k',
                linewidths=0.4)
for c in cs.collections:
    try:
        c_path = c.get_paths()[0]
        # Dump out vertices for plotting later
        filename = 'Exponential_B_M_95per_contour_nsim_%i_nevents_%i.txt' % (n_sim, n_events)
        np.savetxt(filename, c_path.vertices, delimiter=',')
        patch = PathPatch(c_path, transform=ax.transData, facecolor='none',
                          linewidth=0.4, linestyle='--')
    except:
        print('Cannot plot this contour, does not exist')
m = ax.pcolormesh(xi, yi, zi_norm.reshape(xi.shape), shading='gouraud',
                  cmap=plt.cm.Greys, clip_path=patch)


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

# Evaluate a gaussian kde on a regular grid of nbins x nbins over
# data extents
nbins=100
data_ar = np.array([mem, burst])
x, y = data_ar
k = kde.gaussian_kde(data_ar)
xi, yi = np.mgrid[0.98*x.min():1.02*x.max():nbins*1j, 0.98*y.min():1.02*y.max():nbins*1j]
zi = k(np.vstack([xi.flatten(), yi.flatten()]))
# Calculate percentiles
cumsum = np.cumsum(zi)
# Normalise data and calculate bottom 5th percentile for drawing contour
# that contains 95% of the distribution.
zi_norm = zi/max(cumsum)
perc_5th = 0.05*max(zi_norm)
# Plot the data
# Slightly extend the bounds of the data for smoother plotting
cs = ax.contour(xi, yi, zi_norm.reshape(xi.shape), [perc_5th], colors='k',
                linewidths=0.4)
for c in cs.collections:
    try:
        c_path = c.get_paths()[0]
        filename = 'Gamma_B_M_95per_contour_nsim_%i_nevents_%i.txt' % (n_sim, n_events)
        np.savetxt(filename, c_path.vertices, delimiter=',')  
        patch = PathPatch(c_path, transform=ax.transData, facecolor='none',
                          linewidth=0.4, linestyle='--')
    except:
        print('Cannot plot this contour, does not exist')
m = ax.pcolormesh(xi, yi, zi_norm.reshape(xi.shape), shading='gouraud',
                  cmap=plt.cm.Greys, clip_path=patch)

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

# Evaluate a gaussian kde on a regular grid of nbins x nbins over
# data extents
nbins=100
data_ar = np.array([mem, burst])
x, y = data_ar
k = kde.gaussian_kde(data_ar)
xi, yi = np.mgrid[0.98*x.min():1.02*x.max():nbins*1j, 0.98*y.min():1.02*y.max():nbins*1j]
zi = k(np.vstack([xi.flatten(), yi.flatten()]))
# Calculate percentiles
cumsum = np.cumsum(zi)
# Normalise data and calculate bottom 5th percentile for drawing contour
# that contains 95% of the distribution.
zi_norm = zi/max(cumsum)
perc_5th = 0.05*max(zi_norm)
# Plot the data
# Slightly extend the bounds of the data for smoother plotting
cs = ax.contour(xi, yi, zi_norm.reshape(xi.shape), [perc_5th], colors='k',
                linewidths=0.4)
for c in cs.collections:
    try:
        c_path = c.get_paths()[0]
        filename = 'Weibull_B_M_95per_contour_nsim_%i_nevents_%i.txt' % (n_sim, n_events)
        np.savetxt(filename, c_path.vertices, delimiter=',')  
        patch = PathPatch(c_path, transform=ax.transData, facecolor='none',
                          linewidth=0.4, linestyle='--')
    except:
        print('Cannot plot this contour, does not exist')
m = ax.pcolormesh(xi, yi, zi_norm.reshape(xi.shape), shading='gouraud',
                  cmap=plt.cm.Greys, clip_path=patch)


plt.savefig('Weibull_B_M_diagram.png')
