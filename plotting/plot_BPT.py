"""Plot some basic BPT distributions
"""
import os, sys
from matplotlib import pyplot as plt
import numpy as np
from scipy.stats import expon, gamma, weibull_min, invgauss

mu = 100
alphas = [ 0.5, 1., 2., 5., 10.]
#alpha = 1
x_vals = np.arange(0, 4*mu)
# Plot for a range of alpha values
for alpha in alphas:
    bpt = invgauss(alpha, scale=mu)
    pdf_vals = bpt.pdf(x_vals)
    plt.plot(x_vals, pdf_vals, label=alpha)

# Now add exponential
exp_dist = expon(scale = mu)
pdf_vals = exp_dist.pdf(x_vals)
plt.plot(x_vals, pdf_vals, label='Exponential', color='k')

plt.legend()
plt.savefig('BPT_distribution.png')
