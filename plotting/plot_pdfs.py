# Plot pdfs of different distirbutions

import os, sys
from matplotlib import pyplot as plt
import numpy as np
from scipy.stats import expon, gamma, weibull_min, lognorm


mean = 500
exp = expon(scale=mean)

gam = gamma(2.2, scale=500/2.2)

x = np.arange(0.01, 3000, 1)
exp_pdf = exp.pdf(x) 
plt.plot(x, exp_pdf, 'k-', lw=2, label=r'Exponential $(\theta = 500)$')
plt.plot(x, gam.pdf(x), 'b-', lw=2, label=r'Gamma $(\alpha = 2.2, \theta = 500/\alpha)$')
plt.legend()
plt.xlabel('Inter-event time (Years)')
plt.ylabel('Probability density')
# Add mean
plt.xlim([0, 3000])
plt.ylim([0, max(exp_pdf)]) 
plt.plot([mean, mean], [0, max(exp_pdf)], c='0.5', linestyle='--')
plt.tight_layout()
plt.savefig('pdf_examples.png')
