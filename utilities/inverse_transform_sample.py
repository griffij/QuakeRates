"""Generate random samples using inverse transform sampling
so that we can efficiently sample conditional probabilities.

Wikipedia summary of method suffices:
https://en.wikipedia.org/wiki/Inverse_transform_sampling

Jonathan Griffin
University of Otago/Geoscience Australia
May 2020
"""

import os, sys
import numpy as np
from numpy import inf
from scipy.stats import expon, uniform

def ivt_expon(lam, a=0, b=inf, n_samples=1):
    """Generate random samples from an exponential function defined
    by rate lambda and between points a and b.
    """
    a_update = expon.cdf(a, scale=1/lam) # Convert to uniform distribution space
    b_update = expon.cdf(b, scale=1/lam) # and here
    # Get uniform distribution over [a, b] in transformed space
    rv_unif = uniform.rvs(loc=a_update, scale=(b_update-a_update), size=n_samples)
    rv_exp = (-1/lam)*np.log(1-rv_unif)
    return rv_exp


if __name__ == "__main__":
    mean_int = 1000 # Mean interval; lambda = 1/mean_int
    n_samples = 10
    a=10000
    b=inf
    ivt_exp = ivt_expon(1/mean_int, a=a, b=b, n_samples = n_samples)
    print(ivt_exp)
