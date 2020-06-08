"""Functions for bi-linear fitting
With thanks to Trevor Allen, Geoscience Australia
"""

import numpy as np
from numpy import zeros_like
from scipy.odr import Data, Model, ODR, models
import scipy.odr.odrpack as odrpack

# functions to return binary arrays
def highside(x, hx):
    from numpy import zeros_like
    xmod = zeros_like(x)
    
    idx = x >= hx
    xmod[idx] = 1
    return xmod
    
def lowside(x, hx):
    from numpy import zeros_like
    xmod = zeros_like(x)
    
    idx = x <= hx
    xmod[idx] = 1
    return xmod

def bilinear_reg_zero_slope(c, x):
    hx = c[2] # x-hinge
    ans2 = zeros_like(x)
    ans1 = zeros_like(x)
    
    idx1 = x <= hx
    idx2 = x >= hx
    
    modx_lo = lowside(x, hx)
    modx_hi = highside(x, hx)
    
    ans1 = modx_lo * (c[0] * x + c[1])
    yarea = c[0] * hx + c[1]
    ans2 = modx_hi * yarea
    
    return ans1 + ans2  

def bilinear_reg_fix(c, x):
    from numpy import zeros_like
    hxfix = np.log10(2e-4) #4.0 # hinge magnitude
    ans2 = zeros_like(x)
    ans1 = zeros_like(x)

    #idx1 = x <= hx
    #idx2 = x >= hx

    modx_lo = lowside(x, hxfix)
    modx_hi = highside(x, hxfix)

    ans1 = modx_lo * (c[0] * x + c[1])
    yarea = c[0] * hxfix + c[1]
    ans2 = modx_hi * (c[2] * (x-hxfix) + yarea)

    return ans1 + ans2

def bilinear_reg_fix_zero_slope(c, x):
    from numpy import zeros_like
    hxfix = np.log10(1.4e-4) #4.0 # hinge magnitude
    ans2 = zeros_like(x)
    ans1 = zeros_like(x)

    #idx1 = x <= hx
    #idx2 = x >= hx

    modx_lo = lowside(x, hxfix)
    modx_hi = highside(x, hxfix)

    ans1 = modx_lo * (c[0] * x + c[1])
    yarea = c[0] * hxfix + c[1]
#    yarea = c[0] * hx + c[1] 
#    ans2 = modx_hi * (c[2] * (x-hxfix) + yarea)
    ans2 = modx_hi * yarea 
    return ans1 + ans2
