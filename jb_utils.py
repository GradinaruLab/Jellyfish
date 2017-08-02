"""
Set of utilities written by Justin Bois for use in BE/Bi 103 (2014
edition) and beyond.
"""
import os
import glob
import warnings
import matplotlib.path as path
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.cm as cm
import numpy as np

try:
    import numdifftools as nd
except:
    warnings.warn('Unable to import numdifftools.  hess_nd unavailable.',
                  ImportWarning)

try: 
    import pywt
except:
    warnings.warn('Unable to import PyWavelets. visushrink will not work.',
                  ImportWarning)

try:
    import skimage.io
except:
    warnings.warn('Unable to import skimage.  ' \
                  + 'Image processing utils will not work.', ImportWarning)
    
try:
    import Image
except:
    warnings.warn('Unable to import PIL.  ' \
                  + 'Image processing utils will not work.', ImportWarning)

try:
    import skimage.io
    import skimage.measure
except:
    warnings.warn('Unable to import skimage. '\
                  + 'Image processing utils will not work.', ImportWarning)
        
    
# ###############################################
# FOLLOWING ARE UTILITIES DATA SMOOTHING
# ###############################################
# ############################
def epan_kernel(t):
    """
    Epanechnikov kernel.
    """
    return np.logical_and(t > -1.0, t < 1.0) * 3.0 * (1.0 - t**2) / 4.0

# ############################
def tri_cube_kernel(t):
    """
    Tri-cube kernel.
    """
    return np.logical_and(t > -1.0, t < 1.0) * (1.0 - abs(t**3))**3

# ############################
def gauss_kernel(t):
    """
    Gaussian kernel.
    """
    return np.exp(-t**2 / 2.0)

# ############################
def nw_kernel_smooth(x_0, x, y, kernel_fun, lam):
    """
    Gives smoothed data at points x_0 using a Nadaraya-Watson kernel 
    estimator.  The data points are given by NumPy arrays x, y.
        
    kernel_fun must be of the form
        kernel_fun(t), 
    where t = |x - x_0| / lam
    
    This is not a fast way to do it, but it simply implemented!
    """
    
    # Function to give estimate of smoothed curve at single point.
    def single_point_estimate(x_0_single):
        """
        Estimate at a single point x_0_single.
        """
        t = np.abs(x_0_single - x) / lam
        return np.dot(kernel_fun(t), y) / kernel_fun(t).sum()
    
    # If we only want an estimate at a single data point
    if np.isscalar(x_0):
        return single_point_estimate(x_0)
    else:  # Get estimate at all points
        y_smooth = np.empty_like(x_0)
        for i in range(len(x_0)):
            y_smooth[i] = single_point_estimate(x_0[i])
        return y_smooth

