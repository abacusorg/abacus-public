#Misc functions that don't belong anywhere else

import ctypes as ct
import os
import csv
import numpy as np
import contextlib
import re
from glob import glob

import string
import numexpr as ne
from collections import OrderedDict
from .Reglob import reglob

# all values returned in GB
# try to return 'available' RAM, but some systems only have 'total' and 'free', so fallback to those
def get_RAM_info(field=None):
    auto = field is None
    field_names = OrderedDict([('available', 'MemAvailable'), ('total', 'MemTotal'), ('free', 'MemFree')])
    if not auto:
        field = field.lower()
        field_names = {field: field_names[field]}
    
    meminfo = list(csv.reader(open("/proc/meminfo","r"), delimiter = ":"))
    
    for field in field_names:
        for line in meminfo:
            if line[0] == field_names[field]:
                mem = line[1].split()[0]
                return float(mem)/1024**2.
    
    raise RuntimeError('No field(s) named {} in /proc/meminfo'.format(list(field_names.keys())))
    

def getRuntime(timingfile,timetype = "Total Wall Clock Time"):
    timeinfo = csv.reader(open(timingfile,"r"), delimiter = ":")
    for line in timeinfo:
        if line:
            if timetype in line[0]:
                s = line[1].find("s")
                t = 0
                t = eval((line[1])[0:s])
                return t

def cpointer(a):
    if a is None:
        return a  # Sometimes we deliberately want to pass a NULL pointer.  `None` is the way to do that.
    # We always want numpy arrays we pass to C to be contiguous.  This may cause a copy!
    if not a.flags.contiguous:
        raise RuntimeError('Trying to pass non-contiguous array to C function! ' + str(a))
    #a = np.ascontiguousarray(a)
    if a.dtype == np.float32:
        return  a.ctypes.data_as(ct.POINTER(ct.c_float))
    elif a.dtype == np.float64:
        return a.ctypes.data_as(ct.POINTER(ct.c_double))
    elif a.dtype == np.int64:
        return a.ctypes.data_as(ct.POINTER(ct.c_longlong))
    elif a.dtype == np.uint64:
        return a.ctypes.data_as(ct.POINTER(ct.c_ulonglong))
    elif a.dtype == np.int32:
        return a.ctypes.data_as(ct.POINTER(ct.c_int))
    elif a.dtype == np.uint32:
        return a.ctypes.data_as(ct.POINTER(ct.c_uint))
    # Complex types will be passed as float/double arrays
    elif a.dtype == np.complex64:
        return a.ctypes.data_as(ct.POINTER(ct.c_float))
    elif a.dtype == np.complex128:
        return a.ctypes.data_as(ct.POINTER(ct.c_double))
    elif np.issubdtype(a.dtype, np.void):
        return a.ctypes.data_as(ct.c_void_p)
    else:
        raise TypeError(a.dtype)
        
# Use this as a ctypes 'argtype' to transparently pass numpy arrays as C pointers
class ndarray_arg:
    @classmethod
    def from_param(cls, ndarray):
        return cpointer(ndarray)

class asciistring_arg:
    @classmethod
    def from_param(cls, s):
        return ct.c_char_p(s.encode('ascii'))
        
# This is the greatest use of a context manager ever.
@contextlib.contextmanager
def chdir(dirname=None):
  curdir = os.getcwd()
  try:
    if dirname is not None:
      os.chdir(dirname)
    yield
  finally:
    os.chdir(curdir)

# Adds a tick to the given axis with the given value and name
def add_tick(ax, loc, label):
    # Generate ticks and lims
    ax.figure.canvas.draw()
    lims = ax.get_xlim()

    # Modify ticks (this messes up lims... ugh)
    locs = ax.get_xticks()
    labels = ax.get_xticklabels()
    locs = list(locs) + [loc]
    labels = [l.get_text() for l in list(labels)] + [label]
    ax.set_xticks(locs)
    ax.set_xticklabels(labels)
    ax.set_xlim(lims)
    
import matplotlib.pyplot as plt
def matrix_plot(m, fig=None, ax=None, contour=False, subplots_kwargs={}, contour_kwargs={}, imshow_kwargs={}, **kwargs):
    if not fig or not ax:
        fig, ax = plt.subplots(**subplots_kwargs)
    ax.set_xlabel(kwargs.pop('xlabel', ''))
    ax.set_ylabel(kwargs.pop('ylabel', ''))
    ax.set_title(kwargs.pop('title', ''))
    
    if 'extent' in imshow_kwargs:
        xmin,xmax,ymin,ymax = imshow_kwargs['extent']
    else:
        xmin,xmax,ymin,ymax = 0., m.shape[0], 0., m.shape[1]
    x = np.linspace(xmin,xmax,m.shape[0])
    y = np.linspace(ymin,ymax,m.shape[1])
    X, Y = np.meshgrid(y, x)
    
    cax = ax.imshow(m, origin='lower', interpolation='none', **imshow_kwargs)
    fig.colorbar(cax)
    
    if contour:
        colors = contour_kwargs.pop('colors', 'k')
        cs = ax.contour(X, Y, m, colors=colors, **contour_kwargs)
        ax.clabel(cs, fmt='%.2g')
        
    return fig, ax
    
import matplotlib.colors as colors
class MidpointNormalize(colors.Normalize):
    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        colors.Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        # I'm ignoring masked values and all kinds of edge cases to make a
        # simple example...
        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y))

def numexpr_dist(arrays, dist=True):
    """
    To sum an unknown number of arrays in numexpr, you have to build
    a string with all the arrays explicitly given variable names.
    This function does that and evaluates it as a numexpr expression.
    
    This is equivalent to:
    
    dist = sum(a**2. for a in arrays)**.5
    
    Parameters
    ----------
    arrays: iterable of ndarray
        A list of arrays to sum
    dist: bool, optional
        Instead of the plain sum, return the square root of the sum
        of the square of the elements.
        
    Returns
    -------
    dist: ndarray
        The result of the sum
    
    """
    # assign each array a letter
    local_dict = {l:c for l,c in zip(string.ascii_letters,arrays)}
    assert len(local_dict) == len(arrays)
    
    if dist:
        nedist = ne.evaluate('({})**.5'.format('+'.join(s+'**2' for s in local_dict)), local_dict=local_dict)
    else:
        nedist = ne.evaluate('+'.join(local_dict), local_dict=local_dict)
        
    return nedist

def set_global_nthreads(n):
    '''
    Sets environment variables that control the number
    of threads for a variety of applications.
    
    This must be run before numba's first import to have
    any effect on numba.
    '''
    import os
    import sys
    os.environ['OMP_NUM_THREADS'] = str(n)
    if 'numba' in sys.modules:
        import numba
        if numba.config.NUMBA_NUM_THREADS != n:
            print('''Warning: numba has already been imported with {} threads;
                     setting the number of threads will have no effect.'''.format(numba.config.NUMBA_NUM_THREADS))
    os.environ['NUMBA_NUM_THREADS'] = str(n)
    
    import numexpr as ne
    ne.set_num_threads(n)
    ne.set_vml_num_threads(1)  # use numexpr-level parallelism instead
    
    from Abacus.Analysis.PowerSpectrum import PowerSpectrum
    PowerSpectrum.nthreads = n
    #pslib.set_nthreads(n)
    
import contexttimer
class ContextTimer(contexttimer.Timer):
    '''
    A simple extension to the contexttimer lib that adds
    the `cumulative` keyword arg.
    '''
    
    def __init__(self, prefix, *args, cumulative=False, **kwargs):
        self.cumulative = cumulative
        self.cumulative_time = 0.  # all the time before the most recent/current loop
        kwargs['prefix'] = prefix
        super(ContextTimer, self).__init__(*args, **kwargs)
        
    def __enter__(self):
        if self.cumulative and self.end:
            # The time from the most recent loop gets added to the cumulative time
            self.cumulative_time += self.end - self.start
        return super(ContextTimer, self).__enter__()
        
    @property
    def elapsed(self):
        last_time = super(ContextTimer, self).elapsed
        return self.cumulative_time*self.factor + last_time

from matplotlib.mlab import griddata
from sklearn.neighbors.kde import KernelDensity
from scipy.stats import gaussian_kde
def scatter_density(x, y, ax, z=None, size=10., log=False, bw=.03):
    '''
    A common problem in scatter plots is overlapped points.
    One can use 2D histogramming to plot densities instead, but
    one may still be interested in outlier points which
    will not show up very well in histograms.

    An alternative is to still plot points intead of a
    histogram, but color the points by a kernel density
    estimate of the points.  By plotting the "densest" points
    on top, one can get a sense of the density while still
    preserving outliers.

    Parameters
    ----------
    x, y: array_like
        The points to plot
    z: array_like, optional
        If given, does the KDE in 3D.  Only x,y are used for
        plotting though.  In theory, one could use this for
        the zorder as well, but that probably wouldn't look good.
    ax: matplotlib.Axes
        The axes to plot on
    size: float, optional
        The plotting size of the points
    log: bool, optional
        Whether to do the KDE in log-log space
    bw: float, optional
        The KDE bandwidth, or "smoothness", in data units.
        This hugely impacts the runtime. If plotting is taking
        a while, try reducing this.
    '''
    if z is not None:
        xy = np.vstack([x,y,z]).T
    else:
        xy = np.vstack([x,y]).T

    if log:
        pos = (xy > 0).all(axis=-1)
        xy = np.log10(xy[pos])
        x = np.log10(x[pos])
        y = np.log10(y[pos])
    
    #color = gaussian_kde(xy)(xy)
    kde = KernelDensity(kernel='gaussian', bandwidth=bw, rtol=1e-2).fit(xy)
    color = kde.score_samples(xy)
    
    idx = color.argsort()
    x, y, color = x[idx], y[idx], color[idx]
    
    if log:
        sc = ax.scatter(10**x, 10**y, s=size, c=color)
    else:
        sc = ax.scatter(x, y, s=size, c=color)
        
    return x, y, color, kde, sc

import argparse
# Combine argparse mixins to format both the description and defaults
class ArgParseFormatter(argparse.RawDescriptionHelpFormatter, argparse.ArgumentDefaultsHelpFormatter):
    pass
    
if __name__ == '__main__':
    fm = get_RAM_info()
    print(fm) 
    
