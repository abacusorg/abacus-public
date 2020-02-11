#Misc functions that don't belong anywhere else

import ctypes as ct
import os
import csv
import contextlib
import re
from glob import glob
import string
from collections import OrderedDict

import numpy as np
import numexpr as ne

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
    # We always want numpy arrays we pass to C to be contiguous
    if not a.flags.contiguous:
        raise RuntimeError('Trying to pass non-contiguous array to C function! ' + str(a))
    #a = np.ascontiguousarray(a)  This may cause a copy!
    if a.dtype == np.float32:
        return a.ctypes.data_as(ct.POINTER(ct.c_float))
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
def add_tick(ax, loc, label, ha=None):
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

    if ha:
        labels = ax.get_xticklabels()
        labels[-1].set_ha(ha)
    
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

    TODO: rewrite
    '''
    
    def __init__(self, *args, cumulative=False, **kwargs):
        self.cumulative = cumulative
        self.cumulative_time = 0.  # all the time before the most recent/current loop
        args = list(args)
        if len(args) > 0:
            kwargs['prefix'] = args.pop(0)

        kwargs['fmt'] = kwargs.get('fmt', 'time: {:.2f} seconds')
        
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

    def start(self):
        self.__enter__()

    def stop(self, report=True):
        _output = self.output
        self.output=report
        self.__exit__(None,None,None)
        self.output=_output

    def report(self):
        print(" ".join([self.prefix, self.fmt.format(self.elapsed)]))

def scatter_density(x, y, ax, z=None, size=10., log=False, bw=.03, adaptive=False, **scatter_kwargs):
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

    #from matplotlib.mlab import griddata
    from sklearn.neighbors.kde import KernelDensity
    #from scipy.stats import gaussian_kde
    #from KDEpy.TreeKDE import TreeKDE

    if z is not None:
        xy = np.vstack([x,y,z]).T
    else:
        xy = np.vstack([x,y]).T

    if log:
        pos = (xy > 0).all(axis=-1)
        xy = np.log10(xy[pos])
        x = np.log10(x[pos])
        y = np.log10(y[pos])
    
    # Method 1: scipy
    #color = gaussian_kde(xy)(xy)
    
    # Method 2: sklearn
    kde = KernelDensity(kernel='gaussian', bandwidth=bw, rtol=1e-2).fit(xy)
    color = kde.score_samples(xy)

    # Method 3: KDEpy
    #kde = TreeKDE(bw=bw).fit(xy)
    #color = kde.evaluate(xy)
    #color = np.log(color)

    if adaptive:
        # Shrink the bw in high-density regions; grow it in low density regions
        # TODO: very little effect.  Might need to shrink in log space.

        # scale to 0 to 1
        color -= color.min()
        color /= color.max()
        #import pdb; pdb.set_trace()
        alpha = 0.99999
        newbw = bw*(1 - alpha*color)

        print('Starting new KDE...')
        kde = TreeKDE(bw=newbw).fit(xy)
        color = kde.evaluate(xy)
        color = np.log(color)

    
    idx = color.argsort()
    x, y, color = x[idx], y[idx], color[idx]
    
    if log:
        sc = ax.scatter(10**x, 10**y, s=size, c=color, **scatter_kwargs)
    else:
        sc = ax.scatter(x, y, s=size, c=color, **scatter_kwargs)
        
    return x, y, color, kde, sc

import argparse
# Combine argparse mixins to format both the description and defaults
# Use as:
# >>> parser = argparse.ArgumentParser(description='...', formatter_class=Tools.ArgParseFormatter)
class ArgParseFormatter(argparse.RawDescriptionHelpFormatter, argparse.ArgumentDefaultsHelpFormatter):
    pass


import numba
@numba.njit(parallel=True)
def wrap_zero_centered(pos, box):
    '''
    Wraps an array of positions in-place
    to the range [-box/2, box/2).
    
    Parameters
    ----------
    pos: ndarray
        The positions to wrap, already in the same units as box
    box: float
        The box edge length (box is zero-centered)
    '''
    pos = pos.reshape(-1)
    N = len(pos)
    for i in numba.prange(N):
        while pos[i] >= box/2:
            pos[i] -= box
        while pos[i] < -box/2:
            pos[i] += box

@numba.njit(parallel=True)
def wrap_zero_origin(pos, box):
    '''
    Wraps an array of positions in-place
    to the range [0, box).
    
    Parameters
    ----------
    pos: ndarray
        The positions to wrap, already in the same units as box
    box: float
        The box edge length (box has center at L/2)
    '''
    pos = pos.reshape(-1)
    N = len(pos)
    for i in numba.prange(N):
        while pos[i] >= box:
            pos[i] -= box
        while pos[i] < 0:
            pos[i] += box
    
if __name__ == '__main__':
    fm = get_RAM_info()
    print(fm) 
    
