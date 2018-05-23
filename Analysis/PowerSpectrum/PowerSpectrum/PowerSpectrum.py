"""
PowerSpectrum.py

The library interface to compute power spectra with a variety of options
from a variety of file formats.

For a user-friendly command-line interface, see `run_PS.py`.

     Author: Lehman Garrison
     based on dferrer's original 2013 C implementation
"""

import sys
import ctypes as ct
import numpy as np
import scipy
import scipy.interpolate
import os
import glob
import pyfftw
import numexpr as ne


from Abacus.Tools import ndarray_arg
import Abacus.ReadAbacus
from . import misc
from .misc import processargs
from .Histogram import RadialBin, RadialBinGrid, k_bin_edges
from . import TSC

# nthreads and verbosity are set globally
# TODO: add these to processargs
import multiprocessing
nthreads = multiprocessing.cpu_count()

pi = np.pi

verbose = True

# Control the number of OpenMP threads that the underlying C library will use
#import pslib
#from pslib import get_lib
#set_nthreads = pslib.set_nthreads

@processargs
def CalculateFromMem(positions, gridshape, boxsize, *, weights=None, dtype=np.float32, rotate_to=None, bins=-1, log=False, window='window_aliased'):
    """
    Computes the power spectrum of a set of points by binning them into a density field.
    
    Parameters
    ----------
    positions: ndarray of shape (N,3), or list of ndarray
        The points to bin
    gridshape: int or tuple of int
        The FFT mesh size.  Assumes 3D cube if an int.
        Should not include padding for an rfft.
    boxsize: float
        The physical domain size.
        `positions` are assumed to be in a zero-centered box of side length `boxsize`.
    weights: ndarray of shape (N,), or list of ndarray, optional
        Binning weights for the points.
    dtype: np.dtype, optional
        The data type to do the computations in.
    rotate_to: array of floats of shape (3,), optional
        `rotate_to` defines a vector that the cartesian z-hat will be rotated to align with.
        All rotations are periodic "safe", meaning we discard particles at the edges.
        Default: None.
    bins: int or array_like of float, optional
        Number of k-space bins, or list of bin edges.  Default: max(gridshape)//4.
    log: bool, optional
        Whether to bin in log space.  Default: False.
    window: str, optional
        Divide out the TSC window function.
        Options are: None, 'window' and 'window_aliased' (default).
        'window' only divides out the TSC window function.
        'window_aliased' also compensates for aliasing across
        Nyquist, but assumes constant P(k) past Nyquist.
    
    Returns
    -------
    k: ndarray
        Average wavenumber per k-bin
    s: ndarray
        Average power per k-bin
    nb: ndarray
        Total modes per k-bin
    """
    
    # Get the density field
    fielddensity = TSC.BinParticlesFromMem(positions, gridshape, boxsize, weights=weights, dtype=dtype, rotate_to=rotate_to, prep_rfft=True, nthreads=nthreads)
    
    # Do the FFT
    res = FFTAndBin(fielddensity, boxsize, bins=bins, log=log, inplace=True, window=window)
    
    return res


@processargs
def CalculateBySlab(file_pattern, gridshape, boxsize, *, dtype=np.float32, zspace=False, format='pack14', rotate_to=None, bins=-1, log=False, window='window_aliased'):
    '''
    Main entry point for PS computation from files on disk of various formats.
    
    The reading and binning is done slab-by-slab, so we don't need to fit the whole particle set
    in memory.

    Parameters:
    -----------
    file_pattern: str
        A globbing pattern for all the files to read
    boxsize: float
        Size of the box.  Particles (except for RVzel and Gadget format) are assumed to be stored in unit box.
    gridshape: int or tuple of int
        The FFT mesh size.  Assumes 3D cube if an int.
        Should not include padding for an rfft.
    dtype: np.dtype, optional
        Precision for the density array.  Default: np.float32.
    zspace: bool, optional
        Whether to bin the particles in redshift space.  Default: False.
    format: str, optional
        File format for the input files.  Options are 'rvdouble', 'pack14', 'rvzel', 'state', 'gadget'.
        Default: 'pack14'
    rotate_to: array of floats of shape (3,), optional
        `rotate_to` defines a vector that the cartesian z-hat will be rotated to align with.
        All rotations are periodic "safe", meaning we discard particles at the edges.
        Default: None.
    bins: int or array_like of float, optional
        Number of k-space bins, or list of bin edges.  Default: max(gridshape)//4.
    log: bool, optional
        Whether to bin in log space.  Default: False.
    window: str, optional
        Divide out the TSC window function.
        Options are: None, 'window' and 'window_aliased' (default).
        'window' only divides out the TSC window function.
        'window_aliased' also compensates for aliasing across
        Nyquist, but assumes constant P(k) past Nyquist.


    Returns:
    --------
    k: ndarray, shape (nbins,)
        The wavenumbers of the bins
    P: ndarray, shape (nbins,)
        Power in each bin
    nmodes: ndarray, shape (nbins,)
        Number of modes in each bin
        
    If `bins` is None, instead the following is returned:
    power: ndarray, real, shape (gridsize, [gridsize,] gridsize//2 + 1)
        The full 2D/3D power spectrum (Fourier transform of the deltas, squared)
    '''
    
    # Do the binning
    density = TSC.BinParticlesFromFile(file_pattern, boxsize, gridshape, dtype=dtype, zspace=zspace, format=format, rotate_to=rotate_to, prep_rfft=True, nthreads=nthreads)
    
    # Do the FFT
    return FFTAndBin(density, boxsize, bins=bins, log=log, inplace=True, window=window)
    

@processargs
def RebinTheoryPS(ps, gridshape, boxsize, *, bins=-1, log=False, dtype=np.float32):
    """
    If we want to compare a theory power spectrum to one generated by binning
    the density field (as this code does), we need to evaluate the theory PS on
    the same FFT mesh as the density field.  This function does that then returns
    the binned result.
    
    The interpolation to the mesh is done with cubic spline in log space.
    This is the same method our zeldovich code uses, so the results should
    match nearly to machine precision.
    
    Parameters
    ----------
    ps: ndarray, shape (N,2)
        The theory power spectrum to interpolate to the mesh
    boxsize: float
        The domain size
    gridshape: int or tuple of int
        The FFT mesh size.  Assumes 3D cube if an int.
        Should not include padding for an rfft.
    bins: int or array_like of float, optional
        Number of k-space bins, or list of bin edges.  Default: max(gridshape)//4.
    log: bool, optional
        Whether to bin in log space.  Default: False.
        
    Returns
    -------
    k: ndarray, shape (nbins,)
        The wavenumbers of the bins
    s: ndarray, shape (nbins,)
        Power in each bin
    nmodes: ndarray, shape (nbins,)
        Number of modes in each bin
    """
    pkref_interp = misc.loglog_interp(ps[:,0], ps[:,1])
    
    dx = boxsize/np.array(gridshape)
    kgrid = misc.fftfreq_nd(gridshape, dx)
    power = pkref_interp(kgrid)
    
    k, s, nmodes, bin_info = RadialBinGrid(boxsize, power, bins, rfft=True)

    return k, s, nmodes


#@processargs  # need to thin about args for this function
def StackedPower(*args, **kwargs):
    """
    Estimate a power spectrum by averaging many realizations.
    Arguments are the same as for `CalculateFromMem`, but `positions`
    (and optionally `weights`) should be list(s).  Each element in the
    list is interpreted as being from an individual realization of the
    power spectrum.
    
    Note that there's no good way to evaluate the stacked power
    spectrum by stacking the density fields.  You really have to
    do all the FFTs and stack the result.
    
    The stacking is always done on the full 2D/3D power spectrum.
    Any binning is done on the stacked result.
    
    Docstring for `CalculateFromMem` follows:
    """
    gridshape = args[1]
    boxsize = args[2]
    dtype = kwargs.get('dtype', np.float32)
    log = kwargs.get('log', False)
    
    # Remove positions and weights from the args
    # We will only be passing one element at a time from them
    positions = args[0]; args = args[1:]
    weights = kwargs.pop('weights', None)
    bins = kwargs.pop('bins', -1)
    
    assert type(positions) == list
    
    if weights is None:
        weights = [None,]*len(positions)
    else:
        assert type(weights) == list
    
    ps_shape = np.array(gridshape)
    ps_shape[-1] = ps_shape[-1]//2 + 1

    mean_ps = pyfftw.zeros_aligned(ps_shape, dtype=dtype)
    ncat = len(positions)

    for pos,w in zip(positions,weights):
        ps = CalculateFromMem(pos, *args, weights=w, nbins=None, **kwargs)
        ne.evaluate('mean_ps + ps/ncat', out=mean_ps)
    del ps
        
    if bins != None:
        bin_edges = k_bin_edges(gridshape, boxsize, nbins=bins, log=log)
        dx = boxsize/np.array(ps_shape)
        k, s, nmodes, bin_info = RadialBinGrid(boxsize, mean_ps, bin_edges, rfft=True)
        # Effective number of modes from stacking
        nmodes *= ncat
        return k, s, nmodes
    
    return mean_ps
    
StackedPower.__doc__ += CalculateFromMem.__doc__
    

@processargs
def CrossPower(gridshape, boxsize, *, pos1=None, pos2=None, delta1=None, delta2=None, dtype=np.float32, bins=-1, log=False, weights1=None, weights2=None):
    """
    Given two sets of positions or density fields, cross-correlate.

    We define the one-dimensional CC of a 3D delta field is defined as

    < Re( \delta_1(k) \delta_2(k)^* ) > / < ( |\delta_1(k)| |\delta_2(k)| ) >

    where <.> denotes averaging in k-bins.
    
    Parameters
    ----------
    pos1, pos2: ndarrays of shape (np,ndim)
        The particle positions
    weights1, weights2: ndarrays of shape (np,)
        Particle weights (if using particles)
    delta1, delta2: ndarrays of shape `gridshape`
        The density fields
    gridshape: int or tuple of int
        The FFT mesh size.  Assumes 3D cube if an int.
        Should not include padding for an rfft.
    boxsize: float
        The physical domain size
    dtype: np.dtype, optional
        The data type with which to do the binning/FFT
    bins: int or array_like of float, optional
        Number of k-space bins, or list of bin edges.  Default: max(gridshape)//4.
    log: bool, optional
        Whether to bin in log space.  Default: False.
        
    Returns
    -------
    k: np.ndarray, shape (nbin,)
        The wavenumbers of the bins
    P_cross: np.ndarray, shape(nbin,)
        The cross correlation coefficients,
    nmodes: ndarray, shape (nbins,)
        Number of modes in each bin
    """
    
    if type(gridshape) is int:
        gridshape = (gridshape,)*3
    gridshape = np.array(gridshape)
    
    def get_delta(pos, w):
        rho = TSC.BinParticlesFromMem(pos, gridshape, boxsize, prep_rfft=False, weights=w)
        
        rho_av = np.mean(rho,dtype = np.float64).astype(dtype)
        ne.evaluate('rho/rho_av - 1', out=rho)
        return rho
        
    def get_delta_fft(delta):
        delta = _FFT(delta, boxsize, power=False)
        return delta
    
    if delta1 is not None:
        assert pos1 is None
        assert weights1 is None
    else:
        delta1 = get_delta(pos1, weights1)
    if delta2 is not None:
        assert pos2 is None
        assert weights2 is None
    else:
        delta2 = get_delta(pos2, weights2)

    delta1 = get_delta_fft(delta1)
    delta2 = get_delta_fft(delta2)
    
    P_cov = ne.evaluate('(delta1*conj(delta2)).real')
    k_cross, P_cov_binned, nmodes, bin_info = RadialBinGrid(boxsize, P_cov, bins, rfft=True)
    del P_cov

    P_norm = ne.evaluate('(abs(delta1).real*abs(delta2).real)')
    k_cross, P_norm_binned, nmodes, bin_info = RadialBinGrid(boxsize, P_norm, bins, rfft=True, bin_info=bin_info)
    del P_norm
    del delta1, delta2

    P_cross = P_cov_binned/P_norm_binned
    
    return k_cross, P_cross


@processargs(infershape=True)
def FFTAndBin(density, boxsize, *, inplace=False, bins=-1, log=False, window='window_aliased', normalize_dens=True):
    """
    Given a 2D/3D density field, convert to density contrast, FFT, square, normalize, de-window,
    and bin to produce a 1D power spectrum.
    
    Parameters
    ----------
    density: ndarray
        The (unnormalized) density field from e.g. TSC,
        which will be converted to delta = density / density.mean() - 1
    boxsize: float
        The physical domain size
    inplace: bool, optional
        Whether the density field is appropriately sized for an in-place rfft
    nbins: int, optional
        Number of k-space bins.  Default: max(gridshape)//4.
        If None, skip the radial binning and return the full 2D/3D power.
    log: bool, optional
        Whether to bin in log space.  Default: False.
    window: str, optional
        Divide out the TSC window function.
        Options are: None, 'window' and 'window_aliased' (default).
        'window' only divides out the TSC window function.
        'window_aliased' also compensates for aliasing across
        Nyquist, but assumes constant P(k) past Nyquist.
        
    Returns
    -------
    k: ndarray, shape (nbins,)
        The wavenumbers of the bins
    P: ndarray, shape (nbins,)
        Power in each bin
    nmodes: ndarray, shape (nbins,)
        Number of modes in each bin
        
    If `nbins` is None, instead the following is returned:
    power: ndarray, real
        The full 2D/3D power spectrum (Fourier transform of the deltas,
        squared and normalized)
    """
    dtype = density.dtype.type
    gridshape = np.array(density.shape)  # shape of the signal region
    # If the density array has padding, shrink the shape of the valid signal region
    if inplace:
        gridshape[-1] -= 2 if gridshape[-1] % 2 == 0 else 1
    density_real = density[...,:gridshape[-1]]  # The density without the (optional) fft padding
    rho_av = (ne.evaluate('sum(1.*density_real)')/density_real.size).astype(dtype)
    
    # Convert to fractional overdensity
    if normalize_dens:
        assert inplace
        ne.evaluate('density_real/rho_av - 1', out=density_real)
    
    if verbose:
        print("Mean Density (ppc): " + str(rho_av))
        sigma2 = ne.evaluate('sum((1.*density_real)**2)')/density_real.size
        print("RMS density fluctuation: {:.6g}".format(np.sqrt(sigma2)))

    power = _FFT(density, boxsize, window=window, inplace=inplace, power=True)
    if bins is None:
        return power
    
    k, s, nmodes, bin_info = RadialBinGrid(boxsize, power, bins, rfft=True)
    
    return k, s, nmodes


def _FFT(delta, boxsize, *, window='window_aliased', inplace=False, power=True, axes='max3'):
    '''
    Do the RFFT of the given density field.
    
    Parameters
    ----------
    delta: ndarray, 2D or 3D
        The density array.  Should probably be normalized
        already to be a density contrast "delta".
    boxsize: float
        Domain edge length
    window: str, optional
        Divide out the TSC window function.
        Options are: None, 'window' and 'window_aliased' (default).
        'window' only divides out the TSC window function.
        'window_aliased' also compensates for aliasing across
        Nyquist, but assumes constant P(k) past Nyquist.
    inplace: bool, optional
        Do the FFT of `delta` in-place.  `delta` must have
        the standard RFFT padding along the last axis.
        Default: False.
    power: bool, optional
        Return the square of the FFT, normalized to the usual
        power convention in cosmology; i.e. the raw FFTW output
        is multiplied by V/N^2, where N=nx*ny*nz.
        Default: True.
    axes: str or int or array_like of int, optional
        The axes to RFFT over.  Default of 'max3'
        choose the first min(3, input.ndim) axes.
        
    Returns
    -------
    If `power`:
    power: ndarray, real
        `delta_tilde` squared.
    
    else:
    delta_tilde: ndarray, complex
        2D/3D Fourier transform of the input `delta`.
    '''
    dtype = delta.dtype.type

    if axes is 'max3':
        axes = list(range(min(3, delta.ndim)))
    if axes is 'max2':  # for 2D displacement fields
        axes = list(range(min(2, delta.ndim)))

    try:
        axes = np.atleast_1d(np.array(axes))
        if not np.issubdtype(axes.dtype, np.integer):
            raise
    except:
        raise ValueError(axes)
    
    signalshape = np.array(delta.shape)
    if inplace:
        signalshape[axes[-1]] -= 2 if signalshape[axes[-1]]%2 == 0 else 1
        signal = delta[tuple(slice(None,ss) for ss in signalshape)]  # shorten the padded input to the real signal
    else:
        signal = delta
        
    cdtype = {np.float32: np.complex64, np.float64: np.complex128}[dtype]
        
    if inplace:
        output = delta.view(dtype=cdtype)
    else:
        outshape = signalshape.copy()
        outshape[axes[-1]] = outshape[axes[-1]]//2 + 1
        output = pyfftw.empty_aligned(outshape, dtype=cdtype)
    
    # Using FFTW_MEASURE typically results in a 3x speedup, but these FFTs are memory limited, not CPU limited
    fft_obj = pyfftw.FFTW(signal, output, axes=axes, threads=nthreads, flags=['FFTW_ESTIMATE'])
    fft_obj.execute()
    
    if power:
        # Do we want to normalize even if we aren't returning power?
        N = signalshape[axes].prod()
        norm = dtype(boxsize**signal.ndim/N**2.)
        # numexpr insists complex64.real -> float64, so force it to use a smaller array with loser casting rule
        _power = np.empty(output.shape, dtype=dtype)
        ne.evaluate('(output*conj(output)).real*norm', out=_power, casting='same_kind')
        output = _power
        del _power

    if window:
        if window == 'window':
            TSC.tsc_window(output, aliased=False, power=power)
        elif window == 'window_aliased':
            TSC.tsc_window(output, aliased=True, power=power)
        else:
            raise ValueError(window)
    
    # Might be real delta^2, might be complex delta
    return output


@processargs(infershape=True)
def InverseFFTAndBin(input, boxsize, *, bins=-1, log=False, inplace=True):
    """
    Given a 2D/3D complex array, do the inverse Fourier transform and (optionally) bin.
    
    The input array should be of the standard rfft shape.  E.g. for 3D: (gridsize, gridsize, gridsize//2 + 1)
    
    Parameters
    ----------
    input: ndarray, complex, shape (nx, ny, nz//2 + 1)
        The Fourier mesh to transform.  Should be the output of some RFFT routine.
        Real inputs will be copied into a complex array.
    boxsize: float
        The physical domain size
    bins: int or array_like, optional
        Number of radial bins, or a list of bin edges.  Default: gridsize//4.
        If None, skip the radial binning and return the full 2D/3D array.
    log: bool, optional
        Whether to bin in log space.  Default: False.
    inplace: bool, optional
        Do the transform in-place.  Note that the input is always large enough
        to do an in-place transform, and FFTW always does an IRFFTN in-place.
        So this just copies the input before passing to FFTW.  Default: False.
        
    Returns
    -------
    bin_centers: ndarray, shape (nbins,)
        The centers of the bins
    mean: ndarray, shape (nbins,)
        Mean value in each bin
    nsamples: ndarray, shape (nbins,)
        Number of samples in each bin
        
    If `bins` is None, instead the following is returned:
    output: ndarray, real, shape (nx, ny, nz)
        The full 2D/3D result of the inverse transform
    """
    output = _IFFT(input, boxsize, inplace=inplace)
    if bins is None:
        return output
    
    # We want to make an array where each entry (d1,d2,d3) is the distance d from the origin. This will help with radial binning.
    radii = misc.distance_array(input.shape, boxsize, cell_centered=False, periodic=True)
    bin_centers, mean, nsamples, bin_info = RadialBin(radii, input, bins)
    
    return bin_centers, mean, nsamples
   

def _IFFT(input, boxsize, *, inplace=False, axes='max3'):
    '''
    Do the inverse RFFT of the given input.
    
    Parameters
    ----------
    input: ndarray, complex, shape (nx, [ny,] nz//2 + 1)
        The array to IFFT.  Should be the output of some RFFT routine.
        Real inputs will be copied into a complex array.
    boxsize: float
        Domain edge length
    inplace: bool, optional
        Make a copy of the array before transforming (because
        these transforms are always in-place).
        Default: False.
    axes: str or int or array_like of int, optional
        The axes to IRFFT over.  Default of 'max3'
        choose the first min(3, input.ndim) axes.
        
    Returns
    -------
    output: ndarray, real, shape (nx, [ny,] nz)
        The result of the IFFT.
    '''
    # Parse args
    if np.iscomplexobj(input):
        cdtype = input.dtype.type
        dtype = {np.complex64:np.float32, np.complex128:np.float64}[cdtype]
    else:
        dtype = input.dtype.type
        cdtype = {np.float32:np.complex64, np.float64:np.complex128}[dtype]
        
    if axes is 'max3':
        axes = list(range(min(3, input.ndim)))
    if axes is 'max2':  # for 2D displacement fields
        axes = list(range(min(2, input.ndim)))

    try:
        axes = np.atleast_1d(np.array(axes))
        if not np.issubdtype(axes.dtype, np.integer):
            raise
    except:
        raise ValueError(axes)

    signalshape = np.array(input.shape)
    signalshape[axes[-1]] = 2*signalshape[axes[-1]] - 2  # assumes evenness of original dimensions
    paddedsignalshape = signalshape.copy()
    paddedsignalshape[axes[-1]] += 2
        
    # Copy to complex array if real input, or if we want to preserve the input
    # Can FFTW do r2r?
    if not np.iscomplexobj(input) or not inplace:
        input = input.astype(cdtype, copy=True)
    
    # output array is always the input array, but the input may already be a copy
    output = input.view(dtype=dtype).reshape(paddedsignalshape)[tuple(slice(None,ss) for ss in signalshape)]
    
    #import pdb; pdb.set_trace()
    fft_obj = pyfftw.FFTW(input, output, axes=axes, threads=nthreads, direction='FFTW_BACKWARD', flags=['FFTW_ESTIMATE'])
    fft_obj.execute()
    
    # Need to figure out if this is the desired normalization
    #norm = dtype(np.prod(signalshape)/boxsize**input.ndim)
    norm = dtype(1./np.prod(signalshape[axes]))
    ne.evaluate('output*norm', out=output)
    
    return output
