# Copyright 2012-2025 The Abacus Developers
# SPDX-License-Identifier: GPL-3.0-or-later

"""
A collection of utilities related to power spectrum computation.
"""

import numpy as np
import scipy
from Abacus.Tools import numexpr_dist
import numpy.fft as fft


pi = np.pi
    
def fftfreq_nd(shape, dx, rfft=False, mesh=False, kmag=False, dtype=np.float32):
    """
    Uses np.fft.fftfreq to return the sample frequencies for an N-D FFT.
    Returns angular frequencies, instead of the cyclic frequencies that
    the numpy fftfreq gives.
    
    Parameters
    ----------
    shape: array_like
        The shape of the signal; i.e. a.shape for np.fft.fft(a)
    dx: float or array_like
        Sample spacing along each axis (inverse of the sampling rate)
    rfft: bool, optional
        Whether the transform used fft or rfft.
    mesh: bool, optional
        Return a closed, or fleshed-out, grid instead of an open grid.
        Has no effect if using `kmag`.
    kmag: bool, optional
        Return the magnitude instead of vectors
    dtype: np.dtype, optional
        Data type of result.  Default: np.float32
        
    Returns
    -------
    ki: list of length len(shape), or ndarray of shape (nx,[ny],nz[/2+1],[3])
        An open or closed grid of the frequencies along each axis.
        If open:
        The ith entry of ki is an ndarray of len(ki[i].shape) == len(shape),
        with the length of the ith axis being the number of sample
        frequencies given by fft.[r]fftfreq(), and all other axes
        having length 1.
    """
    
    dim = len(shape)
    # Make sure we have dx for all axes
    if type(dx) is float:
        dx = np.array([dx]*dim)
    dx = np.asarray(dx)
    assert len(dx) == dim
    
    ki = []
    for i,this_dx in enumerate(dx):
        s = [1,]*dim
        s[i] = -1
        if i == dim-1 and rfft:  # all but last axis have freqs given by fftfreq
            k = fft.rfftfreq(shape[i], this_dx).reshape(s).astype(dtype)
        else:
            k = fft.fftfreq(shape[i], this_dx).reshape(s).astype(dtype)
        k *= 2*np.pi  # Convert to angular freq
        ki.append(k)
        
    if kmag:
        mag = numexpr_dist(ki)
        return mag
        
    if mesh:
        # Note: this internally converts to double
        ki = np.stack(np.meshgrid(*ki, indexing='ij'), axis=-1).astype(dtype)
        return ki
    
    return ki


def loglog_interp(k, Pk, logx=True, logy=True, **kwargs):
    """
    A convenience function to do cubic spline interpolation in log-log space.
    Useful for iterpolating power spectra.
    
    """
    kind = kwargs.pop('kind', 'cubic')

    target = np.log(Pk) if logy else Pk
    x = np.log(k) if logx else k
    
    loglog_input_interp = scipy.interpolate.interp1d(x, target, kind=kind, **kwargs)
    if logy:
        interp_func = lambda k: np.exp(loglog_input_interp(np.log(k) if logx else k))
    else:
        interp_func = lambda k: loglog_input_interp(np.log(k) if logx else k)
    return interp_func


import numba as nb
def parallel_bcast(ftylist, signature, **kwargs):
    '''
    A decorator to broadcast a function over the first axis of an array in parallel.
    The decorator call signature is identical to `numba.guvectorize()`.
    The decorated function should take an array and a loop index.
    Note that the loop index does not have to appear in the guvectorize args;
    this will be filled in automatically.
    
    The number of threads can only be set by setting the NUMBA_NUM_THREADS
    environment variable before importing numba for the first time.
    
    Inspired by:
    http://stackoverflow.com/questions/35541289/possible-to-use-numba-guvectorize-to-emulate-parallel-forall-prange
    
    Parameters
    ----------
    Same as `numba.guvectorize()`, but with the defaults:
    ```
    nopython = True
    target   = 'parallel'
    ```
    
    Example
    -------
    # Using numba.guvectorize
    # Note this cannot run in parallel, even with target='parallel', because no broadcasting is occurring
    
    @nb.guvectorize([(nb.float32[:,:,:],), (nb.float64[:,:,:],)], '(nx,ny,nz)', nopython=True, target='parallel')
    def jitwindow(power):
        nx,ny,nz = power.shape
        for i in range(nx):
            kx = np.pi*i/nx
            for j in range(ny):
                ky = np.pi*j/ny
                for k in range(nz):
                    kz = np.pi*k/nz
                    power[i,j,k] /= (np.sinc(kx)*np.sinc(ky)*np.sinc(kz))**2.
                
    # Equivalent code using parallel_bcast, which will use multiple CPU cores:
    
    @parallel_bcast([(nb.float32[:,:,:],), (nb.float64[:,:,:],)], '(nx,ny,nz)')
    def window3D_inplace(power, loop_idx):
        i = loop_idx[0]
        nx,ny,nz = power.shape
        kx = np.pi*i/nx
        for j in range(ny):
            ky = np.pi*j/ny
            for k in range(nz):
                kz = np.pi*k/nz
                power[i,j,k] /= (np.sinc(kx)*np.sinc(ky)*np.sinc(kz))**2.
    '''
    
    ftylist = [t + (nb.int64[:],) for t in ftylist]  # could add int32, uints
    
    # Add the loop index to the signature
    signature += ' '
    i = signature.find('->')
    signature = signature[:i] + ',()' + signature[i:]
    
    # Set default kwargs
    if 'nopython' not in kwargs:
        kwargs['nopython'] = True
    if 'target' not in kwargs:
        kwargs['target'] = 'parallel'
        
    def decorator(func):
        bcast_func = nb.guvectorize(ftylist, signature, **kwargs)(func)
        def wrapper(*fargs, **fkwargs):
            loop_inds = list(range(len(fargs[0])))  # could bcast over all arrays
            fargs = fargs + (loop_inds,)
            bcast_func(*fargs, **fkwargs)
        return wrapper
    return decorator

# Load all the functions defined in the cffi library
# TODO: hook numba up to pybind11 functions
#from numba.core.typing import cffi_utils
#from . import _psffilib
#for var in vars(_psffilib.lib):
#    vars()[var] = getattr(_psffilib.lib, var)
#cffi_utils.register_module(_psffilib)
#from _psffilib import ffi as psffi, lib as psffilib

# @parallel_bcast([(nb.float32[:,:,:], nb.float32[:], nb.float32[:], nb.float32[:])], '(nx,ny,nz),(),(),()')
# def shell_fft(out, boxsize, rmax, rmin, loop_idx):
#     '''
#     Analytically evaluate the FFT of the 3D shell-averaging function:
#     4*pi/k*(rmax^2*j_1(rmax*k) - rmin^2*j_1(rmin*k))/(4*pi/3*(rmax^3 - rmin^3))
#     '''
#     nx,ny,nz = out.shape
#     boxsize = boxsize[0]
#     rmax, rmin = rmax[0], rmin[0]
#     i = loop_idx[0]
#     kx2 = (2*np.pi*i/boxsize)**2
#     for j in range(ny):
#         ky2 = (2*np.pi*j/boxsize)**2
#         for k in range(nz):
#             kz2 = (2*np.pi*k/boxsize)**2
#             kmag = np.sqrt(kx2 + ky2 + kz2)
#             out[i,j,k] = 4*np.pi/kmag * (rmax**2*gsl_sf_bessel_j1(rmax*kmag) - rmin**2*gsl_sf_bessel_j1(rmin*kmag)) \
#                          / (4*np.pi/3*(rmax**3 - rmin**3))
#     out[0,0,0] = 1.  # fix the divide by k=0

# @parallel_bcast([(nb.float32[:,:], nb.float32[:], nb.float32[:], nb.float32[:])], '(nx,ny),(),(),()')
# def annulus_fft(out, boxsize, smax, smin, loop_idx):
#     '''
#     Analytically evaluate the FFT of the 2D annulus-averaging function:
#     2*pi/k_perp*(smax*J_1(smax*k) - smin*J_1(smin*k))/(pi*(smax^2 - smin^2))
#     '''
#     nx,ny = out.shape
#     boxsize = boxsize[0]
#     smax, smin = smax[0], smin[0]
#     i = loop_idx[0]
#     kx2 = (2*np.pi*i/boxsize)**2
#     for j in range(ny):
#         ky2 = (2*np.pi*j/boxsize)**2
#         k_perp = np.sqrt(kx2 + ky2)
#         out[i,j] = 2*np.pi/k_perp*(smax*gsl_sf_bessel_J1(smax*k_perp) - smin*gsl_sf_bessel_J1(smin*k_perp)) \
#                    / (pi*(smax**2 - smin**2))
#     out[0,0] = 1.  # fix the divide by k=0


def plot_Pks(ks, Pks=None, cross_Pks=None, deltas=None, deltaks=None, labels=None):
    '''
    Constructs (typically) a three-panel plot of power, ratio of power,
    and cross-correlation.

    TODO: finish this
    '''

    if (Pks is None) != (deltaks is None) != (deltas is None):
        raise ValueError('Must specify exactly one of `Pks`, `deltaks`, `deltas`.')

    nrows = 3 if cross_Pks is not None else 2
    fig, axes = plt.subplots(nrows, 1, sharex=True, figsize=(9,9))

    # broadcast ks
    #if len(ks) != nlines:

    if labels is None:
        labels = [str(i) for i in range(nlines)]

    # power
    l1, = axes[0].loglog(ks[0], Pks[0], label=labels[0])
    l1_disp, = ax1.loglog(k1_disp, Pk1_disp, label='ZA displacement')
    l2, = ax1.loglog(k2, Pk2, label='Reference ($N \geq {}$)'.format(Ncut2), color='k')
    #ax1.axhline(box**3/len(halo_particles1), ls='--', c=l1.get_color())
    ax1.set_ylabel('$P(k)$')
    ax1.tick_params(right=True,top=True)

    # ratio of power
    ax2.plot(k1, np.abs(Pk1/Pk2 - 1), color=l1.get_color())
    ax2.plot(k1, np.abs(Pk1_disp/Pk2 - 1), color=l1_disp.get_color())
    ax2.set_ylabel(r'$|P(k)/P_\mathrm{ref}(k) - 1|$')
    ax2.tick_params(right=True,top=True)
    #ax2.axhline(1.,ls=':',c='k')
    ax2.set_xscale('log'); ax2.set_yscale('log')

    # cross-corr
    ax3.semilogx(k_cross1, P_cross1, color=l1.get_color())
    ax3.semilogx(k_cross1, P_cross1_disp, color=l1_disp.get_color())
    ax3.set_ylim(top=1.0005)
    ax3.set_xlabel(r'$k$ [$h$/Mpc]')
    ax3.set_ylabel('Cross correlation')
    ax3.tick_params(right=True,top=True)
    ax3.axhline(1.,ls=':',c='k')

    # try to plot result with the transfer function applied
    try:
        l1_disp_trans, = ax1.loglog(k1_disp_trans, Pk1_disp_trans, label='ZA disp. + aniso. transfer')
        ax2.plot(k1, np.abs(Pk1_disp_trans/Pk2 - 1), color=l1_disp_trans.get_color())
        ax3.semilogx(k_cross1, P_cross1_disp_trans, color=l1_disp_trans.get_color())
        
        l1_disp_trans_smooth, = ax1.loglog(k1_disp_trans_smooth, Pk1_disp_trans_smooth, \
                                           label='ZA disp. + aniso. transfer + {:g} Mpc/h smooth'.format(smooth))
        ax2.plot(k1, np.abs(Pk1_disp_trans_smooth/Pk2 - 1), color=l1_disp_trans_smooth.get_color())
        ax3.semilogx(k_cross1, P_cross1_disp_trans_smooth, color=l1_disp_trans_smooth.get_color())
    except:
        print('Skipping transfer line')

    ax1.legend()
    fig.tight_layout()
    fig.subplots_adjust(hspace=0.)

    return fig, (ax1, ax2, ax3)

import functools
from . import Histogram
def processargs(original_func=None, infershape=False):
    '''
    Used to decorate functions in the PowerSpectrum class.
    Allows them to flexibly take `bins`/`nbins`.
    '''
    def _processargs(func):
        '''
        - Accept `bins`, `nbins`, but not both.
          `bins` and `nbins` will be accessible from the
          function namespace.
          `bins` will always be an array of bin edges.

        TODO: dens, delta, pos
        TODO: nthreads, verbose
        '''
        defaults = func.__kwdefaults__
        
        @functools.wraps(func)  # so meta
        def wrapper(*args, **_kwargs):
            # Get all kwargs, including defaults if not overridden
            kwargs = defaults.copy()
            kwargs.update(_kwargs)

            # All func sigs must be f(...,gridshape,box,*,...)
            # We need these to automatically populate bin_edges from nbins
            boxsize = args[-1]
            if infershape:
                gridshape = args[-2].shape
            else:
                gridshape = args[-2]

            pass_via_globals = {}

            if 'nbins' in kwargs:
                # Can't specify both kwargs
                assert 'bins' not in kwargs
                # nbins must be int
                assert isinstance(kwargs['nbins'], (np.integer,int))
                
                pass_via_globals['nbins'] = kwargs.pop('nbins')
                
            elif 'bins' in kwargs:
                # Can't specify both kwargs
                assert 'nbins' not in kwargs
                bins = kwargs.pop('bins')
                isint = isinstance(bins, (np.integer,int))
                if isint:
                    pass_via_globals['nbins'] = bins
                elif bins is None:
                    kwargs['bins'] = bins
                else:
                    # wrap lists in arrays
                    bins = np.asanyarray(bins)
                    
                    # Cast integer arrays to float arrays, just to be safe
                    isfloat = issubclass(bins.dtype.type, (float,np.floating))
                    if not isfloat:
                        bins = bins.astype(np.float64)
                    kwargs['bins'] = bins

            if 'bins' not in kwargs and 'nbins' in pass_via_globals:
                bins = Histogram.k_bin_edges(gridshape, boxsize, nbins=pass_via_globals['nbins'], log=kwargs['log'])
                kwargs['bins'] = bins

            func.__globals__.update(pass_via_globals)
            
            return func(*args, **kwargs)

        return wrapper
    if original_func:
        return _processargs(original_func)
    else:
        return _processargs


def process_gridshape(gridshape):
    gridshape = np.atleast_1d(gridshape)
    assert np.issubdtype(gridshape.dtype, np.integer)
    if len(gridshape) == 1:
        gridshape = np.repeat(gridshape, 3)
    gridshape = np.array(gridshape)

    return gridshape
