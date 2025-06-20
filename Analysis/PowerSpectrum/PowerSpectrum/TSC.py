# Copyright 2012-2025 The Abacus Developers
# SPDX-License-Identifier: GPL-3.0-or-later

'''
A Python implementation of triangle-shaped cloud
binning, using Numba for C-like efficiency.
'''
from glob import glob
import warnings
import os.path
from os.path import join as pjoin
import re
import queue
import os
import threading
import gc

import numpy as np
import numba

from .misc import process_gridshape

from Abacus import ReadAbacus
from Abacus.Tools import ContextTimer

from . import pslib

valid_file_formats = list(ReadAbacus.reader_functions.keys())
maxthreads = os.cpu_count()

def BinParticlesFromMem(positions, gridshape, boxsize, weights=None, dtype=np.float32, rotate_to=None, prep_rfft=False, nthreads=-1, inplace=False, norm=False):
    """
    Computes the density field of a set of points by binning them into a grid using TSC.
    
    Parameters
    ----------
    positions: ndarray of shape (N,3), or list of ndarray
        The points to bin
    gridshape: int or tuple of int
        The FFT mesh size.  Assumes 3D cube if an int.
        Should not include padding for an rfft; use `prep_rfft`.
    boxsize: float
        The physical domain size.
        `positions` are assumed to be in a zero-centered box of side length `boxsize`.
    weights: ndarray of shape (N,), or list of ndarray, optional
        Binning weights for the points.
    dtype: np.dtype, optional
        The data type of the density array.
    rotate_to: array of floats of shape (3,), optional
        Rotate the particles before binning.
        Currently not implemented.
        Default: None.
    prep_rfft: bool, optional
        Adjust the input `gridshape` to add padding
        of the correct shape for an in-place rfft.
        Default: False.
    nthreads: int, optional
        Number of threads.  Default of -1 uses all available cores.
    inplace: bool, optional
        Whether we are allowed to modify the `positions` and `weights`
        arrays for performance reasons.
    norm: bool, optional
        Normalize the returned array to a density contrast "delta".
        Default: False.
    
    Returns
    -------
    density: ndarray
        The density field, in units of (weighted) counts per cell.
        Shape is determined by `gridshape` and `prep_rfft`.
    """
    
    # Wrap lone arrays inside lists, so everything will be a list hereafter
    if type(positions) != list:
        positions = [positions]
    if weights is None:
        weights = [None]*len(positions)
    if type(weights) != list:
        weights = [weights]
    assert len(positions) == len(weights)
        
    NPs = [len(p) for p in positions]
        
    # make int gridshape a cube
    gridshape = process_gridshape(gridshape)
    gridshape_real = gridshape.copy()
        
    # add padding beyond what was passed in
    if prep_rfft:
        gridshape[-1] = 2*(gridshape[-1]//2 + 1)
        
    density = np.zeros(gridshape, dtype=dtype)
    density_real = density[...,:gridshape_real[-1]]

    # get the density weighted field grid
    for p,w,NP in zip(positions,weights,NPs):
        TSC(p, density, boxsize, weights=w, prep_rfft=prep_rfft, rotate_to=rotate_to, nthreads=nthreads, inplace=inplace)

    if norm:
        rho_av = density_real.mean(dtype=np.float64).astype(dtype)
        ne.evaluate('density_real/rho_av - 1', out=density_real)
    
    return density
    

def BinParticlesFromFile(file_pattern, boxsize, gridshape, dtype=np.float32, zspace=False, format='pack14', rotate_to=None,
                         prep_rfft=False, nthreads=-1, nreaders=1, readahead=1, verbose=False, return_NP=False):
    '''
    Main entry point for density field computation from files on disk of various formats.
    
    The reading and binning is done slab-by-slab, so we don't need to fit the whole particle set
    in memory.  The reading is overlapped with the binning.

    Parameters:
    -----------
    file_pattern: str
        A globbing pattern for all the files to read
    boxsize: float
        Size of the box.  Particles (except for RVzel and Gadget format) are assumed to be stored in unit box.
    gridshape: int or tuple of int
        The FFT mesh size.  Assumes 3D cube if an int.
        Should not include padding for an rfft; use `prep_rfft`.
    dtype: np.dtype, optional
        Precision for the density array.  Default: np.float32.
    zspace: bool, optional
        Whether to bin the particles in redshift space.  Default: False.
    format: str or callable, optional
        File format for the input files.  Options are 'rvdouble', 'pack14', 'rvzel', 'state', 'gadget'.
        Can also be a function that takes a filename string and returns a (N,d) position array.
        Default: 'pack14'
    rotate_to: array of floats of shape (3,), optional
        `rotate_to` defines a vector that the cartesian z-hat will be rotated to align with.
        All rotations are periodic "safe", meaning we discard particles at the edges.
        Default: None.
    prep_rfft: bool, optional
        Adjust the input `gridshape` to add padding
        of the correct shape for an in-place rfft.
        Default: False.
    nthreads: int, optional
        Number of threads.  Default of -1 uses all available cores.
    nreaders: int, optional
        Number of IO threads. Default: 1
    readahead: int, optional
        How many files to read ahead of what's been binned. -1 allows unbounded readahead.
        Default: 1.
    return_NP: bool, optional
        Return the number of particles read

    Returns
    -------
    density: ndarray
        The density field, in units of (weighted) counts per cell.
        Shape is determined by `gridshape` and `prep_rfft`.
    '''
        
    if type(format) == str:
        format = format.lower()
    if format not in valid_file_formats and not callable(format):
        raise ValueError(format, 'Use one of: {}, or callable'.format(valid_file_formats))
    
    # make int gridshape a cube
    gridshape = process_gridshape(gridshape)
        
    # add padding beyond what was passed in
    if prep_rfft:
        gridshape[-1] = 2*(gridshape[-1]//2 + 1)
        
    density = np.zeros(gridshape, dtype=dtype)
    
    # If rotating the particles, only the inner "rotation safe" region is binned
    # So, shrink the boxsize accordingly
    if rotate_to:
        rotate_to = np.array(rotate_to, dtype=np.float64)
        boxsize /= 3**.5

    if nthreads < 1:
        nthreads = maxthreads
    ReadAbacus_threads = nreaders*ReadAbacus.get_nthreads_by_format(format)
    nthreads = max(1,nthreads-ReadAbacus_threads)
    if verbose:
        print(f'TSC with {nthreads} threads, {ReadAbacus_threads} ReadAbacus threads ({nreaders} IO threads)')


    reader_kwargs = dict(return_pos=True, return_vel=False, zspace=zspace, dtype=dtype, format=format,
                         units=None, readahead=readahead, verbose=verbose, nreaders=nreaders)

    if format == 'rvzel':
        reader_kwargs.update(return_zel=False, add_grid=True, boxsize=boxsize)
    elif format == 'rvdoublezel':
        reader_kwargs.update(return_zel=False, add_grid=True, boxsize=boxsize)
    elif format == 'gadget':
        reader_kwargs.update(boxsize=boxsize)
    elif format == 'state':
        reader_kwargs.update(boxsize=boxsize)
    elif format == 'rvint':
        reader_kwargs.update(boxsize=boxsize)
        
    # Rather than convert units, do the TSC in the native units
    tsc_box = ReadAbacus.get_box_on_disk(file_pattern, format=format)
    if tsc_box == 'box':
        tsc_box = boxsize

    np_read = 0
    tsc_timer = ContextTimer('Cumulative TSC binning', cumulative=True, output=True)
    for data in ReadAbacus.AsyncReader(file_pattern, **reader_kwargs):
        np_read += len(data['pos'])
        with tsc_timer:
            TSC(data['pos'], density, tsc_box, rotate_to=rotate_to, prep_rfft=prep_rfft, nthreads=nthreads, inplace=True)
        del data; gc.collect()
    if verbose:
        print(f'TSC read {np_read} particles = ppd {np_read**(1/3)}')
        print(f'Total TSC binning took {tsc_timer.elapsed:.4g} sec')
    
    if return_NP:
        return density, np_read
    return density


import numexpr as ne
# TODO: wrap modes
def TSC(positions, density, boxsize, weights=None, prep_rfft=False, rotate_to=None, wrap='wrap', nthreads=-1, inplace=False):
    '''
    Main entry point for TSC binning of particles into an existing
    density array.  See `BinParticlesFromMem` to create a density
    array.
    
    This function serves as a wrapper to choose the 2D or 3D Numba
    function, and to do the parallelization.

    Particles are sorted into stripes and threads operate on
    non-adjacent stripes.
    
    Parameters
    ----------
    positions: ndarray of shape (N,3) or (N,2)
        The points to bin
    boxsize: float
        The physical domain size.
        `positions` are assumed to be in a zero-centered box of side length `boxsize`.
    weights: scalar or ndarray of shape (N,), optional
        Binning weights for the points.  Scalar weights are applied too all particles.
    prep_rfft: bool, optional
        Whether the density array has extra padding on the last axis
        for an in-place RFFT.  Particles will not be binned there.
        Default: False.
    rotate_to: array of floats of shape (3,), optional
        Rotate the particles before binning.
        Currently not implemented.
        Default: None.
    nthreads: int, optional
        Number of threads.  Default of -1 uses all available cores.
    inplace: bool, optional
        Whether we are allowed to modify the `positions` and `weights`
        arrays for performance reasons.
    wrap: 
    '''
    if rotate_to:
        raise NotImplementedError('rotate_to')
    if nthreads < 1:
        nthreads = maxthreads
    if weights is not None:
        weights = np.atleast_1d(weights)

    # cast to same dtype
    boxsize = positions.dtype.type(boxsize)
    
    # We may eventually want to provide two wrapping modes
    # 'wrap' and 'clip' (mostly for mocks)
    # Warning: we do this regardless of inplace or not! Probably okay.
    box_wrap(positions, boxsize)
        
    if prep_rfft:
        gz = density.shape[-1]
        gz -= 2 if gz % 2 == 0 else 1
        density = density[...,:gz]

    # Now organize the particles into stripes/chunks
    # chunks must be wider than 3 cells to avoid TSC overlap
    # we want the fewest chunks that is still more than Nthreads*2
    nchunks = len(density)//3  # this is the absolute max
    nchunks = min(nchunks, nthreads*2)

    if nthreads == 1:
        nchunks = 1  # no point

    if nchunks > 1:
        # Partition the particles into chunks
        # The algorithm we really need here is a parallel partition that splits on y-values (not indices)
        # But parallel sort algorithms are more readily available

        ybins = np.linspace(-boxsize/2, boxsize/2, num=nchunks+1, endpoint=True, dtype=positions.dtype)
        positions, weights = sort_pos_and_weight(positions, weights, inplace=inplace)
        splits = np.searchsorted(positions[:,1], ybins[1:-1])
        pchunks = np.array_split(positions, splits)

        if weights is None or len(weights) == 0:
            wchunks = [np.empty(0),]*nchunks
        else:
            if len(weights) > 1:
                wchunks = np.array_split(weights, splits)
            elif len(weights) == 1:
                wchunks = [weights,]*nchunks
            
    else:
        # don't do expensive sorting if we don't need to
        pchunks = [positions]
        wchunks = [weights if weights is not None else np.empty(0)]

    assert len(pchunks) == nchunks
    
    simple = weights is None and density.ndim == 3
    numba_func = numba_tsc_3D_simple if simple else numba_tsc_3D

    if density.ndim == 3:
        spawn_thread = lambda i: threading.Thread(target=numba_func, args=(pchunks[i], density, boxsize), kwargs={'weights':wchunks[i]} if not simple else None)
    elif density.ndim == 2:
        spawn_thread = lambda i: threading.Thread(target=numba_func, args=(pchunks[i], density[...,None], boxsize), kwargs={'weights':wchunks[i]} if not simple else None)
    else:
        raise ValueError(density.ndim)

    # First even, then odd stripes
    set1 = list(range(0,nchunks-1,2))
    set2 = list(range(1,nchunks  ,2))
    last = None if nchunks % 2 == 0 else -1
    
    threads = [spawn_thread(i) for i in set1]
    for thread in threads:
        thread.start()
    for thread in threads:
        thread.join()
        
    threads = [spawn_thread(i) for i in set2]
    for thread in threads:
        thread.start()
    for thread in threads:
        thread.join()
        
    # Now do the last stripe if we had an odd number
    if last is not None:
        thread = spawn_thread(last)
        thread.start()
        thread.join()

# an in-place fast periodic wrap
import numba as nb
@nb.njit(parallel=True)
def box_wrap(p, box):
    p = p.reshape(-1)
    N = len(p)
    for i in nb.prange(N):
        if p[i] >= box/2:
            p[i] -= box
        if p[i] < -box/2:
            p[i] += box
        if p[i] >= box/2:
            p[i] -= box
        if p[i] < -box/2:
            p[i] += box


def sort_pos_and_weight(p, w, inplace):
    # our FFI lib uses 32 bit floats right now
    # our FFI lib assumes 3D right now
    if p.dtype != np.float32 or p.shape[-1] != 3:
        warnings.warn("Warning: particles not float32 in 3D. Falling back to slow sort.")
        order = p[:,1].argsort()
        p = p[order]
        if w is not None and len(w) > 1:
            w = w[order]

    elif w is not None and len(w) > 1:
        assert w.dtype == p.dtype

        # pack into float4
        # possibly a parallel argsort would be better; could try sorting an index array
        # inplace doesn't help here
        pw = np.empty((len(p), p.shape[-1] + 1), dtype=p.dtype)
        pw[:,:3] = p
        pw[:,3] = w
        pslib.y_sorter_weighted(pw)
        # unpack
        # this returns views, but will it make TSC slow?
        p = pw[:,:3]
        w = pw[:,3]
    else:
        if not inplace:
            p = p.copy()
        pslib.y_sorter(p)

    return p, w


@numba.vectorize
def rightwrap(x, L):
    if x >= L:
        return x - L
    return x
    
# This can be called with 2D or 3D positions and densities	
# The 2D density field must have 3 axes, with the last axis flat
# Only the x,y positions will be used in the 2D case
# TODO: rewrite using Numba stencils? They currently don't support periodicity.
@numba.jit(nopython=True, nogil=True, fastmath=False)
def numba_tsc_3D(positions, density, boxsize, weights=np.empty(0)):
    ftype = positions.dtype.type
    itype = np.int16
    gx = itype(density.shape[0])
    gy = itype(density.shape[1])
    gz = itype(density.shape[2])
    
    # because float32 * int32 = float64 in numba...
    fgx = ftype(gx)
    fgy = ftype(gy)
    fgz = ftype(gz)
    
    threeD = gz != 1
    W = ftype(1.)
    Nw = len(weights)
    for n in range(len(positions)):
        # broadcast scalar weights
        if Nw == 1:
            W = ftype(weights[0])  # todo, possibly unsafe
        elif Nw > 1:
            W = ftype(weights[n])
        
        # convert to a position in the grid
        px = (positions[n,0]/boxsize + ftype(.5))*fgx
        py = (positions[n,1]/boxsize + ftype(.5))*fgy
        if threeD:
            pz = (positions[n,2]/boxsize + ftype(.5))*fgz
        
        # round to nearest cell center
        ix = itype(round(px))
        iy = itype(round(py))
        if threeD:
            iz = itype(round(pz))
        
        # calculate distance to cell center
        dx = ftype(ix) - px
        dy = ftype(iy) - py
        if threeD:
            dz = ftype(iz) - pz
        
        # find the tsc weights for each dimension
        wx = ftype(.75) - dx**2
        wxm1 = ftype(.5)*(ftype(.5) + dx)**2
        wxp1 = ftype(.5)*(ftype(.5) - dx)**2
        wy = ftype(.75) - dy**2
        wym1 = ftype(.5)*(ftype(.5) + dy)**2
        wyp1 = ftype(.5)*(ftype(.5) - dy)**2
        if threeD:
            wz = ftype(.75) - dz**2
            wzm1 = ftype(.5)*(ftype(.5) + dz)**2
            wzp1 = ftype(.5)*(ftype(.5) - dz)**2
        else:
            wz = ftype(1.)
        
        # find the wrapped x,y,z grid locations of the points we need to change
        # negative indices will be automatically wrapped
        ixm1 = (ix - itype(1))
        ixw  = rightwrap(ix    , gx)
        ixp1 = rightwrap(ix + itype(1), gx)
        iym1 = (iy - itype(1))
        iyw  = rightwrap(iy    , gy)
        iyp1 = rightwrap(iy + itype(1), gy)
        if threeD:
            izm1 = (iz - itype(1))
            izw  = rightwrap(iz    , gz)
            izp1 = rightwrap(iz + itype(1), gz)
        else:
            izw = itype(0)
        
        # change the 9 or 27 cells that the cloud touches
        density[ixm1, iym1, izw ] += wxm1*wym1*wz  *W
        density[ixm1, iyw , izw ] += wxm1*wy  *wz  *W
        density[ixm1, iyp1, izw ] += wxm1*wyp1*wz  *W
        density[ixw , iym1, izw ] += wx  *wym1*wz  *W
        density[ixw , iyw , izw ] += wx  *wy  *wz  *W
        density[ixw , iyp1, izw ] += wx  *wyp1*wz  *W
        density[ixp1, iym1, izw ] += wxp1*wym1*wz  *W
        density[ixp1, iyw , izw ] += wxp1*wy  *wz  *W
        density[ixp1, iyp1, izw ] += wxp1*wyp1*wz  *W
        
        if threeD:
            density[ixm1, iym1, izm1] += wxm1*wym1*wzm1*W
            density[ixm1, iym1, izp1] += wxm1*wym1*wzp1*W

            density[ixm1, iyw , izm1] += wxm1*wy  *wzm1*W
            density[ixm1, iyw , izp1] += wxm1*wy  *wzp1*W

            density[ixm1, iyp1, izm1] += wxm1*wyp1*wzm1*W
            density[ixm1, iyp1, izp1] += wxm1*wyp1*wzp1*W

            density[ixw , iym1, izm1] += wx  *wym1*wzm1*W
            density[ixw , iym1, izp1] += wx  *wym1*wzp1*W

            density[ixw , iyw , izm1] += wx  *wy  *wzm1*W
            density[ixw , iyw , izp1] += wx  *wy  *wzp1*W

            density[ixw , iyp1, izm1] += wx  *wyp1*wzm1*W
            density[ixw , iyp1, izp1] += wx  *wyp1*wzp1*W

            density[ixp1, iym1, izm1] += wxp1*wym1*wzm1*W
            density[ixp1, iym1, izp1] += wxp1*wym1*wzp1*W

            density[ixp1, iyw , izm1] += wxp1*wy  *wzm1*W
            density[ixp1, iyw , izp1] += wxp1*wy  *wzp1*W

            density[ixp1, iyp1, izm1] += wxp1*wyp1*wzm1*W
            density[ixp1, iyp1, izp1] += wxp1*wyp1*wzp1*W
            

@numba.jit(nopython=True, nogil=True, fastmath=False)
def numba_tsc_3D_simple(positions, density, boxsize):
    ftype = positions.dtype.type
    itype = np.int16
    gx = itype(density.shape[0])
    gy = itype(density.shape[1])
    gz = itype(density.shape[2])
    
    # because float32 * int32 = float64 in numba...
    fgx = ftype(gx)
    fgy = ftype(gy)
    fgz = ftype(gz)
    
    for n in range(len(positions)):
        # convert to a position in the grid
        px = (positions[n,0]/boxsize + ftype(.5))*fgx
        py = (positions[n,1]/boxsize + ftype(.5))*fgy
        pz = (positions[n,2]/boxsize + ftype(.5))*fgz
        
        # round to nearest cell center
        ix = itype(round(px))
        iy = itype(round(py))
        iz = itype(round(pz))
        
        # calculate distance to cell center
        dx = ftype(ix) - px
        dy = ftype(iy) - py
        dz = ftype(iz) - pz
        
        # find the tsc weights for each dimension
        wx = ftype(.75) - dx**2
        wxm1 = ftype(.5)*(ftype(.5) + dx)**2
        wxp1 = ftype(.5)*(ftype(.5) - dx)**2
        wy = ftype(.75) - dy**2
        wym1 = ftype(.5)*(ftype(.5) + dy)**2
        wyp1 = ftype(.5)*(ftype(.5) - dy)**2
        
        wz = ftype(.75) - dz**2
        wzm1 = ftype(.5)*(ftype(.5) + dz)**2
        wzp1 = ftype(.5)*(ftype(.5) - dz)**2
        
        # find the wrapped x,y,z grid locations of the points we need to change
        # negative indices will be automatically wrapped
        ixm1 = (ix - itype(1))
        ixw  = rightwrap(ix    , gx)
        ixp1 = rightwrap(ix + itype(1), gx)
        iym1 = (iy - itype(1))
        iyw  = rightwrap(iy    , gy)
        iyp1 = rightwrap(iy + itype(1), gy)
        
        izm1 = (iz - itype(1))
        izw  = rightwrap(iz    , gz)
        izp1 = rightwrap(iz + itype(1), gz)
        
        # change the 9 or 27 cells that the cloud touches
        density[ixm1, iym1, izm1] += wxm1*wym1*wzm1
        density[ixm1, iym1, izw ] += wxm1*wym1*wz
        density[ixm1, iym1, izp1] += wxm1*wym1*wzp1
        
        density[ixm1, iyw , izm1] += wxm1*wy  *wzm1
        density[ixm1, iyw , izw ] += wxm1*wy  *wz  
        density[ixm1, iyw , izp1] += wxm1*wy  *wzp1
        
        density[ixm1, iyp1, izm1] += wxm1*wyp1*wzm1
        density[ixm1, iyp1, izw ] += wxm1*wyp1*wz  
        density[ixm1, iyp1, izp1] += wxm1*wyp1*wzp1
        
        density[ixw , iym1, izm1] += wx  *wym1*wzm1
        density[ixw , iym1, izw ] += wx  *wym1*wz  
        density[ixw , iym1, izp1] += wx  *wym1*wzp1
        
        density[ixw , iyw , izm1] += wx  *wy  *wzm1
        density[ixw , iyw , izw ] += wx  *wy  *wz  
        density[ixw , iyw , izp1] += wx  *wy  *wzp1
        
        density[ixw , iyp1, izm1] += wx  *wyp1*wzm1
        density[ixw , iyp1, izw ] += wx  *wyp1*wz  
        density[ixw , iyp1, izp1] += wx  *wyp1*wzp1
        
        density[ixp1, iym1, izm1] += wxp1*wym1*wzm1
        density[ixp1, iym1, izw ] += wxp1*wym1*wz  
        density[ixp1, iym1, izp1] += wxp1*wym1*wzp1
        
        density[ixp1, iyw , izm1] += wxp1*wy  *wzm1
        density[ixp1, iyw , izw ] += wxp1*wy  *wz  
        density[ixp1, iyw , izp1] += wxp1*wy  *wzp1
        
        density[ixp1, iyp1, izm1] += wxp1*wyp1*wzm1
        density[ixp1, iyp1, izw ] += wxp1*wyp1*wz          
        density[ixp1, iyp1, izp1] += wxp1*wyp1*wzp1

# The following are window functions for TSC-derived power spectra
# TODO: re-write these with nb.prange

from .misc import parallel_bcast
import numba as nb
@parallel_bcast([(nb.float32[:,:,:],nb.bool_), (nb.float64[:,:,:],nb.bool_)], '(nx,ny,nz),()')
def tsc_window_3d(deltak, power, loop_idx):
    i = loop_idx[0]
    nx,ny,nz = deltak.shape
    nz_signal = 2*nz - 2  # warning: assumes evenness of original dimensions
    # normally these have a factor of pi, but the numpy sinc already includes it
    # also we don't need to worry about reflection for sin, but it matters for sinc
    kx = 1.*i/nx if i < nx//2 else 1.*(i-nx)/nx
    xfac = np.sinc(kx)**6.
    for j in range(ny):
        ky = 1.*j/ny if j < ny//2 else 1.*(j-ny)/ny
        yfac = np.sinc(ky)**6.
        for k in range(nz):
            kz = 1.*k/nz_signal
            zfac = np.sinc(kz)**6.
            if power:
                deltak[i,j,k] /= xfac*yfac*zfac  # could very easily overflow float32 max exponent
            else:
                deltak[i,j,k] /= (xfac*yfac*zfac)**.5
            
@parallel_bcast([(nb.float32[:,:],nb.bool_), (nb.float64[:,:],nb.bool_)], '(nx,ny),()')
def tsc_window_2d(deltak, power, loop_idx):
    i = loop_idx[0]
    nx,ny = deltak.shape
    ny_signal = 2*ny - 2  # warning: assumes evenness of original dimensions
    kx = 1.*i/nx
    xfac = np.sinc(kx)**6.
    for j in range(ny):
        ky = 1.*j/ny_signal
        yfac = np.sinc(ky)**6.
        if power:
            deltak[i,j] /= xfac*yfac
        else:
            deltak[i,j] /= (xfac*yfac)**.5
            
@parallel_bcast([(nb.float32[:,:,:],nb.bool_), (nb.float64[:,:,:],nb.bool_),
                 (nb.complex64[:,:,:],nb.bool_), (nb.complex128[:,:,:],nb.bool_)], '(nx,ny,nz),()')
def tsc_window_aliased_3d(deltak, power, loop_idx):
    i = loop_idx[0]
    nx,ny,nz = deltak.shape
    nz_signal = 2*nz - 2  # warning: assumes evenness of original dimensions
    kx = np.pi*i/nx
    xfac = (1 - np.sin(kx)**2 + 2*np.sin(kx)**4/15)
    for j in range(ny):
        ky = np.pi*j/ny
        yfac = (1 - np.sin(ky)**2 + 2*np.sin(ky)**4/15)
        for k in range(nz):
            kz = np.pi*k/nz_signal
            zfac = (1 - np.sin(kz)**2 + 2*np.sin(kz)**4/15)
            if power:
                deltak[i,j,k] /= xfac*yfac*zfac
            else:
                deltak[i,j,k] /= (xfac*yfac*zfac)**.5

            
@parallel_bcast([(nb.float32[:,:],nb.bool_), (nb.float64[:,:],nb.bool_)], '(nx,ny),()')
def tsc_window_aliased_2d(deltak, power, loop_idx):
    i = loop_idx[0]
    nx,ny = deltak.shape
    ny_signal = 2*ny - 2  # warning: assumes evenness of original dimensions
    kx = np.pi*i/nx
    xfac = (1 - np.sin(kx)**2 + 2*np.sin(kx)**4/15)
    for j in range(ny):
        ky = np.pi*j/ny_signal
        yfac = (1 - np.sin(ky)**2 + 2*np.sin(ky)**4/15)
        if power:
            deltak[i,j] /= xfac*yfac
        else:
            deltak[i,j] /= (xfac*yfac)**.5
            
def tsc_window(deltak, aliased=True, power=True):
    '''
    Selects the right window function to apply
    based on `aliased` and the shape of `deltak`.

    If `power` is True, then deltak is assumed
    to actually be deltak^2, and the square of
    the normalization will be applied.
    '''
    if deltak.ndim == 2:
        if aliased:
            tsc_window_aliased_2d(deltak, power)
        else:
            tsc_window_2d(deltak, power)
    elif deltak.ndim == 3:
        if aliased:
            tsc_window_aliased_3d(deltak, power)
        else:
            tsc_window_3d(deltak, power)
    else:
        raise ValueError(deltak.ndim, power)
