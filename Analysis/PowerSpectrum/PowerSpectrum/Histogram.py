import numpy as np
import scipy
from Abacus.Tools import numexpr_dist
import numpy.fft as fft

from .misc import parallel_bcast
import numba as nb
from .hist_helpers import *

pi = np.pi

def RadialBin(radii, values, bin_edges, bin_info=None):
    """
    Radial histogramming with weights (`values`) that supports
    fast re-evaluation for the same set of `radii` but new weights
    by passing back `bin_info`.
    
    Parameters
    ----------
    radii: ndarray
        The radii of the samples
    values: ndarray
        The values of the samples
    bin_edges: array_like, length nbins+1
        The radial bin edges
    bin_info: 3-tuple of ndarray, optional
        Tuple of `bin_indices`, `count`, and `bin_radii`.
        These will be used to accelerate the binning.  Default: None.
        
    Returns
    -------
    bin_radii: ndarray, shape (nbins,)
        The mean radius of the samples in the bin.
        Empty bins use the bin center.
    mean: ndarray, shape (nbins,)
        The mean value of the samples in each bin.
    count: ndarray, shape (nbins,)
        The number of samples in each bin
    bin_info: 3-tuple of ndarray
        Pass this back to `RadialBin` to accelerate the next binning.
        See `bin_info` parameter.
    """
    radii = radii.reshape(-1)
    values = values.reshape(-1)
    nbins = len(bin_edges) - 1
    
    if bin_info is not None:  # Only need to re-evaluate weights
        bin_indices, count, bin_radii = bin_info
        weight = np.bincount(bin_indices, weights=values, minlength=nbins+2)
        weight = weight[1:-1]
        
        zeros = count == 0
        with np.errstate(divide='ignore'):
            mean = weight/count
        mean[zeros] = 0
        
        bin_info = bin_indices, count, bin_radii
        return bin_radii, mean, count, bin_info
    
    # No bin_info was passed.  Do first (slow) binning.
    
    # Note: digitize seems to convert radii to double internally.
    # Could rewrite this in numba
    bin_indices = np.digitize(radii, bin_edges).astype(np.uint16)
    
    count = np.bincount(bin_indices, minlength=nbins+2)
    bin_radii = np.bincount(bin_indices, weights=radii, minlength=nbins+2)
    
    # The first element of `count` will be the number of occurences of 0
    # which is the number of elements below the first bin.  Discard this.
    # Even though values in the 1st bin are assigned to count[1], discarding
    # count[0] will shift these down so count[0] holds the first bin.
    # Values above the last bin will be assigned `nbins+1`.  Discard these.
    count = count[1:-1]
    bin_radii = bin_radii[1:-1]
    
    zeros = count == 0
    with np.errstate(divide='ignore'):
        bin_radii = bin_radii/count
    bin_centers = (bin_edges[:-1] + bin_edges[1:])/2.
    bin_radii[zeros] = bin_centers[zeros]
    
    bin_info = bin_indices, count, bin_radii
    return RadialBin(radii, values, bin_edges, bin_info=bin_info)


def RadialBinGrid(boxsize, values, bin_edges, pi_bin_edges=None, mu_bin_edges=None, bin_info=None, rfft=False, multipoles=0):
    '''
    Compute the weighted histogram of a set of points on a
    periodic lattice by radial binning.
    
    This is designed to be very memory efficient and never
    contructs temporary arrays.

    The periodicity is in the FFT convention of `x = i*dx`
    for `i < nx/2`, and `x = (nx - i)*dx` for `i >= nx/2`.

    This function uses numba internally and is parallelized.
    
    Parameters
    ----------
    boxsize: float or array_like of float
        Domain edge length(s)
    values: ndarray
        The weights to bin
    bin_edges: array_like of float
        The radial bin edges
    pi_bin_edges: array_like of float
        z-direction bin edges for 2D histogramming
    mu_bin_edges: array_like of float
        mu bin edges for 2D histogramming
    bin_info: 2-tuple of ndarray
        The `bin_counts` and `bin_radii`
    multipoles: int or array_like of int
        The list of  multipoles to return.  If given an iterable,
        the `bin_weight` return value will be a dict keyed by multipole.
        
    Returns
    -------
    bin_radii: ndarray, shape (nbins, [npibins])
        The mean radius of the samples in the bin.
        Empty bins use the bin center.
    bin_weight: ndarray, shape (nbins, [npibins])
        The mean value of the samples in each bin.
    bin_count: ndarray, shape (nbins, [npibins])
        The number of samples in each bin
    bin_info: 2-tuple of ndarray
        Pass this back to `RadialBinGrid` to accelerate the next binning.
        See `bin_info` parameter.
    '''
    assert not (pi_bin_edges is not None and mu_bin_edges is not None), "Can't specify both pi and mu bins!"

    # 2D grids will be binned as a 3D NxMx1 grid
    if values.ndim == 2:
        assert pi_bin_edges is None and mu_bin_edges is None
        values = values[...,None]
    
    gridshape = np.array(values.shape, dtype=np.uint64)

    try:
        len(multipoles)
        flat_multipoles = False
    except TypeError:
        flat_multipoles = True
    multipoles = np.atleast_1d(multipoles)
    assert (np.diff(multipoles) > 0).all()
    
    hist_helper = hist_helper_3D
    rhist_helper = rhist_helper_3D
    whist_helper = whist_helper_3D
    twoD_bin_arg = []
    
    assert values.ndim == 3
    assert (np.diff(bin_edges) > 0).all()
    if pi_bin_edges is not None:
        assert (np.diff(pi_bin_edges) > 0).all()
        
        hist_helper = hist_helper_cyl
        rhist_helper = rhist_helper_cyl
        whist_helper = whist_helper_cyl
        twoD_bin_arg = [pi_bin_edges]

    if mu_bin_edges is not None:
        assert (np.diff(mu_bin_edges) > 0).all()
        
        hist_helper = hist_helper_mu
        rhist_helper = rhist_helper_mu
        whist_helper = whist_helper_mu
        twoD_bin_arg = [mu_bin_edges]
    
    if bin_info is None:
        counts = hist_helper(boxsize, gridshape, bin_edges, rfft, *twoD_bin_arg)
    else:
        counts = bin_info[0]
    
    if bin_info is None:
        radii = rhist_helper(boxsize, gridshape, bin_edges, rfft, *twoD_bin_arg)
        radii /= counts
    else:
        radii = bin_info[1]
    
    weights = whist_helper(boxsize, values, bin_edges, rfft,*twoD_bin_arg, multipoles=multipoles)
    weights /= counts

    if flat_multipoles:
        weights = weights[0]
    else:
        weights = {multipoles[i]:w for (i,w) in enumerate(weights)}
    
    bin_info = (counts, radii)
    
    return radii, weights, counts, bin_info
    

def distance_array(gridshape, box, cell_centered=False, periodic=True, dtype=np.float32):
    """
    Given a grid shape, get the radial distance from the origin
    to each lattice site or grid cell center.  Distances can be periodic wrapped.
    
    This method uses numexpr and is thus extremely efficient.
    
    Parameters
    ----------
    gridshape: array_like of int
        The number of cells in each dimension
    box: float
        The periodic box length
    cell_centered: bool, optional
        Use distances to the centers of cells instead of the lattice sites.
        Default: False
    periodic: bool, optional
        Do a periodic wrap of cell distances along each dimension.
        Default: True
    dtype: np.dtype, optional
        Data type of the result.  Default: np.float32
        
    Returns
    -------
    dist: ndarray of float
        The array of distances to cells, of shape (num_xcells, num_ycells, ...).
    
    """
    gridshape = np.atleast_1d(gridshape)
    steps = box/gridshape
    
    cc = np.ogrid[tuple(slice(0.,box,dx) for dx in steps)]
    cc = [c.astype(dtype) for c in cc]
    if cell_centered:
        cc = [c + 0.5*dx for c,dx in zip(cc,steps)]
    if periodic:
        cc = [c - box*np.round(c/box) for c in cc]
    #dist = sum(c**2. for c in cc)**.5
    dist = numexpr_dist(cc)
    
    assert dist.dtype == dtype, \
            "distance_array does not have the expected dtype! Found '" + str(dist.dtype) + \
            "'; expected '" + str(dtype) + "'"
    
    return dist


def k_bin_edges(gridshape, boxsize, nbins=-1, log=False, bin_like_nfft=0):
    '''
    Produce radial bins in k-space for a given mesh shape.
    Produces spherical bins using the largest bin dimension.
    Bins above Nyquist will thus be sampling corner modes,
    in the case of a cubic grid.
    
    Parameters
    ----------
    gridshape: int or tuple
        Shape of the input signal, not the RFFT mesh
    boxsize: float
        Domain size
    nbins: int, optional
        Number of bins.  Default of -1 uses max(gridshape)/4.
    log: bool, optional
        Bin in log-space.  Default: False.
    '''
    if type(gridshape) is int:
        gridshape = (gridshape,)*3
       
    if bin_like_nfft:
        ngrid = bin_like_nfft
    else:
        ngrid = min(gridshape)  # use min so the k-spheres always fit
    
    kmax = np.pi/(boxsize/ngrid)  # nyquist
    
    if nbins == -1:
        nbins = ngrid//4
        
    if log:
        # Set up uniform bins in log space
        # Histogram.py and logspace are both inclusive on the upper bound of the last bin
        kmin = 2*np.pi/boxsize  # fundamental mode
        bin_edges = np.logspace(np.log10(kmin), np.log10(kmax), num=nbins+1)
    else:
        kmin = 0.  # ensures equal spacing for different ngrid
        bin_edges = np.linspace(kmin, kmax, num=nbins+1)

    # Now clip bins to the real ngrid
    if bin_like_nfft:
        real_kmax = np.pi/(boxsize/min(gridshape))
        bin_edges = bin_edges[bin_edges <= real_kmax]

    return bin_edges
