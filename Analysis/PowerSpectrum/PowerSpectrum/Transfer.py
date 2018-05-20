"""
Compute "transfer functions" (i.e. ratios of power or deltak)
between two sets of particles or density fields.

Typically one is interested in some smooth basis in |k| or
(|k|,mu), so these routines generate smooth fitting
functions.
"""

from . import Histogram, PowerSpectrum
from .misc import fftfreq_nd

import numpy as np
import numpy.polynomial.legendre as npleg

# TODO: could make a decorator that can take pos, delta, or deltak
def AnisotropicTransferMultipoles(deltak1, deltak2, box, nbins=-1, Lmax=4, dtype=np.float32):
    """
    Computes the best-fit T_0, T_2, T_4 to the two deltaks, such that
    in each |k|-bin

    \chi^2 = \sum_i \delta_{2,i}(\mathbf{k}) - delta_{1,i}(\mathbf{k}) (T_0 L_0(\mu_i) + T_2 L_2(\mu_i) + T_4 L_4(\mu_i))

    is minimized.

    Parameters
    ----------
    deltak1, deltak2: ndarrays of shape (ngrid,ngrid,ngrid//2+1)
        The FFT of the normalized delta density field arrays.  `deltak2` is
        the "truth" that `deltak1` will try to fit.

    nbins: int, optional
        Number of linear |k| bins in which to compute T.
        Default of -1 means `ngrid`//4.

    box: float
        The domain size.

    Lmax: int, optional
        The maximum multipole order.  Only even multipoles will
        be computed.

    Returns
    -------
    T: ndarray of shape (`Lmax`+1, `nbins`)
        The transfer function multipole coefficients.
        Only the even multipoles will be non-zero.
        This can be passed to `apply_transfer_function_multipoles`
        to evaluate the Legendre series in each bin.

    korder: ndarray of int of shape (deltak1.size,)
        The indices that sort the `deltak`s into |k|-order.
        Helpful for applying T.
    """
    # Parse args
    shape = np.array(deltak1.shape)
    shape[-1] = 2*shape[-1] - 2  # assumes evenness of original dimensions
    box = float(box)
    if Lmax > 4:
        raise ValueError(Lmax)

    # Build the k-lattice distances and vectors
    kgrid = fftfreq_nd(shape, box/shape, rfft=True, kmag=True, dtype=dtype)
    kvec = fftfreq_nd(shape, box/shape, rfft=True, kmag=False, mesh=False, dtype=dtype)

    mu_k = kvec[-1]/kgrid  # k_z/|k|
    mu_k[0,0,0] = 0.  # this shouldn't matter, since our first k bin doesn't touch k=0

    korder = kgrid.argsort(axis=None)
    korder_kgrid = kgrid.flat[korder]
    korder_mu_k = mu_k.flat[korder]

    # Now every k bin is a contiguous segment of kgrid
    k_bin_edges = Histogram.k_bin_edges(shape, box, nbins=nbins)
    korder_bin_splits = np.searchsorted(korder_kgrid, k_bin_edges)

    korder_deltak1 = deltak1.flat[korder]
    korder_deltak2 = deltak2.flat[korder]

    del kgrid, kvec, mu_k, korder_kgrid

    # Compute the least-squares fit to the anisotropic multipoles
    all_T = []
    for i in range(len(k_bin_edges)-1):
        start = korder_bin_splits[i]
        stop = korder_bin_splits[i+1]
        
        this_delta_pred = korder_deltak1[start:stop]
        this_delta_true = korder_deltak2[start:stop]
        this_mu_k = korder_mu_k[start:stop]
        
        # X is NxL, where each observation is deltak_i*L_j(mu_i)
        X = np.zeros((stop-start,Lmax//2+1), dtype=this_delta_pred.dtype)  # 3 L's
        X[:,0] = npleg.legval(this_mu_k, [1])
        if Lmax >= 2:
            X[:,1] = npleg.legval(this_mu_k, [0,0,1])
        if Lmax >= 4:
            X[:,2] = npleg.legval(this_mu_k, [0,0,0,0,1])
        X *= this_delta_pred.reshape(-1,1)
        
        XH = X.conj().T
        cov = np.matmul(XH, this_delta_true).real
        var = np.linalg.inv(np.matmul(XH,X).real)
        T = np.matmul(var,cov)
        all_T += [T]
        
    all_T = np.array(all_T)
    # Pack L=0,2,4 into columns 0,2,4
    all_T = np.stack([ all_T[:,l//2] if l % 2 == 0 else np.zeros(len(all_T)) for l in range(Lmax+1) ])

    return all_T, korder

import numba as nb
@nb.njit(parallel=True)
def _reorder(arr, neworder):
    '''
    A fast, parallelized function to reorder an array according
    to a new index array.

    All of the following are equivalent:
    >>> newarr = arr[neworder]
    >>> newarr = arr.take(neworder, axis=0)
    >>> newarr = _reorder(arr, neworder)
    '''
    N = len(neworder)
    newarr = np.empty_like(arr[:N])

    for i in nb.prange(N):
        newarr[i] = arr[neworder[i]]
    return newarr
    
def apply_transfer_multipoles(T_multipoles, deltak, box, korder_info=None, unbinned=1.):
    """
    Apply transfer function T_multipoles on deltak.

    Parameters
    ----------
    T_multipoles: ndarray of shape (Lmax+1, nbins)
        The anisotropic multipole coefficients in each bin

    deltak: complex ndarray of shape (ngrid,ngrid,ngrid//2+1[,...])
        The delta(k) field on which to apply the transfer function.
        The first three dimensions are the FFT axes; any axes after
        that are broadcast over.

    box: float
        Box size

    korder: ndarray of int of shape (deltak.size,), optional
        The indices that sort the flattened deltak array into
        |k| order.  Will be computed if not given.

    unbinned: float, optional
        The transfer value to apply to modes that fall outside any k-bin.
        Typical values would be 1 (don't touch them) or 0 (null them out).
        Default: 1.

    Returns
    -------
    deltak_trans: ndarray like `deltak`
        The delta(k) field with the transfer function applied.
    
    korder_info: tuple of ndarray
        A tuple of `(korder, korder_kgrid, korder_mu_k, invkorder)` that can
        be passed back to this function to accelerate the next
        sorting.  Also accepts just `korder`.
    """
    # Parse args
    nbins = T_multipoles.shape[-1]

    orig_shape = deltak.shape
    dispdim_len = 1
    if deltak.ndim == 4:
        dispdim_len = deltak.shape[-1]

    kshape = np.array(deltak.shape[:3])
    shape = kshape.copy()
    shape[-1] = 2*shape[-1] - 2  # assumes evenness of original dimensions

    # tuple of all the reordered arrays
    if type(korder_info) in (list, tuple):
        korder, korder_kgrid, korder_mu_k, invkorder = korder_info
    else:
        # Build the k-lattice distances and vectors
        kgrid = fftfreq_nd(shape, box/shape, rfft=True, kmag=True)
        kvec = fftfreq_nd(shape, box/shape, rfft=True, kmag=False, mesh=False)

        if korder_info is not None:
            korder = korder_info
        else:
            # could use histogram then argpartition
            korder = kgrid.argsort(axis=None)
        
        invkorder = korder.argsort()
        mu_k = kvec[-1]/kgrid  # k_z/|k|
        mu_k[0,0,0] = 0.  # this shouldn't matter, since our first k bin doesn't touch k=0

        korder_kgrid = _reorder(kgrid.reshape(-1), korder)
        korder_mu_k = _reorder(mu_k.reshape(-1), korder)
        del kgrid, kvec, mu_k

    # Now every k bin is a contiguous segment of kgrid
    k_bin_edges = Histogram.k_bin_edges(shape, box, nbins=nbins)
    korder_bin_splits = np.searchsorted(korder_kgrid, k_bin_edges)

    korder_deltak = _reorder(deltak.reshape(-1,dispdim_len), korder)
    if unbinned == 0.:
        korder_deltak_trans = np.zeros_like(korder_deltak)
    else:
        korder_deltak_trans = korder_deltak * unbinned
    del deltak

    for i in range(len(korder_bin_splits)-1):
        start = korder_bin_splits[i]
        stop = korder_bin_splits[i+1]
        #this_deltak = korder_deltak[start:stop]
        this_mu_k = korder_mu_k[start:stop]
        this_T_multipoles = T_multipoles[:,i]
        
        transfer = npleg.legval(this_mu_k, this_T_multipoles)
        korder_deltak_trans[start:stop] = korder_deltak[start:stop]*transfer[...,None]

    del korder_deltak
    # Pack the flattened array back into a cube
    deltak_trans = _reorder(korder_deltak_trans, invkorder)
    del korder_deltak_trans
    deltak_trans = deltak_trans.reshape(orig_shape)

    return deltak_trans, (korder, korder_kgrid, korder_mu_k, invkorder)
