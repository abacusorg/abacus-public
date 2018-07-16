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


# via https://docs.scipy.org/doc/numpy/user/basics.subclassing.html
class TransferMultipoles(np.ndarray):
    '''
    A simple wrapper around ndarray that represents a multipole transfer
    function.  Its main purpose is to carry around metadata that
    represents the binning scheme this function is valid for.
    '''
    def __new__(subtype, bins, box, Lmax=0, kmax=None, ngrid=None, zspace=False, dtype=np.float32):
        if type(bins) is int:
            nbins = bins
            if kmax is None:
                kmax = np.pi/(box/ngrid)
            k_bin_edges = np.linspace(2*np.pi/box, kmax, nbins+1)
        else:
            k_bin_edges = bins
            nbins = len(bins)-1
        
        obj = super(TransferMultipoles, subtype).__new__(subtype, (Lmax+1,nbins), dtype=dtype)
        
        obj.k_bin_edges = k_bin_edges
        obj.nbins = nbins
        obj.box = box
        obj.Lmax = Lmax
        obj.zspace = zspace
        obj[:] = 0.
        
        return obj
    
    def __array_finalize__(self, obj):
        if obj is None:
            return
        
        # Seems strange that we have to duplicate this functionality...
        self.k_bin_edges = obj.k_bin_edges
        self.nbins = obj.nbins
        self.box = obj.box
        self.Lmax = obj.Lmax
        self.zspace = obj.zspace


def apply_transfer_multipoles(Tfer, deltak, unbinned=1.):
    """
    Apply transfer function T_multipoles on deltak.

    Uses a "lookup" algorithm and parallelization over k-planes
    instead of unraveling the density cube as before.

    Parameters
    ----------
    T_multipoles: TransferMultipoles
        A TransferMultipoles object describing the transfer function to apply

    deltak: complex ndarray of shape (ngrid,ngrid,ngrid//2+1[,...])
        The delta(k) field on which to apply the transfer function.
        The first three dimensions are the FFT axes; any axes after
        that are broadcast over.

    unbinned: float, optional
        The transfer value to apply to modes that fall outside any k-bin.
        Typical values would be 1 (don't touch them) or 0 (null them out).
        Default: 1.

    Returns
    -------
    deltak_trans: ndarray like `deltak`
        The delta(k) field with the transfer function applied.
    """
    # Parse args
    # flatten all extra axes
    oshape = deltak.shape
    shape = oshape[:3] + (-1,)
    deltak = deltak.reshape(shape)

    if unbinned == 0.:
        # interestingly, zeros allocates virtual memory, while zeros_like touches the memory
        # we strongly prefer to let each thread have the first touch
        deltak_trans = np.zeros(deltak.shape, dtype=deltak.dtype)
    else:
        deltak_trans = deltak * unbinned

    knorm_bin_edges = Tfer.k_bin_edges * Tfer.box
    _apply_transfer_multipoles(Tfer, knorm_bin_edges, deltak, deltak_trans)

    # unflatten the last axes
    deltak_trans = deltak_trans.reshape(oshape)

    return deltak_trans

from .hist_helpers import legendre

@nb.njit(parallel=True)
def _apply_transfer_multipoles(T_multipoles, knorm_bin_edges, deltak, deltak_trans):
    Nx,Ny,Nz,Ndim = deltak.shape
    Lmaxp1,nbins = T_multipoles.shape

    # use squared edges in units of the fundamental
    k_bin_edges = (knorm_bin_edges/(2*np.pi))**2
    kmin = k_bin_edges.min()
    kmax = k_bin_edges.max()
    assert (np.diff(k_bin_edges) >= 0).all()

    for ix in nb.prange(Nx):
        if ix > Nx/2:
            kx2 = (Nx-ix)**2
        else:
            kx2 = ix**2

        for iy in range(Ny):
            if iy > Ny/2:
                ky2 = (Ny-iy)**2
            else:
                ky2 = iy**2

            t = 0
            for iz in range(Nz):
                kz2 = iz**2
                k2 = kx2 + ky2 + kz2
                mu = (kz2/k2)**0.5
                
                # is this k in some bin?
                if k2 < kmin:
                    continue
                if k2 > kmax:
                    # k will only increase!
                    break

                # look for this kmag bin
                # assumes monotonicity
                # last bin edge is closed, and we know this k is in some bin
                # note that t starts from the spot it stopped last time
                while k_bin_edges[t+1] <= k2 and t < nbins - 1:
                    t += 1

                for d in range(Ndim):
                    deltak_trans[ix,iy,iz,d] = deltak[ix,iy,iz,d] * T_multipoles[0,t]
                    for ell in range(2,Lmaxp1,2):  # only process even multipoles
                        deltak_trans[ix,iy,iz,d] += deltak[ix,iy,iz,d] * T_multipoles[ell,t] * legendre(ell, mu)
