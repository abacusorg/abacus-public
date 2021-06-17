'''
Our numba implementations of histogramming of gridded quantities
are fast, but not very flexible, so there's a fair number of
repititous helper functions.  We store those here.

TODO: we can replace cyl and mu functios with np.prange versions

TODO: if we knew everything was cubic, we could just rescale the bin edges

TODO: we could definitely accelerate the counting/radius histogramming
via the monotonicity of the k vectors, as we do with the value histogramming
'''

import numpy as np
import scipy
from Abacus.Tools import numexpr_dist
import numpy.fft as fft

from .misc import parallel_bcast
import numba as nb

pi = np.pi

# Unweighted histogram
# Presently we only support a cubic domain, but it may be sampled at different spacings in each dimension
@nb.njit(parallel=True)
def hist_helper_3D(boxsize, gridshape, bin_edges, rfft):
    nx,ny,nz = gridshape
    Lx=Ly=Lz = boxsize

    # Lattice spacings
    if rfft:
        _dx = 2*pi/boxsize
        _dy = 2*pi/boxsize
        _dz = 2*pi/boxsize
    else:
        _dx = Lx/nx
        _dy = Ly/ny
        _dz = Lz/nz

    # Do binning with squared distances
    bin_edges2 = bin_edges**2
    nbp1 = len(bin_edges)
    bmin2,bmax2 = bin_edges2.min(), bin_edges2.max()

    # temporary thread workspace
    _hist = np.zeros((nx, nbp1-1), dtype=np.uint64)

    for i in nb.prange(nx):
        if i > nx//2:  # periodic wrap
            dx2 = ((nx-i)*_dx)**2
        else:
            dx2 = (i*_dx)**2
        for j in range(ny):
            if j > ny//2:
                dy2 = ((ny-j)*_dy)**2
            else:
                dy2 = (j*_dy)**2
            for k in range(nz):
                if not rfft and k > nz//2:
                    dz2 = ((nz-k)*_dz)**2
                else:
                    dz2 = (k*_dz)**2
                dist2 = dx2 + dy2 + dz2
                if dist2 < bmin2 or dist2 > bmax2:
                    continue
                for b in range(1,nbp1):
                    if dist2 < bin_edges2[b]:
                        break
                if k > 0 and rfft:
                    _hist[i,b-1] += 2  # double-count values that would be reflected in the full space
                else:
                    _hist[i,b-1] += 1
                #assert b-1 < nbp1-1

    # combine thread results
    histogram = np.sum(_hist, axis=0)
    return histogram

# Radius-weighted histogram
@nb.njit(parallel=True)
def rhist_helper_3D(boxsize, gridshape, bin_edges, rfft):
    nx,ny,nz = gridshape
    Lx=Ly=Lz = boxsize

    # Lattice spacings
    if rfft:
        _dx = 2*pi/boxsize
        _dy = 2*pi/boxsize
        _dz = 2*pi/boxsize
    else:
        _dx = Lx/nx
        _dy = Ly/ny
        _dz = Lz/nz

    # Do binning with squared distances
    bin_edges2 = bin_edges**2
    nbp1 = len(bin_edges)
    bmin2,bmax2 = bin_edges2.min(), bin_edges2.max()

    # temporary thread workspace
    _hist = np.zeros((nx, nbp1-1), dtype=np.float64)

    for i in nb.prange(nx):
        if i > nx//2:  # periodic wrap
            dx2 = ((nx-i)*_dx)**2
        else:
            dx2 = (i*_dx)**2
        for j in range(ny):
            if j > ny//2:
                dy2 = ((ny-j)*_dy)**2
            else:
                dy2 = (j*_dy)**2
            for k in range(nz):
                if not rfft and k > nz//2:
                    dz2 = ((nz-k)*_dz)**2
                else:
                    dz2 = (k*_dz)**2
                dist2 = dx2 + dy2 + dz2
                if dist2 < bmin2 or dist2 > bmax2:
                    continue
                for b in range(1,nbp1):
                    if dist2 < bin_edges2[b]:
                        break
                if k > 0 and rfft or not rfft:
                    _hist[i,b-1] += np.sqrt(dist2)
                else:
                    _hist[i,b-1] += 0.5*np.sqrt(dist2)
                #assert b-1 < nbp1-1

    # combine thread results
    histogram = np.sum(_hist, axis=0)

    # We half-counted everything since we're in the half-plane
    if rfft:
        histogram *= 2

    return histogram

# Weighted histogram
@nb.njit(parallel=True)
def whist_helper_3D(boxsize, values, bin_edges, rfft, multipoles=np.array([0])):
    if not rfft:
        raise NotImplementedError  # the if statement kills numba parallelization 
    nx,ny,nz = values.shape
    Lx=Ly=Lz = boxsize

    # Lattice spacings
    if rfft:
        _dx = 2*pi/boxsize
        _dy = 2*pi/boxsize
        _dz = 2*pi/boxsize
        nznyquist = nz
    else:
        _dx = Lx/nx
        _dy = Ly/ny
        _dz = Lz/nz
        nznyquist = nz//2
        

    # Do binning with squared distances
    bin_edges2 = bin_edges**2
    bmin2,bmax2 = bin_edges2.min(), bin_edges2.max()
    nbp1 = len(bin_edges2)

    # We will compute all multipoles in one pass
    # we will handle the monopole explicitly without mu_k calculation
    n_poles = len(multipoles)
    do_monopole = False
    if multipoles[0] == 0:  # guaranteed sorted
        do_monopole = True

    # "extra poles" are non-monopoles
    have_extra_poles = not do_monopole or n_poles > 1
    extra_pole_start = 1 if do_monopole else 0

    #Lmax = multipoles.max() if len(multipoles) > 0 else 0

    # temporary thread workspace
    _hist = np.zeros((nx, nbp1-1, n_poles), dtype=values.dtype)

    for i in nb.prange(nx):
        dx2 = ((i if i <= nx//2 else nx - i)*_dx)**2
        for j in range(ny):
            dy2 = ((j if j <= ny//2 else ny - j)*_dy)**2
            
            b = 1
            for k in range(nznyquist):
                dz2 = (k*_dz)**2
                dist2 = dx2 + dy2 + dz2

                if dist2 < bmin2:
                    continue
                if dist2 > bmax2:
                    break

                while b < nbp1-1 and dist2 >= bin_edges2[b]:
                    b += 1

                # first pole, probably monopole
                if do_monopole:
                    # We want to double-count values relative to the k=0 plane,
                    # as if we were binning in the whole space and not the 0<=kz<=kny half-space
                    if k == 0 and rfft:
                        _hist[i,b-1,0] += 0.5*values[i,j,k]
                    else:
                        _hist[i,b-1,0] += values[i,j,k]
                    #assert b-1 < nbp1-1

                # higher moments
                # this if statement introduces a performance penalty on the monopole...
                if have_extra_poles:
                    mu_k = np.sqrt(dz2/dist2)  # mu_k = k_z/|k|
                    for p in range(extra_pole_start, n_poles):
                        if k == 0 and rfft:
                            _hist[i,b-1,p] += 0.5*values[i,j,k]*legendre(multipoles[p], mu_k)
                        else:
                            _hist[i,b-1,p] += values[i,j,k]*legendre(multipoles[p], mu_k)

            # # do the second, descending half of the space
            # if not rfft:
            #     raise NotImplementedError  # Need to add poles and kz=0 plane awareness
            #     b = nbp1 - 1
            #     for k in range(nznyquist,nz):
            #         dz = (Lz/nz*(nz-k))
            #         dist2 = dx2 + dy2 + dz**2

            #         if dist2 < bmin2:
            #             break
            #         if dist2 > bmax2:
            #             continue

            #         while dist2 < bin_edges2[b] and b > 0:
            #             b -= 1
                    
            #         _hist[i,b,0] += values[i,j,k]
    
    # combine thread results
    histogram = np.sum(_hist, axis=0)

    # We counted in the half-space; pretend we counted in the full space
    histogram *= 2

    # apply multipoles prefactor
    histogram2 = histogram * (2*multipoles.reshape(-1) + 1)

    # Put multipoles first
    histogram3 = histogram2.T

    return histogram3

@nb.njit
def legendre(n, x):
    if n == 2:
        return 0.5*(3*x**2 - 1)
    elif n == 4:
        return 0.125*(35*x**4 - 30*x**2 + 3)
    elif n == 6:
        return 0.0625*(231*x**6 - 315*x**4 + 105*x**2 - 5)
    else:
        raise NotImplementedError


# Cylindrical histogram
@parallel_bcast([(nb.uint64[:,:], nb.float64[:], nb.uint64[:], nb.float64[:], nb.float64[:], nb.bool_)], '(nb,npib),(),(ndim),(nbp1),(npibp1),()')
def hist_helper_cyl(histogram, boxsize, gridshape, bin_edges, pi_bin_edges, rfft, loop_idx):
    boxsize = boxsize[0]
    i = loop_idx[0]
    nx,ny,nz = gridshape
    Lx=Ly=Lz = boxsize

    if rfft:
        Lx = nx*2*pi/boxsize
        Ly = ny*2*pi/boxsize
        Lz = nz*2*pi/boxsize

    # Do binning with squared distances
    bin_edges = bin_edges**2
    pi_bin_edges = pi_bin_edges**2
    nbp1 = len(bin_edges)
    npibp1 = len(pi_bin_edges)

    dx2 = (Lx/nx*(i if i <= nx//2 else nx - i))**2
    for j in range(ny):
        dy2 = (Ly/ny*(j if j <= ny//2 else ny - j))**2
        dist_perp = dx2 + dy2
        if dist_perp < bin_edges[0] or dist_perp > bin_edges[-1]:
            continue
        for bi in range(1,nbp1):
            if dist_perp < bin_edges[bi]:
                break
        for k in range(nz):
            dz2 = (Lz/nz*(k if k <= nz//2 or rfft else nz - k))**2
            if dz2 < pi_bin_edges[0] or dz2 > pi_bin_edges[-1]:
                continue
            for bj in range(1,npibp1):
                if dz2 < pi_bin_edges[bj]:
                    histogram[bi-1,bj-1] += 1
                    break
            else:  # last bin is closed
                histogram[bi-1,bj-1] += 1
                
# Radius-weighted cylindrical histogram
@parallel_bcast([(nb.float64[:,:], nb.float64[:], nb.uint64[:], nb.float64[:], nb.float64[:], nb.bool_)], '(nb,npib),(),(ndim),(nbp1),(npibp1),()')
def rhist_helper_cyl(histogram, boxsize, gridshape, bin_edges, pi_bin_edges, rfft, loop_idx):
    boxsize = boxsize[0]
    i = loop_idx[0]
    nx,ny,nz = gridshape
    Lx=Ly=Lz = boxsize

    if rfft:
        Lx = nx*2*pi/boxsize
        Ly = ny*2*pi/boxsize
        Lz = nz*2*pi/boxsize

    # Do binning with squared distances
    bin_edges = bin_edges**2
    pi_bin_edges = pi_bin_edges**2
    nbp1 = len(bin_edges)
    npibp1 = len(pi_bin_edges)

    dx2 = (Lx/nx*(i if i <= nx//2 else nx - i))**2
    for j in range(ny):
        dy2 = (Ly/ny*(j if j <= ny//2 else ny - j))**2
        dist_perp = dx2 + dy2
        if dist_perp < bin_edges[0] or dist_perp > bin_edges[-1]:
            continue
        for bi in range(1,nbp1):
            if dist_perp < bin_edges[bi]:
                break
        for k in range(nz):
            dz2 = (Lz/nz*(k if k <= nz//2 or rfft else nz - k))**2
            if dz2 < pi_bin_edges[0] or dz2 > pi_bin_edges[-1]:
                continue
            for bj in range(1,npibp1):
                if dz2 < pi_bin_edges[bj]:
                    histogram[bi-1,bj-1] += np.sqrt(dist_perp + dz2)
                    break
            else:  # last bin is closed
                histogram[bi-1,bj-1] += np.sqrt(dist_perp + dz2)

# Weighted cylindrical histogram
@parallel_bcast([(nb.float64[:,:], nb.float64[:], nb.float32[:,:,:], nb.float64[:], nb.float64[:],nb.bool_)], '(nb,npib),(),(nx,ny,nz),(nbp1),(npibp1),()')
def whist_helper_cyl(histogram, boxsize, values, bin_edges, pi_bin_edges, rfft, loop_idx):
    boxsize = boxsize[0]
    i = loop_idx[0]
    nx,ny,nz = values.shape
    Lx=Ly=Lz = boxsize

    if rfft:
        Lx = nx*2*pi/boxsize
        Ly = ny*2*pi/boxsize
        Lz = nz*2*pi/boxsize

    # Do binning with squared distances
    bin_edges = bin_edges**2
    pi_bin_edges = pi_bin_edges**2
    nbp1 = len(bin_edges)
    npibp1 = len(pi_bin_edges)

    dx2 = (Lx/nx*(i if i <= nx//2 else nx - i))**2
    for j in range(ny):
        dy2 = (Ly/ny*(j if j <= ny//2 else ny - j))**2
        dist_perp = dx2 + dy2
        if dist_perp < bin_edges[0] or dist_perp > bin_edges[-1]:
            continue
        for bi in range(1,nbp1):
            if dist_perp < bin_edges[bi]:
                break
        for k in range(nz):
            dz2 = (Lz/nz*(k if k <= nz//2 or rfft else nz - k))**2
            if dz2 < pi_bin_edges[0] or dz2 > pi_bin_edges[-1]:
                continue
            for bj in range(1,npibp1):
                if dz2 < pi_bin_edges[bj]:
                    histogram[bi-1,bj-1] += values[i,j,k]
                    break
            else:  # last bin is closed
                histogram[bi-1,bj-1] += values[i,j,k]


# Unweighted mu histogram
# Presently we only support a cubic domain, but it may be sampled at different spacings in each dimension
@parallel_bcast([(nb.uint64[:,:], nb.float64[:], nb.uint64[:], nb.float64[:], nb.bool_, nb.float64[:])], '(nb,nmub),(),(ndim),(nbp1),(),(nmubp1)')
def hist_helper_mu(histogram, boxsize, gridshape, bin_edges, rfft, mu_bin_edges, loop_idx):
    i = loop_idx[0]
    boxsize = boxsize[0]
    nx,ny,nz = gridshape
    Lx=Ly=Lz = boxsize

    # When binning an rfft grid in k-space, the lattice is 2*pi*i/L
    # Normally, the lattice is i*L/nx
    if rfft:
        Lx = nx*2*pi/boxsize
        Ly = ny*2*pi/boxsize
        Lz = nz*2*pi/boxsize

    # Do binning with squared distances
    bin_edges = bin_edges**2
    mu_bin_edges = mu_bin_edges**2
    nbp1 = len(bin_edges)
    nmubp1 = len(mu_bin_edges)

    if i > nx//2:  # periodic wrap
        i = nx - i
    dx2 = (i*Lx/nx)**2
    for j in range(ny):
        if j > ny//2:
            j = ny - j
        dy2 = (j*Ly/ny)**2
        for k in range(nz):
            if k > nz//2 and not rfft:
                k = nz - k
            dz2 = (k*Lz/nz)**2
            dist = dx2 + dy2 + dz2
            if dist < bin_edges[0] or dist > bin_edges[-1]:
                continue
            mu2 = dz2/dist
            if mu2 < mu_bin_edges[0] or mu2 > mu_bin_edges[-1]:
                continue
            for bi in range(1,nbp1):
                if dist < bin_edges[bi]:
                    break
            for bj in range(1,nmubp1):
                if mu2 < mu_bin_edges[bj]:
                    break
            histogram[bi-1,bj-1] += 1

# Radius-weighted mu histogram
@parallel_bcast([(nb.float64[:,:], nb.float64[:], nb.uint64[:], nb.float64[:], nb.bool_, nb.float64[:])], '(nb,nmub),(),(ndim),(nbp1),(),(nmubp1)')
def rhist_helper_mu(histogram, boxsize, gridshape, bin_edges, rfft, mu_bin_edges, loop_idx):
    boxsize = boxsize[0]
    i = loop_idx[0]
    nx,ny,nz = gridshape
    Lx=Ly=Lz = boxsize

    if rfft:
        Lx = nx*2*pi/boxsize
        Ly = ny*2*pi/boxsize
        Lz = nz*2*pi/boxsize


    # Do binning with squared distances
    bin_edges = bin_edges**2
    mu_bin_edges = mu_bin_edges**2
    nbp1 = len(bin_edges)
    nmubp1 = len(mu_bin_edges)

    if i > nx//2:  # periodic wrap
        i = nx - i
    dx2 = (i*Lx/nx)**2
    for j in range(ny):
        if j > ny//2:
            j = ny - j
        dy2 = (j*Ly/ny)**2
        for k in range(nz):
            if k > nz//2 and not rfft:
                k = nz - k
            dz2 = (k*Lz/nz)**2
            dist = dx2 + dy2 + dz2
            if dist < bin_edges[0] or dist > bin_edges[-1]:
                continue
            mu2 = dz2/dist
            if mu2 < mu_bin_edges[0] or mu2 > mu_bin_edges[-1]:
                continue
            for bi in range(1,nbp1):
                if dist < bin_edges[bi]:
                    break
            for bj in range(1,nmubp1):
                if mu2 < mu_bin_edges[bj]:
                    break
            histogram[bi-1,bj-1] += np.sqrt(dist)

# Weighted mu histogram
@parallel_bcast([(nb.float64[:,:], nb.float64[:], nb.float32[:,:,:], nb.float64[:], nb.bool_, nb.float64[:]),
                 (nb.float64[:,:], nb.float64[:], nb.float64[:,:,:], nb.float64[:], nb.bool_, nb.float64[:])], '(nb,nmub),(),(nx,ny,nz),(nbp1),(),(nmubp1)')
def whist_helper_mu(histogram, boxsize, values, bin_edges, rfft, mu_bin_edges, loop_idx):
    boxsize = boxsize[0]
    i = loop_idx[0]
    nx,ny,nz = values.shape
    Lx=Ly=Lz = boxsize

    if rfft:
        Lx = nx*2*pi/boxsize
        Ly = ny*2*pi/boxsize
        Lz = nz*2*pi/boxsize


    # Do binning with squared distances
    bin_edges = bin_edges**2
    mu_bin_edges = mu_bin_edges**2
    nbp1 = len(bin_edges)
    nmubp1 = len(mu_bin_edges)

    iwrap = i if i <= nx//2 else nx - i
    dx2 = (iwrap*Lx/nx)**2
    for j in range(ny):
        jwrap = j if j <= ny//2 else ny - j
        dy2 = (jwrap*Ly/ny)**2
        for k in range(nz):
            kwrap = k if k <= nz//2 or rfft else nz - k
            dz2 = (kwrap*Lz/nz)**2
            dist = dx2 + dy2 + dz2
            if dist < bin_edges[0] or dist > bin_edges[-1]:
                continue
            mu2 = dz2/dist
            if mu2 < mu_bin_edges[0] or mu2 > mu_bin_edges[-1]:
                continue
            for bi in range(1,nbp1):
                if dist < bin_edges[bi]:
                    break
            for bj in range(1,nmubp1):
                if mu2 < mu_bin_edges[bj]:
                    break
            histogram[bi-1,bj-1] += values[i,j,k]
