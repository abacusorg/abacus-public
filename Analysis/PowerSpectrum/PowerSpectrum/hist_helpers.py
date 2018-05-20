'''
Our numba implementations of histogramming of gridded quantities
are fast, but not very flexible, so there's a fair number of
repititous helper functions.  We store those here.
'''

import numpy as np
import scipy
from Abacus.Tools import numexpr_dist
import numpy.fft as fft

from .misc import parallel_bcast
import numba as nb

pi = np.pi

# Unweighted histogram
# Could add a 'rfft' kwarg (if supported by numba) to add a factor of pi and make the last axis non-periodic]
# Presently we only support a cubic domain, but it may be sampled at different spacings in each dimension
@parallel_bcast([(nb.uint64[:], nb.float64[:], nb.uint64[:], nb.float64[:], nb.bool_)], '(nb),(),(ndim),(nbp1),()')
def hist_helper_3D(histogram, boxsize, gridshape, bin_edges, rfft, loop_idx):
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
    nbp1 = len(bin_edges)

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
            for b in range(1,nbp1):
                if dist < bin_edges[b]:
                    histogram[b-1] += 1
                    break
            else:  # last bin is closed
                histogram[-1] += 1

# Radius-weighted histogram
@parallel_bcast([(nb.float64[:], nb.float64[:], nb.uint64[:], nb.float64[:], nb.bool_)], '(nb),(),(ndim),(nbp1),()')
def rhist_helper_3D(histogram, boxsize, gridshape, bin_edges, rfft, loop_idx):
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
    nbp1 = len(bin_edges)

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
            for b in range(1,nbp1):
                if dist < bin_edges[b]:
                    histogram[b-1] += np.sqrt(dist)
                    break
            else:  # last bin is closed
                histogram[-1] += np.sqrt(dist)

# Weighted histogram
@parallel_bcast([(nb.float64[:], nb.float64[:], nb.float32[:,:,:], nb.float64[:], nb.bool_),
                 (nb.float64[:], nb.float64[:], nb.float64[:,:,:], nb.float64[:], nb.bool_)], '(nb),(),(nx,ny,nz),(nbp1),()')
def whist_helper_3D(histogram, boxsize, values, bin_edges, rfft, loop_idx):
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
    nbp1 = len(bin_edges)

    dx2 = (Lx/nx*(i if i <= nx//2 else nx - i))**2
    for j in range(ny):
        dy2 = (Ly/ny*(j if j <= ny//2 else ny - j))**2
        for k in range(nz):
            dz2 = (Lz/nz*(k if k <= nz//2 or rfft else nz - k))**2
            dist = dx2 + dy2 + dz2
            if dist < bin_edges[0] or dist > bin_edges[-1]:
                continue
            for b in range(1,nbp1):
                if dist < bin_edges[b]:
                    histogram[b-1] += values[i,j,k]
                    break
            else:  # last bin is closed
                histogram[-1] += values[i,j,k]

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
            kwrap = k if k <= ny//2 or rfft else nz - k
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
