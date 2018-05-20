"""
A library for computing covariance matricies from power spectra;
e.g. computing the covariance of a two-point correlation function,
given the power spectrum that produced it.

Author: Lehman Garrison

"""

import scipy
import numpy as np
from Abacus import Tools
from timeit import default_timer as timer
from . import PowerSpectrum
from . import misc
import numexpr as ne

import pyfftw
from . import Histogram

def Cov_2PCF(ps, boxsize, bin_edges, pi_bin_edges=None, inplace=False, symmetrize=True, shot_noise_np=None, n_cat=None):
    """
    Given a 2D or 3D power spectrum, compute the covariance of the
    binned 2PCF.  The covariance w.r.t. unbinned separations is written as:
    
    C_12 = <Xi(s_1)Xi(s_2)> - <Xi(s_1)><Xi(s_2)>
    
    which is then binned in spherical shells defined by `bin_edges` to get
    C_ab, the covariance w.r.t. bin 'a' and bin 'b'.  We assume Gaussianity
    (independent k-modes) and periodicity (i.e. simulation box).
    
    If `pi_bin_edges` is given, the binning is in cylindrical annuli, with the
    line-of-sight bin edges (z-axis, or r_pi-direction) given by `pi_bin_edges`,
    and the r_perp (projected separation) bins given by `bin_edges`.
    
    If `ps` is 2D, then the binning is in circular annuli, and the Fourier
    transforms are 2D (i.e. computationally cheaper).
    
    Parameters
    ----------
    ps: ndarray of shape (nx,ny,nz/2+1) or (nx,ny/2+1)
        The input 2D or 3D power spectrum (i.e. the square of the FFT of the
        density field).  This should be normalized to the
        standard cosmology PS definition, meaning FFTW or np.fft outputs (squared)
        should be multiplied by V/N^2, where N=nx*ny*nz.
        The Analysis.PowerSpectrum outputs are already normalized.
        nx, ny, and nz should be even.
    boxsize: float
        The physical domain edge length
    bin_edges: array_like of length nbins+1
        The bin edges for either the 3D separations or projected separations,
        depending on the binning mode (see above).
    pi_bin_edges: array_like of length npibins+1, optional
        The bin edges for line-of-sight (z-axis) separations r_pi.
    inplace: bool, optional
        Allow overwriting of `ps`
    symmetrize: bool, optional
        Enforce symmetry of the returned covariance matrix.
        Default: True.
    shot_noise_np: int, optional
        The number of particles (or halos) from which the input power spectrum `ps`
        was measured.  If not None, then the nominal shot noise level
        `P_shot(k) = 1/shot_noise_np` will be subtracted from the input PS.
        The compensatory term will added to the covariance analytically.
        Default: None.
    n_cat: int, optional
        The number of catalogs (realizations) from which the power spectrum
        was measured.  Because the power spectrum is a noisy quantity, squaring
        it gives a result biased by (n_cat + 2)/n_cat.  Providing this parameter
        de-biases the P^2 estimator.  Default: None.
        
    Returns
    -------
    cov: ndarray of shape (ntotbins,ntotbins)
        The 2PCF covariance matrix.
        If both `bin_edges` and `pi_bin_edges` are given, ntotbins = nbins*npibins.
        The i-th perp bin and j-th pi bin is indexed as "i*npibins + j"
    """
    
    # Set up the bins and cov matrix
    bin_edges = np.asarray(bin_edges)
    nbins = len(bin_edges) - 1
    ntotbins = nbins
    if pi_bin_edges is not None:
        pi_bin_edges = np.asarray(pi_bin_edges)
        npibins = len(pi_bin_edges) - 1
        ntotbins = nbins*npibins
    dtype = ps.dtype.type
    Cab = np.zeros((ntotbins, ntotbins), dtype=dtype)
    pi = dtype(np.pi)
    
    # square ps and rename
    ps2_bias = (n_cat + 2.)/n_cat if n_cat else 1.;  ps2_bias = dtype(ps2_bias)
    if inplace:
        ne.evaluate('ps**2/ps2_bias', out=ps)
        ps2 = ps
    else:
        ps2 = ne.evaluate('ps**2/ps2_bias')
    del ps
    
    spherical = pi_bin_edges is None and ps2.ndim == 3
    
    # useful geometric quantities
    if ps2.ndim == 2:
        assert pi_bin_edges is None, "Cannot bin z direction with 2D power spectrum!"
        boxvol = dtype(boxsize**2.)
        bin_areas = np.pi*(bin_edges[1:]**2. - bin_edges[:-1]**2.)
    elif ps2.ndim == 3:
        boxvol = dtype(boxsize**3.)
        
    # subtract shot noise
    if shot_noise_np:
        ne.evaluate('ps2 - (boxvol/shot_noise_np)**2', out=ps2)
    
    # integral prefactor
    ne.evaluate('ps2 * 2/boxvol', out=ps2)
    
    # construct some grid quantities
    gridshape = np.asarray(ps2.shape)
    gridshape[-1] = (gridshape[-1] - 1)*2  # shape of the original binning, only correct for even nz
    assert (gridshape % 2 == 0).all(), 'To remove ifft ambiguity, FFT dimensions must be even'
    
    if spherical:        
        bin_vols = (4./3*np.pi*(bin_edges[1:]**3 - bin_edges[:-1]**3)).astype(dtype)
        
        bin_info, sorted_idx = None, None
        
        for i in range(nbins):
            # Apply the first binning in Fourier space with the j1 factor
            t0 = timer()
            Cab_raw = pyfftw.empty_aligned(ps2.shape, dtype=ps2.dtype)
            misc.shell_fft(Cab_raw, boxsize, bin_edges[i+1], bin_edges[i])
            ne.evaluate('Cab_raw*ps2', out=Cab_raw)
            print("j1 time: {:.2f}s".format(timer() - t0))

            t0 = timer()
            Cab_raw = PowerSpectrum.InverseFFTAndBin(Cab_raw, boxsize, bins=None, inplace=True)
            print("IRFFT time: {:.2f}s".format(timer() - t0))

            # Apply the second binning in real space
            t0 = timer()
            bin_radii, mean_in_bins, count_in_bins, bin_info = Histogram.RadialBinGrid(boxsize, Cab_raw, bin_edges, bin_info=bin_info)
            if (count_in_bins == 0).any():
                print('Warning: no covariance samples in at least one radial bin!  Try coarser bins or a finer mesh?')
            Cab[:,i] = mean_in_bins
            print("histogram time: {:.2f}s".format(timer() - t0))
            
    else:  # Either 3D or 2D cylindrical
        threeD = ps2.ndim == 3
        
        bin_areas = (pi*(bin_edges[1:]**2. - bin_edges[:-1]**2.)).astype(dtype)
        
        if threeD:  # 3D cylindrical
            pi_bin_widths = 2*(pi_bin_edges[1:] - pi_bin_edges[:-1])
            bin_vols = (bin_areas[:,None]*np.atleast_2d(pi_bin_widths)).reshape(-1)
            kpar = misc.fftfreq_nd(gridshape[[-1]], boxsize/gridshape[-1], rfft=True, kmag=True).reshape(1,1,-1)
            
            bin_info = None
            for i in range(nbins):
                # The cylindrical binning (with J1) is applicable for all pi bins
                Cab_perp = pyfftw.empty_aligned(ps2.shape[:-1], dtype=ps2.dtype)
                misc.annulus_fft(Cab_perp, boxsize, bin_edges[i+1], bin_edges[i])
                Cab_perp = Cab_perp[...,None]
                
                for j in range(npibins):
                    # First pi binning (with the sin factor) is in Fourier space
                    piblow,pibhigh = pi_bin_edges[j:j+2].astype(dtype)
                    pi_bin_width = pi_bin_widths[j].astype(dtype)
                    Cab_raw = ne.evaluate('2/kpar * (sin(pibhigh*kpar) - sin(piblow*kpar)) * Cab_perp * ps2 / pi_bin_width')
                    assert Cab_raw.dtype == dtype
                    Cab_raw[...,0] = ps2[...,0]  # fix divide by 0
                    
                    Cab_raw = PowerSpectrum.InverseFFTAndBin(Cab_raw, boxsize, bins=None, inplace=True)

                    # Apply the second binning in real space
                    bin_radii, mean_in_bins, count_in_bins, bin_info = Histogram.RadialBinGrid(boxsize, Cab_raw, bin_edges, pi_bin_edges=pi_bin_edges, bin_info=bin_info)
                    if (count_in_bins == 0).any():
                        print('Warning: no covariance samples in at least one radial bin!  Try coarser bins or a finer mesh?')
                    Cab[:,i*npibins + j] = mean_in_bins.reshape(-1)
            
        else: # 2D annular
            del kpar  # don't need
            
            paired_edges = np.stack((bin_edges[:-1], bin_edges[1:]), axis=-1)
            bin_info = None
            for i in range(nbins):
                # Apply the first binning in Fourier space with the J1 factor
                Cab_raw = pyfftw.empty_aligned(ps2.shape, dtype=ps2.dtype)
                misc.annulus_fft(Cab_raw, boxsize, bin_edges[i+1], bin_edges[i])
                ne.evaluate('Cab_raw*ps2', out=Cab_raw)
                
                Cab_raw = PowerSpectrum.InverseFFTAndBin(Cab_raw, boxsize, bins=None, inplace=True)
                
                # Apply the second binning in real space
                #mean_in_bins, edges, inds = scipy.stats.binned_statistic(dist.reshape(-1), Cab_raw.reshape(-1), bins=bin_edges)
                bin_radii, mean_in_bins, count_in_bins, bin_info = Histogram.RadialBinGrid(boxsize, Cab_raw, bin_edges, bin_info=bin_info)
                if (count_in_bins == 0).any():
                    print('Warning: no covariance samples in at least one radial bin!  Try coarser bins or a finer mesh?')
                Cab[:,i] = mean_in_bins
                
            bin_vols = bin_areas  # only for shot noise

    # Restore subtracted shot noise
    if shot_noise_np:
        shot_noise_diag = np.diag([(boxvol/shot_noise_np)**2.]*len(Cab))
        shot_noise_diag *= 2./boxvol / bin_vols
        Cab += shot_noise_diag
        
    # In most cases, we want to enforce symmetry
    Cab_sym = (Cab + Cab.T)/2.
    frac_diff = np.abs(Cab - Cab.T)/(Cab_sym**2).mean()**.5
    print("max frac diff: ", np.max(frac_diff))
    if (frac_diff > .05).any():
        print('Warning: unusually asymmetric covariance (max frac_diff {:g})'.format(np.max(frac_diff)))
    
    if symmetrize:
        return Cab_sym
    return Cab
    