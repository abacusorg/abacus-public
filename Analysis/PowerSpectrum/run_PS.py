#!/usr/bin/env python3
# Copyright 2012-2025 The Abacus Developers
# SPDX-License-Identifier: GPL-3.0-or-later


'''
A user-friendly command-line interface to PowerSpectrum.

Usage: $ ./run_PS.py --help

Author: Lehman Garrison
'''

import sys
import os
from os.path import join as pjoin, abspath, basename, dirname, getsize
import pdb
import shutil
import argparse
from glob import glob

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from astropy.utils import iers
iers.conf.auto_download = False
import astropy.table

from Abacus.Analysis import common
from Abacus.InputFile import InputFile
from Abacus import Tools

from .PowerSpectrum import TSC, Histogram
from .PowerSpectrum import PowerSpectrum as PS

# A label for the power spectrum based on the properties
def ps_suffix(**kwargs):
    nfft = kwargs['nfft']
    try:
        len(nfft)
    except TypeError:  # convert a plain int to a singlet array
        nfft = [nfft]
        
    # if 2D, add a projected tag
    projected = len(nfft) == 2
    
    # if all dimensions the same, collapse to single int
    if (np.diff(nfft) == 0).all():
        nfft = nfft[0]
    
    ps_fn = ''
    if not kwargs['just_density']:
        if kwargs['track'] is not None:
            ps_fn += '_track{}'.format(kwargs['track'])
        ps_fn += '_nfft{}'.format(kwargs['nfft'])
    else:
        ps_fn += '_ngrid{}'.format(kwargs['nfft'])
    if kwargs['dtype'] != np.float32:
        ps_fn += '_{}bit'.format(np.dtype(kwargs['dtype']).itemsize*8)
    if kwargs['zspace']:
        ps_fn += '_zspace'
    if kwargs['rotate_to']:
        ps_fn += '_rotate({2.1g},{2.1g},{2.1g})'.format(*rotate_to)
    if projected:
        ps_fn += '_projected'
        
    if not kwargs['just_density']:
        ps_fn += '.csv'
    return ps_fn


# If the input is a PS,
# we put the result in the same directory
def get_output_path_ps(input):
    return os.path.dirname(input)

    
def run_PS_on_dir(slicedir, **kwargs):
    # Make a descriptive filename
    this_suffix = ps_suffix(**kwargs)
    kwargs.pop('track')
    
    just_density = kwargs.pop('just_density')
    product_name = 'power' if not just_density else 'density'
    out_parent = kwargs.pop('out_parent')
    if not (outdir := kwargs.pop('outdir')):
        outdir = common.get_output_dir(product_name, slicedir, out_parent=out_parent, **kwargs)
    ps_fn = product_name + this_suffix

    header, header_fn = common.get_header(slicedir, retfn=True)
    
    # box kwargs overrides header
    argbox = kwargs.pop('box')
    if argbox:
        BoxSize = argbox
    else:
        BoxSize = header['BoxSize']
            
    # Make the output dir and store the header
    os.makedirs(outdir, exist_ok=True)
    if header_fn:
        shutil.copy(header_fn, pjoin(outdir, 'header'))

    # Copy the sim-level parameter files
    if not os.path.isdir(pjoin(outdir, os.pardir, 'info')):
        try:
            shutil.copytree(pjoin(slicedir, os.pardir, 'info'), pjoin(outdir, os.pardir, 'info'))
        except:
            pass
    
    print('* Starting {} on {}'.format('PS' if not just_density else 'density', slicedir))
    save_fn = os.path.join(outdir, ps_fn)
    print(f'* and saving to {save_fn}')
    
    nfft = kwargs.pop('nfft')
    if not just_density:
        bins = kwargs['bins']  # needed for kmin,kmax
        raw_results = PS.CalculateBySlab(slicedir, nfft, BoxSize, **kwargs)
        results = to_csv(raw_results, bins, save_fn, header)
    else:
        kwargs.pop('bins', None)
        kwargs.pop('log', None)
        density = TSC.BinParticlesFromFile(slicedir, BoxSize, nfft, **kwargs)
        density.tofile(save_fn)
        results = density

    pardirname = basename(dirname(abspath(save_fn)))
    print(f'* Finished PS ({pardirname}, nfft {nfft})')

    return results, save_fn
            

def run_PS_on_PS(input_ps_fn, **kwargs):    
    outdir = get_output_path_ps(input_ps_fn)
    output_ps_fn = os.path.basename(input_ps_fn)
    
    # Make a descriptive filename
    output_ps_fn += ps_suffix(**kwargs)

    header, header_fn = common.get_header(dirname(input_ps_fn), retfn=True)

    # box kwargs overrides header
    argbox = kwargs.pop('box')
    if argbox:
        BoxSize = argbox
    else:
        BoxSize = header['BoxSize']
        
    # Load the input PS
    input_ps = np.loadtxt(input_ps_fn)
    
    nfft = kwargs.pop('nfft')
    raw_results = PS.RebinTheoryPS(input_ps, nfft, BoxSize, bins=kwargs['bins'], log=kwargs['log'], dtype=kwargs['dtype'])
    outfn = pjoin(outdir, output_ps_fn)
    results = to_csv(raw_results, bins, outfn, header)

    return results, outfn


def run_PS_on_density(input_density_fn, window=None, normalize_dens=False, **kwargs):
    # Make a descriptive filename
    # First have to figure out nfft on disk
    nelem = getsize(input_density_fn)//(np.dtype(kwargs['dtype']).itemsize)
    kwargs.pop('nfft')
    nfft = int(round(nelem**(1/3)))
    if nfft**3 != nelem:
        raise ValueError(f'density file not a cube (nelem {nelem})')
    print(f'Adopting nfft {nfft} corresponding to density cube')

    this_suffix = ps_suffix(nfft=nfft, **kwargs)
    kwargs.pop('track')
    
    inputdir = dirname(input_density_fn)

    out_parent = kwargs.pop('out_parent')
    outdir = common.get_output_dir('power', inputdir, out_parent=out_parent, **kwargs)
    ps_fn = 'power' + this_suffix

    header, header_fn = common.get_header(inputdir, retfn=True)
    
    # box kwargs overrides header
    argbox = kwargs.pop('box')
    if argbox:
        BoxSize = argbox
    else:
        BoxSize = header['BoxSize']
            
    # Make the output dir and store the header
    os.makedirs(outdir, exist_ok=True)
    if header_fn:
        shutil.copy(header_fn, pjoin(outdir, 'header'))

    # Copy the sim-level parameter files
    if not os.path.isdir(pjoin(outdir, os.pardir, 'info')):
        try:
            shutil.copytree(pjoin(inputdir, os.pardir, 'info'), pjoin(outdir, os.pardir, 'info'))
        except:
            pass
    
    print('* Starting PS on {}'.format(input_density_fn))
    save_fn = os.path.join(outdir, ps_fn)
    print(f'* and saving to {save_fn}')
    
    bins = kwargs['bins']  # needed for kmin,kmax
    #bins = Histogram.k_bin_edges(nfft, BoxSize, nbins=-1, dk=dk, log=args['log'])

    density = np.fromfile(input_density_fn, dtype=kwargs['dtype'])
    density = density.reshape(nfft,nfft,nfft)

    # TODO: want window for TSC density fields written to disk, but not straight out of zeldovich
    # TODO: always uses default binning
    raw_results = PS.FFTAndBin(density, BoxSize, bins=bins, log=kwargs['log'],
                                    window=window, normalize_dens=normalize_dens)
    results = to_csv(raw_results, bins, save_fn, header)

    pardirname = basename(dirname(abspath(save_fn)))
    print(f'* Finished PS ({pardirname}, nfft {nfft})')

    return results, save_fn

    
# inputs is a list of slice directories
def run_PS(inputs, **kwargs):
    # Loop over slices
    # Each slice may have multiple binnings/nffts
    # We process all binnings for a given slice in hopes of hitting the filesystem cache
    # TODO: implement our own cache, probably by reworking the power spectrum library into something stateful
    
    binnings = kwargs.pop('binnings')
    maxnfft = np.array([max(t['nfft'] for t in b) for b in binnings])
    order = np.argsort(maxnfft, kind='stable')[::-1]

    nthreads = kwargs.pop('nthreads',-1)
    PS.nthreads = nthreads
    verbose = kwargs.pop('verbose',False)
    PS.verbose = verbose

    # Loop in order of difficulty; i.e. largest nfft first
    # That way we'll fail-fast if we don't fit in memory
    # TODO: reverse stable argsort?
    for i in order:
        input = inputs[i]

        for binning in binnings[i]:  # binning: dict of bins, nfft, track
            if kwargs['format'] == 'density':
                _res = run_PS_on_density(input, **binning, **kwargs)
            # If the input is an output or IC directory
            elif os.path.isdir(input):
                _res = run_PS_on_dir(input, **binning, **kwargs)
            # If the input is a PS file
            elif os.path.isfile(input):
                _res = run_PS_on_PS(input, **binning, **kwargs)
            else:
                raise ValueError(input, "does not exist!")

            make_plot(*_res)


def to_csv(raw_results, bins, fn, header):
    k, Pk, nmodes = raw_results
    if type(Pk) == dict:
        power_cols = {(f'P{ell}' if ell != 0 else 'power'):Pk[ell] for ell in Pk}
    else:
        power_cols = dict(power=Pk)
    t = astropy.table.Table(data=dict(kmin=bins[:-1],
                                      kmax=bins[1:],
                                      kavg=k,
                                      **power_cols,
                                      N_modes=nmodes
                                      ),
                            meta=dict(header)
                            )
    t.write(fn, format='ascii.ecsv')
    return t


def setup_bins(args):
    '''
    -If not scale-free, just pass nbins to k_bin_edges
    -Scale bins on a per-slice basis if scale-free:
        -Set kmin to the fundamental mode of the largest scale factor
        -Set kmax to Nyquist of the smallest scale factor
        -So later slices will have smaller kmin/kmax

    TODO: might want to go past the scale-free bins (either direction)
    if our grid supports it

    TODO: this assumes the presence of a 'header' file

    Sometimes we wish to choose NFFT for each slice as well, so that
    kmax/nfft is roughly constant.

    Returns a 2D list of binning/nfft pairs for every slice.
    Shape is (nslices, ntrack); each binning is a dict with `bins`, `nfft`, and `track`.
    `ntrack` may vary from slice to slice, as not all tracks extend the full z range.
    '''

    ns = args.pop('scalefree_index')
    basea = args.pop('scalefree_base_a')
    nbins = args.pop('nbins')
    dk = args.pop('dk')
    slices = args['input']
    nfft = args.pop('nfft')
    bin_like_nfft = args.pop('bin_like_nfft')
    ntracks = args.pop('ntracks')

    ### end popping args

    if dk is None and nbins is None:
        nbins = -1

    if bool(dk) == bool(nbins):
        raise ValueError((dk,nbins))

    if ns is None:
        ntracks = 1  # tracks only make sense for scale-free
    else:
        if nbins == -1:
            nbins = nfft//4

    nslices = len(slices)

    headers = [common.get_header(p) for p in slices]
    try:
        all_scalefactor = np.array([h['ScaleFactor'] for h in headers])
    except:
        all_scalefactor = np.zeros(len(headers))

    def setup_onetrack(headers=None):
        nslices = len(headers) if headers else 1

        # "box" specified on the command line overrides the header
        if args.get('box') or not headers:
            assert args['box']
            all_boxsize = np.full(nslices, args['box'])
        else:
            all_boxsize = np.asfarray([h['BoxSize'] for h in headers])

        if ns is not None:
            all_scalefactor = np.array([h['ScaleFactor'] for h in headers])
            if basea is None:
                base_scalefactor = all_scalefactor.max()
            else:
                base_scalefactor = basea

            len_rescale = (all_scalefactor/base_scalefactor)**(2./(3+ns))

            # Define kmin and kmax at the base time (latest time)
            kmin = 2*np.pi/all_boxsize[all_scalefactor.argmax()]
            # at the earliest time, then scale to latest time
            kmax = np.pi/(all_boxsize[all_scalefactor.argmin()]/nfft)*len_rescale[all_scalefactor.argmin()]

            if args['log']:
                _bins = np.logspace(np.log10(kmin), np.log10(kmax), num=nbins+1)
            else:
                _bins = np.linspace(kmin, kmax, num=nbins+1)
            all_bins = np.tile(_bins, (nslices,1))

            # these are k bins: the scaling is the inverse of the length rescaling
            # all k will thus get bigger
            all_bins /= len_rescale[:,None]

            # Now choose the smallest even nfft larger than kmax for each time
            all_nfft = np.ceil(all_bins.max(axis=-1)*all_boxsize/np.pi).astype(int)
            all_nfft[all_nfft % 2 == 1] += 1
            #assert (all_nfft <= nfft).all(), all_nfft

            print('--scalefree_index option was specified.  Computed k ranges: ' + str([('{:.3g}'.format(b.min()), '{:.3g}'.format(b.max())) for b in all_bins]))
            print('\tComputed NFFTs: ' + str(all_nfft))
        else:
            # Set up normal bins
            all_bins = np.array([Histogram.k_bin_edges(nfft, all_boxsize[i], nbins=nbins, dk=dk, log=args['log'], bin_like_nfft=bin_like_nfft) \
                                    for i in range(nslices)])
            all_nfft = np.full(nslices, nfft)

        return all_bins, all_nfft

    isort = np.argsort(all_scalefactor, kind='stable')

    binnings = [list() for _ in range(nslices)]  # one per slice

    for t in range(ntracks):
        # Evenly spaced by number of slices
        # Would we prefer a different spacing?
        tlabel = t if ntracks > 1 else None  # for the file name

        # These are the slice indices for this track (increasing towards later times)
        i_thistrack = isort[(ntracks-1-t)*nslices//ntracks:]
        headers_thistrack = [headers[i] for i in i_thistrack]

        all_bins_thistrack, all_nfft_thistrack = setup_onetrack(headers_thistrack)

        for i,itrack in enumerate(i_thistrack):
            binnings[itrack] += [dict(bins=all_bins_thistrack[i], nfft=all_nfft_thistrack[i], track=tlabel)]

    return binnings


def make_plot(results, csvfn):
    '''
    Make a preview plot of the results.
    '''
    # hacky, but consistent with ps_suffix()
    if csvfn.endswith('.csv'):
        pltfn = csvfn[:-4] + '.png'
    else:
        pltfn = csvfn
    projected = 'projected' in csvfn
    
    zspace = 'P2' in results.colnames or 'P4' in results.colnames

    nrow=2 if zspace else 1
    fig, axes = plt.subplots(nrow,1, figsize=(6,4*nrow), sharex=True, squeeze=False)
    axes = axes.reshape(-1)

    valid = results['N_modes'] > 0
    results = results[valid]
    
    if results['power'][0]/results['power'][1] < 0.01:
        results = results[1:]  # trim first bin if small

    label = 'Auto power spectrum'
    axes[0].loglog(results['kavg'], results['power'], label=label)
    if 'P2' in results.colnames:
        axes[1].semilogx(results['kavg'], results['kavg']*results['P2'], label='Quadrupole')
    if 'P4' in results.colnames:
        axes[1].semilogx(results['kavg'], results['kavg']*results['P4'], label='Hexadecapole')
    if zspace:
        axes[1].legend()
        axes[1].set_ylabel(r'$kP(k)$ $[(\mathrm{{Mpc}}/h){exp}]$'.format(exp=('' if projected else '^2')))
    axes[0].legend()
    axes[-1].set_xlabel(r'$k$ [$h$/Mpc]')
    axes[0].set_ylabel(r'$P(k)$ $[(\mathrm{{Mpc}}/h)^{ndim:d}]$'.format(ndim=(2 if projected else 3)))

    fig.savefig(pltfn, bbox_inches='tight')
        
def vector_arg(s):
    try:
        x, y, z = list(map(int, s.strip('()[]').split(',')))
        return x, y, z
    except:
        raise argparse.ArgumentTypeError("Vector must be x,y,z")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Compute power spectra on Abacus outputs or ICs.  Can also evaluate a power spectrum on an FFT mesh.', formatter_class=Tools.ArgParseFormatter)
    parser.add_argument('input', help='The timeslice outputs (or IC directories, or power spectrum file) on which to run PS', nargs='+')
    parser.add_argument('--nfft', help='The size of the FFT (side length of the FFT cube).  Default: 1024', default=1024, type=int)
    parser.add_argument('--format', help='Format of the data to be read', choices=TSC.valid_file_formats + ['density'], type=str.lower)
    parser.add_argument('--rotate-to', help='Rotate the z-axis to the given axis [e.g. (1,2,3)].  Rotations will shrink the FFT domain by sqrt(3) to avoid cutting off particles.', type=vector_arg, metavar='(X,Y,Z)')
    #parser.add_argument('--multipoles', help='Compute these multipoles (e.g. (0,2,4))', type=list)
    parser.add_argument('--projected', help='Project the simulation along the z-axis.  Projections are done after rotations.', action='store_true')
    parser.add_argument('--zspace', help='Displace the particles according to their redshift-space positions.', action='store_true')
    parser.add_argument('--bin-like-nfft', help='Choose bins to be commensurate with this NFFT', default=0, type=int)
    parser.add_argument('--dtype', help='Data type for the binning and FFT.', choices=['float', 'double'], default='float')
    parser.add_argument('--log', help='Do the k-binning in log space.', action='store_true')
    # TODO: move just-density to different util
    parser.add_argument('--just-density', help='Write the density cube out to disk instead of doing an FFT.', action='store_true')
    parser.add_argument('--scalefree_index', help='Automatically scales k_min and k_max according to the scale free cosmology with the given spectral index.  Uses the lowest redshift of all slices as the "base time".', default=None, type=float)
    parser.add_argument('--scalefree_base_a', help='Override the fiducial scale factor for automatic k_min/k_max computation in a scale-free cosmology. Only has an effect with the --scalefree_index option.', default=None, type=float)
    parser.add_argument('--box', help='Override the box size from the header', type=float)
    parser.add_argument('--ntracks', help='The number of "tracks" of power spectra to measure.  Only has an effect with the --scalefree_index option.', type=int, default=1)
    parser.add_argument('--nthreads', help='Number of threads for binning/FFT. Default is all CPUs less ReadAbacus threads.', type=int, default=os.cpu_count())
    parser.add_argument('--nreaders', help='Number of IO threads.', type=int, default=1)
    parser.add_argument('--verbose', help='Verbose status messages', action='store_true', default=True)
    parser.add_argument('--out-parent', help='Parent directory for output')
    parser.add_argument('-o', '--outdir', help='Output directory manual override')

    bingroup = parser.add_mutually_exclusive_group()
    bingroup.add_argument('--nbins', help='Number of k bins.  Default: nfft//4.', type=int)
    bingroup.add_argument('--dk', help='k bin width', type=float)
    # TODO: wisdom option
    
    args = parser.parse_args()
    args = vars(args)

    if args.pop('projected'):
        args['nfft'] = [args['nfft'],]*2
    
    args['binnings'] = setup_bins(args)
    if args['just_density']:
        # No meaning for density cube
        del args['log']
    
    if args['dtype'] == 'float':
        args['dtype'] = np.float32
    elif args['dtype'] == 'double':
        args['dtype'] = np.float64

    with Tools.ContextTimer('* All slices'):
        run_PS(args.pop('input'), **args)
    