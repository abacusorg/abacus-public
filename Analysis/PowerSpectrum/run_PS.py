#!/usr/bin/env python

'''
A user-friendly command-line interface to PowerSpectrum.

Usage: $ ./run_PS.py --help

Author: Lehman Garrison
'''

import sys
import os
from os.path import join as pjoin, abspath, basename, dirname
import pdb
import shutil
import argparse
from glob import glob
from contextlib import ExitStack

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from Abacus.Analysis import common
from Abacus.InputFile import InputFile
from Abacus import Tools

from PowerSpectrum import TSC, Histogram
from PowerSpectrum import PowerSpectrum as PS

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

    
def run_PS_on_dir(folder, **kwargs):
    patterns = [pjoin(folder, '*.dat'), pjoin(folder, 'ic_*'), pjoin(folder, 'position_*'), pjoin(folder, common.gadget_pattern)]
    # Decide which pattern to use
    for pattern in patterns:
        if glob(pattern):
            break
    else:
        assert len(glob(pattern)) > 0, 'Could not find any matches to ' + str(patterns)

    # Make a descriptive filename
    this_suffix = ps_suffix(**kwargs)
    kwargs.pop('track')
    
    just_density = kwargs.pop('just_density')
    product_name = 'power' if not just_density else 'density'
    outdir = common.get_output_dir(product_name, folder)
    ps_fn = product_name + this_suffix

    header, header_fn = common.get_header(folder, retfn=True)
    
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
    if not os.path.isdir(outdir + '/../info'):
        try:
            shutil.copytree(folder+'/../info', outdir + '/../info')
        except:
            pass
    
    print('* Starting {} on {}'.format('PS' if not just_density else 'density', pattern))
    save_fn = os.path.join(outdir, ps_fn)
    print('* and saving to {}'.format(save_fn))
    
    nfft = kwargs.pop('nfft')
    if not just_density:
        bins = kwargs['bins']  # needed for kmin,kmax
        raw_results = PS.CalculateBySlab(pattern, nfft, BoxSize, **kwargs)
        results = to_csv(raw_results, bins, save_fn)

    else:
        density = TSC.BinParticlesFromFile(pattern, BoxSize, nfft, **kwargs)
        density.tofile(save_fn)
        results = density

    pardirname = basename(dirname(abspath(save_fn)))
    print(f'* Finished PS ({pardirname}, nfft {nfft})')

    return results, header, save_fn
            

def run_PS_on_PS(input_ps_fn, **kwargs):    
    outdir = get_output_path_ps(input_ps_fn)
    output_ps_fn = os.path.basename(input_ps_fn)
    
    # Make a descriptive filename
    output_ps_fn += ps_suffix(**kwargs)

    # Read the header
    try:
        header_fn = glob(os.path.dirname(input_ps_fn)+'/*.par')[0]
        header = InputFile(header_fn)
    except IOError:
        print('Could not find "*.par"')
        
    # Load the input PS
    input_ps = np.loadtxt(input_ps_fn)
    
    raw_results = PS.RebinTheoryPS(input_ps, kwargs['nfft'], header.BoxSize, nbins=kwargs['nbins'], log=kwargs['log'], dtype=kwargs['dtype'])
    outfn = pjoin(outdir, output_ps_fn)
    results = to_csv(raw_results, bins, outfn)

    return results, header, outfn

    
# inputs is a list of slice directories
def run_PS(inputs, **kwargs):
    # Loop over slices
    # Each slice may have multiple binnings/nffts
    # We process all binnings for a given slice in hopes of hitting the filesystem cache
    # TODO: implement our own cache, probably by reworking the power spectrum library into something stateful
    
    binnings = kwargs.pop('binnings')
    maxnfft = np.array([max(t['nfft'] for t in b) for b in binnings])
    order = np.argsort(maxnfft, kind='stable')[::-1]

    # Loop in order of difficulty; i.e. largest nfft first
    # That way we'll fail-fast if we don't fit in memory
    # TODO: reverse stable argsort?
    for i in order:
        input = inputs[i]

        for binning in binnings[i]:  # binning: dict of bins, nfft, track
            # If the input is an output or IC directory
            if os.path.isdir(input):
                with ExitStack() as stack:
                    if kwargs.get('format').lower() == 'pack14':
                        stack.enter_context(common.extract_slabs(input))
                    _res = run_PS_on_dir(input, **binning, **kwargs)
            # If the input is a PS file
            elif os.path.isfile(input):
                _res = run_PS_on_PS(input, **kwargs)
            else:
                raise ValueError(input, "does not exist!")

            make_plot(*_res)


def to_csv(raw_results, bins, fn):
    k, Pk, nmodes = raw_results
    results = np.empty(len(k), dtype=[('kmin',float), ('kmax',float), ('kavg',float), ('power',float), ('N_modes',int)])
    results['kmin'] = bins[:-1]
    results['kmax'] = bins[1:]
    results['kavg'] = k
    results['power'] = Pk
    results['N_modes'] = nmodes
    pdresults = pd.DataFrame(results)
    pdresults.to_csv(fn, index=False)

    return results


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
    slices = args['input']
    nfft = args.pop('nfft')
    bin_like_nfft = args.pop('bin_like_nfft')
    ntracks = args.pop('ntracks')

    ### end popping args

    if ns is None:
        ntracks = 1  # tracks only make sense for scale-free
    else:
        if nbins == -1:
            nbins = nfft//4

    nslices = len(slices)

    headers = [common.get_header(p) for p in slices]
    all_scalefactor = np.array([h['ScaleFactor'] for h in headers])

    def setup_onetrack(headers):
        nslices = len(headers)

        # "box" specified on the command line overrides the header
        if args.get('box'):
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
            assert (all_nfft <= nfft).all(), all_nfft

            print('--scalefree_index option was specified.  Computed k ranges: ' + str([('{:.3g}'.format(b.min()), '{:.3g}'.format(b.max())) for b in all_bins]))
            print('\tComputed NFFTs: ' + str(all_nfft))
        else:
            # Set up normal bins
            all_bins = np.array([Histogram.k_bin_edges(nfft, all_boxsize[i], nbins=nbins, log=args['log'], bin_like_nfft=bin_like_nfft) \
                                    for i in range(nslices)])
            all_nfft = np.full(nslices, nfft)

        return all_bins, all_nfft

    isort = np.argsort(all_scalefactor, kind='stable')

    binnings = [list() for _ in range(nslices)]  # one per slice

    for t in range(ntracks):
        # Evenly spaced by number of slices
        # Would we prefer a different spacing?
        tlabel = t if ntracks > 1 else None

        # These are the slice indices for this track (increasing towards later times)
        i_thistrack = isort[t*nslices//ntracks:]
        headers_thistrack = [headers[i] for i in i_thistrack]

        all_bins_thistrack, all_nfft_thistrack = setup_onetrack(headers_thistrack)

        for i,itrack in enumerate(i_thistrack):
            binnings[itrack] += [dict(bins=all_bins_thistrack[i], nfft=all_nfft_thistrack[i], track=tlabel)]

    return binnings


def make_plot(results, header, csvfn):
    '''
    Make a preview plot of the results.
    '''
    # hacky, but consistent with ps_suffix()
    if csvfn.endswith('.csv'):
        pltfn = csvfn[:-4] + '.png'
    else:
        pltfn = csvfn
    projected = 'projected' in csvfn

    fig, ax = plt.subplots()

    valid = results['N_modes'] > 0

    label = 'Auto power spectrum'
    ax.loglog(results[valid]['kavg'], results[valid]['power'], label=label)
    ax.legend()
    ax.set_xlabel(r'$k$ [$h$/Mpc]')
    ax.set_ylabel(r'$P(k)$ $[(\mathrm{{Mpc}}/h)^{ndim:d}]$'.format(ndim=(2 if projected else 3)))

    fig.savefig(pltfn)
        
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
    parser.add_argument('--format', help='Format of the data to be read.  Default: Pack14', default='Pack14', choices=['RVdouble', 'Pack14', 'RVZel', 'state', 'gadget', 'RVTag'])
    parser.add_argument('--rotate-to', help='Rotate the z-axis to the given axis [e.g. (1,2,3)].  Rotations will shrink the FFT domain by sqrt(3) to avoid cutting off particles.', default=None, type=vector_arg, metavar='(X,Y,Z)')
    parser.add_argument('--projected', help='Project the simulation along the z-axis.  Projections are done after rotations.', action='store_true')
    parser.add_argument('--zspace', help='Displace the particles according to their redshift-space positions.', action='store_true')
    parser.add_argument('--nbins', help='Number of k bins.  Default: nfft//4.', default=-1, type=int)
    parser.add_argument('--bin-like-nfft', help='Choose bins to be commensurate with this NFFT', default=0, type=int)
    parser.add_argument('--dtype', help='Data type for the binning and FFT.', choices=['float', 'double'], default='float')
    parser.add_argument('--log', help='Do the k-binning in log space.', action='store_true')
    # TODO: move just-density to different util
    parser.add_argument('--just-density', help='Write the density cube out to disk instead of doing an FFT.', action='store_true')
    parser.add_argument('--scalefree_index', help='Automatically scales k_min and k_max according to the scale free cosmology with the given spectral index.  Uses the lowest redshift of all slices as the "base time".', default=None, type=float)
    parser.add_argument('--scalefree_base_a', help='Override the fiducial scale factor for automatic k_min/k_max computation in a scale-free cosmology. Only has an effect with the --scalefree_index option.', default=None, type=float)
    parser.add_argument('--box', help='Override the box size from the header', type=float)
    parser.add_argument('--ntracks', help='The number of "tracks" of power spectra to measure.  Only has an effect with the --scalefree_index option.', type=int, default=1)
    
    args = parser.parse_args()
    args = vars(args)

    if args.pop('projected'):
        args['nfft'] = [args['nfft'],]*2
    
    if args['just_density']:
        # No meaning for density cube
        del args['nbins'], args['log']
    else:
        args['binnings'] = setup_bins(args)
    
    if args['dtype'] == 'float':
        args['dtype'] = np.float32
    elif args['dtype'] == 'double':
        args['dtype'] = np.float64

    run_PS(args.pop('input'), **args)
    