#!/usr/bin/env python

'''
A user-friendly command-line interface to PowerSpectrum.

Author: Lehman Garrison
'''

import sys
import os
from os.path import join as pjoin
import pdb
import shutil
import argparse
from glob import glob

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pynbody
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
    patterns = [pjoin(folder, '*.dat'), pjoin(folder, 'ic_*'), pjoin(folder, 'position_*'), pjoin(folder, r'*.[0-9]*')]
    # Decide which pattern to use
    for pattern in patterns:
        if glob(pattern):
            break
    else:
        assert len(glob(pattern)) > 0, 'Could not find any matches to ' + str(patterns)

    # Make a descriptive filename
    this_suffix = ps_suffix(**kwargs)
    
    just_density = kwargs.pop('just_density')
    product_name = 'power' if not just_density else 'density'
    outdir = common.get_output_dir(product_name, folder)
    ps_fn = product_name + this_suffix

    # Read the header to get the boxsize
    header_pats = [pjoin(folder, 'header'), pjoin(folder, os.pardir, 'info/*.par'), pjoin(folder, os.pardir, '*.par')]
    headers = sum([glob(h) for h in header_pats], [])
    for header_fn in headers:
        try:
            header = InputFile(header_fn)
            BoxSize = header['BoxSize']
        except IOError:
            continue
        break
    else:
        try:
            del header_fn
        except:
            pass  # Maybe no headers were found 
        # Maybe it's a gadget file?
        try:
            gadget_fn = sorted(glob(pattern))[0]
            f = pynbody.load(gadget_fn)
            BoxSize = float(f.properties['boxsize'])
        except:
            print('Could not find a header in ' + str(header_pats))
            print('or as a gadget file')
            raise
            
    # Make the output dir and store the header
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    try:
        shutil.copy(header_fn, outdir + '/header')
    except:
        # If a normal header doesn't exist, try copying an example gadget file
        shutil.copy(gadget_fn, outdir + '/gadget_header')

    # Copy the sim-level parameter files
    if not os.path.isdir(outdir + '/../info'):
        try:
            shutil.copytree(folder+'/../info', outdir + '/../info')
        except:
            print('Could not copy ../info')
    
    print('Starting {} on {}'.format('PS' if not just_density else 'density', pattern))
    save_fn = os.path.join(outdir, ps_fn)
    print('and saving to {}'.format(save_fn))
    
    if not just_density:
        bins = kwargs['bins']  # needed for kmin,kmax
        raw_results = PS.CalculateBySlab(pattern, kwargs.pop('nfft'), BoxSize, **kwargs)
        results = to_csv(raw_results, bins, save_fn)

    else:
        density = TSC.BinParticlesFromFile(pattern, BoxSize, kwargs.pop('nfft'), **kwargs)
        density.tofile(save_fn)
        results = density

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

    
# folders is a list of slice folders
def run_PS(inputs, **kwargs):
    all_bins = kwargs.pop('all_bins')
    all_nfft = kwargs.pop('all_nfft')

    # Loop in order of difficulty; i.e. largest nfft first
    for i in np.argsort(all_nfft)[::-1]:
        input = inputs[i]
        # If the input is an output or IC directory
        if os.path.isdir(input):
            if kwargs.get('format').lower() == 'pack14':
                with common.extract_slabs(input):
                    _res = run_PS_on_dir(input, bins=all_bins[i], nfft=all_nfft[i], **kwargs)
            else:  # pretty kludgy (TODO: ExitStack)
                _res = run_PS_on_dir(input, **kwargs)
        # If the input is a PS file
        elif os.path.isfile(input):
            _res = run_PS_on_PS(input, **kwargs)
        else:
            raise ValueError(input, "does not exist!")

        make_plot(*_res)
        print('Finished PS.')


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
    -Scale bins on a per-sim basis if scale-free:
        -Set kmin to the fundamental mode of the largest scale factor
        -Set kmax to Nyquist of the smallest scale factor
        -So later slices will have smaller kmin/kmax

    TODO: might want to go past the scale-free bins (either direction)
    if our grid supports it

    TODO: this assumes the presence of a 'header' file

    Sometimes we wish to choose NFFT for each sim as well, so that
    kmax/nfft is roughly constant.

    Returns one bin edges array per sim, and one nfft value per sim
    '''
    ns = args.pop('scalefree_index')
    basea = args.pop('scalefree_base_a')
    nbins = args.pop('nbins')
    nslice = len(args['input'])
    nfft = args.pop('nfft')

    if nbins == -1:
        nbins = nfft//4

    headers = [InputFile(pjoin(p,'header')) for p in args['input']]
    if ns is not None:
        all_scalefactor = np.array([h.ScaleFactor for h in headers])
        if basea is None:
            base_scalefactor = all_scalefactor.max()
        else:
            base_scalefactor = basea

        len_rescale = (all_scalefactor/base_scalefactor)**(2./(3+ns))

        # Define kmin and kmax at the base time (latest time)
        kmin = 2*np.pi/headers[all_scalefactor.argmax()].BoxSize
        # at the earliest time, then scale to latest time
        kmax = np.pi/(headers[all_scalefactor.argmin()].BoxSize/nfft)*len_rescale[all_scalefactor.argmin()]

        if args['log']:
            _bins = np.logspace(np.log10(kmin), np.log10(kmax), num=nbins+1)
        else:
            _bins = np.linspace(kmin, kmax, num=nbins+1)
        all_bins = np.tile(_bins, (len(args['input']),1))

        # these are k bins: the scaling is the inverse of the length rescaling
        # all k will thus get bigger
        all_bins /= len_rescale[:,None]

        # Now choose the smallest even nfft larger than kmax for each time
        all_nfft = np.ceil(all_bins.max(axis=-1)*np.array([h.BoxSize for h in headers])/np.pi).astype(int)
        all_nfft[all_nfft % 2 == 1] += 1
        assert (all_nfft <= nfft).all()

        print('--scalefree_index option was specified.  Computed k ranges: ' + str([('{:.3g}'.format(b.min()), '{:.3g}'.format(b.max())) for b in all_bins]))
        print('\tComputed NFFTs: ' + str(all_nfft))
    else:
        # Set up normal bins
        all_bins = np.array([Histogram.k_bin_edges(nfft, headers[i].BoxSize, nbins=nbins, log=args['log']) \
                                for i in range(nslice)])
        all_nfft = np.full(nslice, nfft)

    return all_bins, all_nfft


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

    label = 'Auto power spectrum'
    ax.loglog(results['kavg'], results['power'], label=label)
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
    parser.add_argument('--format', help='Format of the data to be read.  Default: Pack14', default='Pack14', choices=['RVdouble', 'Pack14', 'RVZel', 'state', 'gadget'])
    parser.add_argument('--rotate-to', help='Rotate the z-axis to the given axis [e.g. (1,2,3)].  Rotations will shrink the FFT domain by sqrt(3) to avoid cutting off particles.', default=None, type=vector_arg, metavar='(X,Y,Z)')
    parser.add_argument('--projected', help='Project the simulation along the z-axis.  Projections are done after rotations.', action='store_true')
    parser.add_argument('--zspace', help='Displace the particles according to their redshift-space positions.', action='store_true')
    parser.add_argument('--nbins', help='Number of k bins.  Default: nfft//4.', default=-1, type=int)
    parser.add_argument('--dtype', help='Data type for the binning and FFT.', choices=['float', 'double'], default='float')
    parser.add_argument('--log', help='Do the k-binning in log space.', action='store_true')
    # TODO: move just-density to different util
    parser.add_argument('--just-density', help='Write the density cube out to disk instead of doing an FFT.', action='store_true')
    parser.add_argument('--scalefree_index', help='Automatically scales k_min and k_max according to the scale free cosmology with the given spectral index.  Uses the lowest redshift of all slices as the "base time".', default=None, type=float)
    parser.add_argument('--scalefree_base_a', help='Override the fiducial scale factor for automatic k_min/k_max computation in a scale-free cosmology. Only has an effect with the --scalefree_index option.', default=None, type=float)
    
    args = parser.parse_args()
    args = vars(args)

    if args.pop('projected'):
        args['nfft'] = [args['nfft'],]*2
    
    if args['just_density']:
        # No meaning for density cube
        del args['nbins'], args['log']
    else:
        args['all_bins'], args['all_nfft'] = setup_bins(args)
    
    if args['dtype'] == 'float':
        args['dtype'] = np.float32
    elif args['dtype'] == 'double':
        args['dtype'] = np.float64

    run_PS(args.pop('input'), **args)
    