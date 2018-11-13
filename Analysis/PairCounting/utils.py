'''
Utilities for the pair counting scripts.
'''

import argparse
import multiprocessing
import os
from os import path
from os.path import join as pjoin
import shutil
from glob import glob

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from Abacus import Tools
from Abacus.InputFile import InputFile
from Abacus.Analysis import common

def setup_bins(args):
    '''
    -Set rmin default to eps/3
    -Scale bins on a per-sim basis if scale-free
    -Choose log/linear binning

    Returns one bin edges array per sim
    '''
    ns = args.pop('scalefree_index')
    basea = args.pop('scalefree_base_a')

    rmin, rmax, nbins = args.pop('rmin'), args.pop('rmax'), args.pop('nbins')

    all_headers = [common.get_header(p) for p in args['primary']]

    # auto defaults to softening/3
    if rmin.lower() == 'auto':
        all_eps = np.array([h['SofteningLength'] for h in all_headers])
        rmin = all_eps.min()/3.

        if ns is not None:
            if (all_eps[0] != all_eps).any():
                raise ValueError("Cannot use auto rmin and specify auto scale-free rescaling if all don't agree.", all_eps)
    else:
        rmin = float(rmin)

    # auto defaults to .002*BoxSize
    if rmax.lower() == 'auto':
        all_box = np.array([h['BoxSize'] for h in all_headers])
        rmax = all_box.max()*0.002

        if ns is not None:
            if (all_box[0] != all_box).any():
                raise ValueError("Cannot use auto rmax and specify auto scale-free rescaling if all don't agree.", all_box)
    else:
        rmax = float(rmax)

    # set up base bins
    if args.pop('linear'):
        bins = np.linspace(rmin, rmax, nbins)
    else:
        bins = np.logspace(np.log10(rmin), np.log10(rmax), nbins)
    all_bins = np.tile(bins, (len(args['primary']),1))

    if ns is not None:
        all_scalefactor = np.array([h['ScaleFactor'] for h in all_headers])
        if basea is None:
            base_scalefactor = all_scalefactor.max()
        else:
            base_scalefactor = basea
        len_rescale = (all_scalefactor/base_scalefactor)**(2./(3+ns))
        all_bins *= len_rescale[:,None]

        print('--scalefree_index option was specified.  Computed rmaxes: ' + str(all_bins.max(axis=-1)))
    
    return all_bins


def default_argparse():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=Tools.ArgParseFormatter)
    parser.add_argument('primary', help='The time slice directory containing the particles.', nargs='+')
    # We want the user to think about rmin and rmax, since any defaults may not be optimal for cross-sim comparison, which is a very common use-case
    parser.add_argument('rmin', help='Minimum pair counting distance (Mpc/h).  Value of "auto" will use 1/3 of softening.')
    parser.add_argument('rmax', help='Maximum pair counting distance (Mpc/h). Value of "auto" will use 1/500 of the box size.')
    parser.add_argument('--secondary', help='The time slice directory containing the secondary particles for cross-correlations.')
    
    parser.add_argument('--format', help='Format of the particle data.', default='Pack14', choices=['RVdouble', 'Pack14', 'RVZel', 'RVTag', 'state', 'gadget'])
    parser.add_argument('--scalefree_index', help='Automatically scales rmin and rmax according to the scale free cosmology with the given spectral index.  Uses the lowest redshift of all primaries as the "base time".', default=None, type=float)
    parser.add_argument('--scalefree_base_a', help='Override the fiducial scale factor for automatic rmin/rmax computation in a scale-free cosmology. Only has an effect with the --scalefree_index option.', default=None, type=float)
    parser.add_argument('--nbins', help='Number of radial bins.', default=100, type=int)
    parser.add_argument('--dtype', help='Data type for internal calculations.', choices=['float', 'double'], default='float')
    parser.add_argument('--linear', help='Do the histogramming in linear space.  Default is log10 space.', default=False, action='store_true')
    parser.add_argument('--out_parent', help='Directory in which to place the data products.', default=None)
    parser.add_argument('--nthreads', help='Number of OpenMP threads for Corrfunc.  NTHREADS <= 0 means use all cores.', default=-1, type=int)
    parser.add_argument('--box', help='Override the box size from the header', type=float)

    #parser.add_argument('--zspace', help='Displace the particles according to their redshift-space positions.', action='store_true')
    # TODO: support non-integer downsampling
    parser.add_argument('--downsample', help='The factor by which to randomly downselect particles before pair counting. Useful for accelerating the computation.', type=int)

    return parser

def process_args(args):
    args = vars(args)
    args['dtype'] = {'float':np.float32, 'double':np.float64}[args['dtype']]

    # set up threads
    nthreads = args.pop('nthreads', -1)
    if nthreads <= 0:
        nthreads = multiprocessing.cpu_count()
    args['nthreads'] = nthreads

    return args


def setup_output_file(primary, out_parent, args, this_rmax):
    output_dir = common.get_output_dir('pair_counting', primary, out_parent=out_parent)
    if args['secondary']:
        simname = path.basename(path.dirname(path.normpath(args['secondary'])))
        output_fn_base = pjoin(output_dir, 'cross_corr_rmax{:.2g}_{}'.format(this_rmax, simname))
    else:
        output_fn_base = pjoin(output_dir, 'auto_corr_rmax{:.2g}'.format(this_rmax))
    output_fn = output_fn_base + '.csv'
    output_fn_plot = output_fn_base + '.png'

    # Make the output dir
    os.makedirs(output_dir, exist_ok=True)

    return output_dir, output_fn, output_fn_plot


def save_header(output_dir, primary, args):
    '''
    Try to save the Abacus header in the output directory,
    but recognize that we also allow processing of non-Abacus
    particles like Gadget that may not have headers.
    '''
    for fn in ['header', 'state']:
        try:
            shutil.copy(pjoin(primary, fn), pjoin(output_dir, 'header'))
            break
        except FileNotFoundError as e:
            pass
    else:
        if args['format'] != 'gadget':
            raise e
            

    if args['secondary']:
        try:
            # the header of the secondary particle set is copied as 'header2'
            shutil.copy(pjoin(args['secondary'], 'header'), pjoin(output_dir, 'header2'))
        except shutil.FileNotFoundError:
            if args['format'] != 'gadget':
                raise

    # Copy the sim-level parameter files
    if not path.isdir(pjoin(output_dir, path.pardir, 'info')):
        try:
            shutil.copytree(pjoin(primary, path.pardir, 'info'), pjoin(output_dir, path.pardir, 'info'))
        except:
            pass


def make_plot(results, fn, headers):
    '''
    Make a preview plot of the results.
    '''
    fig, ax = plt.subplots()
    bin_centers = (results['rmin'] + results['rmax'])/2.

    header = headers[0]
    NP = header['NP']

    bin_vol = 4*np.pi/3.*(results['rmax']**3 - results['rmin']**3)
    if len(headers) > 1:
        header2 = headers[1]
        RR = NP*header2['NP']/header2['BoxSize']**3.*bin_vol
        label = 'cross-correlation'
    else:
        RR = NP*(NP-1)/header['BoxSize']**3.*bin_vol
        label = 'auto-correlation'

    ax.semilogx(bin_centers, results['npairs']/RR - 1, label=label)
    ax.legend()
    ax.set_xlabel(r'$r$ [Mpc/$h$]')
    ax.set_ylabel(r'$\xi(r)$')
    fig.tight_layout()

    fig.savefig(fn)

def read_gadget(dir, downsample=1):
    '''
    Read Gadget positions into an array.
    Really pynbody or nbodykit should be able to do this for us,
    but the former crashes and the latter runs out of memory.
    So we need to be more clever.

    It's tempting to migrate this to ReadAbacus.  Do we really
    want Gadget in there though?
    '''

    import nbodykit
    nbodykit.GlobalCache.resize(0)

    from nbodykit.source import Gadget1Catalog

    globpat = pjoin(dir, common.gadget_pattern)
    #globpat = pjoin(dir, '*.0')
    cdefs = [('Position',('f4',3),'all')]

    wholesnap = Gadget1Catalog(globpat, columndefs=cdefs)

    box = wholesnap.attrs['BoxSize']
    print("* Gadget box size: {}".format(box))
    print("* Found {} Gadget particles".format(len(wholesnap)))

    N = len(wholesnap['Position'])
    dt = wholesnap['Position'].dtype
    # Store as x,y,z columns so corrfunc doesn't need a transpose copy
    allp = [np.empty(N, dtype=dt),
            np.empty(N, dtype=dt),
            np.empty(N, dtype=dt)]

    del wholesnap

    nsaved = 0
    for fn in sorted(glob(globpat)):
        gcat = Gadget1Catalog(fn, columndefs=cdefs)
        thisp = gcat['Position'].compute()
        thisn = len(thisp[::downsample])
        allp[0][nsaved:nsaved+thisn] = thisp[::downsample,0]
        allp[1][nsaved:nsaved+thisn] = thisp[::downsample,1]
        allp[2][nsaved:nsaved+thisn] = thisp[::downsample,2]
        nsaved += thisn
        del gcat, thisp
        print('*   Done', fn)

    # Now shrink the arrays to the downsampled size
    allp = [p[:nsaved] for p in allp]

    return allp, box
