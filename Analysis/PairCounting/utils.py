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
    -Scale bins on a per-slice basis if scale-free
    -Choose log/linear binning

    Returns one bin edges array per slice
    '''
    ns = args.pop('scalefree_index')
    #basea = args.pop('scalefree_base_a', None)

    rmin, rmax, nbin_per_decade = args.pop('rmin'), args.pop('rmax'), args.pop('nbins_decade')

    rmax_mask = args.pop('rmax_mask')
    if rmax_mask is None:
        rmax_mask = rmax

    all_headers = [common.get_header(p) for p in args['primary']]

    # Parse rmin input
    # User can enter 'eps/10' for 10th softening, for example
    try:
        rmin = float(rmin)
    except:
        all_eps = np.array([h['SofteningLength'] for h in all_headers])
        eps = all_eps[0]

        if 'eps' in rmin and (eps != all_eps).any():
            raise ValueError("TODO: Cannot use 'eps' in rmin if all eps don't agree.", all_eps)

        rmin = eval(rmin)

        del eps

    # rmin is defined at the last time
    # scale it back to the first time; that's where our fiducial bins are defined!
    rmin_last = rmin
    if ns is not None:
        all_scalefactor = np.array([h['ScaleFactor'] for h in all_headers])
        
        firsta = all_scalefactor.min()
        lasta = all_scalefactor.max()

        len_rescale = (all_scalefactor/firsta)**(2./(3+ns))  # greater than 1

        rmin /= len_rescale[all_scalefactor.argmax()]

    # Parse rmax input
    # User can enter 'box/500' for .002*box, for example
    try:
        rmax = float(rmax)
    except:
        all_box = np.array([h['BoxSize'] for h in all_headers])
        box = all_box[0]

        rmax = eval(box)

        if (box != all_box).any():
            raise ValueError("TODO: Cannot use 'box' in rmax if all box don't agree.", all_box)

        del box

    # set up base bins (at the earliest time)
    # TODO: should use fixed bin widths, this means bin widths may change with rmax
    nbins = int(np.log10(rmax/rmin)*nbin_per_decade)

    if args.pop('linear'):
        bins = np.linspace(rmin, rmax, nbins+1)
        print('Bin spacing:',bins[1] - bins[0])
    else:
        bins = np.geomspace(rmin, rmax, nbins+1)
        print('Bin log10 spacing:',np.log10(bins[1]/bins[0]))
    all_bins = np.tile(bins, (len(args['primary']),1))

    if ns is not None:
        all_bins *= len_rescale[:,None]

        # Previously, we let rmax grow with the non-linear scale.  This facilitated easy ratios of the results across redshift.
        # But per Michael Joyce, the resolution scale actually *decreases* as 1/a, so it's not actually interesting to go to large rmax at late times
        # So the rmax that the user inputs is now used at the *earliest* time, with later rmaxes decreasing according to an estimate of the resolution scale
        # rmin is set at the latest time, so earlier times probe smaller scales

        all_n1d = np.array([h['NP']**(1/3) for h in all_headers])
        assert all(all_n1d[0] == all_n1d), ("NP must agree across slices in order to determine resolution scale for scale-free binning", all_n1d)

        all_box = np.array([h['BoxSize'] for h in all_headers])
        assert all(all_box[0] == all_box), ("BoxSize must agree across slices in order to determine resolution scale for scale-free binning", all_box)

        resolution_scale = args['resolution_scale']*all_box[0]/all_n1d[0]  # 70*(particle spacing), a kludgy guess!
        res_rescale = (all_scalefactor/firsta)**-1
        resolution_scale *= res_rescale[:,None]

        # Now mask bins above the resolution cut
        all_bins = np.ma.masked_greater(all_bins, np.minimum(rmax_mask, resolution_scale))
        # and bins below rmin *in un-rescaled coords*
        all_bins.mask |= all_bins < rmin_last

        print('--scalefree_index option was specified.  Computed rmaxes: ' + str(all_bins.max(axis=-1)))

        #print(list(zip(all_bins.min(axis=-1), all_bins.max(axis=-1))))
    
    return all_bins


def default_argparse(doc=__doc__):
    parser = argparse.ArgumentParser(description=doc, formatter_class=Tools.ArgParseFormatter)
    parser.add_argument('primary', help='The time slice directory containing the particles.', nargs='+')
    # We want the user to think about rmin and rmax, since any defaults may not be optimal for cross-sim comparison, which is a very common use-case
    parser.add_argument('rmin', help='Minimum pair counting distance (Mpc/h).  Accepts simple math like "eps/10".')
    parser.add_argument('rmax', help='Maximum pair counting distance (Mpc/h). Accepts simple math like "box/500".')
    parser.add_argument('--rmax-mask', help='For scale-free binning, actually only do pair counting up to this RMAX_MASK.  Bins will be set up as if for RMAX though.', type=float, default=None)
    parser.add_argument('--resolution-scale', help='For scale-free binning, mask bins above this effective resolution scale (units of interparticle spacing)', type=float, default=35)
    parser.add_argument('--secondary', help='The time slice directory containing the secondary particles for cross-correlations.')
    
    parser.add_argument('--format', help='Format of the particle data.', default='Pack14', choices=['RVdouble', 'Pack14', 'Pack9', 'RVint', 'RVZel', 'RVTag', 'state', 'gadget'])
    parser.add_argument('--scalefree-index', help='Automatically scales rmin and rmax according to the scale free cosmology with the given spectral index.  Uses the highest redshift for rmax, and the lowest redshift for rmin.', default=None, type=float)
    # TODO: change to amin and amax
    #parser.add_argument('--scalefree_base_a', help='Override the fiducial scale factor for automatic rmin/rmax computation in a scale-free cosmology. Only has an effect with the --scalefree_index option.', default=None, type=float)
    parser.add_argument('--nbins-decade', help='Number of radial bins *per decade*.', default=40, type=int)
    parser.add_argument('--dtype', help='Data type for internal calculations.', choices=['float', 'double'], default='float')
    parser.add_argument('--linear', help='Do the histogramming in linear space.  Default is log10 space.', default=False, action='store_true')
    parser.add_argument('--out-parent', help='Directory in which to place the data products.', default=None)
    parser.add_argument('--nthreads', help='Number of OpenMP threads for Corrfunc.  NTHREADS <= 0 means use all cores.', default=-1, type=int)
    parser.add_argument('--box', help='Override the box size from the header', type=float)

    parser.add_argument('--zspace', help='Displace the particles according to their redshift-space positions.', action='store_true')
    parser.add_argument('--downsample', help='The fraction by which to downselect particles before pair counting. Useful for accelerating the computation.', type=float)

    return parser

def process_args(args):
    args = vars(args)
    args['dtype'] = {'float':np.float32, 'double':np.float64}[args['dtype']]

    # set up threads
    nthreads = args.pop('nthreads', -1)
    if nthreads <= 0:
        nthreads = multiprocessing.cpu_count()
    args['nthreads'] = nthreads

    # homogenize a few args to lowercase
    for k in ['rmin', 'rmax', 'format']:
        if type(args[k]) == str:
            args[k] = args[k].lower()

    return args


def setup_output_file(primary, out_parent, args, this_rmax):
    output_dir = common.get_output_dir('pair_counting', primary, out_parent=out_parent)
    zspace = int(args['zspace'])
    try: 
        ds = float(args['downsample'])
    except TypeError: 
        ds = False 

    if args['secondary']:
        simname = path.basename(path.dirname(path.normpath(args['secondary'])))
        output_fn_base = pjoin(output_dir, 'cross_corr_rmax{:.3g}_{}_zspace_{}'.format(this_rmax, simname, zspace))
    else:
        output_fn_base = pjoin(output_dir, 'auto_corr_rmax{:.3g}_zspace_{}'.format(this_rmax, zspace))

    if ds:
        output_fn_base      += '_ds_{:.3g}'.format(ds)
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

    plt.close(fig)  # why doesn't this happen when fig goes out of scope?

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
