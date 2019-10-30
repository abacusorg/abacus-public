#!/usr/bin/env python3
'''
This is an out-of-core wrapper to Corrfunc, designed to
support pair-counting (auto-correlation or cross-correlation)
on large simulations without loading the whole particle set
in memory.

This is engineered to support Abacus slab-oriented outputs,
in which particles are already spatially sorted into
planes/slabs.

To read the output file as a Numpy record array,
once can use

>>> df = pandas.read_csv(fn)
>>> results = df.to_records(index=False)
'''

'''
TODO: could extend this to accept formats other than
particle slabs (e.g. halos) via a custom reader function.

TODO: wrap functions other than DD

TODO: multiple files per slab, as in our split L0 output
model

TODO: overlap IO and compute
'''

import matplotlib as mpl
mpl.use('Agg')

import Corrfunc
from glob import glob
import os
import os.path as path
from os.path import join as pjoin
import numpy as np
import argparse
import shutil
import pandas as pd
import multiprocessing
import contexttimer
from collections import defaultdict
import re

import utils

from Abacus import ReadAbacus
from Abacus import Tools
from Abacus.InputFile import InputFile
from Abacus.Analysis import common

def DD(primary, secondary=None, pattern='*.dat', format='pack14', dtype=np.float32, chunksize=-1,
    *args, **kwargs):
    """
    Computes DD pair counts on a slab-by-slab basis,
    returning the aggregated result.

    In detail, we may load several slabs at a time if R_max
    demands it.  R_max is maximum bin radius.

    We convert Abacus positions from unit-box to box-sized box
    (so probably Mpc/h).
    
    Parameters
    ----------
    primary: str
        The directory containing the primary particle slab files.
        
    secondary: str, optional
        The directory containing the secondary particle slab files
        for cross-correlation with the primary slabs.
        `pattern2`, `format2`, and `dtype2` can be specified if they
        are different than `directory`'s.
        Default: None.

    pattern: str, optional
        A globbing pattern for the slab files.
        Default: '*.dat'

    format: str, optional
        The format of the slab files.  This is passed
        to the ReadAcacus library.
        Default: 'pack14'

    chunksize: int, optional
        The number of slabs to read at a time.
        `chunksize` <= 0 means choose the smallest
        value allowed by R_max.  This is a tuning
        parameter; larger values will typically
        run faster but use more memory.
        TODO: guess a value for this based on available memory
        Default: -1

    verbose: str, optional
        Print status updates.  Default: True.

    dtype: np.dtype, optional
        The data type to read the particles into.
        Defaut: `np.float32`

    *args, **kwargs: optional
        Additional arguments to pass to `corrfunc.theory.DD`.

    Returns
    -------
    results: Numpy structured array
        A Corrfunc results array.  Pair counting results
        are aggregated across the whole box.

    headers: tuple of InputFile
        The header(s) describing the primary and secondary
        timeslices (the latter only if one was given).
    
    Corrfunc DD documentation
    -------------------------
    """

    verbose = kwargs.get('verbose', True)
    box = float(kwargs.get('boxsize', -1.))
    periodic = kwargs.get('periodic', True)
    try:
        rmax = float(kwargs['binfile'].max())
    except:
        raise NotImplementedError('`binfile` must be array')

    if not periodic:
        raise NotImplementedError('non-periodic not supported')
    autocorr = secondary == None
    # comptue ravg by default
    kwargs['output_ravg'] = kwargs.get('output_ravg', True)

    fns = get_fns(primary, pattern)
    cpd = len(fns)

    ds = kwargs.pop('downsample', None)
    if ds:
        assert ds == 1  # TODO: support downsampling

    # TODO: support directories without headers via file headers or arguments
    header = InputFile(pjoin(primary, 'header'))
    del primary
    assert cpd == header.CPD, "CPD does not match the number of files!"
    if box <= 0.:
        box = header.BoxSize
        kwargs['boxsize'] = box
    
    if secondary:
        pattern2 = kwargs.pop('pattern2', pattern)
        format2 = kwargs.pop('format2', format)
        dtype2 = kwargs.pop('dtype2', dtype)
        fns2 = get_fns(secondary, pattern2)
        # TODO: support incommensurate slabs, including chunksize2
        assert len(fns) == len(fns2)

        header2 = InputFile(pjoin(secondary, 'header'))
        del secondary
        # TODO: maybe we need to check more things?
        assert header.CPD == header2.CPD
        assert header.BoxSize == header2.BoxSize

    # slab width in Mpc
    swidth = len(fns)/box
    # number of neighbor slabs needed
    if chunksize <= 0:
        chunksize = int(np.ceil(rmax/swidth))
    # ensure chunksize divides cpd evenly
    # wastes some space, but simplifies the logic
    while cpd % chunksize != 0:
        chunksize += 1
    assert chunksize > rmax*swidth, "Chunk size {} not large enough for R_max {} (slab width {})?".format(chunksize, rmax, swidth)
    assert cpd//chunksize >= 3, "Chunk size too large; need at least 3 chunks for periodicity."

    if verbose:
        print('Found {} slabs'.format(cpd))
        print('Counting {} pairs in chunks of {:d} slabs'.format('autocorrelation' if autocorr else 'cross-correlation', chunksize))

    # read slabs `start` through `start+n`
    def reader(start, n, two=False):
        start %= cpd
        end = start + n  # this doesn't wrap, but it shouldn't need to
        assert end <= cpd
        _fns = fns2 if two else fns
        fns_to_read = sum((_fns[i] for i in range(start,end)), [])
        if two: p = ReadAbacus.read_many(fns_to_read, return_vel=False, format=format2, dtype=dtype2)['pos']
        else: p = ReadAbacus.read_many(fns_to_read, return_vel=False, format=format, dtype=dtype)['pos']
        # TODO: support formats where the particles are already in `box`
        p *= box
        return p

    # bootstrap the chunk reading
    next = reader(0,chunksize, two=not autocorr)
    if not autocorr:
        secondary = reader(-chunksize,chunksize, two=not autocorr)
    all_results = []

    # now loop over chunks
    for i in range(0,cpd,chunksize):
        if verbose:
            # this actually indicates the position of the "next" chunk
            print('Starting chunk {} -- {}'.format(i, i + chunksize))
            
        if autocorr:
            primary = next
            secondary = primary
        else:
            primary = reader(i,chunksize)
            prev = secondary
            secondary = next
            
            prevcross = Corrfunc.theory.DD(X1=primary[:,0], Y1=primary[:,1], Z1=primary[:,2],
                       X2=prev[:,0], Y2=prev[:,1], Z2=prev[:,2],
                        *args, autocorr=False, copy_particles=False, **kwargs)
            del prev
        del next
            
        # this might be an autocorrelation of primary, or a cross-correlation with the secondary
        auto = Corrfunc.theory.DD(X1=primary[:,0], Y1=primary[:,1], Z1=primary[:,2],
                                  X2=secondary[:,0], Y2=secondary[:,1], Z2=secondary[:,2],
                                    *args, autocorr=autocorr, copy_particles=False, **kwargs)
                     
        if autocorr:
            del secondary
        next = reader(i+chunksize,chunksize, two=not autocorr)
        
        # this will always be a cross with the next
        cross = Corrfunc.theory.DD(X1=primary[:,0], Y1=primary[:,1], Z1=primary[:,2],
                                   X2=next[:,0], Y2=next[:,1], Z2=next[:,2],
                                    *args, autocorr=False, copy_particles=False, **kwargs)
        # Corrfunc doubles autocorrelation pair counts but not cross-correlation
        if autocorr:
            cross['npairs'] *= 2
            prevcross = np.zeros_like(cross)

        # combine the auto and cross results
        assert (auto[['rmin','rmax']] == cross[['rmin','rmax']]).all()
        if not autocorr:
            assert (cross[['rmin','rmax']] == prevcross[['rmin','rmax']]).all()
        result = np.zeros_like(auto)
        result[['rmin','rmax']] = auto[['rmin','rmax']]
        result['npairs'] = auto['npairs'] + cross['npairs'] + prevcross['npairs']
        result['ravg'] = auto['ravg']*auto['npairs'] + cross['ravg']*cross['npairs'] + prevcross['ravg']*prevcross['npairs']
        result['ravg'] /= result['npairs']
        result['weightavg'] = auto['weightavg']*auto['npairs'] + cross['weightavg']*cross['npairs'] + prevcross['weightavg']*prevcross['npairs']
        result['weightavg'] /= result['npairs']

        all_results += [result]

    all_results = np.array(all_results)
    # combine the results across chunks
    assert (all_results['rmin'][:-1] == all_results['rmin'][1:]).all()
    assert (all_results['rmax'][:-1] == all_results['rmax'][1:]).all()

    results = all_results[0].copy()
    results['npairs'] = all_results['npairs'].sum(axis=0)
    results['ravg'] = (all_results['ravg']*all_results['npairs']).sum(axis=0) / results['npairs']
    results['weightavg'] = (all_results['weightavg']*all_results['npairs']).sum(axis=0) / results['npairs']

    if verbose:
        print('Finished.')

    if not autocorr:
        return results, (header, header2)
    return results, (header,)

DD.__doc__ += Corrfunc.theory.DD.__doc__


def get_fns(directory, pattern, slab_regex=r'slab(\d{4})'):
    '''
    Get all files in `directory` that match `pattern`. The
    filenames are parsed to extract the slab number, and 
    the results are returned as a dict with that key.  There
    may be more than one file per slab.
    '''

    # this is all the fns that match the pattern
    _fns = glob(pjoin(directory, pattern))
    # there may be more than one file per slab, so group them by slab number
    fns = defaultdict(list)
    for f in _fns:
        matches = re.findall(slab_regex, f)
        assert len(matches) == 1, "Couldn't determine slab number from filename"
        slabnum = int(matches[0])
        fns[slabnum] += [f]
    return fns


if __name__ == '__main__':
    parser = utils.default_argparse(doc=doc=__doc__)
    parser.add_argument('--chunksize', help='The number of slabs to read at a time. CHUNKSIZE <= 0 means choose the smallest value allowed by RMAX.  This is a tuning parameter; larger values will typically run faster but use more memory.', default=-1, type=int)
    
    args = parser.parse_args()
    args = utils.process_args(args)

    binfiles = utils.setup_bins(args)

    out_parent=args.pop('out_parent')

    for i,primary in enumerate(args.pop('primary')):
        # set up output file
        output_dir = common.get_output_dir('pair_counting', primary, out_parent=out_parent)
        this_rmax = binfiles[i].max()
        if args['secondary']:
            simname = path.basename(path.dirname(path.normpath(args['secondary'])))
            output_fn_base = pjoin(output_dir, 'cross_corr_rmax{:.2g}_{}'.format(this_rmax, simname))
        else:
            output_fn_base = pjoin(output_dir, 'auto_corr_rmax{:.2g}'.format(this_rmax))
        output_fn = output_fn_base + '.csv'
        output_fn_plot = output_fn_base + '.png'

        # Make the output dir
        try:
            os.makedirs(output_dir)
        except:
            if not path.isdir(output_dir):
                raise

        print('* Starting out-of-core pair counting using {} threads with rmax {:.2g}.  Saving result to "{}"'.format(args['nthreads'], this_rmax, output_fn))
        with contexttimer.Timer(fmt='Total time was {:.3g} s', output=True):
            results, headers = DD(primary=primary, verbose=True, binfile=binfiles[i], **args)

        # Save the header and results
        shutil.copy(pjoin(primary, 'header'), output_dir)
        if args['secondary']:
            # the header of the secondary particle set is copied as 'header2'
            shutil.copy(pjoin(args['secondary'], 'header'), pjoin(output_dir, 'header2'))
        # Copy the sim-level parameter files
        if not path.isdir(pjoin(output_dir, path.pardir, 'info')):
            try:
                shutil.copytree(pjoin(primary, path.pardir, 'info'), pjoin(output_dir, path.pardir, 'info'))
            except:
                pass

        pdresults = pd.DataFrame(results)
        pdresults.to_csv(output_fn, index=False)

        utils.make_plot(results, output_fn_plot, headers)
