#!/usr/bin/env python
'''
This is a thin wrapper for running Corrfunc on sets of
particles that can be loaded into memory all at once.

The main purpose of this script is to provide automatic
writing of results to disk.  Also, it handles a number
of different file formats (from Abacus, Gadget, etc).

For particle sets that do not fit into memory, one
can use the `out_of_corrfunc.py` script.
'''

import argparse
import shutil
from os import path
from os.path import join as pjoin
from glob import glob
import timeit

import numpy as np
import Corrfunc
import numexpr as ne
from astropy.table import Table

import utils

from Abacus import ReadAbacus
from Abacus.Tools import ContextTimer
from Abacus.Analysis import common

def run(args):
    binfiles = utils.setup_bins(args)

    out_parent = args.pop('out_parent')

    for i,primary in enumerate(args.pop('primary')):

        # The bin edges are set up as a masked array
        # We want to output the mask, but not pass it to corrfunc!
        bin_edges = binfiles[i]
        bin_edges_valid = np.ma.compressed(bin_edges)

        # set up output file

        this_rmax = bin_edges.max()
        output_dir, output_fn, output_fn_plot = utils.setup_output_file(primary, out_parent, args, this_rmax)

        print('* Starting run_corrfunc.py using {} threads with rmax {:.2g}.  Saving result to "{}"'.format(args['nthreads'], this_rmax, output_fn))

        print('* Starting to read particles...')
        with ContextTimer('Read particles'):
            crosscorr = args['secondary'] is not None
            ds = args['downsample']
            box = args.get('box')

            if args['format'] == 'gadget':
                p,_box = utils.read_gadget(primary, downsample=ds)

                if not box:
                    box = _box
                header = {'NP':len(p[0]), 'BoxSize':box}

                if crosscorr:
                    # TODO: if we have multiple primaries and one secondary this will re-read each time
                    # TODO: support format2
                    ps,_box = utils.read_gadget(args['secondary'], downsample=ds)
                    
                    assert _box == box
                    header2 = {'NP':len(ps[0]), 'BoxSize':box}
            else:
                p, header = ReadAbacus.from_dir(primary, format=args['format'], dtype=args['dtype'], return_vel=False, return_header=True)
                assert header['NP'] == len(p)
                p = p['pos'][::ds]
                
                _box = header.get('BoxSize')
                if not box:
                    box = _box
                header = dict(header)
                header.update({'NP':len(p), 'BoxSize':box})

                ne.evaluate('p*b', out=p, local_dict={'p':p,'b':p.dtype.type(box)})
                p = p.T

                if crosscorr:
                    ps, header2 = ReadAbacus.from_dir(args['secondary'], format=args['format'], dtype=args['dtype'], return_vel=False, return_header=True)
                    _box = header2.get('BoxSize')
                    if _box:
                        assert _box == box
                    ps = ps['pos'][::ds]
                    header2 = dict(header2)
                    header2.update({'NP':len(ps), 'BoxSize':box})
                    ne.evaluate('ps*b', out=ps, local_dict={'ps':ps,'b':ps.dtype.type(box)})
                    ps = ps.T

            headers = [header]
            if crosscorr:
                headers += [header2]

            if crosscorr:
                X2, Y2, Z2 = ps
                del ps
            else:
                X2, Y2, Z2 = None, None, None

        print('* Starting pair counting on {} primary particles...'.format(header['NP']))
        with ContextTimer('Pair counting'):
            results = Corrfunc.theory.DD(autocorr=not crosscorr, nthreads=args['nthreads'], binfile=bin_edges_valid,
                                X1=p[0], Y1=p[1], Z1=p[2],
                                X2=X2, Y2=Y2, Z2=Z2,
                                verbose=True, periodic=True, boxsize=box,
                                output_ravg=False,
                                max_cells_per_dim=500)
        del p, X2, Y2, Z2
        print('* Done.')

        # Save the header and results
        # We measured on the unmasked bins, but now we need to save the results including the masked bins!
        results = Table(results)

        allres = Table(np.zeros(len(bin_edges)-1, dtype=results.dtype))
        allres['rmin'] = bin_edges.data[:-1]
        allres['rmax'] = bin_edges.data[1:]

        mask = (bin_edges.mask[:-1] | bin_edges.mask[1:])
        allres[~mask] = results

        allres['mask'] = mask

        # Save the actual, downsampled number of particles
        allres.meta['N_primary'] = header['NP']
        if crosscorr:
            allres.meta['N_secondary'] = header2['NP']

        allres.write(output_fn, format='ascii.ecsv')

        utils.save_header(output_dir, primary, args)

        utils.make_plot(results, output_fn_plot, headers)

if __name__ == '__main__':
    parser = utils.default_argparse(doc=__doc__)

    args = parser.parse_args()
    args = utils.process_args(args)

    with ContextTimer('All pair counting'):
        run(args)

