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

import numpy as np
import pandas as pd
import Corrfunc
import numexpr as ne

import utils

from Abacus import ReadAbacus
from Abacus.Analysis import common

if __name__ == '__main__':
    parser = utils.default_argparse(doc=__doc__)
    args = parser.parse_args()
    args = utils.process_args(args)

    binfiles = utils.setup_bins(args)

    out_parent = args.pop('out_parent')

    for i,primary in enumerate(args.pop('primary')):
        # set up output file
        this_rmax = binfiles[i].max()
        output_dir, output_fn, output_fn_plot = utils.setup_output_file(primary, out_parent, args, this_rmax)

        print('* Starting run_corrfunc.py using {} threads with rmax {:.2g}.  Saving result to "{}"'.format(args['nthreads'], this_rmax, output_fn))

        print('* Starting to read particles...')
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
        else:
            X2, Y2, Z2 = None, None, None

        print('* Done.  Starting pair counting on {} primary particles...'.format(header['NP']))
        results = Corrfunc.theory.DD(autocorr=not crosscorr, nthreads=args['nthreads'], binfile=binfiles[i],
                            X1=p[0], Y1=p[1], Z1=p[2],
                            X2=X2, Y2=Y2, Z2=Z2,
                            verbose=True, periodic=True, boxsize=box,
                            output_ravg=False)
        del p
        print('* Done.')

        # Save the header and results

        pdresults = pd.DataFrame(results)
        with open(output_fn, 'w') as fp:
            # Save the actual, downsampled number of particles
            fp.write('# N_primary = {}\n'.format(header['NP']))
            if crosscorr:
                fp.write('# N_secondary = {}\n'.format(header2['NP']))
            pdresults.to_csv(fp, index=False)

        utils.save_header(output_dir, primary, args)

        utils.make_plot(results, output_fn_plot, headers)
