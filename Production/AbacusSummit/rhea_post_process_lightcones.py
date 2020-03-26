#!/usr/bin/env python3
'''
Concat lightcone files for a given prefix, packing into ASDF+blosc
This script is launched by disBatch, once per light cone per file type per step.
'''

import argparse
import os
from os.path import dirname, basename, join as pjoin, pardir, getsize, abspath
from glob import glob
import sys

import numpy as np
import asdf
import asdf.compression

from Abacus.Tools import ArgParseFormatter
from Abacus.fast_cksum.cksum_io import CksumWriter, CksumReader
from Abacus.InputFile import InputFile

def process(lcprefix):
    stepdir = dirname(lcprefix)
    lcprefix = basename(lcprefix)
    lcnum, ftype = lcprefix.split('_')

    outdir = abspath(pjoin(pardir, 'lightcones.concat', ftype))
    os.makedirs(outdir, exist_ok=True)

    os.chdir(stepdir)

    fns = sorted(glob(lcprefix + '.*'))

    cksum_reader = CksumReader('checksums.crc32', verbose=True)

    if ftype == 'heal':
        Mway = 4
        pdtype = np.uint32
        colname = 'heal'
        pshape = (-1,)

    elif ftype == 'rv':
        Mway = 12
        pdtype = np.int32
        colname = 'rvint'
        pshape = (-1,3)

    elif ftype == 'pid':
        Mway = 8
        pdtype = np.uint64
        colname = 'packedpid'
        pshape = (-1,)

    nptot = sum(getsize(_fn) for _fn in fns)//pdtype().itemsize
    particles = np.empty(nptot, dtype=pdtype)

    i = 0
    for fn in fns:
        new = np.frombuffer(cksum_reader(fn), dtype=pdtype)
        particles[i:i+len(new)] = new
        i += len(new)
        del new
    assert i == nptot

    particles = particles.reshape(*pshape)

    asdf.compression.set_compression_options(typesize=Mway, asdf_block_size=12*1024**2, blocksize=3*1024**2, shuffle='bitshuffle', nthreads=2)
    tree = {'data': {colname:particles},
            'header': dict(InputFile('header'))}
    asdf_fn = pjoin(outdir, f'{lcprefix}_{stepdir}.asdf')
    with asdf.AsdfFile(tree) as af, CksumWriter(asdf_fn) as fp:
        af.write_to(fp, all_array_compression='blsc')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=ArgParseFormatter)

    parser.add_argument('lcprefix', help='LC prefix like "Step0800/LightCone0_heal"')
    args = parser.parse_args()
    args = vars(args)

    process(**args)

    sys.exit(0)
