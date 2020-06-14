#!/usr/bin/env python3
'''
This script is launched by disBatch, once per chunk per file type per time slice
The argument is a glob pattern to use to find the files for this chunk
'''

import argparse
import os
from os.path import dirname, basename, join as pjoin, pardir, getsize, abspath
from glob import glob
import sys
import re
import io
import gc

import numpy as np
import asdf
import asdf.compression
import blosc

from Abacus.Tools import ArgParseFormatter
from Abacus.fast_cksum.cksum_io import CksumWriter, CksumReader
from Abacus.InputFile import InputFile
from Abacus.ReadAbacus import skip_header

def process(slabpat):
    # pattern looks like slice1.000/*slab001?.L0_pack9.dat
    # Job is launched in $SIM_NAME
    slicedir = dirname(slabpat)  # slice1.000
    slabpat = basename(slabpat)  # *slab001?.L0_pack9.dat
    ftype = slabpat.split('.')[-2]  # L0_pack9
    ftype = ftype.replace('pids','pid')  # de-pluralize
    
    zdir = slicedir.replace('slice','z')  # z1.000
    aggslab = re.search(r'slab\d+', slabpat).group(0)  # slab001
    aggname = '.'.join((aggslab, ftype.replace('_','.')))  # slab001.L0.pack9

    outdir = abspath(pjoin('slices', zdir, ftype))
    os.makedirs(outdir, exist_ok=True)

    fns = sorted(glob(pjoin(slicedir, slabpat)))

    cksum_reader = CksumReader(pjoin(slicedir, 'checksums.crc32'), verbose=True)

    if 'pid' in ftype:
        Mway = 8
        pdtype = np.uint64
        colname = 'pid'
        pshape = (-1,)
        asdf_block_size = 12*1024**2
        blocksize = 3*1024**2

    else:
        Mway = 9
        pdtype = np.byte
        colname = 'pack9'
        pshape = (-1,9)
        asdf_block_size = 9*(1<<21)
        blocksize = 9*(1<<19)

    nptot = sum(getsize(_fn) for _fn in fns)//pdtype().itemsize
    particles = np.empty(nptot, dtype=pdtype)

    i = 0
    for fn in fns:
        with open(fn,'rb') as fp:
            newbytes = cksum_reader(fp)
            if 'pid' not in ftype:
                newbytesfp = io.BytesIO(newbytes)
                skip_header(newbytesfp)  # rvs have header
                newbytes = newbytes[newbytesfp.tell():]
                del newbytesfp
            new = np.frombuffer(newbytes, dtype=pdtype)
        particles[i:i+len(new)] = new
        i += len(new)
        del new, newbytes
        gc.collect()

    if 'pid' in ftype:
        assert i == nptot
    else:
        assert i <= nptot
        particles = particles[:i]

    particles = particles.reshape(*pshape)

    asdf.compression.set_compression_options(typesize=Mway, asdf_block_size=asdf_block_size, blocksize=blocksize, shuffle='bitshuffle', nthreads=4)
    tree = {'data': {colname:particles},
            'header': dict(InputFile(pjoin(slicedir,'header')))}
    asdf_fn = pjoin(outdir, aggname+'.asdf')
    with asdf.AsdfFile(tree) as af, CksumWriter(asdf_fn) as fp:
        af.write_to(fp, all_array_compression='blsc')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=ArgParseFormatter)

    parser.add_argument('slabpat', help='A time slice slab globbing pattern, like slice1.000/*slab001?.L0_pack9.dat')
    args = parser.parse_args()
    args = vars(args)

    process(**args)

    sys.exit(0)
