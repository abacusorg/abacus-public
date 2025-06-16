#!/usr/bin/env python3
# Copyright 2012-2025 The Abacus Developers
# SPDX-License-Identifier: GPL-3.0-or-later

'''
The 2D code needs to read a [x,z] subdomain of the ICs.
But the zeldovich code writes out [slab][x,y,z] ([slab] are separate
files).  This script creates a [slab][z,x,y] copy of the ICs that
the 2D code can read.  It will be automatically invoked by
abacus.py.

We could split the files on z-rank, but then that makes it
harder to change NumZRanks later.  Instead, each z-rank
can just read the contiguous subset of the file it needs.

Note that the 2D ICs are also perfectly valid inputs
to the 1D code.
'''

from pathlib import Path
import multiprocessing.pool

import numpy as np
import tqdm
import click


pbytes = dict(rvzel=32)


@click.command()
@click.argument('icdir')
@click.argument('ppd', type=int)
@click.argument('format')
@click.option('--nthread', type=int)
def command(icdir, ppd, format, nthread=2):
    '''
    Transpose a set of [slab][x,y,z] files from the zeldovich
    code into [slab][z,x,y].  Outputs will be written into
    ``icdir/2D/ic2D_*``.
    '''
    transpose(icdir, ppd**3, format, nthread=nthread)


def transpose(icdir, NP, format, nthread=2):
    '''
    Transpose a set of [slab][x,y,z] files from the zeldovich
    code into [slab][z,x,y].  Outputs will be written into
    ``icdir/2D/ic2D_*``.

    Parameters
    ----------
    icdir : path-like
        Directory containing the IC files

    NP : int
        Number of particles (must be perfect cube)

    format : str
        Format of the IC particles, like 'RVZel'

    nthread : int, optional
        Number of threads for the transpose and IO
    
    '''
    icdir = Path(icdir)
    format = format.lower()

    ppd = int(round(NP**(1/3)))
    if ppd**3 != NP:
        raise ValueError(f"{NP=} expected to be perfect cube "
        "(i.e. an output of the zeldovich IC code) for 2D transpose to work")

    fns = list(sorted(icdir.glob('ic_*')))

    print(f'Transposing {len(fns)} IC files for 2D code...')
    with multiprocessing.pool.ThreadPool(nthread) as pool:
        list(tqdm.tqdm(pool.imap(lambda f: transpose_one(f, ppd, format), fns), total=len(fns)))


def transpose_one(fn, ppd, format):
    b = pbytes[format]
    ic = np.fromfile(fn, dtype=np.byte)
    ic = ic.reshape(-1,ppd,ppd,b)

    outdir = fn.parent / '2D'
    name = fn.name.replace('ic_', 'ic2D_')
    outdir.mkdir(exist_ok=True)

    ic = np.ascontiguousarray(np.moveaxis(ic, 2, 0))
    ic.tofile(outdir / name)


if __name__ == '__main__':
    command()
