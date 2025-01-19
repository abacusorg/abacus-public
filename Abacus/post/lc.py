#!/usr/bin/env python
"""
Concatenate and compress the given lightcone files into ASDF file format,
one per step per file type.
"""

import gc
from pathlib import Path

import asdf
import click
import numpy as np

from Abacus.fast_cksum.cksum_io import CksumReader, CksumWriter
from Abacus.InputFile import InputFile


def check_names(files: list[Path], out: Path) -> str:
    fn0 = files[0]
    parts = fn0.name.split('_')
    lcn = parts[0]
    ftype = parts[1].split('.')[0]
    parent = fn0.parent

    if 'heal' in ftype:
        raise ValueError('Use lc_healpix.py for healpix files')

    for fn in files:
        parts = fn.name.split('_')
        if parts[0] != lcn:
            raise ValueError(f'Lightcone number mismatch: {fn0} vs {fn}')
        if parts[1].split('.')[0] != ftype:
            raise ValueError(f'File type mismatch: {fn0} vs {fn}')
        if fn.parent != parent:
            raise ValueError(f'Parent directory mismatch: {fn0} vs {fn}')

    if ftype not in out.name:
        raise ValueError(f'Output filename {out} does not match file type {ftype}')

    return ftype


@click.command
@click.argument('files', nargs=-1)
@click.option('-o', 'out', required=True, help='Output ASDF filename', metavar='FN')
@click.option('-t', '--nthread', default=4, help='Number of Blosc compression threads')
@click.option('--delete', is_flag=True, help='Delete input files after processing')
def main(files, out, verbose=True, nthread=4, delete=False):
    files = [Path(f) for f in files]
    if not files:
        return
    out = Path(out)
    ftype = check_names(files, out)
    stepdir = files[0].parent

    cksum_reader = CksumReader(stepdir / 'checksums.crc32', verbose=True)
    out.parent.mkdir(parents=True, exist_ok=True)

    if ftype == 'rv':
        Mway = 12
        pdtype = np.int32
        colname = 'rvint'
        pshape = (-1, 3)

    elif ftype == 'pid':
        Mway = 8
        pdtype = np.uint64
        colname = 'packedpid'
        pshape = (-1,)

    nptot = sum(f.stat().st_size for f in files) // pdtype().itemsize
    particles = np.empty(nptot, dtype=pdtype)

    i = 0
    for fn in files:
        new = np.frombuffer(cksum_reader(fn), dtype=pdtype)
        particles[i : i + len(new)] = new
        i += len(new)
        del new
    assert i == nptot
    gc.collect()

    particles = particles.reshape(*pshape)

    compression_kwargs = dict(
        typesize=Mway,
        compression_block_size=12 * 1024**2,
        blosc_block_size=3 * 1024**2,
        shuffle='bitshuffle',
        nthreads=nthread,
    )
    tree = {
        'data': {colname: particles},
        'header': dict(InputFile(stepdir / 'header')),
    }
    with asdf.AsdfFile(tree) as af, CksumWriter(out) as fp:
        af.write_to(
            fp, all_array_compression='blsc', compression_kwargs=compression_kwargs
        )
    print(f'Wrote to {out}')

    if delete:
        for fn in files:
            fn.unlink()


if __name__ == '__main__':
    main()
