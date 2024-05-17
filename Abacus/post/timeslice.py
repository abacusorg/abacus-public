#!/usr/bin/env python
'''
Concatenate and compress the given slabs into a superslab in ASDF file format.
'''

import gc
import io
from pathlib import Path

import asdf
import click
import numpy as np

from Abacus.fast_cksum.cksum_io import CksumReader, CksumWriter
from Abacus.InputFile import InputFile
from Abacus.ReadAbacus import skip_header


@click.command
@click.argument('slabs', nargs=-1)
@click.option('-o', 'out', required=True, help='Output ASDF filename', metavar='FN')
@click.option('-t', '--nthread', default=4, help='Number of Blosc compression threads')
@click.option('--delete', is_flag=True, help='Delete input files after processing')
def main(slabs, out, verbose=True, nthread=4, delete=False):
    slabs = [Path(s) for s in slabs]
    out = Path(out)
    if not slabs:
        return

    slice = slabs[0].parent
    ftype = out.name.split('.')[-2]
    cksum_reader = CksumReader(slice / 'checksums.crc32', verbose=verbose)

    out.parent.mkdir(parents=True, exist_ok=True)

    if ftype == 'pid':
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

    nptot = sum(fn.stat().st_size for fn in slabs)//pdtype().itemsize
    particles = np.empty(nptot, dtype=pdtype)

    i = 0
    for fn in slabs:
        with open(fn,'rb') as fp:
            newbytes = cksum_reader(fp)
            if ftype != 'pid':
                newbytesfp = io.BytesIO(newbytes)
                skip_header(newbytesfp)  # rvs have header
                newbytes = newbytes[newbytesfp.tell():]
                del newbytesfp
            new = np.frombuffer(newbytes, dtype=pdtype)
        particles[i:i+len(new)] = new
        i += len(new)
        del new, newbytes; gc.collect()  # noqa: E702

    if 'pid' in ftype:
        assert i == nptot
    else:
        assert i <= nptot
        particles = particles[:i]

    particles = particles.reshape(*pshape)

    compression_kwargs = dict(
        typesize=Mway,
        compression_block_size=asdf_block_size,
        blosc_block_size=blocksize,
        shuffle='bitshuffle',
        nthreads=nthread,
        )
    tree = {'data': {colname:particles},
            'header': dict(InputFile(slice / 'header'))}
    with asdf.AsdfFile(tree) as af, CksumWriter(out) as fp:
        af.write_to(fp, all_array_compression='blsc', compression_kwargs=compression_kwargs)
    print(f'Wrote to {out}')

    if delete:
        for fn in slabs:
            fn.unlink()


if __name__ == '__main__':
    main()
