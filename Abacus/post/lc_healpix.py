#!/usr/bin/env python
"""
Abacus outputs sparse healpix, which consists of an array of structs like:

{pixel_ID: uint32, N: uint32, {weight: float32, ...}}

The simulation writes these to the `lightcone/Step*/LightCone*_healsparse` files.

This script takes these structs and turns them into a regular healpix map.

TODO: we may also want to add back a mode similar to the old healpix ID outputs
that basically just stores the (sorted) structs in compressed ASDF.
"""

from pathlib import Path

import asdf
import click
import numba
import numpy as np
from timeit import default_timer as timer

from Abacus.fast_cksum.cksum_io import CksumReader, CksumWriter
from Abacus.InputFile import InputFile


def get_healstruct_dtype(HealpixWeightScheme):
    if HealpixWeightScheme == 'vel':
        return np.dtype(
            [
                ('id', np.uint32),
                ('N', np.uint32),
                ('weights', np.float32, 1),
            ],
            align=True,
        )
    elif HealpixWeightScheme == 'vel3':
        return np.dtype(
            [
                ('id', np.uint32),
                ('N', np.uint32),
                ('weights', np.float32, 3),
            ],
            align=True,
        )


def validate_files(files):
    # Check that all files have the same filename prefix, like LightCone0
    prefix = files[0].name.split('_')[0]
    for f in files[1:]:
        if f.name.split('_')[0] != prefix:
            raise ValueError(
                f'All files must have the same filename prefix. Got {f}, expected {files[0]}.'
            )
    return prefix


def get_step_range(headers) -> (int, int):
    steps = [h['FullStepNumber'] for h in headers]
    return min(steps), max(steps)


def read_structs(files: list[Path], healstruct_dtype):
    nstruct = sum(f.stat().st_size for f in files) // healstruct_dtype.itemsize

    structs = np.empty(nstruct, dtype=healstruct_dtype)

    i = 0
    for fn in files:
        cksum_reader = CksumReader(fn.parent / 'checksums.crc32')
        new = np.frombuffer(cksum_reader(fn), dtype=healstruct_dtype)
        structs[i : i + len(new)] = new
        i += len(new)
        del new
    assert i == nstruct

    return structs


@numba.njit(parallel=False, fastmath=True)
def make_map(structs, Nside, maskbits):
    npix = 12 * Nside**2
    density = np.zeros(npix, dtype=np.uint32)
    n_weights = structs['weights'].shape[-1]

    # use float64 to try to avoid underflow in accumulations
    wmaps = [np.zeros(npix, dtype=np.float64) for w in range(n_weights)]

    nstruct = len(structs)

    for i in range(nstruct):
        healid = structs[i]['id']
        density[healid] += structs[i]['N']

        for j in range(n_weights):
            wmaps[j][healid] += structs[i]['weights'][j]

    # TODO: density in float or int? ints compress better
    # density = density.astype(np.float32)
    wmaps = [w.astype(np.float32) for w in wmaps]

    if maskbits:
        # healview = density.view(np.uint32)
        mask = np.uint32(~np.uint32((1 << maskbits) - 1))
        # np.bitwise_and(healview, mask, healview)

        for w in wmaps:
            wview = w.view(np.uint32)
            np.bitwise_and(wview, mask, wview)

    return density, wmaps


def write_to_asdf(healpix, lc_prefix, quantity, step_range, outdir, meta, nthread=4):
    if step_range[0] == step_range[1]:
        steptag = f'Step{step_range[0]:04d}'
    else:
        steptag = 'Step{:04d}-{:04d}'.format(*step_range)

    outfn = outdir / f'heal-{quantity}' / f'{lc_prefix}_heal-{quantity}_{steptag}.asdf'
    meta = meta | {'LCHealpixMapQuantity': quantity}
    tree = meta | {
        'data': {f'heal-{quantity}': healpix},
    }

    # shuffle seems to compress slightly better than bitshuffle
    compression_kwargs = dict(
        typesize=4,
        compression_block_size=12 * 1024**2,
        blosc_block_size=3 * 1024**2,
        shuffle='shuffle',
        nthreads=nthread,
    )

    with asdf.AsdfFile(tree) as af, CksumWriter(outfn) as fp:
        af.write_to(
            fp, all_array_compression='blsc', compression_kwargs=compression_kwargs
        )
    print(f'Wrote to {outfn}')
    print(
        f'IO took {fp.io_timer.elapsed:.4g} sec, {fp.bytes_written/fp.io_timer.elapsed/1e6:.4g} MB/sec'
    )
    print(
        f'Checksum took {fp.checksum_timer.elapsed:.4g} sec, {fp.bytes_written/fp.checksum_timer.elapsed/1e6:.4g} MB/sec'
    )


@click.command
@click.argument('files', nargs=-1)
@click.option('-o', 'outdir', required=True, help='Output directory', metavar='DIR')
@click.option('-t', '--nthread', default=4, help='Number of Blosc compression threads')
@click.option('--delete', is_flag=True, help='Delete input files after processing')
@click.option(
    '-b',
    '--bit-trunc',
    default=0,
    help='Number of low bits to truncate in the float32 maps',
)
@click.option(
    '-d',
    '--dtype',
    default='f4',
    help='Output data type for the float maps',
)
@click.option(
    '-s',
    '--split-halo',
    is_flag=True,
    help='Interpret the high bit of the healpix ID as indicating a halo pixel',
)
def main(files, outdir, verbose=True, nthread=4, delete=False, bit_trunc=0, dtype='f4', split_halo=True):
    files = [Path(f) for f in files]
    if not files:
        return
    outdir = Path(outdir)

    outdir.mkdir(parents=True, exist_ok=True)

    dtype = np.dtype(dtype)

    lc_prefix = validate_files(files)  # LightCone0

    headers = [dict(InputFile(f.parent / 'header')) for f in files]
    healpix_weight_scheme = headers[0]['HealpixWeightScheme']
    healstruct_dtype = get_healstruct_dtype(healpix_weight_scheme)

    quantities = ['density', 'momentum-los']
    if healpix_weight_scheme == 'vel3':
        quantities += ['momentum-theta', 'momentum-phi']
    for q in quantities:
        (outdir / (f'heal-{q}')).mkdir(exist_ok=True)

    structs = read_structs(files, healstruct_dtype)
    step_range = get_step_range(headers)

    t = -timer()
    iord = np.argsort(structs['id'])
    structs = structs.take(iord)
    t += timer()
    print(f'Sort took {t:.4g} sec, {len(structs)/t/1e6:.4g} Mstruct/sec')

    print(f'Max N: {structs["N"].max()}')

    print(structs)

    t = -timer()
    # wouldn't be hard to downgrade nside here if the user wanted
    nside = headers[0]['LCHealpixNside']
    density, wmaps = make_map(structs, nside, maskbits=bit_trunc)

    # numba doesn't support np.float16, so do it here
    # density = density.astype(dtype, copy=False)
    wmaps = [w.astype(dtype, copy=False) for w in wmaps]
    
    for m in [density] + wmaps:
        assert np.all(np.isfinite(m))
    
    t += timer()
    print(f'Map took {t:.4g} sec, {len(structs)/t/1e6:.4g} Mstruct/sec')

    print(density)

    meta = {
        'headers': headers,
        'input_files': [str(f) for f in files],
        'bit_trunc': bit_trunc,
        'healpix_order': 'nest',
    }

    write_to_asdf(
        density, lc_prefix, 'density', step_range, outdir, meta, nthread=nthread
    )
    write_to_asdf(
        wmaps[0], lc_prefix, 'momentum-los', step_range, outdir, meta, nthread=nthread
    )
    
    if healpix_weight_scheme == 'vel3':
        write_to_asdf(
            wmaps[1], lc_prefix, 'momentum-theta', step_range, outdir, meta, nthread=nthread
        )
        write_to_asdf(
            wmaps[2], lc_prefix, 'momentum-phi', step_range, outdir, meta, nthread=nthread
        )

    if delete:
        for fn in files:
            fn.unlink()

if __name__ == '__main__':
    main()
