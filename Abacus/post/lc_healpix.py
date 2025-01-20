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
from timeit import default_timer as timer

import asdf
import click
import numba
import numpy as np

from Abacus.fast_cksum.cksum_io import CksumReader, CksumWriter
from Abacus.InputFile import InputFile

# The bit that indicates a halo pixel
HALO_MASK = np.uint32(0x80000000)

print_ = print
def print(*args, **kwargs):
    kwargs['flush'] = True
    return print_(*args, **kwargs)

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


@numba.njit(parallel=False, fastmath=True)
def _compress_helper(structs):
    nstruct = len(structs)
    j = 0
    for i in range(1, nstruct):
        if structs[i]['id'] == structs[j]['id']:
            structs[j]['N'] += structs[i]['N']
            structs[j]['weights'] += structs[i]['weights']
        else:
            j += 1
            structs[j] = structs[i]

    return j + 1


def compress_healstructs(structs):
    # Sort by healid, co-add identical pixels towards the front.
    structs = structs.take(np.argsort(structs['id']))
    newsize = _compress_helper(structs)
    return structs[:newsize]


def read_structs(files: list[Path], healstruct_dtype):
    nstruct = sum(f.stat().st_size for f in files) // healstruct_dtype.itemsize

    structs = np.empty(nstruct, dtype=healstruct_dtype)

    tio = 0
    tcompress = 0


    i = 0
    # TODO: we could parallelize this, similar to pcopy
    # TODO: if the memory due to IO is the bottleneck, we could read in chunks
    # and compress as we go
    for fn in files:
        cksum_reader = CksumReader(fn.parent / 'checksums.crc32')

        tio -= timer()
        new = np.frombuffer(cksum_reader(fn), dtype=healstruct_dtype)
        tio += timer()

        tcompress -= timer()
        new = compress_healstructs(new)
        tcompress += timer()

        structs[i : i + len(new)] = new
        i += len(new)
        del new
    assert i <= nstruct

    structs = structs[:i]

    # Could do a final compress here, which would also allow us to
    # parallelize make_map

    print(f'IO took {tio:.4g} sec, {nstruct * healstruct_dtype.itemsize / tio / 1e6:.4g} MB/sec')
    print(f'Compression took {tcompress:.4g} sec, {nstruct / tcompress / 1e6:.4g} Mstruct/sec')

    return structs


@numba.njit(parallel=False, fastmath=True)
def make_map(structs, Nside, maskbits, discard_halo_bit):
    npix = 12 * Nside**2
    density = np.zeros(npix, dtype=np.uint32)
    n_weights = structs['weights'].shape[-1]

    # use float64 to try to avoid underflow in accumulations
    wmaps = [np.zeros(npix, dtype=np.float64) for w in range(n_weights)]

    nstruct = len(structs)

    for i in range(nstruct):
        healid = structs[i]['id']
        if discard_halo_bit:
            healid &= ~HALO_MASK
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

    outfn.parent.mkdir(parents=True, exist_ok=True)
    with asdf.AsdfFile(tree) as af, CksumWriter(outfn) as fp:
        af.write_to(
            fp, all_array_compression='blsc', compression_kwargs=compression_kwargs
        )
    print(f'Wrote to {outfn}')
    print(
        f'IO took {fp.io_timer.elapsed:.4g} sec, {fp.bytes_written / fp.io_timer.elapsed / 1e6:.4g} MB/sec'
    )
    print(
        f'Checksum took {fp.checksum_timer.elapsed:.4g} sec, {fp.bytes_written / fp.checksum_timer.elapsed / 1e6:.4g} MB/sec'
    )


@numba.njit(parallel=False, fastmath=True)
def searchsorted_bit(a, mask):
    """Search a sorted array of integers for the first element a[i] where mask & a[i] != 0"""

    lo = 0
    hi = len(a)

    while lo < hi:
        mid = (lo + hi) // 2
        if a[mid] & mask:
            hi = mid
        else:
            lo = mid + 1

    return lo


def process_structs(
    structs,
    headers,
    bit_trunc,
    dtype,
    files,
    totalhalo,
    outdir,
    which_th,
    lc_prefix,
    step_range,
    nthread,
    healpix_weight_scheme,
):
    t = -timer()
    # wouldn't be hard to downgrade nside here if the user wanted
    nside = headers[0]['LCHealpixNside']

    # Currently, Abacus always outputs the halo bit
    discard_halo_bit = True

    density, wmaps = make_map(
        structs,
        nside,
        maskbits=bit_trunc,
        discard_halo_bit=discard_halo_bit,
    )

    # numba doesn't support np.float16, so do it here
    # density = density.astype(dtype, copy=False)
    wmaps = [w.astype(dtype, copy=False) for w in wmaps]

    for m in [density] + wmaps:
        assert np.all(np.isfinite(m))

    t += timer()
    print(f'Map took {t:.4g} sec, {len(structs) / t / 1e6:.4g} Mstruct/sec')

    print(density)

    meta = {
        'headers': headers,
        'input_files': [str(f) for f in files],
        'bit_trunc': bit_trunc,
        'healpix_order': 'nest',
    }

    if len(totalhalo) > 1:
        this_outdir = outdir / which_th
        this_lc_prefix = f'{lc_prefix}_{which_th}'
    else:
        this_outdir = outdir
        this_lc_prefix = lc_prefix

    maps = {'density': density, 'momentum-los': wmaps[0]}
    if healpix_weight_scheme == 'vel3':
        maps |= {'momentum-theta': wmaps[1], 'momentum-phi': wmaps[2]}

    for quantity, mp in maps.items():
        write_to_asdf(
            mp,
            this_lc_prefix,
            quantity,
            step_range,
            this_outdir,
            meta,
            nthread=nthread,
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
    help='Output a healpix map with just the halo pixels',
)
def main(
    files,
    outdir,
    verbose=True,
    nthread=4,
    delete=False,
    bit_trunc=0,
    dtype='f4',
    split_halo=True,
):
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

    structs = read_structs(files, healstruct_dtype)
    N_struct = len(structs)

    step_range = get_step_range(headers)

    t = -timer()
    structs = structs.take(np.argsort(structs['id']))
    t += timer()
    print(f'Sort took {t:.4g} sec, {N_struct / t / 1e6:.4g} Mstruct/sec')

    print(f'Max N: {structs["N"].max()}')

    print(structs)

    totalhalo = ['total']
    if split_halo:
        totalhalo += ['halo']

    for which_th in totalhalo:
        if which_th == 'total':
            start = 0
        else:
            start = searchsorted_bit(structs['id'], HALO_MASK)
            print(
                f'Found {N_struct - start} halo pixels, {(N_struct - start) / N_struct * 100:.4f}%'
            )

        # TODO: calling a function with so many args is unsatisfying, but the other solutions
        # (object-oriented, dataclass config, dict config, etc) all have problems, too.
        process_structs(
            structs[start:],
            headers,
            bit_trunc,
            dtype,
            files,
            totalhalo,
            outdir,
            which_th,
            lc_prefix,
            step_range,
            nthread,
            healpix_weight_scheme,
        )

    if delete:
        for fn in files:
            fn.unlink()


if __name__ == '__main__':
    main()
