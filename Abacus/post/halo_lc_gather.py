# Copyright 2012-2025 The Abacus Developers
# SPDX-License-Identifier: GPL-3.0-or-later

"""
The halo light cone code selects halos from the halo time slice catalog and computes
their intersection with the light cone. In most cases, only some of these halos will
be present in the light cone, so to save space we'd like to discard any halos that
don't appear in the light cone.

The easiest thing to do would be to copy each halo catalog entry into the light cone
catalog, making a standalone data product. However, halos may be repeated many times,
since we tile the box to reach higher redshift. It would be wasteful to repeat the
same halo. So instead, we store the index of the halo in the time slice catalog,
plus some properties of each unique intersection, like the interpolated pos/vel/mass.

The indices refer to the full time slice catalog. But we can save space by getting
rid of any halos in the time slice catalog that aren't referenced in the light cone.
That is what this code does.

This core idea is pretty straightforward; most of the complexity comes from the
cleaning and particle data. We need to tack on the cleaning info, merge the particles
from the original and cleaned catalog, and update the particle indices.

---

Total LC halos: 7059596920
Total unique LC halos: 1322725109
Total snapshot halos: 2260853720

unique is 0.58x of total snapshot

596 GB * 0.58 = 345 GB
  + 152 GB (indirect)
  + ??? A rv
"""

import gc
import shutil
from builtins import print as builtin_print
from pathlib import Path
from timeit import default_timer as timer

import asdf
import click
import numpy as np
from abacusnbody.data.compaso_halo_catalog import CompaSOHaloCatalog
from astropy.table import Table

from Abacus.fast_cksum.cksum_io import CksumWriter

COMPRESSION_KWARGS = dict(
    typesize='auto',
    shuffle='bitshuffle',
    compression_block_size=12 * 1024**2,
    blosc_block_size=3 * 1024**2,
    nthreads=4,
)

DIR_ARG = click.Path(exists=True, file_okay=False, dir_okay=True, path_type=Path)


def print(*args, **kwargs):
    flush = kwargs.pop('flush', True)
    return builtin_print(*args, **kwargs, flush=flush)


@click.command
@click.argument('halo_lc_path', type=DIR_ARG)
@click.argument('sim_path', type=DIR_ARG)
@click.option(
    '--delete',
    '-d',
    is_flag=True,
    help='Delete the original halo light cone catalog after gathering',
)
def main(halo_lc_path, sim_path, delete=False):
    # halo_lc_path is the input
    # sim_path is the input for the time slice catalog, and the output for the gathered catalog

    sim_name = halo_lc_path.parent.name
    zdir = halo_lc_path.name

    assert sim_name == sim_path.name

    snapshot_halo_path = sim_path / 'halos' / zdir

    halo_info_fns = sorted(
        (snapshot_halo_path / 'halo_info').glob('halo_info_*.asdf'),
        key=lambda x: int(x.stem.split('_')[-1]),
    )
    N_superslab = len(halo_info_fns)

    for ss in range(N_superslab):
        # Load all observers into one table
        tread = -timer()
        lc_cat = read_many(halo_lc_path.glob(f'Merger_lc*.{ss:02d}.asdf'))

        # FUTURE: rename these upstream
        lc_cat.rename_columns(
            (
                'InterpolatedPosition',
                'InterpolatedVelocity',
                'InterpolatedN',
                'InterpolatedComovingDist',
                'HaloIndex',
            ),
            (
                'Interpolated_x_L2com',
                'Interpolated_v_L2com',
                'Interpolated_N',
                'Interpolated_ComovingDist',
                'halo_timeslice_index',
            ),
        )

        tread += timer()
        print(f'Read LC catalog in {tread:.4g}s')

        # The halo light cone code seems to be outputting cleaned halos, as well as marking
        # some halos ineligible for the light cone by setting their mass to zero. We'd probably
        # rather make these decisions this upstream, but for now we'll just remove
        # the halos here.
        len_before = len(lc_cat)
        lc_cat = lc_cat[lc_cat['Interpolated_N'] > 0]
        print(f'Removed {len_before - len(lc_cat)} halos with zero mass')

        # Sort on halo_timeslice_index
        tsort = -timer()
        lc_cat.sort('halo_timeslice_index')
        tsort += timer()
        print(f'Sorted LC catalog in {tsort:.4g}s')

        # Get unique halo indices
        tunique = -timer()
        # This could be accelerated because we know the indices are sorted,
        # but it's already fast
        indices = np.unique_values(lc_cat['halo_timeslice_index'])
        tunique += timer()
        print(f'Got unique indices in {tunique:.4g}s')

        # Load the snapshot halo catalog
        tread = -timer()

        # Create a boolean mask with the halo rows to load.
        # CHC will take care of selecting the particles corresponding to the loaded halos,
        # including cleaning
        with asdf.open(halo_info_fns[ss], lazy_load=True) as af:
            N_halo = len(af['data']['N'])
        mask = np.zeros(N_halo, dtype=np.bool)
        mask[indices] = True
        del indices

        snapshot_cat = CompaSOHaloCatalog(
            halo_info_fns[ss],
            subsamples={
                'A': True,
                'B': False,  # only needed for merger trees, not HOD
                'rvint': True,
                'packedpid': True,
                'unpack': False,
            },
            passthrough=True,
            fields='all',
            filter_func=lambda _: mask,
        )

        tread += timer()
        print(f'Read snapshot catalog in {tread:.4g}s')

        tfixup = -timer()

        # No B particles
        snapshot_cat.halos['npstartB'][:] = 0
        snapshot_cat.halos['npoutB'][:] = 0

        # Fix up the light cone catalog indices by making them contiguous, starting from 0.
        # The indices are sorted but not unique.
        lc_cat['halo_timeslice_index'][1:] = (np.diff(lc_cat['halo_timeslice_index']) > 0).cumsum()
        lc_cat['halo_timeslice_index'][0] = 0
        assert lc_cat['halo_timeslice_index'][-1] == len(snapshot_cat.halos) - 1

        tfixup += timer()
        print(f'Fixed up LC catalog in {tfixup:.4g}s')

        # Save the gathered catalog
        tsave = -timer()
        outfn = (
            sim_path / 'lightcone_halos' / zdir / f'lightcone_halo_info_{ss:03d}.asdf'
        )
        particleoutfn = outfn.parent / f'lightcone_rvpid_{ss:03d}.asdf'
        outfn.parent.mkdir(parents=True, exist_ok=True)

        tree = {
            'halo_lightcone': {
                k: np.asarray(v, copy=False) for (k, v) in lc_cat.items()
            },
            'halo_timeslice': {
                k: np.asarray(v, copy=False) for (k, v) in snapshot_cat.halos.items()
            },
            'header': dict(snapshot_cat.header),
            # "header_next" is not really useful, since the shell is centered on a timeslice,
            # rather than bounded by timeslices
            # 'header_next': dict(lc_cat.meta),
        }

        with CksumWriter(outfn) as fp, asdf.AsdfFile(tree) as af:
            af.write_to(
                outfn,
                all_array_compression='blsc',
                compression_kwargs=COMPRESSION_KWARGS,
            )
        del lc_cat, tree, snapshot_cat.halos, af
        gc.collect()

        particle_tree = {
            'data': {
                k: np.asarray(v, copy=False)
                for (k, v) in snapshot_cat.subsamples.items()
            },
            'header': dict(snapshot_cat.header),
        }

        # Save the particle data
        with CksumWriter(particleoutfn) as fp, asdf.AsdfFile(particle_tree) as af:
            af.write_to(
                fp,
                all_array_compression='blsc',
                compression_kwargs=COMPRESSION_KWARGS,
            )
        del snapshot_cat, particle_tree, af
        gc.collect()

        tsave += timer()
        print(f'Saved gathered catalog in {tsave:.4g}s')
        print(f'Wrote to {outfn}')

    if delete:
        shutil.rmtree(halo_lc_path)
        print(f'Deleted {halo_lc_path}')


def read_many(fns):
    """Read many asdf files into a single Astropy table, pre-allocating the columns."""

    # Open all the files and sum up the column lengths
    afs = [asdf.open(fn, memmap=False) for fn in fns]
    first_colname = next(iter(afs[0]['data']))
    Ntot = sum(len(af['data'][first_colname]) for af in afs)
    data0 = afs[0]['data']
    table = Table(
        {
            name: np.empty(
                (Ntot,) + data0[name].shape[1:],
                dtype=data0[name].dtype,
            )
            for name in data0
        },
        copy=False,
        meta=afs[0]['header'],
    )

    # Read the data
    i = 0
    for af in afs:
        n = len(af['data'][first_colname])
        for name in af['data']:
            table[name][i : i + n] = af['data'][name]
        i += n
    assert i == len(table)

    return table


if __name__ == '__main__':
    main()
