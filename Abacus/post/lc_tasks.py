"""
Output a disBatch task file to run lightcone post-processing.

We have two kinds of tasks to issue: one for rv/pid and one for healpix.

For the healpix tasks, this is where we decide how to group the individual
light cone shells (one per timestep) into larger shells.
"""

import sys
from pathlib import Path

import click
import numpy as np

from Abacus.InputFile import InputFile

PARTICLE_TASK = 'python -m Abacus.post.lc'
HEALPIX_TASK = 'python -m Abacus.post.lc_healpix'


def group_by_prefix(fns: list[Path]) -> dict[str, list[Path]]:
    groups = {}
    for fn in fns:  # LightCone0_rv.0
        prefix = fn.name.split('.')[0]  # LightCone0_rv
        if prefix not in groups:
            groups[prefix] = []
        groups[prefix] += [fn]
    return groups


def decide_shells_simple(
    nconcat: int,
    widths: list[float],
    stepnums: list[int],
    redshifts: list[float],
    plot: bool = False,
):
    # Just concatenate steps together in groups of nconcat.

    widths = widths[::-1]
    stepnums = stepnums[::-1]
    redshifts = redshifts[::-1]

    groups = [stepnums[i : i + nconcat] for i in range(0, len(stepnums), nconcat)]

    assert all(len(g) > 0 for g in groups)

    # check that every index was used once
    assert sorted(np.concatenate(groups)) == sorted(stepnums)

    if plot:
        zmap = dict(zip(stepnums, redshifts))
        group_widths = [
            sum(widths[(stepnums == s).argmax()] for s in g) for g in groups
        ]
        plot_widths(groups, group_widths, stepnums, widths, zmap)

    print(
        f'Grouped steps into {len(groups)} groups, mean width of {np.sum(widths) / len(groups):.3g}',
        file=sys.stderr,
    )

    return groups


def decide_shells_greedy(
    width: float,
    widths: list[float],
    stepnums: list[int],
    redshifts: list[float],
    plot: bool = False,
):
    # Group shells using a greedy algorithm

    widths = widths[::-1]
    stepnums = stepnums[::-1]
    redshifts = redshifts[::-1]

    groups = []
    current_group = []
    current_width = 0.0

    for i, (w, s) in enumerate(zip(widths, stepnums)):
        # Start a new group if more than half of the current step would spill the shell.
        if current_width + w > width + w / 2:
            # If we just overfilled a shell, then try to underfill the next shell by carrying
            # over the spillage. This trick borrowed from:
            # https://www.werkema.com/2021/11/01/an-efficient-solution-to-linear-partitioning/
            error = width - current_width

            groups.append(current_group)
            current_group = [s]
            current_width = w - error
        else:
            current_group.append(s)
            current_width += w

    if len(current_group) > 0:
        groups.append(current_group)

    assert all(len(g) > 0 for g in groups)
    # check that every index was used once
    assert sorted(np.concatenate(groups)) == sorted(stepnums)

    group_widths = [sum(widths[(stepnums == s).argmax()] for s in g) for g in groups]

    if plot:
        zmap = dict(zip(stepnums, redshifts))
        plot_widths(groups, group_widths, stepnums, widths, zmap)
    print(
        f'Grouped steps into {len(groups)} groups, '
        f'mean width {np.mean(group_widths):.3g} '
        f'± {np.std(group_widths):.3g} Mpc/h',
        file=sys.stderr,
    )
    return groups


def plot_widths(groups, group_widths, orig_steps, orig_widths, zmap):
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(constrained_layout=True)
    ax.plot([g[0] for g in groups], group_widths, label='concatenated')
    ax.plot(orig_steps, orig_widths, label='raw')
    ax.tick_params(axis='both', right=True, top=True)

    ax.set_xlim(orig_steps[-1], orig_steps[0])

    ax2 = ax.twiny()
    ax2.set_xlim(ax.get_xlim())
    bottomticks = ax.get_xticks()
    ax2.set_xticks(bottomticks)
    ax2.set_xticklabels(
        [f'{zmap[int(s)]:.2f}' if int(s) in zmap else '' for s in bottomticks]
    )
    ax2.set_xlabel('Redshift')

    ax.set_xlabel('Step number')
    ax.set_ylabel('Shell width [Mpc/h]')
    ax.legend()
    fig.savefig(fn := 'shell_widths.png')
    print(f'Saved plot of shell widths to {fn}', file=sys.stderr)


def get_steps_and_headers(lcdir: Path, discard_partialsky: bool = False):
    stepdirs = sorted(lcdir.glob('Step*'))

    # Get the headers that bound all the steps.
    bounding_headers = [
        (dict(InputFile(s / 'header_read')), dict(InputFile(s / 'header_write')))
        for s in stepdirs
    ]

    if discard_partialsky:
        rmax = bounding_headers[0][0]['BoxSize'] * (
            bounding_headers[0][0]['LCBoxRepeats'] + 0.5
        )

        # If the outer edge of the step is within the volume, keep it.
        mask = [
            hpair[0]['CoordinateDistanceHMpc'] <= rmax for hpair in bounding_headers
        ]
        stepdirs = [s for (s, keep) in zip(stepdirs, mask) if keep]
        bounding_headers = [b for (b, keep) in zip(bounding_headers, mask) if keep]

    return stepdirs, bounding_headers


@click.command
@click.argument('lcdir')
@click.option('-l', '--logdir', default='post_log', help='Log directory', metavar='DIR')
@click.option(
    '--delete', is_flag=True, help='Have the tasks delete input files after processing'
)
@click.option(
    '-b',
    '--bit-trunc',
    default=0,
    help='Number of low bits to truncate in the float healpix maps',
)
@click.option(
    '-w',
    '--width',
    default=None,
    help='Target width of each shell in Mpc/h',
    type=float,
)
@click.option(
    '-n', '--nconcat', default=None, help='Number of steps to concatenate', type=int
)
@click.option('-p', '--plot', is_flag=True, help='Plot the shell widths')
@click.option(
    '-d',
    '--dtype',
    default='f4',
    help='Output data type for the float healpix maps',
)
@click.option(
    '-s',
    '--split-halo',
    is_flag=True,
    help='Output a healpix map with just the halo pixels',
)
@click.option(
    '--coarsify-trans-vel',
    default=1,
    help='Coarsify the transverse velocity healpix map Nside by this factor',
)
@click.option(
    '-o',
    '--out',
    'outroot',
    type=click.Path(file_okay=False, dir_okay=True, path_type=Path),
    default=None,
    help='Output directory for the healpix maps. Default: LCDIR/../lightcone-post/',
)
@click.option(
    '--discard-partialsky',
    is_flag=True,
    help='Discard healpix maps at high enough redshift that they do not cover the full sky',
)
@click.option(
    '-a',
    '--algorithm',
    default='greedy',
    type=click.Choice(['simple', 'greedy']),
    help='Algorithm to use for grouping shells',
)
def main(
    lcdir,
    logdir='post_log',
    delete=False,
    bit_trunc=0,
    dtype='f4',
    width=None,
    nconcat=None,
    plot=False,
    split_halo=True,
    coarsify_trans_vel=1,
    outroot=None,
    discard_partialsky=False,
    algorithm='greedy',
):
    lcdir = Path(lcdir).absolute()
    logdir = Path(logdir).absolute()

    if (width is None) == (nconcat is None):
        raise ValueError('Must specify exactly one of --width or --nconcat')

    stepdirs, bounding_headers = get_steps_and_headers(
        lcdir, discard_partialsky=discard_partialsky
    )

    widths = np.array(
        [
            h[0]['CoordinateDistanceHMpc'] - h[1]['CoordinateDistanceHMpc']
            for h in bounding_headers
        ]
    )
    stepnums = np.array([h[1]['FullStepNumber'] for h in bounding_headers], dtype=int)
    redshifts = np.array(
        [h[1]['Redshift'] for h in bounding_headers]
    )  # z at end of step, currently just cosmetic

    if algorithm == 'simple':
        if nconcat is None:
            nconcat = int(np.ceil(width / np.mean(widths)))
        shells = decide_shells_simple(nconcat, widths, stepnums, redshifts, plot=plot)
    elif algorithm == 'greedy':
        if width is None:
            width = np.mean(widths) * nconcat
        shells = decide_shells_greedy(width, widths, stepnums, redshifts, plot=plot)

    for steps in shells:
        stepdirs = [lcdir / f'Step{step:04d}' for step in steps]

        if len(steps) == 1:
            steptag = f'Step{steps[0]:04d}'
        else:
            steptag = f'Step{steps[0]:04d}-{steps[-1]:04d}'

        fns = sum((list(stepdir.glob('LightCone*')) for stepdir in stepdirs), [])
        fns = sorted(f.relative_to(lcdir) for f in fns)
        fns = group_by_prefix(fns)

        print(f'#DISBATCH PREFIX (cd {lcdir}; ')
        print(
            f'#DISBATCH SUFFIX  {"--delete" if delete else ""} ) &> {logdir}/lc-$DISBATCH_TASKID.log'
        )

        if outroot is None:
            outroot = lcdir.parent.resolve(strict=True) / 'lightcone-post'

        for prefix in fns:  # LightCone0_rv
            kind = prefix.split('_')[1]  # rv
            if kind in ('rv', 'pid'):
                # rv/pid accepts an output filename
                out = (
                    outroot / kind / f'{prefix}_{steptag}.asdf'
                )  # SimName/lightcone-post/rv/LightCone0_rv.Step0000.asdf
                task = PARTICLE_TASK
            elif kind == 'healsparse':
                # healpix accepts an output directory because it writes multiple files
                out = outroot
                task = f'{HEALPIX_TASK} -b {bit_trunc} -d {dtype} --coarsify-trans-vel {coarsify_trans_vel}'
                if split_halo:
                    task += ' -s'
            else:
                raise ValueError(f'Unknown kind: {kind}')

            print(' '.join([task, '-o', str(out)] + [str(f) for f in fns[prefix]]))


if __name__ == '__main__':
    main()
