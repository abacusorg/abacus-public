"""
Output a disBatch task file to run lightcone post-processing.

We have two kinds of tasks to issue: one for rv/pid and one for healpix.
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


def decide_shells_simple(nconcat: int, headers: list[dict], plot: bool = False):
    # Just concatenate steps together in groups of nconcat.

    headers = headers[::-1]

    stepnums = np.array([h['FullStepNumber'] for h in headers])
    zmap = {step: h['Redshift'] for step, h in zip(stepnums, headers)}

    groups = [stepnums[i : i + nconcat] for i in range(0, len(stepnums), nconcat)]

    assert all(len(g) > 0 for g in groups)

    # check that every index was used once
    assert sorted(np.concatenate(groups)) == sorted(stepnums)

    if plot:
        widths = np.array(
            [h['FirstHalfEtaKick'] + h['LastHalfEtaKick'] for h in headers]
        )
        widths *= 2997.92  # Mpc/h
        group_widths = [
            sum(widths[(stepnums == s).argmax()] for s in g) for g in groups
        ]
        plot_widths(groups, group_widths, stepnums, widths, zmap)

    print(f'Grouped steps into {len(groups)} groups', file=sys.stderr)

    return groups


def plot_widths(groups, group_widths, orig_steps, orig_widths, zmap):
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots()
    ax.plot([g[0] for g in groups], group_widths, label='concatenated')
    ax.plot(orig_steps, orig_widths, label='raw')
    ax.tick_params(axis='both', right=True, top=True)

    ax.set_xlim(orig_steps[-1], orig_steps[0])

    # add top axis with redshift
    ax2 = ax.twiny()
    ax2.set_xlim(ax.get_xlim())
    bottomticks = ax.get_xticks()
    print(bottomticks)
    ax2.set_xticks(bottomticks)
    ax2.set_xticklabels(
        [f'{zmap[int(s)]:.2f}' if int(s) in zmap else '' for s in bottomticks]
    )
    ax2.set_xlabel('Redshift')

    ax.set_xlabel('Step number')
    ax.set_ylabel('Shell width [Mpc/h]')
    ax.legend()
    fig.savefig('shell_widths.png')


@click.command
@click.argument('lcdir')
@click.option('-l', '--logdir', default='post_log', help='Log directory', metavar='DIR')
@click.option('--delete', is_flag=True, help='Delete input files after processing')
@click.option(
    '-b',
    '--bit-trunc',
    default=0,
    help='Number of low bits to truncate in the float healpix maps',
)
@click.option(
    '-m',
    '--min-width',
    default=None,
    help='Minimum width of each shell in Mpc/h',
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
def main(
    lcdir,
    logdir='post_log',
    delete=False,
    bit_trunc=0,
    dtype='f4',
    min_width=None,
    nconcat=None,
    plot=False,
    split_halo=True,
):
    lcdir = Path(lcdir).absolute()
    logdir = Path(logdir).absolute()

    if (min_width is None) == (nconcat is None):
        raise ValueError('Must specify exactly one of --min-width or --nconcat')

    stepdirs = sorted(lcdir.glob('Step*'))
    headers = [dict(InputFile(s / 'header')) for s in stepdirs]
    
    if nconcat is None:
        nconcat = int(min_width / (headers[-1]['FirstHalfEtaKick'] + headers[-1]['LastHalfEtaKick'])) + 1
        print(f'Using nconcat={nconcat} to get min_width={min_width}', file=sys.stderr)
    
    shells = decide_shells_simple(nconcat, headers, plot=plot)

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
                task = f'{HEALPIX_TASK} -b {bit_trunc} -d {dtype}'
                if split_halo:
                    task += ' -s'
            else:
                raise ValueError(f'Unknown kind: {kind}')

            print(' '.join([task, '-o', str(out)] + [str(f) for f in fns[prefix]]))


if __name__ == '__main__':
    main()
