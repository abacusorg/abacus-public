# Copyright 2012-2025 The Abacus Developers
# SPDX-License-Identifier: GPL-3.0-or-later

'''
Output a disBatch task file to run power spectrum computation.
'''

from pathlib import Path

import click

# TODO: migrate to abacusutils power
TASK = 'python -m Abacus.Analysis.PowerSpectrum.run_PS'


@click.command
@click.argument('inputs', nargs=-1, metavar='INPUT')
@click.option('--nfft', default=256, help='Density grid size')
@click.option('-t', '--nthread', default=1, help='Number of processing threads')
@click.option('-r', '--nreader', default=1, help='Number of IO threads')
@click.option('-l', '--logdir', default='post_log', help='Log directory', metavar='DIR')
def main(inputs, nfft=256, nthread=1, nreader=1, logdir='post_log'):
    inputs = [Path(i).resolve(strict=True) for i in inputs]

    flags = f'--nfft {nfft} --nthreads {nthread} --nreaders {nreader} --nbins {nfft}'

    print(f'#DISBATCH PREFIX logdir={logdir}; (mkdir -p $logdir; {TASK} {flags} ')  # N.B. trailing space
    print('#DISBATCH SUFFIX  ) &> $logdir/power-$DISBATCH_TASKID.log')

    for i in inputs:
        if 'halos' in i.parent.name:
            fmt = 'asdf'
            label = 'AB'
            if not any(True for _ in i.glob('*_rv_*')):
                continue  # no subsample particles
        elif 'slices' in i.parent.name:
            fmt = 'asdf_pack9'
            label = 'pack9'
        else:
            raise ValueError(f'Unknown file format for {i}')
        outdir = i.parents[1] / 'power' / i.name / label
        print(f'-o {outdir} --format {fmt} {i}')
        print(f'-o {outdir} --zspace --format {fmt} {i}')


if __name__ == '__main__':
    main()
