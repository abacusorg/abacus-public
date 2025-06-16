# Copyright 2012-2025 The Abacus Developers
# SPDX-License-Identifier: GPL-3.0-or-later

'''
Output a disBatch task file to run timeslice post-processing.
'''

import warnings
from pathlib import Path

import click
import numpy as np

from Abacus.InputFile import InputFile

TASK = 'python -m Abacus.post.timeslice'

@click.command
@click.argument('slices', nargs=-1)
@click.option('-c', '--chunk', default=10, help='Target number of slabs per superslab')
@click.option('-l', '--logdir', default='post_log', help='Log directory', metavar='DIR')
@click.option('--delete', is_flag=True, help='Delete input files after processing')
def main(slices, chunk=10, logdir='post_log', delete=False):
    slices = [Path(s) for s in slices]

    for slice in slices:
        workdir = slice.parent

        # N.B. trailing space
        print(f'#DISBATCH PREFIX logdir={logdir}; (cd {workdir}; mkdir -p $logdir; {TASK} {"--delete" if delete else ""} ')
        print('#DISBATCH SUFFIX  ) &> $logdir/slice-$DISBATCH_TASKID.log')

        header = dict(InputFile(slice / 'header'))
        cpd = header['CPD']
        NumZRanks = int(header.get('NumZRanks',1))
        n_subslab = cpd*NumZRanks

        kinds = ['L0_pack9', 'field_pack9', 'L0_pack9_pids', 'field_pack9_pids']

        for kind in kinds:
            fns = sorted(slice.glob(f'*.{kind}.dat'))
            fns = [f.relative_to(workdir) for f in fns]
            
            kind_dir = kind[:-1] if 'pid' in kind else kind
            outdir = Path('slices') / slice.name.replace("slice","z") / kind_dir
            
            if len(fns) == 0:
                warnings.warn(
                    f'Skipping {kind} because no files found for "{slice}"',
                    )
                continue

            if len(fns) != n_subslab:
                # TODO: technically this could be a warning instead
                # but we'd need to handle variable-length slabs-per-superslab
                raise ValueError(
                    f'Expected {n_subslab} files, found {len(fns)} for {kind=}, "{slice}"',
                    )

            # build the name, like "slab000.field.pack9.asdf"
            w = 3  # len(str(cpd//chunk))
            name = [
                f'slab{{:0{w}d}}',
                kind.split('_')[0],
                'pack9',
            ]
            if 'pid' in kind:
                name += ['pid']
            name += ['asdf']
            name = '.'.join(name)

            bounds = NumZRanks*np.linspace(0, cpd, cpd//chunk + 1, dtype=int)
            for i in range(len(bounds)-1):
                print('-o ' + str(outdir / name.format(i)) + ' ' +
                      ' '.join([str(f) for f in fns[bounds[i]:bounds[i+1]]])
                      )


if __name__ == '__main__':
    main()
