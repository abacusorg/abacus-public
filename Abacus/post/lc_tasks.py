'''
Output a disBatch task file to run lightcone post-processing.
'''

from pathlib import Path

import click

TASK = 'python -m Abacus.post.lc'


def group_by_prefix(fns: list[Path]) -> dict[str, list[Path]]:
    groups = {}
    for fn in fns:
        prefix = fn.name.split('.')[0]
        if prefix not in groups:
            groups[prefix] = []
        groups[prefix] += [fn]
    return groups


@click.command
@click.argument('lcdir')
@click.option('-l', '--logdir', default='post_log', help='Log directory', metavar='DIR')
@click.option('--delete', is_flag=True, help='Delete input files after processing')
def main(lcdir, logdir='post_log', delete=False):
    lcdir = Path(lcdir)

    for stepdir in sorted(lcdir.glob('Step*')):

        print(f'#DISBATCH PREFIX logdir={logdir}; (cd {stepdir}; mkdir -p $logdir; {TASK} {"--delete" if delete else ""} ')  # N.B. trailing space
        print('#DISBATCH SUFFIX  ) &> $logdir/lc-$DISBATCH_TASKID.log')

        kinds = ['heal', 'pid', 'rv']

        fns = stepdir.glob('LightCone*')
        fns = [f.relative_to(stepdir) for f in fns]
        fns = group_by_prefix(fns)

        outroot = lcdir.parent.resolve(strict=True) / 'lightcone-post'

        for prefix in fns:
            kind = prefix.split('_')[1].split('.')[0]
            if kind not in kinds:
                raise ValueError(f'Unknown file type {kind} in {prefix}')
            outfn = outroot / kind / f'{prefix}_{stepdir.name}.asdf'

            print(' '.join(['-o', str(outfn)] + [str(f) for f in fns[prefix]]))


if __name__ == '__main__':
    main()
