'''
Output a disBatch task file to run group catalog post-processing.
'''

from pathlib import Path

import click

TASK = 'python -m Abacus.post.groups'


@click.command
@click.argument('groups', nargs=-1, metavar='GROUPDIR')
@click.option('-c', '--chunk', default=50, help='Target number of slabs per superslab')
@click.option('--delete', is_flag=True, help='Delete input files on success')
@click.option('-w', '--nworkers', default=1, help='Number of subprocesses per task')
@click.option('-l', '--logdir', default='post_log', help='Log directory')
def main(groups, chunk=10, nworkers=1, delete=False, logdir='post_log'):
    groups = [Path(s) for s in groups]

    flags = [f'-c {chunk} -w {nworkers}']
    if delete:
        flags += ['--delete']
    flags = ' '.join(flags)

    print(f'#DISBATCH PREFIX logdir={logdir}; (mkdir -p $logdir; {TASK} {flags} ')  # N.B. trailing space
    print('#DISBATCH SUFFIX  ) &> $logdir/group-$DISBATCH_TASKID.log')

    for group in groups:
        print(group.resolve())


if __name__ == '__main__':
    main()
