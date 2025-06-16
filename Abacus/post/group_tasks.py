# Copyright 2012-2025 The Abacus Developers
# SPDX-License-Identifier: GPL-3.0-or-later

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
@click.option('-l', '--logdir', default='post_log', help='Log directory', metavar='DIR')
def main(groups, chunk=10, nworkers=1, delete=False, logdir='post_log'):
    groups = [Path(s).absolute() for s in groups]
    logdir = Path(logdir).absolute()

    flags = [f'-c {chunk} -w {nworkers}']
    if delete:
        flags += ['--delete']
    flags = ' '.join(flags)

    print(f'#DISBATCH PREFIX {TASK} {flags} ')  # N.B. trailing space
    print(f'#DISBATCH SUFFIX  &> {logdir}/group-$DISBATCH_TASKID.log')

    for group in groups:
        print(group.absolute())


if __name__ == '__main__':
    main()
