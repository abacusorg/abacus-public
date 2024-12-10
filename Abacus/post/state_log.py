"""
Gather the `log/step*/*.state` files and save them into a single ASDF file.

Usually users are encouraged to get the state from the header of the halo catalog/timeslice/etc.
But occasionally one wants to know the state for every time step, e.g. to examine the
time stepping. It's much faster to have that all in a single file, which this script builds.

One might also want to interpolate chi(z) or other functions from the state. This something of
an anti-pattern, since querying the Abacus.Cosmology module will usually be more accurate.
But maybe one wants to do this offline or not worry about discrepancies between the simulation's
state and results from Abacus.Cosmology.
"""

from pathlib import Path

import asdf
import click
from astropy.table import Table

from Abacus.InputFile import InputFile

@click.command
@click.argument('logdir', type=click.Path(exists=True, file_okay=False, dir_okay=True, path_type=Path))
def main(logdir):
    logdir = Path(logdir)
    state_files = sorted(logdir.glob('step*/*.state'))

    if not state_files:
        raise FileNotFoundError('No state files found in {}'.format(logdir / 'step*'))

    states = []
    for state_file in state_files:
        state = dict(InputFile(state_file))
        # An empty string gets interpreted by InputFile as an empty tuple...
        if state.get('GroupFindingDensitySource') == tuple():
            state['GroupFindingDensitySource'] = ""
        states.append(state)

    states = Table(states)

    outfn = logdir.parent / 'state_log.asdf'
    with asdf.AsdfFile({'data': states}) as af:
        af.write_to(outfn, all_array_compression='blsc')
    
    print(f'Wrote state log to {outfn}')


if __name__ == '__main__':
    main()
