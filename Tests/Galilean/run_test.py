#!/usr/bin/env python3
'''
Set up a uniform grid with a constant velocity and run for a number
of time steps. Then check that we haven't lost any particles!
'''

from pathlib import Path
import argparse

import numpy as np
import numba as nb

from Abacus import abacus
from Abacus import InputFile

PID_MASK = np.uint64(0x7fff7fff7fff)

ci_dtype = np.dtype([('startindex',np.uint64),('startindex_with_ghost',np.uint64),
                    ('count',np.int32), ('active',np.int32),
                    ('mean_square_velocity', np.float32), ('max_component_velocity', np.float32),
                    ('max_component_acceleration', np.float32)],
                    align=True)

def load_ci(params, state):
    lwd = Path(params['LocalWorkingDirectory'])
    cifns = sorted(lwd.parent.glob(lwd.name + '*/read/cellinfo*'))
    
    zsplits = params['NumZRanks'] if params['NumZRanks'] > 1 else 0
    cpd = params['CPD']
    ghost_radius = state['GhostRadius']
    ncell = (cpd + 2*ghost_radius*zsplits)*cpd**2
    ci = []
    for fn in cifns:
        ci += [np.fromfile(fn, dtype=ci_dtype)]

    return ci


def load_aux(params, state):
    lwd = Path(params['LocalWorkingDirectory'])
    auxfns = sorted(lwd.parent.glob(lwd.name + '*/read/aux*'))
    
    aux = np.empty(state['np_with_ghost_state'], dtype=np.uint64)
    i = 0
    slabstart = np.empty(params['CPD']*params['NumZRanks']+1, dtype=np.uint64)
    slabstart[0] = 0
    for x,fn in enumerate(auxfns):
        _aux = np.fromfile(fn, dtype=np.uint64)
        aux[i:i+len(_aux)] = _aux
        i += len(_aux)
        slabstart[x+1] = i
    assert i == len(aux)

    return aux, slabstart


def select_primaries(aux, slabstart, ci, state):
    NP = state['np_state']
    cpd = state['cpd_state']
    ghost_radius = state['GhostRadius']

    primaries = np.empty(NP, dtype=np.uint64)

    i = 0
    for x in range(len(ci)):
        primary_ci = ci[x].reshape(cpd, -1)[:,ghost_radius:-ghost_radius]
        for y in range(len(primary_ci)):
            rowstart = primary_ci[y,0]['startindex_with_ghost'] + slabstart[x]
            rowend = primary_ci[y,-1]['startindex_with_ghost'] + primary_ci[y,-1]['count'].astype(np.uint64) + slabstart[x]
            _p = aux[rowstart : rowend]
            primaries[i:i+len(_p)] = _p
            i += len(_p)
    assert i == NP

    return primaries


@nb.njit(parallel=True)
def check_pid(pid, ppd):
    nbad = 0
    for i in nb.prange(ppd):
        x = i << 32
        for j in range(ppd):
            y = j << 16
            for k in range(ppd):
                z = k
                p = x | y | z
                if pid[i*ppd*ppd + j*ppd + k] != p:
                    nbad += 1  # automatic reduction
    return nbad


def main(nstep):
    abacus.run('abacus.par2', maxsteps=nstep)

    params = InputFile.InputFile('abacus.par')
    state = InputFile.InputFile(Path(params['WorkingDirectory'])/'read'/'state')

    ci = load_ci(params, state)
    aux, slabstart = load_aux(params, state)
    primaries = select_primaries(aux, slabstart, ci, state)

    primary_pid = primaries & PID_MASK
    primary_pid.sort()
    ppd = int(round(state['ppd']))

    print(check_pid(primary_pid, ppd))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('nstep', help='Number of time steps', type=int, default=10, nargs='?')
    
    args = vars(parser.parse_args())
    main(**args)
