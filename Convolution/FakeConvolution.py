#!/usr/bin/env python3

'''
This is a drop-in replacement for the convolution executable.  It generates blank
Taylors and will overwrite multipoles if the parameter file requests it.
This is useful for singlestep performance testing when the convolution
isn't working on a parallel system, for example.

TODO: this only works on the parallel code right now

TODO: do we want to write to last.convlog?

'''

import sys
import argparse
import os
from os.path import join as pjoin

from mpi4py import MPI
import numpy as np

from Abacus import Tools
from Abacus.InputFile import InputFile

_print = print
print = lambda *args,**kwargs: print(*args,**kwargs,file=sys.stderr, flush=True)

def fake_convolution(parfile, verbose=False):
    par = InputFile(parfile)

    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    rank = comm.Get_rank()

    if verbose:
        print('FakeConvolution invoked on rank {} of {}'.format(rank, size))

    try:
        read = par['ReadStateDirectory']
    except KeyError:
        read = pjoin(par['WorkingDirectory'], 'read')

    multipoles = par['MultipoleDirectory']
    taylors = par['TaylorDirectory']
    cpd = par['CPD']

    assert 'MultipoleDirectory2' not in par, 'Not implemented'
    assert 'TaylorDirectory2' not in par, 'Not implemented'

    node_slabs = np.loadtxt(pjoin(read, 'nodeslabs'), dtype=int)

    firstslab = node_slabs[rank]  # inclusive
    lastslab = node_slabs[(rank+1)%size]  # exclusive
    nslab = (lastslab - firstslab)%cpd

    overwrite = par['Conv_IOMode'].lower() == 'overwrite'

    if verbose:
        print('Rank {} will process {} slabs [{},{}) with overwrite = {}'.format(rank, nslab, firstslab, lastslab, overwrite))

    rml = (par['Order']+1)**2
    ncomplex = cpd*(cpd+1)//2*rml
    cdt = np.complex64  # complex128 for double precision

    os.makedirs(taylors, exist_ok=True)

    zero_taylors = np.zeros(ncomplex, dtype=cdt)
    
    for s in range(firstslab,firstslab+nslab):
        s = s % cpd

        if overwrite:
            try:
                os.remove(pjoin(multipoles, 'Multipoles_{:04d}'.format(s)))
            except FileNotFoundError:
                pass  # Multipoles may not exist; that's fine

        zero_taylors.tofile(pjoin(taylors, 'Taylor_{:04d}'.format(s)))

    if verbose:
        print('Rank {} done.'.format(rank))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=Tools.ArgParseFormatter)
    parser.add_argument('parfile', help='Abacus parameter file')
    parser.add_argument('--verbose', help='Print status information', action='store_true')

    args = parser.parse_args()
    args = vars(args)

    fake_convolution(**args)
