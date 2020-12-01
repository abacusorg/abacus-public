#!/usr/bin/env python3
'''
For extremely large CheckForces problems, we sometimes want to run CheckForces,
store the results, and analyze them later/on a different system.  This script
facilitates analysis with MPI, mainly to accelerate reading the TB of acc
data from network file systems.
'''

import argparse
from glob import glob
from os.path import join as pjoin

from mpi4py import MPI
import numpy as np

from Abacus.Tools import ArgParseFormatter
from Abacus.InputFile import InputFile

import CheckForces

DEFAULT_DTYPE = 'f4'

def main(accdir, dtype=DEFAULT_DTYPE):
    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    rank = comm.Get_rank()
    root = 0

    accfns = sorted(glob(pjoin(accdir, 'acc_*')))
    cpd = len(accfns)

    # check that we have all acc slabs
    params = InputFile(pjoin(accdir, 'abacus.par'))
    assert cpd == params['CPD']

    start, end = rank*cpd//size, (rank+1)*cpd//size
    print(f'Rank {rank} processing slabs [{start},{end})')

    myaccfns = [pjoin(accdir, f'acc_{i:04d}') for i in range(start,end)]

    acc_mag, results = CheckForces.analyze_storeforces(params, dtype, slabfns=myaccfns, silent=True, raw=True)

    # Make a dict of nans, will be overwritten with reductions
    all_results = {k:np.full_like(results[k], np.nan) for k in results}
    mpi_ops = {'Max':MPI.MAX,
    		   'Sum':MPI.SUM,
    		   'Min':MPI.MIN,
    		   'Min (nonzero)':MPI.MIN,
    		   'Sum of squares':MPI.SUM}

    for k in results:
	    comm.Reduce(results[k], all_results[k], op=mpi_ops[k], root=root)
    del results

    # Convert 'Sum' into 'Mean' etc
    CheckForces.finalize_results(all_results, params.NP)

    if rank == root:
        print('Global force statistics (equivalent ZA displacement in units of interparticle spacing):')
        for k in all_results:
            print('{k}: {v}'.format(k=k, v=all_results[k]))

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description=__doc__, formatter_class=ArgParseFormatter)
    parser.add_argument('accdir', help='The directory containing the acc slabs')
    parser.add_argument('--dtype', default=DEFAULT_DTYPE, choices=['f4', 'f8'])
    args = parser.parse_args()
    args = vars(args)

    main(**args)
