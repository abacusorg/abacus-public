#!/usr/bin/env python

'''
This script facilitates getting states on and off nodes via MPI.  The most
common uses are when one wants to save and restore a state at the end
and beginning of a queue time allocation.

Another case is when we have a state from the serial code that we want to
distribute to the nodes to run in a parallel fashion.  We can set read state
to be the network disk, but then there's no way to tell the next step to use
the local storage. So it's probably just easier to copy everything to the right
node ahead of time. The '--distribute' mode of this state takes care of that.
It can be triggered by DistributeStateFrom = <source_dir> in the parameter file.

Usage
-----
$ ./move_node_states.py --help

'''


import sys
import argparse
import os
import shutil
from os.path import join as pjoin, basename
from glob import glob
import re

from mpi4py import MPI
import numpy as np

from Abacus import Tools
from Abacus import Reglob
from Abacus.InputFile import InputFile
      

def distribute_from_serial(parfile, source_dir, verbose=True):
    par = InputFile(parfile)

    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    rank = comm.Get_rank()

    if verbose:
        print('Distribute state invoked on rank {} of {}'.format(rank, size), file=sys.stderr)

    cpd = par['CPD']

    # TODO: for now, we'll assume nodeslabs exists.  Could make one here, though
    node_slabs = np.loadtxt(pjoin(source_dir, 'read', 'nodeslabs'), dtype=int)
    assert len(node_slabs) == size

    firstslab = node_slabs[rank]  # inclusive
    lastslab = node_slabs[(rank+1)%size]  # exclusive
    nslab = (lastslab - firstslab)%cpd

    # TODO: for now, we'll just copy over the read state.  Could do multipoles, etc, too
    try:
        localread = par['LocalReadStateDirectory']
    except KeyError:
        localread = pjoin(par['LocalWorkingDirectory'], 'read')

    try:
        read = par['ReadStateDirectory']
    except KeyError:
        read = pjoin(par['WorkingDirectory'], 'read')

    # Must delete existing dirs; not safe to potentially mix state from different states!
    try:
        shutil.rmtree(localread)
    except FileNotFoundError:
        pass
        
    try:
        shutil.rmtree(read)
    except FileNotFoundError:
        pass

    os.makedirs(localread)
    os.makedirs(read)

    for s in range(firstslab,firstslab+nslab):
        s = s % cpd

        fns = glob(pjoin(source_dir, 'read', '*_{:04d}'.format(s)))

        if verbose:
            print('Node {} copying {} files for slab {}'.format(rank, len(fns), s), file=sys.stderr)

        for fn in fns:
            shutil.copy(fn, localread)

    # The first rank will also copy the state file, etc, 
    if rank == 0:
        fns = glob(pjoin(source_dir, 'read', '*'.format(s)))
        names = [basename(f) for f in fns]

        # Some files we know we'll need for consistency
        assert 'state' in names
        assert 'slabsize' in names
        assert 'redlack' in names
        assert 'globaldipole' in names
        assert 'nodeslabs' in names

        for fn in fns:
            # If the filename does not look like 'asdf_1234', copy it
            if re.match(r'^((?!_\d{4}).)*$', basename(fn)):
                shutil.copy(fn, read)

def distribute_to_resume(parfile, verbose=True):
    par = InputFile(parfile)

    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    rank = comm.Get_rank()
    
    if verbose:
        print('Distribute_to_resume state invoked on rank {} of {}'.format(rank, size), file=sys.stderr)
    
    source = pjoin(os.path.dirname(par['WorkingDirectory']), par['SimName'] + '_retrieved_state', 'rank_' + str(rank))
    dest = pjoin(par['LocalWorkingDirectory'])
    
    
    try:
        shutil.rmtree(dest)
    except FileNotFoundError:
        pass
  
    if verbose:
        print('Will copy node {}s state from {} to {}'.format(rank, source, dest), file=sys.stderr)
    
    shutil.copytree(source, dest)
    
    # The first rank will also copy the state file, etc, 
    read_fns_present = True #assume they're present. node 0 will complain if they're not. 
    if rank == 0:
        localread = pjoin(os.path.dirname(par['WorkingDirectory']), par['SimName'], 'read') 
    
        try:
            shutil.rmtree(localread)
        except FileNotFoundError:
            pass
            
        os.mkdir(localread)    
        
        fns = glob(pjoin(os.path.dirname(par['WorkingDirectory']), par['SimName'] + '_retrieved_state', 'read', '*'))
        names = [basename(f) for f in fns]
        
        # Some files we know we'll need for consistency
        read_fns_present = ( set(['state', 'slabsize', 'redlack', 'globaldipole', 'nodeslabs']) <= set(names) )
        
        for fn in fns:
            # If the filename does not look like 'asdf_1234', copy it
            if re.match(r'^((?!_\d{4}).)*$', basename(fn)):
                print('Copying read state file {} to {}'.format(fn, localread), file=sys.stderr)
                shutil.copy(fn, localread)
    
    #wait for read state/multipole recovery if necessary. 
    comm.Barrier()
    read_fns_present = MPI.LAND(read_fns_present)  
    
    if verbose:
        print('Success distributing to resume! Read files present = {}'.format(read_fns_present)))
    
    return read_fns_present
          

def retrieve_state(parfile, verbose=True):
    par = InputFile(parfile)
    
    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    rank = comm.Get_rank()


    if verbose:
        print('Retrieve state invoked on rank {} of {}'.format(rank, size), file=sys.stderr)
    

    dest = pjoin(os.path.dirname(par['WorkingDirectory']), par['SimName'] + '_retrieved_state', 'rank_' + str(rank))
    
    backup_dest_root =  pjoin(os.path.dirname(par['WorkingDirectory']), par['SimName'] + '_retrieved_state_backup')
    backup_dest = pjoin(backup_dest_root, 'rank_' + str(rank))
    
    source = par['LocalWorkingDirectory']
    
    try:
        shutil.rmtree(backup_dest)
    except FileNotFoundError:
        pass
    
     
    try: 
        print('Renaming previous runs retrieved states to backup files')
        shutil.move(dest, backup_dest)
    except FileNotFoundError:
        pass


    if verbose:
        print('Will copy node {}s state from {} to {}'.format(rank, source, dest), file=sys.stderr)
    
    shutil.copytree(source, dest)
    
    #rank 0 will also back up the read directory misc. files, which are about to modified when we requeue during every consecutive timestep. 
    #We want to keep a copy of the set that corresponds to the backed up state we're retrieving. 
    if rank == 0:
        
        dest = pjoin(os.path.dirname(par['WorkingDirectory']), par['SimName'] + '_retrieved_state', 'read')
        
        try:
            shutil.rmtree(dest)
        except FileNotFoundError:
            pass
        
        os.mkdir(dest)
        
        fns = glob(pjoin(os.path.dirname(par['WorkingDirectory']), par['SimName'], 'read', '*'))
        names = [basename(f) for f in fns]

        assert 'state' in names
        assert 'slabsize' in names
        assert 'redlack' in names
        assert 'globaldipole' in names
        assert 'nodeslabs' in names

        for fn in fns:
            # If the filename does not look like 'asdf_1234', copy it
            if re.match(r'^((?!_\d{4}).)*$', basename(fn)): 
                shutil.copy(fn, dest)
    
    if verbose:
        print('Success retrieving state! Removing backup files.')

    try:
        shutil.rmtree(backup_dest_root)
    except FileNotFoundError:
        pass
        
    if verbose:
        print('State retrieval complete.')
    


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=Tools.ArgParseFormatter)
    parser.add_argument('parfile', help='Abacus parameter file')
    parser.add_argument('--source-dir', help='The state directory from which to copy')
    parser.add_argument('--verbose', help='Print status information', action='store_true')
    
    parser.add_argument('--distribute-from-serial', help='Distribute a serial state to nodes', action='store_true')
    parser.add_argument('--retrieve', help="Save all nodes' states to the global disk", action='store_true')
    parser.add_argument('--distribute', help="Distribute nodes' states from the global disk to the nodes", action='store_true')

    args = parser.parse_args()
    args = vars(args)

    if args['distribute_from_serial']:
        distribute_from_serial(args['parfile'], args['sourcedir'], args['verbose'])
    elif args['retrieve']:
        retrieve_state(args['parfile'], args['verbose'])
    elif args['distribute']:
        distribute_to_resume(args['parfile'], args['verbose'])

