#!/usr/bin/env python3

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
from os.path import join as pjoin, basename, dirname, normpath
from glob import glob
import re
from pathlib import Path

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


def distribute_to_resume(parfile, resumedir, verbose=True, allow_new_nnode=False):
    par = InputFile(parfile)

    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    rank = comm.Get_rank()
    
    if verbose:
        print(f'Distribute_to_resume state invoked on rank {rank} of {size}', file=sys.stderr)
    
    source = pjoin(resumedir, 'rank_' + str(rank))
    dest = pjoin(par['LocalWorkingDirectory'] + f'.{rank:04d}')

    try:
        shutil.rmtree(dest)
    except FileNotFoundError:
        pass

    try:
        state = InputFile(pjoin(resumedir, 'read', 'state'))
        old_nnodes = state['NodeSize']
        # For now, we will always trigger recover_multipoles if allow_new_nnode,
        # just because we don't have a good way to coordinate with the upstream Python
        #if old_nnodes == size:
        #    allow_new_nnode = False
    except:
        pass

    if allow_new_nnode:
        # Each node is probably getting slabs from multiple dirs
        state = InputFile(pjoin(resumedir, 'read', 'state'))
        old_nnodes = state['NodeSize']
        cpd = par['CPD']
        old_nodeslabs = np.loadtxt(pjoin(resumedir, 'read', 'nodeslabs'))
        slab_to_oldnode = lambda s: (np.argmin( (s-old_nodeslabs[0])%cpd >= (old_nodeslabs-old_nodeslabs[0])%cpd)-1)%old_nnodes
        new_nodeslabs = np.linspace(0,cpd,num=size+1, endpoint=True, dtype=int)  # we'll write to disk below

        if verbose:
            print(f'Redistributing state to {size} nodes; stored as {len(old_nodeslabs)} on disk.')

        os.makedirs(pjoin(dest,'read'))
        for i in range(*new_nodeslabs[[rank,rank+1]]):
            for fn in glob(pjoin(resumedir, f'rank_{slab_to_oldnode(i)}', 'read', f'*_{i:04d}')):
                #print(f'Copying {fn} to {dest}')
                shutil.copy(fn, pjoin(dest,'read'))

        new_nodeslabs = new_nodeslabs[:-1]  # the cpd entry is implicit on disk
  
    else:
        # Standard operation
        if verbose:
            print('Will copy node {}s state from {} to {}'.format(rank, source, dest), file=sys.stderr)
        
        # write state can be a symlink, must use symlinks=True
        shutil.copytree(source, dest, symlinks=True, copy_function=shutil.copy)
    
    # The first rank will also copy the state file, etc, 
    if rank == 0:
        try:
            read = par['ReadStateDirectory']
        except KeyError:
            read = pjoin(par['WorkingDirectory'], 'read')
    
        try:
            shutil.rmtree(read)
        except FileNotFoundError:
            pass
            
        os.mkdir(read)    
        
        fns = glob(pjoin(resumedir, 'read', '*'))
        names = [basename(f) for f in fns]
        
        # Some files we know we'll need for consistency
        #read_fns_present = int(( set(['state', 'slabsize', 'redlack', 'globaldipole', 'nodeslabs']) <= set(names) ))
        
        # Some files we know we'll need for consistency
        assert 'state' in names
        assert 'slabsize' in names
        assert 'redlack' in names
        assert 'globaldipole' in names

        if allow_new_nnode:
            names = [n for n in names if n != 'nodeslabs']  # don't copy old nodeslabs
            np.savetxt(pjoin(read,'nodeslabs'), new_nodeslabs, fmt='%d')  # write new nodeslabs
        else:
            assert 'nodeslabs' in names
        
        for fn in fns:
            if allow_new_nnode and basename(fn) == 'nodeslabs':
                continue
            # If the filename does not look like 'asdf_1234', copy it
            if re.match(r'^((?!_\d{4}).)*$', basename(fn)):
                if verbose:
                    print('Copying read state file {} to {}'.format(fn, read), file=sys.stderr)
                shutil.copy(fn, read)
    
    if verbose:
        print('Success distributing to resume!')
    
    #wait for read state/multipole recovery if necessary. 
    #comm.Barrier()
    #read_fns_present = comm.Reduce([read_fns_present, MPI.INT], [read_fns_present, MPI.INT], op=MPI.LAND, root=0) 
    #sys.exit(read_fns_present)
    
          

def retrieve_state(parfile, resumedir, verbose=True, delete_ics=False):
    
    par = InputFile(parfile)
    
    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    rank = comm.Get_rank()

    if rank == 0:
        print(f'Retrieve state invoked on {size} ranks', file=sys.stderr)
    
    resumedir = normpath(resumedir)
    dest = pjoin(resumedir, 'rank_' + str(rank))
    
    past =  pjoin(dirname(par['WorkingDirectory']), par['SimName'] + '_retrieved_state_past')
    past_deleting = past + '_DELETING'
    past_node = pjoin(past, 'rank_' + str(rank))
    
    source = par['LocalWorkingDirectory'] + f'.{rank:04d}'
    
    comm.Barrier() 
    
    if rank == 0:

        for p in (past,past_deleting):
            try:
                shutil.rmtree(p)
                print(f'Removed past backup {basename(p)}', file=sys.stderr)
            except FileNotFoundError:
                pass
            
    comm.Barrier() 
    
    if rank == 0:
        try: 
            os.rename(resumedir, past)
            print("Renamed previous run's retrieved states to backup files", file=sys.stderr)
        except FileNotFoundError:
            pass
    
    comm.Barrier()

    if verbose:
        print('Will copy node {}s state from {} to {}'.format(rank, source, dest), file=sys.stderr)
    
    # TODO: make read-only after copy, then read-write after distribution to nodes
    # write state can be a symlink, must use symlinks=True
    shutil.copytree(source, dest, symlinks=True, copy_function=shutil.copy)
    
    #rank 0 will also back up the read directory misc. files, which are about to modified when we requeue during every consecutive timestep. 
    #We want to keep a copy of the set that corresponds to the backed up state we're retrieving. 
    if rank == 0:
        
        dest = pjoin(dirname(par['WorkingDirectory']), par['SimName'] + '_retrieved_state', 'read')
        
        try:
            shutil.rmtree(dest)
        except FileNotFoundError:
            pass
        
        os.mkdir(dest)
        
        fns = glob(pjoin(dirname(par['WorkingDirectory']), par['SimName'], 'read', '*'))
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
    
    comm.Barrier() 

    if rank == 0:
        print('Backup complete. Removing any previous backup.', file=sys.stderr)

        try:
            os.rename(past, past_deleting)
            shutil.rmtree(past_deleting)
        except FileNotFoundError:
            pass

        if delete_ics:
            try: 
                shutil.rmtree(par['InitialConditionsDirectory'])
            except FileNotFoundError:
                pass    

    comm.Barrier() 
    
    if rank == 0:
        print('Success retrieving state!', file=sys.stderr)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=Tools.ArgParseFormatter)
    parser.add_argument('parfile', help='Abacus parameter file')
    parser.add_argument('--source-dir', help='The state directory from which to copy')
    parser.add_argument('--verbose', help='Print status information', action='store_true')
    
    # TODO: git-style subcommands
    parser.add_argument('--distribute-from-serial', help='Distribute a serial state to nodes', action='store_true')
    parser.add_argument('--retrieve', help="Save all nodes' states to the global disk", action='store_true')
    parser.add_argument('--distribute', help="Distribute nodes' states from the global disk to the nodes", action='store_true')
    parser.add_argument('--allow-new-nnode', help="Allow distribution to a different number of nodes than was saved on disk", action='store_true')
    parser.add_argument('--delete-ics', help="Allow deletion of ICs once the first backup has completed", action='store_true')
    parser.add_argument('--safe-cp', help="Use a slower but safer file copy (no os.sendfile())", action='store_true')
    
    parser.add_argument('resumedir', help="Global disk directory to redistribute from")

    args = parser.parse_args()
    args = vars(args)

    if args['safe_cp']:
        shutil._USE_CP_SENDFILE = False

    if args['distribute_from_serial']:
        distribute_from_serial(args['parfile'], args['sourcedir'], args['verbose'])
    elif args['retrieve']:
        retrieve_state(args['parfile'], args['resumedir'], args['verbose'], args['delete_ics'])
    elif args['distribute']:
        distribute_to_resume(args['parfile'], args['resumedir'], args['verbose'], args['allow_new_nnode'])
