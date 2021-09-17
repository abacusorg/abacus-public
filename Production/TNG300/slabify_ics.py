#!/usr/bin/env python
'''Convert Illustris ICs to Abacus.

Normally this would be a job for convert_particles_to_ic.py. But that plumbs the reading through ReadAbacus,
which would need all sorts of special-casing for the conventions of Illustris HDF5 IC files.  So here we are...
'''

import os
from pathlib import Path

import h5py
from Abacus import WriteAbacus

h5fn = Path(os.environ['HOME']) / 'ceph/illustris/ic/TNG300-2/snap_ics.hdf5'
cpd = 405
outdir = h5fn.parent
outformat = 'RVPID'

with h5py.File(h5fn, 'r') as h:
    pos = h['PartType1']['Coordinates'][:]
    vel = h['PartType1']['Velocities'][:]
    pids = h['PartType1']['ParticleIDs'][:]
    
    boxsize = h['Header'].attrs['BoxSize']/1e3
    zinit = h['Header'].attrs['Redshift']
    
    print('Read the following parameters from header:')
    print(dict(h['Header'].attrs))

# Convert from Gadget internal velocities to peculiar km/s
vel *= (1+zinit)**-0.5

writer = WriteAbacus.SlabWriter(len(pos), cpd, boxsize, outdir, outformat, verbose=True)
writer(pos=pos, vel=vel, pid=pids)
del writer
