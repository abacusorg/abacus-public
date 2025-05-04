#!/usr/bin/env python3
"""
Convert a files in a non-Abacus format (like Gadget) to an Abacus format.
Output files will go in a new directory in the input file directory.
Slabs are written in chunks as particles are loaded, so this should
be memory-efficient for large sims.

Particles are arranged into slabs based on their x-positions and 
wrapped to [-BoxSize/2,BoxSize/2)

Caution: if chosing an "Zel"-type output format, the positions will be
later interpreted as displacements by abacus.

Units
=====

Gadget
------
We convert automatically from the Gadget velocity convention to km/s.
Set ICVelocity2Displacement=-1 in the .par file to automatically handle km/s.

Usage
=====
$ ./convert_to_abacus.py --help

"""

import sys
import struct
import os
import os.path as op
import shutil
from glob import glob
import argparse
from collections import defaultdict
import re
from warnings import warn

import numpy as np

from Abacus import Tools
from Abacus import WriteAbacus

supported_formats = ['gadget', 'desi-hdf5']


def main(files, cpd, input_format, output_format, flip=False, input_units=None):
    if flip and 'Zel' not in output_format:
        # TODO: we're not really being careful about displacement vs global
        warn('Specified "flip" but not using a displacement-oriented output format.')

    input_format = input_format.lower()
    assert input_format in supported_formats

    if input_format == 'gadget':
        import pynbody
    elif input_format == 'desi-hdf5':
        import h5py
    
    # If we're dealing with files with numerical suffixes, sort by that suffix
    files = sort_files(files)
        
    # First count the particles in all files
    NP = 0
    for fn in files:
        if input_format == 'gadget':
            f = pynbody.load(fn)
            NP += len(f)
        elif input_format == 'desi-hdf5':
            f = h5py.File(fn, 'r')
            NP += len(f['/Matter/Position'])

    ppd = int(round(NP**(1./3)))
    print(f'Found {NP} particles (ppd {ppd}) in {len(files)} files')
    if ppd**3 != NP:
        warn('not perfect cube!')

    if input_units == 'Quijote':
        print('Loading files assuming Quijote units on disk (pos: a kpc/h; vel: a^0.5 km/s)')

        # For pynbody
        units_system = dict(distance='h^-1 kpc a')
    else:
        units_system = dict(distance='h^-1 Mpc a')
    
    # Get some properties
    if input_format == 'gadget':
        f.set_units_system(**units_system)
        pynbody.snapshot.gadget._do_properties(f) # why do we have to trigger this?

        BoxSize = float(f.properties['boxsize'].in_units('h^-1 Mpc a'))
        print('ppd: {:d}, cpd: {:d}, boxsize: {:f}, redshift: {:f}'.format(ppd, cpd, BoxSize, 1./f.properties['a'] - 1))
        print('All properties:', str(f.properties))
    elif input_format == 'desi-hdf5':
        header = dict(f['/Header'].attrs)
        print('DESI HDF5 header:', str(header))
        if NP != header['NP.Matter']:
            warn(f'Total counted particles {NP} does not match number in header {header["NP.Matter"]}')
        BoxSize = header['BoxSize']
        
    # Format the output dir path
    # TODO: we're assuming all files are in the same dir
    out_path = op.dirname(op.abspath(files[0]))
    if flip:
        out_path += "_flip"
    out_path += f"_{output_format}_cpd{cpd}"
    
    # Make the output dir
    if op.isdir(out_path):
        warn("ICs already exist; removing.")
        shutil.rmtree(out_path)
    os.makedirs(out_path, exist_ok=False)
    print(f'Outputting to: {out_path}')
    
    # Now read the data
    writer = WriteAbacus.SlabWriter(NP, cpd, BoxSize, out_path, output_format)
    for fn in files:
        print(f'Loading particles from {fn}')

        if input_format == 'gadget':
            f = pynbody.load(fn)

            f.set_units_system(**units_system)
            pynbody.snapshot.gadget._do_properties(f) # why do we have to trigger this?

            pos = f['pos']
            vel = f['vel']
            pid = f['iord']

            # in-place conversions
            pos.convert_units('h^-1 Mpc a')
            vel.convert_units('km s^-1')  # Set ICVelocity2Displacement=-1 in the .par file to automatically handle km/s

        elif input_format == 'desi-hdf5':
            f = h5py.File(fn, 'r')
            pos = f['/Matter/Position'][:]
            vel = f['/Matter/Velocity'][:]
            pid = f['/Matter/ParticleID'][:]
        
        # Translate positions and wrap any stragglers
        Tools.wrap_zero_centered(pos, BoxSize)
        #assert (-BoxSize/2 <= pos).all() and (pos < BoxSize/2).all()  # numerical stability issues in the wrapping...

        if flip:
            pos *= -1
            vel *= -1
            
        writer(pos=pos, vel=vel, pid=pid)
        del f
    return 0


def sort_files(files):
    # First, trim a shared suffix
    suffix = longest_common_suffix(files)
    if suffix:
        print(f'Stripping suffix {suffix}')
        trimmed_fns = [fn[:-len(suffix)] for fn in files]
    else:
        trimmed_fns = files

    fnnum = re.compile(r'\d+$')
    fnnums = [fnnum.search(g) for g in trimmed_fns]
    if all(fnnums):
        files = zip(files, (int(n.group(0)) for n in fnnums))
        files = sorted(files, key=lambda f:f[1])
        files = [f[0] for f in files]
    else:
        files = sorted(files)

    return files


def longest_common_suffix(strings):
    s0 = strings[0]
    for i in range(-1,-len(s0)-1,-1):
        if not all(s[i:] == s0[i:] for s in strings):
            break
    else:
        return s0
    
    if i == -1:
        return ''

    return s0[i+1:]


    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=Tools.ArgParseFormatter)
    parser.add_argument('input_format', help='The input file format', choices=supported_formats)
    parser.add_argument('output_format', help='The output file format', choices=['Zeldovich', 'RVdouble', 'RVdoubleZel', 'RVZel', 'RVdoublePID', 'RVPID'])
    parser.add_argument('cpd', help="Cells-per-dimension (number of slabs) to divide the output across", type=int)
    parser.add_argument('files', help='One or more files to convert', nargs='+')
    parser.add_argument('--flip', help='Reverse the displacements relative to their initial Lagrangian grid point', action='store_true')
    parser.add_argument('--input-units', help='Convention for the input units', default=None, choices=['Quijote'])
    
    args = parser.parse_args()
    args = vars(args)

    exit(main(**args))
