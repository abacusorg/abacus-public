#!/usr/bin/env python

'''
This script converts a directory of Abacus slab files
(e.g. initial condition files) from one CPD value to another.

We could add Abacus functionality to do this on-the-fly
(at least for ICs), but I'm not sure such a scheme
would be multi-node friendly.  I think we'd prefer
to just convert the files ahead of time, which this
script allows you to do.

Results will be written in a new directory in the same
directory as the input directory.  The new CPD will be
appended to the output name.

Default command line arguments are chosen for the case
of changing the CPD of standard Abacus ICs.
'''

import argparse
from glob import glob
from os.path import dirname, basename, abspath, join as pjoin
import shutil
import os

from Abacus import Tools
from Abacus import ReadAbacus
from Abacus import WriteAbacus

def main(input_directory, in_pat, in_fmt, out_fmt, new_cpd, box, flip=None, verbose=False):
    # Create the output directory, removing it if it exists.
    input_directory = abspath(input_directory)
    indirname = basename(input_directory)
    outdirname = indirname + '_' + str(new_cpd)
    outdir = pjoin(dirname(input_directory), outdirname)
    shutil.rmtree(outdir, ignore_errors=True)
    os.makedirs(outdir, exist_ok=True)

    writer = WriteAbacus.SlabWriter(NP=0, cpd=new_cpd, boxsize=box, outdir=outdir, format=out_fmt, verbose=True)

    files = sorted(glob(pjoin(input_directory, in_pat)))
    for fn in files:
        print('Processing file "{}"'.format(fn))
        # for now, PID is mandatory
        particles = ReadAbacus.read(fn, format=in_fmt, return_pid=True)
        writer.NP += len(particles)  # writer will check that it wrote NP particles at the end

        writer.ingest(particles=particles)

        del particles


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=Tools.ArgParseFormatter)
    parser.add_argument('input_directory', help='The directory containing the slab files')
    parser.add_argument('new_cpd', help="The cells-per-dimension (number of slabs) to divide the output across", type=int)
    parser.add_argument('--in_fmt', help='Input file format. Can be any ReadAbacus file format.', default='RVZel')
    parser.add_argument('--out_fmt', help='Output file format', choices=['Zeldovich', 'RVdouble', 'RVdoubleZel', 'RVZel', 'RVTag', 'same'], default='RVZel')
    parser.add_argument('--in_pat', help="The file globbing pattern in the input DIRECTORY", type=str, default='ic_*')
    parser.add_argument('--box', help="The box size of the input (and output) particles.  Outputs will be zero-centered.", type=float, default=1.)
    parser.add_argument('--flip', help='Interpret positions as displacements and reverse them.', action='store_true')
    
    args = parser.parse_args()
    args = vars(args)
    exit(main(**args))
