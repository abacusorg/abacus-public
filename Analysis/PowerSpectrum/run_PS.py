#!/usr/bin/env python

'''
A user-friendly command-line interface to PowerSpectrum.

Author: Lehman Garrison
'''

import sys
import os
import PowerSpectrum as PS
from Abacus import InputFile
import pdb
import shutil
import argparse
import numpy as np
from glob import glob

# folders is a list of slice folders
def run_PS(folders, **kwargs):
    for folder in folders:
        split = folder.split('/')
        if len(split[-1]) == 0:
            split = split[:-1]
            
        split[-1] = split[-1].replace('slice', 'z')
        pattern = '{}/{}.{}.slab*.dat'.format(folder, split[-2], split[-1])
        header = InputFile.InputFile(folder+'/header')
        
        # Make a descriptive filename
        ps_fn = 'power'
        ps_fn += '_nfft{}'.format(kwargs['nfft'])
        if kwargs['dtype'] != np.float32:
            ps_fn += '_{}bit'.format(kwargs['dtype'].itemsize)
        if kwargs['zspace']:
            ps_fn += '_zspace'
        if kwargs['rotate_to']:
            ps_fn += '_rotate({2.1g},{2.1g},{2.1g})'.format(*rotate_to)
        if kwargs['projected']:
            ps_fn += '_projected'
        ps_fn += '.csv'
        
        # Change emulator_00 to emulator_00_power
        split[-3] += '_products'
        split.insert(-2,split[-2] + '_products')
        split[-2] += '_power'
        outdir = '/'.join(split) + '/'
        
        # Make the output dir and store the header
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        shutil.copy(folder+'/header', outdir)
        
        # Copy the sim-level parameter files
        for pf in glob(folder+'/../*'):
            if os.path.isfile(pf):
                shutil.copy(pf, outdir+'/..')
        
        print 'Starting PS on {}'.format(pattern)
        print 'and saving to {}/{}'.format(outdir, ps_fn)
        k,s,nmodes = PS.CalculateBySlab(pattern, header.BoxSize, kwargs.pop('nfft'), **kwargs)
        np.savetxt(outdir + ps_fn, zip(k,s), delimiter=',')
        print 'Finished.'
        
        # Touch .ps_done
        with open(folder + '/ps_done', 'a'):
            os.utime(folder + '/ps_done', None)

        
def vector_arg(s):
    try:
        x, y, z = map(int, s.strip('()[]').split(','))
        return x, y, z
    except:
        raise argparse.ArgumentTypeError("Vector must be x,y,z")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Compute power spectra of Abacus outputs.')
    parser.add_argument('slice_folders', help='The timeslice output(s) on which to run PS', nargs='+')
    parser.add_argument('--nfft', help='The size of the FFT (side length of the FFT cube).  Default: 1024', default=1024, type=int)
    parser.add_argument('--format', help='Format of the Abacus timeslice outputs.  Default: Pack14', default='Pack14', choices=['RVdouble', 'LC', 'Pack14'])
    parser.add_argument('--rotate-to', help='Rotate the z-axis to the given axis [e.g. (1,2,3)].  Rotations will shrink the FFT domain by sqrt(3) to avoid cutting off particles.', default=False, type=vector_arg, metavar='(X,Y,Z)')
    parser.add_argument('--projected', help='Project the simulation along the z-axis.  Projections are done after rotations.', action='store_true')
    parser.add_argument('--zspace', help='Displace the particles according to their redshift-space positions.', action='store_true')
    parser.add_argument('--nbins', help='Number of k bins.  Default: nfft/4.', default=-1, type=bool)
    parser.add_argument('--dtype', help='Data type for the binning and FFT.', choices=['float', 'double'], default='float')
    parser.add_argument('--log', help='Do the k-binning in log space.', action='store_true')
    
    args = parser.parse_args()
    args = vars(args)
    
    if args['dtype'] == 'float':
        args['dtype'] = np.float32
    elif args['dtype'] == 'double':
        args['dtype'] = np.float64

    run_PS(args.pop('slice_folders'), **args)
    