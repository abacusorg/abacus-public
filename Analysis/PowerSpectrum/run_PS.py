#!/usr/bin/env python

'''
A user-friendly command-line interface to PowerSpectrum.

Author: Lehman Garrison
'''

import sys
import os
from os.path import join as pjoin
import PowerSpectrum as PS
from Abacus import InputFile
import pdb
import shutil
import argparse
import numpy as np
from glob import glob
from Abacus.Analysis import common
import pynbody
import TSC

# A label for the power spectrum based on the properties
def ps_suffix(**kwargs):
    nfft = kwargs['nfft']
    try:
        len(nfft)
    except TypeError:  # convert a plain int to a singlet array
        nfft = [nfft]
        
    # if 2D, add a projected tag
    projected = len(nfft) == 2
    
    # if all dimensions the same, collapse to single int
    if (np.diff(nfft) == 0).all():
        nfft = nfft[0]
    
    ps_fn = ''
    if not kwargs['just_density']:
        ps_fn += '_nfft{}'.format(kwargs['nfft'])
    else:
        ps_fn += '_ngrid{}'.format(kwargs['nfft'])
    if kwargs['dtype'] != np.float32:
        ps_fn += '_{}bit'.format(np.dtype(kwargs['dtype']).itemsize*8)
    if kwargs['zspace']:
        ps_fn += '_zspace'
    if kwargs['rotate_to']:
        ps_fn += '_rotate({2.1g},{2.1g},{2.1g})'.format(*rotate_to)
    if projected:
        ps_fn += '_projected'
        
    if not kwargs['just_density']:
        ps_fn += '.csv'
    return ps_fn


# If the input is a PS,
# we put the result in the same directory
def get_output_path_ps(input):
    return os.path.dirname(input)
    
    
def run_PS_on_dir(folder, **kwargs):
    patterns = [pjoin(folder, '*.dat'), pjoin(folder, 'ic_*'), pjoin(folder, 'position_*'), pjoin(folder, r'*.[0-9]*')]
    # Decide which pattern to use
    for pattern in patterns:
        if glob(pattern):
            break
    else:
        assert len(glob(pattern)) > 0, 'Could not find any matches to ' + str(patterns)

    # Make a descriptive filename
    this_suffix = ps_suffix(**kwargs)
    
    just_density = kwargs.pop('just_density')
    product_name = 'power' if not just_density else 'density'
    outdir = common.get_output_dir(product_name, folder)
    ps_fn = product_name + this_suffix

    # Read the header to get the boxsize
    header_pats = [pjoin(folder, 'header'), pjoin(folder, os.pardir, 'info/*.par'), pjoin(folder, os.pardir, '*.par')]
    headers = sum([glob(h) for h in header_pats], [])
    for header_fn in headers:
        try:
            header = InputFile.InputFile(header_fn)
            BoxSize = header['BoxSize']
            del header
        except IOError:
            continue
        break
    else:
        try:
            del header_fn
        except:
            pass  # Maybe no headers were found 
        # Maybe it's a gadget file?
        try:
            gadget_fn = sorted(glob(pattern))[0]
            f = pynbody.load(gadget_fn)
            BoxSize = float(f.properties['boxsize'])
        except:
            print 'Could not find a header in ' + str(header_pats)
            print 'or as a gadget file'
            raise
            
    # Make the output dir and store the header
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    try:
        shutil.copy(header_fn, outdir + '/header')
    except:
        # If a normal header doesn't exist, try copying an example gadget file
        shutil.copy(gadget_fn, outdir + '/gadget_header')

    # Copy the sim-level parameter files
    if not os.path.isdir(outdir + '/../info'):
        try:
            shutil.copytree(folder+'/../info', outdir + '/../info')
        except:
            print 'Could not copy ../info'
    
    print 'Starting {} on {}'.format('PS' if not just_density else 'density', pattern)
    save_fn = os.path.join(outdir, ps_fn)
    print 'and saving to {}'.format(save_fn)
    
    if not just_density:
        k,s,nmodes = PS.CalculateBySlab(pattern, BoxSize, kwargs.pop('nfft'), **kwargs)
        np.savetxt(save_fn, zip(k,s,nmodes), delimiter=',', header='k, P(k), N_modes')
    else:
        density = TSC.BinParticlesFromFile(pattern, BoxSize, kwargs.pop('nfft'), **kwargs)
        density.tofile(save_fn)

    # Touch ps_done
    with open(folder + '/ps_done', 'a'):
        os.utime(folder + '/ps_done', None)

            
def run_PS_on_PS(input_ps_fn, **kwargs):    
    outdir = get_output_path_ps(input_ps_fn)
    output_ps_fn = os.path.basename(input_ps_fn)
    
    # Make a descriptive filename
    output_ps_fn += ps_suffix(**kwargs)

    # Read the header
    try:
        header_fn = glob(os.path.dirname(input_ps_fn)+'/*.par')[0]
        header = InputFile.InputFile(header_fn)
    except IOError:
        print 'Could not find "*.par"'
        
    # Load the input PS
    input_ps = np.loadtxt(input_ps_fn)
    
    k,s,nmodes = PS.RebinTheoryPS(input_ps, header.BoxSize, kwargs['nfft'], nbins=kwargs['nbins'], log=kwargs['log'], dtype=kwargs['dtype'])
    np.savetxt(outdir + '/' + output_ps_fn, zip(k,s,nmodes), delimiter=',',header='k, P(k), N_modes')

    
# folders is a list of slice folders
def run_PS(inputs, **kwargs):
    for input in inputs:
        # If the input is an output or IC directory
        if os.path.isdir(input):
            if kwargs.get('format').lower() == 'pack14':
                with common.extract_slabs(input):
                    run_PS_on_dir(input, **kwargs)
            else:  # pretty kludgy...
                run_PS_on_dir(input, **kwargs)
        # If the input is a PS file
        elif os.path.isfile(input):
            run_PS_on_PS(input, **kwargs)
        else:
            raise ValueError(input, "does not exist!")
        print 'Finished PS.'
        
def vector_arg(s):
    try:
        x, y, z = map(int, s.strip('()[]').split(','))
        return x, y, z
    except:
        raise argparse.ArgumentTypeError("Vector must be x,y,z")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Compute power spectra on Abacus outputs or ICs.  Can also evaluate a power spectrum on an FFT mesh.')
    parser.add_argument('input', help='The timeslice outputs (or IC directories, or power spectrum file) on which to run PS', nargs='+')
    parser.add_argument('--nfft', help='The size of the FFT (side length of the FFT cube).  Default: 1024', default=1024, type=int)
    parser.add_argument('--format', help='Format of the data to be read.  Default: Pack14', default='Pack14', choices=['RVdouble', 'Pack14', 'RVZel', 'state', 'gadget'])
    parser.add_argument('--rotate-to', help='Rotate the z-axis to the given axis [e.g. (1,2,3)].  Rotations will shrink the FFT domain by sqrt(3) to avoid cutting off particles.', default=None, type=vector_arg, metavar='(X,Y,Z)')
    parser.add_argument('--projected', help='Project the simulation along the z-axis.  Projections are done after rotations.', action='store_true')
    parser.add_argument('--zspace', help='Displace the particles according to their redshift-space positions.', action='store_true')
    parser.add_argument('--nbins', help='Number of k bins.  Default: nfft/4.', default=-1, type=int)
    parser.add_argument('--dtype', help='Data type for the binning and FFT.', choices=['float', 'double'], default='float')
    parser.add_argument('--log', help='Do the k-binning in log space.', action='store_true')
    parser.add_argument('--just-density', help='Write the density cube out to disk instead of doing an FFT.', action='store_true')
    
    args = parser.parse_args()
    args = vars(args)
    
    if args['just_density']:
        # No meaning for density cube
        args.pop('nbins'); args.pop('log')
    
    if args.pop('projected'):
        args['nfft'] = [args['nfft'],]*2
    
    if args['dtype'] == 'float':
        args['dtype'] = np.float32
    elif args['dtype'] == 'double':
        args['dtype'] = np.float64

    run_PS(args.pop('input'), **args)
    