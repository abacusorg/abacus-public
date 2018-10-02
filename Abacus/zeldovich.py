#!/usr/bin/env python
'''
Wrapper for the zeldovich binary.

Also provides an convenience function to compute sigma8
for a given redshift and cosmology.

usage: ./zeldovich.py <parameter file name>

'''

import os
import os.path as path
from os.path import join as pjoin
import subprocess
import sys
import shutil
import parser
import argparse
from .InputFile import InputFile
from . import GenParam
from . import abacus
from Abacus.Cosmology import AbacusCosmo

zeldovich_dir = pjoin(abacus.abacuspath, 'zeldovich-PLT')
eigmodes_path = pjoin(zeldovich_dir, 'eigmodes128')

# Calculate sigma8 by scaling back params['sigma_8'] from z=0 to the given redshift
def calc_sigma8(params, z='init'):
    if z == 'init':
        z = params['InitialRedshift']
    my_cosmology = AbacusCosmo.MyCosmology()
    my_cosmology.Omega_m = params['Omega_M']
    my_cosmology.Omega_K = params['Omega_K']
    my_cosmology.Omega_DE = params['Omega_DE']
    my_cosmology.H0 = params['H0']
    my_cosmology.w0 = params['w0']
    my_cosmology.wa = params['wa']

    cosmology = AbacusCosmo.Cosmology(1./(1+z), my_cosmology)
    sig8_today = params['sigma_8']
    sig8_at_zinit = sig8_today * cosmology.current.growth / cosmology.today.growth
    
    return sig8_at_zinit


def run(paramfn, allow_eigmodes_fn_override=False):
    '''
    Invokes the zeldovich executable with the given parameter file,
    cleaning up any exisitng output directories first and also
    copying the input power spectrum to the destination.

    If `allow_eigmodes_fn_override` is set, checks if
    the `ZD_PLT_filename` parameter is valid and sets
    it to the current eigmodes file if not.
    '''
    params = InputFile(paramfn)

    if path.exists(params.InitialConditionsDirectory):
        print("Warning: old ICs already exist; removing.")
        shutil.rmtree(params.InitialConditionsDirectory)

    if not path.isdir(params.InitialConditionsDirectory):
        os.makedirs(params.InitialConditionsDirectory)

    if 'ZD_PLT_filename' in params and not path.isfile(params.ZD_PLT_filename):
        # If the filename is the same, then the files should be identical and we can silently override
        if path.basename(params.ZD_PLT_filename) != path.basename(eigmodes_path):
            print('Eigenmodes file "{}" did not exist!'.format(params.ZD_PLT_filename))
            if allow_eigmodes_fn_override:
                print('Overriding to most recent file: "{}".'.format(eigmodes_path))
            else:
                raise ValueError(allow_eigmodes_fn_override)
        params = GenParam.makeInput(paramfn, paramfn, ZD_PLT_filename=eigmodes_path)

    if 'ZD_Pk_filename' in params:
        shutil.copy(params['ZD_Pk_filename'], pjoin(params['InitialConditionsDirectory'], "input.pow"))
    else:
        assert 'ZD_Pk_powerlaw_index' in params
    subprocess.check_call([pjoin(zeldovich_dir, "zeldovich"), paramfn])

    
def run_override_dirs(parfn, out_parent, new_parfn='abacus_ic_fixdir.par'):
    """
    Sometimes we want to regenerate the ICs from a sim run on
    another machine.  Thus, the directories are probably wrong.
    If we have the .par2 file, we can replace the environment
    variables, but if we don't, we need to override the directories,
    which is what this routine is used for.
    
    Note: we generally assume the existence of an info dir,
    which may not be a good assumption.
    """
    
    old_params = InputFile(parfn)
    out_parent = path.abspath(out_parent)
    
    sim_dir = pjoin(out_parent, old_params.SimName)
    new_parfn = pjoin(sim_dir, 'info', new_parfn)
    ic_dir = pjoin(sim_dir, 'ic')
    try:
        os.makedirs(sim_dir)
    except OSError:  # exists
        pass
        
    # If the eigmodes file doesn't exist, look for it in the zeldovich dir
    eigmodes_fn = old_params.ZD_PLT_filename
    if not path.isfile(eigmodes_fn):
        eigmodes_fn = pjoin(zeldovich_dir, path.basename(eigmodes_fn))
        if not path.isfile(eigmodes_fn):
            print('Warning: original eigenmodes filename "{}" not found.  Falling back to most current eigenmodes.'.format(eigmodes_fn))
            eigmodes_fn = eigmodes_path
    
    # Copy over info dir and other important files
    # is parfn in the info dir?
    try: shutil.copytree(pjoin(path.dirname(parfn), path.pardir, 'info'), pjoin(sim_dir, 'info'))
    except: pass
    # is parfn next to the info dir?
    try: shutil.copytree(pjoin(path.dirname(parfn), 'info'), pjoin(sim_dir, 'info'))
    except: pass
    
    if not path.isdir(pjoin(sim_dir, 'info')):
        print("Warning: no info dir found to copy")
        
    # If the Pk file doesn't exist, look for it in the info dir
    kwargs = {}
    if 'ZD_Pk_filename' in old_params:
        pk_fn = old_params.ZD_Pk_filename
        if not path.isfile(pk_fn):
            pk_fn = pjoin(sim_dir, 'info', path.basename(pk_fn))
            assert path.isfile(pk_fn)
            kwargs['ZD_Pk_filename'] = pk_fn
    else:
        assert 'ZD_Pk_powerlaw_index' in old_params
        
    GenParam.makeInput(new_parfn, parfn, InitialConditionsDirectory=ic_dir, ZD_PLT_filename=eigmodes_fn, **kwargs)
    
    run(new_parfn)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Run the zeldovich code (located in {})'.format(zeldovich_dir))
    parser.add_argument('parfile', help='The parameter file.  This is usually the same as the .par file for singlestep.', nargs='+')
    parser.add_argument('--out-parent', help="Overrides the parfile InitialConditionsDirectory (i.e. the zeldovich output directory) with PARENT/SimName/ic.  Create a new abacus_ic.par with the modified parameters: IC dir; and if don't exist: eigmodes, camb_matterpower", metavar='PARENT')
    
    args = parser.parse_args()
    args = vars(args)
    out_parent = args.pop('out_parent')
    
    for parfn in args.pop('parfile'):
        if out_parent:
            run_override_dirs(parfn, out_parent, **args)
        else:
            run(parfn, **args)
