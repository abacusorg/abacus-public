#!/usr/bin/env python
'''
Wrapper for the zeldovich binary.

Also provides an convenience function to compute sigma8
for a given redshift and cosmology.

usage: ./zeldovich.py <parameter file name>

'''

import os
import os.path as path
import subprocess
import sys
import shutil
import parser
import argparse
from .InputFile import InputFile
from . import GenParam
from . import abacus
from Abacus.Cosmology import AbacusCosmo

zeldovich_dir = path.join(abacus.abacuspath, 'zeldovich-PLT')
eigmodes_path = path.join(zeldovich_dir, 'eigmodes128')

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


def run(paramfn):
    params = InputFile(paramfn)

    if path.exists(params.InitialConditionsDirectory):
        print("Warning: old ICs already exist; removing.")
        shutil.rmtree(params.InitialConditionsDirectory)

    if not path.isdir(params.InitialConditionsDirectory):
        os.makedirs(params.InitialConditionsDirectory)

    if 'ZD_Pk_filename' in params:
        shutil.copy(params.ZD_Pk_filename, params.InitialConditionsDirectory + "/input.pow")
    else:
        assert 'ZD_Pk_powerlaw_index' in params
    subprocess.check_call([path.join(zeldovich_dir, "zeldovich"), paramfn])

    
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
    
    sim_dir = path.join(out_parent, old_params.SimName)
    new_parfn = path.join(sim_dir, 'info', new_parfn)
    ic_dir = path.join(sim_dir, 'ic')
    try:
        os.makedirs(sim_dir)
    except OSError:  # exists
        pass
        
    # If the eigmodes file doesn't exist, look for it in the zeldovich dir
    eigmodes_fn = old_params.ZD_PLT_filename
    if not path.isfile(eigmodes_fn):
        eigmodes_fn = path.join(zeldovich_dir, path.basename(eigmodes_fn))
        if not path.isfile(eigmodes_fn):
            print('Warning: original eigenmodes filename "{}" not found.  Falling back to most current eigenmodes.'.format(eigmodes_fn))
            eigmodes_fn = eigmodes_path
    
    # Copy over info dir and other important files
    # is parfn in the info dir?
    try: shutil.copytree(path.join(path.dirname(parfn), path.pardir, 'info'), path.join(sim_dir, 'info'))
    except: pass
    # is parfn next to the info dir?
    try: shutil.copytree(path.join(path.dirname(parfn), 'info'), path.join(sim_dir, 'info'))
    except: pass
    
    if not path.isdir(path.join(sim_dir, 'info')):
        print("Warning: no info dir found to copy")
        
    # If the Pk file doesn't exist, look for it in the info dir
    kwargs = {}
    if 'ZD_Pk_filename' in old_params:
        pk_fn = old_params.ZD_Pk_filename
        if not path.isfile(pk_fn):
            pk_fn = path.join(sim_dir, 'info', path.basename(pk_fn))
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
