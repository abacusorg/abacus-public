#!/usr/bin/env python3
'''
Wrapper for the zeldovich binary.

Also provides an convenience function to compute sigma8
for a given redshift and cosmology.

Valid ways to specify P(k) normalization
========================================
The default zeldovich code behavior is to discard the normalization
of the input power spectrum file and use the value of `ZD_Pk_sigma`
in the parameter file.  From Abacus, one may trigger this behavior
by setting `sigma_8` in the parameter file instead of `ZD_Pk_sigma`
directly.  The `sigma_8` parameter is interpreted as the z=0 value,
and the Python parsing step of the `abacus.par2` file will see this
and call the Abacus cosmology module to compute `ZD_Pk_sigma` via
the growth ratio:

    ZD_Pk_sigma = sigma_8*D(z=InitialRedshift)/D(z=0)

That's all fine when one has a z=0 power spectrum file; then one
doesn't even need to know the file's units of power.  But when your
power spectrum file is at another redshift like z=1, what `sigma_8`
do you specify?  It's not the z=0 value; the shape of the power
spectrum is allowed to change, even aside from the amplitude.

So one would like to know the sigma_8 of the z=1 file, but CAMB
or CLASS may not output that.  So we could extract our sigma_8
integrator from the zeldovich code... or we could just trust the
input file normalization and pass a growth *ratio* instead of
target amplitude.  This parameter is called `ZD_Pk_sigma_ratio`:

    ZD_Pk_sigma_ratio = D(z=InitialRedshift)/D(z=ZD_Pk_file_redshift)

The other new parameter is `ZD_Pk_file_redshift`, which tells us
the redshift of the input power spectrum file.

From the Abacus Python perspective, one may either specify
`ZD_Pk_file_redshift` or `sigma_8`, not both.

From the zeldovich perspective, one may specify `ZD_Pk_sigma` or
`ZD_Pk_sigma_ratio`, not both.

Usage
=====
$ ./zeldovich.py --help

'''

import os
import os.path as path
from os.path import join as pjoin
import shutil
import argparse
import shlex

import numpy as np

from Abacus import GenParam
from Abacus import abacus
from Abacus.Tools import chdir
import Abacus.Cosmology
from Abacus import transpose_ic

zeldovich_dir = pjoin(abacus.abacuspath, 'external', 'zeldovich-PLT')
eigmodes_path = pjoin(zeldovich_dir, 'eigmodes128')

on_the_fly_formats = ['poisson', 'lattice']


def is_on_the_fly_format(fmt):
    return fmt.lower() in on_the_fly_formats


# Calculate sigma8 by scaling back params['sigma_8'] from z=0 to the given redshift
def calc_sigma8(params, z='init'):
    sig8_today = params['sigma_8']
    sig8_at_zinit = sig8_today * calc_growth_ratio(params, z, 0.)
    
    return sig8_at_zinit


# Calculate the growth ratio D(z1)/D(z2)
def calc_growth_ratio(params, z1, z2):
    if z1 == 'init':
        z1 = params['InitialRedshift']
    if z2 == 'init':
        z2 = params['InitialRedshift']

    cosm1 = Abacus.Cosmology.from_params(params, z1)
    cosm2 = Abacus.Cosmology.from_params(params, z2)

    return cosm1.current.growth/cosm2.current.growth


def setup_zeldovich_params(params):
    '''
    Given a parameters dict with a cosmology, set up the power spectrum
    amplitude and growth rate.
    '''
    zd_params = {}

    # If given sigma_8, set ZD_Pk_sigma
    if params.get('sigma_8'):
        if params.get('ZD_Pk_file_redshift'):
            raise ValueError("Must specify one of sigma_8 and ZD_Pk_file_redshift in parameter file")

        sigma8_at_zinit = calc_sigma8(params)

        # If ZD_Pk_sigma was already given, check that it is consistent with what we just computed
        if 'ZD_Pk_sigma' in params:
            if not np.isclose(params['ZD_Pk_sigma'], sigma8_at_zinit, rtol=1e-5):
                raise ValueError(f"ZD_Pk_sigma ({sigma8_at_zinit:f}) calculated from sigma_8 ({params['sigma_8']:f}) conflicts with the provided value of ZD_Pk_sigma ({params['ZD_Pk_sigma']:f})!")
        zd_params['ZD_Pk_sigma'] = sigma8_at_zinit

    elif params.get('ZD_Pk_file_redshift') is not None:
        # If given ZD_Pk_file_redshift, set ZD_Pk_sigma_ratio
        if params.get('sigma_8'):
            raise ValueError("Must specify one of sigma_8 and ZD_Pk_file_redshift in parameter file")

        growth_ratio = calc_growth_ratio(params, 'init', params['ZD_Pk_file_redshift'])

        #If ZD_Pk_sigma_ratio was already given, check that it is consistent with what we just computed
        if 'ZD_Pk_sigma_ratio' in params:
            if not np.isclose(growth_ratio, params['ZD_Pk_sigma_ratio'], rtol=1e-5):
                raise ValueError(f'ZD_Pk_sigma_ratio = {params["ZD_Pk_sigma_ratio"]} in parameter file does not match value {growth_ratio} computed from ZD_Pk_file_redshift = {params["ZD_Pk_file_redshift"]}')
        zd_params['ZD_Pk_sigma_ratio'] = growth_ratio

    else:
        if not params.get('ZD_Pk_sigma') and not params.get('ZD_Pk_sigma_ratio'):
            raise ValueError("Must specify one of ZD_Pk_file_redshift, sigma_8, ZD_Pk_sigma, or ZD_Pk_sigma_ratio in parameter file")

    # Regardless of the growth amplitude method, calculate f_cluster
    cosm_zinit = Abacus.Cosmology.from_params(params, z='init')
    f_cluster_from_smooth = 1 - params.get('Omega_Smooth',0.)/params['Omega_M']
    if 'ZD_f_cluster' in params:
        # If ZD_f_cluster was already given, check that it matches Omega_Smooth
        if not np.isclose(f_cluster_from_smooth, params['ZD_f_cluster'], rtol=1e-5):
            raise ValueError(f'ZD_f_cluster = {params["ZD_f_cluster"]} does not match 1 - Omega_Smooth/Omega_M = {f_cluster_from_smooth}')
    else:
        zd_params['ZD_f_cluster'] = f_cluster_from_smooth

    # Check that f_growth computed in the EdS approximation with f_cluster is consistent with the value from the cosmology module
    f_growth = cosm_zinit.current.f_growth
    fgrowth_from_fcluster = ((1 + 24*f_cluster_from_smooth)**0.5 - 1)/4
    # TODO: when we add a "just density" mode, skip this check
    #if not np.isclose(f_growth, fgrowth_from_fcluster, rtol=2e-4):
    #    raise ValueError(f'fgrowth_from_fcluster = {fgrowth_from_fcluster} from EdS approximation does not match f_growth = {f_growth} computed from cosmology')

    return zd_params


def ensure_2d(param, nthread=2):
    '''
    Create a copy of the ICs suitable for ingestion by the 2D code,
    if one does not already exist.
    '''

    icdir = param['InitialConditionsDirectory']
    if path.isdir(pjoin(icdir, '2D')):
        return

    transpose_ic.transpose(icdir, param['NP'], param['ICFormat'],
        nthread=nthread)


def run(paramfn, allow_eigmodes_fn_override=False, no_parallel=True):
    '''
    Invokes the zeldovich executable with the given parameter file,
    cleaning up any existing output directories first and also
    copying the input power spectrum to the destination.

    If `allow_eigmodes_fn_override` is set, checks if
    the `ZD_PLT_filename` parameter is valid and sets
    it to the current eigmodes file if not.
    '''
    paramfn = os.path.abspath(paramfn)
    params = GenParam.parseInput(paramfn)

    if path.exists(params['InitialConditionsDirectory']):
        print("Warning: old ICs already exist; removing.")
        shutil.rmtree(params['InitialConditionsDirectory'])

    os.makedirs(params['InitialConditionsDirectory'], exist_ok=True)

    if 'ZD_PLT_filename' in params and not path.isfile(params['ZD_PLT_filename']):
        # If the filename is the same, then the files should be identical and we can silently override
        if path.basename(params['ZD_PLT_filename']) != path.basename(eigmodes_path):
            print('Eigenmodes file "{}" did not exist!'.format(params['ZD_PLT_filename']))
            if allow_eigmodes_fn_override:
                print('Overriding to most recent file: "{}".'.format(eigmodes_path))
            else:
                raise ValueError(allow_eigmodes_fn_override)
        params = GenParam.makeInput(paramfn, paramfn, ZD_PLT_filename=eigmodes_path)

    #if 'ZD_Pk_filename' in params:
    #    shutil.copy(params['ZD_Pk_filename'], pjoin(params['InitialConditionsDirectory'], "input.pow"))
    #else:
    #    assert 'ZD_Pk_powerlaw_index' in params

    ZD_cmd = [pjoin(zeldovich_dir, "zeldovich"), paramfn]

    parallel = params.get('Parallel', False) and not no_parallel

    if parallel:
        ZD_cmd = shlex.split(params['ZD_mpirun_cmd']) + ZD_cmd
    with chdir(zeldovich_dir):
        abacus.call_subprocess(ZD_cmd)

    
def run_override_dirs(parfn, out_parent, new_parfn='abacus_ic_fixdir.par', **kwargs):
    """
    Sometimes we want to regenerate the ICs from a sim run on
    another machine.  Thus, the directories are probably wrong.
    If we have the .par2 file, we can replace the environment
    variables, but if we don't, we need to override the directories,
    which is what this routine is used for.
    
    Note: we generally assume the existence of an info dir,
    which may not be a good assumption.
    """
    
    old_params = GenParam.parseInput(parfn)
    out_parent = path.abspath(out_parent)

    run_kwargs = {}
    if 'no_parallel' in kwargs:
        run_kwargs['no_parallel'] = kwargs.pop('no_parallel')
    
    sim_dir = pjoin(out_parent, old_params['SimName'])
    new_parfn = pjoin(sim_dir, new_parfn)
    ic_dir = pjoin(sim_dir, 'ic')
    os.makedirs(sim_dir, exist_ok=True)
    
    # Copy over info dir and other important files
    # is parfn in the info dir?
    #try: shutil.copytree(pjoin(path.dirname(parfn), path.pardir, 'info'), pjoin(sim_dir, 'info'))
    #except: pass
    # is parfn next to the info dir?
    #try: shutil.copytree(pjoin(path.dirname(parfn), 'info'), pjoin(sim_dir, 'info'))
    #except: pass
    
    #if not path.isdir(pjoin(sim_dir, 'info')):
    #    print("Warning: no info dir found to copy")

    # heuristic replacement of abacus directory, somewhat dangerous
    new_params = dict(old_params).copy()
    for k in new_params:
        p = new_params[k]
        if type(p) is str and 'abacus/' in p:
            new_params[k] = pjoin(os.getenv('ABACUS'), p[p.find('abacus/')+len('abacus/'):])

    # If the eigmodes file doesn't exist, look for it in the zeldovich dir
    if new_params['ZD_qPLT']:
        eigmodes_fn = new_params['ZD_PLT_filename']
        if not path.isfile(eigmodes_fn):
            eigmodes_fn = pjoin(zeldovich_dir, path.basename(eigmodes_fn))
            if not path.isfile(eigmodes_fn):
                print('Warning: original eigenmodes filename "{}" not found.  Falling back to most current eigenmodes.'.format(eigmodes_fn))
                eigmodes_fn = eigmodes_path
            new_params['ZD_PLT_filename'] = eigmodes_fn

    ppd = kwargs.pop('ppd', None)
    if ppd:
        new_params['NP'] = ppd**3
        new_params['ZD_NumBlock'] = ppd//4
        print(f'Using ppd {ppd} override')

    new_params['InitialConditionsDirectory'] = ic_dir

    density = kwargs.pop('density', None)
    if density:
        new_params['ZD_qdensity'] = 1

    del_keys = []
    if kwargs.pop('white', None):
        print('Forcing white spectrum with ZD_Pk_powerlaw_index = 0')
        new_params['ZD_Pk_powerlaw_index'] = 0
        new_params['ZD_Pk_norm'] = 0
        new_params['ZD_density_filename'] = "whitenoise%d"
        new_params.pop('ZD_Pk_filename')
        del_keys += ['ZD_Pk_filename']
    else:
        if 'ZD_Pk_filename' in new_params:
            assert path.isfile(new_params['ZD_Pk_filename'])
        else:
            assert 'ZD_Pk_powerlaw_index' in new_params
            
    # if 
    if ('ZD_Pk_file_redshift' in kwargs or 'InitialRedshift' in kwargs) and 'ZD_Pk_sigma_ratio' in new_params:
        del new_params['ZD_Pk_sigma_ratio']

    for k in kwargs:
        try:
            kwargs[k] = eval(kwargs[k])
        except:
            kwargs[k] = kwargs[k]
    new_params.update(kwargs)

    zd_params = setup_zeldovich_params(new_params)
    new_params.update(zd_params)
        
    new_params = GenParam.makeInput(new_parfn, parfn, del_keywords=del_keys, **new_params)
    
    run(new_parfn, **run_kwargs)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=f'Run the zeldovich code (located in {zeldovich_dir})')
    parser.add_argument('parfile', nargs='+',
        help='The parameter file.  This is usually the same as the .par file for singlestep.')
    parser.add_argument('--out-parent', help="Overrides the parfile InitialConditionsDirectory (i.e. the zeldovich output directory) with PARENT/SimName/ic."
                                             "  Create a new abacus_ic.par with the modified parameters: IC dir; and if don't exist: eigmodes, camb_matterpower", metavar='PARENT')
    parser.add_argument('--show-growth', action='store_true',
        help='Just compute the growth factor from z=0 to z_init from the cosmology in the given parameter file. Does not generate ICs.')
    parser.add_argument('--no-parallel', action='store_true',
        help='Do not invoke the ZD_mpirun_cmd to run zeldovich, despite Parallel = 1 in the parameter file', default=True)
    parser.add_argument('--ppd', help='Override ppd with this value', type=int)
    parser.add_argument('--density', help='Enable density output', action='store_true')
    parser.add_argument('--white', help='Ignore power spectrum file and use ZD_Pk_powerlaw_index=0', action='store_true')

    # this solution for arbitrary args borrowed from https://stackoverflow.com/questions/37367331/is-it-possible-to-use-argparse-to-capture-an-arbitrary-set-of-optional-arguments
    parsed, unknown = parser.parse_known_args()
    for arg in unknown:
        if arg.startswith(("-", "--")):
            arg = arg.split('=')[0]
            parser.add_argument(arg, type=str)

    args = parser.parse_args()
    args = vars(args)
    out_parent = args.pop('out_parent')
    ppd = args.pop('ppd')
    density = args.pop('density')
    show_growth = args.pop('show_growth')
    
    for parfn in args.pop('parfile'):
        if show_growth:
            par = GenParam.parseInput(parfn)
            
            print(Abacus.Cosmology.from_params(par, 'init').current.growth)
            
            sigma8_zinit = calc_sigma8(par, z='init')

            # TODO: technically sigma_8 probably doesn't have to be at z=0
            print("sigma_8(z={zinit:g}) = {}\nsigma_8(z=0) = {}\nGrowth ratio D(z={zinit:g})/D(z=0) = {}\n".format(sigma8_zinit, par['sigma_8'], sigma8_zinit/par['sigma_8'], zinit=par['InitialRedshift'],))
        else:
            if out_parent or ppd or density:  # TODO
                run_override_dirs(parfn, out_parent, ppd=ppd, density=density, **args)
            else:
                run(parfn, **args)
