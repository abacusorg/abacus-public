#!/usr/bin/env python
'''
This module is the main entry point for running an Abacus simulation.
For a complete example of a simulation setup that demonstrates how
to use this module, see the `abacus/Production/Example` directory.

This module is mostly a wrapper for the `singlestep` binary. It
handles moving the state folders around, etc.  It can be imported
as a module or directly invoked as a script:

Script usage:
    See abacus/Production/Example

Command line usage:
    abacus.py [-h] [--clean] [--erase-ic] parfile

    Run this sim.

    positional arguments:
      parfile     The parameter file

    optional arguments:
      -h, --help  show this help message and exit
      --clean     Erase the working directory and start over. Otherwise, continue
                  from the existing state. Always preserves the ICs unless
                  --erase-ic.
      --erase-ic  Remove the ICs if they exist.

'''

import os
abacuspath = os.getenv("ABACUS")
if not abacuspath:
    print("Error: Please define $ABACUS to be the absolute path to your abacus distribution")
    sys.exit(1)
    
import os.path as path
from os.path import join as pjoin, normpath, basename, dirname
import subprocess
import sys
import shutil
import ctypes as ct
import argparse
import numpy as np
from glob import glob

from .InputFile import InputFile
from . import Tools
from . import GenParam
from . import zeldovich
from Abacus.Cosmology import AbacusCosmo
import Abacus

EXIT_REQUEUE = 200
site_param_fn = pjoin(abacuspath, 'Production', 'site_files', 'site.def')
directory_param_fn = pjoin(abacuspath, 'Abacus', 'directory.def')

def run(parfn='abacus.par2', config_dir=path.curdir, maxsteps=10000, clean=False, erase_ic=False, output_parfile=None, use_site_overrides=False, override_directories=False, **param_kwargs):
    """
    The main entry point for running a simulation.  Initializes directories,
    copies parameter files, and parses parameter files as appropriate, then
    invokes `singlestep()`.

    Important files should generally be placed in the 'info' directory in
    `config_dir`.  This directory will be copied alongside outputs, data
    products, etc.  `parfn` can also be in the info dir.
    
    Parameters
    ----------
    parfn: str, optional
        The name of the Abacus parameter (.par or .par2) file relative to
        `config_dir`.  Default: abacus.par2
    config_dir: str, optional
        The directory containing the parameter files for this simulation.
        Default: path.curdir (current directory)
    maxsteps: int, optional
        The maximum number of times to invoke singlestep.  Default: 10000.
    clean: bool, optional
        Whether to continue from an existing state or erase it and start over.
        Never erases ICs; only erase_ic does that.  Default: False
    erase_ic: bool, optional
        Erase the initial conditions directory and generate new ICs.
        Default: False
    output_parfile: str, optional
        Filename of the processed parameters file.  Default of None places it
        in the same directory as `parfile`.
    use_site_overrides: bool, optional
        Use the `site_param_fn` file to overwrite tuning parameters
        (like `OMP_NUM_THREADS`) in the given parameter file.  Default: False.
    override_directories: bool, optional
        Discard the directories (like `WorkingDirectory`) listed in the parameter
        file and use the system defaults instead (from environment vars like
        $ABACUS_TMP).  Also allows some zeldovich-PLT paths, like the eigenmodes
        file, to be overwritten.  Default: False.
    param_kwargs: dict, optional
        Extra parameters to pass to the InputFile parser.  Useful for changing
        parameters on-the-fly for e.g. testing different parameters.
        
    Returns
    -------
    retcode: int
        The singlestep return value
    """
    
    parfn = pjoin(config_dir, parfn)
    pardir = dirname(parfn)
    if not output_parfile:
        output_parfile = pjoin(pardir, 'abacus.par')
        
    try:
        maxsteps = int(maxsteps)
        assert(maxsteps >= 0)
    except:
        maxsteps = 10000

    params = preprocess_params(output_parfile, parfn, use_site_overrides=use_site_overrides,
                override_directories=override_directories, **param_kwargs)
    
    basedir = params['WorkingDirectory']
    logdir = params['LogDirectory']
    icdir = params['InitialConditionsDirectory']
    outdir = params['OutputDirectory']

    # If we requested a resume, but there is no state, assume we are starting fresh
    if not clean:
        if path.exists(basedir):
            print('Resuming from existing state.')
        else:
            print('Resume requested but no state exists.  Creating one.')
            clean = True

    if erase_ic and path.exists(icdir):
        print('Erasing IC dir')
        shutil.rmtree(icdir)

    # If not resuming, erase everything but the ICs
    if clean:
        clean_dir(params['WorkingDirectory'], preserve=None if erase_ic else icdir)
        try:
            clean_dir(params['WorkingDirectory2'], preserve=None if erase_ic else icdir)
        except KeyError: pass
        clean_dir(params['MultipoleDirectory'])
        try: clean_dir(params['MultipoleDirectory2'])
        except KeyError: pass
        clean_dir(params['TaylorDirectory'])
        try: clean_dir(params['TaylorDirectory2'])
        except KeyError: pass
        # TODO: deleting the parent ad-hoc is dangerous!  should switch to ConvolutionWorkingDirectory
        clean_dir(dirname(normpath(params['TaylorDirectory'])), preserve=None if erase_ic else icdir)
        clean_dir(dirname(normpath(params['MultipoleDirectory'])), preserve=None if erase_ic else icdir)

        if path.exists(logdir):
            shutil.rmtree(logdir)

        print("Erased previous working directories.")
            
    if not path.exists(basedir):
        os.makedirs(basedir)
    if not path.exists(logdir):
        os.makedirs(logdir)
    
    info_in_path = pjoin(config_dir, 'info')
    try:
        shutil.copy(output_parfile, '{basedir:s}/'.format(basedir=basedir))
        info_out_path = pjoin('{basedir:s}'.format(basedir=basedir), 'info')
        if path.isdir(info_out_path):
            shutil.rmtree(info_out_path)
        shutil.copytree(info_in_path, info_out_path)
    except (shutil.Error, OSError):
        pass  # input == output

    # Copy over all of the input files to the output folder for reference
    if not path.exists(outdir):
        os.makedirs(outdir)
    try:
        shutil.copy(output_parfile, '{outdir:s}/'.format(outdir=outdir))
    except shutil.Error:
        pass  # same file?
    
    info_out_path = pjoin('{outdir:s}'.format(outdir=outdir), 'info')
    try:
        shutil.rmtree(info_out_path)
    except:
        pass
        
    try:
        shutil.copytree(info_in_path, info_out_path)
    except (shutil.Error, OSError):
        pass  # no info dir
        
    with Tools.chdir(basedir):
        # The parfile is directly stored in the basedir
        output_parfile = basename(output_parfile)
        if clean and not path.exists(icdir):
            zeldovich.run(output_parfile, allow_eigmodes_fn_override=override_directories)
        elif clean:
            print('Reusing existing ICs')
        retval = singlestep(output_parfile, maxsteps, AllowIC=clean)

    return retval


def clean_dir(bd, preserve=None, rmdir_ifempty=True):
    '''
    Remove the contents of a directory except for the path
    specified in the "preserve" argument.  Optionally remove
    the directory itself if it's empty.
    '''

    if path.exists(bd):
        # Erase everything in basedir except the ICs
        for d in os.listdir(bd):
            d = pjoin(bd, d)
            try:
                if preserve and path.samefile(d, preserve):
                    continue
            except OSError:
                pass  # icdir may not exist
            try:
                shutil.rmtree(d)
            except OSError:
                os.remove(d)  # may have been a file, not a dir
        if rmdir_ifempty:
            try:
                os.rmdir(bd)
            except OSError:
                pass

def MakeDerivatives(param, floatprec=False): #look for the derivatives required for param and make them if they don't exist
    if not path.exists(param.DerivativesDirectory):
        os.makedirs(param.DerivativesDirectory)
    for ADfn in ['AD32_001.dat','AD32_002.dat','AD32_003.dat','AD32_004.dat','AD32_005.dat','AD32_006.dat','AD32_007.dat','AD32_008.dat', 'AD32_016.dat']: #note! added this in to run Ewald test for far_radius 1-8,16 
        source_ADfn = pjoin(abacuspath, "Derivatives", ADfn)
        if not path.isfile(pjoin(param.DerivativesDirectory, ADfn)):
            shutil.copy(source_ADfn,param.DerivativesDirectory)
    
    fnfmt32 = "fourierspace_float32_%d_%d_%d_%d"%(param.CPD,param.Order,param.NearFieldRadius,param.DerivativeExpansionRadius) + "_%d"
    fnfmt64 = "fourierspace_%d_%d_%d_%d"%(param.CPD,param.Order,param.NearFieldRadius,param.DerivativeExpansionRadius) + "_%d"
            
    derivativenames32 = []; derivativenames64 = []
    for i in range(0,param.CPD//2+1):
        derivativenames32 += [pjoin(param.DerivativesDirectory, fnfmt32%i)]
        derivativenames64 += [pjoin(param.DerivativesDirectory, fnfmt64%i)]
    derivativenames = derivativenames32 if floatprec else derivativenames64
    fnfmt = fnfmt32 if floatprec else fnfmt64

    if not all(map(path.isfile, derivativenames)):
        print('Could not find derivatives in "' + param.DerivativesDirectory + '" ...Creating them')
        print("Error was on file pattern '{}'".format(fnfmt))

        if floatprec:
            # first make sure the derivatives exist in double format
            MakeDerivatives(param, floatprec=False)
            # now make a 32-bit copy
            from . import convert_derivs_to_float32
            for dpath in derivativenames64:
                convert_derivs_to_float32.convert(dpath)
        else:
            with Tools.chdir(param.DerivativesDirectory):
                try:
                    os.remove("./farderivatives")
                except OSError:
                    pass
                subprocess.check_call([pjoin(abacuspath, "Derivatives", "CreateDerivatives"),
                str(param.CPD),str(param.Order),str(param.NearFieldRadius),str(param.DerivativeExpansionRadius)])
    
    
def default_parser():
    """
    A default command-line interface that other scripts may want to use.

    Example
    -------
    parser = abacus.default_parser()
    parser.add_argument(...)  # optional
    args = parser.parse_args()
    retcode = abacus.run('abacus.par2', **vars(args))
    
    """
    parser = argparse.ArgumentParser(description='Run this sim.')
    parser.add_argument('--clean', action='store_true', help="Erase the working directory and start over.  Otherwise, continue from the existing state.  Always preserves the ICs unless --erase-ic.")
    parser.add_argument('--erase-ic', action='store_true', help="Remove the ICs if they exist.")
    parser.add_argument('--maxsteps', type=int, help="Take at most N steps", metavar='N')
    return parser


def preprocess_params(output_parfile, parfn, use_site_overrides=False, override_directories=False, **param_kwargs):
    '''
    There are a number of runtime modifications we want to make
    to the Abacus .par2 file as we parse it.  These typically
    require computation of some sort and can't be done with static
    parsing.

    TODO: probably want an explicit SplitMultipoles parameter

    1) Compute ZD_Pk_sigma from sigma_8
    2) If $ABACUS_SSD2 is defined, define the params {Multipole,Taylor}Directory2
    3) If $ABACUS_TMP2 is defined and SloshState is True, define WorkingDirectory2
    4) If `use_site_overrides` is set, load params from the site.def file
    5) If `override_directories` is set, overwrite directories with defaults from the environment
    '''
    # We will build this dict in order of precedence
    _param_kwargs = param_kwargs.copy()

    if use_site_overrides:
        site = InputFile(site_param_fn)
        for k in site.keys():
            # explicit param_kwargs take precedence over site files
            if k not in param_kwargs:
                param_kwargs[k] = site[k]

    params = GenParam.makeInput(output_parfile, parfn, **param_kwargs)

    if override_directories:
        dirs = GenParam.parseInput(directory_param_fn, varreplace_values=params)
        for k in dirs.keys():
            # explicit param_kwargs take precedence over directory overrides
            if k not in param_kwargs:
                param_kwargs[k] = dirs[k]
    
    if 'sigma_8' in params:
        sigma8_at_zinit = zeldovich.calc_sigma8(params)
        if 'ZD_Pk_sigma' in params:
            assert np.isclose(params['ZD_Pk_sigma'], sigma8_at_zinit, rtol=1e-4),\
                "ZD_Pk_sigma ({:f}) calculated from sigma_8 ({:f}) conflicts with the provided value of ZD_Pk_sigma ({:f})!".format(sigma8_at_zinit, params['sigma_8'], params['ZD_Pk_sigma'])
        param_kwargs['ZD_Pk_sigma'] = sigma8_at_zinit
        params = GenParam.makeInput(output_parfile, parfn, **param_kwargs)

    if params.get('SloshState',False):
        if 'ABACUS_TMP2' in os.environ:
            # WorkingDirectory2 is never processed by Abacus.
            # It will be intercepted by the sloshing logic
            if 'WorkingDirectory2' not in params:
                param_kwargs['WorkingDirectory2'] = pjoin(os.environ['ABACUS_TMP2'], params['SimName'])

        if 'ABACUS_SSD2' in os.environ:
            if 'ConvolutionWorkingDir2' not in params:
                param_kwargs['ConvolutionWorkingDir2'] = pjoin(os.environ['ABACUS_SSD2'], params['SimName'])
    else:
        if 'ABACUS_SSD2' in os.environ:
            if 'MultipoleDirectory2' not in params:
                param_kwargs['MultipoleDirectory2'] = pjoin(os.environ['ABACUS_SSD2'], params['SimName'], 'multipole2')
            if 'TaylorDirectory2' not in params:
                param_kwargs['TaylorDirectory2'] = pjoin(os.environ['ABACUS_SSD2'], params['SimName'], 'taylor2')

    params = GenParam.makeInput(output_parfile, parfn, **param_kwargs)

    return params


def check_multipole_taylor_done(param, state, kind):
    """
    Checks that all the multipoles or Taylors exist and have the expected
    file size. `kind` must be 'Multipole' or 'Taylor'.
    """

    assert kind in ['Multipole', 'Taylor']
    rml = (param.Order+1)**2
    sizeof_Complex = 16 if state.DoublePrecision else 8
    size = param.CPD*(param.CPD+1)/2*rml*sizeof_Complex

    even_dir = param[kind+'Directory']
    odd_dir = param.get(kind+'Directory2', param[kind+'Directory'])
    if kind == 'Multipole':
        kind += 's'  # should consider standardizing our pluralization
    try:
        even = all(path.getsize(pjoin(even_dir, kind+'_{:04d}'.format(i))) == size
            for i in range(0,param.CPD,2))
        odd = all(path.getsize(pjoin(odd_dir, kind+'_{:04d}'.format(i))) == size
            for i in range(1,param.CPD,2))
    except FileNotFoundError:
        return False

    return even and odd


def setup_singlestep_env(param):
    """
    Assign OpenMP threads to cores using the OMP_PLACES
    and OMP_PROC_BIND environment variables.  This has
    to be done in Python instead of singlestep because
    OpenMP offers no C/C++ interfaces to do this.
    
    TODO: Setting intelligent defaults for these would
    require detecting the NUMA architecture of the machine
    """
    
    singlestep_env = os.environ.copy()
    if 'OMP_PLACES' in param:
        singlestep_env['OMP_PLACES'] = param.OMP_PLACES
    if 'OMP_PROC_BIND' in param:
        singlestep_env['OMP_PROC_BIND'] = param.OMP_PROC_BIND
        
    return singlestep_env
    
    
def setup_convolution_env(param):
    """
    Same as singlestep_env, but for the convolution
    """
    
    convolution_env = os.environ.copy()
    if 'Conv_OMP_PLACES' in param:
        convolution_env['OMP_PLACES'] = param.Conv_OMP_PLACES
    if 'Conv_OMP_PROC_BIND' in param:
        convolution_env['OMP_PROC_BIND'] = param.Conv_OMP_PROC_BIND
        
    return convolution_env

    
def setup_state_dirs(paramfn):
    '''
    This function is called once before the singlestep loop begins.
    Its main purpose is to setup the state directories, including
    the symlinks if state sloshing is requested.

    When state sloshing, everything looks normal to the singlestep
    binary, except "read" and "write" are actually symlinks
    to the "state1" and "state2" directories elsewhere.

    When sloshing the state, multipoles & taylors are sloshed
    as well.  The Convolution

    TODO: there will probably be times when we want to slosh
    the main state but split the MT
    '''
    params = GenParam.makeInput(paramfn, paramfn)
    OverwriteState = params.get('OverwriteState',False)
    
    # Normal operation: everything is under the WorkingDirectory or specified individually
    # Can't specify both WorkingDir and ReadDir
    if 'WorkingDirectory' in params:
        assert 'ReadStateDirectory' not in params
        assert 'WriteStateDirectory' not in params
        assert 'PastStateDirectory' not in params
        past = pjoin(params['WorkingDirectory'], "past")
        read = pjoin(params['WorkingDirectory'], "read")
        write = pjoin(params['WorkingDirectory'], "write")
    else:
        past = params['PastStateDirectory']
        read = params['ReadStateDirectory']
        write =params['WriteStateDirectory']

    past = normpath(past)
    read = normpath(read)
    write = normpath(write)

    taylors = normpath(params['TaylorDirectory'])
    multipoles = normpath(params['MultipoleDirectory'])

    if OverwriteState:
        write = read
        past = None

    # Check if the write state already exists
    if path.isfile(pjoin(write, "state")) and not OverwriteState:
        answer = input('\nWrite state "{}" exists and would be overwritten. Do you want to delete it? (y/[n]) '.format(write))
        if answer == "y":
            shutil.rmtree(write)
        else:
            raise ValueError("Cannot proceed if write state exists!")
    
    # Set up symlinks to slosh the state
    if params.get('SloshState',False):
        assert not OverwriteState, 'SloshState is mutually exclusive with OverwriteState!'
        # Should be populated manually or from ABACUS_TMP2
        assert 'WorkingDirectory' in params  # need a working directory to put the symlinks
        assert 'WorkingDirectory2' in params, "Must specify WorkingDirectory2 via $ABACUS_TMP2 or directly in the parameter file!"
        assert 'ConvolutionWorkingDir2' in params, "Must specify ConvolutionWorkingDir2 via $ABACUS_SSD2 or directly in the parameter file!"
        # this signals parity splitting to singlestep
        # TODO: make another parameter to signal this?
        assert 'MultipoleDirectory2' not in params and 'TaylorDirectory2' not in params

        # we use dirname(taylors) as ConvolutionWorkingDir
        assert dirname(taylors) == dirname(multipoles)

        # If the links don't exist, create them and the underlying dirs
        # If they do exist, don't touch them
        def make_link(link, target):
            if path.exists(link):
                # link already existing as a directory insetad of a link is an error!
                assert path.islink(link)
            else:
                try: os.makedirs(target)
                except OSError: pass  # already exists; that's fine
                os.symlink(target, link)

        make_link(read, pjoin(params['WorkingDirectory'], 'state1'))
        make_link(write, pjoin(params['WorkingDirectory2'], 'state2'))
        make_link(taylors, pjoin(dirname(taylors), 'convstate1'))
        make_link(multipoles, pjoin(params['ConvolutionWorkingDir2'], 'convstate2'))

        # Technically we could keep two past states (one in each slosh directory)
        # but in practice this may be too many copies of the state
        past = None
    else:
        if not path.isdir(read):
            os.makedirs(read)
        if not path.isdir(write):
            os.makedirs(write)
        if past and not path.isdir(past):
            os.makedirs(past)

    # Create dirs
    for d in ['MultipoleDirectory', 'MultipoleDirectory2', 'TaylorDirectory', 'TaylorDirectory2']:
        if d in params:
            try:
                os.makedirs(params[d])
            except FileExistsError:
                pass

    return read,write,past,multipoles,taylors


def move_state_dirs(read, write, past, multipoles, taylors):
    '''
    At the end of every timestep, we to move the "write" state
    to the read directory, move "read" to "past", and discard
    "past".

    If we're sloshing the state, we instead want to swap the
    symlinks for `read` and `write`.  We also want to swap
    the symlink for the Taylors (but not the multipoles because
    the convove hasn't happened yet).

    TODO: double check this still works when overwriting the state
    '''

    # move read to past, write to read
    # let's assume sloshing never uses past; otherwise, we need a past in each working dir
    if past is None:
        if not path.samefile(read, write):
            try:
                shutil.rmtree(read)
            except OSError:
                # can't rmtree symlink! but we at least want to delete the contents
                clean_dir(read, preserve=None, rmdir_ifempty=False)
                pass
    else:
        try:
            shutil.rmtree(past)
        except OSError:
            pass  # past doesn't exist
        try:
            os.rename(read, past)
        except OSError:
            pass  # read doesn't exist, must be IC state

    readtmp = read + '_tmp'
    os.rename(write, readtmp)
    try:
        os.rename(read, write)
    except FileNotFoundError:
        # if not a symlink, read was moved to past (or deleted) already
        os.makedirs(write)
    os.rename(readtmp, read)


def remove_MT(md, pattern, rmdir=False):
    # Quick util to clear multipoles/taylors
    for mfn in glob(pjoin(md, pattern)):
        os.remove(mfn)
    if rmdir:
        try:
            os.rmdir(md)
        except OSError:
            pass  # remove dir if empty; we re-create at each step in setup_dirs


def singlestep(paramfn, maxsteps, AllowIC=False, stopbefore=-1):
    """
    Run a number of Abacus timesteps by invoking the `singlestep` executable.

    Parameters
    ----------
    paramfn: str
        File name of the parameter file for this simulation
    maxsteps: int
        Run (at most) `maxsteps` timesteps.
    AllowIC: bool, optional
        Make the first timestep an IC ingestion step.
    stopbefore: int, optional
        Run abacus until it is about to call "singlestep" for step `stopbefore`.
        Useful for attaching a debugger to a certain step without getting
        interference from Python.

    Returns
    -------
    retcode: int
        Returns 0 or 1 for general success or failure;
        or EXIT_REQUEUE if we finish cleanly but have not reached the final z
    """
    finished = False

    if maxsteps == 0:
        return 0

    if not path.exists(paramfn):
        raise ValueError('Parameter file "{:s}" is not accessible'.format(paramfn))
        
    param = InputFile(paramfn)

    # Set up backups
    try:
        backups_enabled = param['BackupStepInterval'] > 0
    except:
        backups_enabled = False  # No BackupIntervalSteps parameter

    #if not backups_enabled:
    #    print 'Warning: automatic state backup not enabled.  Set "BackupStepInterval" and "BackupDirectory" in the parameters file.'
        
    # Set up OpenMP singlestep environment
    singlestep_env = setup_singlestep_env(param)
    convolution_env = setup_convolution_env(param)

    status_log_fn = pjoin(param.OutputDirectory, 'status.log')
    if not os.path.isfile(status_log_fn):
        with open(status_log_fn, 'w') as status_log:
            status_log.write("Step  Redshift\n")

    print("Using parameter file %s and working directory %s."%(paramfn,param.WorkingDirectory))

    # floatprec=True effectively guarantees both precisions are available
    # TODO: how to tell if Abacus was compiled in double precision?
    MakeDerivatives(param, floatprec=True)

    # Create directories
    for d in ['LogDirectory', 'OutputDirectory', 'GroupDirectory']:
        if d in param and param[d] and not path.isdir(param[d]):
            os.makedirs(param[d])

    print("Beginning abacus steps:")

    if AllowIC:
        print("Ingesting IC from "+param.InitialConditionsDirectory+" ... Skipping convolution")
        if path.exists(pjoin(param.InitialConditionsDirectory, "input.pow")):
            shutil.copy(pjoin(param.InitialConditionsDirectory, "input.pow"), pjoin(param.LogDirectory, "input.pow"))

    # Finalize the directories, possibly setting up symlinks for state sloshing
    read,write,past,multipoles,taylors = setup_state_dirs(paramfn)

    for i in range(maxsteps):
        try:
            read_state = InputFile(pjoin(read,"state"))
            ConvDone = check_multipole_taylor_done(param, read_state, kind='Taylor')
            stepnum = read_state.FullStepNumber+1
        except IOError:
            assert AllowIC, "If there is no read state, this must be an IC step!"
            ConvDone = True  # No need to convolve for an IC step
            stepnum = 0

        # Do the convolution
        if not ConvDone:
            # Now check if we have all the multipoles the convolution will need
            if not check_multipole_taylor_done(param, read_state, kind='Multipole'):
                # Invole multipole recovery mode
                print("Warning: missing multipoles! Performing multipole recovery for step {:d}".format(i))
                
                # Build the make_multipoles executable
                with Tools.chdir(pjoin(abacuspath, "singlestep")):
                    subprocess.check_call(['make', 'make_multipoles'])

                # Execute it
                subprocess.check_call([pjoin(abacuspath, "singlestep", "make_multipoles"), paramfn], env=singlestep_env)
                save_log_files(param.LogDirectory, 'step{:04d}.make_multipoles'.format(read_state.FullStepNumber))
                print('\tFinished multipole recovery for read state {}.'.format(read_state.FullStepNumber))

            # Swap the Taylors link.  In effect, this will place the Taylors on the same disk as the multipoles.
            # But that's what we want for sloshing: this was the write disk, so now it will be the upcoming read disk
            # We'll swap the multipoles as soon as we convolve
            if path.islink(taylors):
                taylors_convstate = os.readlink(taylors)
                os.unlink(taylors)
                os.symlink(os.readlink(multipoles), taylors)
            
            print("Performing convolution for step {:d}".format(stepnum))
            subprocess.check_call([pjoin(abacuspath, "Convolution", "ConvolutionDriver"), paramfn], env=convolution_env)
                
            if not check_multipole_taylor_done(param, read_state, kind='Taylor'):
                raise ValueError("Convolution did not complete")
            
            convlogs = glob(pjoin(param.LogDirectory, 'last.conv*'))
            for cl in convlogs:
                shutil.move(cl, cl.replace('last', 'step{:04d}'.format(read_state.FullStepNumber+1)))

            # Warning: Convolution won't work if MultipoleDirectory is the write (or read) state
            # because the states get moved after multipole generation but before convolution.
            remove_MT(param.MultipoleDirectory, 'Multipole*')
            if 'MultipoleDirectory2' in param:
                remove_MT(param.MultipoleDirectory2, 'Multipole*')

            # Now we can move the multipole symlinks
            if path.islink(multipoles):
                os.unlink(multipoles)
                os.symlink(taylors_convstate, multipoles)

            print("\t Finished convolution for state %d"%(read_state.FullStepNumber+1))

        #run singlestep
        print("Running singlestep for step {:d}".format(stepnum))
        if i == stopbefore:
            print('stopbefore = %d was specified; stopping before calling singlestep'%i)
            return 0
        subprocess.check_call([pjoin(abacuspath, "singlestep", "singlestep") , paramfn, str(int(AllowIC))], env=singlestep_env)
        
        # In profiling mode, we don't move the states so we can immediately run the same step again
        if param.get('ProfilingMode', False):
            break
            
        # Check that write/state was written as a test of success
        if not path.exists(pjoin(write, "state")):
            raise ValueError("Singlestep did not complete!")
        write_state = InputFile(pjoin(write, "state"))

        # Update the status log
        with open(status_log_fn, 'a') as status_log:
            status_log.write('{step:4d}  {z:6.3f}\n'.format(step=stepnum, z=write_state.Redshift))

        # save the log and timing files under this step number
        save_log_files(param.LogDirectory, 'step{:04d}'.format(write_state.FullStepNumber))
                        
        shutil.copy(pjoin(write, "state"), pjoin(param.LogDirectory, "step{:04d}.state".format(write_state.FullStepNumber)))
        print(( "\t Finished state {:d}. a = {:.4f}, dlna = {:.3g}, rms velocity = {:.3g}".format(
            write_state.FullStepNumber, write_state.ScaleFactor,
            write_state.DeltaScaleFactor/(write_state.ScaleFactor-write_state.DeltaScaleFactor),
            write_state.RMS_Velocity)))

        # Delete the Taylors right away; want to avoid accidentally
        # running with inconsistent positions and Taylors
        remove_MT(param.TaylorDirectory, 'Taylor_*')
        if 'TaylorDirectory2' in param:
            remove_MT(param.TaylorDirectory2, 'Taylor_*')

        # Check if we need to back up the write state
        if backups_enabled and write_state.FullStepNumber > 0 and write_state.FullStepNumber % param['BackupStepInterval'] == 0:
            backup_dir = pjoin(param['BackupDirectory'], 'backup')
            print('Backing up state to {}'.format(backup_dir))
            tmp_backup = pjoin(param['BackupDirectory'], 'backup_inprogress')
            shutil.copytree(write, tmp_backup)
            if path.isdir(backup_dir):
                shutil.rmtree(backup_dir)
            shutil.move(tmp_backup, backup_dir)
        
        # Make a power spectrum
        dfn = pjoin(read, "density")
        if path.exists(dfn):
            # This import can be slow, so we only do it as necessary
            from Abacus.Analysis.PowerSpectrum import PowerSpectrum as PS
            ngrid = param.PowerSpectrumN1d

            psdtype = np.float64 if read_state.DoublePrecision else np.float32
            density = np.empty((ngrid, ngrid, ngrid + 2), dtype=psdtype)
            density[...,:-2] = np.fromfile(dfn, dtype=psdtype).reshape(ngrid,ngrid,ngrid)
            k,P,nb = PS.FFTAndBin(density, param.BoxSize, inplace=True)
            # we use the write state step number here, even though the particles are positioned at the read state positions
            # this is consistent with the time slice behavior
            shutil.copy(dfn, pjoin(param.LogDirectory, "step{:04d}.density".format(write_state.FullStepNumber)))
            np.savez(pjoin(param.LogDirectory, "step{:04d}.pow".format(write_state.FullStepNumber)), k=k,P=P,nb=nb)
            os.remove(dfn)
        
        # Now shift the states down by one
        taylors_convstate = move_state_dirs(read, write, past, multipoles, taylors)
        
        # check if we've reached the stopping redshift
        # Note that FinalRedshift always trumps,
        # whether it is before or after the TimeSliceZ array
        try:
            finalz = param.FinalRedshift
        except AttributeError:
            # param.FinalRedshift wasn't set, so we'll look for the minimum in the TimeSliceRedshifts array
            if (param.nTimeSlice>0):
                try:
                    finalz = min(param.TimeSliceRedshifts[:param.nTimeSlice])
                except TypeError:  # param.TimeSliceRedshifts is probably just a single number
                    finalz = float(param.TimeSliceRedshifts)
            else:
                finalz = 0.0

        # This logic is deliberately consistent with singlestep.cpp
        # If this is an IC step then we won't have read_state
        if not AllowIC and np.abs(read_state.Redshift - finalz) < 1e-12 and read_state.LPTStepNumber == 0:
            print("Final redshift reached; terminating normally.")
            finished = True
            break
        
        AllowIC = False

        ### end singlestep loop

    # If there is more work to be done, signal that we are ready for requeue
    if not finished:
        return EXIT_REQUEUE

    return 0

def save_log_files(logdir, newprefix, oldprefix='lastrun'):
    for logfn in os.listdir(logdir):
        if logfn.startswith(oldprefix):
            newname = logfn.replace(oldprefix, newprefix, 1)
            shutil.move(pjoin(logdir, logfn), pjoin(logdir, newname))


if __name__ == '__main__':
    parser = default_parser()
    parser.add_argument('parfile', help="The parameter file")
    args = parser.parse_args()
    args = vars(args)
    retcode = run(args.pop('parfile'), **args)
    sys.exit(retcode)
