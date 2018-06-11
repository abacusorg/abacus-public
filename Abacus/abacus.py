#!/usr/bin/env python
'''
Wrapper for the abacus singlestep binary. Handles moving the state folders around, etc.
Can also be imported as a python module.

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
from os.path import join as pjoin
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

def run(parfn='abacus.par2', config_dir=path.curdir, maxsteps=10000, clean=False, erase_ic=False, output_parfile=None, **param_kwargs):
    """
    The main entry point for running a simulation.  Initializes directories,
    copies parameter files, and parses parameter files as appropriate, then
    invokes singlestep.

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
    param_kwargs: dict, optional
        Extra parameters to pass to the InputFile parser.  Useful for changing
        parameters on-the-fly for e.g. testing different parameters.
        
    Returns
    -------
    retcode: int
        The singlestep return value
    """
    
    parfn = pjoin(config_dir, parfn)
    pardir = path.dirname(parfn)
    if not output_parfile:
        output_parfile = pjoin(pardir, 'abacus.par')
        
    try:
        maxsteps = int(maxsteps)
        assert(maxsteps >= 0)
    except:
        maxsteps = 10000

    params = preprocess_params(output_parfile, parfn, **param_kwargs)
    
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

    def clean_basedir(bd, preserve=icdir):
        if path.exists(bd):
            # Erase everything in basedir except the ICs
            for d in os.listdir(bd):
                d = pjoin(bd, d)
                try:
                    if not erase_ic and path.samefile(d, preserve):
                        continue
                except OSError:
                    pass  # icdir may not exist
                try:
                    shutil.rmtree(d)
                except OSError:
                    os.remove(d)  # may have been a file, not a dir
            try:
                os.rmdir(bd)
            except OSError:
                pass
            print("Erased previous working directory")

    # If not resuming, erase everything but the ICs
    if clean:
        clean_basedir(params['WorkingDirectory'])
        try:
            clean_basedir(params['WorkingDirectory2'])
        except KeyError: pass
                
        if path.exists(logdir):
            shutil.rmtree(logdir)
            print("Erased previous log directory.")
            
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
        output_parfile = path.basename(output_parfile)
        if clean and not path.exists(icdir):
            zeldovich.run(output_parfile)
        elif clean:
            print('Reusing existing ICs')
        retval = singlestep(output_parfile, maxsteps, AllowIC=clean)

    return retval


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


def preprocess_params(output_parfile, parfn, **param_kwargs):
    '''
    There are a number of runtime modifications we want to make
    to the Abacus .par2 file as we parse it.  These typically
    require computation of some sort and can't be done with static
    parsing.

    1) Compute ZD_Pk_sigma from sigma_8
    2) If $ABACUS_SSD2 is defined, define the params {Multipole,Taylor}Directory2
    3) If $ABACUS_TMP2 is defined and SloshState is True, define WorkingDirectory2
    '''

    params = GenParam.makeInput(output_parfile, parfn, **param_kwargs)
    
    if 'sigma_8' in params:
        sigma8_at_zinit = zeldovich.calc_sigma8(params)
        if 'ZD_Pk_sigma' in params:
            assert np.isclose(params['ZD_Pk_sigma'], sigma8_at_zinit, rtol=1e-4),\
                "ZD_Pk_sigma ({:f}) calculated from sigma_8 ({:f}) conflicts with the provided value of ZD_Pk_sigma ({:f})!".format(sigma8_at_zinit, params['sigma_8'], params['ZD_Pk_sigma'])
        param_kwargs['ZD_Pk_sigma'] = sigma8_at_zinit
        params = GenParam.makeInput(output_parfile, parfn, **param_kwargs)

    if 'ABACUS_SSD2' in os.environ:
        if 'MultipoleDirectory2' not in params:
            param_kwargs['MultipoleDirectory2'] = pjoin(os.environ['ABACUS_SSD2'], params['SimName'], 'multipole2')
        if 'TaylorDirectory2' not in params:
            param_kwargs['TaylorDirectory2'] = pjoin(os.environ['ABACUS_SSD2'], params['SimName'], 'taylor2')

    if params.get('SloshState',False):
        if 'ABACUS_TMP2' in os.environ:
            # WorkingDirectory2 is never processed by Abacus.
            # They will be intercepted by the sloshing logic
            if 'WorkingDirectory2' not in params:
                param_kwargs['WorkingDirectory2'] = pjoin(os.environ['ABACUS_TMP2'], params['SimName'])

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

    
def setup_state_dirs(paramfn, newparam_suffix='step{step:04d}'):
    '''
    This function is called before every invocation of singlestep.
    Its main purpose is to slosh the state directories if requested
    by swapping ReadStateDirectory and WriteStateDirectory
    each step.

    To do this, paramfn is ingested, and a modified version is
    written back to disk.  The written version doesn't use
    WorkingDirectory and instead directly sets the Read and Write
    dirs.

    When sloshing the state, multipoles & taylors are sloshed
    as well.  MultipoleDirectory2 is interpreted as the slosh
    directory, not the parity-splitting directory.  For that
    reason, MultipoleDirectory2 will not be written to the param
    file.

    TODO: there will probably be times when we want to slosh
    the main state but split the MT

    '''
    params = GenParam.makeInput(paramfn, paramfn)
    OverwriteState = params.get('OverwriteState',False)

    # Parameters to modify in the output file
    # If these remain empty, the input will be unmodified and the filename will be unchanged
    delparams = []
    newparams = {}
    newconvparams = {}
    
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
        write = params['WriteStateDirectory']
    
    # read is the current read state, read2 is the read state of the next singlestep invocation
    read2 = read
    write2 = write

    # Slosh the state
    if params.get('SloshState',False):
        assert not OverwriteState  # Sloshing is mutually exclusive with overwriting!
        # Should be populated manually or from ABACUS_TMP2
        assert 'WorkingDirectory2' in params

        def get_stepnum(dirs):
            stepnum = -1
            for d in dirs:
                if path.isfile(pjoin(d, "state")):
                    assert stepnum < 0  # Should never find more than one read state!
                    stepnum = InputFile(pjoin(d, "state"))['FullStepNumber']
            return stepnum + 1

        # Now decide whether WorkingDirectory/read or WorkingDirectory2/read has the current state
        stepnum = get_stepnum([read, pjoin(params['WorkingDirectory2'], "read")])

        # The MT slosh directories are specified via MTDirectory2 (which might come from $ABACUS_SSD2)
        taylors = params['TaylorDirectory']
        taylors2 = params.get('TaylorDirectory2', taylors)
        multipoles = params['MultipoleDirectory']
        multipoles2 = params.get('MultipoleDirectory2', multipoles)

        if stepnum % 2 == 0:
            write = pjoin(params['WorkingDirectory2'], "write")
            read2 = pjoin(params['WorkingDirectory2'], "read")
            multipoles2, multipoles = multipoles, multipoles2
        else:
            read = pjoin(params['WorkingDirectory2'], "read")
            write2 = pjoin(params['WorkingDirectory2'], "write")
            taylors2, taylors = taylors, taylors2

        # Technically we could keep two past states (one in each slosh directory)
        # but in practice this may be too many copies of the state
        past = None

        # Now we need to remove WorkingDirectory from the parameter file,
        # since we're separately specifying the Read and Write dirs
        # Also: when we restart, need to check both possible read directories for the state file
        delparams += ['WorkingDirectory', 'WorkingDirectory2']
        # For now, slosh multipoles instead of splitting
        delparams += ['MultipoleDirectory2', 'TaylorDirectory2']
        
        newparams['ReadStateDirectory'] = read
        newparams['WriteStateDirectory'] = write
        newparams['TaylorDirectory'] = taylors
        newparams['MultipoleDirectory'] = multipoles

        # The convolution needs to know where the multipoles were put *last* step, not this step
        newconvparams = newparams.copy()
        #newconvparams['TaylorDirectory'] = taylors2
        newconvparams['MultipoleDirectory'] = multipoles2

    if OverwriteState:
        write = read
        past = None

    def insert_suffix(base, suffix=newparam_suffix):
        newparamfn = base.split('.')
        suffix = suffix.format(step=stepnum)
        try:
            i = newparamfn.index('par')
            newparamfn.insert(i, suffix.format(step=stepnum))
        except ValueError:
            newparamfn.append(suffix)
        newparamfn = '.'.join(newparamfn)

        return newparamfn

    # Modifications were requested. Make a par file just for this step.
    if delparams or newparams:
        # Make "abacus.stepN.par"
        newparamfn = insert_suffix(paramfn)
        
        params = GenParam.makeInput(newparamfn, paramfn, del_keywords=delparams, **newparams)
    else:
        newparamfn = paramfn

    if newconvparams:
        convparamfn = insert_suffix(paramfn, newparam_suffix + '.conv')
        convparams = GenParam.makeInput(convparamfn, paramfn, del_keywords=delparams, **newconvparams)
    else:
        convparamfn = paramfn

    # Create dirs
    for d in ['MultipoleDirectory', 'MultipoleDirectory2', 'TaylorDirectory', 'TaylorDirectory2']:
        if d in params:
            try:
                os.makedirs(params[d])
            except FileExistsError:
                pass

    return read,write,read2,write2, past, newparamfn, convparamfn


def remove_MT(md, pattern, rmdir=True):
    # Quick util to clear multipoles/taylors
    for mfn in glob(pjoin(md, pattern)):
        os.remove(mfn)
    if rmdir:
        try:
            os.rmdir(md)
        except OSError:
            pass  # remove dir if empty; we re-create at each step in setup_dirs


def singlestep(paramfn, maxsteps, AllowIC=False, monitor=False, stopbefore=-1):
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
    monitor: bool, optional
        Automated monitoring of abacus
    stopbefore: int, optional
        Run abacus until it is about to call "singlestep" for step `stopbefore`.
        Useful for attaching a debugger to a certain step without getting
        interference from Python.

    Returns
    -------
    retcode: int
        Returns 0 or 1 for general success or failure;
        or the ConvolutionDriver or singlestep failure exit value;
        or EXIT_REQUEUE if we finish cleanly but have not reached the final z
    """
    finished = False

    if maxsteps == 0:
        return 0

    if not path.exists(paramfn):
        raise ValueError('Parameter file "{:s}" is not accessible'.format(paramfn))
        
    param = InputFile(paramfn)
    OverwriteState = param.get('OverwriteState',False)

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
    cmdpipe = None
    respipe = None
    if monitor:
        cmdpipe = open(Monitord.inpipefn,"w")
        cmdpipe.write("WATCH:{name}:{bd}:{paramfn}\n".format(name=param.SimName,bd=param.WorkingDirectory,paramfn=path.basename(paramfn)))
        respipe = os.open(Monitord.outpipefn,os.O_RDONLY)
        response = respipe.readline()
        if response !=  "CONFIRM:{name}\n".format(name=param.SimName):
            raise ValueError(response)

    try:
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

        for i in range(maxsteps):
            # Finalize the parameter file for this step, possibly sloshing the state
            read,write,read2,write2,past,paramfn_thisstep,conv_paramfn_thisstep = setup_state_dirs(paramfn)
            param = InputFile(paramfn_thisstep)
            convparam = InputFile(conv_paramfn_thisstep)

            try:
                read_state = InputFile(pjoin(read,"state"))
                ConvDone = check_multipole_taylor_done(convparam, read_state, kind='Taylor')
                stepnum = read_state.FullStepNumber+1
            except IOError:
                assert AllowIC, "If there is no read state, this must be an IC step!"
                ConvDone = True  # No need to convolve for an IC step
                stepnum = 0

            if path.exists(pjoin(write,"state")) and not OverwriteState:
                print()
                answer = input("Write State exists and would be overwritten. Do you want to delete it? (y/[n]) ")
                if answer == "y":
                    shutil.rmtree(write)
                else:
                    raise ValueError
            if not path.exists(write):
                os.makedirs(write)

            # Do the convolution
            if not ConvDone:
                # Now check if we have all the multipoles the convolution will need
                if not check_multipole_taylor_done(convparam, read_state, kind='Multipole'):
                    # Invole multipole recovery mode
                    print("Warning: missing multipoles! Performing multipole recovery for step {:d}".format(i))
                    
                    # Build the make_multipoles executable
                    with Tools.chdir(pjoin(abacuspath, "singlestep")):
                        retcode = subprocess.call(['make', 'make_multipoles'])
                        if retcode:
                            print('make make_multipoles returned error code {}'.format(retcode))
                            return retcode
                    
                    # Check that the multipole directory exists
                    try:
                        os.makedirs(convparam['MultipoleDirectory'])
                    except FileExistsError:
                        pass

                    # Execute it
                    retcode = subprocess.call([pjoin(abacuspath, "singlestep", "make_multipoles"),conv_paramfn_thisstep], env=singlestep_env)
                    if retcode:
                        print('make_multipoles returned error code {}'.format(retcode))
                        return retcode
                    save_log_files(param.LogDirectory, 'step{:04d}.make_multipoles'.format(read_state.FullStepNumber))
                    print('\tFinished multipole recovery for read state {}.'.format(read_state.FullStepNumber))
                
                print("Performing convolution for step {:d}".format(stepnum))
                retcode = subprocess.call([pjoin(abacuspath, "Convolution", "ConvolutionDriver"),conv_paramfn_thisstep], env=convolution_env)
                if retcode:
                    print('Convolution returned error code {}'.format(retcode))
                    return retcode
                    
                if not check_multipole_taylor_done(convparam, read_state, kind='Taylor'):
                    raise ValueError("Convolution did not complete")
                
                convlogs = glob(pjoin(param.LogDirectory, 'last.conv*'))
                for cl in convlogs:
                    shutil.move(cl, cl.replace('last', 'step{:04d}'.format(read_state.FullStepNumber+1)))
    
                # Warning: Convolution won't work if MultipoleDirectory is the write (or read) state
                # because the states get moved after multipole generation but before convolution.
                # This is contrary to the Taylors behavior and may be worth re-thinking.
                # If we assume singlestep is the more "interrupt-prone" step, then this may be acceptable.
                remove_MT(convparam.MultipoleDirectory, 'Multipole*', rmdir=convparam.MultipoleDirectory != param.MultipoleDirectory)
                if 'MultipoleDirectory2' in convparam:
                    remove_MT(convparam.MultipoleDirectory2, 'Multipole*', rmdir=convparam.MultipoleDirectory2 != param.MultipoleDirectory2)

                print("\t Finished convolution for state %d"%(read_state.FullStepNumber+1))

            #run singlestep
            print("Running singlestep for step {:d}".format(stepnum))
            if i == stopbefore:
                print('stopbefore = %d was specified; stopping before calling singlestep'%i)
                return 0
            retcode = subprocess.call([abacuspath+"/singlestep/singlestep",paramfn_thisstep,str(int(AllowIC))], env=singlestep_env)
            
            # In profiling mode, we don't move the states so we can immediately run the same step again
            if param.get('ProfilingMode', False):
                break
                
            # Check that write/state was written as a test of success
            if not path.exists(pjoin(write, "state")):
                raise ValueError("Singlestep did not complete!")
            if retcode:
                return retcode
            write_state = InputFile(pjoin(write, "state"))

            # Update the status log
            with open(status_log_fn, 'a') as status_log:
                status_log.write('{step:4d}  {z:6.3f}\n'.format(step=stepnum, z=write_state.Redshift))

            # save the log and timing files under this step number
            save_log_files(param.LogDirectory, 'step{:04d}'.format(write_state.FullStepNumber))
                            
            shutil.copy(write+"/state", param.LogDirectory+"/step%.4d.state"%(write_state.FullStepNumber))
            print(( "\t Finished state {:d}. a = {:.4f}, dlna = {:.3g}, rms velocity = {:.3g}".format(
                write_state.FullStepNumber, write_state.ScaleFactor,
                write_state.DeltaScaleFactor/(write_state.ScaleFactor-write_state.DeltaScaleFactor),
                write_state.RMS_Velocity)))
            if monitor: cmdpipe.write("BACKUPOUTPUT:{}\n".format(param.LogDirectory))

            # Delete the Taylors right away; want to avoid accidentally
            # running with inconsistent positions and Taylors
            remove_MT(param.TaylorDirectory, 'Taylor_*')
            if 'TaylorDirectory2' in convparam:
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
                shutil.copy(dfn, param.LogDirectory+"/step%.4d.density"%(write_state.FullStepNumber))
                np.savez(param.LogDirectory+"/step%.4d.pow"%(write_state.FullStepNumber),k=k,P=P,nb=nb)
                os.remove(dfn)
            
            # move read to past, write to read
            if past is not None:
                try:
                    shutil.rmtree(past)
                except OSError:
                    pass  # past doesn't exist
                try:
                    os.rename(read, past)
                except OSError:
                    pass  # read doesn't exist, must be IC state
            try:
                # read2 and write2 are the dirs for the next step
                os.rename(write, read2)
                os.makedirs(write2)
            except OSError:
                pass  # read still exists (want to overwrite)
            
            if read != read2:
                try:
                    shutil.rmtree(read)
                except FileNotFoundError:
                    pass  # When sloshing the state, we don't move read to past, so discard it instead

            # If we created a special parameter file just for this state, delete it
            if paramfn != paramfn_thisstep:
                os.remove(paramfn_thisstep)
            if paramfn != conv_paramfn_thisstep:
                os.remove(conv_paramfn_thisstep)

            
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
                if monitor:
                    cmdpipe.write("FINISH:"+param.SimName)
                print("Final redshift reached; terminating normally.")
                finished = True
                break
            
            AllowIC = False

            if write_state.FullStepNumber%20 ==0 and monitor and path.exists(past):
                shutil.move(past,"backup_wip")
                cmdpipe.write("BACKUPSTATE:{}/backup_wip".format(param.WorkingDirectory))
    except:
        if monitor:
            print(sys.exc_info()[0])
            cmdpipe.write("ERROR: \n")
        raise

    # If there is more work to be done, signal that we are ready for requeue
    if retcode == 0 and not finished:
        return EXIT_REQUEUE

    return retcode

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

