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
import pathlib
import re
import fnmatch
import shlex

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

    Files in the `config_dir` directory will be carried around in the `info`
    directory alongside time slices, data products, etc.
    
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

    params = preprocess_params(output_parfile, parfn, use_site_overrides=use_site_overrides,
                override_directories=override_directories, **param_kwargs)

    parallel = params.get('Parallel', False)

    # These directories are global; i.e. treated the same in the parallel and serial versions
    icdir = params['InitialConditionsDirectory']
    outdir = params['OutputDirectory']
    logdir = params['LogDirectory']
    basedir = params['WorkingDirectory']
    groupdir = params.get('GroupDirectory', '')

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

    # If not resuming, erase old global dirs
    if clean:
        clean_dir(basedir, preserve=icdir if not erase_ic else None)
        clean_dir(outdir, preserve=icdir if not erase_ic else None)
        clean_dir(logdir, preserve=icdir if not erase_ic else None)
        clean_dir(groupdir, preserve=icdir if not erase_ic else None)
            
    os.makedirs(basedir, exist_ok=True)

    for d in ['LogDirectory', 'OutputDirectory', 'GroupDirectory']:
        if d in params and params[d]:
            os.makedirs(params[d], exist_ok=True)
    
    try:
        shutil.copy(output_parfile, basedir)
    except shutil.Error:
        pass  # same file?

    info_out_path = pjoin(basedir, 'info')
    copy_contents(config_dir, info_out_path, clean=True)

    # Copy over all of the input files to the output folder for reference
    os.makedirs(outdir, exist_ok=True)
    try:
        shutil.copy(output_parfile, outdir)
    except shutil.Error:
        pass  # same file?
    
    info_out_path = pjoin(outdir, 'info')
    copy_contents(config_dir, info_out_path, clean=True)
        
    # Note: we are cd-ing into the global working directory if this is a parallel run
    with Tools.chdir(basedir):
        # The parfile is directly stored in the basedir
        output_parfile = basename(output_parfile)
        if clean and not path.exists(icdir) and not zeldovich.is_on_the_fly_format(params['ICFormat']):
            zeldovich.run(output_parfile, allow_eigmodes_fn_override=override_directories)
        elif clean:
            print('Reusing existing ICs')

        retval = singlestep(output_parfile, maxsteps, make_ic=clean)

    return retval


def copy_contents(in_dir, out_dir, clean=True, ignore='*.py'):
    if clean:
        if path.isdir(out_dir):
            shutil.rmtree(out_dir)
    os.makedirs(out_dir, exist_ok=True)

    if type(ignore) == str:
        ignore = [ignore]

    for fn in os.scandir(in_dir):
        if not fn.is_file():
            continue
        if any(fnmatch.fnmatch(fn, ipat) for ipat in ignore):
            continue
        shutil.copy(fn, out_dir)


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

def MakeDerivatives(param, derivs_archive_dir=True, floatprec=False):
    '''
    Look for the derivatives required for param and make them if they don't exist.
    '''

    # We will attempt to copy derivatives from the archive dir if they aren't found
    if derivs_archive_dir == True:
        derivs_archive_dir = pjoin(os.getenv('ABACUS_PERSIST',None), 'Derivatives')
    
    os.makedirs(param.DerivativesDirectory, exist_ok=True)

    #note! added this in to run Ewald test for far_radius 1-8,16 
    for ADnum in list(range(1,9)) + [16]:
        ADfn = 'AD32_{:03d}.dat'.format(ADnum)
        source_ADfn = pjoin(abacuspath, "Derivatives", ADfn)
        if not path.isfile(pjoin(param.DerivativesDirectory, ADfn)):
            shutil.copy(source_ADfn,param.DerivativesDirectory)
    
    suffix = "{0.CPD:d}_{0.Order:d}_{0.NearFieldRadius:d}_{0.DerivativeExpansionRadius:d}".format(param) + "_{slab:d}"
    fnfmt32 = "fourierspace_float32_" + suffix
    fnfmt64 = "fourierspace_" + suffix
            
    derivativenames32, derivativenames64 = [], []
    for i in range(0,param.CPD//2+1):
        derivativenames32 += [fnfmt32.format(slab=i)]
        derivativenames64 += [fnfmt64.format(slab=i)]
    derivativenames = derivativenames32 if floatprec else derivativenames64
    fnfmt = fnfmt32 if floatprec else fnfmt64

    create_derivs_cmd = [pjoin(abacuspath, "Derivatives", "CreateDerivatives"),
                        str(param.CPD),str(param.Order),str(param.NearFieldRadius),str(param.DerivativeExpansionRadius)]

    parallel = param.get('Parallel', False)

    if parallel:
        # TODO: parallel CreateDerivs?
        if 'Conv_mpirun_cmd' in param:
            create_derivs_cmd = shlex.split(param['Conv_mpirun_cmd']) + create_derivs_cmd
        else:
            create_derivs_cmd = shlex.split(param['mpirun_cmd']) + create_derivs_cmd

    # First check if the derivatives are in the DerivativesDirectory
    if not all(path.isfile(pjoin(param.DerivativesDirectory, dn)) for dn in derivativenames):
        
        # If not, check if they are in the canonical $ABACUS_PERSIST/Derivatives directory
        if all(path.isfile(pjoin(derivs_archive_dir, dn)) for dn in derivativenames):
            print('Found derivatives in archive dir "{}". Copying to DerivativesDirectory "{}".'.format(derivs_archive_dir, param.DerivativesDirectory))
            for dn in derivativenames:
                shutil.copy(pjoin(derivs_archive_dir, dn), pjoin(param.DerivativesDirectory, dn))
        else:
            print('Could not find derivatives in "{}" or archive dir "{}". Creating them...'.format(param.DerivativesDirectory, derivs_archive_dir))
            print("Error was on file pattern '{}'".format(fnfmt))

            if floatprec:
                # first make sure the derivatives exist in double format
                MakeDerivatives(param, floatprec=False)
                # now make a 32-bit copy
                from . import convert_derivs_to_float32
                for dpath in derivativenames64:
                    convert_derivs_to_float32.convert(pjoin(param.DerivativesDirectory, dpath))
            else:
                with Tools.chdir(param.DerivativesDirectory):
                    try:
                        os.remove("./farderivatives")
                    except OSError:
                        pass
                    if parallel:
                        print('Dispatching CreateDerivatives with command "{}"'.format(' '.join(create_derivs_cmd)))
                    subprocess.check_call(create_derivs_cmd)
    
    
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
    parser = argparse.ArgumentParser(description='Run this sim.', formatter_class=Tools.ArgParseFormatter)
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

    1) Compute ZD_Pk_sigma from sigma_8
    2) If $ABACUS_SSD2 is defined, define the params {Multipole,Taylor}Directory2
    3) If $ABACUS_TMP2 is defined and StateIOMode is "slosh", define WorkingDirectory2
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

    if params.get('StateIOMode','normal').lower() in ('slosh','stripe'):
        # TODO: state striping not implemented
        if 'ABACUS_TMP2' in os.environ:
            # WorkingDirectory2 is never processed by Abacus if we are sloshing
            # It will be intercepted by the sloshing logic
            if 'WorkingDirectory2' not in params:
                param_kwargs['WorkingDirectory2'] = pjoin(os.environ['ABACUS_TMP2'], params['SimName'])

    Conv_IOMode = params.get('Conv_IOMode','normal').lower()
    if 'ABACUS_SSD2' in os.environ:
        if Conv_IOMode == 'slosh':
            # Don't bother failing with assertions here.  We'll check everything when we actually try to set up the state dirs
            # Similar to WorkingDirectory2, this dir is never used by Abacus proper
            if 'ConvolutionWorkingDir2' not in params:
                param_kwargs['ConvolutionWorkingDir2'] = pjoin(os.environ['ABACUS_SSD2'], params['SimName'])
        elif Conv_IOMode == 'stripe':
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
    the symlinks if state sloshing is requested via `StateIOMode`.

    When state sloshing, everything looks normal to the singlestep
    binary, except "read" and "write" are actually symlinks
    to the "state1" and "state2" directories elsewhere.

    Sloshing/striping of the multipole/taylors can be requested
    with the `Conv_IOMode` parameter.

    When doing a parallel run, we can set up global directories here,
    but not node-local ones.
    '''
    params = GenParam.makeInput(paramfn, paramfn)
    StateIOMode = params.get('StateIOMode','normal').lower()
    Conv_IOMode = params.get('Conv_IOMode','normal').lower()

    parallel = params.get('Parallel', False)

    # Presently, we only support normal operation and overwriting in the parallel version
    if parallel:
        assert StateIOMode in ('normal', 'overwrite')
        assert Conv_IOMode in ('normal', 'overwrite')
    
    # Normal operation: everything is under the WorkingDirectory or specified individually
    try:
        read = params['ReadStateDirectory']
    except KeyError:
        read = pjoin(params['WorkingDirectory'], "read")

    try:
        write = params['WriteStateDirectory']
    except KeyError:
        write = pjoin(params['WorkingDirectory'], "write")

    try:
        past = params['PastStateDirectory']
    except KeyError:
        past = pjoin(params['WorkingDirectory'], "past")

    past = normpath(past)
    read = normpath(read)
    write = normpath(write)

    if not parallel:
        taylors = normpath(params['TaylorDirectory'])
        multipoles = normpath(params['MultipoleDirectory'])
    else:
        # In parallel version, M/T are purely node-local
        taylors = None
        multipoles = None

    if StateIOMode == 'overwrite':
        write = read
        past = None

    # Check if the write state already exists
    if path.isfile(pjoin(write, "state")) and StateIOMode != 'overwrite':
        answer = input('\nWrite state "{}" exists and would be overwritten. Do you want to delete it? (y/[n]) '.format(write))
        if answer == "y":
            shutil.rmtree(write)
        else:
            raise ValueError("Cannot proceed if write state exists!")

    # If the links don't exist, create them and the underlying dirs
    # If they do exist, don't touch them
    def make_link(link, target):
        if path.exists(link):
            # link already existing as a file/directory instead of a link is an error!
            assert path.islink(link)
        else:
            os.makedirs(target, exist_ok=True)
            os.symlink(target, link)
    
    # Set up symlinks to slosh the state
    # This approach ensures that Abacus proper doesn't need to know anything about sloshing
    if StateIOMode == 'slosh':
        # Should be populated manually or from ABACUS_TMP2
        assert 'WorkingDirectory' in params  # need a working directory to put the symlinks
        assert 'WorkingDirectory2' in params, \
            "If StateIOMode = slosh, must specify WorkingDirectory2 via $ABACUS_TMP2 or directly in the parameter file!"

        make_link(read, pjoin(params['WorkingDirectory'], 'state1'))
        make_link(write, pjoin(params['WorkingDirectory2'], 'state2'))

        # Normally we move 'read' to 'past' after singlestep, but that might cross devices if sloshing!
        # Technically we could keep two past states (one in each slosh directory)
        # but in practice this may be too many copies of the state
        past = None
    else:
        os.makedirs(read, exist_ok=True)
        os.makedirs(write, exist_ok=True)
        if past:
            os.makedirs(past, exist_ok=True)

    # Need M/TDirectory2 for either sloshing or striping
    if Conv_IOMode == 'slosh':
        assert 'ConvolutionWorkingDir2' in params, \
            "If Conv_IOMode = slosh, must specify ConvolutionWorkingDir2 via $ABACUS_SSD2 or directly in the parameter file!"

        # we use dirname(taylors) as ConvolutionWorkingDir
        assert dirname(taylors) == dirname(multipoles)

        make_link(taylors, pjoin(dirname(taylors), 'convstate1'))
        make_link(multipoles, pjoin(params['ConvolutionWorkingDir2'], 'convstate2'))

    if Conv_IOMode == 'stripe':
        assert 'MultipoleDirectory2' in params and 'TaylorDirectory2' in params, \
            "If Conv_IOMode = stripe, must specify MultipoleDirectory2 & TaylorDirectory2 via $ABACUS_SSD2 or directly in the parameter file!"

    # Create dirs
    if not parallel:
        for d in ['MultipoleDirectory', 'MultipoleDirectory2', 'TaylorDirectory', 'TaylorDirectory2']:
            if d in params:
                os.makedirs(params[d], exist_ok=True)

    return read,write,past,multipoles,taylors


def move_state_dirs(read, write, past):
    '''
    At the end of every timestep, we to move the "write" state
    to the read directory, move "read" to "past", and discard
    "past".

    If we're sloshing the state, we instead want to swap the
    symlinks for `read` and `write`.  We also want to swap
    the symlink for the Taylors (but not the multipoles because
    the convolve hasn't happened yet).

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

# via https://stackoverflow.com/a/2293793
fp_regex = r'(([1-9][0-9]*\.?[0-9]*)|(\.[0-9]+))([Ee][+-]?[0-9]+)?'

import table_logger
class StatusLogWriter:
    #header = "Step  Redshift  Singlestep  Conv"
    #line_fmt = '{step:4d}  {z:6.1f}  {ss_rate:.1g} Mp/s ({ss_time:.1g} s) {conv_time:.1g}'

    fields = {'Step': '{:4d}',
              'Redshift': '{:.3g}',
              'Singlestep': '{0[0]:.3g} Mp/s ({0[1]:.3g} s)',
              'Conv': '{:.3g} s'}

    topmatter = ['Abacus Status Log',
                 'simname, timestamp',
                 '==================\n\n',]

    def __init__(self, log_fn):
        self.fields = {k:v.format for k,v in self.fields.items()}
        self.log_fp = open(log_fn, 'ab')

        self.log_fp.write('\n'.join(self.topmatter).format().encode('utf-8'))
        self.log_fp.flush()

        self.logger = table_logger.TableLogger(file=self.log_fp,
                                               columns=list(self.fields),
                                               formatters=self.fields)
    def __del__(self):
        self.logger.make_horizontal_border()
        self.log_fp.close()

    def update(self, param, state, ss_time=None, conv_time=None):
        '''
        Update the log with some info about the singlestep and conv
        that just finished.

        The logs have already been moved to their new names, indicated
        by step_num.
        '''

        step_num = state.FullStepNumber

        if not ss_time:
            self.print('Warning: parsing logs to get singlestep time, may miss startup time')
            ss_log_fn = pjoin(param['LogDirectory'], 'step{:04d}.time'.format(step_num))
            try:
                ss_log_txt = pathlib.Path(ss_log_fn).read_text()
            except FileNotFoundError:
                self.print('Warning: could not find timing file "{}"'.format(ss_log_fn))
                ss_time = np.nan

            matches = re.search(r'Total Wall Clock Time\s*:\s*(?P<time>{fp:s})'.format(fp=fp_regex), ss_log_txt)
            ss_time = float(matches.group('time'))
        ss_rate = param['NP']/1e6/ss_time  # Mpart/s

        if not conv_time:
            conv_log_fn = pjoin(param['LogDirectory'], 'step{:04d}.convtime'.format(step_num))
            try:
                conv_log_txt = pathlib.Path(conv_log_fn).read_text()
                matches = re.search(r'ConvolutionWallClock\s*:\s*(?P<time>{fp:s})'.format(fp=fp_regex), conv_log_txt)
                conv_time = float(matches.group('time'))
                self.print('Warning: parsing logs to get convolution time, may miss startup time')
            except:
                conv_time = 0.

        info = dict(Step=step_num, Redshift=state.Redshift, Singlestep=(ss_rate,ss_time), Conv=conv_time)

        self.logger(*(info[k] for k in self.fields))


    def print(self, fmtstring, end='\n', *args, **kwargs):
        '''
        Print a plain statement to the status log
        '''
        self.log_fp.write(('\n * ' + fmtstring.format(*args, **kwargs) + end).encode('utf-8'))
        self.log_fp.flush()
    


def singlestep(paramfn, maxsteps=None, make_ic=False, stopbefore=-1):
    """
    Run a number of Abacus timesteps by invoking the `singlestep` and
    `ConvolutionDriver` executables in a tick-tock fashion.

    If running a parallel job, instead calls `mpirun singlestep` and
    `mpirun ConvolutionDriver` (or whatever MPI launch command has
    been configured).

    Parameters
    ----------
    paramfn: str
        File name of the parameter file for this simulation
    maxsteps: int
        Run (at most) `maxsteps` timesteps.
    make_ic: bool, optional
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

    if maxsteps is None:
        maxsteps = 10000

    maxsteps = int(maxsteps)

    if not path.exists(paramfn):
        raise ValueError('Parameter file "{:s}" is not accessible'.format(paramfn))
        
    param = InputFile(paramfn)

    parallel = param.get('Parallel', False)

    do_fake_convolve = param.get('FakeConvolve', False)

    # Set up backups
    try:
        backups_enabled = param['BackupStepInterval'] > 0
    except:
        backups_enabled = False  # No BackupIntervalSteps parameter
    
    if parallel:
        # TODO: figure out how to signal a backup to the nodes
        backups_enabled = False

    #if not backups_enabled:
    #    print 'Warning: automatic state backup not enabled.  Set "BackupStepInterval" and "BackupDirectory" in the parameters file.'

    ProfilingMode = param.get('ProfilingMode',False)
        
    # Set up OpenMP singlestep environment
    # TODO: check that calling mpirun with this environment is equivalent to calling singlestep with the environment
    singlestep_env = setup_singlestep_env(param)
    convolution_env = setup_convolution_env(param)

    status_log = StatusLogWriter(pjoin(param.OutputDirectory, 'status.log'))
    conv_time = None

    print("Using parameter file {:s} and working directory {:s}.".format(paramfn,param.WorkingDirectory))

    # floatprec=True effectively guarantees both precisions are available
    # TODO: how to tell if Abacus was compiled in double precision?
    if not do_fake_convolve:
        MakeDerivatives(param, floatprec=True)

    if make_ic:
        print("Ingesting IC from "+param.InitialConditionsDirectory+" ... Skipping convolution")
        if path.exists(pjoin(param.InitialConditionsDirectory, "input.pow")):
            shutil.copy(pjoin(param.InitialConditionsDirectory, "input.pow"), pjoin(param.LogDirectory, "input.pow"))

    # Finalize the directories, possibly setting up symlinks for state sloshing
    # Careful: multipoles & taylors can be None if doing a parallel run!
    read,write,past,multipoles,taylors = setup_state_dirs(paramfn)

    if parallel:
        try:
            mpirun_cmd = shlex.split(param['mpirun_cmd'])
            if 'Conv_mpirun_cmd' in param:
                Conv_mpirun_cmd = shlex.split(param['Conv_mpirun_cmd'])
            else:
                Conv_mpirun_cmd = mpirun_cmd
        except KeyError:
            raise RuntimeError('"mpirun_cmd" must be specified in the parameter file if "Parallel = 1"!')

        # Dole out the state to the node if requested (e.g. resuming from a saved state on the network file system)
        distribute_state_from = param.get('DistributeStateFrom', None)
        if distribute_state_from:
            print('Distributing state to nodes...')
            distribute_state_cmd = [pjoin(abacuspath, 'Abacus', 'distribute_state_to_nodes.py'), paramfn, distribute_state_from]
            subprocess.check_call(Conv_mpirun_cmd + distribute_state_cmd)

    print("Beginning abacus steps:")

    for i in range(maxsteps):
        if make_ic:
            ConvDone = True  # No need to convolve for an IC step
            stepnum = 0
        else:
            try:
                read_state = InputFile(pjoin(read,"state"))
                if parallel:
                    ConvDone = False  # Have to take it on faith that the nodes have multipoles and are ready to convolve!
                else:
                    ConvDone = check_multipole_taylor_done(param, read_state, kind='Taylor')
                stepnum = read_state.FullStepNumber + 1
            except:
                # TODO: might have no read state if we're restoring from a backup
                ConvDone = False
                stepnum = 0

        ss_timer, conv_time = None, None

        # Do the convolution
        if not ConvDone:
            if not parallel and not do_fake_convolve:
                # Now check if we have all the multipoles the convolution will need
                if not check_multipole_taylor_done(param, read_state, kind='Multipole'):
                    # Invoke multipole recovery mode
                    print("Warning: missing multipoles! Performing multipole recovery for step {:d}".format(i))
                    
                    # Build the recover_multipoles executable
                    with Tools.chdir(pjoin(abacuspath, "singlestep")):
                        subprocess.check_call(['make', 'recover_multipoles'])

                    # Execute it
                    print("Running recover_multipoles for step {:d}".format(stepnum))
                    subprocess.check_call([pjoin(abacuspath, "singlestep", "recover_multipoles"), paramfn], env=singlestep_env)
                    save_log_files(param.LogDirectory, 'step{:04d}.recover_multipoles'.format(read_state.FullStepNumber))
                    print('\tFinished multipole recovery for read state {}.'.format(read_state.FullStepNumber))

                # Swap the Taylors link.  In effect, this will place the Taylors on the same disk as the multipoles.
                # But that's what we want for sloshing: this was the write disk, so now it will be the upcoming read disk
                # We'll swap the multipoles as soon as we convolve
                if path.islink(taylors):
                    taylors_convstate = os.readlink(taylors)
                    os.unlink(taylors)
                    os.symlink(os.readlink(multipoles), taylors)
            
            if do_fake_convolve:
                ConvolutionDriver_cmd = [pjoin(abacuspath, "Convolution", "FakeConvolution.py"), paramfn]
            else:
                ConvolutionDriver_cmd = [pjoin(abacuspath, "Convolution", "ConvolutionDriver"), paramfn]
            if parallel:
                ConvolutionDriver_cmd = Conv_mpirun_cmd + ConvolutionDriver_cmd
                print('Performing parallel convolution for step {:d} with command "{:s}"'.format(stepnum, ' '.join(ConvolutionDriver_cmd)))
            else:
                print("Performing convolution for step {:d}".format(stepnum))
            with Tools.ContextTimer() as conv_timer:
                subprocess.check_call(ConvolutionDriver_cmd, env=convolution_env)
            conv_time = conv_timer.elapsed

            if ProfilingMode == 2:
                print('\tConvolution step {} finished. ProfilingMode = 2 (Convolution) is active; Abacus will now quit.'.format(stepnum))
                break

            if not parallel and not check_multipole_taylor_done(param, read_state, kind='Taylor'):
                # have to take it on faith in the parallel version!
                raise ValueError("Convolution did not complete")
            
            convlogs = glob(pjoin(param.LogDirectory, 'last.*conv*'))
            for cl in convlogs:
                shutil.move(cl, cl.replace('last', 'step{:04d}'.format(read_state.FullStepNumber+1)))

            if not parallel:
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
        if i == stopbefore:
            print('stopbefore = %d was specified; stopping before calling singlestep'%i)
            return 0
        
        singlestep_cmd = [pjoin(abacuspath, "singlestep", "singlestep"), paramfn, str(int(make_ic))]
        if parallel:
            singlestep_cmd = mpirun_cmd + singlestep_cmd
            print('Running parallel singlestep for step {:d} with command "{:s}"'.format(stepnum, ' '.join(singlestep_cmd)))
        else:
            print("Running singlestep for step {:d}".format(stepnum))
        with Tools.ContextTimer() as ss_timer:
            subprocess.check_call(singlestep_cmd, env=singlestep_env)
        
        # In profiling mode, we don't move the states so we can immediately run the same step again
        if ProfilingMode and ProfilingMode != 2:
            print('\tStep {} finished. ProfilingMode is active; Abacus will now quit and states will not be moved.'.format(stepnum))
            break
            
        # Check that write/state was written as a test of success
        if not path.exists(pjoin(write, "state")):
            raise ValueError("Singlestep did not complete!")
        write_state = InputFile(pjoin(write, "state"))

        # save the log and timing files under this step number
        save_log_files(param.LogDirectory, 'step{:04d}'.format(write_state.FullStepNumber))

        # Update the status log
        status_log.update(param, write_state, ss_timer.elapsed, conv_time)
        
        shutil.copy(pjoin(write, "state"), pjoin(param.LogDirectory, "step{:04d}.state".format(write_state.FullStepNumber)))
        print(( "\t Finished state {:d}. a = {:.4f}, dlna = {:.3g}, rms velocity = {:.3g}".format(
            write_state.FullStepNumber, write_state.ScaleFactor,
            write_state.DeltaScaleFactor/(write_state.ScaleFactor-write_state.DeltaScaleFactor),
            write_state.RMS_Velocity)))

        if not parallel:
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
        
        if not parallel:
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
        move_state_dirs(read, write, past)
        
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
        if not make_ic and np.abs(read_state.Redshift - finalz) < 1e-12 and read_state.LPTStepNumber == 0:
            print("Final redshift of {:g} reached; terminating normally.".format(finalz))
            status_log.print("Final redshift of {:g} reached; terminating normally.".format(finalz))
            finished = True
            break
        
        make_ic = False

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
