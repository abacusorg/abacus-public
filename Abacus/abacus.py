#!/usr/bin/env python3
# Copyright 2012-2025 The Abacus Developers
# SPDX-License-Identifier: GPL-3.0-or-later

'''
This module is the main entry point for running an Abacus simulation.
For a complete example of a simulation setup that demonstrates how
to use this module, see the `abacus/Production/Example` directory.

This module is mostly a wrapper for the `singlestep` binary. It
handles moving the state folders around, etc.  It can be imported
as a module or directly invoked as a script:

Script usage:
    See abacus/Production/Example

'''

import os
import sys
abacuspath = os.getenv("ABACUS")
if not abacuspath:
    print("Error: Please define $ABACUS to be the absolute path to your abacus distribution", file=sys.stderr)
    sys.exit(1)

import os.path as path
from os.path import join as pjoin, normpath, basename, dirname
import subprocess
import shutil
import ctypes as ct
import argparse
from glob import glob
import pathlib
import re
import fnmatch
import shlex
import time
import signal
from tempfile import mktemp, mkstemp
import warnings
from warnings import warn
import datetime

import numpy as np

from .InputFile import InputFile
from . import Tools
from . import GenParam
from . import zeldovich
from . import derivatives
from .Tools import call_subprocess


NEEDS_INTERIM_BACKUP_MINS = 105  #minimum job runtime (in minutes) for which we'd like to do a backup halfway through the job. 
EXIT_REQUEUE = 200
GF_BACKUP_INTERVAL = 2 * 60 * 60  #only backup b/w group finding steps if it's been more than two hours since the last backup.
BACKUP_EVERY_N_MINUTES = 10000


site_param_fn = pjoin(abacuspath, 'Production', 'site_files', 'site.def')
directory_param_fn = pjoin(abacuspath, 'Production', 'directory.def')
wall_timer = time.perf_counter  # monotonic wall clock time

def timestamp():
    return str(datetime.datetime.now().replace(microsecond=0))


def run(parfn='abacus.par2', config_dir=path.curdir, maxsteps=10000, clean=False, erase_ic=False, output_parfile=None, use_site_overrides=False, override_directories=False, just_params=False, param_kwargs=None):
    """
    The main entry point for running a simulation.  Initializes directories,
    copies parameter files, and parses parameter files as appropriate, then
    invokes `singlestep()`.

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
    if not output_parfile:
        # TODO: really we want a way to load the params to memory
        output_parfile = mktemp(prefix='abacus', suffix='.par')

    # TODO: summit doesn't let us write files from the compute nodes.  Is that a problem?
    params = preprocess_params(output_parfile, parfn, use_site_overrides=use_site_overrides,
                override_directories=override_directories, param_kwargs=param_kwargs)
    if just_params:
        return 0

    parallel = params.get('Parallel', False)

    # These directories are global; i.e. treated the same in the parallel and serial versions
    icdir = params['InitialConditionsDirectory']
    outdir = params['OutputDirectory']
    lcdir = params['LCDirectory']
    logdir = params['LogDirectory']
    groupdir  = params.get('GroupDirectory', '')
    basedir = params['WorkingDirectory']
    backupdir = params.get('BackupDirectory', '')

    # might be global; might be local
    multipoledir = params.get('MultipoleDirectory', '')
    taylordir = params.get('TaylorDirectory', '')

    resumedir = ""
    if parallel:
        resumedir = pjoin(dirname(params['WorkingDirectory']), params['SimName'] + '_retrieved_state')
    else:
        resumedir = None

    # If we requested a resume, but there is no state, assume we are starting fresh
    if not clean:
        wd = params.get('WorkingDirectory', basedir)
        if not parallel and not path.exists(wd) and path.exists(pjoin(backupdir, 'backup')):
            print('Restoring state from backup...', end='', flush=True)
            with Tools.ContextTimer() as restore_timer:
                shutil.copytree(pjoin(backupdir, 'backup'), wd, symlinks=True)
            print(f' done in {restore_timer.elapsed:.1f} s', flush=True)

        # TODO: in an interactive parallel session, the state may already exist on the remote nodes. How to check?
        if (parallel and path.exists(resumedir)) or (not parallel and path.exists(basedir)):
            print('Resuming from existing state.')
        else:
            print('Resume requested but no state exists.  Creating one.')
            clean = True

    if erase_ic and path.exists(icdir):
        print('Erasing IC dir')
        shutil.rmtree(icdir)

    # If not resuming, erase old global dirs
    if clean:
        clean_dirs( ( basedir,
                      outdir,
                      logdir,
                      groupdir,
                      lcdir,
                      multipoledir,
                      taylordir ),
                    preserve=icdir if not erase_ic else None
                  )
        #NAM make prettier.
        if parallel and path.exists(resumedir):
            clean_dirs([resumedir],
                preserve=icdir if not erase_ic else None)

    os.makedirs(basedir, exist_ok=True)

    for d in ['LogDirectory', 'OutputDirectory', 'GroupDirectory', 'LCDirectory']:
        if d in params and params[d]:
            os.makedirs(params[d], exist_ok=True)

    shutil.copy(output_parfile, pjoin(basedir, 'abacus.par'))
    output_parfile = pjoin(basedir, 'abacus.par')

    # Copy over all of the input files to the output folder for reference
    os.makedirs(outdir, exist_ok=True)
    try:
        shutil.copy(output_parfile, outdir)
    except shutil.Error:
        pass  # same file?

    # Note: we are cd-ing into the global working directory if this is a parallel run
    with Tools.chdir(basedir):
        # The parfile is directly stored in the basedir
        output_parfile = basename(output_parfile)
        if not zeldovich.is_on_the_fly_format(params['ICFormat']) and not params.get('ExternalICs', False):
            if clean and not zeldovich.have_ics(params):
                zeldovich.run(output_parfile, allow_eigmodes_fn_override=override_directories)
            if params.get('Parallel') and params.get('NumZRanks',1) > 1:
                zeldovich.ensure_2d(params)  # might have switched from 1D to 2D
        elif clean:
            print('Reusing existing ICs')

        retval = singlestep(output_parfile, maxsteps, make_ic=clean, resume_dir = resumedir)

    #if parallel:
    #    handle_requeue(retval, parfn)

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


def clean_dirs(dirs, preserve=None, rmdir_ifempty=True):
    '''
    Remove the contents of a list of directories, except for the path
    specified in the "preserve" argument.  Optionally remove
    the directory itself if it's empty.
    '''

    for bd in dirs:
        if bd == None:
            return 0

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
    parser.add_argument('--just-params', help="Only generate the .par file from the .par2 file", action='store_true')
    return parser


def preprocess_params(output_parfile, parfn, use_site_overrides=False, override_directories=False, param_kwargs=None):
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
    # LHG TODO: not using _param_kwargs, probably a bug?
    if param_kwargs is None:
        param_kwargs = {}
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

    if not params.get('ExternalICs', False) and not zeldovich.is_on_the_fly_format(params['ICFormat']):
        zd_params = zeldovich.setup_zeldovich_params(params)
        if zd_params:
            param_kwargs.update(zd_params)
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

    if not params.get('Parallel',False):
        # In the parallel code, we don't know Multipole/TaylorDirectory until
        # inside singlestep (the path uses the rank). So we should leave those
        # dirs unset in Python. But when we switch to the serial code, we want
        # to set the paths here so we can slosh or do other manipulations.
        if not params.get('TaylorDirectory'):
            param_kwargs['TaylorDirectory'] = pjoin(os.environ['ABACUS_SSD'], params['SimName'], 'taylor')
        if not params.get('MultipoleDirectory'):
            param_kwargs['MultipoleDirectory'] = pjoin(os.environ['ABACUS_SSD'], params['SimName'], 'multipole')


    params = GenParam.makeInput(output_parfile, parfn, **param_kwargs)

    return params


def check_multipole_taylor_done(param, state, kind):
    """
    Checks that all the multipoles or Taylors exist and have the expected
    file size for the serial code. `kind` must be 'Multipole' or 'Taylor'.
    """

    assert kind in ['Multipole', 'Taylor']
    rml = (param.Order+1)**2
    sizeof_Complex = 16 if state.DoublePrecision else 8
    size = param.CPD*(param.CPD+1)//2*rml*sizeof_Complex

    even_dir = param[kind+'Directory']
    odd_dir = param.get(kind+'Directory2', param[kind+'Directory'])
    if kind == 'Multipole':
        kind += 's'  # should consider standardizing our pluralization
    try:
        even = all(path.getsize(pjoin(even_dir, kind+f'_{i:04d}')) == size
            for i in range(0,param.CPD,2))
        odd = all(path.getsize(pjoin(odd_dir, kind+f'_{i:04d}')) == size
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
    if 'OMP_NUM_THREADS' in param and 'OMP_NUM_THREADS' in singlestep_env:
        if int(param['OMP_NUM_THREADS']) != int(singlestep_env['OMP_NUM_THREADS']):
            warn('OMP_NUM_THREADS in the parameter file and the environment do not match. '
                'To avoid confusion, they should be the same (or the environment variable should be unset).')

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


# If the links don't exist, create them and the underlying dirs
# If they do exist, don't touch them
def make_link(link, target):
    if path.exists(link):
        # link already existing as a file/directory instead of a link is an error!
        assert path.islink(link)
        assert path.samefile(os.readlink(link), target)
    else:
        os.makedirs(target, exist_ok=True)
        os.symlink(target, link)


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
    params = InputFile(paramfn)
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
        # make write a symlink to read, using relative paths
        make_link(write, basename(read))
        past = None

    # Check if the write state already exists
    if path.isfile(pjoin(write, "state")) and StateIOMode != 'overwrite':
        answer = input(f'\nWrite state "{write}" exists and would be overwritten. Do you want to delete it? (y/[n]) ')
        if answer == "y":
            shutil.rmtree(write)
        else:
            raise ValueError("Cannot proceed if write state exists!")

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
    '''

    # move read to past, write to read
    # let's assume sloshing never uses past; otherwise, we need a past in each working dir
    if past is None:
        if not path.samefile(read, write):
            try:
                shutil.rmtree(read)
            except OSError:
                # can't rmtree symlink! but we at least want to delete the contents
                clean_dirs(read, preserve=None, rmdir_ifempty=False)
    else:
        try:
            shutil.rmtree(past)
        except OSError:
            pass  # past doesn't exist
        try:
            os.rename(read, past)
        except OSError:
            pass  # read doesn't exist, must be IC state

    if path.exists(read) and path.samefile(read,write):
        return  # overwriting; nothing to do

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
              'Redshift': '{:#.4g}',
              'Rate': '{:#.4g} Mp/s',   #'{0[0]:.4g} Mp/s, {0[1]:.4g}  s)',
              'Elapsed': '{:#.4g} s',  #'{0[0]:.4g} Mp/s, {0[1]:.4g}  s)',
              'Conv': '{:#.4g} s',
              'TauHMpc': '{:#6.1f}',
              'DeltaZ': '{:#.3g}',
              'Time': '{:#.4g}',
              'DeltaT': '{:#.3g}',
              'RMSVel': '{:#.3g}',
              'MaxVel': '{:#7.2f}',
              'DirectPP': '{:6.0f}',
              'RMSCell': '{:#.4g}',
              'MaxL0Sz': '{:7d}',
              'GrpDiam': '{:2d}'
              }

    colwidth = {'Step': 4,
              'Redshift': 8,
              'Rate': 12,
              'Elapsed': 10,
              'Conv': 10,
              'TauHMpc': 7,
              'DeltaZ': 8,
              'Time': 8,
              'DeltaT': 8,
              'RMSVel': 7,
              'MaxVel': 7,
              'DirectPP':8,
              'RMSCell': 8,
              'MaxL0Sz': 7,
              'GrpDiam': 7
              }


    topmatter = ['# Abacus Status Log\n']

    def __init__(self, log_fn):
        self.fields = {k:v.format for k,v in self.fields.items()}
        self.log_fp = open(log_fn, 'ab')

        self.log_fp.write('\n'.join(self.topmatter).format().encode('utf-8'))
        self.log_fp.flush()

        self.logger = table_logger.TableLogger(file=self.log_fp,
                                               columns=list(self.fields),
                                               colwidth=self.colwidth,
                                               formatters=self.fields,
                                               border=False)
    def __del__(self):
        self.logger.make_horizontal_border()
        self.log_fp.close()

    def update(self, param, read_state, write_state, ss_time=None, conv_time=None):
        '''
        Update the log with some info about the singlestep and conv
        that just finished.

        The logs have already been moved to their new names, indicated
        by step_num.

        We gather info from both read_state and write_state.
        At Step N, we are drifting from z_read to z_write.
        Timeslice outputs and groups are coming from z_read, so we output that in the redshift column.
        Lightcones are coming from the annulus between z_read and z_write.
        Velocity stats are from the midway point between the two.
        The timing info is regarding this step, so delta(times) are from write_state.
        '''

        step_num = write_state.FullStepNumber

        if not ss_time:
            self.print('Warning: parsing logs to get singlestep time, may miss startup time')
            ss_log_fn = pjoin(param['LogDirectory'], f'step{step_num:04d}', f'step{step_num:04d}.time')
            try:
                ss_log_txt = pathlib.Path(ss_log_fn).read_text()
            except FileNotFoundError:
                self.print(f'Warning: could not find timing file "{ss_log_fn}"')
                ss_time = np.nan

            matches = re.search(rf'Total Wall Clock Time\s*:\s*(?P<time>{fp_regex:s})', ss_log_txt)
            ss_time = float(matches.group('time'))

        if not conv_time:
            conv_log_fn = pjoin(param['LogDirectory'], f'step{step_num:04d}', f'step{step_num:04d}.convtime')
            try:
                conv_log_txt = pathlib.Path(conv_log_fn).read_text()
                matches = re.search(rf'ConvolutionWallClock\s*:\s*(?P<time>{fp_regex:s})', conv_log_txt)
                conv_time = float(matches.group('time'))
                self.print('Warning: parsing logs to get convolution time, may miss startup time')
            except:
                conv_time = 0.
        ss_rate = param['NP']/1e6/(ss_time+conv_time)  # Mpart/s

        code_to_kms = write_state.VelZSpace_to_kms / write_state.VelZSpace_to_Canonical
        info = dict(Step=step_num, Rate=ss_rate, Elapsed=ss_time+conv_time, Conv=conv_time,
            Redshift=read_state.Redshift, DeltaZ=write_state.DeltaRedshift,
            TauHMpc=read_state.CoordinateDistanceHMpc,
            Time=read_state.Time, DeltaT=write_state.DeltaTime,
            GrpDiam=write_state.MaxGroupDiameter, MaxL0Sz=write_state.MaxL0GroupSize,
            RMSVel=write_state.RMS_Velocity*code_to_kms, MaxVel=write_state.MaxVelocity*code_to_kms,
            DirectPP=write_state.DirectsPerParticle,
            RMSCell=write_state.StdDevCellSize )

        self.logger(*(info[k] for k in self.fields))
        self.log_fp.flush()


    def print(self, fmtstring, end='\n', tee=False, *args, **kwargs):
        '''
        Print a plain statement to the status log
        '''
        out = fmtstring.format(*args, **kwargs) + end
        if tee:
            print(out, end='', flush=True)
        self.log_fp.write(out.encode('utf-8'))
        self.log_fp.flush()


def singlestep(paramfn, maxsteps=None, make_ic=False, stopbefore=-1, resume_dir=None):
    """
    Run a number of Abacus timesteps by invoking the `singlestep` and
    `convolution` executables in a tick-tock fashion.

    If running a parallel job, instead calls `mpirun singlestep` and
    `mpirun convolution` (or whatever MPI launch command has
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
        raise ValueError(f'Parameter file "{paramfn:s}" is not accessible')

    param = InputFile(paramfn)

    parallel = param.get('Parallel', False)
    twoD = param.get('NumZRanks', 0) > 1 and parallel

    do_fake_convolve = param.get('FakeConvolve', False)

    # Set up backups
    try:
        backups_enabled = param['BackupStepInterval'] > 0
    except:
        backups_enabled = False  # No BackupIntervalSteps parameter
    
    run_time_secs = int(os.getenv('SLURM_JOB_END_TIME', '600000')) - int(os.getenv('SLURM_JOB_START_TIME', '0'))
    run_time_secs -= max(min(run_time_secs*0.1, 600), 120)
    run_time_minutes = run_time_secs / 60
    

    start_time = wall_timer()
    print(f"Beginning run at {timestamp()}, running for at most {run_time_minutes} minutes.\n")

    if parallel:
        backups_enabled = False

    #if not backups_enabled:
    #    print 'Warning: automatic state backup not enabled.  Set "BackupStepInterval" and "BackupDirectory" in the parameters file.'

    ProfilingMode = param.get('ProfilingMode',False)

    # Set up OpenMP singlestep environment
    singlestep_env = setup_singlestep_env(param)
    convolution_env = setup_convolution_env(param)

    status_log = StatusLogWriter(pjoin(param.OutputDirectory, 'status.log'))
    conv_time = None

    starting_time = time.time()
    starting_time_str = time.asctime(time.localtime())

    status_log.print(f"# Starting {param.SimName:s} at {starting_time_str:s}")

    print(f"Using parameter file {paramfn:s} and working directory {param.WorkingDirectory:s}.")

    # floatprec=True effectively guarantees both precisions are available
    # TODO: how to tell if Abacus was compiled in double precision?
    if not do_fake_convolve:
        derivatives.make_derivatives(param, floatprec=True, twoD=twoD)

    if make_ic and not zeldovich.is_on_the_fly_format(param.ICFormat):
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
            distribute_state_cmd = [pjoin(abacuspath, 'Abacus', 'move_node_states.py'), paramfn, distribute_state_from, '--distribute-from-serial']
            call_subprocess(Conv_mpirun_cmd + distribute_state_cmd, timeout=param.get('checkpoint_timeout'))

        #check if our previous run was interrupted and saved in the global directory. If yes, redistribute state to nodes.
        #TODO do this by checking if we have a backed-up state available in global directory (instead of looking at param file).

        if not make_ic:  #if this is not a clean run, redistribute the state out to the nodes.
            allow_new_nnode = param.get('AllowResumeToNewNNode',False)
            print('Distributing in order to resume...')
            dist_start = wall_timer()
            distribute_state_cmd = [pjoin(abacuspath, 'Abacus', 'move_node_states.py'), paramfn, '--distribute', resume_dir, '--verbose']
            if param.get('SafeCP', False):
                distribute_state_cmd += ['--safe-cp']
            if allow_new_nnode:
                distribute_state_cmd += ['--allow-new-nnode']
            call_subprocess(Conv_mpirun_cmd + distribute_state_cmd, timeout=param.get('checkpoint_timeout'))
            dist_time = wall_timer() - dist_start
            status_log.print(f"# Distributing state to resume took {dist_time:.1f} seconds.")

            if allow_new_nnode:
                recover_multipoles(paramfn, -1, singlestep_env, mpirun_cmd=Conv_mpirun_cmd)

    print("Beginning abacus steps:")
    if parallel:
        # TODO: this would be useful in the serial code too, especially with ramdisk runs
        emergency_exit_fn = pjoin(param['WorkingDirectory'], 'ABANDON_SHIP')
        print("\n------------------")
        print("To trigger emergency quit safely, create file", emergency_exit_fn)
        print("------------------")

    singlestep_cmd = [pjoin(abacuspath, "singlestep", "singlestep"), paramfn, str(int(make_ic))]
    if parallel:
        singlestep_cmd = mpirun_cmd + singlestep_cmd
        print(f"Using invocation '{' '.join(singlestep_cmd)}'")

        # #if this job is longer than a set time (NEEDS_INTERIM_BACKUP_MINS), we'll need to do a backup halfway through the run.
        interim_backup_complete = run_time_minutes <= NEEDS_INTERIM_BACKUP_MINS

    last_backup = wall_timer()  # time of last backup. Count the very beginning as a "backup"
    for i in range(maxsteps):

        if make_ic:
            ConvDone = True  # No need to convolve for an IC step
            stepnum = 0
        else:
            try:
                read_state = InputFile(pjoin(read,"state"))
                if parallel:
                    # The current merged singlestep + convolution code doesn't need a convolve in between steps
                    ConvDone = True
                else:
                    ConvDone = check_multipole_taylor_done(param, read_state, kind='Taylor')
                stepnum = read_state.FullStepNumber + 1
            except:
                # TODO: might have no read state if we're restoring from a backup
                ConvDone = False
                stepnum = 0

        ss_timer, conv_time = None, None

        make_log_dir_for_step(param['LogDirectory'], stepnum)

        # Do the convolution
        # TODO: do we want to continue to support tick-tock parallel convolve/singlestep?
        if not ConvDone:
            if not do_fake_convolve:
                # Now check if we have all the multipoles the convolution will need
                if not check_multipole_taylor_done(param, read_state, kind='Multipole'):
                    recover_multipoles(paramfn, stepnum, singlestep_env, timeout=param.get('step_timeout'))

                # Swap the Taylors link.  In effect, this will place the Taylors on the same disk as the multipoles.
                # But that's what we want for sloshing: this was the write disk, so now it will be the upcoming read disk
                # We'll swap the multipoles as soon as we convolve
                if path.islink(taylors):
                    taylors_convstate = os.readlink(taylors)
                    os.unlink(taylors)
                    os.symlink(os.readlink(multipoles), taylors)

            if do_fake_convolve:
                convolution_cmd = [pjoin(abacuspath, "Convolution", "FakeConvolution.py"), paramfn]
            else:
                convolution_cmd = [pjoin(abacuspath, "Convolution", "convolution"), paramfn]

            print(f"Performing convolution for step {stepnum:d}")
            with Tools.ContextTimer() as conv_timer:
                call_subprocess(convolution_cmd,
                                env=convolution_env,
                                timeout=param.get('step_timeout'),
                                reset_affinity=param.get('ResetAffinity', True),
                                )
            conv_time = conv_timer.elapsed

            if ProfilingMode == 2:
                print(f'\tConvolution step {stepnum} finished. ProfilingMode = 2 (Convolution) is active; Abacus will now quit.')
                break

            if not check_multipole_taylor_done(param, read_state, kind='Taylor'):
                # have to take it on faith in the parallel version!
                raise ValueError("Convolution did not complete")

            # TODO: pass step num to convolution so it can name the logs appropriately
            convlogs = glob(pjoin(param.LogDirectory, 'last.*conv*'))
            for _cl in convlogs:
                cl = _cl.replace('last', f'step{read_state.FullStepNumber+1:04d}') 
                os.rename(_cl, cl)
                logdir_name = f'step{read_state.FullStepNumber+1:04d}'
                logdir = pjoin(param.LogDirectory, logdir_name)
                shutil.move(cl, pjoin(logdir, basename(cl)))

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
            print(f'Running parallel convolution + singlestep for step {stepnum:d}.')
        else:
            print(f"Running singlestep for step {stepnum:d}")

        with Tools.ContextTimer() as ss_timer:
            try:
                call_subprocess(singlestep_cmd,
                                env=singlestep_env,
                                timeout=param.get('step_timeout'),
                                reset_affinity=param.get('ResetAffinity', True),
                                )
            except subprocess.CalledProcessError as cpe:
                handle_singlestep_error(cpe)
                raise

        # In profiling mode, we don't move the states so we can immediately run the same step again
        if ProfilingMode and ProfilingMode != 2:
            print(f'\tStep {stepnum} finished. ProfilingMode is active; Abacus will now quit and states will not be moved.')
            break

        # Check that write/state was written as a test of success
        write_state_path = pjoin(write, "state")
        if not path.isfile(write_state_path):
            raise ValueError(f'No write state file at "{write_state_path}"; singlestep did not complete!')
        write_state = InputFile(write_state_path)

        # Update the status log
        if stepnum == 0:
            # read_state doesnt exist yet, so pass in the write state instead
            status_log.update(param, write_state, write_state, ss_timer.elapsed, conv_time)
        else:
            status_log.update(param, read_state, write_state, ss_timer.elapsed, conv_time)


        shutil.copy(pjoin(write, "state"),
            pjoin(param.LogDirectory, "step{0:04d}", "step{0:04d}.state").format(write_state.FullStepNumber))
        print("\tFinished state {:d}. a = {:.4f}, dlna = {:.3g}, rms velocity = {:.3g}".format(
            write_state.FullStepNumber, write_state.ScaleFactor,
            write_state.DeltaScaleFactor/(write_state.ScaleFactor-write_state.DeltaScaleFactor),
            write_state.RMS_Velocity))

        if not parallel:
            # Delete the Taylors right away; want to avoid accidentally
            # running with inconsistent positions and Taylors
            remove_MT(param.TaylorDirectory, 'Taylor_*')
            if 'TaylorDirectory2' in param:
                remove_MT(param.TaylorDirectory2, 'Taylor_*')

        # Check if we need to back up the write state
        if backups_enabled and write_state.FullStepNumber > 0 and write_state.FullStepNumber % param['BackupStepInterval'] == 0:
            # TODO: this only works for non-parallel sims and doesn't cover all valid directory structures
            backup_src = param['WorkingDirectory']
            backup_dir = pjoin(param['BackupDirectory'], 'backup')
            status_log.print(f'# Backing up state to {backup_dir}...', tee=True, end='')
            tmp_backup = pjoin(param['BackupDirectory'], 'backup_inprogress')
            shutil.rmtree(tmp_backup, ignore_errors=True)
            backup_concurrency = param.get('BackupConcurrency', 1)
            with Tools.ContextTimer() as backup_timer, Tools.MultithreadedCopier(backup_concurrency) as copier:
                shutil.copytree(backup_src, tmp_backup, symlinks=True, copy_function=copier.copy)
            status_log.print(f' done in {backup_timer.elapsed:.1f} s.', tee=True)
            if path.isdir(backup_dir):
                shutil.rmtree(backup_dir)
            shutil.move(tmp_backup, backup_dir)

            if param.get('DeleteICsAfterFirstBackup', 0) and write_state.FullStepNumber == param['BackupStepInterval']:
                try:
                    shutil.rmtree(pjoin(param['InitialConditionsDirectory']))
                except FileNotFoundError:
                    pass

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
                shutil.copy(dfn, pjoin(param.LogDirectory, "step{0:04d}", "step{0:04d}.density").format(write_state.FullStepNumber))
                np.savez(pjoin(param.LogDirectory, "step{0:04d}", "step{0:04d}.pow").format(write_state.FullStepNumber), k=k,P=P,nb=nb)
                os.remove(dfn)
        
        out_of_time = (wall_timer() - start_time >= run_time_secs)
        if not parallel and out_of_time:
            print('Running out of time in serial job! Exiting.')
            ending_time = time.time()
            ending_time_str = time.asctime(time.localtime())
            ending_time = (ending_time-starting_time)/3600.0    # Elapsed hours
            status_log.print(f"# Serial job terminating prematurely w/ code 0 due to running out of time in job.  {ending_time_str:s} after {ending_time:f} hours.")
            return 0 

        if parallel:
            # Merge all nodes' checksum files into one
            merge_checksum_files(write_state.NodeSizeX * write_state.NodeSizeZ, param)

            # check if we need to gracefully exit and/or backup:
            #    is a group finding step next?          ---> backup to NVME and continue.   (group finding steps are more likely to crash!)
            #    are we ~halfway through our long run?  ---> backup to NVME and continue.   (just in case.)
            #    are we running out of time in the job? ---> backup to global dir and exit.
            #    did the user ask to abandon ship?      ---> backup to global dir and exit.

            abandon_ship = path.exists(emergency_exit_fn)
            
            #if the halfway-there backup hasn't been done already, and we're more than halfway through the job, backup the state now: 
            interim_backup = (not interim_backup_complete) and (wall_timer() - start_time >= run_time_secs / 2) and (run_time_secs > NEEDS_INTERIM_BACKUP_MINS * 60 ) #only relevant for long jobs.
            
            #are we coming up on a group finding step? If yes, backup the state, just in case. 
            pre_gf_backup  = False 
            nGFoutputs = [] 
            output_arrs = [param.get('L1OutputRedshifts'), param.get('TimeSliceRedshifts'), param.get('TimeSliceRedshifts_Subsample')]
            for output_arr in output_arrs:
                try:
                    nGFoutputs.append( len(output_arr) )
                except (AttributeError,TypeError):
                    nGFoutputs.append(0) 

            if (run_time_secs > NEEDS_INTERIM_BACKUP_MINS * 60): 
                for i in range(len(nGFoutputs)):
                    for z in range(nGFoutputs[i]):
                        L1z = output_arrs[i][z] 

            # #if the halfway-there backup hasn't been done already, and we're more than halfway through the job, backup the state now:
            # interim_backup = (not interim_backup_complete) and (wall_timer() - start_time >= run_time_secs / 2) and (run_time_secs > NEEDS_INTERIM_BACKUP_MINS * 60 ) #only relevant for long jobs.

            #are we coming up on a group finding step? If yes, backup the state, just in case.
            pre_gf_backup  = False
            output_arrs = []
            for n in ('L1OutputRedshifts', 'TimeSliceRedshifts', 'TimeSliceRedshifts_Subsample'):
                oa = param.get(n)
                if oa is not None:
                    try:
                        output_arrs += [list(oa)]
                    except TypeError:
                        output_arrs += [[oa]]

            if run_time_secs > NEEDS_INTERIM_BACKUP_MINS * 60:
                for oa in output_arrs:
                    for L1z in oa:
                        if L1z <= -1:
                            continue
                        L1a = 1/(1 + L1z)

                        # here we assume that the next da will be similar to the previous one.
                        # we don't want to backup before every GF step -- that's too expensive.
                        # so let's check if it's been more than two hours since the previous backup.
                        if (write_state.Redshift > L1z + 1e-12) and (1/(1 + write_state.Redshift) + write_state.DeltaScaleFactor > L1a):
                            last_z = L1z == 0.1 or L1z == 0.0  # TODO: magic numbers
                            time_for_gf_backup = (wall_timer() - last_backup) > GF_BACKUP_INTERVAL

                            if last_z or time_for_gf_backup:
                                pre_gf_backup = True
                                break

                    if pre_gf_backup:
                        break

            try:
                intermittent_backup = wall_timer() - last_backup >= BACKUP_EVERY_N_MINUTES*60
            except:
                intermittent_backup = False

            exiting = out_of_time or abandon_ship
            #save = exiting or interim_backup or pre_gf_backup
            save = exiting or pre_gf_backup or intermittent_backup

            if save and param.get('SaveState',True):
                if abandon_ship:
                    exit_message = 'EMERGENCY EXIT REQUESTED'
                if out_of_time:
                    exit_message = 'RUNNING OUT OF JOB TIME'
                # if interim_backup:
                #     exit_message = 'HALFWAY THROUGH A LONG JOB'
                if intermittent_backup:
                    exit_message = f'{BACKUP_EVERY_N_MINUTES} MIN BACKUP'
                if pre_gf_backup:
                    exit_message = 'GROUP FINDING COMING UP'

                print(exit_message, 'backing up node state.')
                status_log.print(f'# {exit_message:s}: backing up node state.')

                restore_time = wall_timer()
                retrieve_state_cmd = [pjoin(abacuspath, 'Abacus', 'move_node_states.py'), paramfn, resume_dir, '--retrieve']
                if param.get('DeleteICsAfterFirstBackup'):
                    retrieve_state_cmd += ['--delete-ics']
                call_subprocess(Conv_mpirun_cmd + retrieve_state_cmd, timeout=param.get('checkpoint_timeout'))

                restore_time = wall_timer() - restore_time

                print(f'Retrieving and storing state took {restore_time} seconds. ', end = '')
                status_log.print(f'# Retrieving and storing state took {restore_time:.1f} seconds.')
                last_backup  = wall_timer()  #log the time of the last backup.


                #checking if path exists explicitly just in case user requested emergency exit while we were retrieveing the state.
                if path.exists(emergency_exit_fn):
                    os.remove(emergency_exit_fn)

        # Now shift the states down by one
        move_state_dirs(read, write, past)

        # check if we've reached the stopping redshift
        # Note that FinalRedshift always trumps,
        # whether it is before or after the TimeSliceZ array
        try:
            finalz = param.FinalRedshift
        except AttributeError:
            # param.FinalRedshift wasn't set, so we'll look for the minimum in the output redshifts arrays.
            redshift_arrs = ['TimeSliceRedshifts', 'TimeSliceRedshifts_Subsample', 'L1OutputRedshifts']
            finalz = None
            for name in redshift_arrs:
                if name in param:
                    _zmin = np.min(param[name])
                    finalz = min(finalz, _zmin) if finalz is not None else _zmin
            if finalz is None:
                finalz = 0.


        # This logic is deliberately consistent with singlestep.cpp
        # If this is an IC step then we won't have read_state
        if (not make_ic and
            np.abs(read_state.Redshift - finalz) < 1e-12 and
            read_state.LPTStepNumber == 0 and
            not (write_state.DoGroupFindingOutput == 1 and write_state.DidGroupFindingOutput == 0)
            ):
            ending_time = time.time()
            ending_time_str = time.asctime(time.localtime())
            ending_time = (ending_time-starting_time)/3600.0    # Elapsed hours
            print(f"Final redshift of {finalz:g} reached; terminating normally after {ending_time:.2f} hours.")
            status_log.print(f"# Final redshift of {finalz:g} reached at {ending_time_str:s}; terminating normally after {ending_time:.2f} hours.")
            finished = True
            break

        if parallel:
            if exiting:
                print('Exiting.')
                if not abandon_ship and maxsteps == 10000:
                    print('Requeueing!')
                    ending_time = time.time()
                    ending_time_str = time.asctime(time.localtime())
                    ending_time = (ending_time-starting_time)/3600.0    # Elapsed hours
                    status_log.print(f"# Terminating normally w/ requeue code.  {ending_time_str:s} after {ending_time:.2f} hours.")
                    return EXIT_REQUEUE
                else:
                    if abandon_ship:
                        print('Requeue disabled because ABANDON_SHIP was requested.')
                    else:
                        print(f'Requeue disabled because maxsteps={maxsteps} was specified.')
                    return 0
            else:
                print('Continuing run.')


        if parallel and abandon_ship:
            ending_time = time.time()
            ending_time_str = time.asctime(time.localtime())
            ending_time = (ending_time-starting_time)/3600.0    # Elapsed hours
            print(f"Emergency exit triggered at {ending_time_str:s}; abandoning ship after {ending_time:.2f} hours.")
            status_log.print(f"# Emergency exit triggered at {ending_time_str:s}; abandoning ship after {ending_time:.2f} hours.")
            os.remove(emergency_exit_fn)
            break

        make_ic = False

    if maxsteps != 10000: #we asked to do only a limited number of steps, and we've successfully completed them. We're done.
        ending_time = time.time()
        ending_time_str = time.asctime(time.localtime())
        ending_time = (ending_time-starting_time)/3600.0    # Elapsed hours
        print(f"Reached maxsteps limit at {ending_time_str:s}; exiting job w/o requeue after {ending_time:.2f} hours.")
        status_log.print(f"# Reached maxsteps limit at {ending_time_str:s}; exiting job w/o requeue after {ending_time:.2f} hours.")
        finished = True

    if finished and parallel and param.get('SaveFinalState',False):
        retrieve_state_cmd = [pjoin(abacuspath, 'Abacus', 'move_node_states.py'), paramfn, resume_dir, '--retrieve']
        if param.get('DeleteICsAfterFirstBackup'):
            retrieve_state_cmd += ['--delete-ics']
        restore_time = wall_timer()
        call_subprocess(Conv_mpirun_cmd + retrieve_state_cmd, timeout=param.get('checkpoint_timeout'))
        restore_time = wall_timer() - restore_time
        print(f'Retrieving and storing final state took {restore_time:.1f} seconds. ', end = '')
        status_log.print(f'# Retrieving and storing final state took {restore_time:.1f} seconds.')

    if finished and not parallel and not param.get('SaveFinalState',False):
        backup_dir = pjoin(param['BackupDirectory'], 'backup')
        if path.isdir(backup_dir):
            shutil.rmtree(backup_dir)

    # If there is more work to be done, signal that we are ready for requeue
    if not finished and not ProfilingMode:
        ending_time = time.time()
        ending_time_str = time.asctime(time.localtime())
        ending_time = (ending_time-starting_time)/3600.0    # Elapsed hours
        status_log.print(f"# Terminating normally.  {ending_time_str:s} after {ending_time:.2f} hours.")
        print(f"About to return EXIT_REQUEUE code {EXIT_REQUEUE}")
        return EXIT_REQUEUE


    return 0


def recover_multipoles(paramfn, stepnum, env, mpirun_cmd=None, timeout=None):
    '''
    If multipoles are missing for any reason, we can invoke a redacted
    singlestep pipeline that just reads in positions and outputs multipoles.
    Then Abacus can proceed as normal.

    First, we have to call `make` on the `recover_multipoles` executable,
    because that's not built as part of normal compilation because we don't
    use it often.
    '''
    # Invoke multipole recovery mode
    print(f"Warning: missing multipoles! Performing multipole recovery for step {stepnum:d}")

    # Build the recover_multipoles executable
    with Tools.chdir(pjoin(abacuspath, "singlestep")):
        call_subprocess(['make', 'recover_multipoles'])

    # Execute it
    print(f"Running recover_multipoles for step {stepnum:d}")

    recov_cmd = [pjoin(abacuspath, "singlestep", "recover_multipoles"), paramfn]

    if mpirun_cmd:
        recov_cmd = mpirun_cmd + [pjoin(abacuspath, "singlestep", "recover_multipoles"), paramfn]

    try:
        call_subprocess(recov_cmd, env=env, timeout=timeout)
    except subprocess.CalledProcessError as cpe:
        handle_singlestep_error(cpe)
        raise

    print(f'\tFinished multipole recovery for step {stepnum}.')


def make_log_dir_for_step(params_logdir, step, symlink_name="last"):
    '''
    We store the log files for each step in their own
    directory.  We also supply a directory symlink that
    points to the most recent log directory.  This function
    creates the new step-numbered log directory and updates
    the symlink to point to it.  Should be called before each
    step.
    '''

    logdir_name = f'step{step:04d}'
    logdir = pjoin(params_logdir, logdir_name)
    # the logdir might exist already
    # if we're profiling and running the same step over and over,
    # for example
    os.makedirs(logdir, exist_ok=True)

    symlink = pjoin(params_logdir, symlink_name)
    if path.exists(symlink):
        assert path.islink(symlink)
    try:
        os.unlink(symlink)
    except FileNotFoundError:
        pass

    os.symlink(logdir_name, symlink)


def merge_checksum_files(nnodes, param=None, dir_globs=None, checksum_fn='checksums.crc32'):
    '''
    Each node computes the CRC32 checksum of its output files.
    Rather than doing a messy MPI merge of the checksums so
    that the root node can write one file, for now each node
    just writes its own file.  So the last step is to merge
    these files into one file, which this function does.

    We basically guess which directories might have checksum
    files.  It's not ideal.  The directories should match
    the types in SlabBuffer::WantChecksum().

    The checksums.crc32 file written by this function
    should match the format of the GNU cksum util.

    We use `nnodes` to construct the filenames of the
    files to check, rather than do a potentially expensive glob.

    Also, if `checksum_fn` already exists, we will assume
    that we have already merged all checksum files in this directory.
    '''

    if not dir_globs:
        dir_globs = [pjoin(param.OutputDirectory, 'slice*'),
                    pjoin(param.GroupDirectory, 'Step*'),
                    pjoin(param.LCDirectory, 'Step*'),
            ]

    for pat in dir_globs:
        for d in glob(pat):
            if path.isfile(checksum_fn):
                # Already merged!
                continue

            cksum_fns = []
            for i in range(nnodes):
                nodefn = pjoin(d,f'checksums.{i:04d}.crc32')
                if path.isfile(nodefn):
                    cksum_fns += [nodefn]

            # Even if empty, process anyway so we write the checksum_fn
            # as a signal of completion

            lines = []
            for fn in cksum_fns:
                with open(fn, 'r') as fp:
                    lines += fp.readlines()

            # Sort on the file name
            lines = [line.split() for line in lines]
            assert(all(len(line) == 3 for line in lines))
            lines = sorted(lines, key=lambda l:l[2])
            lines = [' '.join(line) + '\n' for line in lines]

            with open(pjoin(d,checksum_fn), 'a') as fp:
                fp.writelines(lines)

            for fn in cksum_fns:
                os.remove(fn)



def handle_singlestep_error(error):
    '''
    The singlestep executable may quit with an informative return value,
    such as from a Unix signal.  Let's interpret that for the user.
    The error is not suppressed here and can/should still be raised
    upstream.

    Parameters
    ----------
    error: subprocess.CalledProcessError
        The singlestep exception object
    '''
    if error.returncode == -signal.SIGBUS:
        print('singlestep died with signal SIGBUS! Did the ramdisk run out of memory?', file=sys.stderr)


def showwarning(message, category, filename, lineno, file=None, line=None):
    '''
    A monkey patch for warnings.showwarning that is a little easier to parse.
    '''

    if file is None:
        file = sys.stderr

    print(f'{filename}:{lineno}: {category.__name__}: {message}', file=file)

warnings.showwarning = showwarning
