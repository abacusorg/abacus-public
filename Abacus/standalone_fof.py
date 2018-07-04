#!/usr/bin/env python
'''
Run the Abacus on-the-fly group finder on the given time slice(s).

This invokes the `singlestep/standalone_fof` binary, which contains
a redacted dependency pipeline that reads in pack14 slabs and
does group finding.

The binary takes the slice directory and a parameters file.
The parameters file needs to be consistent with the slice
headers, but there are some parameters that will need to be
different if one is running this analysis on a different
computer than the simulation was run on.  So, this script
will scrape the slice header then replace all the directories
and performance tuning parameters (e.g. OMP_NUM_THREADS) with
the site-local values.

'''
import subprocess
import os
import shutil
from os.path import join as pjoin, isdir, basename, dirname, normpath
import argparse

# This is really silly.  Python doesn't allow relative imports from scripts...
from Abacus import abacus
from Abacus import Tools
from Abacus import GenParam
from Abacus.InputFile import InputFile
import Abacus.Analysis.common

def standalone_fof(slicedir, output_paramfn='standalone_fof.par', use_site_overrides=True, override_directories=True):
    slicedir = normpath(slicedir)

    # This parameters file is special: it acts as both the parameters and read state file
    headerfn = pjoin(slicedir, 'header')
    params = abacus.preprocess_params(output_paramfn, headerfn, use_site_overrides=use_site_overrides, override_directories=override_directories)

    # Make an output directory for group outputs and log files
    for d in ['LogDirectory', 'OutputDirectory', 'GroupDirectory']:
            if d in params and params[d]:
                params[d] = normpath(params[d]) 
                if not isdir(params[d]):
                    os.makedirs(params[d])

    # set up the final output directory
    this_groupdir_name = 'Step{:04d}_z{:5.3f}'.format(params['FullStepNumber'], params['Redshift'])
    product_dir = Abacus.Analysis.common.get_output_dir('groups', slicedir, out_parent=dirname(dirname(params['GroupDirectory'])))
    dest_groupdir_path = pjoin(product_dir, this_groupdir_name)
    if isdir(dest_groupdir_path):
        raise RuntimeError('Destination directory "{}" already exists!'.format(dest_groupdir_path))

    # Copy the parameter file to the simulation dir
    shutil.move(output_paramfn, pjoin(params['GroupDirectory'], output_paramfn))

    with Tools.chdir(pjoin(abacus.abacuspath, 'singlestep')):
        subprocess.check_call(['make', 'standalone_fof'])
        subprocess.check_call(['./standalone_fof', slicedir, pjoin(params['GroupDirectory'], output_paramfn)])

    # Now move the results to a _products directory
    if not isdir(dirname(product_dir)):
        os.makedirs(dirname(product_dir))

    this_groupdir_path = pjoin(params['GroupDirectory'], this_groupdir_name)
    shutil.move(this_groupdir_path, product_dir)

    # copy over the info dir
    dest_info_dir = pjoin(dirname(product_dir), 'info')
    if not isdir(dest_info_dir):
        shutil.copytree(pjoin(dirname(slicedir), 'info'), pjoin(dirname(product_dir), 'info'))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('slice', nargs='+')
    parser.add_argument('--no_site_overrides', action='store_true', default=False)
    parser.add_argument('--no_dir_overrides', action='store_true', default=False)
    args = parser.parse_args()

    use_site_overrides = not args.no_site_overrides
    override_directories = not args.no_dir_overrides

    for sl in args.slice:
        standalone_fof(sl, use_site_overrides=use_site_overrides, override_directories=override_directories)
