#!/usr/bin/env python

"""
A collection of file organization functions
to standardize paths for analysis products.
"""
import os
import os.path as path
from os.path import join as pjoin
import tarfile
import shutil
from glob import glob
from warnings import warn

import contextlib

from Abacus import Tools
from Abacus.Tools import ContextTimer
from Abacus.InputFile import InputFile
from Abacus import ReadAbacus

def get_output_dir(product_name, slice_dir, out_parent=None):
    """
    The directory in which to store outputs for this product.
    """
    # Rename 'slice' to 'z'
    dirname, slicename = path.split(path.abspath(slice_dir))
    slicename = slicename.replace('slice', 'z')
    
    # Turn 'sim' into 'sim_products/sim_[product]'
    dirname, simname = path.split(path.abspath(dirname))
    simname = path.join(simname + '_products', simname + '_' + product_name)
    
    # Replace the parent directory with the given one
    if out_parent:
        dirname = path.abspath(out_parent)
        
    # Join
    outdir = path.join(dirname, simname, slicename)  # no trailing slash
    return outdir
    

def prompt_removal(dir, make_new=False, noprompt=False):
    """
    Ask the user if they want to delete `dir`.  Fail if they do no answer yes.
    Optionally create an empty directory with the same name after the removal.
    """
    if path.exists(dir):
        if not noprompt:
            yn = input('Output directory "{}" exists.  Delete (y/[n])? '.format(dir))
        if noprompt or yn.lower() in ('y', 'yes'):
            print('Removing {}'.format(dir))
            shutil.rmtree(dir)
        else:
            raise RuntimeError('Cannot continue if "{}" exists.'.format(dir))
    
    if make_new:
        os.makedirs(dir)


gadget_pattern = r'*.[0-9]*'


def get_header(dir, retfn=False):
    # Read the header to get the boxsize
    header_pats = [pjoin(dir, 'header'), pjoin(dir, os.pardir, 'info/*.par'), pjoin(dir, os.pardir, '*.par')]
    headers = sum([glob(h) for h in header_pats], [])

    for header_fn in headers:
        try:
            header = InputFile(header_fn)
            BoxSize = header['BoxSize']
        except IOError:
            continue
        break
    else:
        for f in (get_desi_hdf5_header, get_gadget_header):
            try:
                header = f(dir)
                header_fn = None  # safe to return None?
            except:
                raise
                continue
            break
        else:
            warn(f'Could not find a header in {dir}')
            header, header_fn = None, None
            

    if retfn:
        return header, header_fn
    return header


def get_gadget_prefix(dir):
    '''
    Gadget snapshots consist of many files with names {prefix}.0, {prefix}.1, etc.
    This tries to extract just the prefix part, which can be passed to pynbody.load().
    '''
    gadget_prefix = glob(pjoin(dir, gadget_pattern))[0]
    gadget_prefix = gadget_prefix[:gadget_prefix.rfind('.')]
    return gadget_prefix

def get_gadget_header(dir):
    gadget_fn = sorted(glob(pjoin(dir, gadget_pattern)))[0]

    # Extract the relevant properties from the gadget header and return them in a dict
    import pynbody
    f = pynbody.load(gadget_fn)
    header = {'BoxSize': float(f.properties['boxsize']),
              'ScaleFactor': float(f.properties['a'])}
    return header


def get_desi_hdf5_header(dir):
    import h5py

    h5pat = ReadAbacus.default_file_patterns['desi_hdf5']
    h5fn = glob(pjoin(dir, h5pat))[0]

    with h5py.File(h5fn, 'r') as fp:
        header = dict(fp['/AbacusHeader'].attrs)
    return header


import multiprocessing
def make_tar(dir, pattern, tarfn, delete_source=False, nthreads=1, out_parent=''):
    """
    Makes a compressed tar file named `tarfn` inside each
    directory specified in `dirs`.  Files and directories
    that match `pattern` are included in the tarball.
    
    If multiple dirs are given, then `nthreads` can be used
    to run tar in parallel. This can result in a substantial
    speedup, as tar-making is always CPU, not IO, limited.

    Parameters
    ----------
    dirs: str or list of str
        The directories to processes
    pattern: str or list of str
        The globbing pattern(s) to apply inside each dir
    tarfn: str or list of str
        The tar filename(s) to add globbed files to
    delete_source: bool, optional
        Delete source files after finishing writing
        the tar.  Default: False.
    nthreads: int, optional
        Maximum number of threads to use.  `nthreads` <= 0 means use all cores.
        Default: 1.
    out_parent: str, optional
        Place the outputs in a directory rooted at out_parent.
        Default is to make the tarballs in-place.
    
    Returns
    -------
    fns: list of str
        The absolute filenames of any files that were
        added to a tar
    """
    # Make sure all args have the same length
    if type(dir) is str:
        dir = [dir]
    if type(pattern) is str:
        pattern = [pattern]
    if type(tarfn) is str:
        tarfn = [tarfn]
    maxlen = max(len(dir), len(pattern), len(tarfn))
    if len(dir) == 1:
        dir *= maxlen
    if len(pattern) == 1:
        pattern *= maxlen
    if len(tarfn) == 1:
        tarfn *= maxlen
    assert len(dir) == len(pattern) == len(tarfn)
    
    if nthreads <= 0:
        nthreads = multiprocessing.cpu_count()
    
    pool = multiprocessing.Pool(nthreads)
    allfns = pool.map(_make_tar_worker(delete_source, out_parent), list(zip(dir, pattern, tarfn)))
    fns = sum(allfns, [])
    
    return fns

class _make_tar_worker(object):
    def __init__(self, delete_source, out_parent):
        self.delete_source = delete_source
        self.out_parent = out_parent
    
    def __call__(self, args):
        dir, pattern, tarfn = args
        
        # Make the parent directories for the output file
        tarfn = path.abspath(path.join(self.out_parent, dir, tarfn))
        try:
            os.makedirs(os.path.dirname(tarfn))
        except OSError as exc:
            pass
            
        with Tools.chdir(dir):
            fns = glob(pattern)
            if not fns:
                return fns
            with tarfile.open(tarfn, 'w:gz') as tar:
                for fn in fns:
                    tar.add(fn)
            if self.delete_source:
                for fn in fns:
                    try:
                        os.remove(fn)
                    except OSError:  # Maybe a dir?
                        shutil.rmtree(fn)
            # make absolute while still in chdir context
            fns = [path.abspath(fn) for fn in fns]
        return fns

def get_slab_fmt_str(slice_dir, simname):
    slice_z = path.basename(path.abspath(slice_dir)).replace('slice', 'z')
    fmt_str = path.join(slice_dir, '{}.{}.slab'.format(simname, slice_z) + '{slab:04d}.dat')
    return fmt_str


@contextlib.contextmanager
def extract_slabs(dir, verbose=True, tarfn='slabs.tar.gz'):
    '''
    The directory we want to process on may have its
    slabs archived in a tarball.  First check if all
    the files we may need are already available.
    If not, we want to extract these files
    but clean up when leaving the context.
    
    Assumes tarballs do not contain directories.
    
    Parameters
    ----------
    dir: str
        The directory 
    tarfn: str, optional
        The name of the tarfile in `dir` to expand
        Default: 'slabs.tar.gz'
    verbose: bool, optional
        Print the extraction time.  Default: True
    '''
    header = InputFile(path.join(dir, 'header'))
    cpd = int(round(header.CPD))
    fmt_str = get_slab_fmt_str(dir, header.SimName)
    
    # Check that all the files we need exist
    for slab in range(cpd):
        fn = fmt_str.format(slab=slab)
        if not path.isfile(fn):
            break
    else:
        # All files present; do nothing
        yield
        return
    
    if verbose:
        print('Extracting slabs...')
    # If something is missing, extract the tarball
    extract_timer = ContextTimer()
    with extract_timer:
        with tarfile.open(path.join(dir, tarfn)) as tfp:
            tfp.extractall(path=dir)
            fns = tfp.getnames()  # this is smart enough to cache the result from extractall()!
    if verbose:
        print('Extract time: {:.3f}s'.format(extract_timer.elapsed))
    yield

    # Now remove the extracted files
    for fn in fns:
        os.remove(path.join(dir,fn))
        
import argparse
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Parallel tar utility', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('dir', help='The directories to process', nargs='+')
    parser.add_argument('pattern', help='The glob pattern of files to put in the tarball')
    parser.add_argument('tarfn', help='The name of the tar file to create in each directory')
    parser.add_argument('--nthreads', help='Number of threads on this node to divide the problem over', default=1, type=int)
    parser.add_argument('--delete-source', action='store_true', help='Remove the files that were placed the in the tar.')
    parser.add_argument('--out-parent', help='Place the outputs in a directory rooted at out_parent. Default is to make the tarballs in-place.', default='')
    
    args = parser.parse_args()
    
    make_tar(**vars(args))
