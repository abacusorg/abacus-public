"""
A collection of file organization functions
to standardize paths for analysis products.
"""
import os
import os.path as path
import tarfile
import shutil
from glob import glob
from itertools import izip

from Abacus import Tools

def get_output_dir(product_name, slice_dir):
    """
    The directory in which to store outputs for this product.
    """
    # Rename 'slice' to 'z'
    dirname, slicename = path.split(path.abspath(slice_dir))
    slicename = slicename.replace('slice', 'z')
    
    # Turn 'sim' into 'sim_products/sim_[product]'
    dirname, simname = path.split(path.abspath(dirname))
    simname = path.join(simname + '_products', simname + '_' + product_name)
    
    # Join
    outdir = path.join(dirname, simname, slicename, '')  # include trailing slash
    return outdir

import multiprocessing
def make_tar(dir, pattern, tarfn, delete_source=False, nthreads=1):
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
        Maximum number of threads to use.  Default: 1.
    
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
    
    pool = multiprocessing.Pool(nthreads)
    allfns = pool.map(_make_tar_worker(delete_source), zip(dir, pattern, tarfn))
    fns = sum(allfns, [])
    
    return fns

class _make_tar_worker(object):
    def __init__(self, delete_source):
        self.delete_source = delete_source
    
    def __call__(self, args):
        dir, pattern, tarfn = args
        with Tools.chdir(dir):
            fns = glob(pattern)
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
