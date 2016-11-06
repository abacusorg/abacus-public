"""
A collection of file organization functions
to standardize paths for analysis products.
"""
import os
import os.path
import tarfile
import shutil

def get_output_dir(product_name, slice_dir):
    """
    The directory in which to store outputs for this product.
    """
    # Rename 'slice' to 'z'
    dirname, slicename = os.path.split(os.path.abspath(slice_dir))
    slicename = slicename.replace('slice', 'z')
    
    # Turn 'sim' into 'sim_products/sim_[product]'
    dirname, simname = os.path.split(os.path.abspath(dirname))
    simname += '_products/' + simname + '_' + product_name
    
    # Join
    outdir = os.path.join(dirname, simname, slicename)
    return outdir


def make_tar(dirs, pattern, tarfn, delete_source=False):
    """
    Make tar file `tarfn` from files
    that match `pattern` inside any of `dirs`.
    `delete_source` deletes the source files afterwards.
    
    For easy suffix attachment, guaranteed not to have
    a trailing slash.
    """
    for d in dirs:
        with Tools.chdir(d):
            fns = glob(pattern)
            with tarfile.open(tarfn, 'w:gz') as tar:
                for fn in fns:
                    tar.add(fn)
            if delete_source:
                for fn in fns:
                    try:
                        os.remove(fn)
                    except OSError:  # Maybe a dir?
                        shutil.rmtree(fn)

