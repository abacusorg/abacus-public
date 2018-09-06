"""
Tools for efficient reading of various Abacus file formats in Python.
See the Abacus readme for details on the various formats.

The main entry points are the `read()` and `from_dir()` functions
to read a single file or a while directory of files, respectively.
They are just wrappers to call the appropriate reader function
for the given file format, like `read_pack14()` or `read_rvzel()`.

To use the fast pack14 reader, you'll need to build the underlying
library with
`$ make analysis`
from the top level Abacus directory.  Or just build all of Abacus.
The alternative is to use the slow `read_pack14_lite()` function.

Todo:
- Unit conversion options
- native ctypes float3* instead of void* (or CFFI?)
- inline particle downsampling
"""

import os
import os.path as path
from glob import glob
import ctypes as ct
import re

import numpy as np
import numba

from . import abacus
from .InputFile import InputFile
from .Tools import ndarray_arg, asciistring_arg

ralib = ct.cdll.LoadLibrary(path.join(abacus.abacuspath, "clibs", "libreadabacus.so"))

# Set up the arguments and return type for the library functions
# TODO: should consider migrating this to CFFI
for f in (ralib.read_pack14, ralib.read_pack14f):
    f.restype = ct.c_uint64
    f.argtypes = (asciistring_arg, ct.c_size_t, ct.c_int, ct.c_int, ct.c_int, ct.c_int, ndarray_arg)
    

def read(*args, **kwargs):
    """
    A convenience function to read particle data from
    file type `format`.  Simply calls
    the appropriate `read_[format]()` function.

    For usage, see the docstring of the reader
    function for the file format you are using.
    """
    format = kwargs.pop('format')
    types = {'pack14':read_pack14, 'rvdouble':read_rvdouble,
             'rvzel':read_rvzel, 'state':read_state,
             'rvdoubletag':read_rvdoubletag,
             'pack14_lite':read_pack14_lite}
    format = format.lower()
    return types[format](*args, **kwargs)

def from_dir(dir, pattern='*.dat', key=None, **kwargs):
    """
    Read all files from `dir` that match `pattern`.
    
    Parameters
    ----------
    dir: str
        The directory to read files from
    pattern: str, optional
        A bash globbing pattern to find all the files to read
    format: str, optional
        The file format.  Determines which reader function will be used.
    return_header: bool, optional
        Return the header of the last file parsed
    key: func, optional
        A filename sorting key function.  Default: None.
    **kwargs: dict, optional
        Additional args to pass to the reader function
        
    Returns
    -------
    particles: ndarray
        The concatenated particle array
    header: InputFile, optional
        If `return_header` and a header is found, return parsed InputFile
    """
    _key = (lambda k: key(path.basename(k))) if key else None
    files = sorted(glob(path.join(dir, pattern)), key=_key)
    assert files, "No files found matching {}".format(path.join(dir, pattern))

    return read_many(files, **kwargs)


def read_many(files, format='pack14', separate_fields=False, **kwargs):
    """
    Read a list of files into a contiguous array.

    Our files are usually written as something like float6, where
    pos and vel are interwoven.  If one wants to de-interlave them
    into separate arrays, one can use the `separate_fields` arg.
    The reader functions don't know about this, so they'll still
    return one big array.  But then this function can split the
    fields.  This results in some small extra memory usage and
    time spent copying.  But it's granular at the level of the
    individual file sizes.

    Parameters
    ----------
    files: list of str
        The list of file names
    format: str
        The file format. Specifies the reader function.
    separate_fields: bool, optional
        Split the (e.g.) pos and vel fields of the contiguous array returned
        by the reader function into two separate, contiguous arrays.
        Default: False.
    kwargs: dict, optional
        Additional args to pass to the reader function.

    Returns
    -------
    particles: ndarray or dict of ndarray
        The concatenated particle array, or a dict of arrays if `separate_fields`.
    header: InputFile, optional
        If `return_header` and a header is found, return parsed InputFile
    """

    # Allocate enough space to hold the concatenated particle array
    assert files
    total_fsize = sum(path.getsize(fn) for fn in files)
    
    # State has multiple 'psize_on_disk' values
    format = format.lower()
    _format = format
    if format == 'state':
        dtype = kwargs.get('dtype', np.float32)
        if dtype == np.float32:
            _format = format + 'f'
            
    _format = format
    if _format == 'state' and kwargs.pop('dtype_on_disk', np.float32) == np.float64:
        _format = 'state64'
    alloc_NP = get_np_from_fsize(total_fsize, format=_format)
    outdt = output_dtype(**kwargs)
    
    if separate_fields:
        # One array per dtype field
        particle_arrays = {}
        for field in outdt.descr:
            particle_arrays[field[0]] = np.empty(alloc_NP, dtype=field[1:])
    else:
        particles = np.empty(alloc_NP, dtype=outdt)
    
    return_header = kwargs.get('return_header', False)
    
    start = 0
    for fn in files:
        if separate_fields:
            read_into = None
        else:
            read_into = particles[start:]
        
        out = read(fn, format=format, out=read_into, **kwargs)

        if separate_fields:
            if return_header:
                out, header = out
            NP = len(out)
            # Now split the array
            for field in outdt.descr:
                particle_arrays[field[0]][start:start+NP] = out[field[0]]
            del out
        else:
            if return_header:
                NP, header = out
            else:
                NP = out
        start += NP
        
    # Shrink the array to the size that was actually read
    if separate_fields:
        for field in particle_arrays:
            particle_arrays[field] = particle_arrays[field][:start]

        if return_header:
            return particle_arrays, header  # return header associated with last file
    
        return particle_arrays
    else:
        particles = particles[:start]
    
        if return_header:
            return particles, header  # return header associated with last file
    
        return particles

################################
# Begin list of reader functions
################################

def read_pack14(fn, ramdisk=False, return_vel=True, zspace=False, return_pid=False, return_header=False, dtype=np.float32, out=None):
    """
    Read particle data from a file in pack14 format.
    
    Parameters 
    ----------
    fn: str
        The filename to read
    ramdisk: bool, optional
        Whether `fn` resides on a ramdisk or not.  Necessary to know if we can do directIO.
    return_vel: bool, optional
        Return velocities along with other data
    zspace: bool, optional
        Apply redshift-space distortion to particle positions
    return_pid: bool, optional
        Return particle IDs along with other data
    return_header: bool, optional
        If the pack14 file has an ASCII header, return it as a second return value.
    dtype: data-type, optional
        Either np.float32 or np.float64.  Determines the data type the particle data is loaded into
    out: ndarray, optional
        A pre-allocated array into which the particles will be directly loaded.
        
    Returns
    -------
    data: ndarray of length (npart,)
        The particle data.  Positions are in data['pos']; velocities are in data['vel'] (if `return_vel`),
        and PIDs are in data['pid'] (if `return_pid`).
        
    or,
        
    NP: int
        If `out` is given, returns the number of rows read into `out`.
        
    and optionally,
    
    header: InputFile
        If `return_header` and a header is found, return parsed InputFile
    """
    try:
        ralib
    except NameError:
        raise RuntimeError("pack14 C library was not found. Try building Abacus with 'make analysis'? Or use the slower `read_pack14_lite()` function.")
    readers = {np.float32: ralib.read_pack14f,
               np.float64: ralib.read_pack14 }
    assert dtype in readers
    
    with open(fn, 'rb') as fp:
        header = skip_header(fp)
        if header:
            header = InputFile(str_source=header)
        offset = fp.tell()
    
    # Use the given buffer
    if out is not None:
        _out = out
    else:  # or allocate one
        fsize = path.getsize(fn)
        alloc_NP = get_np_from_fsize(fsize, format='pack14')
        ndt = output_dtype(return_vel=return_vel, return_pid=return_pid, dtype=dtype)
        _out = np.empty(alloc_NP, dtype=ndt)
    
    NP = readers[dtype](fn, offset, ramdisk, return_vel, zspace, return_pid, _out.view(dtype=dtype))
    _out = _out[:NP]  # shrink the buffer to the real size
    
    retval = (NP,) if out is not None else (_out,)
    if return_header:
        retval += (header,)
    
    if len(retval) == 1:
        return retval[0]
    return retval
    
    
def read_rvdouble(fn, return_vel=True, zspace=False, dtype=np.float32, out=None, return_header=False):
    """
    Read particle data from a file in rvdouble format.
    
    Usage is the same as `read_pack14`.
    """
    with open(fn, 'rb') as fp:
        header = skip_header(fp)
        if header:
            header = InputFile(str_source=header)
        data = np.fromfile(fp, dtype=(np.float64,6))
    
    if out is not None:
        _out = out
    else:
        ndt = output_dtype(return_vel=return_vel, dtype=dtype)
        _out = np.empty(len(data), dtype=ndt)
    
    _out['pos'][:len(data)] = data[:,:3]
    if zspace:
        _out['pos'][:len(data)] += data[:,3:]
    if return_vel:
        _out['vel'][:len(data)] = data[:,3:]
    
    retval = (_out,) if out is None else (len(data),)
    
    if return_header:
        retval += (header,)
        
    if len(retval) == 1:
        return retval[0]
    return retval


def read_rvdoubletag(fn, return_vel=True, return_pid=False, zspace=False, dtype=np.float32, out=None, return_header=False):
    """
    Read particle data from a file in rvdoubletag format.
    
    Usage is the same as `read_pack14`.
    """
    with open(fn, 'rb') as fp:
        header = skip_header(fp)
        if header:
            header = InputFile(str_source=header)
        data = np.fromfile(fp, dtype=(np.float64,7))
    
    if out is not None:
        _out = out
    else:
        ndt = output_dtype(return_vel=return_vel, return_pid=return_pid, dtype=dtype)
        _out = np.empty(len(data), dtype=ndt)
    
    _out['pos'][:len(data)] = data[:,:3]
    if zspace:
        _out['pos'][:len(data)] += data[:,3:]
    if return_vel:
        _out['vel'][:len(data)] = data[:,3:6]
    if return_pid:
        _out['pid'][:len(data)] = data[:,6].view(dtype=np.uint64)
    
    retval = (_out,) if out is None else (len(data),)
    
    if return_header:
        retval += (header,)
        
    if len(retval) == 1:
        return retval[0]
    return retval
    
    
def read_rvzel(fn, return_vel=True, return_zel=False, return_pid=False, zspace=False, add_grid=False, boxsize=None, dtype=np.float32, out=None, return_header=False):
    """
    Reads RVzel format particle data.  This can be output from our
    zeldovich IC generator or Abacus.  zeldovich outputs don't have
    headers.
    
    Parameters
    ----------
    fn: string
        The filename
    return_vel: bool, optional
        Include the velocity in the returned data
    return_zel: bool, optional
        Include the 3 short ints that identify the initial Lagrangian lattice site
    return_pid: bool, optional
        Include the PID, as computed from the `zel' field with an automatically detected ppd
    zspace: bool, optional
        Apply redshift-space distortions
    add_grid: bool, optional
        Convert displacements to absolute positions by adding the grid locations
        determined by the `zel` field and `boxsize`
    boxsize: float, optional
        Necessary only if `add_grid`.
    dtype: np.dtype, optional
        The dtype to load the particle data into
    out: np.ndarray, optional
        Pre-allocated space to put the particles into
        
    Returns
    -------
    particles: ndarray of shape (npart,)
        dtype of result is determined by the parameters
    or, if `out` is given,
    NP: int
        The number of particles loaded into `out`

    and optionally,
    
    header: InputFile
        If `return_header` and a header is found, return parsed InputFile
    """
    
    rvzel_dt = np.dtype([("zel", np.uint16, 3), ("r", np.float32, 3), ("v", np.float32, 3)], align=True)
    with open(fn, 'rb') as fp:
        header = skip_header(fp)
        raw = np.fromfile(fp, dtype=rvzel_dt)
    
    if out is None:
        return_dt = output_dtype(return_vel=return_vel, return_zel=return_zel, return_pid=return_pid, dtype=dtype)
        particles = np.empty(raw.shape, dtype=return_dt)
    else:
        particles = out

    if len(raw) > 0:
        if add_grid or return_pid:
            ppd = np.array([raw['zel'].max() + 1], dtype=np.uint64)  # necessary for numpy to not truncate the result
        # We are only guaranteed to have a whole number of planes from the zeldovich code, but this might be an Abacus output
        #if add_grid or return_pid:
        #    #make sure we have a whole number of planes; otherwise, ppd might be wrong
        #    planes = len(raw) / ppd**2
        #    assert planes*ppd**2 == len(raw), "ppd {} detected from zel field, but this implies {} particle planes which does not match {} particles".format(ppd, planes, len(raw))
        
        particles['pos'][:len(raw)] = raw['r']
        if add_grid:
            grid = (1.*raw['zel'] / ppd - 0.5)*boxsize
            particles['pos'][:len(raw)] += grid
        if zspace:
            particles['pos'][:len(raw)] += raw['v']
        if return_vel:
            particles['vel'][:len(raw)] = raw['v']
        if return_zel:
            particles['zel'][:len(raw)] = raw['zel']
        if return_pid:
            particles['pid'][:len(raw)] = raw['zel'][:,2] + ppd*raw['zel'][:,1] + ppd**2*raw['zel'][:,0]

    retval = (particles,) if out is None else (len(raw),)
    
    if return_header:
        retval += (header,)
        
    if len(retval) == 1:
        return retval[0]

    return retval
    
    
def read_state(fn, make_global=True, dtype=np.float32, dtype_on_disk=np.float32, return_pid='auto', return_aux=False, return_vel='auto', return_pos='auto', return_header=False, out=None):
    """
    Read an Abacus position or velocity state file (like 'read/position_0000').
    
    The positions will be returned in whatever units they were stored in,
    which is typically zero-centered unit box.
    
    `cpd` is inferred from the cellinfo size, and `slab` is inferred from the filename.
    
    Parameters
    ----------
    fn: str
        Filename to read
    make_global: bool or 'auto', optional
        Whether to convert cell-offset positions into global positions by reading cellinfo
        Default: True`
    dtype_on_disk: np.dtype, optional
        The data type of the positions and cellinfo floats
    dtype: np.dtype, optional
        The output data type
    return_pid: bool, optional
        Return the PID as loaded from the `auxillary_*` slabs
    return_header, bool, optional
        Return the `state` file
        
    Returns
    -------
    data: ndarray of length (npart,)
        The particle data.  Positions are in data['pos']; velocities are in data['vel'] (if `return_vel`),
        and PIDs are in data['pid'] (if `return_pid`).
    or,
    NP: int
        If `out` is given, returns the number of rows read into `out`.
    """
    basename = path.basename(fn)
    if return_pos == 'auto':
        return_pos = 'position' in basename
    if return_vel == 'auto':
        return_vel = 'velocity' in basename
    if return_pid == 'auto':
        return_pid = 'auxillary' in basename 
    if make_global == 'auto':
        make_global = 'position' in basename
        
    # cellinfo dtype
    ci_dtype = np.dtype([('startindex',np.uint64),
                         ('count',np.int32), ('active',np.int32),
                         ('mean_square_velocity', dtype_on_disk), ('max_component_velocity', dtype_on_disk), ('max_component_acceleration', dtype_on_disk)],
                         align=True)
    pos_fn = re.sub(r'\w+(?=_\d{4}$)', r'position', fn)
    vel_fn = re.sub(r'\w+(?=_\d{4}$)', r'velocity', fn)
    ci_fn = re.sub(r'\w+(?=_\d{4}$)', r'cellinfo', fn)
    aux_fn = re.sub(r'\w+(?=_\d{4}$)', r'auxillary', fn)
    slab = int(re.search(r'\d{4}$', fn).group(0))

    if return_header:
        state = InputFile(path.join(path.dirname(fn), 'state'))
    
    # Count the particles to read in
    fsize = path.getsize(pos_fn)
    _format = {np.float32:'state', np.float64:'state64'}[dtype_on_disk]
    NP = get_np_from_fsize(fsize, format=_format)
    
    # Use the given buffer
    if out is not None:
        particles = out
    else:  # or allocate one
        ndt = output_dtype(return_vel=return_vel, return_pid=return_pid, dtype=dtype)
        particles = np.empty(alloc_NP, dtype=ndt)
    del fn
    
    if return_pos:
        particles['pos'][:NP] = np.fromfile(pos_fn, dtype=(dtype_on_disk, 3))
        if make_global:
            cellinfo = np.fromfile(ci_fn, dtype=ci_dtype)
            assert(cellinfo['count'].sum() == NP)
            assert np.all(np.cumsum(cellinfo['count'])[:-1] == cellinfo['startindex'][1:])
            cpd = int(np.round(np.sqrt(len(cellinfo))))
            assert cpd**2 == len(cellinfo)
            cellid = np.arange(len(cellinfo))
            centers = np.empty((cpd**2,3), dtype=dtype)

            centers[:,0] = 1.*slab
            centers[:,1] = cellid//cpd
            centers[:,2] = cellid%cpd
            cpdhalf = (cpd-1)/2.
            centers = (centers - cpdhalf)/cpd
            
            # wastes some space, but should be okay
            particles['pos'][:NP] += np.repeat(centers, cellinfo['count'], axis=0)
                
    if return_vel:
        particles['vel'][:NP] = np.fromfile(vel_fn, dtype=(dtype_on_disk, 3))
    if return_pid:
        particles['pid'][:NP] = np.fromfile(aux_fn, dtype=np.uint64)
        particles['pid'][:NP] &= 0xffffffffff  # PID aux bitmask (lower 40 bits)
    if return_aux:
        particles['aux'][:NP] = np.fromfile(aux_fn, dtype=np.uint64)
            
    if out is not None:
        particles = NP
        
    if return_header:
        return particles, state
    return particles


def read_pack14_lite(fn, return_vel=True, return_pid=False, return_header=False, dtype=np.float32, out=None):
    """
    Read particle data from a file in pack14 format using a Numba-based library instead of
    a C library, as `read_pack14()` uses.

    N.B.: expect this function to be many times slower than read_pack14() (2-3x slower in my tests).
          Someday Numba will catch up...

    Parameters 
    ----------
    fn: str
        The filename to read
    return_vel: bool, optional
        Return velocities along with other data
    return_pid: bool, optional
        Return particle IDs along with other data
    return_header: bool, optional
        If the pack14 file has an ASCII header, return it as a second return value.
    dtype: data-type, optional
        Either np.float32 or np.float64.  Determines the data type the particle data is loaded into
    out: ndarray, optional
        A pre-allocated array into which the particles will be directly loaded.
        
    Returns
    -------
    data: ndarray of length (npart,)
        The particle data.  Positions are in data['pos']; velocities are in data['vel'] (if `return_vel`),
        and PIDs are in data['pid'] (if `return_pid`).
        
    or,
        
    NP: int
        If `out` is given, returns the number of rows read into `out`.
        
    and optionally,
    
    header: InputFile
        If `return_header` and a header is found, return parsed InputFile
    """

    with open(fn, 'rb') as fp:
        header = skip_header(fp)
        if header:
            header = InputFile(str_source=header)
        try:
            raw = np.fromfile(fp, dtype=np.uint8).reshape(-1,14)
        except ValueError as e:
            print('Could not reshape array. Bad file size?')
            raise e

    # Use the given buffer
    if out is not None:
        _out = out
        outdt = out.dtype
    # or allocate one
    else:
        fsize = path.getsize(fn)
        alloc_NP = get_np_from_fsize(fsize, format='pack14')
        outdt = output_dtype(return_vel=return_vel, return_pid=return_pid, dtype=dtype)
        _out = np.empty(alloc_NP, dtype=outdt)
    
    # Now merge the pos and vel fields, if present. It's a little faster in Numba.
    dtlist = []
    if 'pid' in outdt.fields:
        dtlist += [('pid',_out['pid'].dtype)]
        # Numba doesn't like None, so use len 0 arrays instead 
        pid = _out['pid']
    else:
        pid = np.empty(0, dtype=np.uint64)
    
    if 'vel' in outdt.fields:
        dtlist += [('posvel',np.float32,6)]
    else:
        dtlist += [('posvel',np.float32,3)]
        
    merged_dt = np.dtype(dtlist, align=True)
    _out = _out.view(dtype=merged_dt)
    posvel = _out['posvel']
    
    npread = p14lite._read_pack14(raw, posvel, pid)
    _out = _out[:npread]  # shrink the buffer to the real size
    _out = _out.view(dtype=outdt)
    
    retval = (npread,) if out is not None else (_out,)
    if return_header:
        retval += (header,)
    
    if len(retval) == 1:
        return retval[0]
    return retval

##############################
# End list of reader functions
##############################

#############
# Begin utils
#############
    
def skip_header(fp, max_tries=10):
    """
    Some of our files are written with an ASCII header before
    the data section begins.  Given a file pointer, this fast-
    forwards past the header, also returning it as a string if
    it exists.
    
    Parameters
    ----------
    fp: file
        The file pointer to fast-forward
    max_tries: int, optional
        The header is guaranteed to be some multiple of
        4096 bytes.  But if there's no header, then the
        routine will keep jumping forward by 4096 bytes
        looking for the end-of-header, until it reads
        the whole file.  We may want to truncate that
        search early if we suspect there is no header
        (or that it will be small, since our usual
        header is < 4096 bytes).  Default: 10
        
    Returns
    -------
    header: str
        The skipped header, if one is found.
    """
    fsize = path.getsize(fp.name)
    fpstart = fp.tell()
    
    # Specify that the header must end on a 4096 byte boundary
    alignment = 4096
    sentinel = b'\n'  # header ends with ^B\n
    ss = len(sentinel)  # sentinel size = 2

    nomax = type(max_tries) is not int or max_tries < 0
    
    i = 0
    while nomax or i < max_tries:
        fp.seek(alignment - ss, 1)  # jump forward by 4094
        headerend = fp.read(ss)
        if len(headerend) < ss:
            # EoF reached without finding a header.  Reset fp and return.
            fp.seek(fpstart)
            return
        if headerend == sentinel:
            break  # success!
        i += 1
    else:
        # We didn't find the header; reset the pointer to where it was
        fp.seek(fpstart)
        return
    
    # Found a header.  Jump back to beginning and read    
    headerend = fp.tell()
    fp.seek(fpstart)
    headersize = headerend - fp.tell()
    header = fp.read(headersize)
    assert len(header) == headersize, "Header unexpectedly small!"
    header = header[:-2]  # trim the last two bytes
    return str(header)


# These defaults have to be consistent with the reader function defaults
def output_dtype(return_vel=True, return_pid=False, return_zel=False, return_aux=False, dtype=np.float32, **kwargs):
    """
    Construct the dtype of the output array.
    """
    ndt_list = []
    if return_pid:
        ndt_list += [('pid', np.uint64)]
    if return_aux:
        ndt_list += [('aux', np.uint64)]
    ndt_list += [('pos', dtype, 3)]
    if return_vel:
        ndt_list += [('vel', dtype, 3)]
    if return_zel:
        ndt_list += [("zel", np.uint16, 3)]
    ndt = np.dtype(ndt_list, align=True)

    return ndt

# Set up a few library utils
try:
    from . import abacus
    ralib = ct.cdll.LoadLibrary(path.join(abacus.abacuspath, "clibs", "libreadabacus.so"))

    # Set up the arguments and return type for the library functions
    # TODO: switch our C library to use CFFI
    for f in (ralib.read_pack14, ralib.read_pack14f):
        f.restype = ct.c_uint64
        f.argtypes = (asciistring_arg, ct.c_size_t, ct.c_int, ct.c_int, ct.c_int, ct.c_int, ndarray_arg)
except (OSError, ImportError):
    pass  # no pack14 library found

# An upper limit on the number of particles in a file, based on its size on disk
get_np_from_fsize = lambda fsize, format, downsample=1: int(np.ceil(float(fsize)/psize_on_disk[format]/downsample**3.))
psize_on_disk = {'pack14': 14, 'pack14_lite':14, 'rvdouble': 6*8, 'state64':3*8, 'state':3*4, 'rvzel':32, 'rvdoubletag':7*8}
