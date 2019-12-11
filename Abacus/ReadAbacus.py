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
- automatic globbing patterns for a given file format
- standardize return signatures
"""

import os
import os.path as ppath
from os.path import join as pjoin, basename, dirname
from glob import glob
import ctypes as ct
import re
import threading
import queue
from warnings import warn

import numpy as np
import numba

import ctypes

from .InputFile import InputFile
from .Tools import ndarray_arg, asciistring_arg
from .Tools import wrap_zero_centered, wrap_zero_origin


def read(*args, **kwargs):
    """
    A convenience function to read particle data from
    file type `format`.  Simply calls
    the appropriate `read_[format]()` function.
    `format` may also be a function that accepts a filename
    and returns the particle data.
    For usage, see the docstring of the reader
    function for the file format you are using.
    """
    format = kwargs.pop('format')
    if type(format) == str:
        format = format.lower()

    units = kwargs.pop('units', None)
    box_on_disk = kwargs.pop('box_on_disk', None)
    
    # Readers should return None if a header isn't found
    original_return_header = kwargs.get('return_header')
    kwargs['return_header'] = True

    try:
        ret = reader_functions[format](*args, **kwargs)
    except KeyError:
        # Try calling format as a function
        ret = format(*args, **kwargs)

    if type(ret) is tuple:  # unambiguous way to unpack this?
        data, header = ret
    else:
        data = ret

    # todo: compare str and float
    if units != None:
        if box_on_disk is None:
            fn = args[0]
            box_on_disk = get_box_on_disk(fn, format)
        if units != box_on_disk:
            try:
                box = kwargs['boxsize']
            except:
                box = header['BoxSize']
            data['pos'] *= eval(str(units))/eval(str(box_on_disk))
            # vel conversion?

    if original_return_header:
        return data, header
    return data


def from_dir(dir, pattern=None, key=None, **kwargs):
    """
    Read all files from `dir` that match `pattern`.
    
    Parameters
    ----------
    dir: str
        The directory to read files from
    pattern: str or None, optional
        A bash globbing pattern to find all the files to read.
        If None, use the default pattern for the given format.
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
    if pattern is None:
        if is_ic_path(path):
            pattern = 'ic_*'
        else:
            format = kwargs.get('format')
            pattern = get_file_pattern(format)

    _key = (lambda k: key(ppath.basename(k))) if key else None
    files = []

    if type(pattern) is str:
        pattern = (pattern,)
        
    for p in pattern:
        files += glob(pjoin(dir, p))
    files = sorted(files, key=_key)

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
        TODO: if format is a function, need to guess alloc_NP
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
    total_fsize = sum(ppath.getsize(fn) for fn in files)
    
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
    alloc_NP = get_np_from_fsize(total_fsize, format=_format, downsample=kwargs.get('downsample'))
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
            print('pa',  particle_arrays[field][:start])

            particle_arrays[field] = particle_arrays[field][:start]

        if return_header:
            return particle_arrays, header  # return header associated with last file
        return particle_arrays
    else:
        print('p', particles, particles[:start])

        particles = particles[:start]
        if return_header:
            return particles, header  # return header associated with last file
    
        return particles


def AsyncReader(path, readahead=1, chunksize=1, key=None, verbose=False, return_fn=False, **kwargs):
    '''
    This is a generator that reads files in a separate thread in the background,
    yielding them one at a time as they become available.
    The intended usage is:
    ```
    for chunk_of_particles in AsyncReader('/slice/directory'):
        process(chunk_of_particles)
    ```
    The next file will be read in the background while `process` is running.  For code
    that only needs a small chunk of the simulation in memory at a time, this will run
    faster and use less memory.
    The amount of readahead is configurable to ensure one doesn't run out of memory.
    In order to overlap compute and IO, the IO code must release the GIL.  Since np.fromfile()
    does that, most of our readers do that as well.  pack14 doesn't use fromfile, but ctypes
    also releases the GIL.
    Parameters
    ----------
    path: str or list of str
        A directory or file pattern, or a list of files.
    readahead: int, optional
        The number of files to read ahead of the last-yielded file.
        Default of 1 is fine when compute is slower than IO; a value
        too big might run out of memory.
    chunksize: int, optional
        The number of files to read at a time (into a contiguous array).
        Default: 1.
    key: callable, optional
        A filename sorting key.  Default: None
    verbose: bool, optional
        Print status.  Default: False
    return_fn: bool, optional
        Yield the filename alongside the data.  Default: False
    **kwargs: dict, optional
        Extra arguments to be passed to `read` or `read_many`.
    '''

    assert chunksize == 1, "chunksize > 1 not yet implemented"

    _key = (lambda k: key(ppath.basename(k))) if key else None


    # First, determine the file names to read

    if type(path) is str:
        if ppath.isdir(path):
            if is_ic_path(path):
                pattern = 'ic_*'
            else:
                pattern = get_file_pattern(kwargs.get('format'))
            files = sorted(glob(pjoin(path, pattern)), key=_key)
        elif ppath.isfile(path):
            files = [path]
        else:
            # Not a dir and not a file; interpret as a glob pattern
            files = sorted(glob(path), key=_key)
    else:
        # Already a list?
        files = path

    assert files, f'No files found matching "{path}"!'

    NP = 0
    if readahead < 0:
        readahead = len(files)
    file_queue = queue.Queue(maxsize=readahead+1)

    def reader_loop():
        nonlocal NP
        reader_kwargs = kwargs.copy()
        format = reader_kwargs.pop('format')

        Nfn = len(files)

        # Read and bin the particles
        for i,filename in enumerate(files):
            if verbose:
                print(f'Reading {i+1}/{Nfn} ' +'("{}")... '.format(basename(files[i])), end='', flush=True)
            data = read(filename, format=format, **reader_kwargs)
            if verbose:
                print('done.', flush=True)
            NP += len(data)

            file_queue.put((data,filename))  # blocks until free slot available
        file_queue.put(None)  # signal termination

    io_thread = threading.Thread(target=reader_loop)
    io_thread.start()
    
    while True:
        data = file_queue.get()
        if data is None:
            break
        data,fn = data
        if return_fn:
            yield data, fn
        else:
            yield data

    io_thread.join()
    assert file_queue.empty()

    

################################
# Begin list of reader functions
################################

def read_pack14(fn, ramdisk=False, return_vel=True, zspace=False, return_pid=False, return_header=False, dtype=np.float32, boxsize=None, downsample=None, out=None):
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
    boxsize: optional
        Ignored; included for compatibility
    downsample: float, optional
        The downsample fraction.  Downsampling is performed using the same PID hash as Abacus proper.
        Default of None means no downsampling.
        
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

    if zspace:
        return_vel = True

    try:
        ralib
    except NameError:
        raise RuntimeError("pack14 C library was not found. Try building Abacus with 'make analysis'? Or use the slower `read_pack14_lite()` function.")
    readers = {np.float32: ralib.read_pack14f,
               np.float64: ralib.read_pack14 }
    dtype = np.dtype(dtype).type
    assert dtype in readers, dtype
    

    with open(fn, 'rb') as fp:
        header = skip_header(fp)
        if header:
            header = InputFile(str_source=header)
        offset = fp.tell()
    
    # Use the given buffer
    if out is not None:
        _out = out
    else:  # or allocate one
        fsize = ppath.getsize(fn)
        alloc_NP = get_np_from_fsize(fsize, format='pack14', downsample=downsample)
        ndt = output_dtype(return_vel=return_vel, return_pid=return_pid, dtype=dtype)
        _out = np.empty(alloc_NP, dtype=ndt)

    try:
        if downsample > 1:
            warn(f'Downsample factor {downsample} is greater than 1!  A fraction less than 1 is expected.')
    except TypeError: 
        if downsample is None:
            downsample = 1.1  # any number larger than 1 will take all particles
    
    
    NP = readers[dtype](fn, offset, ramdisk, return_vel, False, return_pid, downsample, _out.view(dtype=dtype))


    # shrink the buffer to the real size
    if out is None:
        _out.resize(NP, refcheck=False)
    elif not return_vel:
        _out = _out[:NP]
    else:
        _out = _out[:2*NP] 

    if return_vel:
        data = np.array([list(triplet) for line in _out for triplet in line])
        pos = data[::2]
        vel = data[1::2]
        if zspace:
            vel = vel  * (header['VelZSpace_to_kms']/header['VelZSpace_to_Canonical']); 
            pos += vel # add velocities. Don't do box wrap. 
        _out = np.array(pos)

    retval = (NP,) if out is not None else (_out,)
    if return_header:
        retval += (header,)
    
    if len(retval) == 1:
        return retval[0]
    return retval

def read_pack9(fn, ramdisk=False, return_vel=True, zspace=False, return_pid=False, return_header=False, dtype=np.float32, boxsize=None, downsample=None, out=None):
    """
    Read particle data from a file in pack9 format.
    
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
    boxsize: optional
        Ignored; included for compatibility
    downsample: float, optional
        The downsample fraction.  Downsampling is performed using the same PID hash as Abacus proper.
        Default of None means no downsampling.
        
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
    if zspace: 
        return_vel = True
    try:
        ralib
    except NameError:
        raise RuntimeError("pack9 C library was not found. Try building Abacus with 'make analysis'?")
    readers = {np.float32: ralib.read_pack9f,
               np.float64: ralib.read_pack9 }
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
        fsize = ppath.getsize(fn)
        alloc_NP = get_np_from_fsize(fsize, format='pack9', downsample=downsample)
        ndt = output_dtype(return_vel=return_vel, return_pid=return_pid, dtype=dtype)
        _out = np.empty(alloc_NP, dtype=ndt)

    try:
        if downsample > 1:
            warn(f'Downsample factor {downsample} is greater than 1!  A fraction less than 1 is expected.')
    except TypeError: 
        if downsample is None:
            downsample = 1.1  # any number larger than 1 will take all particles

    readers[dtype].argtypes = [ctypes.c_char_p, ctypes.c_size_t, ctypes.c_int, ctypes.c_int, ctypes.c_int, ctypes.c_int, ctypes.c_double, ctypes.c_double, ctypes.c_void_p] 
    NP = readers[dtype](bytes(fn, encoding='utf-8'), offset, ramdisk, return_vel, False, return_pid, downsample, 0, _out.view(dtype=dtype).ctypes.data_as(ctypes.POINTER(ctypes.c_void_p)))

    # shrink the buffer to the real size
    if out is None:
        _out.resize(NP, refcheck=False)
    elif not return_vel:
        _out = _out[:NP]
    else:
        _out = _out[:2*NP] 

    if return_vel:
        data = np.array([list(triplet) for line in _out for triplet in line])
        pos = data[::2]
        vel = data[1::2]
        if zspace:
            vel = vel  * (header['VelZSpace_to_kms']/header['VelZSpace_to_Canonical']); 

            pos += vel # add velocities. Don't do box wrap. 
        _out = np.array(pos)

    
    retval = (NP,) if out is not None else (_out,)
    if return_header:
        retval += (header,)
    
    if len(retval) == 1:
        return retval[0]
    return retval    


def rvint_unpack(data):
    velscale = 6000.0/2048.0
    posscale = 2.0**(-32.0)

    iv = (data&0xfff).astype(int)
    ix = (data-iv).astype(int)
    pos = posscale*ix 
    vel = velscale*(iv - 2048)
  
    return pos, vel


def read_rvint(fn, return_vel = True, return_pid=False, zspace=False, dtype=np.float32, out=None, return_header=False, double=False, tag=False,  downsample=None):

    disk_base_dt = np.int32 
    disk_dt = [('pv',disk_base_dt,3)]

    print("Reading rvint from ", fn)

    with open(fn, 'rb') as fp:
        header = skip_header(fp)
        if header:
            header = InputFile(str_source=header)
        data = np.fromfile(fp, dtype=disk_dt)
    
    data = np.array([list(triplet) for line in data for triplet in line])

    state = InputFile(pjoin(ppath.dirname(fn), 'header'))

    # Use the given buffer
    if out is not None:
        _out = out  

    pos = np.zeros([len(data), 3]) 
    vel = np.zeros([len(data), 3])

    pos, vel = rvint_unpack(data) 

    _out['pos'][:len(data)] = pos

    # vel = vel  / (state['VelZSpace_to_kms']/state['VelZSpace_to_Canonical']); 

    if zspace:
        print(pos, vel) 
        _out['pos'][:len(data)] += vel
    if return_vel:
        _out['vel'][:len(data)] = vel

    retval = (_out,) if out is None else (len(data),)
    if return_header:
        retval += (state,)

    if len(retval) == 1:
        return retval[0]

    return retval

    
def read_rvtag(*args,**kwargs):
    return read_rv(*args, tag=True, **kwargs)

def read_rvdouble(*args,**kwargs):
    return read_rv(*args, double=True, dtype=np.float64, **kwargs)

def read_rvdoubletag(*args,**kwargs):
    return read_rv(*args, double=True, tag=True, dtype=np.float64, **kwargs)

def read_rv(fn, return_vel=True, return_pid=False, zspace=False, dtype=np.float32, out=None, return_header=False, double=False, tag=False):
    """
    Read particle data from a file in rvdoubletag format.
    
    Usage is the same as `read_pack14`.
    Parameters
    ----------
    double: bool, optional
        Whether the format on disk is RVdoubleTag or just RVTag.
        Default: False
    tag: bool, optional
        Whether the format on disk is RVTag or just RV.
        Default: False
    """
    if return_pid:
        assert tag

    disk_base_dt = np.float64 if double else np.float32
    disk_dt = [('pos',disk_base_dt,3), ('vel',disk_base_dt,3)]
    if tag:
        disk_dt += [('pid',np.uint64)]

    with open(fn, 'rb') as fp:
        header = skip_header(fp)
        if header:
            header = InputFile(str_source=header)
        data = np.fromfile(fp, dtype=disk_dt)
    
    if out is not None:
        _out = out
    else:
        # special case: can we return the data raw?
        if tag == return_pid and return_vel and dtype == disk_base_dt and not zspace:
            if return_header:
                return data, header
            else:
                return data
        else:
            ndt = output_dtype(return_vel=return_vel, return_pid=return_pid, dtype=dtype)
            _out = np.empty(len(data), dtype=ndt)
    
    _out['pos'][:len(data)] = data['pos']

    if zspace:
        _out['pos'][:len(data)] += data['vel']
    if return_vel:
        _out['vel'][:len(data)] = data['vel']
    if return_pid:
        _out['pid'][:len(data)] = data['pid']
    
    retval = (_out,) if out is None else (len(data),)
    
    if return_header:
        retval += (header,)
        
    if len(retval) == 1:
        return retval[0]
    return retval
    

def read_rvdoublezel(*args, **kwargs):
    return read_rvzel(*args, double=True, **kwargs)
    
def read_rvzel(fn, return_vel=True, return_zel=False, return_pid=False, zspace=False, add_grid=False, boxsize=None, dtype=np.float32, out=None, return_header=False, double=False):
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
    double: bool, optional
        Whether the input file format is RVdoubleZel
        
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
    
    dtype_on_disk = np.float64 if double else np.float32

    rvzel_dt = np.dtype([("zel", np.uint16, 3), ("r", dtype_on_disk, 3), ("v", dtype_on_disk, 3)], align=True)
    with open(fn, 'rb') as fp:
        header = skip_header(fp)
        raw = np.fromfile(fp, dtype=rvzel_dt)

    if header:
        header = InputFile(str_source=header)
    else:
        # Look for a file named 'header'
        try:
            header = InputFile(fn=pjoin(ppath.dirname(fn), 'header'))
        except:
            header = None

    if header and not boxsize:
        output_type = header.get('OutputType', 'ic')
        if output_type == 'TimeSlice':
            boxsize = 1.
        elif output_type == 'ic':
            boxsize = header['BoxSize']
        else:
            raise ValueError(output_type)
    
    if out is None:
        return_dt = output_dtype(return_vel=return_vel, return_zel=return_zel, return_pid=return_pid, dtype=dtype)
        particles = np.empty(raw.shape, dtype=return_dt)
    else:
        particles = out

    if len(raw) > 0:
        if add_grid or return_pid:
            if header:
                ppd = np.array([np.round(header['NP']**(1./3))], dtype=np.int64)
            else:
                # This will only work for ICs, where we don't have a header
                ppd = np.array([raw['zel'].max() + 1], dtype=np.int64)  # necessary for numpy to not truncate the result

        # We are only guaranteed to have a whole number of planes from the zeldovich code, but this might be an Abacus output
        #if add_grid or return_pid:
        #    #make sure we have a whole number of planes; otherwise, ppd might be wrong
        #    planes = len(raw) / ppd**2
        #    assert planes*ppd**2 == len(raw), "ppd {} detected from zel field, but this implies {} particle planes which does not match {} particles".format(ppd, planes, len(raw))
        
        particles['pos'][:len(raw)] = raw['r']
        if add_grid:
            if not boxsize:
                raise ValueError("Need to specify boxsize if using add_grid (or a header must be present from which the box size can be read)")
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
    basename = ppath.basename(fn)
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
        state = InputFile(pjoin(ppath.dirname(fn), 'state'))
    
    # Count the particles to read in
    fsize = ppath.getsize(pos_fn)
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
        fsize = ppath.getsize(fn)
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


def read_gadget(fn, **kwargs):
    '''
    Bare-bones Gadget reader using pynbody.
    '''
    dtype = kwargs['dtype']

    import pynbody
    f = pynbody.load(fn)
    data = f['pos'].astype(dtype=dtype)
    #data = np.array(f['pos'] - box_on_disk/2., dtype=('pos',dtype,3))  # Shift to zero-centered

    return {'pos':data}


def read_desi_hdf5(fn, **kwargs):
    '''
    HDF5 format used in the DESI CosmoSim WG.
    '''
    import h5py

    data = {}
    with h5py.File(fn, 'r') as fp:
        data['pos'] = fp['/Matter/Position'][:]

    return data


##############################
# End list of reader functions
##############################

#############
# Begin utils
#############
    
def skip_header(fp, max_tries=10, encoding='utf-8'):
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
    fsize = ppath.getsize(fp.name)
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
    return header.decode(encoding)


# These defaults have to be consistent with the reader function defaults
def output_dtype(return_vel=True, return_pid=False, return_zel=False, return_aux=False, dtype=np.float32, **kwargs):
    """
    Construct the dtype of the output array.
    """
    ndt_list = []
    if return_pid:
        ndt_list += [('pid', np.int64)]
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
    ralib = ct.cdll.LoadLibrary(pjoin(abacus.abacuspath, "clibs", "libreadabacus.so"))

    # Set up the arguments and return type for the library functions
    # TODO: switch our C library to use CFFI
    for f in (ralib.read_pack14, ralib.read_pack14f):
        f.restype = ct.c_uint64
        f.argtypes = (asciistring_arg, ct.c_size_t, ct.c_int, ct.c_int, ct.c_int, ct.c_int, ct.c_double, ndarray_arg)
except (OSError, ImportError):
    raise
    pass  # no pack14 library found


# An UPPER LIMIT on the number of particles in a file, based on its size on disk and downsampling rate
def get_np_from_fsize(fsize, format, downsample=None):

    ds = 1. if downsample is None else downsample

    ds = max(min(1.,ds),0.)  # clip to [0,1]

    N = float(fsize)/psize_on_disk[format]*ds

    # No downsampling uncertainty
    if downsample is None:
        return int(np.ceil(N))

    # Add 6-sigma downsampling uncertainty
    N += 6*np.sqrt(N)

    return int(np.ceil(N))

psize_on_disk = {'pack14': 14, 'pack14_lite':14, 'pack9': 9, 'rvint': 3*4, 'rvdouble': 6*8, 'state64':3*8, 'state':3*4, 'rvzel':32, 'rvtag':32, 'rvdoubletag':7*8}
reader_functions = {'pack14':read_pack14, 'pack9':read_pack9, 'rvint': read_rvint, 'rvdouble':read_rvdouble,

                    'rvzel':read_rvzel, 'state':read_state,
                    'rvdoublezel':read_rvdoublezel,
                    'rvdoubletag':read_rvdoubletag,
                    'rvtag':read_rvtag, 'rv':read_rvtag,
                    'pack14_lite':read_pack14_lite,
                    'gadget':read_gadget,
                    'desi_hdf5':read_desi_hdf5}
default_file_patterns = {'pack14':('*.dat',),
                         'pack9':('*L0_pack9.dat', '*field_pack9.dat'),
                         'rvint' : ('*rv_A*', '*rv_B*'),
                         'state':('position_*',),
                         'desi_hdf5':('*.hdf5',),
                         'gadget':('*.*',)
                         }
fallback_file_pattern = ('*.dat',)
                    
default_box_on_disk = {'desi_hdf5':'box',
                'pack14':1.,
                'rvtag':'box',
                'gadget':'box',
}

def get_file_pattern(format):
    try:
        return default_file_patterns[format]
    except KeyError:
        return fallback_file_pattern

def is_ic_path(path):
    '''
    Our initial conditions files are stored on disk in a BoxSize box,
    not unit box, so it's convenient to detect if this is one of our
    standard IC paths.
    '''
    # Does the directory name contain 'ic'?
    abspattern = os.path.abspath(path)
    abspattern = abspattern.split(os.sep)
    isic = re.search(r'\bic(\b|(?=_))', os.sep.join(abspattern))

    return isic

def get_box_on_disk(fn, format):
    if format == 'gadget':
        import pynbody
        f = pynbody.load(fn)
        box_on_disk = float(f.properties['boxsize'])
        return box_on_disk

    isic = is_ic_path(fn)

    # We assume IC data is stored with a BoxSize box, not unit box
    box_on_disk = 'box' if isic else default_box_on_disk[format]

    return box_on_disk
