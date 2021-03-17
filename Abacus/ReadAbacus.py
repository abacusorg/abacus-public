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
The alternative is to use the slower `read_pack14_lite()` function.
Todo:
- Unit conversion options
- native ctypes float3* instead of void* (or CFFI?)
- inline particle downsampling
- standardize return signatures
"""

import os
from os.path import join as pjoin, basename, dirname, isdir, isfile, getsize
from glob import glob
import ctypes
import ctypes as ct
import re
import threading
import queue
from warnings import warn
import gc

from astropy.utils import iers
iers.conf.auto_download = False
import numpy as np
import numba
import asdf
from astropy.table import Table

from .InputFile import InputFile
from .Tools import ndarray_arg, asciistring_arg
from .Tools import wrap_zero_centered, wrap_zero_origin
from .Tools import ContextTimer
from .abacus_halo_catalog import unpack_rvint

# blosc is the decompression underlying our ASDF files
# No gain is likely beyond 4 threads
BLOSC_THREADS = 4

def read(*args, **kwargs):
    """
    A convenience function to read one file's particle data
    from file type `format`.  Simply calls
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

    if callable(format):
        reader_func = format
    else:
        try:
            reader_func = reader_functions[format]
        except KeyError:
            raise ValueError(f'Unknown format "{format}"')

    ret = reader_func(*args, **kwargs)

    if type(ret) is tuple:  # unambiguous way to unpack this?
        data, header = ret
    else:
        data = ret

    if hasattr(data, 'meta'):
        header = data.meta

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
            with ContextTimer('read() box scale'):
                data['pos'] *= eval(str(units))/eval(str(box_on_disk))
            # vel conversion?

    if original_return_header:
        return data, header
    return data


def from_dir(dir, pattern=None, key=None, **kwargs):
    """
    Read all files from `dir` that match `pattern`.  If a pattern
    is not given, one will be inferred from `format`.
    
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

    if not isdir(dir):
        raise ValueError(f'Path "{dir}" was passed to from_dir() but it is not a directory')

    files = get_files_from_path(dir, pattern=pattern, key=key, **kwargs)

    return read_many(files, **kwargs)


def read_many(files, format='pack14', **kwargs):
    """
    Read a list of files into a contiguous array.

    Parameters
    ----------
    files: list of str
        The list of file names
    format: str
        The file format. Specifies the reader function.
        TODO: if format is a function, need to guess alloc_NP
    kwargs: dict, optional
        Additional args to pass to the reader function.

    Returns
    -------
    particles: Astropy table
        The concatenated particle table
    """

    # Allocate enough space to hold the concatenated particle array
    if not files:
        raise ValueError("No files passed to read_many()!")

    verbose = kwargs.get('verbose',False)
    
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
    alloc_NP = get_alloc_np(files, format=_format, downsample=kwargs.get('downsample'))

    particles = allocate_table(alloc_NP, **kwargs)
    
    return_header = kwargs.get('return_header', False)
    header = None

    tot_read_time = 0
    tot_unpack_time = 0
    start = 0
    for fn in files:
        offsetparticles = particles[start:]
        out = read(fn, format=format, out=offsetparticles, **kwargs)

        if hasattr(offsetparticles,'meta'):
            header = offsetparticles.meta
            NP = out
            #tot_read_time += offsetparticles.meta['read_time']
            #tot_unpack_time += offsetparticles.meta['unpack_time']
        elif return_header:
            NP, header = out
        else:
            NP = out
        start += NP

    if verbose:
        # TODO: rates
        print(f'Total ReadAbacus read time: {tot_read_time:.4g} sec')
        print(f'Total ReadAbacus unpack time: {tot_unpack_time:.4g} sec')

    # Shrink the array to the size that was actually read
    particles = particles[:start]
    if header:
        particles.meta.update(header)

    return particles


def AsyncReader(path, readahead=1, chunksize=1, key=None, verbose=False, return_fn=False,
                nreaders=1, **kwargs):
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
    
    Multiple IO threads may be used, in which case they'll all push to the same
    results queue.  So this is a multiple-producer, single consumer model.  Be
    aware that multiple threads will increase memory usage, because nreaders*(readahead+1)
    results may be sitting in the queue.

    Parameters
    ----------
    path: str or list of str
        A directory or file pattern, or a list of files.
    readahead: int, optional
        The number of files that each thread is allowed to read ahead.
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
    nreaders: int, optional
        Number of IO threads. Default: 1
    **kwargs: dict, optional
        Extra arguments to be passed to `read` or `read_many`.
    '''

    if chunksize != 1:
        raise NotImplementedError("chunksize > 1 not yet implemented")

    format = kwargs.pop('format')
    files = get_files_from_path(path, format=format)
    Nfn = len(files)

    NP = 0
    if readahead < 0:
        readahead = len(files)
        
    fn_queue = queue.Queue()
    data_queue = queue.Queue(maxsize=nreaders*(readahead+1))
    
    for i,fn in enumerate(files):
        fn_queue.put((i,fn))

    def reader_loop():
        reader_kwargs = kwargs.copy()

        # Read a file and push to queue
        while True:
            try:
                i,filename = fn_queue.get_nowait()
            except queue.Empty:
                break
            
            if verbose:
                print(f'Reading {i+1}/{Nfn}... ', end='', flush=True)
            data = read(filename, format=format, **reader_kwargs)
            if verbose:
                print('done.', flush=True)

            data_queue.put((data,filename))  # blocks until free slot available
            del data

    io_threads = [threading.Thread(target=reader_loop) for _ in range(nreaders)]
    for it in io_threads:
        it.start()
    
    tot_read_time, tot_unpack_time = 0, 0
    for _ in range(Nfn):
        data = data_queue.get()
        data,fn = data
        NP += len(data)

        if hasattr(data, 'meta'):
            tot_read_time += data.meta['read_time']
            tot_unpack_time += data.meta['unpack_time']
        if return_fn:
            yield data, fn
        else:
            yield data
        del data
        gc.collect()
    
    for it in io_threads:
        it.join()
    assert data_queue.empty()
    assert fn_queue.empty()

    if verbose:
        print(f'AsyncReader read {NP} particles')
        print(f'Total ReadAbacus read time: {tot_read_time:.4g} sec')
        print(f'Total ReadAbacus unpack time: {tot_unpack_time:.4g} sec')

    

################################
# Begin list of reader functions
################################

def read_packN(N, fn, return_pos=True, return_vel=True, zspace=False, return_pid=False, return_header=False, dtype=np.float32, boxsize=None, downsample=None, out=None):
    """
    Read particle data from a file in pack9 or pack14 format.
    
    Parameters 
    ----------
    N: int
        The packN format to read, 9 or 14
    fn: str
        The filename to read
    return_vel: bool, optional
        Return velocities along with other data
    zspace: bool, optional
        Apply redshift-space distortion to particle positions
    return_pid: bool, optional
        Return particle IDs along with other data
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
    data: astropy Table of length (npart,)
        The particle data.  Positions are in data['pos']; velocities are in data['vel'] (if `return_vel`),
        and PIDs are in data['pid'] (if `return_pid`).
        
    or,
        
    NP: int
        If `out` is given, returns the number of rows read into `out`.
    """

    if N not in (9, 14):
        raise ValueError("N must be 9 or 14!")

    try:
        ralib
    except NameError:
        raise RuntimeError("packN C library was not found. Try building Abacus with 'make analysis'? Or use the slower `read_pack14_lite()` function.")
    dtype = np.dtype(dtype).type
    reader = packN_readers[N][dtype]
    
    try:
        if downsample > 1:
            warn(f'Downsample factor {downsample} is greater than 1!  A fraction less than 1 is expected.')
    except TypeError: 
        if downsample is None:
            downsample = 1.1  # any number larger than 1 will take all particles

    with open(fn, 'rb') as fp, ContextTimer('Read', output=False) as timer:
        header = skip_header(fp)
        if header:
            header = InputFile(str_source=header)
        data = np.fromfile(fp, dtype=np.byte).reshape(-1,N)
    
    # Use the given buffer
    if out is not None:
        _out = out
    else:  # or allocate one
        alloc_NP = get_alloc_np(fn, format=f'pack{N}', downsample=downsample)
        _out = allocate_table(alloc_NP, return_vel=return_vel, return_pid=return_pid, dtype=dtype)
    _out.meta.update(header)

    posout = _out['pos'] if return_pos else None
    velout = _out['vel'] if return_vel else None
    pidout = _out['pid'] if return_pid else None
    
    nthread = 1
    with ContextTimer('Unpack',output-False) as unpack_timer:
        NP = reader(data, data.nbytes, nthread, zspace, downsample, posout, velout, pidout)

    _out.meta['read_time'] = timer.elapsed + _out.meta.get('read_time',0.)
    _out.meta['unpack_time'] = unpack_timer.elapsed + _out.meta.get('unpack_time',0.)

    # shrink the buffer to the real size
    if out is None:
        # TODO: can we call resize on each table column safely? are there any deep refs?
        #_out.resize(NP, refcheck=False)
        _out = _out[:NP]
    else:
        _out = _out[:NP]
    
    if out:
        return NP
    return _out


read_pack14 = lambda *args,**kwargs: read_packN(14, *args, **kwargs)
read_pack9 = lambda *args,**kwargs: read_packN(9, *args, **kwargs)
read_pack14.__doc__ = read_pack9.__doc__ = read_packN.__doc__


def read_rvint(fn, return_vel = True, return_pid=False, zspace=False, dtype=np.float32, out=None, return_header=False, double=False, tag=False,  downsample=None):
    if return_pid:
        raise NotImplementedError  # comes from different file

    disk_dt = np.int32

    with open(fn, 'rb') as fp:
        header = skip_header(fp)
        if header:
            header = InputFile(str_source=header)
        data = np.fromfile(fp, dtype=disk_dt).reshape(-1,3)

    state = InputFile(pjoin(dirname(fn), 'header'))

    # Use the given buffer
    if out is not None:
        _out = out
    else:
        alloc_NP = get_alloc_np(fn, format='rvint', downsample=downsample)
        _out = allocate_table(alloc_NP, return_vel=return_vel, return_pid=return_pid, dtype=dtype)
    # Know the final size right away
    _out = _out[:len(data)]

    unpack_rvint(data, header['BoxSize'], float_dtype=dtype, posout=_out['pos'], velout=_out['vel'] if return_vel else None)

    if zspace:
        _out['pos'][:,0] += _out['vel'][:,0] / state['VelZSpace_to_kms']  # km/s * s Mpc/km /Mpc --> dimensionless box units. 

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
        Whether the format on disk is RVdoublePID or just RVPID.
        Default: False
    tag: bool, optional
        Whether the format on disk is RVPID or just RV.
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
            _out = allocate_table(len(data), return_vel=return_vel, return_pid=return_pid, dtype=dtype)
    
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
    
def read_rvzel(fn, return_pos=True, return_vel=True, return_zel=False, return_pid=False, zspace=False, add_grid=False, boxsize=None, dtype=np.float32, out=None, return_header=False, double=False):
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
    with open(fn, 'rb') as fp, ContextTimer('Read', output=False) as read_timer:
        header = skip_header(fp)
        raw = np.fromfile(fp, dtype=rvzel_dt)

    if header:
        header = InputFile(str_source=header)
    else:
        # Look for a file named 'header'
        try:
            header = InputFile(fn=pjoin(dirname(fn), 'header'))
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
        particles = allocate_table(len(raw), return_pos=return_pos, return_vel=return_vel, return_zel=return_zel, return_pid=return_pid, dtype=dtype)
    else:
        particles = out
        
    if header:
        particles.meta.update(header)
        
    unpack_timer = ContextTimer('unpack')
    unpack_timer.Start()
    
    from . import zeldovich
    D = zeldovich.calc_growth_ratio(header, 1.4, 'init')
    D = 1.
    print(f'Scaling displacements by D={D:.3g}')

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
        
        if return_pos:
            particles['pos'][:len(raw)] = raw['r']*D
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
            
    unpack_timer.stop(report=False)
    
    particles.meta['read_time'] = particles.meta.get('read_time',0.) + read_timer.elapsed
    particles.meta['unpack_time'] = particles.meta.get('unpack_time',0.) + unpack_timer.elapsed

    retval = (particles,) if out is None else (len(raw),)
    
    if return_header:
        retval += (header,)
        
    if len(retval) == 1:
        return retval[0]

    return retval
    
    
def read_state(fn, make_global=True, dtype=np.float32, dtype_on_disk=np.float32,
                return_pid='auto', return_aux=False, return_vel='auto', return_pos='auto', return_header=False,
                pid_bitmask=0x7fff7fff7fff, out=None, zspace=False, boxsize=None):
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
    pid_bitmask: int, optional
        Bitmask to apply to the aux field to get the PID.
        Pre-AbacusSummit sims used 0xffffffffff.
        
    Returns
    -------
    data: ndarray of length (npart,)
        The particle data.  Positions are in data['pos']; velocities are in data['vel'] (if `return_vel`),
        and PIDs are in data['pid'] (if `return_pid`).
    or,
    NP: int
        If `out` is given, returns the number of rows read into `out`.
    """
    base = basename(fn)
    if return_pos == 'auto':
        return_pos = 'position' in base
    if return_vel == 'auto':
        return_vel = 'velocity' in base
    if return_pid == 'auto':
        return_pid = 'auxillary' in base 
    if make_global == 'auto':
        make_global = 'position' in base

    if zspace:
        raise NotImplementedError
        
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
        try:
            state = InputFile(pjoin(dirname(fn), 'state'))
        except FileNotFoundError:
            state = None
    
    # Count the particles to read in
    _format = {np.float32:'state', np.float64:'state64'}[dtype_on_disk]
    NP = get_alloc_np(pos_fn, format=_format)
    
    # Use the given buffer
    if out is not None:
        particles = out
    else:  # or allocate one
        particles = allocate_table(NP, return_vel=return_vel, return_pid=return_pid, return_aux=return_aux, dtype=dtype)
    del fn

    read_timer = ContextTimer('read', cumulative=True, output=False)
    unpack_timer = ContextTimer('unpack', cumulative=True, output=False)
    
    if return_pos:
        with read_timer:
            particles['pos'][:NP] = np.fromfile(pos_fn, dtype=(dtype_on_disk, 3))
        if make_global:
            with read_timer:
                cellinfo = np.fromfile(ci_fn, dtype=ci_dtype)
            unpack_timer.Start()
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
            unpack_timer.stop(report=False)
                
    if return_vel:
        with read_timer:
            particles['vel'][:NP] = np.fromfile(vel_fn, dtype=(dtype_on_disk, 3))
    if return_aux:
        with read_timer:
            particles['aux'][:NP] = np.fromfile(aux_fn, dtype=np.uint64)
    if return_pid:
        if return_aux:
            with unpack_timer:
                particles['pid'][:NP] = particles['aux'][:NP] & pid_bitmask
        else:
            with read_timer:  # TODO: technically unpacking as well
                particles['pid'][:NP] = np.fromfile(aux_fn, dtype=np.uint64) & pid_bitmask

    # TODO: particles could be out, which is probably a new table object, so meta won't get propagated...
    particles.meta['read_time'] = read_timer.elapsed + particles.meta.get('read_time',0.)
    particles.meta['unpack_time'] = unpack_timer.elapsed + particles.meta.get('unpack_time',0.)
            
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
    # or allocate one
    else:
        alloc_NP = get_alloc_np(fn, format='pack14')
        _out = allocate_table(alloc_NP, return_vel=return_vel, return_pid=return_pid, dtype=dtype)
    
    # Now merge the pos and vel fields, if present. It's a little faster in Numba.
    dtlist = []
    if 'pid' in _out.colnames:
        dtlist += [('pid',_out['pid'].dtype)]
        # Numba doesn't like None, so use len 0 arrays instead 
        pid = _out['pid']
    else:
        pid = np.empty(0, dtype=np.uint64)
    
    if 'vel' in _out.colnames:
        dtlist += [('posvel',np.float32,6)]
    else:
        dtlist += [('posvel',np.float32,3)]
        
    merged_dt = np.dtype(dtlist, align=True)
    _out = _out.view(dtype=merged_dt)
    posvel = _out['posvel']
    
    raise NotImplementedError('finish pack14lite conversion to tables')
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


def read_asdf(fn, colname=None, out=None, return_pos='auto', return_vel='auto', return_pid='auto', dtype=np.float32,
                load_header=True, verbose=True, blosc_threads=None, zspace=False, **kwargs):
    '''
    ASDF format used with AbacusSummit.  This interface is designed for the scenario where
    the distinction between halo and field is unimportant.  For halo-ortiented access, use
    abacus_halo_catalog.

    TODO: maybe this can be our "backdoor pilot" to converting ReadAbacus to astropy tables
    '''

    base = basename(fn)
    if return_pos == 'auto':
        return_pos = 'rv' in base
    if return_vel == 'auto':
        return_vel = 'rv' in base
    if return_pid == 'auto':
        return_pid = 'pid' in base

    if return_pid:
        raise NotImplementedError('pid reading not finished')

    asdf_data_key = kwargs.pop('asdf_data_key','data')
    asdf_header_key = kwargs.pop('asdf_data_key','header')
    #if len(kwargs) != 0:
    #    raise ValueError(f'Unknown args! {kwargs}')

    if blosc_threads is None:
        blosc_threads = get_nthreads_by_format('asdf')

    import asdf.compression
    try:
        asdf.compression.validate('blsc')
    except:
        # Note: this is a temporary solution until blosc is integrated into ASDF, or until we package a pluggable decompressor
        raise RuntimeError('Error: your ASDF installation does not support Blosc compression.  Please run "pip install git+https://github.com/lgarrison/asdf.git"')
    asdf.compression.set_decompression_options(nthreads=blosc_threads)

    from . import abacus_halo_catalog
    import astropy.table

    with asdf.open(fn, lazy_load=True, copy_arrays=True) as af:
        if colname is None:
            _colnames = ['rvint', 'pack9', 'packedpid']
            for cn in _colnames:
                if cn in af.tree[asdf_data_key]:
                    if colname is not None:
                        raise ValueError(f"More than one key of {_colnames} found in asdf file {fn}. Need to specify colname!")
                    colname = cn
            if colname is None:
                raise ValueError(f"Could not find any of {_colnames} in asdf file {fn}. Need to specify colname!")

        header = af.tree[asdf_header_key]

        with ContextTimer('Read ASDF', output=False) as timer:
            data = af.tree['data'][colname][:]
        if verbose:
            print(f'Read {data.nbytes/1e6:.4g} MB from {basename(fn)} in {timer.elapsed:.4g} sec at {data.nbytes/1e6/timer.elapsed:.4g} MB/s')

        maxN = len(data)
        if out is not None:
            _out = out
            if load_header:
                out.meta.update(header)
        else:
            _out = allocate_table(maxN, return_pos=return_pos, return_vel=return_vel)
            _out.meta.update(header)

        _out.meta['read_time'] = timer.elapsed + _out.meta.get('read_time',0.)

        timer = ContextTimer('Unpack bits')
        timer.Start()
        _posout = _out['pos'] if return_pos else False
        _velout = _out['vel'] if return_vel else False
        if colname == 'rvint':
            npos,nvel = abacus_halo_catalog.unpack_rvint(data, header['BoxSize'], float_dtype=dtype, posout=_posout, velout=_velout,
                                                         zspace=zspace, VelZSpace_to_kms=header['VelZSpace_to_kms'])
            nread = max(npos,nvel)
        elif colname == 'pack9':
            p9_unpacker = packN_readers[9][dtype]
            if _posout is False:
                _posout = None  # ctypes takes None as the null pointer
            if _velout is False:
                _velout = None
            # unpack using the same number of threads as blosc
            nread = p9_unpacker(data, data.nbytes, blosc_threads, zspace, 1.1, _posout, _velout, None)
        elif colname == 'packedpid':
            justpid, lagr_pos, tagged, density = abacus_halo_catalog.unpack_pids(data, header['BoxSize'], header['ppd'])
        timer.stop(report=False)
        if verbose:
            totalbytes = sum(_out[n].nbytes for n in _out.colnames)
            print(f'Unpacked {totalbytes/1e6:.4g} MB from {basename(fn)} in {timer.elapsed:.4g} sec at {totalbytes/1e6/timer.elapsed:.4g} MB/s')

        _out.meta['unpack_time'] = timer.elapsed + _out.meta.get('unpack_time',0.)

        if out is not None:
            return nread
        else:
            return _out[:nread]


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
def allocate_table(N, return_pos=True, return_vel=True, return_pid=False, return_zel=False, return_aux=False, dtype=np.float32):
    """
    Construct an empty Astropy table of length N to hold particle information
    """
    ndt_list = []
    if return_pid:
        ndt_list += [('pid', np.int64)]
    if return_aux:
        ndt_list += [('aux', np.uint64)]
    if return_pos:
        ndt_list += [('pos', (dtype, 3))]
    if return_vel:
        ndt_list += [('vel', (dtype, 3))]
    if return_zel:
        ndt_list += [("zel", (np.uint16, 3))]
    ndt = np.dtype(ndt_list, align=True)

    particles = Table()
    for field in ndt_list:
        particles.add_column(np.empty(N, dtype=field[1]), copy=False, name=field[0])

    return particles

# Set up a few library utils
try:
    from . import abacus
    ralib = ct.cdll.LoadLibrary(pjoin(abacus.abacuspath, "clibs", "libreadabacus.so"))

    # Set up the arguments and return type for the library functions
    # TODO: switch our C library to use CFFI
    for f in (ralib.unpack_pack14, ralib.unpack_pack9, ralib.unpack_pack14f, ralib.unpack_pack9f):
        f.restype = ct.c_uint64
        f.argtypes = (ndarray_arg, ct.c_size_t, ct.c_int, ct.c_int, ct.c_double, ndarray_arg, ndarray_arg, ndarray_arg)
        
    packN_readers = {9: {np.float32:ralib.unpack_pack9f,
                         np.float64: ralib.unpack_pack9 },
                    14: {np.float32:ralib.unpack_pack14f,
                         np.float64: ralib.unpack_pack14}
                    }
except (OSError, ImportError):
    #raise
    packN_readers = None
    pass  # no pack14 library found




def get_alloc_np(fns, format, downsample=None):
    '''An upper limit on the number of particles in a file,
        based on its header information or size on disk,
        accounting for downsampling'''
    format = format.lower()

    if type(fns) is str:
        fns = (fns,)

    # If reading ASDF files, can use the YAML info
    if 'asdf' in format:
        if downsample:
            raise NotImplementedError("ASDF downsampling not yet implemented")
        asdf_data_key = 'data'
        Ntot = 0
        for fn in fns:
            with asdf.open(fn, lazy_load=True) as af:
                N = None
                for col in af.tree[asdf_data_key]:
                    _N = len(af.tree[asdf_data_key][col])
                    if N is None:
                        N = _N
                    else:
                        if _N != N:
                            # Could allow colname specification, but this should suffice for now
                            raise ValueError(f'Columns of different length found in ASDF data tree: {_N} and {N}')
            Ntot += N
        return Ntot

    # For ordinary files, need to guess based on filesize
    if format not in psize_on_disk:
        raise ValueError(f'Particle size on disk not known for {format}')

    fsize = sum(getsize(fn) for fn in fns)
    ds = 1. if downsample is None else downsample

    ds = max(min(1.,ds),0.)  # clip to [0,1]

    N = fsize/psize_on_disk[format]*ds

    # No downsampling uncertainty
    if downsample is None:
        return int(np.ceil(N))

    # Add 6-sigma downsampling uncertainty
    N += 6*np.sqrt(N)

    return int(np.ceil(N))



psize_on_disk = {'pack14': 14, 'pack14_lite':14, 'pack9': 9, 'rvint': 3*4, 'rvdouble': 6*8, 'state64':3*8, 'state':3*4, 'rvzel':32, 'rvtag':32, 'rvdoubletag':7*8}
reader_functions = {'pack14':read_pack14, 'pack9':read_pack9, 'pack14_lite':read_pack14_lite,
                    'rvdouble':read_rvdouble, 'rvzel':read_rvzel, 'rvdoublezel':read_rvdoublezel, 'rvdoubletag':read_rvdoubletag, 'rvtag':read_rvtag, 'rv':read_rvtag,
                    'gadget':read_gadget,
                    'desi_hdf5':read_desi_hdf5,
                    'state':read_state,
                    'rvint':read_rvint,
                    'asdf':read_asdf, 'asdf_b':read_asdf, 'asdf_a':read_asdf, 'asdf_pack9':read_asdf}
default_file_patterns = {'pack14':'*.dat',
                         'pack9':('*L0_pack9.dat', '*field_pack9.dat'),
                         'rvint' : ('*rv_A*', '*rv_B*'),
                         'state':'position_*',
                         'desi_hdf5':'*.hdf5',
                         'gadget':'*.*',
                         'asdf':('*.asdf','field_*/*.asdf','halo_*/*.asdf'),
                         'asdf_a':('*_A_*.asdf','*_A/*.asdf',),
                         'asdf_b':('*_B_*.asdf','*_B/*.asdf',),
                         'asdf_pack9':('field_pack9/*.asdf','L0_pack9/*.asdf'),
                         }
fallback_file_pattern = '*.dat'
                    
default_box_on_disk = {'desi_hdf5':'box',
                'pack14':1.,
                'pack9':1.,
                'rvtag':'box',
                'gadget':'box',
                'state':1.,
                'asdf_a':'box',
                'asdf_b':'box',
                'asdf':'box',
                'asdf_pack9':1.,
}

_nthreads_by_format = dict(asdf=BLOSC_THREADS)

def get_nthreads_by_format(format, default_nthreads=1):
    format = format.lower()
    if 'asdf' in format:  # capture 'asdf_A' and 'asdf_B'
        return _nthreads_by_format['asdf']
    elif format in _nthreads_by_format:
        return _nthreads_by_format[format]
    return default_nthreads


def get_file_patterns(format, return_pos=True, return_vel=True, return_pid=False, **kwargs):
    # TODO
    format = format.lower()

    # For ASDF, we want to select the subdirectories based on the fields requested
    if format in ('asdf','asdf_a','asdf_b'):
        pats = []
        AB = re.match(r'asdf(?:_(?P<AB>\w))?', format).group('AB')
        if AB:
            AB = AB.upper()
        else:
            AB = '*'
        if return_pos or return_vel:
            pats += [f'field_rv_{AB}_*.asdf', f'field_rv_{AB}/field_rv_{AB}_*.asdf',
                     f'halo_rv_{AB}_*.asdf', f'halo_rv_{AB}/halo_rv_{AB}_*.asdf']
        if return_pid:
            # TODO: read_asdf is designed to do this for us. But what if only one of rv/pid exist?
            pats += [f'field_pid_{AB}_*.asdf', f'field_pid_{AB}/field_pid_{AB}_*.asdf',
                     f'halo_pid_{AB}_*.asdf', f'halo_pid_{AB}/halo_pid_{AB}_*.asdf']
        return pats

    try:
        pats = default_file_patterns[format]
    except KeyError:
        pats = fallback_file_pattern

    if type(pats) is str:
        pats = (pats,)
    return pats


def get_files_from_path(path, format=None, pattern=None, key=None, **kwargs):
    '''
    Given a path (which can be a single file, globbing patern, dir,
    or list of these), return the corresponding list of files.  If passing
    directories, must provide a format or pattern.

    key is an optional sorting key.
    '''

    _key = (lambda k: key(basename(k))) if key else None

    if type(path) is str:
        paths = (path,)
    else:
        # If not a string, must already be a list
        try:
            len(path)
        except:
            raise ValueError(f'path must be a string or iterable, not {type(path)}!')
        paths = path

    if type(pattern) is str:
        patterns = (pattern,)

    files = []
    for path in paths:
        # Determine the file names to read
        if type(path) is not str:
            raise ValueError(f'All path values must be str. Found {type(path)}')

        if isdir(path):
            if pattern is None:
                if is_ic_path(path):
                    patterns = ('ic_*',)
                else:
                    if format is None:
                        raise ValueError("Must give file format if not giving pattern!")
                    patterns = get_file_patterns(format, **kwargs)

            for p in patterns:
                files += sorted(glob(pjoin(path, p)), key=_key)
        elif isfile(path):
            files += [path]
        else:
            # Not a dir and not a file; interpret as a glob pattern
            files += glob(path)

    files = sorted(files, key=_key, reverse=True)

    if not files:
        raise ValueError(f'No files found matching "{path}" for format "{format}"!')

    return files


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
