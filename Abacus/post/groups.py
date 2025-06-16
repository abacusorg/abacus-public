#!/usr/bin/env python
# Copyright 2012-2025 The Abacus Developers
# SPDX-License-Identifier: GPL-3.0-or-later


"""
Convert raw binary outputs from Abacus's on-the-fly halo/group finder
into ASDF file format.  This is done once after completion of a
simulation; thereafter, users will use abacusutils.data.compaso_halo_catalog
to read the ASDF catalogs.  Thus, end users should not have to use this file.

For developers examining the raw binary catalogs, the Numpy dtype
"halo_dt" defined in this file may be useful as well.

This file used to be called `convert_raw_groups_to_asdf.py`.

Usage
=====
$ python -m Abacus.post.groups --help

Details
=======
Abacus writes halo information to the 'group' directory.
Typically, it writes 9 files per slab in the GroupDirectory:
    - `halo_info_SSSS`
        The L1 groups (aka halos).  Contains stats like
        CoM positions and velocities and moments of the particles.
        Also indicates the number of subsampled particles in the
        `halo_pids_A/B` and `halo_rv_A/B` files.
        
    - `halo_rv_A_SSSS` and `halo_rv_B_SSSS` (`halo_rv_A/B_SSSS` for short)
        A subsample of particles in L1 groups, in `RVint` format.
        The subsampling is done as a function of the PID;
        "subsample-able" particles are the same as taggable particles.
        The particles are written in the same order as `halos`.
        
    - `halo_pids_A_SSSS` and `halo_pids_B_SSSS` (`halo_pids_A/B_SSSS` for short)
        The corresponding 64-bit particle IDs, which also contain
        information about the Lagrangian position of the particles,
        whether they are tagged and their local density.
    
    - `field_rv_A_SSSS` and `field_rv_B_SSSS` (`field_rv_A/B_SSSS` for short)
        Same as `halo_rv_A/B_SSSS`, but only for the field particles.

    - `field_pids_A_SSSS` and`field_pids_B_SSSS` (`field_pids_A/B_SSSS` for short)
        Same as `halo_pids_A/B_SSSS`, but only for the field particles.
        
"""

import gc
import multiprocessing
import os
import shutil
import sys
import tempfile
import warnings
from pathlib import Path

import asdf
import cffi
import click
import numpy as np

from Abacus.fast_cksum.cksum_io import CksumReader, CksumWriter
from Abacus.InputFile import InputFile
from Abacus.Tools import ContextTimer

# The path to the halostat C struct definition
halostat_cstruct_h = Path(os.getenv('ABACUS')) / 'singlestep' / 'FOF' / 'halostat_cstruct.h'

# Eventually these will be external dependencies
fast_cksum_store = Path(os.getenv('ABACUS')) / 'external' / 'fast-cksum' / 'bin' / 'fast_cksum_store'
fast_cksum_cat = Path(os.getenv('ABACUS')) / 'external' / 'fast-cksum' / 'bin' / 'fast_cksum_cat'

@click.command
@click.argument('group_dir_or_file')
@click.option('-c', '--chunk', default=50,
              help='Target number of slabs per superslab',
              )
@click.option('--delete', is_flag=True,
              help='Delete the input files on success',
              )
@click.option('-v', '--verbose', is_flag=True,
              help='Print status messages',
              )
@click.option('-w', '--nworkers', default=1,
              help='Number of workers/subprocesses to use. Parallelism is over superslabs.',
              )
@click.option('--empty-ok', is_flag=True, default=True,
              help='Whether a catalog without any halo_info files will raise an error',
              )
def convert(group_dir_or_file, chunk=50, asdf_data_key='data', asdf_header_key='header', asdf_compression='blsc', delete=False,
            verbose=True, nworkers=1, empty_ok=True):
    """
    Convert the raw binary halo_info files in the given group directory into ASDF format.
    The ASDF files will consist of "super slabs" of approximately `chunk` slabs each.

    Particle PIDs and RVs will be concatenated by the `chunk` factor as well, but they
    stay in raw binary, not ASDF.  The concatenation of these files requires reindexing
    the particle indices in the halo_info files, which this function also does.

    This function may also be invoked on a single halo_info file instead of a group directory;
    this would mainly be used for debugging.

    Parameters
    ----------
    group_dir_or_file: str
        The path to a group directory or a single halo_info file

    chunk: int, optional
        The number of slabs to concatenate into each chunk.
        Default: 50

    asdf_data_key: str, optional
        The name of the key to use for storing the halo catalog in the output ASDF file.
        The default choice of 'data' allows Astropy to read in the catalog without
        additional user input.  This choice needs to be synchronized with abacusutils.

    asdf_header_key: str, optional
        The name of the key to use for storing the Abacus header in the output ASDF file.
        Astropy Tables use the name 'meta' instead of 'header', but this technically isn't
        an Table ASDF file, so they won't automatically read this field regardless.  So the
        choice just needs to be synchronized with abacusutils.
        Default: 'header'

    asdf_compression: str, optional
        The compression algorithm to use for the ASDF halo table columns.
        Note that one cannot use memory-mapped IO on the resulting catalogs.
        More details here: https://asdf.readthedocs.io/en/latest/asdf/arrays.html#compression
        Default is 'blsc', which is Abacus's tuned Blosc compression

    verbose: bool, optional
        Print extra status information.

    nworkers: int, optional
        Number of simultaneous workers/subprocesses for chunk processing. Default: 1

    empty_ok: bool, optional
        Whether a totally empty catalog (no halo_info files) will raise an error or not.
        Default: False
            
    """

    group_dir_or_file = Path(group_dir_or_file).resolve(strict=True)

    if group_dir_or_file.is_dir():
        groupdir = group_dir_or_file
        oneslab = False
    else:
        if not group_dir_or_file.name.startswith('halo_info_'):
            raise ValueError('If passing file, must be halo_info file')
        groupdir = group_dir_or_file.parent
        oneslab = int(group_dir_or_file.name.split('.')[0].split('_')[-1])
    del group_dir_or_file  # use groupdir

    # Number of slabs
    header = dict(InputFile(groupdir / 'header'))
    cpd = int(header['CPD'])

    # twoD?
    NodeSizeZ = header.get('NodeSizeZ', 1)
    twoD = NodeSizeZ > 1

    # Find the closest redshift in the header and use that in the output dir name
    # Want to be consistent across sims, and cat z's are best-effort!
    groupdir_asdf = (groupdir.parents[1] /
        groupdir.parent.name.replace('group','halos') /
        'z{:.3f}'.format(find_closest_output_redshift(header, header['Redshift']))
    )
    groupdir_asdf.mkdir(parents=True, exist_ok=False)
    shutil.copy2(groupdir / 'header', groupdir_asdf)

    # When computing chunk boundaries, don't skip empty slabs, just to maintain synchronicity among sims
    # Note we're spreading out any unevenness in the chunk divisions here, rather than sticking it all at the end
    if oneslab is False:
        nchunks = max(cpd//chunk,1)
        chunk_boundaries = np.linspace(0,cpd, num=nchunks+1, endpoint=True, dtype=int)
    else:
        nchunks = 1
        chunk_boundaries = np.array([oneslab,oneslab+1])
        cpd = 1
    del chunk  # now use chunk_boundaries

    if verbose:
        print(f'\n==== Processing {cpd} slabs from {groupdir} in {nchunks} chunks')

    file_types = ['halo_info', 'halo_pid_A', 'halo_pid_B', 'field_pid_A', 'field_pid_B',
                    'halo_rv_A', 'halo_rv_B', 'field_rv_A', 'field_rv_B']
    source_type_overrides = {'halo_pid_A':'halo_pids_A',
                             'halo_pid_B':'halo_pids_B',
                             'field_pid_A':'field_pids_A',
                             'field_pid_B':'field_pids_B'}

    file_patterns = {}
    for name in file_types:
        fnstem = source_type_overrides.get(name, name)  # de-pluralize pids from here on out
        if twoD:
            file_patterns[name] = groupdir / (fnstem + '_{:04d}.{:03d}')
        else:
            file_patterns[name] = groupdir / (fnstem + '_{:04d}')

    index_files_timer = ContextTimer('Indexing files')

    index_files_timer.Start()
    all_fns = {}  # all_fns[file_type][chunk_index] is a list of files for that type and chunk
    # why not just use flat arrays? because different files are allowed to be missing for different types
    for t in file_types:
        all_fns[t] = []
        for c in range(nchunks):
            lo,hi = chunk_boundaries[[c,c+1]]
            if twoD:
                _fns = [Path(str(file_patterns[t]).format(i,k)) for i in range(lo,hi) for k in range(NodeSizeZ)]
            else:
                _fns = [Path(str(file_patterns[t]).format(i)) for i in range(lo,hi)]

            if not empty_ok or t != 'halo_info':
                _fns = [f for f in _fns if f.is_file()]
            all_fns[t] += [_fns]

        if any(all_fns[t]):
            (groupdir_asdf / t).mkdir(exist_ok=False)

        # Empty particle files are not output (?)
        #assert len(all_fns[t]) <= len(all_fns['halo_info'])
        # Other consistency checks?
    index_files_timer.stop(report=True)
        
    if not empty_ok and not sum(len(chunk_fns) for chunk_fns in all_fns['halo_info']):
        raise FileNotFoundError(f'Could not find any halo_info files matching {file_patterns["halo_info"]}')

    # We need to know if we have field subsamples because the halo files will cotain both L0 and L1 if so (TODO: is this right?)
    # TODO: more robust way to indicate this?  Field in header?
    have_field_subsamples = any(any(chunk_fns for chunk_fns in all_fns[name]) for name in file_types if name.startswith('field'))

    if verbose:
        if have_field_subsamples:
            print('This catalog has field subsamples')
        else:
            print('This catalog does not have field subsamples')

    # A list of dict, one per chunk
    chunk_args = []
    for c in range(nchunks):
        chunk_fns = {name:all_fns[name][c] for name in file_types}
        if len(chunk_fns['halo_info']) == 0 and not have_field_subsamples:
            print('Found completely empty chunk; skipping')
            continue
        chunk_args += [dict(c=c, chunk_fns=chunk_fns)]

    # function object for passing state
    # TODO: getting pretty verbose, maybe the convert() function should just be the function object
    chunk_converter = ChunkConverter(verbose=verbose, have_field_subsamples=have_field_subsamples,
                                    groupdir_asdf=groupdir_asdf,
                                    asdf_header_key=asdf_header_key, asdf_data_key=asdf_data_key,
                                    header=header, asdf_compression=asdf_compression)

    if verbose:
        print(f'About to launch {nworkers} subprocesses to process {len(chunk_args)} chunks')

    with multiprocessing.Pool(processes=nworkers) as pool:
        pool.map(chunk_converter, chunk_args)

    # All chunks done; do any epilogue work

    # Merge all checksum files
    with ContextTimer("Merge checksum files"):
        for name in file_types:
            d = groupdir_asdf / name
            if not d.is_dir():
                continue
            lines = []
            fns = list(d.glob('*.crc32'))
            assert 'checksums.crc32' not in [fn.name for fn in fns]
            for fn in fns:
                with open(fn,'r') as fp:
                    lines += fp.readlines()
            lines.sort(key=lambda L:L.split()[-1])
            with (d / 'checksums.crc32').open('w') as fp:
                fp.writelines(lines)
            for fn in fns:
                fn.unlink()

    # Delete original files
    # We wait until all chunks are done to facilitate easy recovery
    # It does mean higher peak disk usage, though
    if delete:
        # just touch a file so we notice if we stopped partway through deletion for any reason
        # although the presence of the merged checksum files is a good indicator too
        with (groupdir / 'READY_FOR_DELETION').open('w'):
            pass
        for chunk_arg in chunk_args:  # for chunk
            for file_type in chunk_arg['chunk_fns']:  # for file type
                for fn in chunk_arg['chunk_fns'][file_type]:  # for file
                    fn.unlink(missing_ok=True)

        (groupdir / 'checksums.crc32').unlink()
        (groupdir / 'header').unlink()
        (groupdir / 'READY_FOR_DELETION').unlink()
        groupdir.rmdir()  # should be empty!  If not, we missed a file

        try:
            # The last job removes the group dir
            groupdir.parent.rmdir()
        except OSError:
            pass  # not empty


class ChunkConverter:
    def __init__(self, verbose, have_field_subsamples, groupdir_asdf, asdf_header_key, asdf_data_key, header, asdf_compression):
        self.verbose = verbose
        self.have_field_subsamples = have_field_subsamples
        self.groupdir_asdf = groupdir_asdf
        self.asdf_header_key = asdf_header_key
        self.asdf_data_key = asdf_data_key
        self.header = header
        self.asdf_compression = asdf_compression


    def __call__(self, kwargs):
        return self.convert_chunk(**kwargs)


    def convert_chunk(self, c, chunk_fns):  # c is chunk num
        print()  # newline

        # Amount of memory stored for each halo
        per_halo = halo_dt.itemsize

        file_types = list(chunk_fns.keys())

        n_halo_info_files = len(chunk_fns['halo_info'])
        if chunk_fns['halo_info']:
            groupdir = chunk_fns['halo_info'][0].parent
        else:
            groupdir = chunk_fns['field_rv_B'][0].parent
        cksum_fn = groupdir / 'checksums.crc32'

        # make a column container for the whole chunk
        nhalo_eachslab = np.array([get_nelem_from_filesize(fn, per_halo) for fn in chunk_fns['halo_info']], dtype=int)
        nhalo_thischunk = nhalo_eachslab.sum()
        halo_start = np.zeros(n_halo_info_files+1,dtype=int)
        halo_start[1:] = np.cumsum(nhalo_eachslab)
        chunkhalos = {col:np.empty(nhalo_thischunk, dtype=halo_dt[col]) for col in halo_dt.fields}

        if self.verbose:
            print(f'Found {nhalo_thischunk} halos ({nhalo_thischunk:.3g}) in chunk {c}')

        transpose_timer = ContextTimer('Transpose halos', cumulative=True, output=False)
        reindex_timer = ContextTimer('Reindex halos')

        if nhalo_thischunk:
            # now read each slab in the chunk, doing the transpose as we go
            cksum_reader = CksumReader(cksum_fn, verbose=True)  # could do outside chunk loop
            for i,hfn in enumerate(chunk_fns['halo_info']):
                if not hfn.is_file():
                    continue  # might have to skip empties if empty_ok
                hstart,hend = halo_start[[i,i+1]]
                slabhalos = np.frombuffer(cksum_reader(hfn), dtype=halo_dt)
                with transpose_timer:
                    for col in chunkhalos:
                        chunkhalos[col][hstart:hend] = slabhalos[col]
                del slabhalos
                gc.collect()
            print('Read halo info time: {:.1f} seconds, {:.1f} MB/s IO, {:.2f} GB/s CRC32'.format(
                cksum_reader.io_timer.elapsed, cksum_reader.bytes_read/1e6/cksum_reader.io_timer.elapsed,
                cksum_reader.bytes_read/1e9/cksum_reader.checksum_timer.elapsed))
            del cksum_reader
            transpose_timer.report()

            # reindex the chunk
            reindex_timer.Start()
            if not self.have_field_subsamples:
                njump1 = reindex_contiguous_subsamples(chunkhalos['npstartA'], chunkhalos['npoutA'])  # noqa: F841
                njump2 = reindex_contiguous_subsamples(chunkhalos['npstartB'], chunkhalos['npoutB'])  # noqa: F841
                #assert njump1 == n_halo_info_files-1, (njump1, n_halo_info_files)  # might fail for very small slices
                #assert njump2 == n_halo_info_files-1, (njump2, n_halo_info_files)
            else:
                #TODO: if we have L0 but no halos, this will break
                reindex_subsamples_from_filesize(chunkhalos, chunk_fns['halo_pid_A'], 'A', nhalo_eachslab, np.int64().itemsize)
                reindex_subsamples_from_filesize(chunkhalos, chunk_fns['halo_pid_B'], 'B', nhalo_eachslab, np.int64().itemsize)
            reindex_timer.stop(report=True)

            # All done with halo records; write them to ASDF
            tree = {self.asdf_header_key:self.header,
                    self.asdf_data_key:chunkhalos}
            asdf_fn = self.groupdir_asdf / 'halo_info' / f'halo_info_{c:03d}.asdf'
            # TODO: per-column tunings?
            compression_kwargs = dict(typesize='auto', shuffle='shuffle', compression_block_size=12*1024**2, blosc_block_size=3*1024**2, nthreads=4)
            with asdf.AsdfFile(tree) as af, CksumWriter(asdf_fn) as fp:
                af.write_to(fp, all_array_compression=self.asdf_compression, compression_kwargs=compression_kwargs)
            #asdf_size = os.stat(asdf_fn).st_size  # do we want to report compressed or uncompressed IO?

            tot_time = fp.io_timer.elapsed + fp.checksum_timer.elapsed
            print('Write halo info into ASDF time: {:.1f} seconds, {:.1f} MB/s IO, {:.2f} GB/s CRC32 (compression: {})'.format(tot_time,
                    fp.bytes_written/1e6/fp.io_timer.elapsed, fp.bytes_written/1e9/fp.checksum_timer.elapsed, self.asdf_compression))
            del chunkhalos, tree, af
            gc.collect()

        # Concatenate the corresponding particle slabs
        cat_timer = ContextTimer('Concatenate particle files', cumulative=True, output=False)
        cat_timer.nbytes = 0
        particle_write_time = 0
        particle_cksum_time = 0
        particle_bytes_written = 0
        cksum_reader = CksumReader(cksum_fn, verbose=True)
        for name in file_types:
            # halo_info is the asdf, already did that
            if name == 'halo_info':
                continue

            if not chunk_fns[name]:
                # Probably just don't have subsamples
                continue
                # However, could be a sparse chunk, in which case we might want to write a blank file(?)
                # raise NotImplementedError

            outfn = self.groupdir_asdf / name / (name + f'_{c:03d}.asdf')
            #assert not isfile(outfn), f'Concatenated file {outfn} already exists?'

            if '_rv_' in name:
                pdtype = np.int32  # keep signedness, makes later conversion mildly more convenient
                pshape = (-1,3)
                Mway = 12
                colname = 'rvint'
            else:
                pdtype = np.uint64
                pshape = (-1,)
                Mway = 8
                colname = 'packedpid'

            compression_kwargs = dict(typesize=Mway, shuffle='bitshuffle', compression_block_size=12*1024*1024, blosc_block_size=3*1024*1024, nthreads=4)
            npart_catted = sum(_cfn.stat().st_size for _cfn in chunk_fns[name])/pdtype().itemsize
            assert int(npart_catted) == npart_catted
            npart_catted = int(npart_catted)
            catted = np.empty(npart_catted, dtype=pdtype)
            i = 0
            for cfn in chunk_fns[name]:
                new = np.frombuffer(cksum_reader(cfn), dtype=pdtype)
                with cat_timer:
                    catted[i:i+len(new)] = new
                i += len(new)
                del new
                gc.collect()
            assert i == len(catted)
            cat_timer.nbytes += catted.nbytes
            catted = catted.reshape(*pshape)

            tree = {'data': {colname:catted},
                    'header': self.header}
            with asdf.AsdfFile(tree) as af, CksumWriter(outfn) as fp:
                af.write_to(fp, all_array_compression=self.asdf_compression, compression_kwargs=compression_kwargs)
            particle_write_time += fp.io_timer.elapsed
            particle_cksum_time += fp.checksum_timer.elapsed
            particle_bytes_written += fp.bytes_written
            del catted, tree, af
            gc.collect()

        tot_time = particle_write_time + particle_cksum_time
        if tot_time:
            print('Read particles time: {:.1f} seconds, {:.1f} MB/s IO, {:.2f} GB/s CRC32'.format(
                cksum_reader.io_timer.elapsed, cksum_reader.bytes_read/1e6/cksum_reader.io_timer.elapsed,
                cksum_reader.bytes_read/1e9/cksum_reader.checksum_timer.elapsed))
            print('Concat particles: {:.1f} seconds, {:.1f} MB/s'.format(cat_timer.elapsed, cat_timer.nbytes/1e6/cat_timer.elapsed))
            print('Write + compress particles into ASDF: {:.1f} seconds, {:.1f} MB/s compress+IO, {:.2f} GB/s CRC32 (compression: {})'.format(tot_time,
                    particle_bytes_written/1e6/particle_write_time, particle_bytes_written/1e9/particle_cksum_time, self.asdf_compression))

        # Encourage the print messages to flush to disk
        sys.stdout.flush()
        sys.stderr.flush()


def find_closest_output_redshift(header, targetz):
    '''
    Catalog redshifts are "best effort", meaning we will not shorten the timestep
    to hit them exactly.  But we want to give the directories consistent redshift
    names.  This function tells us what redshift to use for that.  Downstream users
    are responsible for using the redshift in the header, not the directory name.
    '''
    keys = ('TimeSliceRedshifts', 'L1OutputRedshifts', 'TimeSliceRedshifts_Subsample')
    allz = []
    for k in keys:
        try:
            allz += list(header.get(k,[]))
        except TypeError:
           allz += [header.get(k)]
    allz = np.array(allz)
    if allz.size <= 0:
        warnings.warn(f'Failed to get output redshifts from header for consistent naming purposes! Will fall back to actual redshift {targetz}.')
        return targetz
    return allz[np.abs(allz-targetz).argmin()]


def reindex_contiguous_subsamples(subsamp_start, subsamp_len):
    """
    If we concatenate halos and particles into big files/arrays, the "subsample start"
    indices in the halos table no longer correspond to the concatenated particle array.
    But we can easily reconstruct the correct indices as the cumulative sum of the
    subsample sizes.

    This only works if there are no un-indexed particles in the subsample files, such
    as L0 particles that are not indexed as part of any L1 halo.  The indexed particles
    must be contiguous, hence the name of this function.
    
    Parameters
    ----------
    subsamp_start: ndarray
        The concatenated subsample start array. Must be in its original order (i.e. the order
        on disk, because the subsamples on disk are in the same order).  Will be updated
        in-place.

    subsamp_len: ndarray
        The concatenated subsample lengths corresponding to `subsamp_start`.  Will not
        be modified.
        
    Returns
    -------
    njump: int
        The number of discontinuities that were fixed in the subsample indexing.
        Should be equal to the number of file splits.
    """
    
    # The number of "discontinuities" in particle indexing should equal the number of files
    # This helps us make sure the halos were not reordered
    njump = (subsamp_start[:-1] + subsamp_len[:-1] != subsamp_start[1:]).sum()
    
    # Now reindex the halo records
    subsamp_start[0] = 0
    subsamp_start[1:] = subsamp_len.cumsum()[:-1]
    
    return njump



def reindex_subsamples_from_filesize(halos, halo_particle_fns, AB, N_halo_per_file, field_size):
    '''
    For subsample redshifts where we have L1s followed by L0s in the halo_pids files,
    we need to reindex using the total number of PIDs in the file, not the npout fields,
    which only have the L1s.
    '''

    nh = 0
    for k,fn in enumerate(halo_particle_fns):
        nh += N_halo_per_file[k]
        np_thisfile = get_nelem_from_filesize(fn, field_size)
        halos['npstart'+AB][nh:] += np_thisfile


# Get the number of binary records a file contains based on the file size and record size
# Raises an exception if the file size is not divisible by the record size
def get_nelem_from_filesize(filename, elem_size):
    try:
        fsize = Path(filename).stat().st_size
    except FileNotFoundError:  # allow missing halo_info files if empty_ok
        return 0
    nelem = fsize//elem_size
    if nelem*elem_size != fsize:
        raise RuntimeError(f'File size {fsize} of {filename} not divisible by record size {elem_size}')
    return nelem


##################################################################################
# The rest of the file concerns automatically generating the halo_dt Numpy dtype
# from the C struct defintion
##################################################################################

def generate_halostat_dtype(halostat_cstruct_fn, add_int16_suffix=True):
    '''
    Get the Numpy dtype corresponding to the `HaloStat` struct.
    Uses CFFI to generate the dtype automatically by parsing the
    `halostat_cstruct_fn` file.
    
    Note that this is the "raw" dtype and that one probably
    needs to unpack the compressed integer fields into floats.
    To facilitate this, the `add_int16_suffix` will append
    '_i16' or '_u16' to the names of the int16_t or uint16_t fields.
    
    Parameters
    ----------
    halostat_cstruct_fn: str
        The path to the "halostat_cstruct.h" file that contains
        the "struct HaloStat" definition in plain C (no C++).
        Usually in "abacus/singlestep/FOF/".
        
    add_int16_suffix: bool, optional
        Append "_i16" or "_u16" to the int16_t and uint16_t
        field names.  Default: True.
    
    '''
    with open(halostat_cstruct_fn, 'r') as fp:
        s = fp.read()
    
    ffibuilder = cffi.FFI()
    
    # Provide the struct in a CFFI-parseable format
    ffibuilder.cdef(s)
    
    # Trigger the compiler to make sure CFFI got it right
    with tempfile.TemporaryDirectory() as tmpdir:
        ffibuilder.verify(s, tmpdir=tmpdir)
    
    ct = ffibuilder.typeof('struct HaloStat')
    dtype = ctype_to_dtype(ffibuilder, ct)
    
    if add_int16_suffix:
        newnames = []
        for name in dtype.names:
            if dtype[name].base == np.int16:
                newnames += [name + '_i16']
            elif dtype[name].base == np.uint16:
                newnames += [name + '_u16']
            else:
                newnames += [name]
        dtype.names = newnames
    
    # Do some sanity checks
    assert ffibuilder.sizeof('struct HaloStat') == dtype.itemsize
    assert ffibuilder.alignof('struct HaloStat') == dtype.alignment
    
    return dtype

_FLOAT_TYPES = set(['float', 'double', 'long double'])

def _sub_overrides(overrides, prefix):
    out = {}
    for key, value in overrides.items():
        if key.startswith(prefix):
            out[key[len(prefix):]] = value
            
    return out

def ctype_to_dtype(ffi, ctype, overrides=None):
    """
    Convert a CFFI ctype representing a struct to a numpy dtype.
    This does not necessarily handle all possible cases, but should correctly
    account for things like padding.
    
    Adapted from: https://gist.github.com/bmerry/71d67ae0dae8612585f7
    
    Parameters
    ----------
    ffi
        cffi object imported from the library module
    ctype : `CType`
        Type object created from CFFI
    overrides : dict, optional
        Map elements of the type to specified numpy types. The keys are
        strings. To specify an element of a structure, use ``.name``. To
        specify the elements of an array type, use ``[]``. These strings can
        be concatenated (with no whitespace) to select sub-elements.
    """
    if overrides is None:
        overrides = {}
    try:
        return overrides['']
    except KeyError:
        pass

    if ctype.kind == 'primitive':
        if ctype.cname in _FLOAT_TYPES:
            return np.dtype('f' + str(ffi.sizeof(ctype)))
        elif ctype.cname == 'char':
            return np.dtype('c')
        elif ctype.cname == '_Bool':
            return np.dtype(np.bool_)
        else:
            test = int(ffi.cast(ctype, -1))
            if test == -1:
                return np.dtype('i' + str(ffi.sizeof(ctype)))
            else:
                return np.dtype('u' + str(ffi.sizeof(ctype)))
            
    elif ctype.kind == 'struct':
        names = []
        formats = []
        offsets = []
        for field in ctype.fields:
            if field[1].bitsize != -1:
                raise ValueError('bitfields are not supported')
            names.append(field[0])
            sub_overrides = _sub_overrides(overrides, '.' + field[0])
            formats.append(ctype_to_dtype(ffi, field[1].type, sub_overrides))
            offsets.append(field[1].offset)
        
        # LHG: We use align=True here so that dtype.alignment == 8,
        # but it should have no effect on the itemsize or offsets
        return np.dtype(dict(names=names, formats=formats, offsets=offsets, itemsize=ffi.sizeof(ctype)), align=True)
    
    elif ctype.kind == 'array':
        shape = []
        prefix = ''
        while ctype.kind == 'array' and prefix not in overrides:
            shape.append(ctype.length)
            ctype = ctype.item
            prefix += '[]'
        sub_overrides = _sub_overrides(overrides, prefix)
        return np.dtype((ctype_to_dtype(ffi, ctype, sub_overrides), tuple(shape)))
    
    else:
        raise ValueError(f'Unhandled kind {ctype.kind}')

    
halo_dt = generate_halostat_dtype(halostat_cstruct_h)

if __name__ == '__main__':
    convert()
