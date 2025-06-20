#!/usr/bin/env python3

# Warning!  h5py 2.6 has a few bugs that may affect the following code.

import argparse
from glob import glob
from pkg_resources import parse_version
import os
from os.path import join as pjoin, dirname, basename
import psutil
import multiprocessing
import sys
import gc

import numpy as np
import numpy.lib.recfunctions as rfn

from Abacus.InputFile import InputFile
from Abacus.Tools import ArgParseFormatter

# Binary data types

# from io_internal.h
'''
struct binary_output_header {
    uint64_t magic;
    int64_t snap, chunk;
    float scale, Om, Ol, h0;
    float bounds[6];
    int64_t num_halos, num_particles;
    float box_size, particle_mass;
    int64_t particle_type;
    int32_t format_revision;
    char rockstar_version[VERSION_MAX_SIZE];
    char unused[BINARY_HEADER_SIZE - (sizeof(char)*VERSION_MAX_SIZE) - (sizeof(float)*12) - sizeof(int32_t) - (sizeof(int64_t)*6)];
};
'''
header_dt = np.dtype([('magic', np.uint64),
                        ('snap', np.int64), ('chunk', np.int64),
                        ('scale', np.float32), ('Om', np.float32), ('Ol', np.float32), ('h0', np.float32),
                        ('bounds', np.float32, 6),
                        ('num_halos', np.int64), ('num_particles', np.int64),
                        ('box_size', np.float32), ('particle_mass', np.float32),
                        ('particle_type', np.int64),
                        ('format_revision', np.int32),
                        ('rockstar_version', 'S12'),
                        ('unused', np.uint8, 256-112)], align=True)
ROCKSTAR_MAGIC = 0xfadedacec0c0d0d0
PARTICLE_TYPE_FULL = 2
assert header_dt.itemsize == 256


# From halo.h
'''
struct halo {
  int64_t id, parent_id;
  float pos[6], corevel[3], bulkvel[3];
  float m, m_SO, r, child_r, vmax_r, mgrav, vmax, rvmax, rs, klypin_rs, vrms,
    J[3], energy, spin, alt_m[4], alt_m_SO[4], Xoff, Voff, b_to_a, c_to_a, A[3],
    b_to_a2, c_to_a2, A2[3],
    bullock_spin, kin_to_pot, m_pe_b, m_pe_d, halfmass_radius;
  int64_t num_p, num_child_particles, p_start, desc, flags, n_core, subsamp_start, subsamp_len;
  float min_pos_err, min_vel_err, min_bulkvel_err;
};
'''

halo_dt = np.dtype([('id',np.int64), ('parent_id',np.int64),
                   ('pos',np.float32,3), ('vel',np.float32,3), ('corevel',np.float32,3), ('bulkvel',np.float32,3),
                   ('m',np.float32), ('m_SO',np.float32), ('r',np.float32), ('child_r',np.float32), ('vmax_r',np.float32), ('mgrav',np.float32), ('vmax',np.float32), ('rvmax',np.float32), ('rs',np.float32), ('klypin_rs',np.float32), ('vrms',np.float32),
                   ('J',np.float32,3), ('energy',np.float32), ('spin',np.float32), ('alt_m',np.float32,4), ('alt_m_SO',np.float32,4), ('Xoff',np.float32), ('Voff',np.float32), ('b_to_a',np.float32), ('c_to_a',np.float32), ('A',np.float32,3), ('b_to_a2',np.float32), ('c_to_a2',np.float32), ('A2',np.float32,3),
                   ('bullock_spin',np.float32), ('kin_to_pot',np.float32), ('m_pe_b',np.float32), ('m_pe_d',np.float32), ('halfmass_radius',np.float32),
                   ('num_p',np.int64), ('num_child_particles',np.int64), ('p_start',np.int64), ('desc',np.int64), ('flags',np.int64), ('n_core',np.int64), ('subsamp_start',np.int64), ('subsamp_len',np.int64),
                   ('min_pos_err',np.float32), ('min_vel_err',np.float32), ('min_bulkvel_err',np.float32)],
                     align=True)

'''
struct particle {
  int64_t id;
  float pos[6];
};
'''
particle_dt = np.dtype([('id', np.int64), ('pos', np.float32, 3), ('vel', np.float32, 3)], align=True)

# Convert one halo catalog .bin file and its associated particle file
def convert_binfile(binfn):
    gc.collect()
    import h5py
    assert parse_version(h5py.__version__) > parse_version('2.6')

    assert(type(binfn) == str)

    halo_dir = dirname(binfn)
    abacus_header = InputFile(pjoin(halo_dir, 'header'))
    rockstar_cfg = InputFile(pjoin(halo_dir, 'rockstar.cfg'))

    # Don't bother with the particle files if they're empty!
    no_particles = rockstar_cfg['SUBSAMPLE_FRAC'] == 0.

    with open(binfn, 'rb') as bfp:
        # Read the binary header
        bin_header = np.fromfile(bfp, dtype=header_dt, count=1)
        assert bin_header['magic'] == ROCKSTAR_MAGIC
        assert bin_header['particle_type'] == PARTICLE_TYPE_FULL
        
        # read the binary halos
        bin_halos = np.fromfile(bfp, dtype=halo_dt)
        
    # halo sanity checks
    assert len(bin_halos) == bin_header['num_halos']
    for field in bin_halos.dtype.fields:
        if field not in ['m_pe_b', 'rs']:  # these are sometimes inf; that's probably something we can ignore and just pass into the catalogs
            assert np.isfinite(bin_halos[field]).all(), (field, binfn)

    # Count the particles
    particle_mass = abacus_header['ParticleMassHMsun']

    h5halos_N = np.empty(bin_halos.shape, dtype=np.dtype([('N', np.int32), ('alt_N', np.int32,4), ('N_SO', np.int32), ('alt_N_SO', np.int32,4)]))
    h5halos_N['N'] = np.round(bin_halos['m'] / particle_mass)
    h5halos_N['alt_N'] = np.round(bin_halos['alt_m'] / particle_mass)
    h5halos_N['N_SO'] = np.round(bin_halos['m_SO'] / particle_mass)
    h5halos_N['alt_N_SO'] = np.round(bin_halos['alt_m_SO'] / particle_mass)
    bin_halos = rfn.merge_arrays([bin_halos, h5halos_N], flatten=True, usemask=False)

    # Check that the mass is a multiple of the count
    assert np.allclose(bin_halos['N'] * particle_mass, bin_halos['m'], rtol=1e-4)
    assert np.allclose(bin_halos['N_SO'] * particle_mass, bin_halos['m_SO'], rtol=1e-4)

    # write the halos as hdf5
    hfn = binfn[:-3] + 'h5'
    with h5py.File(hfn, 'w') as hfp:
        dset = hfp.create_dataset('halos', data=bin_halos)
        dset.attrs.update(vars(abacus_header))
        
    # Now read the particles
    if not no_particles:
        pfn = binfn.replace('/halos_', '/particles_')
        with open(pfn, 'rb') as pfp:
            particles = np.fromfile(pfp, dtype=particle_dt)
            
        # particle sanity checks
        assert len(particles) == bin_header['num_particles']
        assert len(particles) == bin_halos['subsamp_len'].sum(), (len(particles), bin_halos['subsamp_len'].sum())
        for field in particles.dtype.fields:
            good = np.isfinite(particles[field])
            try:
                assert good.all()
            except:
                print(field, pfn)
                bad = ~(good.all(axis=-1))
                print(particles[bad])
                raise
        
        # write the particles as hdf5
        hfn = pfn[:-3] + 'h5'
        with h5py.File(hfn, 'w') as hfp:
            pdset = hfp.create_dataset('particles', data=particles)
            pdset.attrs.update(vars(abacus_header))

        # Remove the particle binary file
        os.remove(pfn)
    
    # Remove the halo binary file
    os.remove(binfn)
    gc.collect()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Convert Rockstar binary halo catalogs and particle subsamples to HDF5',
      formatter_class=ArgParseFormatter)
    parser.add_argument('halo_dirs', help='Directory containing *.bin halo catalogs and (optionally) *.particles particle subsamples',
      nargs='+', metavar='HALO_DIR')
    parser.add_argument('--nproc', help='Number of processes to use to process files.',
      default=len(psutil.Process().cpu_affinity()), type=int)

    args = parser.parse_args()

    binfiles = []
    for hd in args.halo_dirs:
        binfiles += glob(pjoin(hd, 'halos_*.bin'))

    print(f'Processing {len(binfiles)} files with {args.nproc} processes')
    # The Pool context manager calls terminate() instead of close().  Sad!
    pool = multiprocessing.get_context("fork").Pool(args.nproc)
    pool.map(convert_binfile, binfiles)  #, chunksize=len(binfiles)//args.nproc)
    pool.close()
    pool.join()
