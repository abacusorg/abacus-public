#!/usr/bin/env python

# Warning!  h5py 2.6 has a few bugs that may affect the following code.

import h5py
from pkg_resources import parse_version
assert parse_version(h5py.__version__) > parse_version('2.6')

import argparse
from glob import glob
import numpy as np
import pdb
import os
from os.path import join as pjoin
from Abacus.InputFile import InputFile

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

def convert_cats(halo_dirs):
    if type(halo_dirs) == str:
        halo_dirs = [halo_dirs]
        
    for halo_dir in halo_dirs:
        abacus_header = InputFile(pjoin(halo_dir, 'header'))
        for binfn in glob(pjoin(halo_dir, 'halos_*.bin')):
            pfn = binfn.replace('/halos_', '/particles_')
            
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
                assert np.isfinite(bin_halos[field]).all()

            # write the halos as hdf5
            hfn = binfn[:-3] + 'h5'
            with h5py.File(hfn, 'w') as hfp:
                dset = hfp.create_dataset('halos', data=bin_halos)
                dset.attrs.update(vars(abacus_header))
                
            # Now read the particles
            with open(pfn, 'rb') as pfp:
                particles = np.fromfile(pfp, dtype=particle_dt)
                
            # particle sanity checks
            assert len(particles) == bin_header['num_particles']
            assert len(particles) == bin_halos['subsamp_len'].sum(), (len(particles), bin_halos['subsamp_len'].sum())
            for field in particles.dtype.fields:
                assert np.isfinite(particles[field]).all()
            
            # write the particles as hdf5
            hfn = pfn[:-3] + 'h5'
            with h5py.File(hfn, 'w') as hfp:
                pdset = hfp.create_dataset('particles', data=particles)
                pdset.attrs.update(vars(abacus_header))
            
            # Remove the halo and particle binary files
            os.remove(binfn)
            os.remove(pfn)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Convert Rockstar binary halo catalogs and particle subsamples to HDF5')
    parser.add_argument('halo_dirs', help='Directory containing *.bin halo catalogs and *.particles particle subsamples', nargs='+', metavar='HALO_DIR')
    
    args = parser.parse_args()

    convert_cats(args.halo_dirs)
    