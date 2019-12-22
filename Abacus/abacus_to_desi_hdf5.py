#!/usr/bin/env python3
'''
This script converts Abacus particle data to the HDF5 format in the 
schema set by the DESI cosmo-sim working group:

https://desi.lbl.gov/trac/wiki/CosmoSimsWG/IC

'''


import argparse
import os
from os.path import dirname, basename, abspath, join as pjoin
import re
import gc
import multiprocessing

from Abacus import Tools
from Abacus import ReadAbacus
from Abacus.InputFile import InputFile
import Abacus.Cosmology

import numpy as np
import h5py

class Converter:
    def __init__(self, **kwargs):
        self.kwargs = kwargs

    def __call__(self, fns):
        kwargs = self.kwargs
        format = kwargs['format']
        ds = kwargs['ds']
        dtype = kwargs.pop('dtype')

        read_args = dict(format=format, dtype=dtype, return_vel=True, return_pid=True, return_header=True, return_fn=True, verbose=True)
        if ds and ds > 1:
            read_args['return_zel'] = True
        if format == 'rvzel':
            read_args['add_grid'] = True

        for (p,header),fn in ReadAbacus.AsyncReader(fns, **read_args):
            convert(p, header, fn, **kwargs)
        
        print('All done.')


def convert(p, header, fn, format, ds=None, out_parent=None):
    ppd = np.int64(np.round(header.NP**(1./3)))
    assert ppd**3 == header.NP

    if ds and ds > 1:
        assert int(ds) == ds
        ppd_ds = ppd//ds
        assert ppd/ds == ppd_ds
        p = p.reshape(-1,ppd,ppd)

        # Check that the first plane all agrees on its index
        izel0 = int(p['zel'][0,0,0,0])
        assert (p['zel'][0,...,0] == izel0).all()
        
        s = -izel0 % ds

        # Downsample by the given factor. We have to check the index of the first plane to decide whether to include it
        p = p[s::ds, ::ds, ::ds]
        p = p.reshape(-1)

        header['NP'] = np.int64(header['NP']/ds**3)
        header['ParticleMassHMsun'] *= ds**3
        header['ParticleMassMsun'] *= ds**3
    else:
        ppd_ds = ppd

    # better way to know the box on disk?
    if format == 'rvzel':
        p['vel'] *= header.VelZSpace_to_kms/header.BoxSize
    else:
        p['pos'] *= header.BoxSize
        p['vel'] *= header.VelZSpace_to_kms

    # Wrap the particles from [0,L), without changing the origin
    pos_wrapped = np.ascontiguousarray(p['pos'])
    Tools.wrap_zero_origin(pos_wrapped, header.BoxSize)

    # Copy over header values from Abacus to the DESI cosmo sim names
    h5header = {}
    for k in ['BoxSize', 'InitialRedshift', ('NP', 'NP.Matter'), 'Omega_M', 'Omega_DE', 'H0',
                ('ParticleMassHMsun', 'ParticleMass.Matter'), ('CPD', 'NumFilesPerSnapshot'),
                'f_growth', 'Redshift']:
        
        if type(k) is tuple:
            k,h5k = k
        else:
            h5k = k

        h5header[h5k] = header[k]

    # 'Growth' in Abacus is normalized to D~a at early times, not D~1 at z=0
    cosmo = Abacus.Cosmology.from_params(header, z=header['Redshift'])
    h5header['GrowthRatio'] = cosmo.current.growth/cosmo.today.growth
    #h5header['GrowthRatio'] = header.ZD_Pk_sigma/header['sigma_8']

    # `HubbleNow` in Abacus is in units of H0
    h5header['HubbleNow'] = header['HubbleNow']*header['H0']
    # `VelZSpace_to_kms` in Abacus is the inverse of RSDFactor
    h5header['RSDFactor'] = (header['VelZSpace_to_kms']/header['BoxSize'])**-1

    if out_parent is None:
        out_parent = dirname(dirname(fn))
    h5dir = pjoin(out_parent, basename(dirname(fn)) + '_desi_hdf5_' + str(ppd_ds))
    os.makedirs(h5dir, exist_ok=True)

    name_stem = basename(fn)
    # TODO: more general extension detection
    if name_stem.endswith('.dat'):
        name_stem = name_stem[:-4]

    h5fn = pjoin(h5dir, name_stem + '.hdf5')

    with h5py.File(h5fn, 'w') as fp:
        fp.create_dataset("/Matter/Position", data=pos_wrapped)
        fp.create_dataset("/Matter/Velocity", data=p['vel'])
        fp.create_dataset("/Matter/ParticleID", data=p['pid'])

        # should this be a group or empty dataset?
        dset = fp.create_group("/AbacusHeader")
        for k in vars(header):
            dset.attrs[k] = header[k]

        dset = fp.create_group("/Header")
        for k in h5header:
            dset.attrs[k] = h5header[k]

    del p, pos_wrapped
    gc.collect()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=Tools.ArgParseFormatter)
    parser.add_argument('abacus-file', help='One or more Abacus files', nargs='+')
    parser.add_argument('--format', help='The Abacus file format', choices=ReadAbacus.reader_functions.keys(), default='rvzel')
    parser.add_argument('--ds', help='Downsample-per-dimension factor', type=int)
    parser.add_argument('--out-parent', help='Directory in which to create a directory for the outputs')
    parser.add_argument('--dtype', help='Precision in which to read and write the outputs', choices=['float32', 'float64'], default='float32')
    parser.add_argument('--nproc', help='Number of processes to use to process files.', default=multiprocessing.cpu_count(), type=int)
    
    args = parser.parse_args()
    args = vars(args)

    abacus_files = args.pop('abacus-file')
    nproc = args.pop('nproc')
    print(f'Processing {len(abacus_files)} files with {nproc} processes', flush=True)

    converter = Converter(**args)

    chunks = [abacus_files[s::nproc] for s in range(nproc)]  # sad that chunksize won't do this!
    with multiprocessing.get_context("fork").Pool(nproc) as pool:
        pool.map(converter, chunks)
