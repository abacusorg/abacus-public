#!/usr/bin/env python
'''
This script converts Abacus particle data to the HDF5 format in the 
schema set by the DESI cosmo-sim working group:

https://desi.lbl.gov/trac/wiki/CosmoSimsWG/IC

'''


import argparse
import os
from os.path import dirname, basename, abspath, join as pjoin
import re

from Abacus import Tools
from Abacus import ReadAbacus
from Abacus.InputFile import InputFile
import Abacus.Cosmology

import pynbody
import numpy as np
import h5py

def convert(fn, format, ds=None):
    header = InputFile(pjoin(dirname(fn), 'header'))
    ppd = np.int64(np.round(header.NP**(1./3)))
    assert ppd**3 == header.NP

    p = ReadAbacus.read(fn, format=format, dtype=np.float64,
                        return_vel=True, return_pid=True, return_zel=True, add_grid=True, boxsize=header.BoxSize)

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

    p['vel'] *= header.VelZSpace_to_kms/header.BoxSize

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
    h5header['GrowthRatio'] = header.ZD_Pk_sigma/header['sigma_8']
    # `HubbleNow` in Abacus is in units of H0
    h5header['HubbleNow'] = header['HubbleNow']*header['H0']
    # `VelZSpace_to_kms` in Abacus is the inverse of RSDFactor
    h5header['RSDFactor'] = (header['VelZSpace_to_kms']/header['BoxSize'])**-1
    
    h5dir = dirname(fn) + '_desi_hdf5_' + str(ppd_ds)
    os.makedirs(h5dir, exist_ok=True)

    h5fn = pjoin(h5dir, basename(fn) + '.hdf5')

    with h5py.File(h5fn, 'w') as fp:
        fp.create_dataset("/Matter/Position", data=p['pos'])
        fp.create_dataset("/Matter/Velocity", data=p['vel'])
        fp.create_dataset("/Matter/ParticleID", data=p['pid'])

        # should this be a group or empty dataset?
        dset = fp.create_group("/AbacusHeader")
        for k in vars(header):
            dset.attrs[k] = header[k]

        dset = fp.create_group("/Header")
        for k in h5header:
            dset.attrs[k] = h5header[k]


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=Tools.ArgParseFormatter)
    parser.add_argument('abacus_file', help='One or more Abacus files', nargs='+')
    parser.add_argument('--format', help='The Abacus file format', choices=ReadAbacus.reader_functions.keys(), default='rvzel')
    parser.add_argument('--ds', help='Downsample-per-dimension factor', type=int)
    
    args = parser.parse_args()
    args = vars(args)

    abacus_files = args.pop('abacus_file')

    for i,af in enumerate(abacus_files):
        print('Converting {}/{} (\"{}\")... '.format(i+1,len(abacus_files), basename(af)), end='', flush=True)
        convert(af, **args)
        print('Done.')
