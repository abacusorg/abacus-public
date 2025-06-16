#!/usr/bin/env python
# Copyright 2012-2025 The Abacus Developers
# SPDX-License-Identifier: GPL-3.0-or-later

'''
Helper routines for the merger tree modules.

TODO: this file was copied over from abacus_mergertree and has not yet been modernized

Original filename: merger_tree_library.py
'''

import sys

import numba
import numpy as np
from abacusnbody.data.compaso_halo_catalog import CompaSOHaloCatalog
from numba import jit


def read_halo_catalogue(step_now, return_header = False):
    print("Reading halo catalogue from step %s"%(step_now))
    sys.stdout.flush()

    cat      = CompaSOHaloCatalog(
                step_now,
                load_subsamples="AB_halo_pid",
                convert_units=True,
                fields = "merger",
                unpack_bits = "merger",
                )
    pids     = cat.subsamples[:]["pid"].data
    rho      = cat.subsamples[:]["density"].data
    nstartA  = cat.halos[:]["npstartA"].data
    ntagA    = cat.halos[:]["npoutA"].data
    nstartB  = cat.halos[:]["npstartB"].data
    ntagB    = cat.halos[:]["npoutB"].data
    pos      = cat.halos[:]["x_L2com"].data
    vmax     = cat.halos[:]["vcirc_max_L2com"].data
    nphalo   = cat.halos[:]["N"].data
    numhalos = cat.numhalos
    ntag     = ntagA + ntagB

    nslice = cat.header["FullStepNumber"]
    mpart  = cat.header["ParticleMassHMsun"]
    z      = cat.header["Redshift"]
    box    = cat.header["BoxSizeHMpc"]
    mhalo  = nphalo * mpart

    if return_header:
        return cat.header, box, nslice, z, numhalos, nphalo, mhalo, pos, vmax, nstartA, ntagA, nstartB, ntagB, ntag, pids, rho
    else:
        return box, nslice, z, numhalos, nphalo, mhalo, pos, vmax, nstartA, ntagA, nstartB, ntagB, ntag, pids, rho

# It's probably worth giving these haloes a unique identifier

def indxxHalo(Nslice, Nhaloes):
    indxx  = Nslice * 1e12 + np.arange(Nhaloes)
    indxx  = indxx.astype(int)
    return indxx

@jit(nopython=True)
def indxxHaloSlabwise(Nslice, numhaloArray, filextArray):
    indxxArray = np.zeros(np.sum(numhaloArray), dtype = numba.int64)
    for nn in range(len(numhaloArray)):
        if nn == 0:
            startindex = 0
            endindex = startindex + numhaloArray[nn]
        else:
            startindex = np.sum(numhaloArray[:nn])
            endindex = startindex + numhaloArray[nn]
        indxxArray[startindex : endindex] = Nslice*1e12 + filextArray[nn]*1e9 + np.arange(numhaloArray[nn])
    return indxxArray

# Create a utility for resizing/appending to an existing HDF5 array

def writeDset(dataset, appendArray):
    oldLen = dataset.len()
    newLen = oldLen + len(appendArray)
    try:
        dataset.resize((newLen,appendArray.shape[1]))
    except IndexError:
        dataset.resize((newLen,))
    dataset[oldLen:newLen] = appendArray

# Function to associate particle IDs with a unique halo ID
def filled_array(i, start, end, length):
    out = np.zeros((length), dtype=int)
    np.add.at(out,start,i)
    np.add.at(out,end,-i)
    return out.cumsum()
