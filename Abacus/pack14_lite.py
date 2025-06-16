# Copyright 2012-2025 The Abacus Developers
# SPDX-License-Identifier: GPL-3.0-or-later

'''
This file contains a number of functions for reading our
pack14 file format in Python.  Our standard way to do this
is with a C library, but it's often easier to distribute
Python code than C.  The following functions are implemented
in Numba in an attempt to get them to be not horribly
slower than C.  Currently, they're "only" 2-3x slower
in my tests.

This file is not intended to be used on its own; the
main entry point is the `ReadAbacus` library (specifically
the "pack14_lite" file format will call these functions).

Note: the "lite" in "pack14_lite" refers to the way
we read it (i.e. with Python/Numba), not file format
itself, which is ordinary pack14.
'''

import numpy as np
import numba as nb

cell_header_spec = [ ('vscale', nb.int16),
                     ('cpd', nb.int16),
                    ('ijk', nb.int16[:]) ]
@nb.jitclass(cell_header_spec)
class cell_header:
    def __init__(self, vscale, cpd, ijk):
        self.vscale = vscale
        self.cpd = cpd
        self.ijk = ijk
        
    def islegal(self):
        return self.vscale > 0

@nb.njit
def iscell(c):
    return c[0] == 0xff

@nb.njit
def unpack_cell(c, ch, s):
    #expand_to_short(c, s)  # numba is bad at inlining this
    s1 = np.int16(c[1])
    s4 = np.int16(c[4])
    s7 = np.int16(c[7])
    
    low = np.int16(0xf)
    high = np.int16(0xf0)
    four = np.int16(4)
    
    s[1] = (((s1&high)<<four)  | np.int16(c[2])) - np.int16(2048)
    s[2] = ((s4&low) | (np.int16(c[3])<<four)) - np.int16(2048)
    s[3] = (((s4&high)<<four)  | np.int16(c[5])) - np.int16(2048)
    s[4] = ((s7&low) | (np.int16(c[6])<<four)) - np.int16(2048)
    s[5] = (((s7&high)<<four)  | np.int16(c[8])) - np.int16(2048)

    ch.cpd = s[1] + np.int16(2000)
    ch.vscale = s[2] + np.int16(2000)
    ch.ijk = s[3:6] + np.int16(2000)
    
    assert ch.islegal()

@nb.njit
def expand_to_short(c, s):    
    s1 = np.int16(c[1])
    s4 = np.int16(c[4])
    s7 = np.int16(c[7])
    
    low = np.int16(0xf)
    high = np.int16(0xf0)
    four = np.int16(4)

    s[0] = ((s1&low) | (np.int16(c[0])<<four)) - np.int16(2048)
    s[1] = (((s1&high)<<four)  | np.int16(c[2])) - np.int16(2048)
    s[2] = ((s4&low) | (np.int16(c[3])<<four)) - np.int16(2048)
    s[3] = (((s4&high)<<four)  | np.int16(c[5])) - np.int16(2048)
    s[4] = ((s7&low) | (np.int16(c[6])<<four)) - np.int16(2048)
    s[5] = (((s7&high)<<four)  | np.int16(c[8])) - np.int16(2048)

@nb.njit
def unpack_id(c):
    id = (np.uint64(c[9])<<np.uint64(32)) | (c[10:14].view(np.uint32)[0])
    return id

@nb.njit
def unpack(c, header, s, posvel, pid, n):
    # maybe we would prefer to only unpack the relevant bytes
    #expand_to_short(c, s)
    s1 = np.int16(c[1])
    s4 = np.int16(c[4])
    s7 = np.int16(c[7])
    
    low = np.int16(0xf)
    high = np.int16(0xf0)
    four = np.int16(4)

    invcpd = 1./header.cpd

    # maybe we would prefer separate functions instead of `if` statements?
    # it barely makes a difference in the timings
    if len(pid) > 0:
        pid[n] = unpack_id(c)
    s[0] = ((s1&low) | (np.int16(c[0])<<four)) - np.int16(2048)
    s[1] = (((s1&high)<<four)  | np.int16(c[2])) - np.int16(2048)
    s[2] = ((s4&low) | (np.int16(c[3])<<four)) - np.int16(2048)
    for i in range(3):
        posvel[n,i] = (s[i]/2000. + header.ijk[i] + 0.5)*invcpd - 0.5
    if posvel.shape[1] > 3:
        s[3] = (((s4&high)<<four)  | np.int16(c[5])) - np.int16(2048)
        s[4] = ((s7&low) | (np.int16(c[6])<<four)) - np.int16(2048)
        s[5] = (((s7&high)<<four)  | np.int16(c[8])) - np.int16(2048)
        for i in range(3,6):
            posvel[n,i] = s[i]/2000*header.vscale*invcpd

@nb.njit
def _read_pack14(raw, posvel, pid):
    n = 0
    maxnp = len(raw)
    
    header = cell_header(0, 0, np.zeros(3, dtype=np.int16))
    _short = np.empty(6, dtype=np.int16)
    
    for i in range(maxnp):
        if iscell(raw[i]):
            unpack_cell(raw[i], header, _short)
        else:
            unpack(raw[i], header, _short, posvel, pid, n)
            n += 1
    return n