#!/usr/bin/env python
# Copyright 2012-2025 The Abacus Developers
# SPDX-License-Identifier: GPL-3.0-or-later


'''
This utility converts particle files, usually time slice outputs, into
IC files that Abacus can read.

While we could teach Abacus to read its own packN files, the L0 output
scheme means the slab domain could be pretty wide. This way, we can
at least provide a proper slab segmentation.

Note that if the particle files use the Abacus "spread out" PID
convention, you will need to use the `--repack_pid` flag to repack
these as contiguous for Abacus IC ingestion.
'''

import argparse
from pathlib import Path
import gc
import re

import numpy as np

from Abacus.Tools import ArgParseFormatter, ContextTimer
from Abacus import ReadAbacus, WriteAbacus

del __builtins__.format, __builtins__.input  # safety...

AUXXPID = 0x7fff
AUXYPID = 0x7fff0000
AUXZPID = 0x7fff00000000
NUMPIDBITS = 45

def convert_to_ic(input, input_format, output=None, output_format='RVPID', cpd=None, box=None, repack_pid=False):
    '''
    Parameters
    ----------
    input: path-like
        Dir containing the particles to convert
    input_format: str
        File format to the input. Accepts any of `ReadAbacus.formats`.
    output: path-like, optional
        Output dir. Will write to a subdir: '{output}/{input}_{output_format}'.
        Default is same as dir as input.
    output_format: str, optional
        The output format. Accepts any of `WriteAbacus.formats`.
    cpd: int, optional
        The CPD of the output. Default is to leave unchanged from input.
    box: float, optional
        The boxsize of the input.  Default is to guess.
    '''
    input = Path(input)
    if output is None:
        output = input.parent
    output = Path(output) / f'{input}_{output_format}'
    output.mkdir(parents=True)
    
    header = ReadAbacus.get_any_header(input, format=input_format)
    
    if cpd is None:
        cpd = header['CPD']
    if box is None:
        box = ReadAbacus.get_box_on_disk('', input_format)
    print(f'Assuming particles on disk have box {box}')
    
    writer = WriteAbacus.SlabWriter(header['NP'], cpd, box, output, output_format, verbose=True)
    
    for data in ReadAbacus.AsyncReader(input, format=input_format,
                                verbose=True, return_pos=True, return_vel=True, return_pid=True):
        if repack_pid:
            data['pid'] = make_pids_contiguous(data['pid'])
        writer(**data)
    del writer
    
def make_pids_contiguous(pid):
    '''PIDs in our output formats use the distributed/non-contiguous convention.
    But the IC reader code assumes contiguous (normal) bit packing, because the ICs
    may be shared between codes, for one thing. So we should just repack the
    non-contiguous PIDs as contiguous in this code.
    '''

    # nothing to do for x
    pid = ((pid & AUXYPID) >> 1) | (pid & ~(AUXYPID >> 1))
    pid = ((pid & AUXZPID) >> 2) | (pid & ~(AUXZPID >> 2))
    pid &= (1<<NUMPIDBITS) - 1  # zero all but the PID bits
    
    return pid

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=ArgParseFormatter)
    parser.add_argument('input', help='The input dir, like slice1.000', nargs='+')
    parser.add_argument('-i', help='Format of the particles being input to this script', choices=ReadAbacus.formats, required=True)
    parser.add_argument('-o', help='The format of the particles being output from this script', choices=WriteAbacus.formats)
    parser.add_argument('--repack-pid', help='Repack the our "spread out" PID bit convention into bit-contiguous PIDs', action='store_true')
    
    args = vars(parser.parse_args())
    
    for indir in args.pop('input'):
        convert_to_ic(indir, **args)
