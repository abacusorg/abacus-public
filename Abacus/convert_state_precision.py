#!/usr/bin/env python3
'''
Convert state positions, velocities, and cellinfos from one
float precision to another; e.g. float to double.  This is
mostly useful in the context of checking the effect of floating
point precision on Abacus forces.

This routine does nothing to the multipoles & Taylors; they will
have to be re-generated.

The input state will be renamed with its original dtype, and the
new state will have the original name.  E.g.
    `read` --> `read_float32`
    `read` will be float64
'''

import argparse
import os.path
from os.path import join as pjoin

import numpy as np
import shutil

from Abacus import Tools


class Converter:

    def __init__(self, dir, from_dtype, to_dtype):
        dir = os.path.normpath(dir)

        self.dir = dir
        self.from_dtype = from_dtype
        self.to_dtype = to_dtype

        self.from_ci_dtype = np.dtype([('startindex',np.uint64),
                             ('count',np.int32), ('active',np.int32),
                             ('mean_square_velocity', from_dtype), ('max_component_velocity', from_dtype), ('max_component_acceleration', from_dtype)],
                             align=True)
        self.to_ci_dtype = np.dtype([('startindex',np.uint64),
                             ('count',np.int32), ('active',np.int32),
                             ('mean_square_velocity', to_dtype), ('max_component_velocity', to_dtype), ('max_component_acceleration', to_dtype)],
                             align=True)

        dtype_names = {np.float32:'float32', np.float64:'float64'}
        self.outdir = dir + '_' + dtype_names[to_dtype]
        self.dirbackup = dir + '_' + dtype_names[from_dtype]

        os.makedirs(self.outdir, exist_ok=True)

    def convert(self, fn):
        orig = np.fromfile(pjoin(self.dir, fn), dtype=self.from_dtype)
        new = orig.astype(self.to_dtype)

        assert np.isfinite(new).all()

        new.tofile(pjoin(self.outdir, fn))

    def copy(self, fn):
        shutil.copy(pjoin(self.dir, fn), pjoin(self.outdir, fn))

    def convert_cellinfo(self, fn):
        orig = np.fromfile(pjoin(self.dir, fn), dtype=self.from_ci_dtype)
        new = orig.astype(self.to_ci_dtype)

        for field in self.to_ci_dtype.fields:
            assert np.isfinite(new[field]).all()

        new.tofile(pjoin(self.outdir, fn))

    def finish(self):
        for fn in ['state', 'slabsize', 'redlack', 'globaldipole'] :
            try:
                shutil.copy(pjoin(self.dir, fn), pjoin(self.outdir, fn))
            except:
                raise
                pass  # not required to exist

        # Move read --> read_float32
        # Move read_float64 --> read
        os.rename(self.dir, self.dirbackup)
        os.rename(self.outdir, self.dir)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=Tools.ArgParseFormatter)
    parser.add_argument('state_dir', help='The directory containing the state')
    parser.add_argument('--from-prec', help='Input state precision', type=str, choices=['float','double'], default='float')
    parser.add_argument('--to-prec', help='Output state precision', type=str, choices=['float','double'], default='double')
    parser.add_argument('--nslab', help='Number of slabs to convert.  <=0 means all slabs.', type=int, default=-1)

    args = parser.parse_args()
    args = vars(args)

    dtype_dict = {'float':np.float32, 'double':np.float64}
    args['from_prec'] = dtype_dict[args['from_prec']]
    args['to_prec'] = dtype_dict[args['to_prec']]

    c = Converter(args['state_dir'], args['from_prec'], args['to_prec'])

    if args['nslab'] <= 0:
        nslab = len(glob(pjoin(args['state_dir'], 'position_*')))
    else:
        nslab = args['nslab']

    for i in range(nslab):
        c.convert('position_{:04d}'.format(i))
        c.convert('velocity_{:04d}'.format(i))
        c.convert_cellinfo('cellinfo_{:04d}'.format(i))
        c.copy('auxillary_{:04d}'.format(i))

    c.finish()