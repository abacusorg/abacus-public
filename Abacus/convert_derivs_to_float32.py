#!/usr/bin/env python3
# Copyright 2012-2025 The Abacus Developers
# SPDX-License-Identifier: GPL-3.0-or-later

"""
usage: convert_derivs_to_float32.py [-h]
                                    fourierspace_file [fourierspace_file ...]

Takes 64-bit derivatives files and makes a 32-bit copy of them. Does not
modify the original derivatives. The new files will be named
"float32/fourierspace_float32_*".

This script is primarily invoked automatically through abacus.py.

positional arguments:
  fourierspace_file  The files to convert. Usually something like
                     "fourierspace_1485_8_2_8_*".

optional arguments:
  -h, --help         show this help message and exit
"""

import argparse
import sys
from pathlib import Path

import numpy as np


def convert(dpath, from_dtype=np.float64, to_dtype=np.float32, tag='float32'):
  dpath = Path(dpath)
  ddir = dpath.parent
  dfn = dpath.name

  if not dfn.startswith('fourierspace'):
    print(
      "Error: file '{}' does not begin with fourierspace!  Wrong file?".format(
        dpath
      )
    )
    sys.exit(1)

  derivs = np.fromfile(dpath, dtype=from_dtype)

  # Create the output directory if it doesn't exist
  output_dir = ddir / tag
  output_dir.mkdir(exist_ok=True)

  new_fn = output_dir / dfn.replace('fourierspace', 'fourierspace_' + tag)

  new_derivs = derivs.astype(to_dtype)
  assert np.isfinite(new_derivs).all(), (
    'Derivatives not finite after conversion.  Did an exponent overfloat float32?'
  )

  new_derivs.tofile(new_fn)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Takes 64-bit derivatives files and makes a 32-bit copy of them.  Does not modify the original derivatives.  The new files will be named "fourierspace_float32_*".'
    )
    parser.add_argument(
        'fourierspace_file',
        help='The files to convert. Usually something like "fourierspace_1485_8_2_8_*".',
        nargs='+',
    )

    args = parser.parse_args()

    for dfile in args.fourierspace_file:
        convert(dfile)
