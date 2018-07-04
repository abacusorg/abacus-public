# Script to build the python module "AbacusCosmo"
# which is an interface to Abacus' Cosmology module.

# Build with the command:
# python setup.py build_ext --inplace

# See README for usage.

from distutils.core import setup, Extension
import os

abacus = os.path.expandvars('$ABACUS')

setup(ext_modules=[Extension('_AbacusCosmo',
                    sources=['{abacus}/include/Cosmology.cpp'.format(abacus=abacus), 'AbacusCosmo.i'],
                    include_dirs=['{abacus}/include'.format(abacus=abacus)],
                    swig_opts=['-I{abacus}/include'.format(abacus=abacus), '-c++'],
                    )
                  ]
      )
      
