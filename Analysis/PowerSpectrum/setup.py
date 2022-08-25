import os
from setuptools import setup
from pathlib import Path

from pybind11.setup_helpers import Pybind11Extension, build_ext

ABACUS = Path(os.getenv('ABACUS'))

ext_modules = [
    Pybind11Extension(
        "PowerSpectrum.pslib",
        ['PowerSpectrum/bind.cpp'],
        include_dirs=[],
        libraries = ['gsl', 'gslcblas', 'tbb'],
    ),
]

setup(ext_modules=ext_modules, cmdclass={"build_ext": build_ext})
