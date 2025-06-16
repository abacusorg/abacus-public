# Copyright 2012-2025 The Abacus Developers
# SPDX-License-Identifier: GPL-3.0-or-later

# Build the AbacusCosmo extension using pybind11

import os
from setuptools import setup
from pathlib import Path

from pybind11.setup_helpers import Pybind11Extension, build_ext

ABACUS = Path(os.getenv('ABACUS'))

ext_modules = [
    Pybind11Extension(
        "AbacusCosmo",
        ['bind.cpp'],
        include_dirs=[ABACUS/'include'],
    ),
]

setup(ext_modules=ext_modules, cmdclass={"build_ext": build_ext})
