#!/usr/bin/env python3

'''
Compile the fast_cksum library.
'''

import os
from os.path import join as pjoin

from cffi import FFI
ffibuilder = FFI()

ffibuilder.cdef("""
	uint32_t crc32_fast(const void* data, size_t length);
	uint32_t crc32_fast_partial(const void* data, size_t length);
	uint32_t crc32_fast_finalize(size_t total_length, uint32_t previousCrc32);
	""")
ffibuilder.cdef("""
	uint32_t crc32_fast_partial(const void* data, size_t length, uint32_t previousCrc32);
	""", override=True)

ffibuilder.set_source("_crc32_fast",  # name of the output C extension
"""
    #include "crc32_fast.h"
""",
	include_dirs=[pjoin(os.pardir, os.pardir, 'include')],
    sources=[pjoin(os.pardir, os.pardir, 'include', 'crc32_fast.cpp')],
    extra_compile_args=[])

if __name__ == "__main__":
    ffibuilder.compile(verbose=True)
