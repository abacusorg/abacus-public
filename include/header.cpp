// Copyright 2012-2025 The Abacus Developers
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef __HEADER__
#define __HEADER__

// First include global definitions from ./configure
#include "config.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <cstring>
#include <cmath>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <vector>
#include <algorithm>
#include <float.h>
#include <complex>
#include <omp.h>
#include <sys/mman.h>
#include <sched.h>
#include <pthread.h>
#include <filesystem>

#include <sys/time.h>
#include <errno.h>
#include <string.h>
#include <fcntl.h>

#include <ctype.h>
#include <sys/select.h>
#include <stdlib.h>

#include <stdint.h>
#include <unistd.h>

#include "threevector.hh"

// TODO: use compiled lib rather than header-only
#include <fmt/std.h>
#include <fmt/ostream.h>

namespace fs = std::filesystem;

using uint32 = uint32_t;
using int32 = int32_t;
using uint64 = uint64_t;
using int64 = int64_t;

#define __STDC_FORMAT_MACROS
#include <inttypes.h>

#include "AbacusComplex.hh"
using Complex = AbacusComplex<double>;
using FFTComplex = AbacusComplex<double>;

#include "AbacusVector.hh"

// So that we can easily adjust to double precision for kinematics
#ifdef DOUBLEPRECISION
using FLOAT = double;
using FLOAT3 = double3;
#else 
using FLOAT = float;
using FLOAT3 = float3;
#endif

// dtype for Multipoles/Taylors and Derivatives on disk
// might eventually be used for computations as well
#ifdef DOUBLEPRECISION
using MTCOMPLEX = AbacusComplex<double>;
using DFLOAT = double ;
#else
using MTCOMPLEX = AbacusComplex<float>;
using DFLOAT = float;
#endif

#define LIKELY(x) __builtin_expect((x),1)
#define UNLIKELY(x) __builtin_expect((x),0)

template<class C, typename T>
bool contains(const C &a, const C &b, const T e) { return std::find(a, b, e) != b; };

#endif
