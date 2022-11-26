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

#include <sys/time.h>
#include <errno.h>
#include <string.h>
#include <fcntl.h>

#include <ctype.h>
#include <sys/select.h>
#include <stdlib.h>

#include <stdint.h>
#include <unistd.h>

#define uint32 uint32_t
#define int32 int32_t
#define uint64 uint64_t
#define int64 int64_t

#define __STDC_FORMAT_MACROS
#include <inttypes.h>

typedef std::complex<double> Complex;

// So that we can easily adjust to double precision for kinematics
#ifdef DOUBLEPRECISION
#define FLOAT double
#define FLOAT3 double3
#else 
#define FLOAT float 
#define FLOAT3 float3
#endif

// dtype for Multipoles/Taylors and Derivatives on disk
// might eventually be used for computations as well
#ifdef DOUBLEPRECISION
typedef std::complex<double> MTCOMPLEX;
typedef double DFLOAT;
#else
typedef std::complex<float> MTCOMPLEX;
typedef float DFLOAT;
#endif

#define LIKELY(x) __builtin_expect((x),1)
#define UNLIKELY(x) __builtin_expect((x),0)

using namespace std;

