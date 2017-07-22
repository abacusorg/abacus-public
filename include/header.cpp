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

#include <sys/time.h>
#include <errno.h>
#include <string.h>
#include <fcntl.h>

#include <ctype.h>
#include <sys/select.h>
#include <stdlib.h>

#include <stdint.h>
#include <unistd.h>

#define uint64 uint64_t

#define __STDC_FORMAT_MACROS
#include <inttypes.h>

typedef std::complex<double> Complex;

// dtype for Multipoles/Taylors and Derivatives on disk
// might eventually be used for computations as well
#ifdef DOUBLEPRECISION
typedef std::complex<double> MTCOMPLEX;
typedef double DFLOAT;
#else
typedef std::complex<float> MTCOMPLEX;
typedef float DFLOAT;
#endif

using namespace std;

