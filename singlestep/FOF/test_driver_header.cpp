// This is a code to provide a test for the group finding.

// There are a couple of other unit tests:
// g++ -DTEST -fopenmp -lgomp -O2 fof_sublist.cpp 
// g++ -DTEST -fopenmp -lgomp -O2 slab_accum.cpp

#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <algorithm>
#include <vector>
#include <cstdint>
#include <assert.h>

#include <sys/time.h>
#include "STimer.cc"
#include "PTimer.cc"

#ifdef OMP
    #include <omp.h>
#else
    int omp_get_max_threads() { return 1; }
    int omp_get_num_threads() { return 1; }
    int omp_get_thread_num() { return 0; }
#endif

#include <unistd.h>
#include <string.h>
#include "immintrin.h"

#include "promote_numeric.h"
#include "threevector.hh"
#include "pprint.cc"

#define STDLOG(...) { int a=0; }

#define assertf(_mytest,...) do { \
    if (!(_mytest)) { \
        STDLOG(0,"Failed Assertion: %s\n", #_mytest); STDLOG(1,__VA_ARGS__); \
        fprintf(stderr,"Failed Assertion: %s\n", #_mytest); \
        fpprint(std::cerr, __VA_ARGS__); \
        assert(0==99); \
    }} while(0)




typedef float Float;
typedef unsigned long long uint64;

#include "parallel.partition.cpp"
#include "multiappendlist.cpp"

typedef float3 posstruct;
typedef float3 velstruct;

class auxstruct {
  public:
    uint64 val;
    uint64 pid() { return val; }
};

class Cell {
  public:
    posstruct *pos;
    velstruct *vel;
    auxstruct *aux;
    int n;
    int count() { return n; }
};

class grid {
    public:
    int cpd;
    grid(int &_cpd) { cpd = _cpd;}
    int WrapSlab(int s) {
    	while (s<0) s+=cpd;
	while (s>=cpd) s-=cpd;
	return s;
    }
    integer3 WrapCell(int x, int y, int z) {
        integer3 val;
	val.x = WrapSlab(x);
	val.y = WrapSlab(y);
	val.z = WrapSlab(x);
	return val;
    }
};



