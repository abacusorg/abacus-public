// This is a code to provide a test for the group finding.

// There are a couple of other unit tests:
// g++ -DTEST -fopenmp -lgomp -O2 fof_sublist.cpp 
// g++ -DTEST -fopenmp -lgomp -O2 slab_accum.cpp

#define STANDALONE_FOF

#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <algorithm>
#include <vector>
#include <cstdint>
#include <assert.h>
#include "stdlog.cc"

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

/* #define STDLOG(...) { int a=0; }

#define assertf(_mytest,...) do { \
    if (!(_mytest)) { \
        STDLOG(0,"Failed Assertion: %s\n", #_mytest); STDLOG(1,__VA_ARGS__); \
        fprintf(stderr,"Failed Assertion: %s\n", #_mytest); \
        fpprint(std::cerr, __VA_ARGS__); \
        assert(0==99); \
    }} while(0)

*/



typedef float Float;
typedef float3 Float3;
typedef unsigned long long uint64;

#include "parallel.partition.cpp"
#include "multiappendlist.cpp"

typedef float3 posstruct;
typedef float3 velstruct;
typedef float3 accstruct;

#define AUXTAGGABLEBIT 48llu //Can the particle be tagged.
#define AUXTAGGEDBIT 49llu //Has the particle been tagged
#define AUXL0BIT 50llu //Is the particle in a level 0 group
#define AUXL1BIT 51llu //Is the particle in a levl 1 group
class auxstruct {
  public:
    uint64 val;
    uint64 pid() { return val; }

    // Group and subsample related bits
    inline void set_taggable() {
        // The TAGGABLE bit should be set at the beginning of the sim and not changed.
        val |= ((uint64)1 << AUXTAGGABLEBIT);
    }
    inline bool is_taggable() {
        return val & ((uint64)1 << AUXTAGGABLEBIT);
    }
    inline void set_tagged() {
        // The TAGGED bit is a lasting tag, once set.
        val |= ((uint64)1 << AUXTAGGEDBIT);
    }
    inline bool is_tagged() {
        return val & ((uint64)1 << AUXTAGGEDBIT);
    }

    inline void reset_L01_bits() {
        // We need to be able to unset these bits each time we run groupfinding
        uint64 mask = ((uint64)1 << AUXL0BIT) + ((uint64)1 << AUXL1BIT);
        val &= ~mask;
    }
    inline void set_L0() {
        val |= ((uint64)1 << AUXL0BIT);
    }
    inline bool is_L0() {
        return val & ((uint64)1 << AUXL0BIT);
    }
    inline void set_L1() {
        val |= ((uint64)1 << AUXL1BIT);
    }
    inline bool is_L1() {
        return val & ((uint64)1 << AUXL1BIT);
    }
};

class Cell {
  public:
    posstruct *pos;
    velstruct *vel;
    auxstruct *aux;
    accstruct *acc;
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

void setup_log() {
    stdlog_threshold_global = 1;
    stdlog.open("/tmp/fof.log");
    STDLOG_TIMESTAMP;
}
