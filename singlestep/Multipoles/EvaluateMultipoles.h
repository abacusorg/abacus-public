#ifndef MULTIPOLES
#define MULTIPOLES

#include "basemultipoles.h"

#ifdef AVX512MULTIPOLES
#include "avx512_calls.c"
#endif

#include "externalmultipoles.h"

class Multipoles : public basemultipoles { 
public:
    Multipoles(int order);
    ~Multipoles();

#ifdef AVXMULTIPOLES
    d4 *ip1x[128];
    d4 *ip2x[128];

    d4 *ip1y[128];
    d4 *ip2y[128];

    d4 *ip1z[128];
    d4 *ip2z[128];

    d4 *cx[128];
    d4 *cy[128];
    d4 *cz[128];

    d4 *masses1[128];
    d4 *masses2[128];

    d4 *globalM[128];
#endif

    void EvaluateCartesianMultipoles(FLOAT3 *p, int n, FLOAT3 center, double *CM);
    void AnalyticCartesianMultipoles(FLOAT3 *p, int n, FLOAT3 center, double *CM);
    void      ASMCartesianMultipoles(FLOAT3 *p, int n, FLOAT3 center, double *CM);
    void   AVX512CartesianMultipoles(FLOAT3 *p, int n, FLOAT3 center, double *CM);
    void UnrolledCartesianMultipoles(FLOAT3 *p, int n, FLOAT3 center, double *CM);
    void      VSXCartesianMultipoles(FLOAT3 *p, int n, FLOAT3 center, double *CM);
};

// The following are the function pointers to the automatically generated kernels
// e.g. the AVX512 kernels from generateCartesianAVX512.py

#define REP16(K) { \
     K<1>,  K<2>,  K<3>,  K<4>, \
     K<5>,  K<6>,  K<7>,  K<8>, \
     K<9>,  K<10>, K<11>, K<12>, \
     K<13>, K<14>, K<15>, K<16> \
};

#ifdef AVXMULTIPOLES 
void (*CMptr[24])( d4 *ip1x, d4 *ip2x, d4 *ip1y, d4 *ip2y, d4 *ip1z, d4 *ip2z, 
                   d4 *cx, d4 *cy, d4 *cz, d4 *globalM, 
                   d4 *mass1, d4 *mass2) = {
     MultipoleKernel1,  MultipoleKernel2,  MultipoleKernel3,  MultipoleKernel4,
     MultipoleKernel5,  MultipoleKernel6,  MultipoleKernel7,  MultipoleKernel8,
     MultipoleKernel9,  MultipoleKernel10, MultipoleKernel11, MultipoleKernel12,
     MultipoleKernel13, MultipoleKernel14, MultipoleKernel15, MultipoleKernel16
};
#endif

#ifdef AVX512MULTIPOLES 
void (*CM512ptr[24])( AVX512_DOUBLES &px, AVX512_DOUBLES &py, AVX512_DOUBLES &pz,
                      AVX512_DOUBLES &cx, AVX512_DOUBLES &cy, AVX512_DOUBLES &cz,
                      AVX512_DOUBLES *CM ) = REP16(Multipole512Kernel);
#endif

#ifdef UNROLLEDMULTIPOLES 
void (*CM_unrolled_ptr[24])(FLOAT3 *p, int n, double3 center, double *CM) = REP16(MultipoleUnrolledKernel);
#endif

#ifdef VSXMULTIPOLES 
void (*CM_VSX_ptr[24])(FLOAT3 *p, int n, double3 center, double *CM) = REP16(MultipoleVSXKernel);
#endif

#undef REP16

#endif
