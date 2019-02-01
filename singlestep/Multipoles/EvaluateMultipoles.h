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

    void EvaluateCartesianMultipoles(FLOAT3 *p, int n, FLOAT3 expansioncenter, double *cm);
    void AnalyticCartesianMultipoles(FLOAT3 *p, int n, FLOAT3 expansioncenter, double *cm);
    void      ASMCartesianMultipoles(FLOAT3 *p, int n, FLOAT3 expansioncenter, double *cm);
    void   AVX512CartesianMultipoles(FLOAT3 *p, int n, FLOAT3 expansioncenter, double *cm);
    void UnrolledCartesianMultipoles(FLOAT3 *p, int n, FLOAT3 expansioncenter, double *cm);

};

// The following are the function pointers to the automatically generated kernels
// e.g. the AVX512 kernels from generateCartesianAVX512.py

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
                      AVX512_DOUBLES *CM ) = {
     Multipole512Kernel1,  Multipole512Kernel2,  Multipole512Kernel3,  Multipole512Kernel4,
     Multipole512Kernel5,  Multipole512Kernel6,  Multipole512Kernel7,  Multipole512Kernel8,
     Multipole512Kernel9,  Multipole512Kernel10, Multipole512Kernel11, Multipole512Kernel12,
     Multipole512Kernel13, Multipole512Kernel14, Multipole512Kernel15, Multipole512Kernel16
};
#endif

#ifdef UNROLLEDMULTIPOLES 
void (*CM_unrolled_ptr[24])(FLOAT3 p, double3 center, double *CM) = {
     MultipoleUnrolledKernel1,  MultipoleUnrolledKernel2,  MultipoleUnrolledKernel3,  MultipoleUnrolledKernel4,
     MultipoleUnrolledKernel5,  MultipoleUnrolledKernel6,  MultipoleUnrolledKernel7,  MultipoleUnrolledKernel8,
     MultipoleUnrolledKernel9,  MultipoleUnrolledKernel10, MultipoleUnrolledKernel11, MultipoleUnrolledKernel12,
     MultipoleUnrolledKernel13, MultipoleUnrolledKernel14, MultipoleUnrolledKernel15, MultipoleUnrolledKernel16
};
#endif

#endif
