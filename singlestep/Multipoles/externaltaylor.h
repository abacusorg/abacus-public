#ifndef EXTERNALTAYLOR
#define EXTERNALTAYLOR

#ifdef AVXMULTIPOLES
#ifndef __D4DECL__
#define __D4DECL__
typedef struct { double v[4]; } d4; 
#endif

extern "C" void DispatchTaylorAVXKernel(int order, d4 * px, d4 * py, d4 * pz, d4 * cx, d4 * cy, d4 * cz, d4 * tx, d4 * ty, d4 * tz, d4 * ax, d4 * ay, d4 * az ); 
#endif

#ifdef AVX512MULTIPOLES
#include "avx512_calls.h"

void DispatchTaylor512Kernel(int order, AVX512_DOUBLES &px, AVX512_DOUBLES &py, AVX512_DOUBLES &pz, AVX512_DOUBLES &cx, AVX512_DOUBLES &cy, AVX512_DOUBLES &cz,
                                 AVX512_DOUBLES *tx, AVX512_DOUBLES *ty, AVX512_DOUBLES *tz, AVX512_DOUBLES &ax, AVX512_DOUBLES &ay, AVX512_DOUBLES &az )
#endif

#ifdef UNROLLEDMULTIPOLES
void DispatchTaylorUnrolledKernel(int order, FLOAT3 *p, int n, double3 center, double3 *Q, float3 *acc);
#endif

#ifdef VSXMULTIPOLES
void DispatchTaylorVSXKernel(int order, FLOAT3 *p, int n, double3 center, double3 *Q, float3 *acc);
#endif


#endif
