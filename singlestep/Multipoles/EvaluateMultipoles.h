#ifndef EVALUATEMULTIPOLES
#define EVALUATEMULTIPOLES

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
    void      AVXCartesianMultipoles(FLOAT3 *p, int n, FLOAT3 center, double *CM);
    void   AVX512CartesianMultipoles(FLOAT3 *p, int n, FLOAT3 center, double *CM);
    void UnrolledCartesianMultipoles(FLOAT3 *p, int n, FLOAT3 center, double *CM);
    void      VSXCartesianMultipoles(FLOAT3 *p, int n, FLOAT3 center, double *CM);
};

#endif
