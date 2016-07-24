#ifndef MULTIPOLES
#define MULTIPOLES

#include "basemultipoles.h"

#ifdef AVXMULTIPOLES
#include "externalmultipoles.h"
#endif


class Multipoles : public  basemultipoles { 
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

    void AnalyticCartesianMultipoles(double3 *p, int n, double3 expansioncenter, double *cm);
    void      ASMCartesianMultipoles(double3 *p, int n, double3 expansioncenter, double *cm);
    
    void AnalyticCartesianMultipoles( float3 *p, int n,  float3 expansioncenter, double *cm);
    void      ASMCartesianMultipoles( float3 *p, int n,  float3 expansioncenter, double *cm);

};

#endif
