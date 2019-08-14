#ifndef TAYLOR
#define TAYLOR

#ifdef AVX512MULTIPOLES
#include "avx512_calls.h"
#endif

#include "externaltaylor.h"
#include "basemultipoles.h"


class Taylor : public basemultipoles {
public:
    Taylor(int order);
    ~Taylor(void);

    void         EvaluateTaylor(double *CT, FLOAT3 expansioncenter, int np, FLOAT3 *ps, FLOAT3 *a);
    void AnalyticEvaluateTaylor(double *CT, FLOAT3 expansioncenter, int np, FLOAT3 *ps, FLOAT3 *a);
    void      AVXEvaluateTaylor(double *CT, FLOAT3 expansioncenter, int np, FLOAT3 *ps, FLOAT3 *a);
    void   AVX512EvaluateTaylor(double *CT, FLOAT3 expansioncenter, int np, FLOAT3 *ps, FLOAT3 *a);
    void UnrolledEvaluateTaylor(double *CT, FLOAT3 expansioncenter, int np, FLOAT3 *ps, FLOAT3 *a);
    void      VSXEvaluateTaylor(double *CT, FLOAT3 expansioncenter, int np, FLOAT3 *ps, FLOAT3 *a);

#define MAXTHREADS 128
#ifdef AVXMULTIPOLES
    d4 *cx[MAXTHREADS];
    d4 *cy[MAXTHREADS];
    d4 *cz[MAXTHREADS];

    d4 *ax[MAXTHREADS];
    d4 *ay[MAXTHREADS];
    d4 *az[MAXTHREADS];

    d4 *px[MAXTHREADS];
    d4 *py[MAXTHREADS];
    d4 *pz[MAXTHREADS];

    d4 *Qx[MAXTHREADS];
    d4 *Qy[MAXTHREADS];
    d4 *Qz[MAXTHREADS];
#endif

};

#endif
