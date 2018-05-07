#ifndef TAYLOR
#define TAYLOR

#ifdef AVX512MULTIPOLES
#include "avx512_calls.h"
#endif

#ifdef AVXMULTIPOLES
#include "externaltaylor.h"
#endif

#include "basemultipoles.h"


class Taylor : public basemultipoles {
public:
    Taylor(int order);
    ~Taylor(void);

    void         EvaluateTaylor(double *CT, FLOAT3 expansioncenter, int np, FLOAT3 *ps, FLOAT3 *a);
    void AnalyticEvaluateTaylor(double *CT, FLOAT3 expansioncenter, int np, FLOAT3 *ps, FLOAT3 *a);
    void      ASMEvaluateTaylor(double *CT, FLOAT3 expansioncenter, int np, FLOAT3 *ps, FLOAT3 *a);
    void   AVX512EvaluateTaylor(double *CT, FLOAT3 expansioncenter, int np, FLOAT3 *ps, FLOAT3 *a);

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

/*
#ifdef AVX512MULTIPOLES
    AVX512_DOUBLES *cx512[MAXTHREADS];
    AVX512_DOUBLES *cy512[MAXTHREADS];
    AVX512_DOUBLES *cz512[MAXTHREADS];

    AVX512_DOUBLES *ax512[MAXTHREADS];
    AVX512_DOUBLES *ay512[MAXTHREADS];
    AVX512_DOUBLES *az512[MAXTHREADS];

    AVX512_DOUBLES *px512[MAXTHREADS];
    AVX512_DOUBLES *py512[MAXTHREADS];
    AVX512_DOUBLES *pz512[MAXTHREADS];

    AVX512_DOUBLES *Qx512[MAXTHREADS];
    AVX512_DOUBLES *Qy512[MAXTHREADS];
    AVX512_DOUBLES *Qz512[MAXTHREADS];
#endif
*/
};

#endif
