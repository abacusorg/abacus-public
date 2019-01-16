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
    void      ASMEvaluateTaylor(double *CT, FLOAT3 expansioncenter, int np, FLOAT3 *ps, FLOAT3 *a);
    void   AVX512EvaluateTaylor(double *CT, FLOAT3 expansioncenter, int np, FLOAT3 *ps, FLOAT3 *a);
    void UnrolledEvaluateTaylor(double *CT, FLOAT3 expansioncenter, int np, FLOAT3 *ps, FLOAT3 *a);

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

#ifdef AVXMULTIPOLES
void (*Tptr[17])(d4 * px, d4 * py, d4 * pz, d4 * cx, d4 * cy, d4 * cz, 
                 d4 * tx, d4 * ty, d4 * tz, d4 * ax, d4 * ay, d4 * az) = { 
 TaylorKernel0, TaylorKernel1, TaylorKernel2, TaylorKernel3, TaylorKernel4, 
 TaylorKernel5, TaylorKernel6, TaylorKernel7, TaylorKernel8, 
 TaylorKernel9, TaylorKernel10, TaylorKernel11, TaylorKernel12, 
 TaylorKernel13, TaylorKernel14, TaylorKernel15, TaylorKernel16 }; 
#endif


// We are not using the manually unrolled versions right now
#ifdef AVX512MULTIPOLES
void (*Tptr512[17])(AVX512_DOUBLES &px, AVX512_DOUBLES &py, AVX512_DOUBLES &pz, AVX512_DOUBLES &cx, AVX512_DOUBLES &cy, AVX512_DOUBLES &cz, 
                    AVX512_DOUBLES *tx, AVX512_DOUBLES *ty, AVX512_DOUBLES *tz, AVX512_DOUBLES &ax, AVX512_DOUBLES &ay, AVX512_DOUBLES &az) = { 
 Taylor512Kernel0, Taylor512Kernel1, Taylor512Kernel2, Taylor512Kernel3, Taylor512Kernel4, 
 Taylor512Kernel5, Taylor512Kernel6, Taylor512Kernel7, Taylor512Kernel8, 
 Taylor512Kernel9, Taylor512Kernel10, Taylor512Kernel11, Taylor512Kernel12, 
 Taylor512Kernel13, Taylor512Kernel14, Taylor512Kernel15, Taylor512Kernel16 };
#endif

#ifdef UNROLLEDMULTIPOLES
void (*Tptr_unrolled[17])(FLOAT3 *p, int n, double3 center, double3 *Q, float3 *acc) = { 
 TaylorUnrolledKernel0, TaylorUnrolledKernel1, TaylorUnrolledKernel2, TaylorUnrolledKernel3, TaylorUnrolledKernel4, 
 TaylorUnrolledKernel5, TaylorUnrolledKernel6, TaylorUnrolledKernel7, TaylorUnrolledKernel8, 
 TaylorUnrolledKernel9, TaylorUnrolledKernel10, TaylorUnrolledKernel11, TaylorUnrolledKernel12, 
 TaylorUnrolledKernel13, TaylorUnrolledKernel14, TaylorUnrolledKernel15, TaylorUnrolledKernel16 };
#endif

#endif
