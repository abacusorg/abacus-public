#ifndef SLABTAYLOR
#define SLABTAYLOR

#include "fftw3.h"
#include "EvaluateTaylor.h"

class SlabTaylor : public Taylor {
public:
    SlabTaylor(int order, int cpd);
    ~SlabTaylor(void);

    void InverseFFTY_Taylor( MTCOMPLEX *STT);
    void InverseFFTZ_Taylor( int y, MTCOMPLEX *STT);

    STimer FFTTaylor;
    STimer TaylorR2C;
    STimer TaylorKernel;
    STimer ConstructOffsets;

    void CellTaylorFromPencilTaylor(int z,   double *T);
    void EvaluateSlabTaylor(int x,  FLOAT3 *FA,  FLOAT3 *spos, int *count, int *offset,  FLOAT3 *cc, MTCOMPLEX *TaylorCoefficients);

    int *count; 	// Will be allocated to [cpd*cpd].  Number of particles
    int *offset; 	// Will be allocated to [cpd*cpd].  Start of particles
    FLOAT3 *cc; 	// [cpd*cpd] to hold the cell centers

private:

    double **TaylorPencil; double *tfactor; double **cartesian;
    Complex **in_1d; Complex **out_1d;
    Complex **in_c2r; double  **out_c2r;

    fftw_plan plan_backward_c2r_1d[128];
    fftw_plan plan_backward_c2c_1d[128];
    int cpd;
    int nprocs;
    double invcpd3;
};

#endif
