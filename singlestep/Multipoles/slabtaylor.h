#ifndef SLABTAYLOR
#define SLABTAYLOR

#include "fftw3.h"
#include "EvaluateTaylor.h"

class SlabTaylor : public Taylor {
public:
    SlabTaylor(int order, int cpd);
    virtual ~SlabTaylor(void);

    virtual void EvaluateSlabTaylor(int x, FLOAT3 *FA, const FLOAT3 *spos,
                            const int *count, const int *offset, const int *ghost_offsets,
                           const FLOAT3 *cc, const MTCOMPLEX *TaylorCoefficients) = 0;
    
    // no-ops in 1D
    virtual void ComputeIFFTZAndMPI(int x, MTCOMPLEX *outslab) { };
    virtual void LaunchAllToAll() { };
    
    virtual int CheckAnyMPIDone() { return 0; };
    virtual int IsMPIDone(int slab) { return 1; };


    STimer FFTTaylor;
    STimer TaylorR2C;
    STimer TaylorKernel;
    STimer ConstructOffsets;

    // 2D
    STimer UnpackRecvBuf;
    STimer FFTZTaylor;
    STimer FillMPIBufs;
    STimer AllToAll;
    STimer CheckMPI;

    int *count; 	// Will be allocated to [cpd*cpd].  Number of particles
    int *offset; 	// Will be allocated to [cpd*cpd].  Start of particles
    FLOAT3 *cc; 	// [cpd*cpd] to hold the cell centers

protected:
    int cpd;
    int nprocs;
    double invcpd3;
    int cpdp1half;
    int rml_cpd_pad;

    double **TaylorPencil;
    double *tfactor;
    double **cartesian;
    
    Complex **in_1d;
    Complex **out_1d;
    Complex **in_c2r;
    double  **out_c2r;

    fftw_plan plan_backward_c2r_1d[1024];
    fftw_plan plan_backward_c2c_1d[1024];

    void CellTaylorFromPencil(int z,  double *T, const double *pencil);
};

class SlabTaylorLocal : public SlabTaylor {
public:
    SlabTaylorLocal(int order, int cpd);
    ~SlabTaylorLocal(void);

    void EvaluateSlabTaylor(int x, FLOAT3 *FA, const FLOAT3 *spos,
                            const int *count, const int *offset, const int *ghost_offsets,
                            const FLOAT3 *cc, const MTCOMPLEX *TaylorCoefficients);

private:
    void InverseFFTY(Complex *out, const MTCOMPLEX *in);
    void InverseFFTZ(int y, double *out, const Complex *in);

    Complex *transposetmp;
};

#endif
