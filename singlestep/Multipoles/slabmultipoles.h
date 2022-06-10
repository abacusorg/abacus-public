#ifndef SLABMULTIPOLES
#define SLABMULTIPOLES

#include "fftw3.h"
#include "EvaluateMultipoles.h"
#include "redlack.h"

class SlabMultipoles : public Multipoles, public Redlack { 
public:
    
    SlabMultipoles(int order, int cpd); 
    virtual ~SlabMultipoles(void) = 0;

    virtual void ComputeMultipoleFFT( int x,  FLOAT3 *spos, int *count, int *offset,  FLOAT3 *cc, MTCOMPLEX *tmp)
        = 0;
    
    // no-ops in 1D
    virtual void ComputeFFTZ(int x, MTCOMPLEX *outslab) { };
    virtual void CheckAnyMPIDone() { };
    virtual int IsMPIDone(int slab) { return 1; };

    STimer FFTMultipole;
    STimer MultipoleC2R;
    STimer MultipoleKernel;
    STimer ConstructOffsets;
    PTimer _kernel, _c2r, _ffty, _cellmultipoles;

    // 2D
    STimer FFTZMultipole;
    STimer FFTZTranspose;
    STimer UnpackRecvBuf;
    STimer FillMPIBufs;
    STimer AllToAll;
    STimer CheckMPI;

    int *count;         // Will be allocated to [cpd*cpd].  Number of particles
    int *offset;        // Will be allocated to [cpd*cpd].  Start of particles
    FLOAT3 *cc;         // [cpd*cpd] to hold the cell centers

protected:

    fftw_plan plan_forward_c2c_1d[1024];
    fftw_plan plan_forward_r2c_1d[1024];

    Complex **in_1d;
    Complex **out_1d;
    double **in_r2c;
    Complex **out_r2c;
    double **cartesian;
    Complex *transposetmp;
    int nprocs;

    int rml_cpdp1half_pad;
};

class SlabMultipolesLocal : public SlabMultipoles { 
public:
    
    SlabMultipolesLocal(int order, int cpd); 
    ~SlabMultipolesLocal(void);

    void ComputeMultipoleFFT( int x,  FLOAT3 *spos, int *count, int *offset,  FLOAT3 *cc, MTCOMPLEX *tmp);

private:
    double **rowmultipoles;

    void FFTY(MTCOMPLEX *out, const Complex *in);
    void FFTZ(Complex *out, const double *in);
};

#endif
