#include "header.cpp"
#include "threevector.hh"
#include "STimer.h"
#include "PTimer.h"
#include "basemultipoles.h"
#include "slabtaylor.h"
#include "EvaluateTaylor.cpp"
#include "factorial.h"

SlabTaylorLocal::~SlabTaylorLocal(void){
    free(transposetmp);

}

SlabTaylorLocal::SlabTaylorLocal(int order, int _cpd)
    : SlabTaylor(order, _cpd) {

    assert(posix_memalign((void **) &transposetmp, PAGE_SIZE, sizeof(Complex)*rml_cpd_pad*cpdp1half) == 0);

}

SlabTaylor::SlabTaylor(int order, int _cpd) : Taylor(order) {
    cpd = _cpd;
    invcpd3 = 1.0/cpd/cpd/cpd;
    nprocs = omp_get_max_threads();
    cpdp1half = (cpd+1)/2;
    rml_cpd_pad = ((cpd*rml + PAGE_SIZE-1)/PAGE_SIZE)*PAGE_SIZE;

    assert(posix_memalign((void **) &tfactor, PAGE_SIZE, sizeof(double)*completemultipolelength) == 0);
    
    int i,j,k;
    FORALL_COMPLETE_MULTIPOLES_BOUND(i,j,k,order)  
        tfactor[cmap(i,j,k)] = pow(-1.0,i+j+k)/(fact(i)*fact(j)*fact(k));

    in_1d           = (Complex **)  malloc(sizeof(fftw_complex*) * nprocs);
    out_1d          = (Complex **)  malloc(sizeof(Complex*) * nprocs);
    in_c2r          = (Complex **)  malloc(sizeof(Complex*) * nprocs);
    out_c2r         = (double **)   malloc(sizeof(double*) * nprocs);
    
    TaylorPencil    = (double **) malloc(sizeof(double*) * nprocs);
    cartesian       = (double **) malloc(sizeof(double*) * nprocs);

    assert(posix_memalign((void **) &cc, PAGE_SIZE, cpd*node_z_size*sizeof(FLOAT3)) == 0);
    assert(posix_memalign((void **) &count, PAGE_SIZE, cpd*node_z_size*sizeof(int)) == 0);
    assert(posix_memalign((void **) &offset, PAGE_SIZE, cpd*node_z_size*sizeof(int)) == 0);

    #pragma omp parallel for schedule(static)
    for(int g = 0; g < nprocs; g++){

        // Normally, one would allocate these with fftw_malloc
        // But we wish to impose a stricter alignment so that each array can land on a separate NUMA node
        // To enforce this, we go ahead and touch the memory with memset from a parallel region
        assert(posix_memalign((void **) &in_1d[g], PAGE_SIZE, sizeof(fftw_complex) * cpd) == 0);
        memset(in_1d[g], 0, sizeof(fftw_complex) * cpd);

        assert(posix_memalign((void **) &out_1d[g], PAGE_SIZE, sizeof(Complex) * cpd) == 0);
        memset(out_1d[g], 0, sizeof(Complex) * cpd);
        
        assert(posix_memalign((void **) &in_c2r[g], PAGE_SIZE, sizeof(Complex) * (cpd+1)/2) == 0);
        memset(in_c2r[g], 0, sizeof(Complex) * (cpd+1)/2);

        assert(posix_memalign((void **) &out_c2r[g], PAGE_SIZE, sizeof(double) * cpd) == 0);
        memset(out_c2r[g], 0, sizeof(double) * cpd);
        
        assert(posix_memalign((void **) &(TaylorPencil[g]), PAGE_SIZE, sizeof(double)*rml*cpd) == 0);
        memset(TaylorPencil[g], 0, sizeof(double)*rml*cpd);

        assert(posix_memalign((void **) &(cartesian[g]), PAGE_SIZE, sizeof(double)*cml) == 0);
        memset(cartesian[g], 0, sizeof(double)*cml);

    }
        
    // FFTW planning is not thread safe; this loop must remain serial!
    for(int g=0;g<nprocs;g++) {
         plan_backward_c2r_1d[g] = fftw_plan_dft_c2r_1d(cpd, 
                                      (fftw_complex *) in_c2r[g], 
                                      (double *) out_c2r[g], 
                                      FFTW_PATIENT);
        
         plan_backward_c2c_1d[g] = fftw_plan_dft_1d(cpd, 
                                      (fftw_complex *) in_1d[g], 
                                      (fftw_complex *) out_1d[g], 
                                      FFTW_BACKWARD, FFTW_PATIENT);
    }
}

SlabTaylor::~SlabTaylor(void) {
    for(int g=0;g<nprocs;g++) {
        fftw_destroy_plan( plan_backward_c2r_1d[g] );
        fftw_destroy_plan( plan_backward_c2c_1d[g] );

        free(in_1d[g]);
        free(out_1d[g]);
        free(in_c2r[g]);
        free(out_c2r[g]);
        
        free(TaylorPencil[g]); free(cartesian[g]);
    }
    free(in_1d); free(out_1d);
    free(in_c2r); free(out_c2r);
    free(TaylorPencil); free(cartesian);
    free(tfactor);
    free(count);
    free(offset);
    free(cc);
}

void SlabTaylorLocal::InverseFFTY(Complex *out, const MTCOMPLEX *in) {
    // in: [cpdp1half, rml, cpd]
    // out: [cpdp1half, rml, cpd]

    FFTTaylor.Start();
    #pragma omp parallel for schedule(static)
    for(int z=0;z<cpdp1half;z++) {
        int g = omp_get_thread_num();
        for(int m=0;m<rml;m++) {
            for(int y=0;y<cpd;y++) in_1d[g][y] = (Complex) in[z*cpd*rml + m*cpd + y];
            fftw_execute( plan_backward_c2c_1d[g] );
            for(int y=0;y<cpd;y++) out[z*rml_cpd_pad + m*cpd + y] = out_1d[g][y];
        }
    }
    FFTTaylor.Stop();
}

void SlabTaylorLocal::InverseFFTZ(int y, double *out, const Complex *in) {
    // in: [cpdp1half, rml, cpd]
    // out: [rml, cpd]

    int g = omp_get_thread_num();
    for(int m=0;m<rml;m++) {
        for(int z=0;z<cpdp1half;z++) 
            in_c2r[g][z] = in[z*rml_cpd_pad + m*cpd + y];
        fftw_execute( plan_backward_c2r_1d[g] );
        for(int z=0;z<cpd;z++) 
            out[m*cpd + z] = out_c2r[g][z];
    }
}

void SlabTaylor::CellTaylorFromPencil(int z,  double *T, const double *pencil) {
    // pencil: [rml, cpd]
    // T: [cml]
    
    int a,b,c;
    FORALL_REDUCED_MULTIPOLES_BOUND(a,b,c,order) 
        T[cmap(a,b,c)] = pencil[ rmap(a,b,c)*cpd + z ];
    TRACEFREERECURSION(a,b,c,order) 
        T[cmap(a,b,c)] = -T[cmap(a+2,b,c-2)] - T[cmap(a,b+2,c-2)];
    for(int m=0;m<completemultipolelength;m++) T[m] *= tfactor[m];
}


void SlabTaylorLocal::EvaluateSlabTaylor(int x, FLOAT3 *FA, const FLOAT3 *spos,
                                        const int *count, const int *offset, const int *_ghost_offset_unused,
                                        const FLOAT3 *cc, const MTCOMPLEX *TaylorCoefficients){
    // TaylorCoefficients: [(cpd+1)/2,m,cpd]
    // FA: particle accelerations

    STimer wc;
    PTimer _r2c, _tkernel, _zfft, _kernel_r2c;

    InverseFFTY(transposetmp, TaylorCoefficients);
    wc.Start();

    NUMA_FOR(y,0,cpd, NO_CLAUSE, FALLBACK_DYNAMIC){
        int g = omp_get_thread_num();
        _zfft.Start();
        InverseFFTZ( y, TaylorPencil[g], transposetmp);
        _zfft.Stop();
        
        _kernel_r2c.Start();
        for(int z=0;z<cpd;z++) {
            int i = y*cpd + z;
            //_r2c.Start();
            CellTaylorFromPencil(z, cartesian[g], TaylorPencil[g]);
            //_r2c.Stop();
            
            //_tkernel.Start();
            FLOAT3 *aa = &(FA[offset[i]]);
            memset(aa, 0, sizeof(FLOAT3)*count[i]);

            EvaluateTaylor( cartesian[g], 
                               cc[i], count[i], (float3*) &spos[offset[i]], aa);
            //_tkernel.Stop();
        }
        _kernel_r2c.Stop();
    }
    NUMA_FOR_END;
    wc.Stop();
    
    /*double seq = _r2c.Elapsed() + _tkernel.Elapsed() + _zfft.Elapsed();
    double f_r2c = _r2c.Elapsed()/seq;
    double f_kernel = _tkernel.Elapsed()/seq;
    double f_zfft = _zfft.Elapsed()/seq;

    struct timespec  seq_r2c = scale_timer(f_r2c, wc.get_timer() );
    struct timespec  seq_tkernel = scale_timer(f_kernel, wc.get_timer() );
    struct timespec  seq_zfft = scale_timer(f_zfft, wc.get_timer() );

    TaylorR2C.increment( seq_r2c  );
    TaylorKernel.increment( seq_tkernel );
    FFTTaylor.increment( seq_zfft );*/

    double seq = _kernel_r2c.Elapsed() + _zfft.Elapsed();
    double f_kernel_r2c = _kernel_r2c.Elapsed()/seq;
    double f_zfft = _zfft.Elapsed()/seq;

    struct timespec  seq_kernel_r2c = scale_timer(f_kernel_r2c, wc.get_timer() );
    struct timespec  seq_zfft = scale_timer(f_zfft, wc.get_timer() );

    TaylorKernel.increment(seq_kernel_r2c);
    FFTTaylor.increment( seq_zfft );
}
