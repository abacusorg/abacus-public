#include "header.cpp"
#include "threevector.hh"
#include "STimer.h"
#include "PTimer.h"
#include "basemultipoles.h"
#include "slabmultipoles.h"
#include "EvaluateMultipoles.cpp"
#include "redlack.cpp"

SlabMultipolesLocal::~SlabMultipolesLocal(void){
    free(transposetmp);

    for(int g = 0; g < nprocs; g++)
        free(rowmultipoles[g]);
    delete[] rowmultipoles;
}

SlabMultipolesLocal::SlabMultipolesLocal(int order, int cpd)
    : SlabMultipoles(order, cpd) {
    
    rowmultipoles   = new double*[nprocs];
    assert(posix_memalign((void **) &transposetmp, PAGE_SIZE, sizeof(Complex) * rml_cpdp1half_pad*cpd) == 0);

    #pragma omp parallel for schedule(static)
    for(int g = 0; g < nprocs; g++){
        assert(posix_memalign((void **) &(rowmultipoles[g]), PAGE_SIZE, sizeof(double)*rml*cpd) == 0);
        memset(rowmultipoles[g], 0, sizeof(double)*rml*cpd);
    }
}

SlabMultipoles::~SlabMultipoles(void) {
    for(int g=0;g<nprocs;g++) {
        fftw_destroy_plan(plan_forward_c2c_1d[g]);
        fftw_destroy_plan(plan_forward_r2c_1d[g]);
        
        free(in_1d[g]);
        free(out_1d[g]);
        free(in_r2c[g]);
        free(out_r2c[g]);
        
        free(cartesian[g]);
    }


    free(in_1d);
    free(out_1d);
    free(in_r2c);
    free(out_r2c);

    free(cartesian);
    free(count);
    free(offset);
    free(cc);

}

SlabMultipoles::SlabMultipoles(int order, int cpd) : Multipoles(order), 
                                                     Redlack(cpd)  {

    nprocs      = omp_get_max_threads();
    // rml*cpdp1half, padded to page
    rml_cpdp1half_pad = ((cpdp1half*rml + PAGE_SIZE-1)/PAGE_SIZE)*PAGE_SIZE;
    
    in_1d           = (Complex **) malloc(sizeof(Complex*) * nprocs);
    out_1d          = (Complex **) malloc(sizeof(fftw_complex*) * nprocs);
    in_r2c          = (double **)  malloc(sizeof(double*) * nprocs);
    out_r2c         = (Complex **) malloc(sizeof(Complex*) * nprocs);
    
    cartesian       = (double **)  malloc(sizeof(double*)*nprocs);

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

        assert(posix_memalign((void **) &out_1d[g], PAGE_SIZE, sizeof(fftw_complex) * cpd) == 0);
        memset(out_1d[g], 0, sizeof(fftw_complex) * cpd);
        
        assert(posix_memalign((void **) &in_r2c[g], PAGE_SIZE, sizeof(double) * cpd) == 0);
        memset(in_r2c[g], 0, sizeof(double) * cpd);

        assert(posix_memalign((void **) &out_r2c[g], PAGE_SIZE, sizeof(Complex) * (cpd+1)/2) == 0);
        memset(out_r2c[g], 0, sizeof(Complex) * (cpd+1)/2);

        assert(posix_memalign((void **) &(cartesian[g]), PAGE_SIZE, sizeof(double)*cml) == 0);
        memset(cartesian[g], 0, sizeof(double)*cml);
    }

    // FFTW planning is not thread safe; this loop must remain serial!
    for(int g=0;g<nprocs;g++) {
        plan_forward_c2c_1d[g] = 
            fftw_plan_dft_1d(cpd, (fftw_complex *) in_1d[g], 
                                  (fftw_complex *) out_1d[g], 
                                  FFTW_FORWARD, FFTW_PATIENT);
        
        plan_forward_r2c_1d[g] = 
            fftw_plan_dft_r2c_1d(cpd, (double *) in_r2c[g], 
                                      (fftw_complex *) out_r2c[g], 
                                      FFTW_PATIENT);
    }
}

void SlabMultipolesLocal::FFTY(MTCOMPLEX *out, const Complex *in) {
    // out: [(cpd+1)/2, rml, cpd]
    // in: [cpd, rml*(cpd+1)/2 + pad]
    
    FFTMultipole.Start();
    #pragma omp parallel for schedule(static)
    for(int z=0;z<cpdp1half;z++){
        int g = omp_get_thread_num();
        for(int m=0;m<rml;m++) {
            for(int y=0;y<cpd;y++) 
                in_1d[g][y] = in[y*rml_cpdp1half_pad + m*cpdp1half + z];
            fftw_execute(plan_forward_c2c_1d[g]);

            // Here we demote complex<double> to complex<float>
            for(int y=0;y<cpd;y++) 
                out[ z*cpd*rml + m*cpd + y] = (MTCOMPLEX) out_1d[g][y];
        }
    }
    FFTMultipole.Stop();
}

void SlabMultipolesLocal::FFTZ(Complex *out, const double *in) {
    // Execute the FFT along the z dimension (for a single y).
    // Intended to be called from an OpenMP parallel region.
    // out: shape [rml,(cpd+1)/2]
    // in: shape [cpd, rml]

    int g = omp_get_thread_num();
    for(int m=0;m<rml;m++) {
        for(int z=0; z<cpd; z++)
            in_r2c[g][z] = in[ z*rml + m ];
        
        fftw_execute( plan_forward_r2c_1d[g] );

        // We store the data non-transposed in this out-of-place buffer,
        // and then execute the transpose in FFTY
        for(int z=0;z<cpdp1half;z++) 
            out[m*cpdp1half + z] = out_r2c[g][z];
    }
}

void SlabMultipolesLocal::ComputeMultipoleFFT( int x, FLOAT3 *spos, 
                     int *count, int *offset, FLOAT3 *cc, MTCOMPLEX *out) {
    STimer wc;
    PTimer _kernel, _c2r, _fftz;
    pdouble localMassSlabX[nprocs];
    padded<double3> localdipole[nprocs];
    double *localMassSlabZ = new double[nprocs*cpd];
    
    // Init with OpenMP to preserve thread locality
    #pragma omp parallel for schedule(static)
    for(int g = 0; g < nprocs; g++){
        localMassSlabX[g] = 0.;
        localdipole[g] = double3(0.);
        for(int z = 0; z < cpd; z++)
            localMassSlabZ[g*cpd + z] = 0.;
    }
    
    wc.Start();
    NUMA_FOR(y,0,cpd)
        int g = omp_get_thread_num();
        for(int z=0;z<cpd;z++) {
            _kernel.Start();
            int i = y*cpd + z;
            EvaluateCartesianMultipoles( &(spos[offset[i]]),
                                    count[i], cc[i], cartesian[g] );
            _kernel.Stop();
            
            _c2r.Start();
            DispatchCartesian2Reduced(order, cartesian[g], &(rowmultipoles[g][z*rml]));
            double Mxyz = rowmultipoles[g][z*rml];
            localMassSlabX[g] += Mxyz;
            localMassSlabZ[g*cpd + z] += Mxyz;
            MassSlabY[y] += Mxyz;  // thread dimension, no race TODO: still false sharing
            
            localdipole[g] += double3(rowmultipoles[g][z*rml + rmap(1,0,0) ],
                                rowmultipoles[g][z*rml + rmap(0,1,0) ],
                                rowmultipoles[g][z*rml + rmap(0,0,1) ] );
            _c2r.Stop();
        }
        _fftz.Start();
        FFTZ( &(transposetmp[y*rml_cpdp1half_pad]), rowmultipoles[g]);
        _fftz.Stop();
    }
    
    for(int g = 0; g < nprocs; g++){
        MassSlabX[x] += localMassSlabX[g];
        globaldipole += localdipole[g];
        for(int z = 0; z < cpd; z++)
            MassSlabZ[z] += localMassSlabZ[g*cpd + z];
    }
    delete[] localMassSlabZ;
    
    wc.Stop();
    FFTY(out, transposetmp);
    
    double seq = _kernel.Elapsed() + _c2r.Elapsed() + _fftz.Elapsed();
    
    struct timespec  seq_kernel = scale_timer(_kernel.Elapsed()/seq, wc.get_timer() );
    struct timespec  seq_c2r = scale_timer(_c2r.Elapsed()/seq, wc.get_timer() );
    struct timespec  seq_fftz = scale_timer(_fftz.Elapsed()/seq, wc.get_timer() );

    MultipoleKernel.increment(seq_kernel);
    MultipoleC2R.increment(seq_c2r);
    FFTMultipole.increment(seq_fftz);
}
