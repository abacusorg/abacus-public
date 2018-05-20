/* This is a macro template file meant to be loaded via #include "fft.h"
 *  Must define:
 *  FFT_FUNC_NAME:
 *       The function name
 *  FLOAT, FLOAT3:
 *       The data types
 *  den, denR, den_2D, denR_2D; arr*:
 *       Ways of accessing arrays
 *
 *  Can define:
 *  WINDOW:
 *      Apply the window function
 *  INPLACE:
 *      Do the FFT in-place
 *  TWO_D:
 *      Do a 2D FFT instead of 3D
 *  DOUBLEPRECISION:
 *      Use the fftw library instead of fftwf
 */

// Two sets of functions: fftw, and fftwf
#ifdef DOUBLEPRECISION
#define FFTW_INIT_THREADS fftw_init_threads
#define FFTW_PLAN_WITH_NTHREADS fftw_plan_with_nthreads
#define FFTW_COMPLEX fftw_complex
#define FFTW_PLAN fftw_plan
#define FFTW_PLAN_DFT_R2C_2D fftw_plan_dft_r2c_2d
#define FFTW_PLAN_DFT_R2C_3D fftw_plan_dft_r2c_3d
#define FFTW_MALLOC fftw_malloc
#define FFTW_EXECUTE fftw_execute
#define FFTW_FREE fftw_free
#else
#define FFTW_INIT_THREADS fftwf_init_threads
#define FFTW_PLAN_WITH_NTHREADS fftwf_plan_with_nthreads
#define FFTW_COMPLEX fftwf_complex
#define FFTW_PLAN fftwf_plan
#define FFTW_PLAN_DFT_R2C_2D fftwf_plan_dft_r2c_2d
#define FFTW_PLAN_DFT_R2C_3D fftwf_plan_dft_r2c_3d
#define FFTW_MALLOC fftwf_malloc
#define FFTW_EXECUTE fftwf_execute
#define FFTW_FREE fftwf_free
#endif

#ifdef TWO_D
#define COMPLEX_OUT_SHAPE gridN1D*(gridN1D/2+1)
#else
#define COMPLEX_OUT_SHAPE gridN1D*gridN1D*(gridN1D/2+1)
#endif

// On return, density[...,:gridN1D/2+1] will be filled with the square amplitudes
// of the complex scalar field (the "delta tilde's")
void FFT_FUNC_NAME(FLOAT *density, uint64_t gridN1D, double boxsize){
    // Set up the FFT
    FFTW_INIT_THREADS();
    FFTW_PLAN_WITH_NTHREADS(omp_get_max_threads());
    
    // Need to allocate memory for FFTW planning
    // Could remove this and switch to FFTW_ESTIMATE if we're running out of memory
    FLOAT *wisdom_scratch_in = (FLOAT *) malloc(sizeof(FLOAT)*2*COMPLEX_OUT_SHAPE);  // Same size as density
    assert(wisdom_scratch_in != NULL);
    
#ifdef INPLACE
    // If doing in-place transform, the real density will already be padded to (n1d/2 + 1)*2 on the last axis
    FFTW_COMPLEX *output = (FFTW_COMPLEX *) density;
    FFTW_COMPLEX *wisdom_scratch_out = (FFTW_COMPLEX *) wisdom_scratch_in;
#else
    // If doing out-of-place, need to allocate a complex output of size (n1d/2 + 1) on the last axis
    FFTW_COMPLEX *output = (FFTW_COMPLEX *) FFTW_MALLOC(sizeof(FFTW_COMPLEX)*COMPLEX_OUT_SHAPE);
    FFTW_COMPLEX *wisdom_scratch_out = output;
    assert(output != NULL);
#endif
    
#ifdef TWO_D
    // Create wisdom if it doesn't exist
    FFTW_PLAN_DFT_R2C_2D(gridN1D,gridN1D,wisdom_scratch_in,wisdom_scratch_out,FFTW_MEASURE);
    FFTW_PLAN plan = FFTW_PLAN_DFT_R2C_2D(gridN1D,gridN1D,density,output,FFTW_WISDOM_ONLY);
#else
    // Create wisdom if it doesn't exist
    //FFTW_PLAN_DFT_R2C_3D(gridN1D,gridN1D,gridN1D,wisdom_scratch_in,wisdom_scratch_out,FFTW_MEASURE);
    FFTW_PLAN plan = FFTW_PLAN_DFT_R2C_3D(gridN1D,gridN1D,gridN1D,density,output,FFTW_ESTIMATE);
#endif
    assert(plan != NULL);  // Wisdom creation failed?

    // Do the FFT
    //feenableexcept(FE_INVALID | FE_DIVBYZERO);
    FFTW_EXECUTE(plan);

    //copy back the results and divide out the window function
#ifdef TWO_D
    FLOAT norm = boxsize*boxsize/(1.0*gridN1D*gridN1D*gridN1D*gridN1D);
    #ifdef INPLACE
        #define DEN denR_2D
    #else
        #define DEN den_2D
    #endif
    #define OUT(i,j) reinterpret_cast<__complex__ FLOAT*>(output)[(gridN1D/2+1)*i + j]

    #pragma omp parallel for schedule(dynamic,1)
    for(uint64_t i = 0; i < gridN1D; i++){
        for(uint64_t j = 0; j < gridN1D/2+1; j++){
            DEN(i,j) = squ(CABS(OUT(i,j)));  // Square norm of the complex density

            //normalize the power spectrum
            #ifdef WINDOW
            FLOAT3 kv;
            kv.x = (FLOAT)(i)/gridN1D * PI;
            kv.y = (FLOAT)(j)/gridN1D * PI;
            
            FLOAT win = window(kv.x)*window(kv.y);
            DEN(i,j) *= norm/win;
            #else
            DEN(i,j) *= norm;
            #endif
        }
    }
    
#else
    FLOAT norm = boxsize*boxsize*boxsize/(1.0*gridN1D*gridN1D*gridN1D*gridN1D*gridN1D*gridN1D);
    #ifdef INPLACE
        #define DEN denR
    #else
        #define DEN den
    #endif
    #define OUT(i,j,k) reinterpret_cast<__complex__ FLOAT*>(output)[gridN1D*(gridN1D/2+1)*i + (gridN1D/2+1)*j + k]
    
    #pragma omp parallel for schedule(dynamic,1)
    for(uint64_t i = 0; i < gridN1D; i++){
        for(uint64_t j = 0; j < gridN1D; j++){
            for(uint64_t k = 0; k < gridN1D/2+1; k++){
                DEN(i,j,k) = squ(CABS(OUT(i,j,k)));

                //normalize the power spectrum
                #ifdef WINDOW
                FLOAT3 kv;
                kv.x = (FLOAT)(i)/gridN1D * PI;
                kv.y = (FLOAT)(j)/gridN1D * PI;
                kv.z = (FLOAT)(k)/gridN1D * PI;
                
                FLOAT win = window(kv.x)*window(kv.y)*window(kv.z);
                DEN(i,j,k) *= norm/win;
                #else
                DEN(i,j,k) *= norm;
                #endif
            }
        }
    }
#endif
    
free(wisdom_scratch_in);

#ifndef INPLACE
    FFTW_FREE(output);
#endif
}

#undef DEN
#undef OUT
#undef COMPLEX_OUT_SHAPE

#undef FFTW_INIT_THREADS 
#undef FFTW_PLAN_WITH_NTHREADS 
#undef FFTW_COMPLEX 
#undef FFTW_PLAN 
#undef FFTW_PLAN_DFT_R2C_2D 
#undef FFTW_PLAN_DFT_R2C_3D 
#undef FFTW_MALLOC 
#undef FFTW_EXECUTE 
#undef FFTW_FREE 

