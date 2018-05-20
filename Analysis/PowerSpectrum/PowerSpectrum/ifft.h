/* This is a macro template file meant to be loaded via #include "ifft.h"
 *  Must define:
 *  IFFT_FUNC_NAME:
 *       The function name
 *  FLOAT, FLOAT3:
 *       The data types
 *  arr*:
 *       Ways of accessing arrays
 *
 *  Can define:
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
#define FFTW_PLAN_DFT_C2R_2D fftw_plan_dft_c2r_2d
#define FFTW_PLAN_DFT_C2R_3D fftw_plan_dft_c2r_3d
#define FFTW_MALLOC fftw_malloc
#define FFTW_EXECUTE fftw_execute
#define FFTW_FREE fftw_free
#else
#define FFTW_INIT_THREADS fftwf_init_threads
#define FFTW_PLAN_WITH_NTHREADS fftwf_plan_with_nthreads
#define FFTW_COMPLEX fftwf_complex
#define FFTW_PLAN fftwf_plan
#define FFTW_PLAN_DFT_C2R_2D fftwf_plan_dft_c2r_2d
#define FFTW_PLAN_DFT_C2R_3D fftwf_plan_dft_c2r_3d
#define FFTW_MALLOC fftwf_malloc
#define FFTW_EXECUTE fftwf_execute
#define FFTW_FREE fftwf_free
#endif

#ifdef TWO_D
#define COMPLEX_IN_SHAPE gridN1D*(gridN1D/2 + 1)
#define REAL_OUT_SHAPE gridN1D*gridN1D
#else
#define COMPLEX_IN_SHAPE gridN1D*gridN1D*(gridN1D/2 + 1)
#define REAL_OUT_SHAPE gridN1D*gridN1D*gridN1D
#endif

// On return, input.view(dtype=float)[...,:gridN1D] will be filled with the real transform result
// The transform is always in-place, because the complex arrays always have enough space to do so
void IFFT_FUNC_NAME(FFTW_COMPLEX *input, uint64_t gridN1D){
    // Set up the FFT
    FFTW_INIT_THREADS();
    FFTW_PLAN_WITH_NTHREADS(omp_get_max_threads());
    
    // Need to allocate memory for FFTW planning
    // Could remove this and switch to FFTW_ESTIMATE if we're running out of memory
    FFTW_COMPLEX *wisdom_scratch_in = (FFTW_COMPLEX *) malloc(sizeof(FFTW_COMPLEX)*COMPLEX_IN_SHAPE);  // At least as big as 'input'
    assert(wisdom_scratch_in != NULL);
    
    FLOAT *output = (FLOAT *) input;
    FLOAT *wisdom_scratch_out = (FLOAT *) wisdom_scratch_in;
    
#ifdef TWO_D
    // Create wisdom if it doesn't exist
    FFTW_PLAN_DFT_C2R_2D(gridN1D,gridN1D,wisdom_scratch_in,wisdom_scratch_out,FFTW_MEASURE);
    FFTW_PLAN plan = FFTW_PLAN_DFT_C2R_2D(gridN1D,gridN1D,input,output,FFTW_WISDOM_ONLY);
#else
    // Create wisdom if it doesn't exist
    FFTW_PLAN_DFT_C2R_3D(gridN1D,gridN1D,gridN1D,wisdom_scratch_in,wisdom_scratch_out,FFTW_MEASURE);
    FFTW_PLAN plan = FFTW_PLAN_DFT_C2R_3D(gridN1D,gridN1D,gridN1D,input,output,FFTW_WISDOM_ONLY);
#endif
    assert(plan != NULL);  // Wisdom creation failed?
    free(wisdom_scratch_in);

    // Do the FFT
    //feenableexcept(FE_INVALID | FE_DIVBYZERO);
    FFTW_EXECUTE(plan);

    return;
    //copy back the results and divide out the window function
#ifdef TWO_D
    FLOAT norm = 1./(1.0*gridN1D*gridN1D);
    #define OUT(i,j) output[2*(gridN1D/2 + 1)*i + j]

    #pragma omp parallel for schedule(dynamic,1)
    for(uint64_t i = 0; i < gridN1D; i++){
        for(uint64_t j = 0; j < gridN1D; j++){
            OUT(i,j) *= norm;
        }
    }
#else
    FLOAT norm = 1./(1.0*gridN1D*gridN1D*gridN1D);
    #define OUT(i,j,k) output[gridN1D*2*(gridN1D/2 + 1)*i + 2*(gridN1D/2 + 1)*j + k]
    
    #pragma omp parallel for schedule(dynamic,1)
    for(uint64_t i = 0; i < gridN1D; i++){
        for(uint64_t j = 0; j < gridN1D; j++){
            for(uint64_t k = 0; k < gridN1D; k++){
                OUT(i,j,k) *= norm;
            }
        }
    }
#endif
}

#undef OUT
#undef COMPLEX_IN_SHAPE
#undef REAL_OUT_SHAPE

#undef FFTW_INIT_THREADS 
#undef FFTW_PLAN_WITH_NTHREADS 
#undef FFTW_COMPLEX 
#undef FFTW_PLAN 
#undef FFTW_PLAN_DFT_C2R_2D 
#undef FFTW_PLAN_DFT_C2R_3D 
#undef FFTW_MALLOC 
#undef FFTW_EXECUTE 
#undef FFTW_FREE 

