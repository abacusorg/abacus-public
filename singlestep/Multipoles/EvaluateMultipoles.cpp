#ifdef TEST

// Compile test driver with:
// $ make EvaluateMultipolesTest

// Can use these to override ./configure for testing
//#define AVXMULTIPOLES
//#define AVX512MULTIPOLES
#define UNROLLEDMULTIPOLES

#include "threevector.hh"
#define FLOAT float
#define FLOAT3 float3

#include "basemultipoles.cpp"
#endif  // #ifdef TEST

#include "EvaluateMultipoles.h"

Multipoles::Multipoles(int order) : basemultipoles(order) {
#ifdef AVXMULTIPOLES
for(int g=0;g<omp_get_max_threads();g++) {
    int rv;
    rv = posix_memalign( (void **) &(ip1x[g]), 256, 512 ); assert(rv==0);
    rv = posix_memalign( (void **) &(ip2x[g]), 256, 512 ); assert(rv==0);
    rv = posix_memalign( (void **) &(ip1y[g]), 256, 512 ); assert(rv==0);
    rv = posix_memalign( (void **) &(ip2y[g]), 256, 512 ); assert(rv==0);
    rv = posix_memalign( (void **) &(ip1z[g]), 256, 512 ); assert(rv==0);
    rv = posix_memalign( (void **) &(ip2z[g]), 256, 512 ); assert(rv==0);

    rv = posix_memalign( (void **) &(cx[g]),   256, 512 ); assert(rv==0);
    rv = posix_memalign( (void **) &(cy[g]),   256, 512 ); assert(rv==0);
    rv = posix_memalign( (void **) &(cz[g]),   256, 512 ); assert(rv==0);

    rv = posix_memalign( (void **) &(masses1[g]), 256, 512  ); assert(rv==0);
    rv = posix_memalign( (void **) &(masses2[g]), 256, 512  ); assert(rv==0);
    rv = posix_memalign( (void **) &(globalM[g]), 256, 32768); assert(rv==0);
}
#endif
}

Multipoles::~Multipoles(){
#ifdef AVXMULTIPOLES
    for(int g=0;g<omp_get_max_threads();g++) {
        free(ip1x[g]);
        free(ip2x[g]);
        free(ip1y[g]);
        free(ip2y[g]);
        free(ip1z[g]);
        free(ip2z[g]);

        free(cx[g]);
        free(cy[g]);
        free(cz[g]);

        free(masses1[g]);
        free(masses2[g]);
        free(globalM[g]);
    }
#endif
}

// Dispatch function
void Multipoles::EvaluateCartesianMultipoles(FLOAT3 *p, int n, FLOAT3 center, double *cm) {
    #ifdef AVX512MULTIPOLES
    AVX512CartesianMultipoles(p, n, center, cm);
    #elif defined(AVXMULTIPOLES)
    ASMCartesianMultipoles(p, n, center, cm);
    #elif defined(UNROLLEDMULTIPOLES)
    UnrolledCartesianMultipoles(p, n, center, cm);
    #else
    AnalyticCartesianMultipoles(p, n, center, cm);
    #endif
}

void Multipoles::AnalyticCartesianMultipoles(FLOAT3 *p, int n, FLOAT3 center, 
                                             double *cm) {

    for(int m=0;m<completemultipolelength;m++) cm[m] = 0;

    // We up-cast the positions in the AVX versions, so for consistency do that here
    double3 dcenter = double3(center);
    for(int q=0;q<n;q++) {
        double3 r =  p[q] - dcenter;

        double fi,fij,fijk;
        fi = 1.0;

        for(int i=0;i<=order;i++) {
            fij = fi;
            for(int j=0;j<=order-i;j++) {
                fijk = fij;
                for(int k=0;k<=order-i-j;k++) {
                    cm[ cmap(i,j,k) ] += fijk;
                    fijk *= r.z;
                }
                fij *= r.y;
            }
            fi *= r.x;
        }
    }
}

#ifdef UNROLLEDMULTIPOLES

void Multipoles::UnrolledCartesianMultipoles(FLOAT3 *p, int n, FLOAT3 center, 
                                             double *CM) {

    for(int m = 0; m < completemultipolelength; m++)
        CM[m] = 0;

    // We up-cast the positions in the AVX versions, so for consistency do that here
    double3 dcenter = double3(center);
    for(int q=0;q<n;q++) {
        CM_unrolled_ptr[order-1](p[q], dcenter, CM);
    }
}

#endif

#ifdef AVXMULTIPOLES
void Multipoles::ASMCartesianMultipoles(FLOAT3 *xyz, int n, FLOAT3 center, 
                                        double *CM) {
    int g = omp_get_thread_num();

    // TODO: off by one error?
    for(int k=0;k<completemultipolelength;k++) CM[k] = 0;

    for(int k=0;k<completemultipolelength;k++) 
        for(int j=0;j<4;j++) globalM[g][k].v[j] = 0;

    for(int j=0;j<4;j++) cx[g][0].v[j] = center.x;
    for(int j=0;j<4;j++) cy[g][0].v[j] = center.y;
    for(int j=0;j<4;j++) cz[g][0].v[j] = center.z;

    int end = n-(n%8);

    for(int i=end;i<n;i++) {    
        double cp[completemultipolelength];
        FLOAT3 p = xyz[i];
        AnalyticCartesianMultipoles( &p, 1, center,  &(cp[0]) );
        for(int k=0;k<completemultipolelength;k++) CM[k] += cp[k];
        // do these with analytic and add to CM
    }

    for(int j=0;j<4;j++) masses1[g][0].v[j]  = 1;
    for(int j=0;j<4;j++) masses2[g][0].v[j]  = 1;

    for(int k=0;k<end;k+=8) {

        for(int j=0;j<4;j++) ip1x[g][0].v[j] = xyz[k+j].x; 
        for(int j=0;j<4;j++) ip2x[g][0].v[j] = xyz[k+4+j].x;

        for(int j=0;j<4;j++) ip1y[g][0].v[j] = xyz[k+j].y;
        for(int j=0;j<4;j++) ip2y[g][0].v[j] = xyz[k+4+j].y;

        for(int j=0;j<4;j++) ip1z[g][0].v[j] = xyz[k+j].z;
        for(int j=0;j<4;j++) ip2z[g][0].v[j] = xyz[k+4+j].z;

        (CMptr[order-1])( &(ip1x[g][0]), &(ip2x[g][0]), 
                        &(ip1y[g][0]), &(ip2y[g][0]),
                        &(ip1z[g][0]), &(ip2z[g][0]),
                        &(cx[g][0]), &(cy[g][0]), &(cz[g][0]), 
                        &(globalM[g][0]),  &(masses1[g][0]), &(masses2[g][0]) );

    }

    for(int k=0;k<completemultipolelength;k++) 
       for(int j=0;j<4;j++) CM[k] += globalM[g][k].v[j];
}
#endif

void Multipoles::AVX512CartesianMultipoles(FLOAT3 *xyz, int n, FLOAT3 center, double *CM) {
#ifdef AVX512MULTIPOLES
    AVX512_DOUBLES cx512 = AVX512_SET_DOUBLE(center.x);
    AVX512_DOUBLES cy512 = AVX512_SET_DOUBLE(center.y);
    AVX512_DOUBLES cz512 = AVX512_SET_DOUBLE(center.z);

    AVX512_DOUBLES CM512[cml];
    for(int i = 0; i < cml; i++)
        CM512[i] = AVX512_SETZERO_DOUBLE();

    //AVX512_DOUBLES zk[order+1];
    //zk[0] = AVX512_SET_DOUBLE(1.0);

    int n_aligned = n - (n % AVX512_NVEC_DOUBLE);
    int nleft = n - n_aligned;

    for(int k=0; k <= n_aligned-AVX512_NVEC_DOUBLE; k += AVX512_NVEC_DOUBLE) {
        // Load 8 DOUBLE3s as List3s
        AVX512_DOUBLES px512, py512, pz512;
        for(int j = 0; j < AVX512_NVEC_DOUBLE; j++){
            px512[j] = xyz[k+j].x;
            py512[j] = xyz[k+j].y;
            pz512[j] = xyz[k+j].z;
        }

        // This function calls a manually unrolled version of the commented-out code below
        CM512ptr[order-1](px512, py512, pz512,
                        cx512, cy512, cz512,
                        CM512);

        /*AVX512_DOUBLES deltax = AVX512_SUBTRACT_DOUBLES(px512, cx512);
        AVX512_DOUBLES deltay = AVX512_SUBTRACT_DOUBLES(py512, cy512);
        AVX512_DOUBLES deltaz = AVX512_SUBTRACT_DOUBLES(pz512, cz512);

        // Evaluate these 8 particles
        AVX512_DOUBLES fi,fij; //,fijk;
        fi = AVX512_SET_DOUBLE(1.0);

        // Precompute z^k
        for(int i = 1; i <= order; i++)
            zk[i] = AVX512_MULTIPLY_DOUBLES(zk[i-1], deltaz);

        int i = 0;
        for(int a=0;a<=order;a++) {
            fij = fi;
            for(int b=0;b<=order-a;b++) {
                //fijk = fij;
                for(int c=0;c<=order-a-b;c++) {
                    //CM512[i] = AVX512_ADD_DOUBLES(CM512[i], fijk);
                    //fijk = AVX512_MULTIPLY_DOUBLES(fijk, deltaz);
                    CM512[i] = AVX512_FMA_ADD_DOUBLES(fij, zk[c], CM512[i]);
                    i++;
                }
                fij = AVX512_MULTIPLY_DOUBLES(fij, deltay);
            }
            fi = AVX512_MULTIPLY_DOUBLES(fi, deltax);
        }*/
    }

    // We could manually unroll this masked version too.  Might be faster for small cells, which could be important
    if(n_aligned < n){
        // Load nleft DOUBLE3s as List3s
        AVX512_DOUBLES px512 = AVX512_SETZERO_DOUBLE();
        AVX512_DOUBLES py512 = AVX512_SETZERO_DOUBLE();
        AVX512_DOUBLES pz512 = AVX512_SETZERO_DOUBLE();
        for(int j = 0; j < nleft; j++){
            px512[j] = xyz[n_aligned+j].x;
            py512[j] = xyz[n_aligned+j].y;
            pz512[j] = xyz[n_aligned+j].z;
        }

        AVX512_DOUBLES deltax = AVX512_SUBTRACT_DOUBLES(px512, cx512);
        AVX512_DOUBLES deltay = AVX512_SUBTRACT_DOUBLES(py512, cy512);
        AVX512_DOUBLES deltaz = AVX512_SUBTRACT_DOUBLES(pz512, cz512);

        // Evaluate these 8 particles
        AVX512_DOUBLES fi,fij,fijk;
        fi = AVX512_SET_DOUBLE(1.0);
        AVX512_MASK_DOUBLE left_mask = masks_per_misalignment_value_double[nleft];

        int i = 0;
        for(int a=0;a<=order;a++) {
            fij = fi;
            for(int b=0;b<=order-a;b++) {
                fijk = fij;
                for(int c=0;c<=order-a-b;c++) {
                    // Using the FMA version here isn't faster
                    CM512[i] = AVX512_MASK_ADD_DOUBLES(CM512[i], left_mask, CM512[i], fijk);
                    i++;
                    fijk = AVX512_MULTIPLY_DOUBLES(fijk, deltaz);
                }
                fij = AVX512_MULTIPLY_DOUBLES(fij, deltay);
            }
            fi = AVX512_MULTIPLY_DOUBLES(fi, deltax);
        }
    }

    int i = 0;
    for(int a=0;a<=order;a++)
        for(int b=0;b<=order-a;b++)
            for(int c=0;c<=order-a-b;c++)
                // TODO: I think this can just be i
                CM[cmap(a,b,c)] = AVX512_HORIZONTAL_SUM_DOUBLES(CM512[i++]);
#endif
}

#ifdef TEST
#include <cstdlib>
#include <cstring>
#include <chrono>
#include <gsl/gsl_rng.h>

void compare_multipoles(double *cm1, double* cm2, int64_t n, double rtol){
    int64_t nbad = 0;
    double max_frac_diff = 0;
    for(int64_t i = 0; i < n; i++){
        double m1 = cm1[i];
        double m2 = cm2[i];
        if (m1 == m2)
            continue;
        double frac_diff = std::abs(m1 - m2)/(std::abs(m1) + std::abs(m2));
        max_frac_diff = std::max(max_frac_diff, frac_diff);
        if(frac_diff > rtol || !std::isfinite(frac_diff)){
            nbad++;
            //std::cout << cm1[i] <<std::endl;
            //std::cout << cm2[i] <<std::endl;
        }
    }
    printf("\t>>> %lld (%.2f%%) mismatched multipoles\n", nbad, (FLOAT) nbad/n*100);
    printf("\t>>> Max frac error: %.2g \n", max_frac_diff);
}

void report(const char* prefix, int64_t npart, std::chrono::duration<double> elapsed){
    std::cout << prefix << " time: " << elapsed.count() << " sec" << std::endl;
    std::cout << "\t" << npart/1e6/elapsed.count() << " Mpart per second" << std::endl;
    std::cout.flush();
}

// Usage: ./EvaluateMultipolesTest [ppc]
int main(int argc, char **argv){
    Multipoles MP(8);

    int64_t ncell = 10*1875*1875;
    int64_t ppc = 44;
    if (argc > 1)
        ppc = atoi(argv[1]);
    float rtol=1e-6;
    int64_t npart = (int64_t)ncell*ppc;

    double *cartesian1, *cartesian2, *cartesian3, *cartesian4;
    FLOAT3 center(0.1,0.2,0.3);
    FLOAT3 *xyz;

    // ========================== //
    // Set up RNG
    int nthread = omp_get_max_threads();
    gsl_rng *rng[nthread];

    for(int i = 0; i < nthread; i++){
        rng[i] = gsl_rng_alloc(gsl_rng_mt19937);
        gsl_rng_set(rng[i], 999 + i);  // Seed the RNG
    }

    assert(posix_memalign((void **) &cartesian1, 4096, sizeof(double)*MP.cml*ncell) == 0);
    assert(posix_memalign((void **) &cartesian2, 4096, sizeof(double)*MP.cml*ncell) == 0);
    assert(posix_memalign((void **) &cartesian3, 4096, sizeof(double)*MP.cml*ncell) == 0);
    assert(posix_memalign((void **) &cartesian4, 4096, sizeof(double)*MP.cml*ncell) == 0);

    assert(posix_memalign((void **) &xyz, 4096, sizeof(FLOAT3)*npart) == 0);
    #pragma omp parallel for schedule(static)
    for(uint64_t i = 0; i < npart; i++){
        int t = omp_get_thread_num();
        xyz[i].x = gsl_rng_uniform(rng[t]);
        xyz[i].y = gsl_rng_uniform(rng[t]);
        xyz[i].z = gsl_rng_uniform(rng[t]);
    }

    /****************************************/

    auto begin = std::chrono::steady_clock::now();
    auto end = std::chrono::steady_clock::now();

#ifdef AVXMULTIPOLES
    // ASM Multipoles
    begin = std::chrono::steady_clock::now();
    #pragma omp parallel for schedule(static)
    for(int64_t k = 0; k < ncell; k++){
        FLOAT3 *thisxyz = xyz + k*ppc;
        double *thisct = cartesian1 + k*MP.cml;
        
        MP.ASMCartesianMultipoles(thisxyz, ppc, center, thisct);
    }
    end = std::chrono::steady_clock::now();
    report("ASM Multipoles", npart, end-begin);
#endif

#ifdef AVX512MULTIPOLES
    // AVX-512 Multipoles
    begin = std::chrono::steady_clock::now();
    #pragma omp parallel for schedule(static)
    for(int64_t k = 0; k < ncell; k++){
        FLOAT3 *thisxyz = xyz + k*ppc;
        double *thisct = cartesian3 + k*MP.cml;
        
        MP.AVX512CartesianMultipoles(thisxyz, ppc, center, thisct);
    }
    end = std::chrono::steady_clock::now();
    report("AVX-512 Multipoles", npart, end-begin);
    compare_multipoles(cartesian1, cartesian3, MP.cml*ncell, rtol);
#endif

    // Analytic Multipoles
    begin = std::chrono::steady_clock::now();
    #pragma omp parallel for schedule(static)
    for(int64_t k = 0; k < ncell; k++){
        FLOAT3 *thisxyz = xyz + k*ppc;
        double *thisct = cartesian2 + k*MP.cml;
        
        MP.AnalyticCartesianMultipoles(thisxyz, ppc, center, thisct);
    }
    end = std::chrono::steady_clock::now();
    report("Analytic Multipoles", npart, end-begin);
    //compare_multipoles(cartesian1, cartesian2, MP.cml*ncell, rtol);

#ifdef UNROLLEDMULTIPOLES
    // Unrolled Multipoles
    begin = std::chrono::steady_clock::now();
    #pragma omp parallel for schedule(static)
    for(int64_t k = 0; k < ncell; k++){
        FLOAT3 *thisxyz = xyz + k*ppc;
        double *thisct = cartesian4 + k*MP.cml;
        
        MP.UnrolledCartesianMultipoles(thisxyz, ppc, center, thisct);
    }
    end = std::chrono::steady_clock::now();
    report("Unrolled Multipoles", npart, end-begin);
    compare_multipoles(cartesian2, cartesian4, MP.cml*ncell, rtol);
#endif

    free(cartesian1); free(cartesian2); free(cartesian3);
    free(xyz);
    return 0;
}
#endif
