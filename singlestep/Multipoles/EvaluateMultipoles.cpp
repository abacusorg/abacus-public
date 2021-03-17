#ifdef TEST

#include "config.h"

// Compile test driver with:
// $ make EvaluateMultipolesTest

// Can use these to override ./configure for testing

#ifdef HAVE_AVX
//#define AVXMULTIPOLES
#endif

#ifdef HAVE_AVX512_F
#define AVX512MULTIPOLES
#endif

#define UNROLLEDMULTIPOLES

#ifdef HAVE_VSX
//#define VSXMULTIPOLES
#endif

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
    AVXCartesianMultipoles(p, n, center, cm);
    #elif defined(VSXMULTIPOLES)
    VSXCartesianMultipoles(p, n, center, cm);
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

#ifdef VSXMULTIPOLES

void Multipoles::VSXCartesianMultipoles(FLOAT3 *p, int n, FLOAT3 center, 
                                             double *CM) {

    for(int m = 0; m < completemultipolelength; m++)
        CM[m] = 0;

    // We up-cast the positions in the AVX versions, so for consistency do that here
    double3 dcenter = double3(center);
    DispatchMultipoleVSXKernel(order, p, n, dcenter, CM);
}

#endif

#ifdef UNROLLEDMULTIPOLES

void Multipoles::UnrolledCartesianMultipoles(FLOAT3 *p, int n, FLOAT3 center, 
                                             double *CM) {

    for(int m = 0; m < completemultipolelength; m++)
        CM[m] = 0;

    // We up-cast the positions in the AVX versions, so for consistency do that here
    double3 dcenter = double3(center);
    DispatchMultipoleUnrolledKernel(order, p, n, dcenter, CM);
}

#endif

#ifdef AVXMULTIPOLES
void Multipoles::AVXCartesianMultipoles(FLOAT3 *xyz, int n, FLOAT3 center, 
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

        DispatchMultipoleAVXKernel(order, &(ip1x[g][0]), &(ip2x[g][0]), 
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
    DispatchMultipole512Kernel(order, xyz, n, center, CM);
#endif
}

/****************************************************************************************/

#ifdef TEST
#include <cstdlib>
#include <cstring>
#include <chrono>
#include <gsl/gsl_rng.h>

int have_any_results = 0;
void compare_multipoles(double *cm1, double* cm2, int64_t n, double rtol){
    if(!have_any_results){
        have_any_results = 1;
        return;
    }

    int64_t nbad = 0;
    double max_frac_diff = 0;
    #pragma omp parallel for schedule(static) reduction(+:nbad) reduction(max:max_frac_diff)
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
    printf("\t>>> %zd (%.2f%%) mismatched multipoles\n", nbad, (FLOAT) nbad/n*100);
    printf("\t>>> Max frac error: %.2g \n", max_frac_diff);
    fflush(stdout);
}

void report(const char* prefix, int64_t npart, std::chrono::duration<double> elapsed, int cml, int nthread){
    //double nflop = 332; //  order 2: 22; order 8: 332 (per particle)
    double nflop = 385; //  FMA version: order 2: ?; order 8: 385 (per particle)
    auto t = elapsed.count();

    std::cout << prefix << " time: " << t << " sec" << std::endl;
    printf("\t%.3f Mpart per second (%.2g DP-GFLOPS per thread)\n", npart/1e6/t, nflop*npart/1e9/t/nthread);
    std::cout.flush();
}

// Usage: ./EvaluateMultipolesTest [ppc]
int main(int argc, char **argv){
    Multipoles MP(8);

    int64_t cpd = 1875;
    int64_t ncell = 1*cpd*cpd;
    int64_t ppc = 52;   
    if (argc > 1)
        ppc = atoi(argv[1]);
    float rtol=1e-6;
    int64_t npart = (int64_t)ncell*ppc;

    double *current_cartesian = NULL, *last_cartesian = NULL;
    FLOAT3 center(1e-3,2e-3,3e-3);
    FLOAT3 *xyz;

    // ========================== //
    // Set up RNG
    int nthread = omp_get_max_threads();
    printf("Running with %d threads, ppc %zd\n", nthread, ppc);
    gsl_rng *rng[nthread];

    for(int i = 0; i < nthread; i++){
        rng[i] = gsl_rng_alloc(gsl_rng_mt19937);
        gsl_rng_set(rng[i], 999 + i);  // Seed the RNG
    }

    assert(posix_memalign((void **) &xyz, PAGE_SIZE, sizeof(FLOAT3)*npart) == 0);
    #pragma omp parallel for schedule(static)
    for(int64_t i = 0; i < npart; i++){
        int t = omp_get_thread_num();
        xyz[i].x = gsl_rng_uniform(rng[t])/cpd;
        xyz[i].y = gsl_rng_uniform(rng[t])/cpd;
        xyz[i].z = gsl_rng_uniform(rng[t])/cpd;
    }

    assert(posix_memalign((void **) &current_cartesian, PAGE_SIZE, sizeof(double)*MP.cml*ncell) == 0);
    assert(posix_memalign((void **) &last_cartesian, PAGE_SIZE, sizeof(double)*MP.cml*ncell) == 0);

    /****************************************/

    auto begin = std::chrono::steady_clock::now();
    auto end = std::chrono::steady_clock::now();

#ifdef AVXMULTIPOLES
    // AVX Multipoles

    // zero the outputs
    std::swap(current_cartesian, last_cartesian);
    #pragma omp parallel for schedule(static)
    for(int64_t k = 0; k < ncell; k++){
        double *thisct = current_cartesian + k*MP.cml;
        memset(thisct, 0, MP.cml*sizeof(double));
    }

    begin = std::chrono::steady_clock::now();
    #pragma omp parallel for schedule(static)
    for(int64_t k = 0; k < ncell; k++){
        FLOAT3 *thisxyz = xyz + k*ppc;
        double *thisct = current_cartesian + k*MP.cml;
        
        MP.AVXCartesianMultipoles(thisxyz, ppc, center, thisct);
    }
    end = std::chrono::steady_clock::now();
    report("AVX Multipoles", npart, end-begin, MP.cml, nthread);
    compare_multipoles(current_cartesian, last_cartesian, MP.cml*ncell, rtol);
#endif

#ifdef AVX512MULTIPOLES
    // AVX-512 Multipoles

    // zero the outputs
    std::swap(current_cartesian, last_cartesian);
    #pragma omp parallel for schedule(static)
    for(int64_t k = 0; k < ncell; k++){
        double *thisct = current_cartesian + k*MP.cml;
        memset(thisct, 0, MP.cml*sizeof(double));
    }

    begin = std::chrono::steady_clock::now();
    #pragma omp parallel for schedule(static)
    for(int64_t k = 0; k < ncell; k++){
        FLOAT3 *thisxyz = xyz + k*ppc;
        double *thisct = current_cartesian + k*MP.cml;
        
        MP.AVX512CartesianMultipoles(thisxyz, ppc, center, thisct);
    }
    end = std::chrono::steady_clock::now();
    report("AVX-512 Multipoles", npart, end-begin, MP.cml, nthread);
    compare_multipoles(current_cartesian, last_cartesian, MP.cml*ncell, rtol);
#endif

#ifdef VSXMULTIPOLES
    // VSX Multipoles
    
    // zero the outputs
    std::swap(current_cartesian, last_cartesian);
    #pragma omp parallel for schedule(static)
    for(int64_t k = 0; k < ncell; k++){
        double *thisct = current_cartesian + k*MP.cml;
        memset(thisct, 0, MP.cml*sizeof(double));
    }
    
    begin = std::chrono::steady_clock::now();
    #pragma omp parallel for schedule(static)
    for(int64_t k = 0; k < ncell; k++){
        FLOAT3 *thisxyz = xyz + k*ppc;
        double *thisct = current_cartesian + k*MP.cml;
        
        MP.VSXCartesianMultipoles(thisxyz, ppc, center, thisct);
    }
    end = std::chrono::steady_clock::now();
    report("VSX Multipoles", npart, end-begin, MP.cml, nthread);
    //return 0;
    compare_multipoles(current_cartesian, last_cartesian, MP.cml*ncell, rtol);
#endif

#ifdef UNROLLEDMULTIPOLES
    // Unrolled Multipoles
    
    // zero the outputs
    std::swap(current_cartesian, last_cartesian);
    #pragma omp parallel for schedule(static)
    for(int64_t k = 0; k < ncell; k++){
        double *thisct = current_cartesian + k*MP.cml;
        memset(thisct, 0, MP.cml*sizeof(double));
    }

    begin = std::chrono::steady_clock::now();
    #pragma omp parallel for schedule(static)
    for(int64_t k = 0; k < ncell; k++){
        FLOAT3 *thisxyz = xyz + k*ppc;
        double *thisct = current_cartesian + k*MP.cml;
        
        MP.UnrolledCartesianMultipoles(thisxyz, ppc, center, thisct);
    }
    end = std::chrono::steady_clock::now();
    report("Unrolled Multipoles", npart, end-begin, MP.cml, nthread);
    compare_multipoles(current_cartesian, last_cartesian, MP.cml*ncell, rtol);
#endif

    // Analytic Multipoles
    
    // zero the outputs
    std::swap(current_cartesian, last_cartesian);
    #pragma omp parallel for schedule(static)
    for(int64_t k = 0; k < ncell; k++){
        double *thisct = current_cartesian + k*MP.cml;
        memset(thisct, 0, MP.cml*sizeof(double));
    }

    begin = std::chrono::steady_clock::now();
    #pragma omp parallel for schedule(static)
    for(int64_t k = 0; k < ncell; k++){
        FLOAT3 *thisxyz = xyz + k*ppc;
        double *thisct = current_cartesian + k*MP.cml;
        
        MP.AnalyticCartesianMultipoles(thisxyz, ppc, center, thisct);
    }
    end = std::chrono::steady_clock::now();
    report("Analytic Multipoles", npart, end-begin, MP.cml, nthread);
    compare_multipoles(current_cartesian, last_cartesian, MP.cml*ncell, rtol);

    free(current_cartesian); free(last_cartesian);
    free(xyz);
    return 0;
}
#endif
