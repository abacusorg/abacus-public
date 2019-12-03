#ifdef TEST

// Compile test driver with:
// $ make EvaluateTaylorTest

// Can use these to override ./configure for testing
//#define AVXMULTIPOLES
//#define AVX512MULTIPOLES
#define UNROLLEDMULTIPOLES

#include "threevector.hh"
#define FLOAT float
#define FLOAT3 float3

#include "basemultipoles.cpp"

#endif

#include "EvaluateTaylor.h"


Taylor::~Taylor(void) {
#ifdef AVXMULTIPOLES
    for(int g=0;g<omp_get_max_threads();g++) {
        free(cx[g]);
        free(cy[g]);
        free(cz[g]);
        
        free(ax[g]);
        free(ay[g]);
        free(az[g]);
        
        free(px[g]);
        free(py[g]);
        free(pz[g]);

        free(Qx[g]);
        free(Qy[g]);
        free(Qz[g]);
    }
#endif
}

Taylor::Taylor(int order) : basemultipoles(order) {
    assert(omp_get_max_threads() < MAXTHREADS);
    
#ifdef AVXMULTIPOLES
    for(int g=0;g<omp_get_max_threads();g++) {
        int rv;
        rv = posix_memalign( (void **) &(cx[g]), 256, 512 ); assert(rv==0);
        rv = posix_memalign( (void **) &(cy[g]), 256, 512 ); assert(rv==0);
        rv = posix_memalign( (void **) &(cz[g]), 256, 512 ); assert(rv==0);

        rv = posix_memalign( (void **) &(ax[g]), 256, 512 ); assert(rv==0);
        rv = posix_memalign( (void **) &(ay[g]), 256, 512 ); assert(rv==0);
        rv = posix_memalign( (void **) &(az[g]), 256, 512 ); assert(rv==0);

        rv = posix_memalign( (void **) &(px[g]), 256, 512 ); assert(rv==0);
        rv = posix_memalign( (void **) &(py[g]), 256, 512 ); assert(rv==0);
        rv = posix_memalign( (void **) &(pz[g]), 256, 512 ); assert(rv==0);

        rv = posix_memalign( (void **) &(Qx[g]), 256, 32768 ); assert(rv==0);
        rv = posix_memalign( (void **) &(Qy[g]), 256, 32768 ); assert(rv==0);
        rv = posix_memalign( (void **) &(Qz[g]), 256, 32768 ); assert(rv==0);
        assert(rv == 0);
    }
#endif

}

// This is the dispatch function
void Taylor::EvaluateTaylor(double *CT, FLOAT3 expansioncenter, int np,
                                    FLOAT3 *ps, FLOAT3 *acc){
    #ifdef AVX512MULTIPOLES
    AVX512EvaluateTaylor(CT, expansioncenter, np, ps, acc);
    #elif defined(AVXMULTIPOLES)
    AVXEvaluateTaylor(CT, expansioncenter, np, ps, acc);
    #elif defined(UNROLLEDMULTIPOLES)
    UnrolledEvaluateTaylor(CT, expansioncenter, np, ps, acc);
    #elif defined(VSXMULTIPOLES)
    VSXEvaluateTaylor(CT, expansioncenter, np, ps, acc);
    #else
    AnalyticEvaluateTaylor(CT, expansioncenter, np, ps, acc);
    #endif
}

void Taylor::AnalyticEvaluateTaylor(double *CT, FLOAT3 expansioncenter, int np,
                                    FLOAT3 *ps, FLOAT3 *acc) {
    assert(np>=0);

    for(int p=0;p<np;p++) {

        int i,j,k;
        double fi,fij,fijk;
        double3 delta = ps[p] - expansioncenter;
        double3 a = double3(0);

        fi = 1.0;
        FOR(i,0,order-1) {
            fij = fi;
            FOR(j,0,order-1-i) {
                fijk = fij;
                FOR(k,0,order-1-i-j) {
                    a.x -= (i+1) * CT[ cmap(i+1,j  ,k  ) ] * fijk;
                    a.y -= (j+1) * CT[ cmap(i  ,j+1,k  ) ] * fijk;
                    a.z -= (k+1) * CT[ cmap(i  ,j  ,k+1) ] * fijk;
                    fijk *= delta.z;
                }
                fij *= delta.y;
            }
            fi *= delta.x;
        }
        acc[p] += a;
    }
}

#ifdef UNROLLEDMULTIPOLES

void Taylor::UnrolledEvaluateTaylor(double *CT, FLOAT3 center, int n, FLOAT3 *p, FLOAT3 *acc) {
    // We up-cast the positions in the AVX versions, so for consistency do that here
    double3 dcenter = double3(center);

    // This function call contains the loop over particles; allows vectorization
    DispatchTaylorUnrolledKernel(order, p, n, dcenter, CT, acc);
}

#endif

#ifdef VSXMULTIPOLES

void Taylor::VSXEvaluateTaylor(double *CT, FLOAT3 center, int n, FLOAT3 *p, FLOAT3 *acc) {

    // We up-cast the positions in the AVX versions, so for consistency do that here
    double3 dcenter = double3(center);

    // This function call contains the loop over particles; allows vectorization
    DispatchTaylorVSXKernel(order, p, n, dcenter, CT, acc);
}

#endif

#ifdef AVXMULTIPOLES
void Taylor::AVXEvaluateTaylor( double *CT, FLOAT3 center, int n, FLOAT3 *xyz,
                                FLOAT3 *acc) {
    int g = omp_get_thread_num();

    for(int j=0;j<4;j++) {  
        cx[g][0].v[j] = center.x;
        cy[g][0].v[j] = center.y;
        cz[g][0].v[j] = center.z;
    }

    int i,a,b,c;
    i = 0;
    FOR(a,0,order-1)
        FOR(b,0,order-1-a)
            FOR(c,0,order-1-a-b) {
                for(int j=0;j<4;j++) {
                    Qx[g][i].v[j] = (a+1)*CT[ cmap(a+1,b  ,c  ) ];
                    Qy[g][i].v[j] = (b+1)*CT[ cmap(a  ,b+1,c  ) ];
                    Qz[g][i].v[j] = (c+1)*CT[ cmap(a  ,b  ,c+1) ];
                }
                i++;
            }

    int l = 0;
    int m = n-1;

    for(int k=l;k<=m;k+=4) {
        int end=k+3; if(k+3>m) end = m;
        
        for(int j=k;j<=end;j++) {
            px[g][0].v[j-k] = xyz[j].x;
            py[g][0].v[j-k] = xyz[j].y;
            pz[g][0].v[j-k] = xyz[j].z;
        }

        for(int j=0;j<4;j++) ax[g][0].v[j] = 1;

        DispatchTaylorAVXKernel(order, &(px[g][0]),&(py[g][0]),&(pz[g][0]),
                       &(cx[g][0]),&(cy[g][0]),&(cz[g][0]),
                       &(Qx[g][0]),&(Qy[g][0]),&(Qz[g][0]),
            &(ax[g][0]),&(ay[g][0]),&(az[g][0]) );

        for(int j=k;j<=end;j++) {
            acc[j].x -= ax[g][0].v[j-k];
            acc[j].y -= ay[g][0].v[j-k];
            acc[j].z -= az[g][0].v[j-k];
        }
    }
}
#endif

void Taylor::AVX512EvaluateTaylor( double *CT, FLOAT3 center, int n, FLOAT3 *xyz,
                                FLOAT3 *acc) {
#ifdef AVX512MULTIPOLES
    DispatchTaylor512Kernel(order, CT, center, n, xyz, acc);
#endif
}

#ifdef TEST
#include <cstdlib>
#include <cstring>
#include <chrono>
#include <gsl/gsl_rng.h>

int have_any_results = 0;
void compare_acc(FLOAT3 *acc1, FLOAT3* acc2, int nacc, double rtol){
    if(!have_any_results){
        have_any_results = 1;
        return;
    }

    int nbad = 0;
    double max_frac_diff = 0;
    #pragma omp parallel for schedule(static) reduction(+:nbad) reduction(max:max_frac_diff)
    for(int i = 0; i < nacc; i++){
        double3 a1(acc1[i]);
        double3 a2(acc2[i]);
        if (a1 == a2)
            continue;
        double frac_diff = (a1 - a2).norm()/(a1 + a2).norm();
        max_frac_diff = std::max(max_frac_diff, frac_diff);
        if(frac_diff > rtol || !std::isfinite(frac_diff)){
            nbad++;
            //std::cout << acc1[i] <<std::endl;
            //std::cout << acc2[i] <<std::endl;
        }
    }
    printf("\t>>> %d (%.2f%%) mismatched accels\n", nbad, (FLOAT) nbad/nacc*100);
    printf("\t>>> Max frac error: %.2g \n", max_frac_diff);
    fflush(stdout);
}

void report(const char* prefix, int64_t npart, std::chrono::duration<double> elapsed, int nthread){
    //double nflop = 842*npart + 360;  // Q version
    double nflop = 1199*npart;  // no Q version

    auto t = elapsed.count();

    std::cout << prefix << " time: " << t << " sec" << std::endl;
    printf("\t %.3f Mpart per second (%.3g DP-GFLOPS per thread)\n", npart/1e6/t, nflop/1e9/t/nthread);
    fflush(stdout);
}

int main(int argc, char **argv){
    Taylor TY(8);

    int64_t cpd = 65;
    int64_t ncell = 1*cpd*cpd;
    int64_t ppc = 52;
    if (argc > 1)
        ppc = atoi(argv[1]);
    float rtol=1e-6;
    uint64_t npart = ncell*ppc;

    double *cartesian;  // re-use cartesians for all cells (?)
    FLOAT3 center(1e-3,2e-3,3e-3);
    FLOAT3 *xyz;
    FLOAT3 *current_acc, *last_acc;

    // ========================== //
    // Set up RNG
    int nthread = omp_get_max_threads();
    printf("Running with %d threads, ppc %zd\n", nthread, ppc);
    gsl_rng *rng[nthread];

    for(int i = 0; i < nthread; i++){
        rng[i] = gsl_rng_alloc(gsl_rng_mt19937);
        gsl_rng_set(rng[i], 999 + i);  // Seed the RNG
    }


    assert(posix_memalign((void **) &cartesian, 4096, sizeof(double)*TY.cml) == 0);
    for (int i = 0; i < TY.cml; i++)
        cartesian[i] = (double)rand()/RAND_MAX;

    assert(posix_memalign((void **) &xyz, 4096, sizeof(FLOAT3)*npart) == 0);
    #pragma omp parallel for schedule(static)
    for(uint64_t i = 0; i < npart; i++){
        int t = omp_get_thread_num();
        xyz[i].x = gsl_rng_uniform(rng[t])/cpd;
        xyz[i].y = gsl_rng_uniform(rng[t])/cpd;
        xyz[i].z = gsl_rng_uniform(rng[t])/cpd;
    }

    assert(posix_memalign((void **) &last_acc, 4096, sizeof(FLOAT3)*npart) == 0);
    assert(posix_memalign((void **) &current_acc, 4096, sizeof(FLOAT3)*npart) == 0);

    auto begin = std::chrono::steady_clock::now();
    auto end = std::chrono::steady_clock::now();

#ifdef AVXMULTIPOLES
    // AVX Taylors

    // zero the outputs
    std::swap(current_acc, last_acc);
    #pragma omp parallel for schedule(static)
    for(int64_t k = 0; k < ncell; k++){
        FLOAT3 *thisacc = current_acc + k*ppc;
        memset(thisacc, 0, ppc*sizeof(FLOAT3));
    }

    begin = std::chrono::steady_clock::now();
    #pragma omp parallel for schedule(static)
    for(int64_t k = 0; k < ncell; k++){
        FLOAT3 *thisxyz = xyz + k*ppc;
        FLOAT3 *thisacc = current_acc + k*ppc;
        
        TY.AVXEvaluateTaylor(cartesian, center, ppc, thisxyz, thisacc);
    }
    end = std::chrono::steady_clock::now();
    report("AVX Taylors", npart, end-begin, nthread);
    compare_acc(current_acc, last_acc, npart, rtol);
#endif

#ifdef AVX512MULTIPOLES
    // AVX-512 Taylors

    // zero the outputs
    std::swap(current_acc, last_acc);
    #pragma omp parallel for schedule(static)
    for(int64_t k = 0; k < ncell; k++){
        FLOAT3 *thisacc = current_acc + k*ppc;
        memset(thisacc, 0, ppc*sizeof(FLOAT3));
    }

    begin = std::chrono::steady_clock::now();
    #pragma omp parallel for schedule(static)
    for(int64_t k = 0; k < ncell; k++){
        FLOAT3 *thisxyz = xyz + k*ppc;
        FLOAT3 *thisacc = current_acc + k*ppc;
        
        TY.AVX512EvaluateTaylor(cartesian, center, ppc, thisxyz, thisacc);
    }
    end = std::chrono::steady_clock::now();
    report("AVX-512 Taylors", npart, end-begin, nthread);
    compare_acc(current_acc, last_acc, npart, rtol);
#endif

#ifdef VSXMULTIPOLES
    // Analytic Taylors

    // zero the outputs
    std::swap(current_acc, last_acc);
    #pragma omp parallel for schedule(static)
    for(int64_t k = 0; k < ncell; k++){
        FLOAT3 *thisacc = current_acc + k*ppc;
        memset(thisacc, 0, ppc*sizeof(FLOAT3));
    }

    begin = std::chrono::steady_clock::now();
    #pragma omp parallel for schedule(static)
    for(int64_t k = 0; k < ncell; k++){
        FLOAT3 *thisxyz = xyz + k*ppc;
        FLOAT3 *thisacc = current_acc + k*ppc;

        TY.VSXEvaluateTaylor(cartesian, center, ppc, thisxyz, thisacc);
    }
    end = std::chrono::steady_clock::now();
    report("VSX Taylors", npart, end-begin, nthread);
    compare_acc(current_acc, last_acc, npart, rtol);
#endif


#ifdef UNROLLEDMULTIPOLES
    // Unrolled Taylors

    // zero the outputs
    std::swap(current_acc, last_acc);
    #pragma omp parallel for schedule(static)
    for(int64_t k = 0; k < ncell; k++){
        FLOAT3 *thisacc = current_acc + k*ppc;
        memset(thisacc, 0, ppc*sizeof(FLOAT3));
    }

    begin = std::chrono::steady_clock::now();
    #pragma omp parallel for schedule(static)
    for(int64_t k = 0; k < ncell; k++){
        FLOAT3 *thisxyz = xyz + k*ppc;
        FLOAT3 *thisacc = current_acc + k*ppc;
        
        TY.UnrolledEvaluateTaylor(cartesian, center, ppc, thisxyz, thisacc);
    }
    end = std::chrono::steady_clock::now();
    report("Unrolled Taylors", npart, end-begin, nthread);
    compare_acc(current_acc, last_acc, npart, rtol);
#endif

    // Analytic Taylors

    // zero the outputs
    std::swap(current_acc, last_acc);
    #pragma omp parallel for schedule(static)
    for(int64_t k = 0; k < ncell; k++){
        FLOAT3 *thisacc = current_acc + k*ppc;
        memset(thisacc, 0, ppc*sizeof(FLOAT3));
    }
    
    begin = std::chrono::steady_clock::now();
    #pragma omp parallel for schedule(static)
    for(int64_t k = 0; k < ncell; k++){
        FLOAT3 *thisxyz = xyz + k*ppc;
        FLOAT3 *thisacc = current_acc + k*ppc;
        
        TY.AnalyticEvaluateTaylor(cartesian, center, ppc, thisxyz, thisacc);
    }
    end = std::chrono::steady_clock::now();
    report("Analytic Taylors", npart, end-begin, nthread);
    compare_acc(current_acc, last_acc, npart, rtol);

    return 0;
}
#endif
