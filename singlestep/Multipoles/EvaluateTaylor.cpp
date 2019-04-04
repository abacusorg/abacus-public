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

/*
// For now we are using stack-allocated vectors
#ifdef AVX512MULTIPOLES
    for(int g=0;g<omp_get_max_threads();g++) {
        int rv;
        // 64-byte alignment or larger is fine for AVX-512
        // Why 512 bytes?  That's 16 AVX vectors.  maxorder vectors?  or just padding?
        rv = posix_memalign( (void **) &(cx512[g]), 256, 1024 ); assert(rv==0);
        rv = posix_memalign( (void **) &(cy512[g]), 256, 1024 ); assert(rv==0);
        rv = posix_memalign( (void **) &(cz512[g]), 256, 1024 ); assert(rv==0);

        rv = posix_memalign( (void **) &(ax512[g]), 256, 1024 ); assert(rv==0);
        rv = posix_memalign( (void **) &(ay512[g]), 256, 1024 ); assert(rv==0);
        rv = posix_memalign( (void **) &(az512[g]), 256, 1024 ); assert(rv==0);

        rv = posix_memalign( (void **) &(px512[g]), 256, 1024 ); assert(rv==0);
        rv = posix_memalign( (void **) &(py512[g]), 256, 1024 ); assert(rv==0);
        rv = posix_memalign( (void **) &(pz512[g]), 256, 1024 ); assert(rv==0);

        rv = posix_memalign( (void **) &(Qx512[g]), 256, 65536 ); assert(rv==0);
        rv = posix_memalign( (void **) &(Qy512[g]), 256, 65536 ); assert(rv==0);
        rv = posix_memalign( (void **) &(Qz512[g]), 256, 65536 ); assert(rv==0);
        assert(rv == 0);
    }
#endif
*/
}

// This is the dispatch function
void Taylor::EvaluateTaylor(double *CT, FLOAT3 expansioncenter, int np,
                                    FLOAT3 *ps, FLOAT3 *acc){
    #ifdef AVX512MULTIPOLES
    AVX512EvaluateTaylor(CT, expansioncenter, np, ps, acc);
    #elif defined(AVXMULTIPOLES)
    ASMEvaluateTaylor(CT, expansioncenter, np, ps, acc);
    #elif defined(UNROLLEDMULTIPOLES)
    UnrolledEvaluateTaylor(CT, expansioncenter, np, ps, acc);
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
    int cml_orderm1 = (order)*(order+1)*(order+2)/6;
    double3 Q[cml_orderm1];


    // We could unroll the below loop for more speed, but it already runs at ~60 Mpart/s/core
    int i,a,b,c;
    i = 0;
    FOR(a,0,order-1)
        FOR(b,0,order-1-a)
            FOR(c,0,order-1-a-b) {
                Q[i].x = (a+1)*CT[ cmap(a+1,b  ,c  ) ];
                Q[i].y = (b+1)*CT[ cmap(a  ,b+1,c  ) ];
                Q[i].z = (c+1)*CT[ cmap(a  ,b  ,c+1) ];
                i++;
            }

    // We up-cast the positions in the AVX versions, so for consistency do that here
    double3 dcenter = double3(center);

    // This function call contains the loop over particles; allows vectorization
    Tptr_unrolled[order](p, n, dcenter, Q, acc);
}

#endif

#ifdef AVXMULTIPOLES
void Taylor::ASMEvaluateTaylor( double *CT, FLOAT3 center, int n, FLOAT3 *xyz,
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

        (Tptr[order])( &(px[g][0]),&(py[g][0]),&(pz[g][0]),
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
    AVX512_DOUBLES _cx512 = AVX512_SET_DOUBLE(center.x);
    AVX512_DOUBLES _cy512 = AVX512_SET_DOUBLE(center.y);
    AVX512_DOUBLES _cz512 = AVX512_SET_DOUBLE(center.z);

    // I think this has to be big enough to enumerate all compresed taylors in the loop below
    int cml_orderm1 = (order)*(order+1)*(order+2)/6;
    AVX512_DOUBLES _Qx512[cml_orderm1];
    AVX512_DOUBLES _Qy512[cml_orderm1];
    AVX512_DOUBLES _Qz512[cml_orderm1];

    int i,a,b,c;
    i = 0;
    FOR(a,0,order-1)
        FOR(b,0,order-1-a)
            FOR(c,0,order-1-a-b) {
                _Qx512[i] = AVX512_SET_DOUBLE((a+1)*CT[ cmap(a+1,b  ,c  ) ]);
                _Qy512[i] = AVX512_SET_DOUBLE((b+1)*CT[ cmap(a  ,b+1,c  ) ]);
                _Qz512[i] = AVX512_SET_DOUBLE((c+1)*CT[ cmap(a  ,b  ,c+1) ]);
                i++;
            }

    int n_aligned = n - (n % AVX512_NVEC_DOUBLE);

    for(int k=0; k <= n_aligned - AVX512_NVEC_DOUBLE; k += AVX512_NVEC_DOUBLE) {
        // Load 8 DOUBLE3s as List3s
        AVX512_DOUBLES _px512, _py512, _pz512;
        for(int j = 0; j < AVX512_NVEC_DOUBLE; j++){
            _px512[j] = xyz[k+j].x;
            _py512[j] = xyz[k+j].y;
            _pz512[j] = xyz[k+j].z;
        }

        AVX512_DOUBLES _ax512 = AVX512_SET_DOUBLE(0.);
        AVX512_DOUBLES _ay512 = AVX512_SET_DOUBLE(0.);
        AVX512_DOUBLES _az512 = AVX512_SET_DOUBLE(0.);

        // It doesn't seem faster to use this manually unrolled version
        /*(Tptr512[order])( _px512, _py512, _pz512,
                          _cx512, _cy512, _cz512,
                          _Qx512, _Qy512, _Qz512,
                          _ax512, _ay512, _az512 );*/

        // Evaluate these 8 particles
        AVX512_DOUBLES fi,fij,fijk;

        AVX512_DOUBLES deltax = AVX512_SUBTRACT_DOUBLES(_px512, _cx512);
        AVX512_DOUBLES deltay = AVX512_SUBTRACT_DOUBLES(_py512, _cy512);
        AVX512_DOUBLES deltaz = AVX512_SUBTRACT_DOUBLES(_pz512, _cz512);

        fi = AVX512_SET_DOUBLE(1.0);
        i = 0;  // nested counter
        FOR(a,0,order-1) {
            AVX512_DOUBLES fij = fi;
            FOR(b,0,order-1-a) {
                AVX512_DOUBLES fijk = fij;
                //#pragma unroll (2)  // not faster
                FOR(c,0,order-1-a-b) {
                    _ax512 = AVX512_FMA_ADD_DOUBLES(_Qx512[i], fijk, _ax512);
                    _ay512 = AVX512_FMA_ADD_DOUBLES(_Qy512[i], fijk, _ay512);
                    _az512 = AVX512_FMA_ADD_DOUBLES(_Qz512[i], fijk, _az512);
                    i++;
                    fijk = AVX512_MULTIPLY_DOUBLES(fijk, deltaz);
                }
                fij = AVX512_MULTIPLY_DOUBLES(fij, deltay);
            }
            fi = AVX512_MULTIPLY_DOUBLES(fi, deltax);
        }

        // Unpack the List3s into FLOAT3s
        for(int j = 0; j < AVX512_NVEC_DOUBLE; j++){
            acc[k+j].x -= _ax512[j];
            acc[k+j].y -= _ay512[j];
            acc[k+j].z -= _az512[j];
        }
    }

    if(n_aligned < n){
        int nleft = n - n_aligned;
        // Load nleft DOUBLE3s as List3s
        AVX512_DOUBLES _px512 = AVX512_SET_DOUBLE(0.);
        AVX512_DOUBLES _py512 = AVX512_SET_DOUBLE(0.);
        AVX512_DOUBLES _pz512 = AVX512_SET_DOUBLE(0.);
        for(int j = 0; j < nleft; j++){
            _px512[j] = xyz[n_aligned+j].x;
            _py512[j] = xyz[n_aligned+j].y;
            _pz512[j] = xyz[n_aligned+j].z;
        }

        // Evaluate these nleft particles
        AVX512_DOUBLES fi,fij,fijk;
        AVX512_DOUBLES _ax512 = AVX512_SET_DOUBLE(0.);
        AVX512_DOUBLES _ay512 = AVX512_SET_DOUBLE(0.);
        AVX512_DOUBLES _az512 = AVX512_SET_DOUBLE(0.);

        AVX512_DOUBLES deltax = AVX512_SUBTRACT_DOUBLES(_px512, _cx512);
        AVX512_DOUBLES deltay = AVX512_SUBTRACT_DOUBLES(_py512, _cy512);
        AVX512_DOUBLES deltaz = AVX512_SUBTRACT_DOUBLES(_pz512, _cz512);

        fi = AVX512_SET_DOUBLE(1.0);
        i = 0;  // nested counter
        FOR(a,0,order-1) {
            AVX512_DOUBLES fij = fi;
            FOR(b,0,order-1-a) {
                AVX512_DOUBLES fijk = fij;
                FOR(c,0,order-1-a-b) {
                    _ax512 = AVX512_FMA_ADD_DOUBLES(_Qx512[i], fijk, _ax512);
                    _ay512 = AVX512_FMA_ADD_DOUBLES(_Qy512[i], fijk, _ay512);
                    _az512 = AVX512_FMA_ADD_DOUBLES(_Qz512[i], fijk, _az512);
                    i++;
                    fijk = AVX512_MULTIPLY_DOUBLES(fijk, deltaz);
                }
                fij = AVX512_MULTIPLY_DOUBLES(fij, deltay);
            }
            fi = AVX512_MULTIPLY_DOUBLES(fi, deltax);
        }

        // Unpack the List3s into FLOAT3s
        for(int j = 0; j < nleft; j++){
            acc[n_aligned+j].x -= _ax512[j];
            acc[n_aligned+j].y -= _ay512[j];
            acc[n_aligned+j].z -= _az512[j];
        }
    }
#endif
}

#ifdef TEST
#include <cstdlib>
#include <cstring>
#include <chrono>
#include <gsl/gsl_rng.h>

void compare_acc(FLOAT3 *acc1, FLOAT3* acc2, int nacc, double rtol){
    int nbad = 0;
    double max_frac_diff = 0;
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
}

void report(const char* prefix, int64_t npart, std::chrono::duration<double> elapsed){
    std::cout << prefix << " time: " << elapsed.count() << " sec" << std::endl;
    std::cout << "\t" << npart/1e6/elapsed.count() << " Mpart per second" << std::endl;
}

int main(int argc, char **argv){
    Taylor TY(8);

    int ncell = 10*1875*1875;
    int ppc = 44;
    if (argc > 1)
        ppc = atoi(argv[1]);
    float rtol=1e-6;
    uint64_t npart = ncell*ppc;

    double *cartesian;  // re-use cartesians for all cells (?)
    FLOAT3 center(0.1,0.2,0.3);
    FLOAT3 *xyz;
    FLOAT3 *acc1, *acc2, *acc3, *acc4;

    // ========================== //
    // Set up RNG
    int nthread = omp_get_max_threads();
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
        xyz[i].x = gsl_rng_uniform(rng[t]);
        xyz[i].y = gsl_rng_uniform(rng[t]);
        xyz[i].z = gsl_rng_uniform(rng[t]);
    }

    assert(posix_memalign((void **) &acc1, 4096, sizeof(FLOAT3)*npart) == 0);
    assert(posix_memalign((void **) &acc2, 4096, sizeof(FLOAT3)*npart) == 0);
    assert(posix_memalign((void **) &acc3, 4096, sizeof(FLOAT3)*npart) == 0);
    assert(posix_memalign((void **) &acc4, 4096, sizeof(FLOAT3)*npart) == 0);

    auto begin = std::chrono::steady_clock::now();
    auto end = std::chrono::steady_clock::now();

#ifdef AVXMULTIPOLES
    // ASM Taylors
    begin = std::chrono::steady_clock::now();
    #pragma omp parallel for schedule(static)
    for(int k = 0; k < ncell; k++){
        FLOAT3 *thisxyz = xyz + k*ppc;
        FLOAT3 *thisacc = acc1 + k*ppc;
        
        memset(thisacc, 0, sizeof(FLOAT3)*ppc);
        TY.ASMEvaluateTaylor(cartesian, center, ppc, thisxyz, thisacc);
    }
    end = std::chrono::steady_clock::now();
    report("ASM Taylors", npart, end-begin);
#endif

#ifdef AVX512MULTIPOLES
    // AVX-512 Taylors
    begin = std::chrono::steady_clock::now();
    #pragma omp parallel for schedule(static)
    for(int k = 0; k < ncell; k++){
        FLOAT3 *thisxyz = xyz + k*ppc;
        FLOAT3 *thisacc = acc3 + k*ppc;
        
        memset(thisacc, 0, sizeof(FLOAT3)*ppc);
        TY.AVX512EvaluateTaylor(cartesian, center, ppc, thisxyz, thisacc);
    }
    end = std::chrono::steady_clock::now();
    report("AVX-512 Taylors", npart, end-begin);
    compare_acc(acc1, acc3, npart, rtol);
#endif


    // Analytic Taylors
    begin = std::chrono::steady_clock::now();
    #pragma omp parallel for schedule(static)
    for(int k = 0; k < ncell; k++){
        FLOAT3 *thisxyz = xyz + k*ppc;
        FLOAT3 *thisacc = acc2 + k*ppc;
        
        memset(thisacc, 0, sizeof(FLOAT3)*ppc);
        TY.AnalyticEvaluateTaylor(cartesian, center, ppc, thisxyz, thisacc);
    }
    end = std::chrono::steady_clock::now();
    report("Analytic Taylors", npart, end-begin);
    //compare_acc(acc1, acc2, npart, rtol);

#ifdef UNROLLEDMULTIPOLES
    // Analytic Taylors
    begin = std::chrono::steady_clock::now();
    #pragma omp parallel for schedule(static)
    for(int k = 0; k < ncell; k++){
        FLOAT3 *thisxyz = xyz + k*ppc;
        FLOAT3 *thisacc = acc4 + k*ppc;
        
        memset(thisacc, 0, sizeof(FLOAT3)*ppc);
        TY.UnrolledEvaluateTaylor(cartesian, center, ppc, thisxyz, thisacc);
    }
    end = std::chrono::steady_clock::now();
    report("Unrolled Taylors", npart, end-begin);
    compare_acc(acc4, acc2, npart, rtol);
#endif

    return 0;
}
#endif
