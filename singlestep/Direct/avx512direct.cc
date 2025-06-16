// Copyright 2012-2025 The Abacus Developers
// SPDX-License-Identifier: GPL-3.0-or-later

/* AVX-512 implementations of our directs kernels
   Some of the AVX-512 ideas are based on Manodeep Sinha's AVX-512 pair counting implementation for Corrfunc.

    Compile test driver with:
    $ make avx512direct
*/


#ifdef TEST
#define FLOAT float
#define FLOAT3 float3
#define RSQRT 1.0f/sqrtf
#endif

#include "threevector.hh"
#include "StructureOfLists.cc"
#include "avx512_calls.c"

class AVX512Direct {
public:
    AVX512Direct(){};
    ~AVX512Direct(){};
    void compute(List3<FLOAT> &psrc, List3<FLOAT> &psink, FLOAT3 &delta, float _inv_eps2, FLOAT3 *pacc);
    void compute_fallback(List3<FLOAT> &psrc, List3<FLOAT> &psink, FLOAT3 &delta, float _inv_eps2, FLOAT3 *pacc);
};

#if defined(DIRECTSINGLESPLINE) || defined(TEST)

// TODO: If we're just using this for microstepping, no need for particle set offset
// Could make a templated version includes offsets for pencils
// TODO: Should sinks == sources?
// The particle arrays must be 64-byte aligned
void AVX512Direct::compute(List3<FLOAT> &psrc, List3<FLOAT> &psink, FLOAT3 &delta, float _inv_eps2, FLOAT3 *pacc){
    assert(delta.x == 0);
    assert(delta.y == 0);
    assert(delta.z == 0);
    int nsink = psink.N;
    int nsrc = psrc.N;

    // Constants
    AVX512_FLOATS inv_eps2 = AVX512_SET_FLOAT(_inv_eps2);
    AVX512_FLOATS small = AVX512_SET_FLOAT((FLOAT) 1e-32);
    AVX512_FLOATS m_fifteen = AVX512_SET_FLOAT((FLOAT) -15);
    AVX512_FLOATS six = AVX512_SET_FLOAT((FLOAT) 6);
    AVX512_FLOATS ten = AVX512_SET_FLOAT((FLOAT) 10);

    int nsrc_aligned = nsrc - (nsrc % AVX512_NVEC);
    
    for(int i = 0; i < nsink; i++){
        // Load a single sink: broadcast x,y,z into three vectors
        AVX512_FLOATS xsink = AVX512_SET_FLOAT(psink.X[i]);
        AVX512_FLOATS ysink = AVX512_SET_FLOAT(psink.Y[i]);
        AVX512_FLOATS zsink = AVX512_SET_FLOAT(psink.Z[i]);

        // TODO: do we need to read the accs and accumulate into that, or is setting them sufficient?
        AVX512_FLOATS xacc = AVX512_SET_FLOAT((FLOAT) 0);
        AVX512_FLOATS yacc = AVX512_SET_FLOAT((FLOAT) 0);
        AVX512_FLOATS zacc = AVX512_SET_FLOAT((FLOAT) 0);

        // Mask-off loop
        #pragma unroll (16)
        for(int j = 0; j <= nsrc_aligned - AVX512_NVEC; j += AVX512_NVEC){
            // Load exactly 16 sources
            AVX512_FLOATS xsource = AVX512_LOAD_FLOATS_ALIGNED(&(psrc.X[j]));
            AVX512_FLOATS ysource = AVX512_LOAD_FLOATS_ALIGNED(&(psrc.Y[j]));
            AVX512_FLOATS zsource = AVX512_LOAD_FLOATS_ALIGNED(&(psrc.Z[j]));

            // Compute dx,dy,dz = source - sink
            AVX512_FLOATS dx = AVX512_SUBTRACT_FLOATS(xsource, xsink);
            AVX512_FLOATS dy = AVX512_SUBTRACT_FLOATS(ysource, ysink);
            AVX512_FLOATS dz = AVX512_SUBTRACT_FLOATS(zsource, zsink);

            // dr2 = (dx*dx + dy*dy + dz*dz)*inv_eps2 + (FLOAT)1e-32;
            AVX512_FLOATS dr2 = AVX512_SQUARE_FLOAT(dx);
            dr2 = AVX512_FMA_ADD_FLOATS(dy, dy, dr2);
            dr2 = AVX512_FMA_ADD_FLOATS(dz, dz, dr2);
	    // TODO: Would do the FOF check here
            dr2 = AVX512_FMA_ADD_FLOATS(dr2, inv_eps2, small);
            //dr2 = AVX512_MULTIPLY_FLOATS(dr2, inv_eps2);

            // f = RSQRT(dr2);
            AVX512_FLOATS f = AVX512_RSQRT14_FLOAT(dr2);  // approximate rsqrt
            // In the future, we may want to add a NR iteration here
            // Uncomment to use slow, accurate rsqrt
            //for(int k = 0; k < AVX512_NVEC; k++)
            //    f[k] = 1.0f/sqrtf(dr2[k]);

            /* Compute 1/r^3 and the softened form; take the min with mask blends
            TODO: could change to (1/r^2)(1/r) via an inv and sqrt
            if(R <= 1){
                f = MIN(f*f*f, (-15*f + 6)*dr2 + 10);
            }
            else{
               f = f*f*f;
            }
             */
            AVX512_FLOATS f3 = AVX512_MULTIPLY_FLOATS(f, AVX512_SQUARE_FLOAT(f));
            AVX512_FLOATS softened = AVX512_FMA_ADD_FLOATS(m_fifteen, f, six);
            softened = AVX512_FMA_ADD_FLOATS(softened, dr2, ten);
            AVX512_MASK min_mask = AVX512_COMPARE_FLOATS(f3, softened, _CMP_LT_OS);  // ordered, signaling NaNs
            f = AVX512_BLEND_FLOATS_WITH_MASK(min_mask, softened, f3);

            // a.x -= f*dx;
            xacc = AVX512_FNMA_ADD_FLOATS(f, dx, xacc);
            yacc = AVX512_FNMA_ADD_FLOATS(f, dy, yacc);
            zacc = AVX512_FNMA_ADD_FLOATS(f, dz, zacc);
        }

        // Now do the last iteration with masks
        if(nsrc_aligned < nsrc){
            // Load up to 16 sources, masking extras as zeros
            AVX512_MASK m_mask_left = masks_per_misalignment_value_FLOAT[nsrc-nsrc_aligned];
            AVX512_FLOATS xsource = AVX512_MASKZ_LOAD_FLOATS_ALIGNED(m_mask_left, &(psrc.X[nsrc_aligned]));
            AVX512_FLOATS ysource = AVX512_MASKZ_LOAD_FLOATS_ALIGNED(m_mask_left, &(psrc.Y[nsrc_aligned]));
            AVX512_FLOATS zsource = AVX512_MASKZ_LOAD_FLOATS_ALIGNED(m_mask_left, &(psrc.Z[nsrc_aligned]));

            AVX512_FLOATS dx = AVX512_SUBTRACT_FLOATS(xsource, xsink);
            AVX512_FLOATS dy = AVX512_SUBTRACT_FLOATS(ysource, ysink);
            AVX512_FLOATS dz = AVX512_SUBTRACT_FLOATS(zsource, zsink);

            AVX512_FLOATS dr2 = AVX512_SQUARE_FLOAT(dx);
            dr2 = AVX512_FMA_ADD_FLOATS(dy, dy, dr2);
            dr2 = AVX512_FMA_ADD_FLOATS(dz, dz, dr2);
            dr2 = AVX512_FMA_ADD_FLOATS(dr2, inv_eps2, small);

            AVX512_FLOATS f = AVX512_RSQRT14_FLOAT(dr2);
            
            AVX512_FLOATS f3 = AVX512_MULTIPLY_FLOATS(f, AVX512_SQUARE_FLOAT(f));
            AVX512_FLOATS softened = AVX512_FMA_ADD_FLOATS(m_fifteen, f, six);
            softened = AVX512_FMA_ADD_FLOATS(softened, dr2, ten);
            
            AVX512_MASK min_mask = AVX512_COMPARE_FLOATS(f3, softened, _CMP_LT_OS);
            f = AVX512_BLEND_FLOATS_WITH_MASK(min_mask, softened, f3);

            // We're about to discard the current sink; start prefetching the next one
            PREFETCH(psink.X[i+1]);
            PREFETCH(psink.Y[i+1]);
            PREFETCH(psink.Z[i+1]);

            xacc = AVX512_MASK3_FNMA_ADD_FLOATS(f, dx, xacc, m_mask_left);
            yacc = AVX512_MASK3_FNMA_ADD_FLOATS(f, dy, yacc, m_mask_left);
            zacc = AVX512_MASK3_FNMA_ADD_FLOATS(f, dz, zacc, m_mask_left);
        }

        // Sum the accelerations for this sink and store
        pacc[i].x = AVX512_HORIZONTAL_SUM_FLOATS(xacc);
        pacc[i].y = AVX512_HORIZONTAL_SUM_FLOATS(yacc);
        pacc[i].z = AVX512_HORIZONTAL_SUM_FLOATS(zacc);
    }
}

// A plain C++ List3 implementation
void AVX512Direct::compute_fallback(List3<FLOAT> &psrc, List3<FLOAT> &psink, FLOAT3 &delta, float inv_eps2, FLOAT3 *pacc){
    for(int i = 0; i < psink.N; i++){
        FLOAT xsink = psink.X[i];
        FLOAT ysink = psink.Y[i];
        FLOAT zsink = psink.Z[i];

        for(int j = 0; j < psrc.N; j++){
            FLOAT dx = psrc.X[j] - xsink;
            FLOAT dy = psrc.Y[j] - ysink;
            FLOAT dz = psrc.Z[j] - zsink;
            
            FLOAT dr2 = (dx*dx + dy*dy + dz*dz)*inv_eps2 + (FLOAT)1e-32;
            FLOAT f = RSQRT(dr2);
            
            f = std::min(f*f*f, (-15*f + 6)*dr2 + 10);

            pacc[i].x -= f*dx;
            pacc[i].y -= f*dy;
            pacc[i].z -= f*dz;
        }
    }
}

#else
#error "AVX-512 directs currently only implemented for spline softening"
#endif

#ifdef TEST
#include <cstdlib>
#include <cstring>
#include <chrono>

#undef RSQRT
#define DIRECTSINGLESPLINE
#include "DirectCPUKernels.cc"

#include "avxdirectfloatNR.cc"

FLOAT3* cpu_directs(int nsrc, FLOAT3 *psrc, int nsink, FLOAT3 *psink, FLOAT3 &delta, float _inv_eps2){
    FLOAT3 *acc;
    posix_memalign((void **) &acc, 64, sizeof(FLOAT3)*nsink);

    for(int i = 0; i < nsink; i++){
        acc[i] = FLOAT3(0.);
        for(int j = 0; j < nsrc; j++){
            directkernel<0>(psink[i], psrc[j], acc[i], _inv_eps2);
        }
    }
    return acc;
}

void compare_acc(FLOAT3 *acc1, FLOAT3* acc2, int nacc, double rtol){
    int nbad = 0;
    double max_frac_diff = -1;
    for(int i = 0; i < nacc; i++){
        double3 a1(acc1[i]);
        double3 a2(acc2[i]);
        if (a1 == a2)
            continue;
        double frac_diff = (a1 - a2).norm()/(a1 + a2).norm();
        max_frac_diff = std::max(max_frac_diff, frac_diff);
        if(frac_diff > rtol || !std::isfinite(frac_diff)){
            nbad++;
            // fmt::print("{}\n", acc1[i]);
            // fmt::print("{}\n", acc2[i]);
        }
    }
    fmt::print("\t>>> {:d} ({:.2f}%) mismatched accels\n", nbad, (FLOAT) nbad/nacc*100);
    fmt::print("\t>>> Max frac error: {:.2g} \n", max_frac_diff);
}

void report(const std::string &prefix, int64_t ndirect, std::chrono::duration<double> elapsed){
    fmt::print("{} directs time: {} sec\n", prefix, elapsed.count());
    fmt::print("\t{} billion directs per second\n", ndirect/1e9/elapsed.count());
}

int main(int argc, char **argv){
    double rtol = 1e-4;  // acc comparison tolerance
    FLOAT eps = 1e-2;  // softening length

    int nsrc = argc < 2 ? 30000 : atoi(argv[1]);
    List3<FLOAT> sources(nsrc);
    for (uint64_t i = 0; i < nsrc; i++){
	    sources.X[i] = (FLOAT) rand()/RAND_MAX;
	    sources.Y[i] = (FLOAT) rand()/RAND_MAX;
	    sources.Z[i] = (FLOAT) rand()/RAND_MAX;

        /*sources.X[i] = 0.;
        sources.Y[i] = 0.;
        sources.Z[i] = 1.;*/
	}
    
    int nsink = nsrc;
    List3<FLOAT> sinks = sources;
    sinks.owndata = false;
    /*List3<FLOAT> sinks(nsink);
    for(uint64_t j = 0; j < nsink; j++){
	    sinks.X[j] = rand()/RAND_MAX;
	    sinks.Y[j] = rand()/RAND_MAX;
	    sinks.Z[j] = rand()/RAND_MAX;

        //sinks.X[j] = 0.;
        //sinks.Y[j] = 0.;
        //sinks.Z[j] = 0.;
	}*/

    fmt::print("eps: {:g}; rtol: {:g}\n", eps, rtol);
    fmt::print("Finished particle generation.  Starting directs.\n");
    AVX512Direct dd;

    /****************************************/
    // AVX-512 version

    // Do we want acc to be a float3 or list3?
    FLOAT3 *acc;
    posix_memalign((void **) &acc, 64, sizeof(FLOAT3)*nsink);
    memset(acc, 0, sizeof(FLOAT3)*nsink);

    FLOAT3 delta(0,0,0);

    auto begin = std::chrono::steady_clock::now();
    dd.compute(sources, sources, delta, 1./(eps*eps), acc);
    auto end= std::chrono::steady_clock::now();
    report("AVX-512", (int64_t) nsink*nsrc, end-begin);

    // multiply the final prefactor
    //for(int i = 0; i < nsink; i++)
    //    acc[i] /= std::pow(eps, 3);

    //for(int i = 0; i < 10; i++)
    //    fmt::print("({:g}, {:g}, {:g})\n", acc[i].x, acc[i].y, acc[i].z);

    // CPU FLOAT3 version
    FLOAT3 *f3sources;
    posix_memalign((void **) &f3sources, 64, sizeof(FLOAT3)*nsrc);
    FLOAT3 *f3sinks = f3sources;
    //posix_memalign((void **) &f3sinks, 64, sizeof(FLOAT3)*nsink);
    for(int i = 0; i < nsrc; i++){
        f3sources[i].x = sources.X[i];
        f3sources[i].y = sources.Y[i];
        f3sources[i].z = sources.Z[i];
    }

    begin = std::chrono::steady_clock::now();
    FLOAT3 *cpuacc = cpu_directs(nsrc, f3sources, nsink, f3sinks, delta, 1./(eps*eps));
    end = std::chrono::steady_clock::now();
    report("CPU FLOAT3", (int64_t) nsink*nsrc, end-begin);
    compare_acc(acc, cpuacc, nsink, rtol);

    /****************************************/
    // CPU AVX version (FLOAT3)
    FLOAT3 *acc_avx_plummer;
    posix_memalign((void **) &acc_avx_plummer, 64, sizeof(FLOAT3)*nsink);
    memset(acc_avx_plummer, 0, sizeof(FLOAT3)*nsink);

    AVXDirectFloatNR directfloatnr(nsrc+4096);
    begin = std::chrono::steady_clock::now();
    directfloatnr.compute(nsrc, f3sources, nsink, f3sinks, delta, eps, acc_avx_plummer);
    end = std::chrono::steady_clock::now();
    report("Plummer AVX", (int64_t) nsink*nsrc, end-begin);
    compare_acc(acc, acc_avx_plummer, nsink, rtol);

    /****************************************/
    // CPU List3 version
    FLOAT3 *acc_list3;
    posix_memalign((void **) &acc_list3, 64, sizeof(FLOAT3)*nsink);
    memset(acc_list3, 0, sizeof(FLOAT3)*nsink);

    begin = std::chrono::steady_clock::now();
    dd.compute_fallback(sources, sources, delta, 1./(eps*eps), acc_list3);
    end = std::chrono::steady_clock::now();
    report("CPU List3", (int64_t) nsink*nsrc, end-begin);
    compare_acc(acc, acc_list3, nsink, rtol);

    return 0;
}
#endif
