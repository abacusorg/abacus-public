// AVX-512 implementations of our directs kernels
// This implementation is based on Manodeep Sinha's AVX-512 pair counting implementation for Corrfunc

#ifdef TEST
#define FLOAT float
#define FLOAT3 float3
#endif

#include "threevector.hh"
#include "StructureOfLists.cc"
#include "avx512_calls.c"

class AVX512Direct {
public:
    AVX512Direct(){};
    ~AVX512Direct(){};
    void compute(int nsrc, List3<FLOAT> &psrc, int nsink, List3<FLOAT> psink, FLOAT3 &delta, float _inv_eps2,
                 FLOAT3 *pacc);
};

#if defined(DIRECTSINGLESPLINE) || defined(TEST)

// TODO: If we're just using this for microstepping, no need for particle set offset
// Could make a templated version includes offsets for pencils
// TODO: could make a version that takes the position transpose instead of FLOAT3s, at least for the sources
void AVX512Direct::compute(int nsrc, List3<FLOAT> &psrc, int nsink, List3<FLOAT> psink, FLOAT3 &delta, float _inv_eps2,
                 FLOAT3 *pacc){
    assert(delta.x == 0);
    assert(delta.y == 0);
    assert(delta.z == 0);

    // Constants
    AVX512_FLOATS inv_eps2 = AVX512_SET_FLOAT(_inv_eps2);
    AVX512_FLOATS small = AVX512_SET_FLOAT((FLOAT) 1e-32);
    AVX512_FLOATS m_fifteen = AVX512_SET_FLOAT((FLOAT) -15);
    AVX512_FLOATS six = AVX512_SET_FLOAT((FLOAT) 6);
    AVX512_FLOATS ten = AVX512_SET_FLOAT((FLOAT) 10);
    
    for(int i = 0; i < nsink; i++){
        // Load a single sink: broadcast x,y,z into three vectors
        AVX512_FLOATS xsink = AVX512_SET_FLOAT(psink.X[i]);
        AVX512_FLOATS ysink = AVX512_SET_FLOAT(psink.Y[i]);
        AVX512_FLOATS zsink = AVX512_SET_FLOAT(psink.Z[i]);

        // TODO: do we need to read the accs and accumulate into that, or is setting them sufficient?
        AVX512_FLOATS xacc = AVX512_SET_FLOAT((FLOAT) 0);
        AVX512_FLOATS yacc = AVX512_SET_FLOAT((FLOAT) 0);
        AVX512_FLOATS zacc = AVX512_SET_FLOAT((FLOAT) 0);

        for(int j = 0; j < nsrc; j += AVX512_NVEC){
            // Load up to 16 sources, masking extras as zeros
            AVX512_MASK m_mask_left = (nsrc - j) >= AVX512_NVEC ? ~0:masks_per_misalignment_value_FLOAT[nsrc-j];
            // TODO: maybe we can guarantee alignment?
            // TODO: Is it faster to load with masks or do a "cleanup" loop after?
            AVX512_FLOATS xsource = AVX512_MASKZ_LOAD_FLOATS_ALIGNED(m_mask_left, &(psrc.X[j]));
            AVX512_FLOATS ysource = AVX512_MASKZ_LOAD_FLOATS_ALIGNED(m_mask_left, &(psrc.Y[j]));
            AVX512_FLOATS zsource = AVX512_MASKZ_LOAD_FLOATS_ALIGNED(m_mask_left, &(psrc.Z[j]));

            // Compute dx,dy,dz = source - sink
            AVX512_FLOATS dx = AVX512_SUBTRACT_FLOATS(xsink, xsource);
            AVX512_FLOATS dy = AVX512_SUBTRACT_FLOATS(ysink, ysource);
            AVX512_FLOATS dz = AVX512_SUBTRACT_FLOATS(zsink, zsource);

            // dr2 = (dx*dx + dy*dy + dz*dz)*inv_eps2 + (FLOAT)1e-32;
            // TODO: combine these?
            AVX512_FLOATS dr2 = AVX512_SQUARE_FLOAT(dx);
            dr2 = AVX512_FMA_ADD_FLOATS(dy, dy, dr2);
            dr2 = AVX512_FMA_ADD_FLOATS(dz, dz, dr2);
            dr2 = AVX512_FMA_ADD_FLOATS(dr2, inv_eps2, small);

            // f = RSQRT(dr2);
            AVX512_FLOATS f = AVX512_RSQRT14_FLOAT(dr2);  // approximate rsqrt
            printf("%f \n", f[0]);

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
            // TODO: is this inefficient because we're "creating" a new vector instead of specifying a destination?
            xacc = AVX512_MASK3_FNMA_ADD_FLOATS(f, dx, xacc, m_mask_left);
            yacc = AVX512_MASK3_FNMA_ADD_FLOATS(f, dy, yacc, m_mask_left);
            zacc = AVX512_MASK3_FNMA_ADD_FLOATS(f, dz, zacc, m_mask_left);
        }

        // Sum the accelerations for this sink and store
        // TODO: can we make this aligned? 
        pacc[i].x = AVX512_HORIZONTAL_SUM_FLOATS(xacc);
        pacc[i].y = AVX512_HORIZONTAL_SUM_FLOATS(yacc);
        pacc[i].z = AVX512_HORIZONTAL_SUM_FLOATS(zacc);
    }
}

#else
#error "AVX-512 directs currently only implemented for spline softening"
#endif

#ifdef TEST
int main(void){
    AVX512Direct dd;

    int nsrc = 1;
    List3<FLOAT> sources(nsrc);
    sources.X[0] = 0;
    sources.Y[0] = 0;
    sources.Z[0] = 1;
    
    int nsink = 1;
    List3<FLOAT> sinks(nsink);
    sinks.X[0] = 0;
    sinks.Y[0] = 0;
    sinks.Z[0] = 0;
    
    // Do we want acc to be a float3 or list3?
    FLOAT3 acc[1] = {FLOAT3(0.)};

    FLOAT eps = .5;
    FLOAT3 delta(0,0,0);

    dd.compute(nsrc, sources, nsink, sinks, delta, 1./(eps*eps), acc);

    // multiply the final prefactor
    for(int i = 0; i < nsink; i++)
        acc[i] /= std::pow(eps, 3);

    printf("(%f, %f, %f)\n", acc[0].x, acc[0].y, acc[0].z);

    return 0;
}
#endif