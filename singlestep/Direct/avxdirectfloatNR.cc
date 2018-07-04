#include "avxdirectfloatNR.h"

AVXDirectFloatNR::AVXDirectFloatNR(int maxsrc) : maxsrc(maxsrc) {
    assert( posix_memalign((void **)&jpdata,  64, sizeof(jpstruct<float>)*maxsrc)     == 0 );
    assert( posix_memalign((void **)&ipdata,  64, sizeof(ipstruct<float,4>))  == 0 );
    assert( posix_memalign((void **)&deltas,  64, sizeof(ipstruct<float,4>))  == 0 );
    assert( posix_memalign((void **)&accdata, 64, sizeof(apstruct<float,4>)) == 0 );

    for(int i=0;i<maxsrc;i++) {
        jpdata[i].x = 0;
        jpdata[i].y = 0;
        jpdata[i].z = 0;
        jpdata[i].m = 0;
    }
}

AVXDirectFloatNR::~AVXDirectFloatNR() {
    free(jpdata); free(ipdata); free(deltas); free(accdata);
}


void AVXDirectFloatNR::storejpdata(int nsrc, ThreeVector<float> *psrc) {
    for(int j=0; j<nsrc; j++) {
        __m128 pf = {psrc[j].x, psrc[j].y, psrc[j].z, (float)1};
        *(__m128 *)(jpdata+j) = pf;
    }
    int e = 4-nsrc%4;
    assert((nsrc+e)%4==0);
    for(int j=nsrc;j<nsrc+e;j++) {
        jpdata[j].x =0;
        jpdata[j].y = 0;
        jpdata[j].z = 0;
        jpdata[j].m = 0;
    }
    
}

/*
  Apply a particleset src of length nsrc to a particleset sink of length nsink with coordinate systems separated
  by delta, to produce a particleset of acceleration and potential.
 */
void AVXDirectFloatNR::compute(int nsrc, ThreeVector<float> *psrc,
                           int nsink, ThreeVector<float> *psink,
                           float3 &delta, float eps2,
                           ThreeVector<float> *pacc) {
    
    float conpot = -0.5;
    //assert(nsrc + 4 <= maxsrc);

    float conacc = conpot*conpot*conpot;
    
    for(int i=0; i<4; i++) {
        deltas->x[i] = (float)delta.x;
        deltas->y[i] = (float)delta.y;
        deltas->z[i] = (float)delta.z;
    }

#if 0
    for(int i=0;i<4;i++) {
        jpdata[i].x = 0;
        jpdata[i].y = 0;
        jpdata[i].z = 0;
        jpdata[i].m = 0;
    }
#endif

    storejpdata(nsrc,psrc); 
 
    
    for(int p=0; p<nsink; p+=4) {

        int e = 4; if(p+e>nsink) e = nsink-p;
        for(int j=0; j<e; j++) {
            ipdata->x[j] = psink[p+j].x;
            ipdata->y[j] = psink[p+j].y;
            ipdata->z[j] = psink[p+j].z;
            ipdata->eps2[j] = eps2;
        }
        for(int j =e; j <4; j++){
            ipdata->x[j] = 0.0f;
            ipdata->y[j] = 0.0f;
            ipdata->z[j] = 0.0f;
            ipdata->eps2[j] = 1.0f;
        }

        // applied all of the sources to these sinks
        KernelAccPot(ipdata, jpdata, nsrc, deltas, eps2, accdata);

        for(int j=0; j<e; j++) {
            pacc[p+j].x -= accdata->x[j]*conacc;
            pacc[p+j].y -= accdata->y[j]*conacc;
            pacc[p+j].z -= accdata->z[j]*conacc;
        }
    }
}

#define AX    YMM08
#define AY    YMM09
#define AZ    YMM10
#define PHI   YMM11

#define DX    YMM12
#define DY    YMM13
#define DZ    YMM14
#define MJ    YMM07

#define J1    YMM00
#define J2    YMM01
#define X2    YMM02
#define xX2   XMM02
#define Y2    YMM03

#define XI    YMM04
#define YI    YMM05
#define ZI    YMM06
#define EPS2  YMM15

#include "avxsseabrev.h"
#define ALIGN64 __attribute__ ((aligned(64)))

void AVXDirectFloatNR::KernelAccPot(ipstruct<float,4> *ipdata, jpstruct<float> *jpdata, int nsrc, 
                                    ipstruct<float,4> *deltas, float eps2,
                                    apstruct<float,4> *accdata) {
    int j;

    float three ALIGN64;
    three = -3.0;

    PREFETCH(jpdata[0]);

    VZEROALL; // to zero out acc registers for accumulation

    LOADPS(*ipdata->x, XMM04);    // load 4 xi's into XMM04 (aliased into XI)
    LOADPS(*ipdata->y, XMM05);
    LOADPS(*ipdata->z, XMM06);
    LOADPS(*ipdata->eps2, XMM15);

    VADDPS_M(*deltas->x, XMM04, XMM04); // add deltas to sink position
    VADDPS_M(*deltas->y, XMM05, XMM05); 
    VADDPS_M(*deltas->z, XMM06, XMM06); 

    VPERM2F128(XI, XI, XI, 0x00); // copy low 4 sp into high 4 sp for xi, yi, zi, and eps
    VPERM2F128(YI, YI, YI, 0x00);
    VPERM2F128(ZI, ZI, ZI, 0x00);
    VPERM2F128(EPS2, EPS2, EPS2, 0x00);
    // so that, e.g., XI = (xi0, xi1, xi2, xi3, xi0, xi1, xi2, xi3) and so on...


    VLOADPS(*(jpdata), J1);
    VLOADPS(*(jpdata+2), J2);

    jpdata += 4;

    VSHUFPS(J1, J1, X2, 0x00);
    VSHUFPS(J1, J1, Y2, 0x55);
    VSHUFPS(J1, J1, MJ, 0xff);
    VSHUFPS(J1, J1, J1, 0xaa);

    for(j = 0 ; j < nsrc; j += 4) {

        VSUBPS(XI, X2, DX);
        VSUBPS(YI, Y2, DY);
        VSUBPS(ZI, J1, DZ);

        VMULPS(DX, DX, X2);
        VMULPS(DZ, DZ, J1);
        VMULPS(DY, DY, Y2);

        VADDPS(J1, X2, X2);
        VADDPS(EPS2, Y2, Y2);
        VADDPS(Y2, X2, Y2);

        VBROADCASTSS(three, J1);
        VRSQRTPS(Y2, X2);

        VMULPS(X2, Y2, Y2);
        VMULPS(X2, Y2, Y2);
        VADDPS(J1, Y2, Y2);
        VMULPS(Y2, X2, X2);

        VLOADPS(*(jpdata), J1);

        VMULPS(X2, MJ, MJ);
        VMULPS(X2, X2, Y2);

        VMULPS(MJ, Y2, Y2);
        //        VSUBPS(MJ, PHI, PHI);

        VMULPS(Y2, DX, DX);
        VMULPS(Y2, DY, DY);
        VMULPS(Y2, DZ, DZ);

        VSHUFPS(J2, J2, X2, 0x00);
        VSHUFPS(J2, J2, MJ, 0xff);
        VSHUFPS(J2, J2, Y2, 0x55);
        VSHUFPS(J2, J2, J2, 0xaa);

        VADDPS(DX, AX, AX);
        VADDPS(DY, AY, AY);
        VADDPS(DZ, AZ, AZ);

        VSUBPS(XI, X2, DX);
        VSUBPS(YI, Y2, DY);
        VSUBPS(ZI, J2, DZ);

        VMULPS(DX, DX, X2);
        VMULPS(DZ, DZ, J2);
        VMULPS(DY, DY, Y2);

        VADDPS(J2, X2, X2);
        VADDPS(EPS2, Y2, Y2);
        VADDPS(Y2, X2, Y2);

        VBROADCASTSS(three, J2);
        VRSQRTPS(Y2, X2);

        VMULPS(X2, Y2, Y2);
        VMULPS(X2, Y2, Y2);
        VADDPS(J2, Y2, Y2);
        VMULPS(Y2, X2, X2);

        VLOADPS(*(jpdata+2), J2);    

        VMULPS(X2, MJ, MJ);
        VMULPS(X2, X2, Y2);

        jpdata += 4;
        PREFETCH(*(jpdata));

        VMULPS(MJ, Y2, Y2);
        VSUBPS(MJ, PHI, PHI);

        VMULPS(Y2, DX, DX);
        VMULPS(Y2, DY, DY);
        VMULPS(Y2, DZ, DZ);

        VSHUFPS(J1, J1, X2, 0x00);
        VSHUFPS(J1, J1, MJ, 0xff);
        VSHUFPS(J1, J1, Y2, 0x55);
        VSHUFPS(J1, J1, J1, 0xaa);

        VADDPS(DX, AX, AX);
        VADDPS(DY, AY, AY);
        VADDPS(DZ, AZ, AZ);
    }

    VEXTRACTF128(AX, XMM00, 0x01);   // extract upper 128 from AX into XMM00
    VADDPS(AX, YMM00, AX);           // add to lower 128
    VEXTRACTF128(AY, XMM01, 0x01);
    VADDPS(AY, YMM01, AY);
    VEXTRACTF128(AZ, XMM02, 0x01);
    VADDPS(AZ, YMM02, AZ);
    VEXTRACTF128(PHI, XMM03, 0x01);
    VADDPS(PHI, YMM03, PHI);

    STORPS(XMM08, *accdata->x);
    STORPS(XMM09, *accdata->y);
    STORPS(XMM10, *accdata->z);
    //    STORPS(XMM11, *accdata->pot);

}
#undef AX
#undef AY
#undef AZ
#undef PHI
#undef DX
#undef DY
#undef DZ
#undef MJ
#undef J1
#undef J2
#undef X2
#undef Y2
#undef XI
#undef YI
#undef ZI
#undef EPS2
