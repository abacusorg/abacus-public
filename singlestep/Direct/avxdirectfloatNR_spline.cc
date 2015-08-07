AVXDirectFloatNR::AVXDirectFloatNR(int maxsrc) : maxsrc(maxsrc) {
    assert( posix_memalign((void **)&jpdata,  64, sizeof(jpstruct<float>)*maxsrc)     == 0 );
    assert( posix_memalign((void **)&ipdata,  64, sizeof(ipstruct<float,8>))  == 0 );
    assert( posix_memalign((void **)&deltas,  64, sizeof(ipstruct<float,8>))  == 0 );
    assert( posix_memalign((void **)&accdata, 64, sizeof(apstruct<float,8>)) == 0 );
}

AVXDirectFloatNR::~AVXDirectFloatNR() {
    free(jpdata); free(ipdata); free(deltas); free(accdata);
}


void AVXDirectFloatNR::storejpdata(int nsrc, ThreeVector<float> *psrc) {
    for(int j=0; j<nsrc; j++) {
        jpdata[j].x = psrc[j].x;
        jpdata[j].y = psrc[j].y;
        jpdata[j].z = psrc[j].z;
        jpdata[j].m = 1.0f;
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
    
    float conpot = sqrt(2.0);


    float conacc = conpot*conpot*conpot;
    
    eps2 = sqrt(eps2);

    for(int i=0; i<8; i++) {
        deltas->x[i] = delta.x;
        deltas->y[i] = delta.y;
        deltas->z[i] = delta.z;
    }

    storejpdata(nsrc, psrc);
    
    for(int p=0; p<nsink; p+=8) {

        int e = 8; if(p+e>nsink) e = nsink-p;
        for(int j=0; j<e; j++) {
            ipdata->x[j] = psink[p+j].x;
            ipdata->y[j] = psink[p+j].y;
            ipdata->z[j] = psink[p+j].z;
            ipdata->eps2[j] = eps2;
        }
        for(int j = e; j < 8; j++)
            ipdata->eps2[j] = eps2;

        // applied all of the sources to these sinks
        KernelAccPot(ipdata, jpdata, nsrc, deltas, eps2, accdata);

        for(int j=0; j<e; j++) {
            pacc[p+j].x -= accdata->x[j];
            pacc[p+j].y -= accdata->y[j];
            pacc[p+j].z -= accdata->z[j];
        }
    }
}

#define AX    YMM08
#define AY    YMM09
#define AZ    YMM10
#define UN   YMM11

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

// do one xj on eight xi's per loop
void AVXDirectFloatNR::KernelAccPot(ipstruct<float,8> *ipdata, jpstruct<float> *jpdata, int nsrc, 
                                     ipstruct<float,8> *deltas, float eps2, 
                                     apstruct<float,8> *accdata) {
    int j;

    
    float spline_r1_0 ALIGN64;
    float spline_r1_1 ALIGN64;
    float spline_r1_2 ALIGN64;
    float spline_r1_3 ALIGN64;
    float spline_r1_4 ALIGN64;

    float spline_r0_0 ALIGN64;
    float spline_r0_1 ALIGN64;

    float d_one ALIGN64;
    float d_two ALIGN64;
    float d_three ALIGN64;

    d_one = 1.0f;
    d_two = 2.0f;
    d_three = 3.0f;

    spline_r1_0 = -1.0f/15.0f;
    spline_r1_1 = 8.0f/3.0f;
    spline_r1_2 = -3.0f;
    spline_r1_3 = 6.0f/5.0f;
    spline_r1_4 = -1.0f/6.0f;

    spline_r0_0 = 4.0f/3.0f;
    spline_r0_1 = 1.0f/2.0f;

    PREFETCH(*jpdata);

    VZEROALL; // to zero out acc registers for accumulation


    VLOADPS(*ipdata->x, XI);    // load 4 xi's 
    VLOADPS(*ipdata->y, YI);
    VLOADPS(*ipdata->z, ZI);

    VADDPS_M(*deltas->x, XI, XI); // add deltas to sink position
    VADDPS_M(*deltas->y, YI, YI); 
    VADDPS_M(*deltas->z, ZI, ZI); 

    VBROADCASTSS(jpdata->x, X2); // load 1 xj into four copies
    VBROADCASTSS(jpdata->y, Y2);
    VBROADCASTSS(jpdata->z, J2);
    VBROADCASTSS(jpdata->m, MJ);
    jpdata++;

    for(j=0; j<nsrc; j++) {
        VLOADPS(*ipdata->eps2, EPS2);

        PREFETCH(*jpdata);  // does this do anything -- unlikely!

        VSUBPS(XI, X2, DX);
        VSUBPS(ZI, J2, DZ);
        VSUBPS(YI, Y2, DY);

        VMULPS(DX, DX, X2);      // X2 = DX^2
        VMULPS(DZ, DZ, J2);
        VMULPS(DY, DY, Y2);

        VADDPS(X2, J2, J2);   // J2 = X2 + J2 = DX^2 + DZ^2
        VADDPS(J2, Y2, Y2);   // Y2 = Y2 + J2 = DX^2 + DY^2 + DZ^2 = R^2

        VSQRTPS(Y2,Y2); //Y2 = sqrt(R^2) = R
        VDIVPS(EPS2,Y2,X2); //x2 = R/eps = u
        
        //prune out any u > 3
        VBROADCASTSS(d_three,UN); //UN = 3
        VCMPPS(X2,UN,J1,0x11);
        //X2 = (u >= 3): 3 : u;
        VBLENDVPS(J1,UN,X2,X2); 

        
        //we will now accumulate the two spline factors into J1 and J2
        //we count up through the powers of u in UN
        VMULPS(X2,X2,UN);//UN = u**2
        
        VBROADCASTSS(spline_r1_0,J1); //J1 = -1/14
        VBROADCASTSS(spline_r1_3,EPS2); //eps2 = 6/5
        VBROADCASTSS(spline_r0_0,J2); //J2 = 4/3

        VMULPS(UN,EPS2,EPS2); // EPS2 = 6/5 u**2 
        VSUBPS(EPS2,J2,J2); //J2 = 4/3 - 6/5 u**2 

        VMULPS(UN,X2,UN); //UN = u**3
        VBROADCASTSS(spline_r1_1,EPS2); //eps2 = 8/3
        VMULPS(EPS2,UN,EPS2); //eps2 = 8/3 u**3
        VADDPS(EPS2,J1,J1); //J1 = -1/14 + 8/3 u**3

        VBROADCASTSS(spline_r0_1,EPS2);//eps2 = 1/2
        VMULPS(EPS2,UN,EPS2);//eps2 = 1/2 u**3
        VADDPS(EPS2,J2,J2); //J2 = 4/3 - 6/5 u**2 + 1/2 u**3
        

        VMULPS(UN,X2,UN);//UN = u**4
        VBROADCASTSS(spline_r1_2,EPS2);//EPS2 = -3
        VMULPS(UN,EPS2,EPS2); //EPS2 = -3 u**4
        VADDPS(EPS2,J1,J1);// J1 = -1/14 + 8/3 u**3 - 3 u**4
        
        VMULPS(UN,X2,UN); //UN = u**5
        VBROADCASTSS(spline_r1_3,EPS2);//EPS2 = 6/5
        VMULPS(UN,EPS2,EPS2); //EPS2 = 6/5 u**5
        VADDPS(EPS2,J1,J1); //J1 = -1/14 + 8/3 u**3 - 3 u**4 +  6/5 u**5
        
        VMULPS(UN,X2,UN); //UN = u**5
        VBROADCASTSS(spline_r1_4,EPS2);//EPS2 = 1/6
        VMULPS(UN,EPS2,EPS2); //EPS2 = 1/6 u**6
        VADDPS(EPS2,J1,J1); //J1 = -1/14 + 8/3 u**3 - 3 u**4 +  6/5 u**5 - 1/6 u**6

        VBROADCASTSS(d_one,EPS2); //EPS2 = 1.0;
        VCMPPS(EPS2, X2,UN, 0x11); //UN = (u < 1.0)

        //choose whether to use J1 or J2 as the softening factor;
        VBLENDVPS(UN,J2,J1,J1);

        VBROADCASTSS(d_two, J2);//J2 = 2.0

        VCMPPS(J2,X2,UN,0x11); //UN = (u < 2.0)

        VBLENDVPS(UN, J1, EPS2, J1); //J1 = UN ? J1 : 1.0;

        //J1 is now correctly set to be the softening factor
        
        VSUBPS(J2,J2,J2); //J2 = 0;

        VCMPPS(EPS2,X2,UN,0x11); //UN = (u <1.0)

        VLOADPS(*ipdata->eps2, EPS2);
        VBLENDVPS(UN,EPS2,Y2,J2); //J2 = (Y2 < 1.0) ? Y2 : eps;



        VDIVPS(J2,J1,J1); //J1 = spf/r
        VDIVPS(J2,J1,J1); //J1 = spf/r**2
        VDIVPS(J2,J1,Y2); //J1 = spf/r**3

        VMULPS(MJ, Y2, Y2);   // Y2 = m/R^3

        VMULPS(Y2, DX, DX);   // DX = DX m / R^3 = dx m / R^3
        VMULPS(Y2, DY, DY);
        VMULPS(Y2, DZ, DZ);

        VBROADCASTSS(jpdata->x, X2);
        VBROADCASTSS(jpdata->y, Y2);
        VBROADCASTSS(jpdata->z, J2);
        VBROADCASTSS(jpdata->m, MJ);
        jpdata++;

        VADDPS(DX, AX, AX);
        VADDPS(DY, AY, AY);
        VADDPS(DZ, AZ, AZ);
    }

    VSTORPS(AX, *accdata->x);
    assert(isfinite(*accdata->x));
    VSTORPS(AY, *accdata->y);
    assert(isfinite(*accdata->y));
    VSTORPS(AZ, *accdata->z);
    assert(isfinite(*accdata->z));
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
