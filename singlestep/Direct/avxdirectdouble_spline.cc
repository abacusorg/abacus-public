// Copyright 2012-2025 The Abacus Developers
// SPDX-License-Identifier: GPL-3.0-or-later

AVXDirectDouble::AVXDirectDouble(int maxsrc) : maxsrc(maxsrc) {
    assert( posix_memalign((void **)&jpdata,  64, sizeof(jpstruct<double>)*maxsrc)     == 0 );
    assert( posix_memalign((void **)&ipdata,  64, sizeof(ipstruct<double,4>))  == 0 );
    assert( posix_memalign((void **)&deltas,  64, sizeof(ipstruct<double,4>))  == 0 );
    assert( posix_memalign((void **)&accdata, 64, sizeof(apstruct<double,4>)) == 0 );
}

AVXDirectDouble::~AVXDirectDouble() {
    free(jpdata); free(ipdata); free(deltas); free(accdata);
}


void AVXDirectDouble::storejpdata(int nsrc, ThreeVector<double> *psrc) {
    for(int j=0; j<nsrc; j++) {
        __m256d pd = {psrc[j].x, psrc[j].y, psrc[j].z, (double)1};
        *(__m256d *)(jpdata+j) = pd;
    }
}


/*
  Apply a particleset src of length nsrc to a particleset sink of length nsink with coordinate systems separated
  by delta, to produce a particleset of acceleration and potential.
 */

void AVXDirectDouble::compute(int nsrc, ThreeVector<double> *psrc, 
                           int nsink, ThreeVector<double> *psink, 
                           double3 &delta, double eps_inv, 
                           ThreeVector<double> *pacc) {
    
    double conpot = sqrt(2.0);
    double eps = 1./eps_inv;  // Should change kernels to accept eps_inv

    double conacc = conpot*conpot*conpot;
    
    for(int i=0; i<4; i++) {
        deltas->x[i] = delta.x;
        deltas->y[i] = delta.y;
        deltas->z[i] = delta.z;
    }

    storejpdata(nsrc, psrc);
    
    for(int p=0; p<nsink; p+=4) {

        int e = 4; if(p+e>nsink) e = nsink-p;
        for(int j=0; j<e; j++) {
            ipdata->x[j] = psink[p+j].x;
            ipdata->y[j] = psink[p+j].y;
            ipdata->z[j] = psink[p+j].z;
            ipdata->eps2[j] = eps;
        }
        for(int j = e; j < 4; j++)
            ipdata->eps2[j] = eps;

        // applied all of the sources to these sinks
        KernelAccPot(ipdata, jpdata, nsrc, deltas, eps, accdata);

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

// do one xj on four xi's per loop
void AVXDirectDouble::KernelAccPot(ipstruct<double,4> *ipdata, jpstruct<double> *jpdata, int nsrc, 
                                     ipstruct<double,4> *deltas, double eps2, 
                                     apstruct<double,4> *accdata) {
    int j;

    
    double spline_r1_0 ALIGN64;
    double spline_r1_1 ALIGN64;
    double spline_r1_2 ALIGN64;
    double spline_r1_3 ALIGN64;
    double spline_r1_4 ALIGN64;

    double spline_r0_0 ALIGN64;
    double spline_r0_1 ALIGN64;

    double d_one ALIGN64;
    double d_two ALIGN64;
    double d_three ALIGN64;
    d_one = 1.0;
    d_two = 2.0;
    d_three = 3.0;

    spline_r1_0 = -1.0/15.0;
    spline_r1_1 = 8.0/3.0;
    spline_r1_2 = -3.0;
    spline_r1_3 = 6.0/5.0;
    spline_r1_4 = -1.0/6.0;

    spline_r0_0 = 4.0/3.0;
    spline_r0_1 = 1.0/2.0;

    PREFETCH(*jpdata);

    VZEROALL; // to zero out acc registers for accumulation


    VLOADPD(*ipdata->x, XI);    // load 4 xi's 
    VLOADPD(*ipdata->y, YI);
    VLOADPD(*ipdata->z, ZI);

    VADDPD_M(*deltas->x, XI, XI); // add deltas to sink position
    VADDPD_M(*deltas->y, YI, YI); 
    VADDPD_M(*deltas->z, ZI, ZI); 

    VBROADCASTSD(jpdata->x, X2); // load 1 xj into four copies
    VBROADCASTSD(jpdata->y, Y2);
    VBROADCASTSD(jpdata->z, J2);
    VBROADCASTSD(jpdata->m, MJ);
    jpdata++;

    for(j=0; j<nsrc; j++) {
        VLOADPD(*ipdata->eps2, EPS2);

        PREFETCH(*jpdata);  // does this do anything -- unlikely!

        VSUBPD(XI, X2, DX);
        VSUBPD(ZI, J2, DZ);
        VSUBPD(YI, Y2, DY);

        VMULPD(DX, DX, X2);      // X2 = DX^2
        VMULPD(DZ, DZ, J2);
        VMULPD(DY, DY, Y2);

        VADDPD(X2, J2, J2);   // J2 = X2 + J2 = DX^2 + DZ^2
        VADDPD(J2, Y2, Y2);   // Y2 = Y2 + J2 = DX^2 + DY^2 + DZ^2 = R^2

        VSQRTPD(Y2,Y2); //Y2 = sqrt(R^2) = R
        VDIVPD(EPS2,Y2,X2); //x2 = R/eps = u
        //prune out any u > 3
        VBROADCASTSD(d_three,UN); //UN = 3
        VCMPPD(X2,UN,J1,0x11);
        //X2 = (u >= 3): 3 : u;
        VBLENDVPD(J1,UN,X2,X2);
        
        //we will now accumulate the two spline factors into J1 and J2
        //we count up through the powers of u in UN
        VMULPD(X2,X2,UN);//UN = u**2
        
        VBROADCASTSD(spline_r1_0,J1); //J1 = -1/14
        VBROADCASTSD(spline_r1_3,EPS2); //eps2 = 6/5
        VBROADCASTSD(spline_r0_0,J2); //J2 = 4/3

        VMULPD(UN,EPS2,EPS2); // EPS2 = 6/5 u**2 
        VSUBPD(EPS2,J2,J2); //J2 = 4/3 - 6/5 u**2 

        VMULPD(UN,X2,UN); //UN = u**3
        VBROADCASTSD(spline_r1_1,EPS2); //eps2 = 8/3
        VMULPD(EPS2,UN,EPS2); //eps2 = 8/3 u**3
        VADDPD(EPS2,J1,J1); //J1 = -1/14 + 8/3 u**3

        VBROADCASTSD(spline_r0_1,EPS2);//eps2 = 1/2
        VMULPD(EPS2,UN,EPS2);//eps2 = 1/2 u**3
        VADDPD(EPS2,J2,J2); //J2 = 4/3 - 6/5 u**2 + 1/2 u**3
        

        VMULPD(UN,X2,UN);//UN = u**4
        VBROADCASTSD(spline_r1_2,EPS2);//EPS2 = -3
        VMULPD(UN,EPS2,EPS2); //EPS2 = -3 u**4
        VADDPD(EPS2,J1,J1);// J1 = -1/14 + 8/3 u**3 - 3 u**4
        
        VMULPD(UN,X2,UN); //UN = u**5
        VBROADCASTSD(spline_r1_3,EPS2);//EPS2 = 6/5
        VMULPD(UN,EPS2,EPS2); //EPS2 = 6/5 u**5
        VADDPD(EPS2,J1,J1); //J1 = -1/14 + 8/3 u**3 - 3 u**4 +  6/5 u**5
        
        VMULPD(UN,X2,UN); //UN = u**5
        VBROADCASTSD(spline_r1_4,EPS2);//EPS2 = 1/6
        VMULPD(UN,EPS2,EPS2); //EPS2 = 1/6 u**6
        VADDPD(EPS2,J1,J1); //J1 = -1/14 + 8/3 u**3 - 3 u**4 +  6/5 u**5 - 1/6 u**6

        VBROADCASTSD(d_one,EPS2); //EPS2 = 1.0;
        VCMPPD(EPS2, X2,UN, 0x11); //UN = (u < 1.0)

        //choose whether to use J1 or J2 as the softening factor;
        VBLENDVPD(UN,J2,J1,J1);

        VBROADCASTSD(d_two, J2);//J2 = 2.0

        VCMPPD(J2,X2,UN,0x11); //UN = (u < 2.0)

        VBLENDVPD(UN, J1, EPS2, J1); //J1 = UN ? J1 : 1.0;

        //J1 is now correctly set to be the softening factor
        
        VSUBPD(J2,J2,J2); //J2 = 0;

        VCMPPD(EPS2,X2,UN,0x11); //UN = (u <1.0)

        VLOADPD(*ipdata->eps2, EPS2);
        VBLENDVPD(UN,EPS2,Y2,J2); //J2 = (Y2 < 1.0) ? Y2 : eps;



        VDIVPD(J2,J1,J1); //J1 = spf/r
        VDIVPD(J2,J1,J1); //J1 = spf/r**2
        VDIVPD(J2,J1,Y2); //J1 = spf/r**3

        VMULPD(MJ, Y2, Y2);   // Y2 = m/R^3

        VMULPD(Y2, DX, DX);   // DX = DX m / R^3 = dx m / R^3
        VMULPD(Y2, DY, DY);
        VMULPD(Y2, DZ, DZ);

        VBROADCASTSD(jpdata->x, X2);
        VBROADCASTSD(jpdata->y, Y2);
        VBROADCASTSD(jpdata->z, J2);
        VBROADCASTSD(jpdata->m, MJ);
        jpdata++;

        VADDPD(DX, AX, AX);
        VADDPD(DY, AY, AY);
        VADDPD(DZ, AZ, AZ);
    }

    VSTORPD(AX, *accdata->x);
    assert(isfinite(*accdata->x));
    VSTORPD(AY, *accdata->y);
    assert(isfinite(*accdata->y));
    VSTORPD(AZ, *accdata->z);
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
