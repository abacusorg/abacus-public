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
                           double3 &delta, double eps2, 
                           ThreeVector<double> *pacc) {
    
    double conpot = sqrt(2.0);


    double conacc = conpot*conpot*conpot;
    
    //double eps2 = eps*eps;

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
            ipdata->eps2[j] = eps2;
        }
        for(int j = e; j < 4; j++)
            ipdata->eps2[j] = eps2;

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

// do one xj on four xi's per loop
void AVXDirectDouble::KernelAccPot(ipstruct<double,4> *ipdata, jpstruct<double> *jpdata, int nsrc, 
                                     ipstruct<double,4> *deltas, double eps2, 
                                     apstruct<double,4> *accdata) {
    int j;

    double threehalf ALIGN64;
    threehalf = 1.5;

    PREFETCH(*jpdata);

    VZEROALL; // to zero out acc registers for accumulation

    VBROADCASTSD(threehalf, J1);

    VLOADPD(*ipdata->x, XI);    // load 4 xi's 
    VLOADPD(*ipdata->y, YI);
    VLOADPD(*ipdata->z, ZI);
    VLOADPD(*ipdata->eps2, EPS2);

    VADDPD_M(*deltas->x, XI, XI); // add deltas to sink position
    VADDPD_M(*deltas->y, YI, YI); 
    VADDPD_M(*deltas->z, ZI, ZI); 

    VBROADCASTSD(jpdata->x, X2); // load 1 xj into four copies
    VBROADCASTSD(jpdata->y, Y2);
    VBROADCASTSD(jpdata->z, J2);
    VBROADCASTSD(jpdata->m, MJ);
    jpdata++;

    for(j=0; j<nsrc; j++) {

        PREFETCH(*jpdata);  // does this do anything -- unlikely!

        VSUBPD(XI, X2, DX);
        VSUBPD(ZI, J2, DZ);
        VSUBPD(YI, Y2, DY);

        VMULPD(DX, DX, X2);      // X2 = DX^2
        VMULPD(DZ, DZ, J2);
        VMULPD(DY, DY, Y2);

        VADDPD(X2, J2, J2);   // J2 = X2 + J2 = DX^2 + DZ^2
        VADDPD(EPS2, Y2, Y2); // Y2 = Y2 + EPS2 = DY^2 + eps^2
        VADDPD(J2, Y2, Y2);   // Y2 = Y2 + J2 = DX^2 + DY^2 + DZ^2 + eps^2 = R^2

        VADDPD(Y2, Y2, X2);
        VCVTPD2PS(X2, xX2);   // convert to float to use rsqrt

        VRSQRTPS(xX2, xX2);   // 1/sqrt(2*R^2)
        VCVTPS2PD(xX2, X2);

        VMULPD(X2, X2, J2);   // J2 = x0^2
        VMULPD(J2, Y2, J2);   // J2 = x0^2 R^2
        VSUBPD(J2, J1, J2);   // J2 = 1.5 - x0^2 R^2
        VMULPD(J2, X2, X2);   // X2 = 1.5 x0 - x0^3 R^2

        VMULPD(X2, X2, J2);   // J2 = x0^2
        VMULPD(J2, Y2, J2);   // J2 = x0^2 R^2
        VSUBPD(J2, J1, J2);   // J2 = 1.5 - x0^2 R^2
        VMULPD(J2, X2, X2);   // X2 = 1.5 x0 - x0^3 R^2

        VMULPD(X2, MJ, MJ);   // MJ = MJ/R = m/R
        //        VSUBPD(MJ, PHI, PHI); // PHI = m/R - PHI

        VMULPD(X2, X2, Y2);   // Y2 = 1/R^2
        VMULPD(MJ, Y2, Y2);   // Y2 = m/R * 1/R^2 = m/R^3

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
    VSTORPD(AY, *accdata->y);
    VSTORPD(AZ, *accdata->z);
    //    VSTORPD(PHI, *accdata->pot);

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
