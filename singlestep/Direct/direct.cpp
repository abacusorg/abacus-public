// Copyright 2012-2025 The Abacus Developers
// SPDX-License-Identifier: GPL-3.0-or-later

#ifdef AVXDIRECT
    #ifdef DIRECTCUBICSPLINE
    #include "avxdirectdouble_spline.cc"
    #include "avxdirectfloatNR_spline.cc"
    #elif defined DIRECTCUBICPLUMMER
    #error "AVX cubic plummer softening not implemented yet"
    #elif defined DIRECTSINGLESPLINE
    #error "AVX single spline not implemented yet"
    #else
    #include "avxdirectdouble.cc"
    #include "avxdirectfloatNR.cc"
    #endif
#endif

#include "direct.h"

Direct::Direct(void) {
#ifdef AVXDIRECT
    directdouble = new AVXDirectDouble(MAXSOURCELENGTH);
    directfloatnr = new AVXDirectFloatNR(MAXSOURCELENGTH);
#endif
}

Direct::~Direct(void) {
#ifdef AVXDIRECT
    delete directdouble;
    delete directfloatnr;
#endif
}

void Direct::AVXExecute(FLOAT3 *sinks, FLOAT3 *sources, int nsinks, int nsources, FLOAT3 delta, FLOAT eps, FLOAT3 *SA) {
    #ifdef AVXDIRECT
        #ifdef DOUBLEPRECISION
            directdouble->compute(nsources, sources, nsinks, sinks, delta, eps, SA);
        #else  
            directfloatnr->compute(nsources, sources, nsinks, sinks, delta, eps, SA);
        #endif
    
    #else 
        Execute(sinks, sources, nsinks, nsources, delta, eps, SA); 
    #endif
}

#include "DirectCPUKernels.cc"

void Direct::Execute(FLOAT3 *sinks, FLOAT3 *sources, int nsinks, int nsources, FLOAT3 delta, FLOAT eps, FLOAT3 *SA) {
    for(int q=0;q<nsinks;q++) {
        FLOAT3 acc = FLOAT3(0);
        FLOAT3 sink = delta+sinks[q];
        //#pragma simd assert
        for(int p=0;p<nsources;p++) {
            directkernel<0>(sink,sources[p],acc,eps);
        }     
        SA[q] += acc;
    }
}
