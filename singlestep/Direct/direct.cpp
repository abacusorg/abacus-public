#ifdef AVXDIRECT
    #ifdef DIRECTSPLINE
    #include "avxdirectdouble_spline.cc"
    #include "avxdirectfloatNR_spline.cc"
    #else
    #include "avxdirectdouble.cc"
    #include "avxdirectfloatNR.cc"
    #endif
#endif


Direct::Direct(void) {
#ifdef AVXDIRECT
    directdouble = new AVXDirectDouble(MAXSOURCELENGTH);
    directfloatnr = new AVXDirectFloatNR(MAXSOURCELENGTH);
#endif
}

Direct::~Direct(void) { }

void Direct::AVXExecute(FLOAT3 *sinks, FLOAT3 *sources, int nsinks, int nsources, FLOAT3 delta, FLOAT eps2, FLOAT3 *SA) {
    #ifdef AVXDIRECT
        #ifdef DOUBLEPRECISION
            directdouble->compute(nsources, sources, nsinks, sinks, delta, eps2, SA);
        #else  
            directfloatnr->compute(nsources, sources, nsinks, sinks, delta, eps2, SA);
        #endif
    
    #else 
        Execute(sinks, sources, nsinks, nsources, delta, eps2, SA); 
    #endif
}

#include "DirectCPUKernels.cc"

void Direct::Execute(FLOAT3 *sinks, FLOAT3 *sources, int nsinks, int nsources, FLOAT3 delta, FLOAT eps2, FLOAT3 *SA) {
    for(int q=0;q<nsinks;q++) {
        FLOAT3 acc = FLOAT3(0);
        FLOAT3 sink = delta+sinks[q];
        #pragma simd
        for(int p=0;p<nsources;p++) {
            directkernel<0>(sink,sources[p],acc,eps2);
        }     
        SA[q] += acc;
    }
}
