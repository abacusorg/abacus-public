#ifdef AVXDIRECT
    #include "avxdirectdouble.cc"
    #include "avxdirectfloatNR.cc"
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
        #ifdef DIRECTSPLINE
            #error "Spline softening is not implemented in AVX yet."
        #endif
    
        #ifdef DOUBLEPRECISION
            directdouble->compute(nsources, sources, nsinks, sinks, delta, eps2, SA);
        #else  
            directfloatnr->compute(nsources, sources, nsinks, sinks, delta, eps2, SA);
        #endif
    
    #else 
        Execute(sinks, sources, nsinks, nsources, delta, eps2, SA); 
    #endif
}

void Direct::Execute(FLOAT3 *sinks, FLOAT3 *sources, int nsinks, int nsources, FLOAT3 delta, FLOAT eps2, FLOAT3 *SA) {
    for(int q=0;q<nsinks;q++) {
        FLOAT3 acc = FLOAT3(0);
        FLOAT3 x = delta+sinks[q];
        for(int p=0;p<nsources;p++) {
            FLOAT3 rr = sources[p] - x;
            double ir = 1.0/sqrt(rr.norm2() + eps2);
            double ir3 = ir*ir*ir;
            acc -= ir3*rr;
        }
        SA[q] += acc;
    }
}
