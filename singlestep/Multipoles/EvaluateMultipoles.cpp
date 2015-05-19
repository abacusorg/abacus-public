#ifdef AVXMULTIPOLES 
void (*CMptr[24])( d4 *ip1x, d4 *ip2x, d4 *ip1y, d4 *ip2y, d4 *ip1z, d4 *ip2z, 
                   d4 *cx, d4 *cy, d4 *cz, d4 *globalM, 
                   d4 *mass1, d4 *mass2) = {
     MultipoleKernel1,  MultipoleKernel2,  MultipoleKernel3,  MultipoleKernel4,
     MultipoleKernel5,  MultipoleKernel6,  MultipoleKernel7,  MultipoleKernel8,
     MultipoleKernel9,  MultipoleKernel10, MultipoleKernel11, MultipoleKernel12,
     MultipoleKernel13, MultipoleKernel14, MultipoleKernel15, MultipoleKernel16
};
#endif

Multipoles::Multipoles(int order) : basemultipoles(order) {
#ifdef AVXMULTIPOLES
for(int g=0;g<omp_get_num_procs();g++) {
    int rv;
    rv = posix_memalign( (void **) &(ip1x[g]), 256, 512 ); assert(rv==0);
    rv = posix_memalign( (void **) &(ip2x[g]), 256, 512 ); assert(rv==0);
    rv = posix_memalign( (void **) &(ip1y[g]), 256, 512 ); assert(rv==0);
    rv = posix_memalign( (void **) &(ip2y[g]), 256, 512 ); assert(rv==0);
    rv = posix_memalign( (void **) &(ip1z[g]), 256, 512 ); assert(rv==0);
    rv = posix_memalign( (void **) &(ip2z[g]), 256, 512 ); assert(rv==0);

    rv = posix_memalign( (void **) &(cx[g]),   256, 512 ); assert(rv==0);
    rv = posix_memalign( (void **) &(cy[g]),   256, 512 ); assert(rv==0);
    rv = posix_memalign( (void **) &(cz[g]),   256, 512 ); assert(rv==0);

    rv = posix_memalign( (void **) &(masses1[g]), 256, 512  ); assert(rv==0);
    rv = posix_memalign( (void **) &(masses2[g]), 256, 512  ); assert(rv==0);
    rv = posix_memalign( (void **) &(globalM[g]), 256, 32768); assert(rv==0);
}
#endif
}

void Multipoles::AnalyticCartesianMultipoles(double3 *p, int n, double3 center, 
                                             double *cm) {

    for(int m=0;m<completemultipolelength;m++) cm[m] = 0;

    for(int q=0;q<n;q++) {
        double3 r =  p[q] - center;

        double fi,fij,fijk;
        fi = 1.0;

        for(int i=0;i<=order;i++) {
            fij = fi;
            for(int j=0;j<=order-i;j++) {
                fijk = fij;
                for(int k=0;k<=order-i-j;k++) {
                    cm[ cmap(i,j,k) ] += fijk; fijk *= r.z;
                }
                fij *= r.y;
            }
            fi *= r.x;
        }
    }
}

void Multipoles::AnalyticCartesianMultipoles(float3 *p, int n, float3 center, 
                                             double *cm) {

    for(int m=0;m<completemultipolelength;m++) cm[m] = 0;

    for(int q=0;q<n;q++) {
        double3 r =  p[q] - center;

        double fi,fij,fijk;
        fi = 1.0;

        for(int i=0;i<=order;i++) {
            fij = fi;
            for(int j=0;j<=order-i;j++) {
                fijk = fij;
                for(int k=0;k<=order-i-j;k++) {
                    cm[ cmap(i,j,k) ] += fijk; fijk *= r.z;
                }
                fij *= r.y;
            }
            fi *= r.x;
        }
    }
}

void Multipoles::ASMCartesianMultipoles(double3 *xyz, int n, double3 center, 
                                        double *CM) {


#ifdef AVXMULTIPOLES

    int g = omp_get_thread_num();

    for(int k=0;k<completemultipolelength;k++) CM[k] = 0;

    for(int k=0;k<completemultipolelength;k++) 
        for(int j=0;j<4;j++) globalM[g][k].v[j] = 0;

    for(int j=0;j<4;j++) cx[g][0].v[j] = center.x;
    for(int j=0;j<4;j++) cy[g][0].v[j] = center.y;
    for(int j=0;j<4;j++) cz[g][0].v[j] = center.z;

    int end = n-(n%8);

    for(int i=end;i<n;i++) {    
        double cp[completemultipolelength];
        double3 p = xyz[i];
        AnalyticCartesianMultipoles( &p, 1, center,  &(cp[0]) );
        for(int k=0;k<completemultipolelength;k++) CM[k] += cp[k];
        // do these with analytic and add to CM
    }

    for(int j=0;j<4;j++) masses1[g][0].v[j]  = 1;
    for(int j=0;j<4;j++) masses2[g][0].v[j]  = 1;

    for(int k=0;k<end;k+=8) {

        for(int j=0;j<4;j++) ip1x[g][0].v[j] = xyz[k+j].x; 
        for(int j=0;j<4;j++) ip2x[g][0].v[j] = xyz[k+4+j].x;

        for(int j=0;j<4;j++) ip1y[g][0].v[j] = xyz[k+j].y;
        for(int j=0;j<4;j++) ip2y[g][0].v[j] = xyz[k+4+j].y;

        for(int j=0;j<4;j++) ip1z[g][0].v[j] = xyz[k+j].z;
        for(int j=0;j<4;j++) ip2z[g][0].v[j] = xyz[k+4+j].z;

        (CMptr[order-1])( &(ip1x[g][0]), &(ip2x[g][0]), 
                        &(ip1y[g][0]), &(ip2y[g][0]),
                        &(ip1z[g][0]), &(ip2z[g][0]),
                        &(cx[g][0]), &(cy[g][0]), &(cz[g][0]), 
                        &(globalM[g][0]),  &(masses1[g][0]), &(masses2[g][0]) );

    }

    for(int k=0;k<completemultipolelength;k++) 
       for(int j=0;j<4;j++) CM[k] += globalM[g][k].v[j];
#else
     AnalyticCartesianMultipoles(xyz, n, center, CM);
#endif

}

void Multipoles::ASMCartesianMultipoles(float3 *xyz, int n, float3 center, 
                                        double *CM) {


#ifdef AVXMULTIPOLES

    int g = omp_get_thread_num();
    //assert(g < 16);

    for(int k=0;k<completemultipolelength;k++) CM[k] = 0;

    for(int k=0;k<completemultipolelength;k++) 
        for(int j=0;j<4;j++) globalM[g][k].v[j] = 0;

    for(int j=0;j<4;j++) cx[g][0].v[j] = center.x;
    for(int j=0;j<4;j++) cy[g][0].v[j] = center.y;
    for(int j=0;j<4;j++) cz[g][0].v[j] = center.z;

    int end = n-(n%8);

    for(int i=end;i<n;i++) {    
        double cp[completemultipolelength];
        float3 p = xyz[i];
        AnalyticCartesianMultipoles( &p, 1, center,  &(cp[0]) );
        for(int k=0;k<completemultipolelength;k++) CM[k] += cp[k];
        // do these with analytic and add to CM
    }

    for(int j=0;j<4;j++) masses1[g][0].v[j]  = 1;
    for(int j=0;j<4;j++) masses2[g][0].v[j]  = 1;

    for(int k=0;k<end;k+=8) {

        for(int j=0;j<4;j++) ip1x[g][0].v[j] = xyz[k+j].x; 
        for(int j=0;j<4;j++) ip2x[g][0].v[j] = xyz[k+4+j].x;

        for(int j=0;j<4;j++) ip1y[g][0].v[j] = xyz[k+j].y;
        for(int j=0;j<4;j++) ip2y[g][0].v[j] = xyz[k+4+j].y;

        for(int j=0;j<4;j++) ip1z[g][0].v[j] = xyz[k+j].z;
        for(int j=0;j<4;j++) ip2z[g][0].v[j] = xyz[k+4+j].z;

        (CMptr[order-1])( &(ip1x[g][0]), &(ip2x[g][0]), 
                        &(ip1y[g][0]), &(ip2y[g][0]),
                        &(ip1z[g][0]), &(ip2z[g][0]),
                        &(cx[g][0]), &(cy[g][0]), &(cz[g][0]), 
                        &(globalM[g][0]),  &(masses1[g][0]), &(masses2[g][0]) );

    }

    for(int k=0;k<completemultipolelength;k++) 
       for(int j=0;j<4;j++) CM[k] += globalM[g][k].v[j];
#else
     AnalyticCartesianMultipoles(xyz, n, center, CM);
#endif

}
