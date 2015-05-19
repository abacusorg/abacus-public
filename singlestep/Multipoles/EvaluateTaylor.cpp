#ifdef AVXMULTIPOLES
void (*Tptr[17])(d4 * px, d4 * py, d4 * pz, d4 * cx, d4 * cy, d4 * cz, 
                 d4 * tx, d4 * ty, d4 * tz, d4 * ax, d4 * ay, d4 * az) = { 
 TaylorKernel0, TaylorKernel1, TaylorKernel2, TaylorKernel3, TaylorKernel4, 
 TaylorKernel5, TaylorKernel6, TaylorKernel7, TaylorKernel8, 
 TaylorKernel9, TaylorKernel10, TaylorKernel11, TaylorKernel12, 
 TaylorKernel13, TaylorKernel14, TaylorKernel15, TaylorKernel16 }; 
#endif


Taylor::~Taylor(void) { }

Taylor::Taylor(int order) : basemultipoles(order) {
#ifdef AVXMULTIPOLES
    for(int g=0;g<omp_get_max_threads();g++) {
        int rv;
        rv = posix_memalign( (void **) &(cx[g]), 256, 512 ); assert(rv==0);
        rv = posix_memalign( (void **) &(cy[g]), 256, 512 ); assert(rv==0);
        rv = posix_memalign( (void **) &(cz[g]), 256, 512 ); assert(rv==0);

        rv = posix_memalign( (void **) &(ax[g]), 256, 512 ); assert(rv==0);
        rv = posix_memalign( (void **) &(ay[g]), 256, 512 ); assert(rv==0);
        rv = posix_memalign( (void **) &(az[g]), 256, 512 ); assert(rv==0);

        rv = posix_memalign( (void **) &(px[g]), 256, 512 ); assert(rv==0);
        rv = posix_memalign( (void **) &(py[g]), 256, 512 ); assert(rv==0);
        rv = posix_memalign( (void **) &(pz[g]), 256, 512 ); assert(rv==0);

        rv = posix_memalign( (void **) &(Qx[g]), 256, 32768 ); assert(rv==0);
        rv = posix_memalign( (void **) &(Qy[g]), 256, 32768 ); assert(rv==0);
        rv = posix_memalign( (void **) &(Qz[g]), 256, 32768 ); assert(rv==0);
    }
#endif
}


void Taylor::AnalyticEvaluateTaylor(double *CT, double3 expansioncenter, int np,
                                    double3 *ps, double3 *acc) {
    assert(np>=0);

    for(int p=0;p<np;p++) {

        int i,j,k;
        double fi,fij,fijk;
        double3 delta = ps[p] - expansioncenter;
        double3 a = double3(0);

        fi = 1.0;
        FOR(i,0,order-1) {
            fij = fi;
            FOR(j,0,order-1-i) {
                fijk = fij;
                FOR(k,0,order-1-i-j) {
                    a.x -= (i+1) * CT[ cmap(i+1,j  ,k  ) ] * fijk;
                    a.y -= (j+1) * CT[ cmap(i  ,j+1,k  ) ] * fijk;
                    a.z -= (k+1) * CT[ cmap(i  ,j  ,k+1) ] * fijk;
                    fijk *= delta.z;
                }
                fij *= delta.y;
            }
            fi *= delta.x;
        }
        acc[p] += a;
    }
}

void Taylor::AnalyticEvaluateTaylor(double *CT, float3 expansioncenter, int np,
                                    float3 *ps, float3 *acc) {
    assert(np>=0);

    for(int p=0;p<np;p++) {

        int i,j,k;
        double fi,fij,fijk;
        double3 delta = ps[p] - expansioncenter;
        double3 a = double3(0);

        fi = 1.0;
        FOR(i,0,order-1) {
            fij = fi;
            FOR(j,0,order-1-i) {
                fijk = fij;
                FOR(k,0,order-1-i-j) {
                    a.x -= (i+1) * CT[ cmap(i+1,j  ,k  ) ] * fijk;
                    a.y -= (j+1) * CT[ cmap(i  ,j+1,k  ) ] * fijk;
                    a.z -= (k+1) * CT[ cmap(i  ,j  ,k+1) ] * fijk;
                    fijk *= delta.z;
                }
                fij *= delta.y;
            }
            fi *= delta.x;
        }
        acc[p] += a;
    }
}

void Taylor::ASMEvaluateTaylor( double *CT, double3 center, int n, double3 *xyz,
                                double3 *acc) {

#ifdef AVXMULTIPOLES

    int g = omp_get_thread_num();

    for(int j=0;j<4;j++) {  
        cx[g][0].v[j] = center.x;
        cy[g][0].v[j] = center.y;
        cz[g][0].v[j] = center.z;
    }

    int i,a,b,c;
    i = 0;
    FOR(a,0,order-1)
        FOR(b,0,order-1-a)
            FOR(c,0,order-1-a-b) {
                for(int j=0;j<4;j++) {
                    Qx[g][i].v[j] = (a+1)*CT[ cmap(a+1,b  ,c  ) ];
                    Qy[g][i].v[j] = (b+1)*CT[ cmap(a  ,b+1,c  ) ];
                    Qz[g][i].v[j] = (c+1)*CT[ cmap(a  ,b  ,c+1) ];
                }
                i++;
            }

    int l = 0;
    int m = n-1;

    for(int k=l;k<=m;k+=4) {
        int end=k+3; if(k+3>m) end = m;
        
        for(int j=k;j<=end;j++) {
            px[g][0].v[j-k] = xyz[j].x;
            py[g][0].v[j-k] = xyz[j].y;
            pz[g][0].v[j-k] = xyz[j].z;
        }

        for(int j=0;j<4;j++) ax[g][0].v[j] = 1;

        (Tptr[order])( &(px[g][0]),&(py[g][0]),&(pz[g][0]),
                       &(cx[g][0]),&(cy[g][0]),&(cz[g][0]),
                       &(Qx[g][0]),&(Qy[g][0]),&(Qz[g][0]),
            &(ax[g][0]),&(ay[g][0]),&(az[g][0]) );

        for(int j=k;j<=end;j++) {
            acc[j].x -= ax[g][0].v[j-k];
            acc[j].y -= ay[g][0].v[j-k];
            acc[j].z -= az[g][0].v[j-k];
        }
    }

#else 

    AnalyticEvaluateTaylor(CT, center, n, xyz, acc);

#endif

}

void Taylor::ASMEvaluateTaylor( double *CT, float3 center, int n, float3 *xyz, 
                                float3 *acc) {

#ifdef AVXMULTIPOLES

    int g = omp_get_thread_num();

    for(int j=0;j<4;j++) {  
        cx[g][0].v[j] = center.x;
        cy[g][0].v[j] = center.y;
        cz[g][0].v[j] = center.z;
    }

    int i,a,b,c;
    i = 0;
    FOR(a,0,order-1)
        FOR(b,0,order-1-a)
            FOR(c,0,order-1-a-b) {
                for(int j=0;j<4;j++) {
                    Qx[g][i].v[j] = (a+1)*CT[ cmap(a+1,b  ,c  ) ];
                    Qy[g][i].v[j] = (b+1)*CT[ cmap(a  ,b+1,c  ) ];
                    Qz[g][i].v[j] = (c+1)*CT[ cmap(a  ,b  ,c+1) ];
                }
                i++;
            }

    int l = 0;
    int m = n-1;

    for(int k=l;k<=m;k+=4) {
        int end=k+3; if(k+3>m) end = m;
        
        for(int j=k;j<=end;j++) {
            px[g][0].v[j-k] = xyz[j].x;
            py[g][0].v[j-k] = xyz[j].y;
            pz[g][0].v[j-k] = xyz[j].z;
        }

        for(int j=0;j<4;j++) ax[g][0].v[j] = 1;

        (Tptr[order])( &(px[g][0]),&(py[g][0]),&(pz[g][0]),
                       &(cx[g][0]),&(cy[g][0]),&(cz[g][0]),
                       &(Qx[g][0]),&(Qy[g][0]),&(Qz[g][0]),
            &(ax[g][0]),&(ay[g][0]),&(az[g][0]) );

        for(int j=k;j<=end;j++) {
            acc[j].x -= ax[g][0].v[j-k];
            acc[j].y -= ay[g][0].v[j-k];
            acc[j].z -= az[g][0].v[j-k];
        }
    }

#else 

    AnalyticEvaluateTaylor(CT, center, n, xyz, acc);

#endif

}

