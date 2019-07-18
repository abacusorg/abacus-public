class InCoreConvolution : public basemultipoles {
private:
    Complex *_mcache,*_tcache;
    double *_dcache,*mfactor;
    int cpd,blocksize,nblocks,cpdhalf,CompressedMultipoleLengthXY;
    const int thread_padding = 1024;
    int bufsize, padblocksize;

public:
    InCoreConvolution(int order, int cpd, int blocksize) : 
            basemultipoles(order), cpd(cpd), blocksize(blocksize) {

        CompressedMultipoleLengthXY  = ((1+cpd)*(3+cpd))/8;
        nblocks = (cpd*cpd)/blocksize;
        cpdhalf = (cpd-1)/2;

        /* We're going to allocate a bunch of work arrays, three per thread.
        We've been given a blocksize, and at present it must divide CPD^2
        evenly, which means that it must be odd.  
        However, we'd like to pad each multipole out to a multiple of 4,
        so that even doubles are 32-byte aligned.
        */
        padblocksize = (blocksize/4+1)*4;

        bufsize = completemultipolelength*padblocksize;
        bufsize = 512*(bufsize/512+1);   
            // Round up to a multiple of 512 to avoid thread contention


        // TODO: Actually, we only need tcache for one blocksize
        // per core, since we compute each multipole and then save it.
        // That could increase the blocksize.  But initial testing
        // showed this slowed the code down, so DJE didn't implement.

        int cs = omp_get_max_threads() * bufsize;
        int ret;
        ret = posix_memalign((void **)&_mcache, 4096, cs*sizeof(Complex));
        assert(ret==0);
        ret = posix_memalign((void **)&_tcache, 4096, cs*sizeof(Complex));
        assert(ret==0);
        ret = posix_memalign((void **)&_dcache, 4096, cs*sizeof(double));
        // ret = posix_memalign((void **)&_dcache, 4096, cs*sizeof(Complex));
        assert(ret==0);

        /// _mcache = new Complex[cs]; 
        /// _dcache = new  double[cs]; 
        /// _tcache = new Complex[cs];
        mfactor = new double[completemultipolelength];
        int a,b,c;
        FORALL_COMPLETE_MULTIPOLES_BOUND(a,b,c,order) {
            mfactor[ cmap(a,b,c) ] = 
                pow(-1.0,a+b+c)/(fact2(2*(a+b+c)-1)*fact(a)*fact(b)*fact(c));
        }
    }
    ~InCoreConvolution(void) {
        /// delete[] _mcache; delete[] _dcache; delete[] _tcache; 
        free(_mcache); free(_tcache); free(_dcache);
        delete[] mfactor;
    }
    void InCoreConvolve(Complex *FFTM, DFLOAT *CompressedDerivatives);
    unsigned long long int  ConvolutionArithmeticCount(void);
};

unsigned long long int InCoreConvolution::ConvolutionArithmeticCount(void) {
    int a,b,c,i,j,k;

    unsigned long long int l = 0;

    /// This appears to be a bug
    unsigned long long int cpd3 = ((unsigned long long int ) cpd)*cpd*cpd;
    /// unsigned long long int cpd3 = ((unsigned long long int ) cpd)*cpd*(1+cpd)/2;

    FORALL_REDUCED_MULTIPOLES_BOUND(a,b,c,order) {

            FORALL_COMPLETE_MULTIPOLES_BOUND(i,j,k,order-a-b-c) {
                if( ((a+i)+(b+j)+(c+k))%2 ==0 ) {
                    l += cpd3 * 4;
                }
            }

            FORALL_COMPLETE_MULTIPOLES_BOUND(i,j,k,order-a-b-c) {
                if( ((a+i)+(b+j)+(c+k))%2 == 1 ) {
                    l += cpd3* 4;
                }
            }

            l += cpd3*4; // complex multiply
    }

    return l;

}

inline void MVMpair(Complex *t, Complex *m, double *d, int n) { 
    // We require here that d[] is already broadcast into a list of pairs
    // of doubles, so that we can do the FMA math as a pure vector of doubles.
    double *td = (double *)t;
    double *md = (double *)m;
    n*=2;
    #pragma simd assert
    for(int i=0;i<n;i++) td[i] += md[i] * d[i]; 
}

inline void MVM(Complex *t, Complex *m, double *d, int n) { 
    //#pragma simd
    // __builtin_assume_aligned
    //#pragma vector aligned
    // for(int i=0;i<n;i++) t[i] += m[i] * d[i]; 

    // An SSE equivalent
    for (int i=0; i<n; i++) {
        __m128d xx = _mm_load1_pd((double *)(d+i));
        __m128d yy = _mm_load_pd((double *)(m+i));
        __m128d zz = _mm_load_pd((double *)(t+i));
        __m128d xy = _mm_mul_pd(xx,yy);
        zz = _mm_add_pd(zz,xy);
        _mm_store_pd((double *)(t+i),zz);
    }
}

void InCoreConvolution::InCoreConvolve(Complex *FFTM, DFLOAT *CompressedD) {
    // why is this 50% faster with dynamic?
    #pragma omp parallel for schedule(dynamic)
    for(int block=0;block<nblocks;block++) {
            
        Complex *FM, *FMabc, *FMap2bcm2, *FMabp2cm2;
        double  *FD, *FDabc, *FDap2bcm2, *FDabp2cm2, *FDd;
        Complex *mcache, *tcache;
        double  *dcache;
        int xyz,a,b,c,i,j,k,xyzbegin;

        int gindex = omp_get_thread_num()*bufsize;
        dcache = &(_dcache[gindex]); 
        mcache = &(_mcache[gindex]);
        tcache = &(_tcache[gindex]);

        xyzbegin = block*blocksize;
        // Load the reduced multipoles into the cache
        FORALL_REDUCED_MULTIPOLES_BOUND(a,b,c,order) {
            FM = &(FFTM[rmap(a,b,c) * cpd*cpd + xyzbegin]);
            int m = cmap(a,b,c)*padblocksize;
            FOR(xyz,0,blocksize-1) mcache[m + xyz] = FM[xyz];
        }

        // Apply the trace-free recursion to generate complete multipoles
        TRACEFREERECURSION(a,b,c,order) {
            FMabc       = &(mcache[cmap(a  ,b  ,c  ) * padblocksize]);
            FMap2bcm2   = &(mcache[cmap(a+2,b  ,c-2) * padblocksize]);
            FMabp2cm2   = &(mcache[cmap(a  ,b+2,c-2) * padblocksize]);
            FOR(xyz,0,blocksize-1) 
                FMabc[xyz] = - FMap2bcm2[xyz] - FMabp2cm2[xyz];
        }

        // Apply the scaling to the complete multipoles
        FORALL_COMPLETE_MULTIPOLES_BOUND(a,b,c,order)  {
            FM = &(mcache[cmap(a,b,c) * padblocksize]);
            double mf = mfactor[ cmap(a,b,c) ];
            FOR(xyz,0,blocksize-1) FM[xyz] *= mf;
        }

        // Load the derivatives.  This requires some care.
        // The derivatives were only stored in one half-quadrant.
        // Swapping (x,y) and (i,j) makes no change; that completes the quadrant.
        // Then the other quadrants are related by parity.
        
        FORALL_REDUCED_MULTIPOLES_BOUND(a,b,c,order) {

            int m = cmap(a,b,c)*padblocksize;
            // Rather than unpacking the (x,y) by mod, we'll do 
            // one and then increment
            int y = xyzbegin%cpd;
            int x = (xyzbegin-y)/cpd;
            FOR(xyz,0,blocksize-1) {
                /// int l = xyzbegin + xyz;
                /// int y = l%cpd;
                /// int x = (l-y)/cpd;
                int xp = x; if(xp>cpdhalf) xp = cpd - x;
                int yp = y; if(yp>cpdhalf) yp = cpd - y;

                double xR;
                if(yp>=xp)  {
                    int i = rmap(b,a,c)*CompressedMultipoleLengthXY;
                    xR = (double) CompressedD[ i + (yp*(yp+1))/2 + xp ];
                }
                else {
                    int i = rmap(a,b,c)*CompressedMultipoleLengthXY;
                    xR = (double) CompressedD[ i + (xp*(xp+1))/2 + yp ];
                }

                if( (yp!=y) && ((b&0x1)==1) ) xR *= -1;
                if( (xp!=x) && ((a&0x1)==1) ) xR *= -1;
                dcache[m + xyz] = xR;
                y++;
                if (y==cpd) { y = 0; x++; }  // Accomplish the mod wrap
            }
        }
        /* TODO: In the above code, the completion of the first quadrant
        results in a poor access pattern.  Is this really worth the speed
        penalty to halve the size of the derivative file?  It would be 
        1/8 of the Multipoles even if we didn't.  And even if we do
        want that small size, is it wise to do the half-quadrant 
        completion only on demand?  Perhaps a block-tranpose would be faster.
        At present, we're filling half the plane with badly-ordered
        accesses, whereas we could get away with only a half-quadrant.
        */

        // Recurse to restore the complete derivatives
        TRACEFREERECURSION(a,b,c,order) {
            FDabc       = &(dcache[cmap(a  ,b  ,c  ) * padblocksize]);
            FDap2bcm2   = &(dcache[cmap(a+2,b  ,c-2) * padblocksize]);
            FDabp2cm2   = &(dcache[cmap(a  ,b+2,c-2) * padblocksize]);
            FOR(xyz,0,blocksize-1) 
                FDabc[xyz] = - FDap2bcm2[xyz] - FDabp2cm2[xyz];
        }

/*
        // Now we want to broadcast each double into two, so that
        // the MVM math is trivial double-FMA.  Until now, dcache
        // has only used the first half of its buffer, length
        // completemultipolelength*blocksize.  Now we want to double
        // that.
        FDd = dcache+completemultipolelength*blocksize*2-2;
        for (xyz=completemultipolelength*blocksize-1; xyz>=0; xyz--, FDd-=2) {
            FDd[1] = FDd[0] = dcache[xyz];
        }
*/

        // Now compute the Taylors from M & D
        // We do this in two parts, since the D's are stored as real,
        // but are either pure real or pure imaginary, according
        // to their i+j+k parity.
        FORALL_REDUCED_MULTIPOLES_BOUND(a,b,c,order) {

            FOR(xyz,0,blocksize-1) tcache[xyz] = 0;
            FORALL_COMPLETE_MULTIPOLES_BOUND(i,j,k,order-a-b-c) {
                // Odd parity only
                if( (((a+i)+(b+j)+(c+k))&0x1) == 1 ) {
                    FD  = &(dcache[cmap(a+i,b+j,c+k) * padblocksize]);
                    FM  = &(mcache[cmap(i  ,j  ,k  ) * padblocksize]);
                    MVM(tcache, FM, FD, blocksize);
                }
            }
            // Now multiply by i, before we accumulate the rest
            double *tcache_dbl = (double *)tcache;
            for (xyz=0; xyz<2*blocksize; xyz+=2) {
                double tmp = tcache_dbl[xyz];
                tcache_dbl[xyz] = -tcache_dbl[xyz+1];
                tcache_dbl[xyz+1] = tmp;
            }

            FORALL_COMPLETE_MULTIPOLES_BOUND(i,j,k,order-a-b-c) {
                // Even parity only 
                if( (((a+i)+(b+j)+(c+k))&0x1) ==0 ) {
                    FD  = &(dcache[cmap(a+i,b+j,c+k) * padblocksize]);
                    FM  = &(mcache[cmap(i  ,j  ,k  ) * padblocksize]);
                    MVM(tcache,FM, FD, blocksize); 
                }
            }
            Complex *FT = &(FFTM[rmap(a,b,c) * cpd*cpd + xyzbegin]);
            FOR(xyz,0,blocksize-1) FT[xyz] = tcache[xyz];
        }
    }
}
