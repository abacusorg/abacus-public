// Copyright 2012-2025 The Abacus Developers
// SPDX-License-Identifier: GPL-3.0-or-later

class InCoreConvolution : public basemultipoles {
private:
    Complex *_mcache,*_tcache;
    double *_dcache,*mfactor;
    int cpd,cpdky_pad,blocksize,nblocks,cpdhalf,CompressedMultipoleLengthXY;
    const int thread_padding = 1024;
    int bufsize, padblocksize, blocksize_even;
    int twoD;  // Load the derivatives in the format expected by the 2D code
    int kysize;  // 2D: each z-slab has [x: 0..cpd][y: 0..kysize]

public:
    InCoreConvolution(int order, int cpd, int _kysize, int blocksize, int _twoD, int _cpdky_pad = 0) : 
            basemultipoles(order), cpd(cpd), blocksize(blocksize), twoD(_twoD), kysize(_kysize) {
				
		STDLOG(3," Entering ICC constructor\n");
		// In the serial case, a set of x-y multipoles has cpd*cpd cells.
        // In the parallel case, we pad out y to _cpdky_pad -- some number slightly bigger than cpd*node_ky_size.
		
		if (_cpdky_pad == 0) cpdky_pad = cpd * kysize; 
		else cpdky_pad = _cpdky_pad;
		
        CompressedMultipoleLengthXY  = ((1+cpd)*(3+cpd))/8;  // only 1D
        nblocks = (cpd*kysize)/blocksize;
        cpdhalf = (cpd-1)/2;

        /* We're going to allocate a bunch of work arrays, three per thread.
        We've been given a blocksize, and at present it must divide CPD^2
        evenly, which means that it must be odd.  
        However, we'd like to pad each multipole out to a multiple of 2,
        so that Complex are 32-byte aligned and can use AVX.
        */

        blocksize_even = blocksize;
        if ((blocksize&1)==1) blocksize_even += 1; 
        assert((blocksize_even&1)==0);   // n must be even
        padblocksize = blocksize_even;

        bufsize = completemultipolelength*padblocksize;
        bufsize = 512*(bufsize/512+1);   
            // Round up to a multiple of 512 to avoid thread contention


        // TODO: Actually, we only need tcache for one blocksize
        // per core, since we compute each multipole and then save it.
        // That could increase the blocksize.  But initial testing
        // showed this slowed the code down, so DJE didn't implement.

        int cs = omp_get_max_threads() * bufsize;
        int ret;
        ret = posix_memalign((void **)&_mcache, PAGE_SIZE, cs*sizeof(Complex));
        assert(ret==0);
        ret = posix_memalign((void **)&_tcache, PAGE_SIZE, cs*sizeof(Complex));
        assert(ret==0);
        ret = posix_memalign((void **)&_dcache, PAGE_SIZE, cs*sizeof(double));
        assert(ret==0);

        mfactor = new double[completemultipolelength];
				
		
        int a,b,c;
        FORALL_COMPLETE_MULTIPOLES_BOUND(a,b,c,order) {
            mfactor[ cmap(a,b,c) ] = 
                pow(-1.0,a+b+c)/(fact2(2*(a+b+c)-1)*fact(a)*fact(b)*fact(c));
        }
				
    }
    ~InCoreConvolution(void) {
        free(_mcache); free(_tcache); free(_dcache);
        delete[] mfactor;
    }
    void InCoreConvolve(Complex *FFTM, DFLOAT *CompressedDerivatives);
    unsigned long long int  ConvolutionArithmeticCount(void);
};


unsigned long long int InCoreConvolution::ConvolutionArithmeticCount(void) {
    int a,b,c,i,j,k;

    // Not counting math to unwrap the derivatives during loading.
    // Not counting math to set up the loops and pointers.
    unsigned long long int l = 0;
    unsigned long long int cpd3 = ((unsigned long long int ) cpd)*cpd*(1+cpd)/2;

    TRACEFREERECURSION(a,b,c,order) { 
        l+=4;   // Multipole tracefree recursion
        l+=2;   // Derivative tracefree recursion
    }
    FORALL_COMPLETE_MULTIPOLES_BOUND(a,b,c,order)  {
        l+=2;   // Multipole rescaling
    }

    FORALL_REDUCED_MULTIPOLES_BOUND(a,b,c,order) {
        FORALL_COMPLETE_MULTIPOLES_BOUND(i,j,k,order-a-b-c) {
            l += 4;  // FMA work
        }
        l += 2; // multiply by I
    }

    return l*cpd3;
}

// =======================================================
// Some intrinsics for Complex vector operations, as we have 
// some evidence that gcc isn't figuring these out.

/* The following routines are trying to speed up calculations on the
complex<double> vectors.  It seems that gcc isn't vectorizing effectively;
we're seeing notable improvements with simple intrinsics. */

inline void FMAvector(Complex *t, Complex *m, double *d, int n) { 
    // Multiply a complex vector by a scalar and add to a complex accumulator
    // We are doing this in sets of two cells, since that fits AVX
    // and may also suit the memory access pattern of Power9 slices.

    #ifdef HAVE_AVX
    // An AVX equivalent -- this does two at once, so we have to pad to 2!
    // Will compute off the end.  
    // We have to broadcast two consecutive doubles into (g,g,h,h)
    // First broadcast to = (g,h,g,h), then shuffle to (g,g,h,h)
    // TODO: Could FMA here, but testing showed little gain 
    #define AVX_COMPLEX_DOUBLE_FMA \
        {__m256d _xx = _mm256_broadcast_pd((__m128d *)(d)); d+=2; \
        __m256d xx = _mm256_shuffle_pd(_xx,_xx, 0xc); \
        __m256d yy = _mm256_load_pd((double *)(m)); m+=2; \
        __m256d zz = _mm256_load_pd((double *)(t)); \
        zz = _mm256_add_pd(_mm256_mul_pd(xx,yy),zz); \
        _mm256_store_pd((double *)(t), zz); t+=2;}

    n = (n>>1);   // Divide by 2; this is how many cell pairs
    for (int i=0;i<n;i++) AVX_COMPLEX_DOUBLE_FMA
    return;
    #endif

    #ifdef HAVE_VSX
    // In VSX, let's try typecasting from Complex
    vector<double> *_t = (vector<double> *)t;
    vector<double> *_m = (vector<double> *)m;
    #pragma GCC ivdep
    for (int i=0; i<n;i++) t[i] += m[i] * d[i];
    return;
    #endif

    int i=0;
    /* This is the base version */
    // #pragma omp simd aligned(t,m:16) aligned(d:8)
    // for(i=0; i<n;i++) t[i] += m[i] * d[i];

    /* // A simple unrolled version; will compute off the end in 2's 
    #pragma omp simd aligned(t,m:32) aligned(d:16)
    for(i=0; i<n;i+=2) {
        t[i] += m[i] * d[i]; 
        t[i+1] += m[i+1] * d[i+1]; 
    }
    */ 
    // An alternate, with pointer movement
    #pragma GCC ivdep
    for(i=0; i<(n>>1);i++) {
        *t += (*m)*(*d); t++; m++; d++;
        *t += (*m)*(*d); t++; m++; d++;
    }

    /*
    // An SSE equivalent, as an example for VSX
    for (int i=0; i<n; i++) {
        __m128d xx = _mm_load1_pd((double *)(d+i));
        __m128d yy = _mm_load_pd((double *)(m+i));
        __m128d zz = _mm_load_pd((double *)(t+i));
        __m128d xy = _mm_mul_pd(xx,yy);
        zz = _mm_add_pd(zz,xy);
        _mm_store_pd((double *)(t+i),zz);
    }
    */

}

inline void Multiply_by_scalar(Complex *m, double &s, int n) {
    // Multiply a complex vector by a constant scalar
    // The AVX code requires n to be divisible by 2
    int xyz;
    #ifndef HAVE_AVX
        #pragma omp simd aligned(m:16) 
        FOR(xyz,0,n-1) m[xyz] *= s;
    #else
        __m256d scalar = _mm256_broadcast_sd((double *)&s);
        n = (n>>1);
        for (xyz=0; xyz<n; xyz++) {
            __m256d mm = _mm256_load_pd((double *)(m));
            mm = _mm256_mul_pd(mm,scalar);
            _mm256_store_pd((double *)(m),mm); m+=2;
        }
    #endif
}

inline void Multiply_by_I(Complex *t, int blocksize) {
    // Multiply a complex vector by i really means swapping to (-im, re)
    int xyz;
    #ifndef HAVE_AVX
        double *tcache_dbl = (double *)t;
		#pragma omp simd aligned(t:16) 		
        for (xyz=0; xyz<2*blocksize; xyz+=2) {
            double tmp = tcache_dbl[xyz];
            tcache_dbl[xyz] = -tcache_dbl[xyz+1];
            tcache_dbl[xyz+1] = tmp;
        }
    #else
        // Consider using SSE intrinsics to do this flip
        __m128d neg = _mm_set_pd(0.0, -0.0);
        for (xyz=0; xyz<blocksize; xyz++) {
            __m128d tswap = _mm_load_pd((double *)(t+xyz));
            tswap = _mm_shuffle_pd(tswap,tswap,0x1);
            tswap = _mm_xor_pd(tswap, neg);
            _mm_store_pd((double *)(t+xyz),tswap);
        }
    #endif
}

inline void Set_Vector_to_Zero(Complex *t, int n) {
    // The AVX code requires t to be 32-byte aligned and n to be even
    #ifndef HAVE_AVX
        #pragma omp simd aligned(t:16) 
        for (int i=0; i<n; i++) t[i] = 0.0;
    #else
        __m256d zero = _mm256_set1_pd(0.0);
        n = (n>>1);
        for (int i=0; i<n; i++) { _mm256_store_pd((double *)t, zero); t+=2; }
    #endif
}

// =============================================================
// And now here's the main routine

void InCoreConvolution::InCoreConvolve(Complex *FFTM, DFLOAT *CompressedD) {
	
	STDLOG(3, "Entering InCoreConvolve\n");

	#pragma omp parallel for schedule(static)
    for(int block=0;block<nblocks;block++) {
            
        Complex *FM, *FMabc, *FMap2bcm2, *FMabp2cm2;
        double  *FD, *FDabc, *FDap2bcm2, *FDabp2cm2;
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
            FM = &(FFTM[rmap(a,b,c) * cpdky_pad + xyzbegin]);
            int m = cmap(a,b,c)*padblocksize;
            FOR(xyz,0,blocksize-1) mcache[m + xyz] = FM[xyz];
        }

        // Apply the trace-free recursion to generate complete multipoles
        TRACEFREERECURSION(a,b,c,order) {
            FMabc       = &(mcache[cmap(a  ,b  ,c  ) * padblocksize]);
            FMap2bcm2   = &(mcache[cmap(a+2,b  ,c-2) * padblocksize]);
            FMabp2cm2   = &(mcache[cmap(a  ,b+2,c-2) * padblocksize]);
            // TODO: Could vectorize this, at some cost of obscurity
            FOR(xyz,0,blocksize-1) 
                FMabc[xyz] = - FMap2bcm2[xyz] - FMabp2cm2[xyz];
        }

        // Apply the scaling to the complete multipoles
        FORALL_COMPLETE_MULTIPOLES_BOUND(a,b,c,order)  {
            int cmap_abc = cmap(a,b,c);
            FM = &(mcache[cmap_abc * padblocksize]);
            Multiply_by_scalar(FM, mfactor[cmap_abc], blocksize_even);
            // double mf = mfactor[ cmap_abc ];
            // FOR(xyz,0,blocksize-1) FM[xyz] *= mf;
        }

        if(twoD){
            // Load the derivatives, stored as [ky: 0..kysize][m: 0..rml][kx: 0..cpdp1half]
            // ky and kx are stored as a full quadrant, rather than half

            FORALL_REDUCED_MULTIPOLES_BOUND(a,b,c,order) {

                int m = cmap(a,b,c)*padblocksize;
                // Rather than unpacking the (x,y) by mod, we'll do 
                // one and then increment
                int y = xyzbegin%kysize;
                int x = xyzbegin/kysize;
                FOR(xyz,0,blocksize-1) {
                    int xp = x; if(xp>cpdhalf) xp = cpd - x;
                    int yp = y; if(yp>cpdhalf) yp = cpd - y;  // kysize never exceeds cpdhalf in 2D

                    double xR = (double) CompressedD[ yp*rml*(cpd+1)/2 + rmap(a,b,c)*(cpd+1)/2 + xp ];

                    if( (yp!=y) && ((b&0x1)==1) ) xR *= -1;  // not triggered in 2D
                    if( (xp!=x) && ((a&0x1)==1) ) xR *= -1;
                    dcache[m + xyz] = xR;
                    y++;
                    if (y==kysize) { y = 0; x++; }  // Accomplish the mod wrap
                }
            }

        } else {  // 1D
            
            // Load the derivatives.  This requires some care.
            // The derivatives were only stored in one half-quadrant.
            // Swapping (x,y) and (i,j) makes no change; that completes the quadrant.
            // Then the other quadrants are related by parity.
            
            FORALL_REDUCED_MULTIPOLES_BOUND(a,b,c,order) {

                int m = cmap(a,b,c)*padblocksize;
                // Rather than unpacking the (x,y) by mod, we'll do 
                // one and then increment
                int y = xyzbegin%kysize;
                int x = xyzbegin/kysize;
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
                    if (y==kysize) { y = 0; x++; }  // Accomplish the mod wrap
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
        }

        // Recurse to restore the complete derivatives
        TRACEFREERECURSION(a,b,c,order) {
            FDabc       = &(dcache[cmap(a  ,b  ,c  ) * padblocksize]);
            FDap2bcm2   = &(dcache[cmap(a+2,b  ,c-2) * padblocksize]);
            FDabp2cm2   = &(dcache[cmap(a  ,b+2,c-2) * padblocksize]);
            FOR(xyz,0,blocksize-1) 
                FDabc[xyz] = - FDap2bcm2[xyz] - FDabp2cm2[xyz];
        }

        // Now compute the Taylors from M & D
        // We do this in two parts, since the D's are stored as real,
        // but are either pure real or pure imaginary, according
        // to their i+j+k parity.
        FORALL_REDUCED_MULTIPOLES_BOUND(a,b,c,order) {
            // FOR(xyz,0,blocksize-1) tcache[xyz] = 0;
            Set_Vector_to_Zero(tcache, blocksize_even);
            FORALL_COMPLETE_MULTIPOLES_BOUND(i,j,k,order-a-b-c) {
                // Odd parity only
                if( (((a+i)+(b+j)+(c+k))&0x1) == 1 ) {
                    FD  = &(dcache[cmap(a+i,b+j,c+k) * padblocksize]);
                    FM  = &(mcache[cmap(i  ,j  ,k  ) * padblocksize]);
                    // We know the rest of the k's are simple
                    for (; k<=order-a-b-c-i-j; k+=2, FD+=padblocksize*2, FM+=padblocksize*2)
                    FMAvector(tcache, FM, FD, blocksize_even);
                }
            }
            // Now multiply by i, before we accumulate the rest
            Multiply_by_I(tcache, blocksize_even);

            FORALL_COMPLETE_MULTIPOLES_BOUND(i,j,k,order-a-b-c) {
                // Even parity only 
                if( (((a+i)+(b+j)+(c+k))&0x1) ==0 ) {
                    FD  = &(dcache[cmap(a+i,b+j,c+k) * padblocksize]);
                    FM  = &(mcache[cmap(i  ,j  ,k  ) * padblocksize]);
                    // We know the rest of the k's are simple
                    for (; k<=order-a-b-c-i-j; k+=2, FD+=padblocksize*2, FM+=padblocksize*2)
                    FMAvector(tcache,FM, FD, blocksize_even); 
                }
            }
            Complex *FT = &(FFTM[rmap(a,b,c) * cpdky_pad + xyzbegin]);
            // FOR(xyz,0,blocksize-1) FT[xyz] = tcache[xyz];
            memcpy(FT, tcache, blocksize*sizeof(Complex));
        }
    }
	
	STDLOG(3, "Exiting InCoreConvolve\n");
	
}
