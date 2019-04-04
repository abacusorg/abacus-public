class InCoreConvolution : public basemultipoles {
public:
    InCoreConvolution(int order, int cpd, int blocksize) : 
            basemultipoles(order), cpd(cpd), blocksize(blocksize) {

        CompressedMultipoleLengthXY  = ((1+cpd)*(3+cpd))/8;
        nblocks = (cpd*cpd)/blocksize;
        cpdhalf = (cpd-1)/2;
        int cs = omp_get_max_threads() * (completemultipolelength * blocksize + thread_padding);  // extra 1K padding between threads
        _mcache = new Complex[cs]; 
        _dcache = new  double[cs]; 
        _tcache = new Complex[cs];
        mfactor = new double[completemultipolelength];
        int a,b,c;
        FORALL_COMPLETE_MULTIPOLES_BOUND(a,b,c,order) {
            mfactor[ cmap(a,b,c) ] = 
                pow(-1.0,a+b+c)/(fact2(2*(a+b+c)-1)*fact(a)*fact(b)*fact(c));
        }
    }
    ~InCoreConvolution(void) {
        delete[] _mcache; delete[] _dcache; delete[] _tcache; delete[] mfactor;
    }
    void InCoreConvolve(Complex *FFTM, DFLOAT *CompressedDerivatives);
    unsigned long long int  ConvolutionArithmeticCount(void);

private:
    Complex *_mcache,*_tcache;
    double *_dcache,*mfactor;
    int cpd,blocksize,nblocks,cpdhalf,CompressedMultipoleLengthXY;
    const int thread_padding = 1024;
};

unsigned long long int InCoreConvolution::ConvolutionArithmeticCount(void) {
    int a,b,c,i,j,k;

    unsigned long long int l = 0;

    unsigned long long int cpd3 = ((unsigned long long int ) cpd)*cpd*cpd;

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

void MVM(Complex *t, Complex *m, double *d, int n) { 
    //#pragma simd
    // __builtin_assume_aligned
    //#pragma vector aligned
    #pragma simd assert
    for(int i=0;i<n;i++) t[i] += m[i] * d[i]; 
}

void InCoreConvolution::InCoreConvolve(Complex *FFTM, DFLOAT *CompressedD) {
    // why is this 50% faster with dynamic?
    #pragma omp parallel for schedule(dynamic)
    for(int block=0;block<nblocks;block++) {
            
        Complex *FM, *FMabc, *FMap2bcm2, *FMabp2cm2;
        double  *FD, *FDabc, *FDap2bcm2, *FDabp2cm2;
        Complex *mcache, *tcache;
        double  *dcache;
        int xyz,a,b,c,i,j,k,xyzbegin;

        int gindex = omp_get_thread_num()*(blocksize*completemultipolelength + thread_padding);
        dcache = &(_dcache[gindex]); 
        mcache = &(_mcache[gindex]);
        tcache = &(_tcache[gindex]);

        xyzbegin = block*blocksize;
        FORALL_REDUCED_MULTIPOLES_BOUND(a,b,c,order) {
            FM = &(FFTM[rmap(a,b,c) * cpd*cpd + xyzbegin]);
            int m = cmap(a,b,c);
            FOR(xyz,0,blocksize-1) mcache[m*blocksize + xyz] = FM[xyz];
        }

        TRACEFREERECURSION(a,b,c,order) {
            FMabc       = &(mcache[cmap(a  ,b  ,c  ) * blocksize]);
            FMap2bcm2   = &(mcache[cmap(a+2,b  ,c-2) * blocksize]);
            FMabp2cm2   = &(mcache[cmap(a  ,b+2,c-2) * blocksize]);
            FOR(xyz,0,blocksize-1) 
                FMabc[xyz] = - FMap2bcm2[xyz] - FMabp2cm2[xyz];
        }

        FORALL_COMPLETE_MULTIPOLES_BOUND(a,b,c,order)  {
            FM = &(mcache[cmap(a,b,c) * blocksize]);
            double mf = mfactor[ cmap(a,b,c) ];
            FOR(xyz,0,blocksize-1) FM[xyz] *= mf;
        }

        FORALL_REDUCED_MULTIPOLES_BOUND(a,b,c,order) {

            int m = cmap(a,b,c);
            FOR(xyz,0,blocksize-1) {
                int l = xyzbegin + xyz;
                int y = l%cpd;
                int x = (l-y)/cpd;
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

                if( (yp!=y) && (b%2==1) ) xR *= -1;
                if( (xp!=x) && (a%2==1) ) xR *= -1;
                dcache[m*blocksize + xyz] = xR;
            }
        }

        TRACEFREERECURSION(a,b,c,order) {
            FDabc       = &(dcache[cmap(a  ,b  ,c  ) * blocksize]);
            FDap2bcm2   = &(dcache[cmap(a+2,b  ,c-2) * blocksize]);
            FDabp2cm2   = &(dcache[cmap(a  ,b+2,c-2) * blocksize]);
            FOR(xyz,0,blocksize-1) 
                FDabc[xyz] = - FDap2bcm2[xyz] - FDabp2cm2[xyz];
        }

        FORALL_REDUCED_MULTIPOLES_BOUND(a,b,c,order) {
            Complex *FT = &(FFTM[rmap(a,b,c) * cpd*cpd + xyzbegin]);

            FOR(xyz,0,blocksize-1) tcache[xyz] = 0;
            FORALL_COMPLETE_MULTIPOLES_BOUND(i,j,k,order-a-b-c) {
                if( ((a+i)+(b+j)+(c+k))%2 ==0 ) {
                    FD  = &(dcache[cmap(a+i,b+j,c+k) * blocksize]);
                    FM  = &(mcache[cmap(i  ,j  ,k  ) * blocksize]);
                    MVM(tcache, FM, FD, blocksize);
                }
            }
            FOR(xyz,0,blocksize-1) FT[xyz] = tcache[xyz];

            FOR(xyz,0,blocksize-1) tcache[xyz] = 0;
            FORALL_COMPLETE_MULTIPOLES_BOUND(i,j,k,order-a-b-c) {
                if( ((a+i)+(b+j)+(c+k))%2 == 1 ) {
                    FD  = &(dcache[cmap(a+i,b+j,c+k) * blocksize]);
                    FM  = &(mcache[cmap(i  ,j  ,k  ) * blocksize]);
                    MVM(tcache,FM, FD, blocksize); 
                }
            }
            Complex I(0,1); FOR(xyz,0,blocksize-1) FT[xyz] += tcache[xyz]*I;
        }
    }
}
