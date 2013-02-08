class OutofCoreConvolution { 
public: 
    
    OutofCoreConvolution(void) { } 
    ~OutofCoreConvolution(void) { } 

    int runtime_cpd;

    int runtime_order;
    int runtime_NearFieldRadius;
    int runtime_DerivativeExpansionRadius;
    int runtime_IsRamDisk;
    int runtime_DiskBufferSizeKB;

    int runtime_ConvolutionCacheSizeMB;
    int runtime_MaxConvolutionRAMMB;

    char runtime_DerivativesDirectory[1024];
    char runtime_MultipoleDirectory[1024];
    char runtime_TaylorDirectory[1024];

    char runtime_MultipolePrefix[1024];
    char runtime_TaylorPrefix[1024];

    void Convolve( int _runtime_cpd, int _runtime_order, int _runtime_NearFieldRadius, int _runtime_DerivativeExpansionRadius,
                   int _runtime_IsRamDisk, int _runtime_DiskBufferSizeKB, int _runtime_ConvolutionCacheSizeMB, int _runtime_MaxConvolutionRAMMB,
                   char *_runtime_DerivativesDirectory, char *_runtime_MultipoleDirectory, char *_runtime_TaylorDirectory,
                   char *_runtime_MultipolePrefix, char *_runtime_TaylorPrefix );

    int blocksize, zwidth;

    STimer ReadDerivatives;
    STimer ReadMultipoles;
    STimer WriteTaylor;

    STimer ForwardZFFTMultipoles;
    STimer InverseZFFTTaylor;

    STimer ConvolutionArithmetic;
    STimer ArraySwizzle;
    
    STimer ConvolveWallClock;

    unsigned long long int ReadDerivativesBytes;
    unsigned long long int ReadMultipolesBytes;
    unsigned long long int WriteTaylorBytes;
    unsigned long long int ops;
    unsigned long long int totalMemoryAllocated;

private:

    int cpd,order,rml,CompressedMultipoleLengthXY;

    void BlockConvolve(void);
    void WriteDiskTaylor(int z);
    void ReadDiskDerivatives(int z);
    void ReadDiskMultipoles(int z);

    int mapM[8192];
    int remap[8192];

    ReadDirect *RD_RDD;
    ReadDirect *RD_RDM;
    WriteDirect *WD_WDT;

    Complex *DiskBuffer;
    double *CompressedDerivatives;
    Complex *TemporarySpace;
};


void OutofCoreConvolution::ReadDiskMultipoles(int z) { 

    for(int x=0;x<cpd;x++) {
        char fn[1024];
        sprintf(fn,"%s/%s_%d", runtime_MultipoleDirectory, runtime_MultipolePrefix, mapM[(x +(cpd-1)/2)%cpd]  );

        ReadMultipoles.Start();
        size_t s = sizeof(Complex); s *= zwidth; s *= cpd*rml;
        size_t offset = z; offset *= cpd; offset *= rml; offset *= sizeof(Complex);
        RD_RDM->BlockingRead( fn, (char *) &(TemporarySpace[0]), s, offset );
        ReadMultipolesBytes += s;
        ReadMultipoles.Stop();

        ArraySwizzle.Start();

        for(int zb=0;zb<zwidth;zb++)
            for(int m=0;m<rml;m++)
                for(int y=0;y<cpd;y++)
                    DiskBuffer[zb*rml*cpd*cpd + m*cpd*cpd + x*cpd + y] = 
                        TemporarySpace[ zb*cpd*rml + y*rml + m ];

        ArraySwizzle.Stop();
    }
}

void OutofCoreConvolution::WriteDiskTaylor(int z) {

    for(int x=0;x<cpd;x++) {
        ArraySwizzle.Start();

        for(int zb=0;zb<zwidth;zb++)
            for(int m=0;m<rml;m++)
                for(int y=0;y<cpd;y++)
                    TemporarySpace[ zb*cpd*rml + y*rml + m ] = 
                        DiskBuffer[  zb*rml*cpd*cpd + m*cpd*cpd + x*cpd + y];

        ArraySwizzle.Stop();

        char fn[1024];
        sprintf(fn,"%s/%s_%d", runtime_TaylorDirectory, runtime_TaylorPrefix, remap[ (x + (cpd-1)/2)%cpd] );

        WriteTaylor.Start();
        size_t s = sizeof(Complex); s *= zwidth; s *= cpd*rml; 
        WD_WDT->BlockingAppend(fn, (char *) &(TemporarySpace[0]), s);
        WriteTaylorBytes += s;
        WriteTaylor.Stop(); 
    }
}

void OutofCoreConvolution::ReadDiskDerivatives(int z) { 
    char fn[1024];
    sprintf(fn,"%s/fourierspace_%d_%d_%d_%d_%d",
           runtime_DerivativesDirectory, cpd,order, runtime_NearFieldRadius, runtime_DerivativeExpansionRadius, z);

    ReadDerivatives.Start();
    size_t s;
    s = sizeof(double);
    s *= rml;
    s *= CompressedMultipoleLengthXY;
    printf("CompressedMultipoleLengthXY = %d \n", CompressedMultipoleLengthXY );
    RD_RDD->BlockingRead( fn, (char *) &(CompressedDerivatives[0]), s, 0 );
    ReadDerivativesBytes += s;
    ReadDerivatives.Stop();
}


void OutofCoreConvolution::BlockConvolve(void) {
    ConvolveWallClock.Start();

    int cml = ((order+1)*(order+2)*(order+3))/6; 
    int nprocs = omp_get_num_procs();
    int cacherambytes = runtime_ConvolutionCacheSizeMB*(1024*1024);

    blocksize = 0;
    for(blocksize=cpd*cpd;blocksize>=2;blocksize--) 
        if((cpd*cpd)%blocksize==0)  // 2.5 = 2 Complex (mcache,tcache) 1 double dcache
            if( nprocs*2.5*cml*blocksize*sizeof(Complex) < cacherambytes) break;

    fftw_plan plan_forward_1d[16];
    fftw_plan plan_backward_1d[16];

    InCoreConvolution ICC(order,cpd,blocksize);

    ops = ICC.ConvolutionArithmeticCount();

    Complex *in_1d = new Complex[nprocs*cpd];
    Complex *out_1d = new Complex[nprocs*cpd];

    for(int g=0;g<nprocs;g++) {
        plan_forward_1d[g] = fftw_plan_dft_1d(cpd, (fftw_complex *) &(in_1d[g*cpd]), (fftw_complex *) &(out_1d[g*cpd]), FFTW_FORWARD, FFTW_MEASURE);
        plan_backward_1d[g] = fftw_plan_dft_1d(cpd, (fftw_complex *) &(in_1d[g*cpd]), (fftw_complex *) &(out_1d[g*cpd]), FFTW_BACKWARD, FFTW_MEASURE);
    }

    for(int zblock=0;zblock<(cpd+1)/2;zblock+=zwidth) {
    
        ReadDiskMultipoles(zblock);

        for(int z=zblock;z<zblock+zwidth; z++) {
       
            ReadDiskDerivatives(z);

            Complex *Mtmp = &( DiskBuffer[ (z-zblock)*rml*cpd*cpd ] );


            ForwardZFFTMultipoles.Start();

            #pragma omp parallel for schedule(dynamic,1) 
            for(int m=0;m<rml;m++) {
                int g = omp_get_thread_num();
                for(int y=0;y<cpd;y++) {
                    for(int x=0;x<cpd;x++) 
                        in_1d[g*cpd + x] = Mtmp[m*cpd*cpd + x*cpd + y];
                    fftw_execute(plan_forward_1d[g]);
                    for(int x=0;x<cpd;x++) 
                        Mtmp[m*cpd*cpd + x*cpd +  y] = out_1d[g*cpd + x];
                }
            }

            ForwardZFFTMultipoles.Stop();


                ConvolutionArithmetic.Start();
            ICC.InCoreConvolve(Mtmp, CompressedDerivatives);
                ConvolutionArithmetic.Stop();

            InverseZFFTTaylor.Start();

            #pragma omp parallel for schedule(dynamic,1)
            for(int m=0;m<rml;m++) {
                int g = omp_get_thread_num();
                for(int y=0;y<cpd;y++) {
                    for(int x=0;x<cpd;x++) 
                        in_1d[g*cpd + x] = Mtmp[ m*cpd*cpd + x*cpd + y];
                    fftw_execute(plan_backward_1d[g]);
                    for(int x=0;x<cpd;x++) 
                        Mtmp[  m*cpd*cpd + x*cpd + y ]  = out_1d[g*cpd + x];
                }
            }

            InverseZFFTTaylor.Stop();
        }

        WriteDiskTaylor(zblock);
    }

    delete[] in_1d;
    delete[] out_1d;

    ConvolveWallClock.Stop();
}


void OutofCoreConvolution::Convolve( int _runtime_cpd, int _runtime_order, int _runtime_NearFieldRadius, int _runtime_DerivativeExpansionRadius, 
                                     int _runtime_IsRamDisk, int _runtime_DiskBufferSizeKB, int _runtime_ConvolutionCacheSizeMB, int _runtime_MaxConvolutionRAMMB,
                                     char *_runtime_DerivativesDirectory,  char *_runtime_MultipoleDirectory, char *_runtime_TaylorDirectory,
                                     char *_runtime_MultipolePrefix, char *_runtime_TaylorPrefix ) {

    assert(_runtime_cpd%2==1); assert(_runtime_cpd>0);
    assert(_runtime_order<=16); assert(_runtime_order>=1);
    assert(_runtime_NearFieldRadius>=1);
    assert(_runtime_NearFieldRadius<=(_runtime_cpd-1)/2);
    assert( (_runtime_DerivativeExpansionRadius==8) || (_runtime_DerivativeExpansionRadius==16) );
    assert( (_runtime_IsRamDisk == 1) || (_runtime_IsRamDisk==0) );
    assert(_runtime_DiskBufferSizeKB>=1);
    assert( _runtime_ConvolutionCacheSizeMB >= 1);
    assert( _runtime_MaxConvolutionRAMMB >= 1);

    ReadDerivatives.Clear();
    ReadMultipoles.Clear();
    WriteTaylor.Clear();
    ForwardZFFTMultipoles.Clear();
    InverseZFFTTaylor.Clear();
    ConvolutionArithmetic.Clear();
    ArraySwizzle.Clear();
    ConvolveWallClock.Clear();
    

    ReadDerivativesBytes=0;
    ReadMultipolesBytes=0;
    WriteTaylorBytes=0;
    ops=0;
    totalMemoryAllocated=0;

    runtime_cpd                         = _runtime_cpd;
    runtime_order                       = _runtime_order;
    runtime_NearFieldRadius             = _runtime_NearFieldRadius;
    runtime_DerivativeExpansionRadius   = _runtime_DerivativeExpansionRadius;
    runtime_IsRamDisk                   = _runtime_IsRamDisk;
    runtime_DiskBufferSizeKB            = _runtime_DiskBufferSizeKB;
    runtime_ConvolutionCacheSizeMB      = _runtime_ConvolutionCacheSizeMB; 
    runtime_MaxConvolutionRAMMB         = _runtime_MaxConvolutionRAMMB;

    sprintf(runtime_DerivativesDirectory, "%s", _runtime_DerivativesDirectory);
    sprintf(runtime_MultipoleDirectory,   "%s", _runtime_MultipoleDirectory);
    sprintf(runtime_TaylorDirectory,      "%s", _runtime_TaylorDirectory);

    sprintf(runtime_MultipolePrefix,      "%s", _runtime_MultipolePrefix);
    sprintf(runtime_TaylorPrefix,         "%s", _runtime_TaylorPrefix);

    CheckDirectoryExists(runtime_TaylorDirectory);
    CheckDirectoryExists(runtime_MultipoleDirectory);
    CheckDirectoryExists(runtime_DerivativesDirectory);

    long int  sdb = runtime_DiskBufferSizeKB;
    sdb *= 1024;

    int direct = 1; if(runtime_IsRamDisk==1) direct = 0; 

    RD_RDD = new ReadDirect(direct,sdb);
    RD_RDM = new ReadDirect(direct,sdb);
    WD_WDT = new WriteDirect(direct,sdb);

    cpd = runtime_cpd;
    order = runtime_order;    

    CompressedMultipoleLengthXY  = ((1+cpd)*(3+cpd))/8;

    for(int i=0;i<cpd;i++) mapM[ (i+(cpd+1)/2)%cpd ] = i;

    unsigned long long int rambytes = runtime_MaxConvolutionRAMMB;
    rambytes *= (1024*1024);

    rml = (order+1)*(order+1);
    unsigned long long int zslabbytes = rml*cpd*cpd*sizeof(Complex);

    if(rambytes<zslabbytes) { 
        printf("Each slab requires      %lld bytes\n",zslabbytes);
        printf("You allow a maximum of  %lld bytes\n",rambytes);
        printf("[ERROR] rambytes<zlabbtes\n");
        exit(1);
    }

    for(int i=0;i<cpd;i++)    {
        int j = (i + (cpd+1)/2)%cpd;
        int k = (j+ (cpd-1)/2)%cpd;
        remap[j] = k;
    }

    int n = (cpd+1)/2;
    for(zwidth=n;zwidth>=2;zwidth--) {
        if(n%zwidth==0) if( zwidth*zslabbytes < rambytes) break;
    }

    ssize_t s;

    s = sizeof(Complex);
    s *= zwidth * rml * cpd * cpd;
    DiskBuffer = (Complex *) malloc(s);
    assert(DiskBuffer != NULL);
    totalMemoryAllocated += s;

    s = sizeof(double);
    s *= rml*CompressedMultipoleLengthXY;
    CompressedDerivatives  = (double *) malloc(s);
    assert(CompressedDerivatives != NULL);
    totalMemoryAllocated += s;

    s = sizeof(Complex);
    s *= zwidth * rml;
    s *= cpd;
    TemporarySpace  = (Complex *) malloc(s);
    assert( TemporarySpace != NULL);
    totalMemoryAllocated += s;

    for(int i=0;i<cpd;i++) {
        char cmd[1024];
        sprintf(cmd,"touch %s/Taylor_%d", runtime_TaylorDirectory, i);
        int rv = system(cmd);
        assert(rv!=-1);
    }

    BlockConvolve();

    delete RD_RDD;
    delete RD_RDM;
    delete WD_WDT;

    free(CompressedDerivatives);
    free(TemporarySpace);
    free(DiskBuffer);
}
