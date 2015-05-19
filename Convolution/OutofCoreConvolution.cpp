#include "InCoreConvolution.cpp"
#include "OutofCoreConvolution.h"

void OutofCoreConvolution::ReadDiskMultipoles(int z) { 
    for(int x=0;x<cpd;x++) {
        char fn[1024];
        sprintf(fn,"%s/%s_%d", CP.runtime_MultipoleDirectory, 
                               CP.runtime_MultipolePrefix, 
                               mapM[(x +(cpd-1)/2)%cpd]  );

        ReadMultipoles.Start();
        size_t s = sizeof(Complex); s *= zwidth; s *= cpd*rml;
        size_t offset = z; 
        offset *= cpd; offset *= rml; offset *= sizeof(Complex);
        RD_RDM->BlockingRead( fn, (char *) &(TemporarySpace[0]), s, offset );
        CS.ReadMultipolesBytes += s;
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
        sprintf(fn,"%s/%s_%d",  CP.runtime_TaylorDirectory, 
                                CP.runtime_TaylorPrefix, 
                                remap[ (x + (cpd-1)/2)%cpd] );

        WriteTaylor.Start();
        size_t s = sizeof(Complex); s *= zwidth; s *= cpd*rml; 
        WD_WDT->BlockingAppend(fn, (char *) &(TemporarySpace[0]), s);
        CS.WriteTaylorBytes += s;
        WriteTaylor.Stop(); 
    }
}


void OutofCoreConvolution::ReadDiskDerivatives(int z) { 
    char fn[1024];
    sprintf(fn,"%s/fourierspace_%d_%d_%d_%d_%d",
            CP.runtime_DerivativesDirectory, 
            cpd,order,CP.runtime_NearFieldRadius, 
            CP.runtime_DerivativeExpansionRadius, z);

    ReadDerivatives.Start();
    size_t s;
    s = sizeof(double);
    s *= rml;
    s *= CompressedMultipoleLengthXY;
    RD_RDD->BlockingRead( fn, (char *) &(CompressedDerivatives[0]), s, 0 );
    CS.ReadDerivativesBytes += s;
    ReadDerivatives.Stop();
}

void OutofCoreConvolution::BlockConvolve(void) {
    ConvolveWallClock.Start();

    int cml = ((order+1)*(order+2)*(order+3))/6; 
    int nprocs = omp_get_max_threads();
    size_t cacherambytes = CP.runtime_ConvolutionCacheSizeMB*(1024LL*1024LL);

    blocksize = 0;
    for(blocksize=cpd*cpd;blocksize>=2;blocksize--) 
        if((cpd*cpd)%blocksize==0)
            if( nprocs*2.5*cml*blocksize*sizeof(Complex) < cacherambytes) break;
                // 2.5 = 2 Complex (mcache,tcache) 1 double dcache
    #ifdef GPUFFT
    int ngpu;
    int iCPD = cpd;
    checkCudaErrors(cuda::cudaGetDeviceCount(&ngpu));
    cuda::cufftHandle plan_forward[ngpu];
    cuda::cufftHandle plan_backward[ngpu];
    
    cuda::cufftDoubleComplex * in_1d[ngpu];
    cuda::cufftDoubleComplex * out_1d[ngpu];
    cuda::cudaStream_t dev_stream[ngpu];

    for(int g = 0; g <ngpu; g++){
        checkCudaErrors(cuda::cudaSetDevice(g));
        checkCudaErrors(cuda::cudaMalloc((void **)(&in_1d[g]),cpd*cpd*sizeof(cuda::cufftDoubleComplex)));
        checkCudaErrors(cuda::cudaMalloc((void **)(&out_1d[g]),cpd*cpd*sizeof(cuda::cufftDoubleComplex)));
        checkCudaErrors(cuda::cudaStreamCreate(&dev_stream[g]));
    }

    #else
    fftw_plan plan_forward_1d[nprocs];
    fftw_plan plan_backward_1d[nprocs];

    Complex *in_1d = new Complex[nprocs*cpd];
    Complex *out_1d = new Complex[nprocs*cpd];

    #endif

    InCoreConvolution ICC(order,cpd,blocksize);

    CS.ops = ICC.ConvolutionArithmeticCount();

    #ifdef GPUFFT
    for(int g = 0; g <ngpu; g++){
        checkCudaErrors(cuda::cudaSetDevice(g));
        int icpd2 = cpd*cpd;
        checkCufftErrors(cuda::cufftPlanMany(&(plan_forward[g]) ,1,&iCPD,&icpd2,cpd,1,&icpd2,cpd,1,cuda::CUFFT_Z2Z,cpd));
        checkCufftErrors(cuda::cufftPlanMany(&(plan_backward[g]),1,&iCPD,&icpd2,cpd,1,&icpd2,cpd,1,cuda::CUFFT_Z2Z,cpd));
        checkCufftErrors(cuda::cufftSetStream(plan_forward[g],dev_stream[g]));
        checkCufftErrors(cuda::cufftSetStream(plan_backward[g],dev_stream[g]));
    }

    #else
    for(int g=0;g<nprocs;g++) {
        plan_forward_1d[g] = fftw_plan_dft_1d(cpd, 
                                (fftw_complex *) &(in_1d[g*cpd]), 
                                (fftw_complex *) &(out_1d[g*cpd]), 
                                FFTW_FORWARD, FFTW_MEASURE);
        plan_backward_1d[g] = fftw_plan_dft_1d(cpd, 
                                (fftw_complex *) &(in_1d[g*cpd]), 
                                (fftw_complex *) &(out_1d[g*cpd]), 
                                FFTW_BACKWARD, FFTW_MEASURE);
    }
    #endif

    for(int zblock=0;zblock<(cpd+1)/2;zblock+=zwidth) {
    	if (zblock +zwidth >= (cpd+1)/2) zwidth = (cpd+1)/2 -zblock;
        ReadDiskMultipoles(zblock);
        for(int z=zblock;z<zblock+zwidth; z++) {
	           
            ReadDiskDerivatives(z);

            Complex *Mtmp = &( DiskBuffer[ (z-zblock)*rml*cpd*cpd ] );

            ForwardZFFTMultipoles.Start();
            #ifdef GPUFFT
            for(int m=0;m<rml;m++) {
                int g = m%ngpu;
                checkCudaErrors(cuda::cudaSetDevice(g));
                checkCudaErrors(cuda::cudaMemcpyAsync(in_1d[g],&(Mtmp[m*cpd*cpd]),
                            cpd*cpd*sizeof(cuda::cufftDoubleComplex),cuda::cudaMemcpyHostToDevice,dev_stream[g] ));
                checkCufftErrors(cuda::cufftExecZ2Z(plan_forward[g],in_1d[g],out_1d[g],CUFFT_FORWARD));
                checkCudaErrors(cuda::cudaMemcpyAsync(&(Mtmp[m*cpd*cpd]),out_1d[g],
                            cpd*cpd*sizeof(cuda::cufftDoubleComplex),cuda::cudaMemcpyDeviceToHost,dev_stream[g] ));
            }
            for(int g = 0; g <ngpu; g++){
                checkCudaErrors(cuda::cudaSetDevice(g));
                cudaStreamSynchronize(dev_stream[g]);
            }
            #else
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
            #endif
            ForwardZFFTMultipoles.Stop();

            ConvolutionArithmetic.Start();
            ICC.InCoreConvolve(Mtmp, CompressedDerivatives);
            ConvolutionArithmetic.Stop();

            InverseZFFTTaylor.Start();
            #ifdef GPUFFT
            for(int m=0;m<rml;m++) {
                int g = m%ngpu;
                checkCudaErrors(cuda::cudaSetDevice(g));
                checkCudaErrors(cuda::cudaMemcpyAsync(in_1d[g],&(Mtmp[m*cpd*cpd]),
                            cpd*cpd*sizeof(cuda::cufftDoubleComplex),cuda::cudaMemcpyHostToDevice,dev_stream[g] ));
                checkCufftErrors(cuda::cufftExecZ2Z(plan_forward[g],in_1d[g],out_1d[g],CUFFT_INVERSE));
                checkCudaErrors(cuda::cudaMemcpyAsync(&(Mtmp[m*cpd*cpd]),out_1d[g],
                            cpd*cpd*sizeof(cuda::cufftDoubleComplex),cuda::cudaMemcpyDeviceToHost,dev_stream[g] ));
            }
            for(int g = 0; g <ngpu; g++){
                checkCudaErrors(cuda::cudaSetDevice(g));
                cudaStreamSynchronize(dev_stream[g]);
            }
            #else
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
            #endif
            InverseZFFTTaylor.Stop();
        }

        WriteDiskTaylor(zblock);
    }
    #ifdef GPUFFT
    for(int g = 0; g <ngpu; g++){
        checkCudaErrors(cuda::cudaSetDevice(g));
        cuda::cufftDestroy(plan_forward[g]);
        cuda::cufftDestroy(plan_backward[g]);
        cuda::cudaFree(in_1d[g]);
        cuda::cudaFree(out_1d[g]);
        cuda::cudaStreamDestroy(dev_stream[g]);
        cuda::cudaDeviceReset();
    }
    #else
    delete[] in_1d;
    delete[] out_1d;
    for(int g=0;g<nprocs;g++) {
        fftw_destroy_plan(plan_forward_1d[g]);
        fftw_destroy_plan(plan_backward_1d[g]);
    }
    #endif
    

    ConvolveWallClock.Stop();
}

void OutofCoreConvolution::Convolve( ConvolutionParameters _CP ) {

    CP = _CP;

    assert(CP.runtime_cpd%2==1); assert(CP.runtime_cpd>0);
    assert(CP.runtime_order<=16); assert(CP.runtime_order>=1);
    assert(CP.runtime_NearFieldRadius>=1);
    assert(CP.runtime_NearFieldRadius<=(CP.runtime_cpd-1)/2);

    assert( (CP.runtime_DerivativeExpansionRadius==8) || 
            (CP.runtime_DerivativeExpansionRadius==16) );
    assert( (CP.runtime_IsRamDisk == 1) || (CP.runtime_IsRamDisk==0) );
    assert(CP.runtime_DiskBufferSizeKB>=1);
    assert(CP.runtime_ConvolutionCacheSizeMB >= 1);
    assert(CP.runtime_MaxConvolutionRAMMB >= 1);

    ReadDerivatives.Clear();
    ReadMultipoles.Clear();
    WriteTaylor.Clear();
    ForwardZFFTMultipoles.Clear();
    InverseZFFTTaylor.Clear();
    ConvolutionArithmetic.Clear();
    ArraySwizzle.Clear();
    ConvolveWallClock.Clear();
    
    CS.ReadDerivativesBytes=0;
    CS.ReadMultipolesBytes=0;
    CS.WriteTaylorBytes=0;
    CS.ops=0;
    CS.totalMemoryAllocated=0;

    CS.runtime_ConvolutionCacheSizeMB = CP.runtime_ConvolutionCacheSizeMB;
    CS.runtime_cpd      = CP.runtime_cpd;
    CS.runtime_order    = CP.runtime_order;

    CheckDirectoryExists(CP.runtime_TaylorDirectory);
    CheckDirectoryExists(CP.runtime_MultipoleDirectory);
    CheckDirectoryExists(CP.runtime_DerivativesDirectory);

    size_t  sdb = CP.runtime_DiskBufferSizeKB;
    sdb *= 1024LLU;

    int direct = CP.runtime_IsRamDisk; 

    RD_RDD = new ReadDirect(direct,sdb);
    RD_RDM = new ReadDirect(direct,sdb);
    WD_WDT = new WriteDirect(direct,sdb);

    cpd = CP.runtime_cpd;
    order = CP.runtime_order;    

    CompressedMultipoleLengthXY  = ((1+cpd)*(3+cpd))/8;

    for(int i=0;i<cpd;i++) mapM[ (i+(cpd+1)/2)%cpd ] = i;

    unsigned long long int rambytes = CP.runtime_MaxConvolutionRAMMB;
    rambytes *= (1024LLU*1024LLU);

    rml = (order+1)*(order+1);
    unsigned long long int zslabbytes = rml*cpd*cpd*sizeof(Complex);
    printf("Each slab requires      %lld bytes\n",zslabbytes);
    printf("You allow a maximum of  %lld bytes\n",rambytes);
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
        /*if(n%zwidth==0)*/ if( zwidth*zslabbytes < rambytes) break;
    }
    printf("Resulting zwidth: %llu \n",zwidth);
    ssize_t s;

    s = sizeof(Complex);
    s *= zwidth * rml * cpd * cpd;
    posix_memalign((void **) &DiskBuffer, 4096,s);
    //DiskBuffer = (Complex *) malloc(s);
    if(DiskBuffer == NULL) printf("Tried to alloc aligned %llu bytes\n",s);
    assert(DiskBuffer != NULL);
    CS.totalMemoryAllocated += s;

    s = sizeof(double);
    s *= rml*CompressedMultipoleLengthXY;
    posix_memalign((void **) &CompressedDerivatives,4096,s);
    //CompressedDerivatives  = (double *) malloc(s);
    assert(CompressedDerivatives != NULL);
    CS.totalMemoryAllocated += s;

    s = sizeof(Complex);
    s *= zwidth * rml;
    s *= cpd;
    posix_memalign((void **) &TemporarySpace,4096,s);
    //TemporarySpace  = (Complex *) malloc(s);
    assert( TemporarySpace != NULL);
    CS.totalMemoryAllocated += s;

    for(int i=0;i<cpd;i++) {
        char cmd[1024];
        sprintf(cmd,"touch %s/Taylor_%d", CP.runtime_TaylorDirectory, i);
        int rv = system(cmd);
        assert(rv!=-1);
    }

    BlockConvolve();

    CS.blocksize             = blocksize;
    CS.zwidth                = zwidth;
    CS.ReadDerivatives       = ReadDerivatives.Elapsed();
    CS.ReadMultipoles        = ReadMultipoles.Elapsed();
    CS.WriteTaylor           = WriteTaylor.Elapsed();
    CS.ForwardZFFTMultipoles = ForwardZFFTMultipoles.Elapsed();
    CS.InverseZFFTTaylor     = InverseZFFTTaylor.Elapsed();
    CS.ConvolutionArithmetic = ConvolutionArithmetic.Elapsed();
    CS.ArraySwizzle          = ArraySwizzle.Elapsed();
    CS.ConvolveWallClock     = ConvolveWallClock.Elapsed();

    double accountedtime  = CS.ConvolutionArithmetic;
           accountedtime += CS.ForwardZFFTMultipoles + CS.InverseZFFTTaylor;
           accountedtime += CS.ReadDerivatives + CS.ReadMultipoles + CS.WriteTaylor;
           accountedtime += CS.ArraySwizzle;

    CS.Discrepency = CS.ConvolveWallClock - accountedtime;

    delete RD_RDD;
    delete RD_RDM;
    delete WD_WDT;

    free(CompressedDerivatives);
    free(TemporarySpace);
    free(DiskBuffer);
}
