#include "InCoreConvolution.cpp"
#include "OutofCoreConvolution.h"

void OutofCoreConvolution::ReadDiskMultipolesAndDerivs(int zstart) {
#ifdef CONVIOTHREADED
    // We have to wait for all IO threads, because they are all reading into the same block
    for(int i = 0; i < CP.niothreads; i++){
        WaitForIO.Start();
        // All the threads return the same block; doesn't matter that we re-pop
        iothreads[i]->pop(CurrentBlock);
        WaitForIO.Stop();
    }
    
    PlaneBuffer = CurrentBlock->mtblock;
    CompressedDerivatives = CurrentBlock->dblock;
#else
    CurrentBlock->read(zstart, zwidth, 0);
    CurrentBlock->read_derivs(zstart, zwidth, 0);
#endif
}


// Select a single z-plane out of the x-plane-ordered multipoles block
void OutofCoreConvolution::SwizzleMultipoles(int z){
    ArraySwizzle.Start();
    #pragma omp parallel for schedule(static)
    for(int m=0;m<rml;m++)
        for(int x=0;x<cpd;x++)
            for(int y=0;y<cpd;y++)
                DiskBuffer[m*cpd*cpd + x*cpd + y] = 
                    PlaneBuffer[x][z*cpd*rml + m*cpd + y ];
    ArraySwizzle.Stop();
}


// Pack a single z-plane into the x-plane-ordered Taylors block
void OutofCoreConvolution::SwizzleTaylors(int z){
    // Array swizzle selects the plane belonging to a single slab/file from a block
    // Swizzle does actually benefit from multithreading, so we don't want to do that in the IO thread
    ArraySwizzle.Start();
    // These loops go from outer to inner for the destination array
    #pragma omp parallel for schedule(static)
    for(int x=0;x<cpd;x++) {
        for(int m=0;m<rml;m++)
            for(int y=0;y<cpd;y++)
                PlaneBuffer[x][z*cpd*rml + m*cpd + y ] = 
                    DiskBuffer[ m*cpd*cpd + x*cpd + y]*invcpd3;
    }
    ArraySwizzle.Stop();
}


void OutofCoreConvolution::WriteDiskTaylor(int z) {
#ifdef CONVIOTHREADED
    // Push to all threads.  Each will take the part of the block it needs.
    for(int i = 0; i < CP.niothreads; i++)
        iothreads[i]->push(CurrentBlock);
#else
    CurrentBlock->write(zwidth, 0);
#endif
}

void OutofCoreConvolution::BlockConvolve(void) {
    int nprocs = omp_get_max_threads();
    size_t cacherambytes = CP.runtime_ConvolutionCacheSizeMB*(1024LL*1024LL);

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
    fftw_plan *plan_forward_1d = new fftw_plan[nprocs];
    fftw_plan *plan_backward_1d = new fftw_plan[nprocs];

    Complex **in_1d = new Complex*[nprocs];
    Complex **out_1d = new Complex*[nprocs];
    for(int g = 0; g < nprocs; g++){
        in_1d[g] = new Complex[cpd];
        out_1d[g] = new Complex[cpd];
    }

    #endif

    InCoreConvolution *ICC = new InCoreConvolution(order,cpd,blocksize);

    CS.ops = ICC->ConvolutionArithmeticCount();

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
                                (fftw_complex *) in_1d[g], 
                                (fftw_complex *) out_1d[g], 
                                FFTW_FORWARD, FFTW_PATIENT);
        
        plan_backward_1d[g] = fftw_plan_dft_1d(cpd, 
                                (fftw_complex *) in_1d[g], 
                                (fftw_complex *) out_1d[g], 
                                FFTW_BACKWARD, FFTW_PATIENT);
    }
    #endif

    for(int zblock = 0; zblock < (cpd + 1)/2; zblock += zwidth) {
    	if (zblock + zwidth >= (cpd + 1)/2)
            zwidth = (cpd+1)/2 - zblock;
        STDLOG(1, "Starting z %d to %d\n", zblock, zblock + zwidth);
        
        ReadDiskMultipolesAndDerivs(zblock);
        
        for(int z = zblock; z < zblock + zwidth; z++) {
            SwizzleMultipoles(z - zblock);

            Complex *Mtmp = &( DiskBuffer[0] );

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
            #pragma omp parallel for schedule(static)
            for(int m=0;m<rml;m++) {
                int g = omp_get_thread_num();
                for(int y=0;y<cpd;y++) {
                    for(int x=0;x<cpd;x++) 
                        in_1d[g][x] = Mtmp[m*cpd*cpd + x*cpd + y];
                    fftw_execute(plan_forward_1d[g]);
                    for(int x=0;x<cpd;x++) 
                        Mtmp[m*cpd*cpd + x*cpd + y] = out_1d[g][x];
                }
            }
            
            #endif
            ForwardZFFTMultipoles.Stop();

            ConvolutionArithmetic.Start();
            ICC->InCoreConvolve(Mtmp, CompressedDerivatives[z-zblock]);
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
            #pragma omp parallel for schedule(static)
            for(int m=0;m<rml;m++) {
                int g = omp_get_thread_num();
                for(int y=0;y<cpd;y++) {
                    for(int x=0;x<cpd;x++) 
                        in_1d[g][x] = Mtmp[m*cpd*cpd + x*cpd + y];
                    fftw_execute(plan_backward_1d[g]);
                    for(int x=0;x<cpd;x++) 
                        Mtmp[m*cpd*cpd + x*cpd + y] = out_1d[g][x];
                }
            }
            #endif
            InverseZFFTTaylor.Stop();
            
            SwizzleTaylors(z - zblock);
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
    for(int g=0;g<nprocs;g++) {
        fftw_destroy_plan(plan_forward_1d[g]);
        fftw_destroy_plan(plan_backward_1d[g]);
        delete[] in_1d[g];
        delete[] out_1d[g];
    }
    delete[] in_1d;
    delete[] out_1d;
    
    delete[] plan_forward_1d;
    delete[] plan_backward_1d;
    #endif
    
    delete ICC;

}

void OutofCoreConvolution::Convolve( ConvolutionParameters _CP ) {

    CP = _CP;

    assert(CP.runtime_cpd%2==1); assert(CP.runtime_cpd>0);
    assert(CP.runtime_order<=16); assert(CP.runtime_order>=1);
    assert(CP.runtime_NearFieldRadius>=1);
    assert(CP.runtime_NearFieldRadius<=(CP.runtime_cpd-1)/2);

    assert(	(CP.runtime_DerivativeExpansionRadius==1) ||
			(CP.runtime_DerivativeExpansionRadius==2) ||
			(CP.runtime_DerivativeExpansionRadius==3) ||
			(CP.runtime_DerivativeExpansionRadius==4) || 
			(CP.runtime_DerivativeExpansionRadius==5) || //note! added to change far radius 
			(CP.runtime_DerivativeExpansionRadius==6) || 
			(CP.runtime_DerivativeExpansionRadius==7) ||
			(CP.runtime_DerivativeExpansionRadius==8) || 
            (CP.runtime_DerivativeExpansionRadius==16) );
    assert( (CP.runtime_IsRamDisk == 1) || (CP.runtime_IsRamDisk==0) );
    assert(CP.runtime_DiskBufferSizeKB>=1);
    assert(CP.runtime_ConvolutionCacheSizeMB >= 1);
    assert(CP.runtime_MaxConvolutionRAMMB >= 1);

    ForwardZFFTMultipoles.Clear();
    InverseZFFTTaylor.Clear();
    ConvolutionArithmetic.Clear();
    ArraySwizzle.Clear();
    
    CS.ReadDerivativesBytes=0;
    CS.ReadMultipolesBytes=0;
    CS.WriteTaylorBytes=0;
    CS.ops=0;
    CS.totalMemoryAllocated=0;

    CS.runtime_ConvolutionCacheSizeMB = CP.runtime_ConvolutionCacheSizeMB;

    CheckDirectoryExists(CP.runtime_TaylorDirectory);
    CheckDirectoryExists(CP.runtime_MultipoleDirectory);
    CheckDirectoryExists(CP.runtime_DerivativesDirectory);
    
    cpd = CP.runtime_cpd;
    order = CP.runtime_order;
    invcpd3 = 1./cpd/cpd/cpd;
    blocksize = CP.blocksize;
    CompressedMultipoleLengthXY = CP.CompressedMultipoleLengthXY;
    rml = CP.rml;
    zwidth = CP.zwidth;
    
    // Check that all the multipole files exists
    // This is to prevent losing Taylors by accidentally convolving after an interrupted singlestep
    for(int i=0;i<cpd;i++) {
        char fn[1024];
        CP.MultipoleDirectory(i, fn);
        CheckFileExists(fn);
    }

    size_t sdb = CP.runtime_DiskBufferSizeKB;
    sdb *= 1024LLU;

    int direct = CP.runtime_IsRamDisk; 

    RD_RDD = new ReadDirect(direct,sdb);
    RD_RDM = new ReadDirect(direct,sdb);
    WD_WDT = new WriteDirect(direct,sdb);
    
    STDLOG(0,"Resulting zwidth: %d \n",zwidth);
    size_t s;

    // Disk buffer will only hold one z-slab at a time
    s = sizeof(Complex) * rml * cpd * cpd;
    int memalign_ret = posix_memalign((void **) &DiskBuffer, 4096,s);
    if(DiskBuffer == NULL) printf("Tried to alloc aligned %.2f MB\n", (double) s/1024/1024.);
    assert(memalign_ret == 0);
    assert(DiskBuffer != NULL);
    CS.totalMemoryAllocated += s;

#ifdef CONVIOTHREADED
    { // scope to avoid leakage
    // Allocate blocks and pass them to threads.
    // Blocks are shared among threads.
    int nblocks = (int) ceil((CP.runtime_cpd+1)/2./CP.zwidth);
    int read_ahead = min(2, nblocks);
    Block *blocks[read_ahead];
    for(int i = 0; i < read_ahead; i++){
        blocks[i] = new Block(CP);
        CS.totalMemoryAllocated += blocks[i]->alloc_bytes;
    }

    iothreads = new ConvIOThread*[CP.niothreads];
    for(int i = 0; i < CP.niothreads; i++){
        iothreads[i] = new ConvIOThread(CP, read_ahead, blocks, i);  // starts io thread
    }
    }
#else
    CurrentBlock = new Block(CP);
    CS.totalMemoryAllocated += CurrentBlock->alloc_bytes;
    
    PlaneBuffer = CurrentBlock->mtblock;
    CompressedDerivatives = CurrentBlock->dblock;
#endif

    // Create the Taylor files, erasing them if they existed
    for(int i=0;i<cpd;i++) {
        char fn[1024];
        CP.TaylorDirectory(i, fn);
        FILE *f = fopen(fn, "wb");
        assert(f != NULL);
        fclose(f);
    }
    
    BlockConvolve();

#ifdef CONVIOTHREADED
    for(int i = 0; i < CP.niothreads; i++){
        ConvIOThread *iothread = iothreads[i];
        WaitForIO.Start();
        iothread->join();  // wait for io to terminate
        WaitForIO.Stop();
        
        CS.ReadDerivatives       += iothread->get_deriv_read_time();
        CS.ReadMultipoles        += iothread->get_multipole_read_time();
        CS.WriteTaylor           += iothread->get_taylor_write_time();
        CS.ReadMultipolesBytes   += iothread->get_multipole_bytes_read();
        CS.WriteTaylorBytes      += iothread->get_taylor_bytes_written();
        CS.ReadDerivativesBytes  += iothread->get_derivative_bytes_read();
        delete iothread;
    }
    delete iothreads;
    CS.WaitForIO += WaitForIO.Elapsed();
#else
    CS.WriteTaylorBytes      = CurrentBlock->WriteTaylorBytes;
    CS.ReadMultipolesBytes   = CurrentBlock->ReadMultipoleBytes;
    CS.ReadDerivativesBytes  = CurrentBlock->ReadDerivativeBytes;
    CS.WriteTaylor           = CurrentBlock->WriteTaylor.Elapsed();
    CS.ReadMultipoles        = CurrentBlock->ReadMultipoles.Elapsed();
    CS.ReadDerivatives       = CurrentBlock->ReadDerivatives.Elapsed();
    delete CurrentBlock;
#endif

    CS.ComputeCores          = omp_get_max_threads();
    CS.ForwardZFFTMultipoles = ForwardZFFTMultipoles.Elapsed();
    CS.InverseZFFTTaylor     = InverseZFFTTaylor.Elapsed();
    CS.ConvolutionArithmetic = ConvolutionArithmetic.Elapsed();
    CS.ArraySwizzle          = ArraySwizzle.Elapsed();

    delete RD_RDD;
    delete RD_RDM;
    delete WD_WDT;

    free(DiskBuffer);
}
