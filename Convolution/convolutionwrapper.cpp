#include "header.cpp"
#include "threevector.hh"

#ifdef IOTHREADED
#define CONVIOTHREADED
#endif

#include "STimer.cc"
#include "PTimer.cc"

STimer TotalWallClock;
STimer Setup;
STimer ConvolutionWallClock;

#include "file.cpp"
#include "factorial.cpp"
#include "iolib.cpp"

#include "stdlog.cc"

#include "Parameters.cpp"

#ifdef GPUFFT
namespace cuda{
    #include <cuda.h>
    #include <cuda_runtime_api.h>
    #include <cufft.h>
}
#include "CufftErrors.h"
#else
#include <fftw3.h>
#endif

#include "threadaffinity.h"
#include "ConvolutionLibrary.cpp"
#include "basemultipoles.cpp"

#include <fenv.h>

void dumpstats(OutofCoreConvolution *OCC, char *fn) {

    FILE *fp;
    fp = fopen(fn,"w");
    assert(fp!=NULL);


    double accountedtime  = OCC->CS.ConvolutionArithmetic;
           accountedtime += OCC->CS.ForwardZFFTMultipoles + OCC->CS.InverseZFFTTaylor;
#ifdef CONVIOTHREADED
           accountedtime += OCC->CS.WaitForIO;
#else
           accountedtime += OCC->CS.ReadDerivatives + OCC->CS.ReadMultipoles + OCC->CS.WriteTaylor;
#endif
           accountedtime += OCC->CS.ArraySwizzle;
    double discrepency = OCC->CS.ConvolveWallClock - accountedtime;

    int computecores = OCC->CS.ComputeCores;
    fprintf(fp,"Convolution parameters:  RamAllocated = %dMB CacheSizeMB = %dMB nreal_cores=%d blocksize=%d zwidth=%d cpd=%d order=%d\n\n",
        (int) (OCC->CS.totalMemoryAllocated/(1<<20)), OCC->CS.runtime_ConvolutionCacheSizeMB, computecores, (int) OCC->CP.blocksize, (int) OCC->CP.zwidth, OCC->CP.runtime_cpd, OCC->CP.runtime_order);

    fprintf(fp,"\t ConvolutionWallClock:  %2.2e seconds \n", OCC->CS.ConvolveWallClock );
    fprintf(fp,"\t \t %50s : %1.1e seconds\n", "Array Swizzling", OCC->CS.ArraySwizzle );
    
#ifdef CONVIOTHREADED
    double e = OCC->CS.ReadDerivativesBytes/OCC->CS.ReadDerivatives/(1.0e+6);
    fprintf(fp,"\t \t %50s : %1.1e seconds --> rate was %4.0f MB/s\n", "ReadDiskDerivatives [per thread]", OCC->CS.ReadDerivatives, e );
    
    e = OCC->CS.ReadMultipolesBytes/OCC->CS.ReadMultipoles/(1.0e+6);
    fprintf(fp,"\t \t %50s : %1.1e seconds --> rate was %4.0f MB/s\n", "ReadDiskMultipoles [per thread]", OCC->CS.ReadMultipoles, e );
    
    e = OCC->CS.WriteTaylorBytes/OCC->CS.WriteTaylor/(1.0e+6);
    fprintf(fp,"\t \t %50s : %1.1e seconds --> rate was %4.0f MB/s\n", "WriteDiskTaylor [per thread]", OCC->CS.WriteTaylor, e );
    
    fprintf(fp,"\t \t %50s : %1.1e seconds\n", "Waiting for IO thread", OCC->CS.WaitForIO);
#else
    double e = OCC->CS.ReadDerivativesBytes/OCC->CS.ReadDerivatives/(1.0e+6);
    fprintf(fp,"\t \t %50s : %1.1e seconds --> rate was %4.0f MB/s\n", "ReadDiskDerivatives", OCC->CS.ReadDerivatives, e );
    
    e = OCC->CS.ReadMultipolesBytes/OCC->CS.ReadMultipoles/(1.0e+6);
    fprintf(fp,"\t \t %50s : %1.1e seconds --> rate was %4.0f MB/s\n", "ReadDiskMultipoles", OCC->CS.ReadMultipoles, e );
    
    e = OCC->CS.WriteTaylorBytes/OCC->CS.WriteTaylor/(1.0e+6);
    fprintf(fp,"\t \t %50s : %1.1e seconds --> rate was %4.0f MB/s\n", "WriteDiskTaylor", OCC->CS.WriteTaylor, e );
#endif
    
    double Gops = ((double) OCC->CS.ops)/(1.0e+9);
    fprintf(fp,"\t \t %50s : %1.1e seconds for %5.3f billion double precision operations\n", "Convolution Arithmetic", OCC->CS.ConvolutionArithmetic, Gops );
    
    fprintf(fp,"\t \t %50s : %1.1e seconds\n", "Forward FFT Z Multipoles", OCC->CS.ForwardZFFTMultipoles );
    fprintf(fp,"\t \t %50s : %1.1e seconds\n", "Inverse FFT Z Taylor",         OCC->CS.InverseZFFTTaylor );
    
    fprintf(fp,"\t %50s : %1.1e seconds which is %d%% \n", "Unaccounted remaining wallclock time", discrepency, (int) (discrepency/OCC->CS.ConvolveWallClock*100) );

    double cae = OCC->CS.ConvolutionArithmetic;
    double farithp   = cae/OCC->CS.ConvolveWallClock*100;
    double ffftp     = (OCC->CS.ForwardZFFTMultipoles + OCC->CS.InverseZFFTTaylor)/OCC->CS.ConvolveWallClock*100;
    double fiop      = (OCC->CS.ReadDerivatives + OCC->CS.ReadMultipoles + OCC->CS.WriteTaylor )/OCC->CS.ConvolveWallClock*100;
    double swzp      = OCC->CS.ArraySwizzle/OCC->CS.ConvolveWallClock*100;

    fprintf(fp,"\n \t Summary: Fourier Transforms = %2.0f%%     Convolution Arithmetic = %2.0f%%     Array Swizzle = %2.0f%%", ffftp, farithp, swzp );
#ifdef CONVIOTHREADED
    double fiow      = OCC->CS.WaitForIO/OCC->CS.ConvolveWallClock*100;
    fprintf(fp,"\n \t                                  Non-blocking Disk IO = %2.0f%%    Waiting for IO Thread = %2.0f%% \n", fiop, fiow);
#else
    fprintf(fp,"    Disk IO = %2.0f%% \n", fiop);
#endif
    fprintf(fp,"\t          Arithmetic rate = %2.0f DGOPS --> rate per core = %1.1f DGOPS\n", Gops/cae, Gops/cae/computecores );
    fprintf(fp,"\t          [DGOPS == Double Precision Billion operations per second]\n");
    fprintf(fp,"\n");

    fclose(fp);
}

void setup_openmp(){
    int max_threads = omp_get_max_threads();
    int ncores = omp_get_num_procs();
    int nthreads = P.Conv_OMP_NUM_THREADS > 0 ? P.Conv_OMP_NUM_THREADS : max_threads + P.Conv_OMP_NUM_THREADS;
#ifdef CONVIOTHREADED
    for (int i = 0; i < MAX_IO_THREADS; i++){
        P.Conv_IOCores[i] = P.Conv_IOCores[i];
        //STDLOG(2, "IO thread %d assigned to core %d\n", i, P.Conv_IOCores[i]);
    }
#endif
    
    assertf(nthreads <= max_threads, "Trying to use more OMP threads (%d) than omp_get_max_threads() (%d)!  This will cause global objects that have already used omp_get_max_threads() to allocate thread workspace (like PTimer) to fail.\n");
    assertf(nthreads <= ncores, "Trying to use more threads (%d) than cores (%d).  This will probably be very slow.\n", nthreads, ncores);
    
    omp_set_num_threads(nthreads);
    STDLOG(1, "Initializing OpenMP with %d threads (system max is %d; P.Conv_OMP_NUM_THREADS is %d)\n", nthreads, max_threads, P.Conv_OMP_NUM_THREADS);

    // If threads are bound to cores via OMP_PROC_BIND,
    // then identify free cores for use by IO thread
    if(omp_get_proc_bind() == omp_proc_bind_false){
        STDLOG(1, "OMP_PROC_BIND = false; threads will not be bound to cores\n");
    }
    else{
        int core_assignments[nthreads];
        #pragma omp parallel for schedule(static)
        for(int g = 0; g < nthreads; g++){
            assertf(g == omp_get_thread_num(), "OpenMP thread %d is executing wrong loop ieration (%d)\n", omp_get_thread_num(), g);
            core_assignments[g] = sched_getcpu();
        }
        std::ostringstream core_log;
        core_log << "Thread->core assignments:";
        for(int g = 0; g < nthreads; g++)
            core_log << " " << g << "->" << core_assignments[g];
        core_log << "\n";
        STDLOG(1, core_log.str().c_str());
        
        // Assign the main CPU thread to core 0 to avoid the IO thread during serial parts of the code
        /*int main_thread_core = P.Conv_IOCore != 0 ? 0 : 1;
        set_core_affinity(main_thread_core);
        STDLOG(1, "Assigning main convolution thread to core %d\n", main_thread_core);*/
    }

}

int main(int argc, char ** argv){
	TotalWallClock.Start();
	Setup.Start();

	if (argc!=2) {
	       // Can't use assertf() or QUIT here: stdlog not yet defined!
	       fprintf(stderr, "Error: command line must have 1 parameter given, not %d.\nLegal usage: %s PARAM_FILE\n", argc-1, argv[0]);
	       assert(0==99);
	    }

	    P.ReadParameters(argv[1],1);

	    // Setup the log
	    stdlog_threshold_global = P.LogVerbosity;
	    char logfn[1050];
	    sprintf(logfn,"%s/last.convlog", P.LogDirectory);
	    stdlog.open(logfn);
	    OutofCoreConvolution OCC;
	    STDLOG(1,"Read parameter file\n");
	    Setup.Stop();
        
        setup_openmp();

	    ConvolutionParameters p;
	    p.runtime_ConvolutionCacheSizeMB = P.ConvolutionCacheSizeMB;
        STDLOG(1, "Using cache size %d MB\n", p.runtime_ConvolutionCacheSizeMB);
	    p.runtime_DerivativeExpansionRadius = P.DerivativeExpansionRadius;
	    strcpy(p.runtime_DerivativesDirectory,P.DerivativesDirectory);
	    p.runtime_DiskBufferSizeKB = 1LL<<21;
	    p.runtime_IsRamDisk = P.RamDisk;
	    p.runtime_MaxConvolutionRAMMB = P.MAXRAMMB;
	    strcpy(p.runtime_MultipoleDirectory,P.MultipoleDirectory);

        // Multipole/TaylorDirectory2 will default to the primary directory
        p.niothreads = 1;
        if(strcmp(P.MultipoleDirectory2,STRUNDEF) == 0)
            strcpy(p.runtime_MultipoleDirectory2,P.MultipoleDirectory);
        else{
            strcpy(p.runtime_MultipoleDirectory2,P.MultipoleDirectory2);
#ifdef CONVIOTHREADED
            // Two IO threads if we were given two Multipole directories
            p.niothreads = 2;
#endif
        }
        if(strcmp(P.TaylorDirectory2,STRUNDEF) == 0)
            strcpy(p.runtime_TaylorDirectory2,P.TaylorDirectory);
        else{
            strcpy(p.runtime_TaylorDirectory2,P.TaylorDirectory2);
#ifdef CONVIOTHREADED
            p.niothreads = 2;
#endif
        }

	    sprintf(p.runtime_MultipolePrefix, "Multipoles");
	    p.runtime_NearFieldRadius = P.NearFieldRadius;
	    strcpy(p.runtime_TaylorDirectory,P.TaylorDirectory);
	    p.runtime_cpd = P.cpd;
	    p.runtime_order = P.order;
        p.rml = (P.order+1)*(P.order+1);
        p.CompressedMultipoleLengthXY = ((1+P.cpd)*(3+P.cpd))/8;
	    sprintf(p.runtime_TaylorPrefix, "Taylor");
        p.delete_multipoles_after_read = P.OverwriteState;
        
        int cml = ((P.order+1)*(P.order+2)*(P.order+3))/6;
        int nprocs = omp_get_max_threads();
        size_t cacherambytes = p.runtime_ConvolutionCacheSizeMB*(1024LL*1024LL);

        int blocksize = 0;
        for(blocksize=P.cpd*P.cpd;blocksize>=2;blocksize--) 
            if((P.cpd*P.cpd)%blocksize==0)
                if(nprocs*2.5*cml*blocksize*sizeof(Complex) < cacherambytes) break;
                    // 2.5 = 2 Complex (mcache,tcache) 1 double dcache
        p.blocksize = blocksize;
        
        uint64_t rambytes = p.runtime_MaxConvolutionRAMMB;
        rambytes *= 1024*1024;
        
        // If doing IO in parallel, need one block for reading, one for compute
        // The swizzle "block" can just be a single z-plane, but is usually twice the precision
#ifdef CONVIOTHREADED
        int n_alloc_block = 2;
#else
        int n_alloc_block = 1;
#endif
        uint64_t zslabbytes = p.rml*P.cpd*P.cpd*n_alloc_block*sizeof(MTCOMPLEX);
        zslabbytes += n_alloc_block*sizeof(double)*p.rml*p.CompressedMultipoleLengthXY;  // derivatives block
        STDLOG(0,"Each slab requires      %.2f MB\n",zslabbytes/1024/1024.);
        STDLOG(0,"You allow a maximum of  %.2f MB\n",rambytes/1024/1024.);
        uint64_t swizzlebytes = p.rml*P.cpd*P.cpd*sizeof(Complex);
        if(rambytes < zslabbytes) { 
            fprintf(stderr, "Each slab requires      %.2f MB\n", zslabbytes/1024/1024.);
            fprintf(stderr, "You allow a maximum of  %.2f MB\n", rambytes/1024/1024.);
            fprintf(stderr, "[ERROR] rambytes<zslabbytes\n");
            exit(1);
        }

        int zwidth = 0;
        for(zwidth=(P.cpd+1)/2;zwidth >= 1;zwidth--) {
            if( zwidth*zslabbytes + swizzlebytes < rambytes) break;
        }

        // If we're allocating more than one block, make sure we have at least one z-split
        if(n_alloc_block > 1)
            zwidth = min((P.cpd+1)/4, zwidth);

        p.zwidth = zwidth;
        STDLOG(0,"Resulting zwidth: %d \n",zwidth);
        
        for (int i = 0; i < MAX_IO_THREADS; i++)
            p.io_cores[i] = P.Conv_IOCores[i];
    
        STDLOG(2, "MTCOMPLEX (multipole/taylor) dtype width: %d\n", (int) sizeof(MTCOMPLEX));
        STDLOG(2, "DFLOAT (derivatives)         dtype width: %d\n", (int) sizeof(DFLOAT));

	    ConvolutionWallClock.Start();
	    STDLOG(1,"Starting Convolution\n");
	    OCC.Convolve(p);
	    STDLOG(1,"Convolution Complete\n");
	    ConvolutionWallClock.Stop();
	    TotalWallClock.Stop();
        
        OCC.CS.ConvolveWallClock = ConvolutionWallClock.Elapsed();
        
	    char timingfn[1050];
	    sprintf(timingfn,"%s/last.convtime",P.LogDirectory);
	    dumpstats(&OCC,timingfn);
	    stdlog.close();
        exit(0);
}
