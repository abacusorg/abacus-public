#include "mpi_header.cpp"

#include "header.cpp"
#include "threevector.hh"

#include "stdlog.cc"

#ifdef IOTHREADED
#define CONVIOTHREADED
#endif

#include "STimer.cc"
#include "PTimer.cc"

STimer TotalWallClock;
STimer Setup;
STimer ConvolutionWallClock;


#include "factorial.cpp"
#ifndef PARALLEL
#include "iolib.cpp"
#endif

int node_zstart = -1, node_zwidth = 1;
int first_slab_on_node = 0, first_slab_finished = -1, total_slabs_on_node = -1;


#include "Parameters.cpp"

#ifdef PARALLEL
	MPI_Datatype MPI_skewer;
#endif

/*
Each node will read multipole files for:
	-- all m
	-- all y
	-- range of z current BlockConvolve() is doing. 
    -- range of x each node is responsible for (set in singlestep). 

Each node will then send multipole files to its buddy nodes 
such that at the end of the day, each node has multipoles (and derivs) for:
	-- all m
	-- all y
	-- range of z NODE is responsible for (a subset of the z's the current BlockConvolve() is doing)
    -- all x.  

Each node will then perform convolution and get Taylors for the range of z the NODE is responsible for.

Each node will then receive Taylors from its buddy nodes 
such that at the end of the day, each node has Taylors for:
	-- range of z current BlockConvolve() is doing. 
    -- range of x each node is responsible for (set in singlestep).  
 
*/



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
    fprintf(fp,"Convolution parameters:  RamAllocated = %dMB CacheSizeMB = %4.1fMB nreal_cores=%d blocksize=%d zwidth=%d cpd=%d order=%d",
        (int) (OCC->CS.totalMemoryAllocated/(1<<20)), OCC->CS.runtime_ConvolutionCacheSizeMB, computecores, (int) OCC->CP.blocksize, (int) OCC->CP.zwidth, OCC->CP.runtime_cpd, OCC->CP.runtime_order);

#ifdef CONVIOTHREADED
    fprintf(fp, " niothread=%d", OCC->CP.niothreads);
#endif
    fprintf(fp,"\n\n");

    fprintf(fp,"\t ConvolutionWallClock:  %2.2e seconds \n", OCC->CS.ConvolveWallClock );
    fprintf(fp,"\t \t %50s : %1.2e seconds\n", "Array Swizzling", OCC->CS.ArraySwizzle );
    
#ifdef CONVIOTHREADED
    double e = OCC->CS.ReadDerivativesBytes/OCC->CS.ReadDerivatives/(1.0e+6);
    fprintf(fp,"\t \t %50s : %1.2e seconds --> rate was %4.0f MB/s\n", "ReadDiskDerivatives [per thread]", OCC->CS.ReadDerivatives, e );
    
    e = OCC->CS.ReadMultipolesBytes/OCC->CS.ReadMultipoles/(1.0e+6);
    fprintf(fp,"\t \t %50s : %1.2e seconds --> rate was %4.0f MB/s\n", "ReadDiskMultipoles [per thread]", OCC->CS.ReadMultipoles, e );
	
#ifdef PARALLEL
    e = OCC->CS.TransposeBufferingBytes/OCC->CS.TransposeBuffering/(1.0e+6);
    fprintf(fp,"\t \t %50s : %1.2e seconds --> rate was %4.0f MB/s\n", "Transpose Buffering", OCC->CS.TransposeBuffering, e );
	
    e = OCC->CS.TransposeAlltoAllvBytes/OCC->CS.TransposeAlltoAllv/(1.0e+6);
    fprintf(fp,"\t \t %50s : %1.2e seconds --> rate was %4.0f MB/s\n", "Transpose MPI AlltoAllv", OCC->CS.TransposeAlltoAllv, e );
#endif
    
    e = OCC->CS.WriteTaylorBytes/OCC->CS.WriteTaylor/(1.0e+6);
    fprintf(fp,"\t \t %50s : %1.2e seconds --> rate was %4.0f MB/s\n", "WriteDiskTaylor [per thread]", OCC->CS.WriteTaylor, e );
    
    fprintf(fp,"\t \t %50s : %1.2e seconds\n", "Waiting for IO thread", OCC->CS.WaitForIO);
#else
    double e = OCC->CS.ReadDerivativesBytes/OCC->CS.ReadDerivatives/(1.0e+6);
    fprintf(fp,"\t \t %50s : %1.2e seconds --> rate was %4.0f MB/s\n", "ReadDiskDerivatives", OCC->CS.ReadDerivatives, e );
    
    e = OCC->CS.ReadMultipolesBytes/OCC->CS.ReadMultipoles/(1.0e+6);
    fprintf(fp,"\t \t %50s : %1.2e seconds --> rate was %4.0f MB/s\n", "ReadDiskMultipoles", OCC->CS.ReadMultipoles, e );
	
    e = OCC->CS.WriteTaylorBytes/OCC->CS.WriteTaylor/(1.0e+6);
    fprintf(fp,"\t \t %50s : %1.2e seconds --> rate was %4.0f MB/s\n", "WriteDiskTaylor", OCC->CS.WriteTaylor, e );
#endif
    
    double Gops = ((double) OCC->CS.ops)/(1.0e+9);
    fprintf(fp,"\t \t %50s : %1.2e seconds for %5.3f billion double precision operations\n", "Convolution Arithmetic", OCC->CS.ConvolutionArithmetic, Gops );
    
    fprintf(fp,"\t \t %50s : %1.2e seconds\n", "Forward FFT Z Multipoles", OCC->CS.ForwardZFFTMultipoles );
    fprintf(fp,"\t \t %50s : %1.2e seconds\n", "Inverse FFT Z Taylor",         OCC->CS.InverseZFFTTaylor );
    
    fprintf(fp,"\t %50s : %1.2e seconds which is %d%% \n", "Unaccounted remaining wallclock time", discrepency, (int) (discrepency/OCC->CS.ConvolveWallClock*100) );

    double cae = OCC->CS.ConvolutionArithmetic;
    double farithp   = cae/OCC->CS.ConvolveWallClock*100;
    double ffftp     = (OCC->CS.ForwardZFFTMultipoles + OCC->CS.InverseZFFTTaylor)/OCC->CS.ConvolveWallClock*100;
    double fiop      = (OCC->CS.ReadDerivatives + OCC->CS.ReadMultipoles + OCC->CS.WriteTaylor )/OCC->CS.ConvolveWallClock*100;
#ifdef PARALLEL
	fiop += (OCC->CS.TransposeBuffering + OCC->CS.TransposeAlltoAllv)/OCC->CS.ConvolveWallClock*100;
#endif
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
        //STDLOG(3, "IO thread {:d} assigned to core {:d}\n", i, P.Conv_IOCores[i]);
    }
#endif
    
    assertf(nthreads <= max_threads,
        "Trying to use more OMP threads ({:d}) than omp_get_max_threads() ({:d})!  This will cause global objects that have already used omp_get_max_threads() to allocate thread workspace (like PTimer) to fail.\n",
        nthreads, max_threads);
    assertf(nthreads <= ncores,
        "Trying to use more threads ({:d}) than cores ({:d}).  This will probably be very slow.\n",
        nthreads, ncores);
    
    omp_set_num_threads(nthreads);
    STDLOG(0, "Initializing OpenMP with {:d} threads (system max is {:d}; P.Conv_OMP_NUM_THREADS is {:d})\n", nthreads, max_threads, P.Conv_OMP_NUM_THREADS);

    // If threads are bound to cores via OMP_PROC_BIND,
    // then identify free cores for use by IO thread
    if(omp_get_proc_bind() == omp_proc_bind_false){
        STDLOG(1, "OMP_PROC_BIND = false; threads will not be bound to cores\n");
    }
    else{
        int core_assignments[nthreads];
        #pragma omp parallel for schedule(static)
        for(int g = 0; g < nthreads; g++){
            assertf(g == omp_get_thread_num(), "OpenMP thread {:d} is executing wrong loop iteration ({:d})\n", omp_get_thread_num(), g);
            core_assignments[g] = sched_getcpu();
        }
        std::string core_log;
        core_log += "Thread->core assignments:";
        for(int g = 0; g < nthreads; g++)
            core_log += fmt::format(" {:d}->{:d}", g, core_assignments[g]);
        STDLOG(1, "{}\n", core_log);
        
        // Assign the main CPU thread to core 0 to avoid the IO thread during serial parts of the code
        /*int main_thread_core = P.Conv_IOCore != 0 ? 0 : 1;
        set_core_affinity(main_thread_core);
        STDLOG(1, "Assigning main convolution thread to core {:d}\n", main_thread_core);*/

        for(int g = 0; g < nthreads; g++)
            for(int h = 0; h < g; h++)
                assertf(core_assignments[g] != core_assignments[h], "Two OpenMP threads were assigned to the same core! This will probably be very slow. Check OMP_NUM_THREADS and OMP_PLACES?\n");
    }

}

int choose_zwidth(int Conv_zwidth, int cpd, ConvolutionParameters &CP){
    // Choose the convolution zwidth (number of planes to operate on at once)
    // Use the zwidth in the parameter file if given
    // Negative (default) means choose one based on available RAM
    // Zero means set zwidth to max

    uint64_t rambytes = CP.runtime_MaxConvolutionRAMMB;
    rambytes *= 1024*1024;
    
    // If doing IO in parallel, need one block for reading, one for compute
    // The swizzle "block" can just be a single z-plane, but is usually twice the precision
#ifdef CONVIOTHREADED
    int n_alloc_block = 2;
#else
    int n_alloc_block = 1;
#endif
    uint64_t zslabbytes = CP.rml*cpd*cpd*n_alloc_block*sizeof(MTCOMPLEX);
    zslabbytes += n_alloc_block*sizeof(DFLOAT)*CP.rml*CP.CompressedMultipoleLengthXY;  // derivatives block
    STDLOG(0,"Each slab requires      {:.2f} MB\n",zslabbytes/1024/1024.);
    STDLOG(0,"You allow a maximum of  {:.2f} MB\n",rambytes/1024/1024.);
    uint64_t swizzlebytes = CP.rml*cpd*cpd*sizeof(Complex);
    if(rambytes < zslabbytes) { 
        fprintf(stderr, "Each slab requires      %.2f MB\n", zslabbytes/1024/1024.);
        fprintf(stderr, "You allow a maximum of  %.2f MB\n", rambytes/1024/1024.);
        fprintf(stderr, "[ERROR] rambytes<zslabbytes\n");
        exit(1);
    }



    // If we are on the ramdisk, then we know the problem fits in memory! Just do the whole thing at once
    // If we aren't overwriting, there might be a small efficiency gain from smaller zwidth since reading requires a memcpy()
    // TODO: need to support ramdisk offsets if we want to support zwidth < max
    if(CP.is_ramdisk()){
		CP.z_slabs_per_node = (int)((cpd + 1)/2 / MPI_size) + 1; 
        STDLOG(0, "Forcing zwidth = full and z_slabs_per_node = {:d} since we are using ramdisk\n", CP.z_slabs_per_node);
        return (cpd + 1)/2;  // full width
    }
	
	
	//If we are doing a multi-node Convolve, set zwidth = 1. NAM TODO: May want to extend this to have multiple z per node later. 
#ifdef PARALLEL
	// STDLOG(0, "Forcing zwidth = {:d} in multi-node convolve, where each node does {:d} slab(s).\n", CP.z_slabs_per_node * MPI_size, CP.z_slabs_per_node);
// 	return CP.z_slabs_per_node * MPI_size;
	
	CP.z_slabs_per_node = (int)((cpd + 1)/2 / MPI_size) + 1; 
    STDLOG(0, "Forcing zwidth = full and z_slabs_per_node = {:d} since we are using ramdisk\n", CP.z_slabs_per_node);
    return (cpd + 1)/2;  // full width
#endif
	

    if(Conv_zwidth > 0){
        return Conv_zwidth;
    }
    else if(Conv_zwidth == 0){
        return (cpd + 1)/2;  // full width
    }
    else {
        int zwidth = 0;
        for(zwidth = (cpd+1)/2; zwidth >= 1;zwidth--) {
            if(zwidth*zslabbytes + swizzlebytes < rambytes)
                break;
        }

        // If we're allocating more than one block, make sure we have at least one z-split
        if(n_alloc_block > 1)
            zwidth = min((cpd+1)/4, zwidth);
        return zwidth;
    }

    return -1;
}


void InitializeParallel(int &size, int &rank, int rml, int CPD) {
    #ifdef PARALLEL
         // Start up MPI
         int ret;

         MPI_Init_thread(NULL, NULL, MPI_THREAD_SERIALIZED, &ret);  // TODO can we move this into the IO thread and change to funneled?
         assertf(ret>=MPI_THREAD_FUNNELED, "MPI_Init_thread() claims not to support MPI_THREAD_FUNNELED.\n");
         MPI_Comm_size(MPI_COMM_WORLD, &size);
         MPI_Comm_rank(MPI_COMM_WORLD, &rank);
         sprintf(NodeString,".%04d",rank);
		 
  		 MPI_Type_contiguous(rml*CPD, MPI_MTCOMPLEX, &MPI_skewer);
  		 MPI_Type_commit(&MPI_skewer);
		 
    #else
    #endif
    return;
}

void FinalizeParallel() {
    #ifdef PARALLEL
         // Finalize MPI
		 MPI_Type_free(&MPI_skewer);
         MPI_Finalize();
         STDLOG(0,"Calling MPI_Finalize()");
    #else
    #endif
}

int main(int argc, char ** argv){
	TotalWallClock.Start();
	Setup.Start();

	if (argc!=2) {
	       // Can't use assertf() or QUIT here: stdlog not yet defined!
	       fprintf(stderr, "Error: command line must have 1 parameter given, not %d.\nLegal usage: %s PARAM_FILE\n", argc-1, argv[0]);
	       assert(0==99);
	    }
		
	// Set up MPI

	    P.ReadParameters(argv[1],1);
		
		InitializeParallel(MPI_size, MPI_rank, (P.order+1)*(P.order+1), P.cpd);
		

	    // Setup the log
	    stdlog_threshold_global = P.LogVerbosity;
        // TODO: put log in step-numbered dir
        // Conv needs to learn the step number first
	    char logfn[1050];
	    sprintf(logfn,"%s/last%s.convlog", P.LogDirectory, NodeString);
	    stdlog.open(logfn);
	    STDLOG(1,"Read parameter file\n");
	    Setup.Stop();
        
        setup_openmp();

	    ConvolutionParameters CP(MPI_size);
	    CP.runtime_ConvolutionCacheSizeMB = P.ConvolutionCacheSizeMB;
	    CP.runtime_ConvolutionL1CacheSizeMB = P.ConvolutionL1CacheSizeMB;
        STDLOG(1, "Using L3 cache size {:f} MB\n", CP.runtime_ConvolutionCacheSizeMB);
        STDLOG(1, "Using L1 cache size {:f} MB\n", CP.runtime_ConvolutionL1CacheSizeMB);
	    CP.runtime_DerivativeExpansionRadius = P.DerivativeExpansionRadius;
	    strcpy(CP.runtime_DerivativesDirectory,P.DerivativesDirectory);
	    CP.runtime_DIOBufferSizeKB = 1LL<<11;
	    CP.runtime_AllowDIO = P.AllowDirectIO;
	    CP.runtime_MaxConvolutionRAMMB = P.MAXRAMMB;
	    strcpy(CP.runtime_MultipoleDirectory, P.MultipoleDirectory);

        CP.ProfilingMode = P.ProfilingMode;

	    sprintf(CP.runtime_MultipolePrefix, "Multipoles");
	    CP.runtime_NearFieldRadius = P.NearFieldRadius;
	    strcpy(CP.runtime_TaylorDirectory, P.TaylorDirectory);
	    CP.runtime_cpd = P.cpd;
	    CP.runtime_order = P.order;
        CP.rml = (P.order+1)*(P.order+1);
        CP.CompressedMultipoleLengthXY = ((1+P.cpd)*(3+P.cpd))/8;
	    sprintf(CP.runtime_TaylorPrefix, "Taylor");
        
        CP.StripeConvState = strcmp(P.Conv_IOMode, "stripe") == 0;
        CP.OverwriteConvState = strcmp(P.Conv_IOMode, "overwrite") == 0;

        // Determine number of IO threads
        CP.niothreads = 1;
        if(strcmp(P.MultipoleDirectory2,STRUNDEF) != 0){
            strcpy(CP.runtime_MultipoleDirectory2,P.MultipoleDirectory2);
#ifdef CONVIOTHREADED
            // Two IO threads if we were given two Multipole directories
            CP.niothreads = 2;
#endif
        }
        if(strcmp(P.TaylorDirectory2,STRUNDEF) != 0){
            strcpy(CP.runtime_TaylorDirectory2,P.TaylorDirectory2);
#ifdef CONVIOTHREADED
            CP.niothreads = 2;
#endif
        }
        
        int cml = ((P.order+1)*(P.order+2)*(P.order+3))/6;
        int nprocs = omp_get_max_threads();
        size_t cacherambytes = CP.runtime_ConvolutionCacheSizeMB*(1024LL*1024LL);
        size_t L1cacherambytes = CP.runtime_ConvolutionL1CacheSizeMB*(1024LL*1024LL);

        int blocksize = 1;
        for (int b=2; b<P.cpd*P.cpd; b++) {
            if (2.5*b*sizeof(Complex)>=L1cacherambytes) break;
                // Too much L1 memory: can't hold one example of D,M,T
            if (nprocs*2.5*cml*b*sizeof(Complex)>=cacherambytes) break;
                // Too much L3 memory, can't hold all D,M,T, so stop looking
            if ((P.cpd*P.cpd)%b == 0) blocksize = b;  // Could use this value
        }
            // 1.5 = 1 Complex (mcache) 1 double dcache
            // 2.5 = 2 Complex (mcache,tcache) 1 double dcache
            // 3.0 = 3 Complex (mcache,tcache,dcache)


        // TODO: This algorithm is silly.  We should start at the bottom!
        /*
        int blocksize = 0;
        for(blocksize=P.cpd*P.cpd;blocksize>=2;blocksize--) 
            if((P.cpd*P.cpd)%blocksize==0)
                if(nprocs*3.0*cml*blocksize*sizeof(Complex) < cacherambytes) break;
        */
                    
        CP.blocksize = blocksize;
		
#ifdef PARALLEL
		if (CP.niothreads != 1){

			QUIT("Multi-node (parallel, MPI) implementation does not support more than one IO thread! Crashing...\n");
		} 
		CP.z_slabs_per_node = 2; //may want to automate this choice or put an assert here to make sure we choose a value that uses all nodes. but the overkill version works too, where some nodes don't do anything at all (i.e. cpd = 33, num nodes = 5, zslabspernode = 8. )
#endif

        // Find out which slabs reside on which nodes
        // This will initialize first_slab_on_node and total_slabs_on_node even in the non-parallel version
        int read_all_nodes = 1;
        ReadNodeSlabs(read_all_nodes, CP.first_slabs_all, CP.total_slabs_all);
        
        CP.zwidth = choose_zwidth(P.Conv_zwidth, P.cpd, CP);
		
        for (int i = 0; i < MAX_IO_THREADS; i++)
            CP.io_cores[i] = P.Conv_IOCores[i];
    
        STDLOG(3, "MTCOMPLEX (multipole/taylor) dtype width: {:d}\n", (int) sizeof(MTCOMPLEX));
        STDLOG(3, "DFLOAT (derivatives)         dtype width: {:d}\n", (int) sizeof(DFLOAT));

        // Finished determining parameters; start the main work object
        OutofCoreConvolution OCC(CP);

	    ConvolutionWallClock.Start();
	    STDLOG(1,"Starting Convolution\n");
	    OCC.Convolve();
	    STDLOG(1,"Convolution Complete\n");
	    ConvolutionWallClock.Stop();
	    TotalWallClock.Stop();
        
        OCC.CS.ConvolveWallClock = ConvolutionWallClock.Elapsed();
        
	    char timingfn[1050];
	    sprintf(timingfn,"%s/last%s.convtime",P.LogDirectory,NodeString);
		
		
	    FinalizeParallel();  // This may be the last synchronization point?
		
		
	    dumpstats(&OCC,timingfn);
	    stdlog.close();

        // Delete the Taylors if this was profiling mode
        if(CP.ProfilingMode == 2){
            char cmd[1050];
            sprintf(cmd, "rm -f %s/Taylor_????", CP.runtime_TaylorDirectory);
            int ret = system(cmd);
            assertf(ret == 0, "Command \"%s\" failed\n", cmd);
            if(CP.StripeConvState){
                sprintf(cmd, "rm -f %s/Taylor_????", CP.runtime_TaylorDirectory2);
                ret = system(cmd);
                assertf(ret == 0, "Command \"%s\" failed\n", cmd);
            }
        }
        exit(0);
}
