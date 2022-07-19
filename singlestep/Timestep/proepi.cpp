/*! \brief Setup and teardown code for singlestep
 * proepi.cpp has several responsibilities:
 *  1.It is the root inclusion point for the entirety of Singlestep. All compile-time dependencies
 *      must be included here.
 *  2. Declaration of most global objects that are used between modules. Global objects used
 *      only within a module are typically defined in that module. Some global declarations are
 *      also non-conforming, and are dispersed throughout the code. TODO: This should be rectified
 *      eventually.
 *      Conventionally, objects should be declared like
 *              #include "SomeObjectDefinition.cpp"
 *              SomeObject * GlobalInstanceOfSomeObject;
 *  3. The prologue function, which initializes global objects like the Parameters class that
 *      must procede any useful work. All global objects declared in the beginning of the file
 *      should be initialized here, preferrably by a simple call to their constructor.
 *  4. The epilogue function, which handles teardown of most global objects. Conventionally, this
 *      should be done by deleting the object, with specific teardown code in the destructor
 *      for that module.
*/

#include <bitset>

#include <sys/time.h>
#include <sys/resource.h>

#include "mpi_header.cpp"
#include "header.cpp"
#include "threevector.hh"
#include "float3p1.cpp"    // This will define FLOAT3p1
        // Note that float3p1 and double3p1 are never defined.
        // Must call after FLOAT and FLOAT3 and float3 are defined

#include "StructureOfLists.cc"

#include "stdlog.cc"

//The following section declares a variety of global timers for several steps in the code
//TODO: This should probably go in some sort of reporting class to clean up this section.
#include "STimer.cc"
// #define PTIMER_DUMMY   // Uncommenting this will cause all PTimers to no-op and return Elapsed() = 1e-12 sec.
#include "PTimer.cc"
//#include "Limiter.cpp" 

#include "numa_for.h"

STimer FinishFreeSlabs; 
STimer FinishPartition;
STimer FinishSort;
STimer FinishCellIndex;
STimer FinishMerge;
STimer ComputeMultipoles;
STimer WriteMergeSlab;
STimer WriteMultipoleSlab;
STimer QueueMultipoleMPI;
STimer ParallelConvolveDestructor;

STimer OutputLightConeSearch;
STimer OutputTimeSlice;
STimer OutputLightCone;
STimer OutputBin;

STimer NearForceDirects;
STimer NearForceSetup;
STimer NearForceDirectsOnly;

STimer *SlabForceTime;
STimer *SlabForceLatency;
STimer *SlabFarForceTime;

STimer OutputGroup,GroupExecute;

STimer TaylorCompute;

STimer AddAccel;
STimer KickCellTimer;

STimer DriftMove, DriftRebin;

STimer prologue;
STimer epilogue;
STimer WallClockDirect;
STimer SingleStepSetup;
STimer SingleStepTearDown;
STimer IOFinish;

STimer SlabAccumFree;
STimer ReleaseFreeMemoryTime;

uint64 naive_directinteractions = 0;

//********************************************************************************

#include "file.h"
#include "grid.cpp"
grid *Grid;

#include "Parameters.cpp"
#include "statestructure.cpp"
State ReadState, WriteState;

#include "particlestruct.cpp"

// #include "ParticleCellInfoStructure.cpp"
// #include "maxcellsize.cpp"

#include "slabsize.cpp"
SlabSize *SS;

#include "IC_base.h"
#include "slabbuffer.cpp"
SlabBuffer *SB;

#include "slab_accum.cpp"
    // Code to establish templated slab-based storage of flexible size 
    // that is cell indexed and multi-threaded by pencil

#include "halostat.hh"

// Two quick functions so that the I/O routines don't need to know 
// about the SB object. TODO: Move these to an io specific file
void IO_SetIOCompleted(int arenatype, int arenaslab) {
	SB->SetIOCompleted(arenatype, arenaslab); }
void IO_DeleteArena(int arenatype, int arenaslab)    {
	SB->DeAllocate(arenatype, arenaslab); }


#include "threadaffinity.h"

// TODO: do we need to support non-threaded io? threaded io still has the blocking option
#ifdef IOTHREADED
#include "io_thread.cpp"
#else
#include "io_dio.cpp"
//#include "io_fopen.cpp"
#endif

#include "cellparticles.cpp"
CellParticles *CP;

#include "dependency.cpp"

// Need this for both insert.cpp and timestep.cpp.
int FINISH_WAIT_RADIUS = 1;

// Forward-declare GFC
class GroupFindingControl;
GroupFindingControl *GFC = NULL;

#include "multiappendlist.cpp"
#include "insert.cpp"
#include "drift.cpp"
#include "merge.cpp"

#include "direct.h"
#include "direct.cpp"
#include "directdriver.cpp"

#include "kick.cpp"

#include "basemultipoles.cpp"
#include "redlack.cpp"
#include "slabmultipoles.cpp"
#include "factorial.cpp"
#include "slabtaylor.cpp"
SlabTaylor *TY;
SlabMultipoles *MF;
Redlack *RL;

#ifdef PARALLEL
#include "slabmultipoles_mpi.cpp"
#include "slabtaylor_mpi.cpp"
#endif

#include "multipole_taylor.cpp"

#include "check.cpp"

#include "Cosmology.cpp"
Cosmology *cosm;
#include "lpt.cpp"
#include "loadIC.cpp"

#include "output_timeslice.cpp"
#include "LightCones.cpp"

#include "binning.cpp"
FLOAT * density; //!< Array to accumulate gridded densities in for low resolution inline power-spectra.

#include "groupfinding.cpp"
#include "microstep.cpp"
#include "output_field.cpp"    // Field particle subsample output

int first_slab_on_node, total_slabs_on_node, first_slab_finished;
int * first_slabs_all = NULL;
int * total_slabs_all = NULL;

	// The first read slab to be executed by this nodes,
	// as well as the total number and the first finished
	// In the single node code, this is simply 0 and CPD.
#include "node_slabs.cpp"

// FFTW Wisdom
char wisdom_file[1024];
int wisdom_exists;
void init_fftw();
void finish_fftw();

#include "timestep.cpp"
#include "parallel.cpp"
#include "reporting.cpp"

#include <fenv.h>

/*! \brief Initializes global objects
 *
 */
void Prologue(Parameters &P, int MakeIC, int NoForces) {
    STDLOG(1,"Entering Prologue()\n");
    STDLOG(2,"Size of accstruct is %d bytes\n", sizeof(accstruct));
    prologue.Clear();
    prologue.Start();

    int cpd = P.cpd;
    int order = P.order;
    long long int np = P.np;
    assert(np>0);

    init_fftw();  // wisdom import, etc, before SlabMultipoles or anything that uses FFTW

    // Look in ReadState to see what PosSlab files are available
    // TODO: Haven't implemented this yet
    first_slab_on_node = 0;
    // first_slab_on_node = ReadState.FullStepNumber; // A fun test
    total_slabs_on_node = cpd;
    // TODO: This fails for Spiral with first!=0 because the IC have
    // put all particles into input slab 0.

    // Call this to setup the Manifests
    SetupManifest(2*P.GroupRadius+1);

    Grid = new grid(cpd);
    SB = new SlabBuffer(cpd, order);
    CP = new CellParticles(cpd, SB);

    STDLOG(2,"Initializing Multipoles\n");
    if(MPI_size_z > 1){
        #ifdef PARALLEL
        MF = new SlabMultipolesMPI(order, cpd);
        #endif
    } else {
        MF = new SlabMultipolesLocal(order, cpd);
    }

    STDLOG(2,"Setting up insert list\n");
    uint64 maxILsize = P.np+1;
    // IC steps and LPT steps may need more IL slabs.  Their pipelines are not as long as full (i.e. group finding) steps
    if (MakeIC || LPTStepNumber() > 0){
        if (P.NumSlabsInsertListIC>0) maxILsize =(maxILsize* P.NumSlabsInsertListIC)*node_z_size/P.cpd/P.cpd+1;
    } else {
        if (P.NumSlabsInsertList>0) maxILsize   =(maxILsize* P.NumSlabsInsertList)*node_z_size/P.cpd/P.cpd+1;
    }
    IL = new InsertList(cpd, maxILsize);
    STDLOG(2,"Maximum insert list size = %l, size %l MB\n", maxILsize, maxILsize*sizeof(ilstruct)/1024/1024);

    STDLOG(2,"Setting up IO\n");

    char logfn[1050];
    sprintf(logfn,"%s/step%04d%s.iolog",
        WriteState.LogDirectory, WriteState.FullStepNumber, NodeString);
    allow_directio_global = P.AllowDirectIO;
    STDLOG(1,"Setting global AllowDirectIO = %d\n", P.AllowDirectIO);
    IO_Initialize(logfn);

    SS = new SlabSize(P.cpd, MPI_size_z, MPI_rank_z);
    ReadNodeSlabs();

    if(!MakeIC) {
            // ReadMaxCellSize(P);
        SS->load_from_params(P);
        
        if(!NoForces){
            if(MPI_size_z > 1){
                #ifdef PARALLEL
                TY = new SlabTaylorMPI(order,cpd);
                #endif
            } else {
                TY = new SlabTaylorLocal(order,cpd);
            }
        } else {
            TY = NULL;
        }
        
        RL = new Redlack(cpd);

        SlabForceTime = new STimer[cpd];
        SlabForceLatency = new STimer[cpd];
        SlabFarForceTime = new STimer[cpd];

		RL->ReadInAuxiallaryVariables(P.ReadStateDirectory);

        NFD = NoForces ? NULL : new NearFieldDriver(P.NearFieldRadius);
    } else {
        TY = NULL;
        RL = NULL;
        NFD = NULL;

        SlabForceTime = NULL;
        SlabForceLatency = NULL;
        SlabFarForceTime = NULL;
    }

    prologue.Stop();
    STDLOG(1,"Leaving Prologue()\n");
}

/*! \brief Tears down global objects
 *
 */
void Epilogue(Parameters &P, bool MakeIC) {
    STDLOG(1,"Entering Epilogue()\n");
    epilogue.Clear();
    epilogue.Start();

    // IO_Terminate();

	if(IL->length!=0) { IL->DumpParticles(); assert(IL->length==0); }

    if(SS != NULL){
        SS->store_from_params(P);
        WriteState.np_with_ghost_state = SS->total_new_size_with_ghost();
   	}


    if(MF != NULL){ // Some pipelines, like standalone_fof, don't use multipoles
        MF->GatherRedlack();    // For the parallel code, we have to coadd the inputs
        STDLOG(3, "globaldipole = (%g,%g,%g)\n", MF->globaldipole.x, MF->globaldipole.y, MF->globaldipole.z);

        if (MPI_rank==0) {
            MF->ComputeRedlack();
            MF->WriteOutAuxiallaryVariables(P.WriteStateDirectory);
        }
        delete MF;
        MF = NULL;
    }

    if(ReadState.DoBinning){
        STDLOG(1,"Outputting Binned Density\n");
        char denfn[1024];
        // TODO: Should this be going to ReadState or WriteState or Output?
        int ret = snprintf(denfn, 1024, "%s/density%s",P.ReadStateDirectory, NodeString);
        assert(ret >= 0 && ret < 1024);
        FILE * densout = fopen(denfn,"wb");
        fwrite(density,sizeof(FLOAT),P.PowerSpectrumN1d*P.PowerSpectrumN1d*P.PowerSpectrumN1d,densout);
        fclose(densout);
        delete density; density = 0;
    }

    SB->report_peak();
    delete SB;
    SB = NULL;
    STDLOG(2,"Deleted SB\n");
    delete CP;
    CP = NULL;
    delete IL;
    IL = NULL;
    delete SS;
    SS = NULL;
    delete Grid;
    Grid = NULL;
	
	FreeManifest();
    TeardownNeighborExchange();

    if(0 and P.ForceOutputDebug){
        #ifdef CUDADIRECT
        STDLOG(1,"Direct Interactions: CPU (%lu) and GPU (%lu)\n",
                    NFD->DirectInteractions_CPU,NFD->TotalDirectInteractions_GPU);
        if(!(NFD->DirectInteractions_CPU == NFD->TotalDirectInteractions_GPU)){
            printf("Error:\n\tDirect Interactions differ between CPU (%lu) and GPU (%lu)\n",
                    NFD->DirectInteractions_CPU,NFD->TotalDirectInteractions_GPU);
            //assert(NFD->DirectInteractions_CPU == NFD->TotalDirectInteractions_GPU);
        }
        #endif
    }
    
    delete TY;
    TY = NULL;
    STDLOG(2,"Deleted TY\n");
    delete RL;
    RL = NULL;
    delete[] SlabForceLatency;
    SlabForceLatency = NULL;
    delete[] SlabForceTime;
    SlabForceTime = NULL;
    delete[] SlabFarForceTime;
    SlabFarForceTime = NULL;
    delete GFC;
    GFC = NULL;
    if(NFD) {
        WriteState.DirectsPerParticle = (double)1.0e9*NFD->gdi_gpu/P.np;
        STDLOG(2,"About to delete NFD/kill the GPUs\n");
        delete NFD;
        NFD = NULL;
    } else {
        WriteState.DirectsPerParticle = 0;
    }

    finish_fftw();

    finish_numa_for();  // can't use NUMA_FOR after this

    // Report peak memory usage
    struct rusage rusage;
    assert(getrusage(RUSAGE_SELF, &rusage) == 0);
    STDLOG(0, "Peak resident memory usage was %.3g GB\n", (double) rusage.ru_maxrss / 1024 / 1024);
    
    epilogue.Stop();
    // This timing does not get written to the timing log, so it had better be small!
    STDLOG(1,"Leaving Epilogue(). Epilogue took %.2g sec.\n", epilogue.Elapsed());
}

void init_fftw(){
    // Import FFTW wisdom, before SlabMultipoles or anything that does FFT planning
    sprintf(wisdom_file, "%s/fftw_%d.wisdom", P.WorkingDirectory, P.cpd);
    wisdom_exists = fftw_import_wisdom_from_filename(wisdom_file);
    STDLOG(1, "Wisdom import from file \"%s\" returned %d (%s).\n", wisdom_file, wisdom_exists, wisdom_exists == 1 ? "success" : "failure");
}

void finish_fftw(){
    if(MPI_rank == 0){
        int ret = fftw_export_wisdom_to_filename(wisdom_file);
        STDLOG(1, "Wisdom export to file \"%s\" returned %d.\n", wisdom_file, ret);
    }
   fftw_cleanup();  // better not call this before exporting wisdom!
}

std::vector<std::vector<int>> free_cores;  // list of cores on each socket that are not assigned a thread (openmp, gpu, io, etc)
void init_openmp(){
    // First report the CPU affinity bitmask of the master thread; might be useful in diagnosing OMP_PLACES problems
    cpu_set_t mask;
    std::ostringstream affinitylog;

    assertf(pthread_getaffinity_np(pthread_self(), sizeof(cpu_set_t), &mask) == 0, "pthread_getaffinity_np failed\n");
    affinitylog << "Core affinity for the master thread: ";
    for (int i = 0; i < CPU_SETSIZE; i++) {
        if(CPU_ISSET(i, &mask))
            affinitylog << i << " ";
    }
    affinitylog << "\n";
    STDLOG(2, affinitylog.str().c_str());

    // Now tell OpenMP singlestep to use the desired number of threads
    int max_threads = omp_get_max_threads();
    int ncores = omp_get_num_procs();
    int nthreads = P.OMP_NUM_THREADS > 0 ? P.OMP_NUM_THREADS : max_threads + P.OMP_NUM_THREADS;

    // On summitdev (which has 160 thread slices), singlestep crashes.
    // I suspect this is because some of our stack-allocated arrays of size nthreads get too big
    // In practice, we will probably not use this many cores because there aren't nearly that many physical cores
    nthreads = min(128, nthreads);

    assertf(nthreads <= max_threads, "Trying to use more OMP threads (%d) than omp_get_max_threads() (%d)!  This will cause global objects that have already used omp_get_max_threads() to allocate thread workspace (like PTimer) to fail.\n",
        nthreads, max_threads);
    assertf(nthreads <= ncores, "Trying to use more threads (%d) than cores (%d).  This will probably be very slow.\n", nthreads, ncores);

    omp_set_num_threads(nthreads);
    STDLOG(1, "Initializing OpenMP with %d threads (system max is %d; P.OMP_NUM_THREADS is %d)\n", nthreads, max_threads, P.OMP_NUM_THREADS);

    // If threads are bound to cores via OMP_PROC_BIND,
    // then identify free cores for use by GPU and IO threads
    int core_assignments[nthreads];

    if(omp_get_proc_bind() == omp_proc_bind_false){
        //free_cores = NULL;  // signal that cores are not bound to threads
        STDLOG(1, "OMP_PROC_BIND = false; threads will not be bound to cores\n");

        // no need to init core_assignments here; init_numa_for will ignore it if omp_proc_bind_false
    }
    else{
        //bool is_core_free[ncores] = {true};
        #pragma omp parallel for schedule(static)
        for(int g = 0; g < nthreads; g++){
            assertf(g == omp_get_thread_num(), "OpenMP thread %d is executing wrong loop iteration (%d)\n", omp_get_thread_num(), g);
            core_assignments[g] = sched_getcpu();
        }
        std::ostringstream core_log;
        core_log << "Thread->core assignments:";
        for(int g = 0; g < nthreads; g++)
            core_log << " " << g << "->" << core_assignments[g];
        core_log << "\n";
        STDLOG(1, core_log.str().c_str());

        for(int g = 0; g < nthreads; g++)
            for(int h = 0; h < g; h++)
                assertf(core_assignments[g] != core_assignments[h], "Two OpenMP threads were assigned to the same core! This will probably be very slow. Check OMP_NUM_THREADS and OMP_PLACES?\n");

        // Assign the main CPU thread to core 0 to avoid the GPU/IO threads during serial parts of the code
        // Actually, this is unnecessary because the first OpenMP thread *is* the main CPU thread
        //int main_thread_core = 0;
        //set_core_affinity(main_thread_core);
        //STDLOG(1, "Assigning main singlestep thread to core %d\n", main_thread_core);
    }

    // Initialize the helper variables needed for "NUMA For"
    init_numa_for(nthreads, core_assignments);
}

void setup_log(){
    std::setvbuf(stdout,(char *)_IONBF,0,0);
    std::setvbuf(stderr,(char *)_IONBF,0,0);

    stdlog_threshold_global = P.LogVerbosity;
    char logfn[1050];
    sprintf(logfn,"%s/step%04d%s.log",
        WriteState.LogDirectory, WriteState.FullStepNumber, NodeString);
    stdlog.open(logfn);
    STDLOG_TIMESTAMP;
    STDLOG(0, "Log established with verbosity %d.\n", stdlog_threshold_global);
}


void load_read_state(int MakeIC){
    // Do an initial load of the read state,
    // basically just so we can get the step number and open the log with the right filename.
    // We are also going to cheat and set the WriteState step number here, but we'll verify
    // it below in InitWriteState.

    // Note we don't have STDLOG yet here; any reporting about ReadState should go in check_read_state()

    char rstatefn[1050];
    sprintf(rstatefn, "%s/state", P.ReadStateDirectory);

    if(MakeIC){
        // By this point, we should have cleaned up any old state directories
        if(access(rstatefn,0) != -1){
            QUIT("Read state file \"%s\" was found, but this is supposed to be an IC step. Terminating.\n", rstatefn);
        }

        // So that this number is the number of times forces have been computed.
        // The IC construction will yield a WriteState that is number 0,
        // so our first time computing forces will read from 0 and write to 1.
        ReadState.FullStepNumber = -1;

        // We have to fill in a few items, just to bootstrap the rest of the code.
        ReadState.ScaleFactor = 1.0/(1+P.InitialRedshift);

        // Probably not used, but this declares that the ReadState slabs have no ghosts
        ReadState.GhostRadius = 0;
    } else {
        // We're doing a normal step
        // Check that the read state file exists
        if(access(rstatefn,0) == -1){
            QUIT("Read state file \"%s\" is inaccessible and this is not an IC step. Terminating.\n", rstatefn);
        }

        ReadState.read_from_file(P.ReadStateDirectory);
    }

    // InitWriteState wants to use STDLOG, but we need WriteState.FullStepNumber to set up the log filename
    // So bootstrap that here.
    WriteState.FullStepNumber = ReadState.FullStepNumber + 1;
    snprintf(WriteState.LogDirectory, 1024, "%s/step%04d", P.LogDirectory, WriteState.FullStepNumber);
}

void check_read_state(const int MakeIC){
    // Check if ReadStateDirectory is accessible, or if we should
    // build a new state from the IC file

    if(MakeIC){
        STDLOG(0,"Generating initial State from initial conditions\n");
    } else {
        STDLOG(0,"Read ReadState from %s\n",P.ReadStateDirectory);
        ReadState.AssertStateLegal(P);
    }
}

// A few actions that we need to do before choosing the timestep
void InitWriteState(int MakeIC, const char *pipeline, const char *parfn){
    // Even though we do this in BuildWriteState, we want to have the step number
    // available when we choose the time step.
    assert(WriteState.FullStepNumber == ReadState.FullStepNumber+1);  // already did this in load_read_state()
    WriteState.LPTStepNumber = LPTStepNumber();

    // We generally want to do re-reading on the last LPT step
    WriteState.Do2LPTVelocityRereading = 0;
    if (WriteState.FullStepNumber == P.LagrangianPTOrder-1
        && (strcmp(P.ICFormat, "RVdoubleZel") == 0 || strcmp(P.ICFormat, "RVZel") == 0)){
        WriteState.Do2LPTVelocityRereading = 1;
    }

    // Decrease the softening length if we are doing a 2LPT step
    // This helps ensure that we are using the true 1/r^2 force
    /*if(LPTStepNumber()>0){
        WriteState.SofteningLengthNow = P.SofteningLength / 1e4;  // This might not be in the growing mode for this choice of softening, though
        STDLOG(0,"Reducing softening length from %f to %f because this is a 2LPT step.\n", P.SofteningLength, WriteState.SofteningLengthNow);
    }
    else{
        WriteState.SofteningLengthNow = P.SofteningLength;
    }*/

    // Is the softening fixed in proper coordinates?
    if(P.ProperSoftening){
        // TODO: use the ReadState or WriteState ScaleFactor?  We haven't chosen the timestep yet.
        WriteState.SofteningLengthNow = min(P.SofteningLength/ReadState.ScaleFactor, P.SofteningMax);
        STDLOG(1, "Adopting a comoving softening of %d, fixed in proper coordinates\n", WriteState.SofteningLengthNow);
    }
    else{
        WriteState.SofteningLengthNow = P.SofteningLength;
        STDLOG(1, "Adopting a comoving softening of %d, fixed in comoving coordinates\n", WriteState.SofteningLengthNow);
    }

    // Now scale the softening to match the minimum Plummer orbital period
#if defined DIRECTCUBICSPLINE
    strcpy(WriteState.SofteningType, "cubic_spline");
    WriteState.SofteningLengthNowInternal = WriteState.SofteningLengthNow * 1.10064;
#elif defined DIRECTSINGLESPLINE
    strcpy(WriteState.SofteningType, "single_spline");
    WriteState.SofteningLengthNowInternal = WriteState.SofteningLengthNow * 2.15517;
#elif defined DIRECTCUBICPLUMMER
    strcpy(WriteState.SofteningType, "cubic_plummer");
    WriteState.SofteningLengthNowInternal = WriteState.SofteningLengthNow * 1.;
#else
    strcpy(WriteState.SofteningType, "plummer");
    WriteState.SofteningLengthNowInternal = WriteState.SofteningLengthNow;
#endif

    if(strcmp(P.StateIOMode, "overwrite") == 0){
        WriteState.OverwriteState = 1;
        STDLOG(1, "StateIOMode = \"overwrite\"; write state will overwrite read state\n");
    }
    if(strcmp(P.StateIOMode, "stripe") == 0){
        WriteState.StripeState = 1;
        assertf(0, "State striping currently not implemented\n");
    }

    if(strcmp(P.Conv_IOMode, "stripe") == 0){
        WriteState.StripeConvState = 1;
        STDLOG(1,"Striping multipoles and taylors\n");
    }
    else if(strcmp(P.Conv_IOMode, "overwrite") == 0){
        WriteState.OverwriteConvState = 1;
        STDLOG(1,"Overwriting multipoles and taylors\n");
    }

    strcpy(WriteState.Pipeline, pipeline);
    strcpy(WriteState.ParameterFileName, parfn);

    // Check if WriteStateDirectory/state exists, and fail if it does
    char wstatefn[1050];
    int ret = snprintf(wstatefn, 1050, "%s/state", P.WriteStateDirectory);
    assert(ret >= 0 && ret < 1050);
    if(access(wstatefn,0) !=-1 && !WriteState.OverwriteState)
        QUIT("WriteState \"%s\" exists and would be overwritten. Please move or delete it to continue.\n", wstatefn);
}


void InitKernelDensity(){
    /* N.B. we are opting to store these concepts in the WriteState,
    even though they relate to the ReadState positions and cosmology.
    We may revisit this choice, but at least it allows us to easily
    log these values as part of the write state file, and avoids the
    unseemly practice of modifying the ReadState.
    */

    #ifdef COMPUTE_FOF_DENSITY
    #ifdef CUDADIRECT   // For now, the CPU doesn't compute FOF densities, so signal this by leaving Rad2=0.
    if (P.DensityKernelRad==0) {
        // Default to the L0 linking length
        if (GFC != NULL) WriteState.DensityKernelRad2 = GFC->linking_length;
        else WriteState.DensityKernelRad2 = P.FoFLinkingLength[0]/pow(P.np,1./3); 
        WriteState.DensityKernelRad2 *= WriteState.DensityKernelRad2*(1.0+1.0e-5);

        // We use square radii.  The radius is padded just a little
        // bit so we don't risk underflow with 1 particle at r=b
        // in comparison to the self-count.
        WriteState.L0DensityThreshold = 0.0;
        // Use this as a signal to use DensityKernalRad2 (in code units,
        // not cosmic units) as the threshold,
        // which means that a particle is L0 eligible if there is any
        // non-self particle within the L0 linking length
    } else {
        WriteState.DensityKernelRad2 = P.DensityKernelRad/pow(P.np,1./3);
        WriteState.DensityKernelRad2 *= WriteState.DensityKernelRad2;
        WriteState.L0DensityThreshold = P.L0DensityThreshold;
    }

    WriteState.FOFunitdensity    = P.np*4.0*M_PI*2.0/15.0*pow(WriteState.DensityKernelRad2,2.5)+1e-30;
    WriteState.invFOFunitdensity = 1.0/WriteState.FOFunitdensity;
    #endif
    #endif
}

template<class C, typename T>
bool contains(C *a, C *b, T e) { return std::find(a, b, e) != b; };


void InitGroupFinding(int MakeIC){
    /*
    Request output of L1 groups and halo/field subsamples if:
    - this redshift is a member of L1OutputRedshift; or
    - this redshift is a member of TimeSliceRedshifts_Subsample
    - this redshift is a member of TimeSliceRedshifts.

    If L1OutputRedshifts is not set, then we instead fall back to L1Output_dlna.
    In that case, we check if we are crossing a L1Output_dlna checkpoint by going from ReadState to WriteState.

    We may not end up outputting groups if group finding is not enabled.  This is signaled by GFC = NULL.

    We need to enable group finding if:
    - We are doing microstepping
    - We are outputting groups
    But we can't enable it if:
    - AllowGroupFinding is disabled
    - ForceOutputDebug is enabled
    - This is an IC or 2LPT step
    ForceOutputDebug outputs accelerations as soon as we compute them
    i.e. before GroupFinding has a chance to rearrange them

    2D:
    timestep (1): shorten timestep
    timestep (2): do a da=0 timestep, with forces, to half-kick velocities to synchronicity (also send ghosts for group finding)
    timestep (3): do a da=0, no-force, group-finding timestep (send ghosts)
    timestep (4): resume timestepping with forces

    // TODO: is there anyway to recycle ghosts? the only bits that need updating are groupfinding tags
    // TODO: is it okay to locally tag particles whose groups we de-dup to another node?
    */

    // Decide if the next step will do group finding
    WriteState.DoGroupFindingOutput = 0;

    for(int i = 0; i < P.nTimeSliceL1; i++){
        double L1z = P.L1OutputRedshifts[i];
        double L1a = 1.0/(1+L1z);

        if(contains(P.TimeSliceRedshifts, P.TimeSliceRedshifts + P.nTimeSlice, L1z) || 
            contains(P.TimeSliceRedshifts_Subsample, P.TimeSliceRedshifts_Subsample + P.nTimeSliceSubsample, L1z)
         ) continue;


        if(ReadState.DoGroupFindingOutput != (i+2) &&
            ReadState.ScaleFactor < L1a && WriteState.ScaleFactor >= L1a){

            // We don't shorten our timestep to land exactly on a merger tree redshift.
            // Is the next time step going to take us further from the target z?
            if(fabs(ReadState.ScaleFactor - L1a) > fabs(WriteState.ScaleFactor - L1a)){
                WriteState.DoGroupFindingOutput = i+2;
                STDLOG(0,"Group finding on the write state requested by L1OutputRedshifts[%d]\n", i);
            } else {
                ReadState.DoGroupFindingOutput = i+2;
                STDLOG(0,"Group finding on the read state requested by L1OutputRedshifts[%d]\n", i);
            }
            
            break;
        }
    }

    if(WriteState.DoTimeSliceOutput && !WriteState.DoSubsampleOutput) WriteState.DoSubsampleOutput = 1;

    if(WriteState.DoSubsampleOutput && !WriteState.DoGroupFindingOutput) WriteState.DoGroupFindingOutput = 1;

    if(P.L1Output_dlna >= 0 && !WriteState.DoGroupFindingOutput &&
        !(ReadState.DoGroupFindingOutput && WriteState.DeltaScaleFactor == 0)){
        WriteState.DoGroupFindingOutput = P.L1Output_dlna == 0 ||
                    ( log(WriteState.ScaleFactor) - log(ReadState.ScaleFactor) >= P.L1Output_dlna ) ||
                    ( fmod(log(WriteState.ScaleFactor), P.L1Output_dlna) < fmod(log(ReadState.ScaleFactor), P.L1Output_dlna) );
        STDLOG(0,"Group finding at this redshift requested by L1Output_dlna\n");
    }
    
    if(P.OutputAllHaloParticles && !WriteState.DoSubsampleOutput){
        WriteState.DoSubsampleOutput = 1;
        STDLOG(0,"OutputAllHaloParticles = 1; forcing subsample A = 100%% and B = 0%%\n");
    }

    if(MPI_size_z > 1 && ReadState.DoGroupFindingOutput && !ReadState.VelIsSynchronous){
        WriteState.DoGroupFindingOutput = ReadState.DoGroupFindingOutput;
        WriteState.DoSubsampleOutput = ReadState.DoSubsampleOutput;
        ReadState.DoGroupFindingOutput = 0;
        ReadState.DoSubsampleOutput = 0;
        STDLOG(0, "Deferring 2D group finding to next step when vel is synchronous\n");
    }

    if(MPI_size_z > 1 && ReadState.DoGroupFindingOutput && !ReadState.HaveAuxDensity){
        WriteState.DoGroupFindingOutput = ReadState.DoGroupFindingOutput;
        WriteState.DoSubsampleOutput = ReadState.DoSubsampleOutput;
        ReadState.DoGroupFindingOutput = 0;
        ReadState.DoSubsampleOutput = 0;
        STDLOG(0, "Deferring 2D group finding to next step when we have aux densities\n");
    }

    if(!P.AllowGroupFinding){
        ReadState.DoGroupFindingOutput = 0;
        WriteState.DoGroupFindingOutput = 0;
        STDLOG(0, "P.AllowGroupFinding does not allow group finding.\n");
    }

    // Are we planning to store densities in aux?
    WriteState.HaveAuxDensity = MPI_size_z > 1 &&
        WriteState.DoGroupFindingOutput && WriteState.VelIsSynchronous && !ReadState.HaveAuxDensity && !MakeIC;

    // Or we might be eligible to reuse the densities if nothing is going to move
    if (ReadState.HaveAuxDensity && std::abs(WriteState.DeltaEtaDrift) < 1e-12) WriteState.HaveAuxDensity = 1;

    // Set up the density kernel
    // This used to inform the group finding, but also might be output as part of lightcones even on non-group-finding steps

    WriteState.DensityKernelRad2 = 0.0;   // Don't compute densities.  This is set in InitKernelDensity if we decide to compute densities
    WriteState.L0DensityThreshold = 0.0;  // Only has any effect if doing group finding

    // done planning for WriteState. Does ReadState tell us to find groups?

    // Can we enable group finding?
    if((P.MicrostepTimeStep > 0 || ReadState.DoGroupFindingOutput) &&
        !(!P.AllowGroupFinding || P.ForceOutputDebug || MakeIC || LPTStepNumber())){
        STDLOG(1, "Setting up group finding\n");

        STDLOG(2, "Group finding: %d, subsample output: %d, timeslice output: %d.\n",
            ReadState.DoGroupFindingOutput, ReadState.DoSubsampleOutput, ReadState.DoTimeSliceOutput);

        GFC = new GroupFindingControl(P.FoFLinkingLength[0]/pow(P.np,1./3),
                    #ifdef SPHERICAL_OVERDENSITY
                    P.SODensity[0], P.SODensity[1],  //by this point, the SODensity and L0DensityThresholds have been rescaled with redshift, in PlanOutput.
                    #else
                    P.FoFLinkingLength[1]/pow(P.np,1./3),
                    P.FoFLinkingLength[2]/pow(P.np,1./3),
                    #endif
                    P.cpd, node_z_start_ghost, node_z_size_with_ghost,
                    P.GroupRadius, P.MinL1HaloNP, P.np);

        #ifdef SPHERICAL_OVERDENSITY
        WriteState.SODensityL1 = P.SODensity[0];
        WriteState.SODensityL2 = P.SODensity[1];
        #endif

        InitKernelDensity();

        if (WriteState.L0DensityThreshold==0) {
            STDLOG(1,"Passing L0DensityThreshold = 0 to signal to use anything with a neighbor\n");
        } else {
            STDLOG(1,"Using L0DensityThreshold = %f\n", WriteState.L0DensityThreshold);
        }

        if(MPI_size_z > 1){
            assertf(GHOST_RADIUS >= GROUP_RADIUS,
                    "GHOST_RADIUS=%d not big enough for GROUP_RADIUS=%d\n",
                    GHOST_RADIUS, GROUP_RADIUS);
        }
    } else{
        GFC = NULL;  // be explicit
        ReadState.DoGroupFindingOutput = 0;
        ReadState.DoSubsampleOutput = 0;  // We currently do not support subsample outputs without group finding
        STDLOG(1, "Group finding not enabled for this step.\n");

        // We aren't doing group finding, but we may be doing output, so init the kernel density anyway
        InitKernelDensity();
    }

    STDLOG(1,"Using DensityKernelRad2 = %f (%f of interparticle)\n", WriteState.DensityKernelRad2, sqrt(WriteState.DensityKernelRad2)*pow(P.np,1./3.));

}

// Check whether "d" is actually a global directory, and thus not eligible for deletion
int IsTrueLocalDirectory(const char* d){
    if(samefile(d, P.WorkingDirectory) ||
        samefile(d, P.ReadStateDirectory) ||
        samefile(d, P.WriteStateDirectory) ||
        samefile(d, P.InitialConditionsDirectory)
        ) {
        return 0;
    }

    return 1;
}


// This function creates all the "local" directories for singlestep;
// i.e. all the directories that the Python code didn't create.
// In the parallel code, that means this function is responsible for creating all node-local directories
// This also deletes existing state directories if MakeIC is invoked
void SetupLocalDirectories(const int MakeIC){
    /* TODO: probably deprecated
    // Resume from a backed-up state
    if(!MakeIC){
        // If BackupDirectory exists and LocalReadStateDirectory does not
        // then we should set [Local]ReadStateDirectory to BackupDirectory
        // Need to set both because can't risk an inconsistent Read and LocalRead state!

        if(CheckFileExists(P.BackupDirectory) == 2
            && CheckFileExists(P.LocalReadStateDirectory) != 2){
            fprintf(stderr, "Local read state dir \"%s\" not found; reading from BackupDirectory \"%s\"\n",
                P.LocalReadStateDirectory, P.BackupDirectory);

            sprintf(P.LocalReadStateDirectory, "%s/read", P.BackupDirectory);
            sprintf(P.ReadStateDirectory, "%s/read", P.BackupDirectory);
        }
    }*/

    // TODO: might want to delete old derivatives directory here,
    // but the risk of accidentally deleting the global derivatives is very high
    char *dirs[] = {P.LocalWorkingDirectory,
                    P.LocalReadStateDirectory,
                    P.LocalWriteStateDirectory,
                    P.TaylorDirectory,
                    P.MultipoleDirectory,
                    P.TaylorDirectory2,
                    P.MultipoleDirectory2
                };

    for(int i = 0; i < sizeof(dirs)/sizeof(char*); i++){
        const char *d = dirs[i];

        if(strcmp(d, STRUNDEF) != 0 && strlen(d) > 0){
            // The following functions don't care if the directory already exists or not
            if(MakeIC && IsTrueLocalDirectory(d)){
                RemoveDirectories(d);
                STDLOG(1, "Removed directory \"%s\"\n", d);
            }

            if(d == P.LocalWriteStateDirectory && strcmp(P.StateIOMode, "overwrite") == 0) {
                CreateSymlink(P.LocalReadStateDirectory, P.LocalWriteStateDirectory);
            } else {
                int res = CreateDirectories(d);
                assertf(res == 0, "Creating directory \"%s\" failed for reason %s!\n", d, strerror(errno));
                STDLOG(1, "Created directory \"%s\"\n", d);
            }
        }
    }
}

// Move the state directories at the end of a timestep.
// Deletes "read", and moves "write" to "read"
// singlestep doesn't know about "past" directory.
void MoveLocalDirectories(){
    if(WriteState.OverwriteState){
        // Nothing to do!
        return;
    }

    if(IsTrueLocalDirectory(P.LocalReadStateDirectory)){
        STDLOG(1, "Removing read directory\n");
        int res = RemoveDirectories(P.LocalReadStateDirectory);
        assertf(res == 0, "Failed to remove read directory!\n");
    }

    if(IsTrueLocalDirectory(P.LocalWriteStateDirectory)){
        STDLOG(1, "Moving write directory to read\n");
        int res = rename(P.LocalWriteStateDirectory, P.LocalReadStateDirectory);
        assertf(res == 0, "Failed to rename write to read!\n");
    }
}


void FinalizeWriteState() {
    WriteNodeSlabs();  // We do this here because it will need a MPI Barrier

    #ifdef PARALLEL
        STDLOG(1,"Node MinCellSize = %d, MaxCellSize = %d\n",
            WriteState.MinCellSize, WriteState.MaxCellSize);
        STDLOG(1,"Maximum v_j in node is %f.\n", WriteState.MaxVelocity);
        STDLOG(1,"Maximum a_j in node is %f.\n", WriteState.MaxAcceleration);
        STDLOG(1,"Minimum cell Vrms/Amax in node is %f.\n", WriteState.MinVrmsOnAmax);
        STDLOG(1,"Unnormalized node RMS_Velocity = %f.\n", WriteState.RMS_Velocity);
        STDLOG(1,"Unnormalized node StdDevCellSize = %f.\n", WriteState.StdDevCellSize);
        STDLOG(1,"Maximum group diameter in node is %d.\n", WriteState.MaxGroupDiameter); 
        STDLOG(1,"Maximum L0 group size in node is %d.\n", WriteState.MaxL0GroupSize); 
    
        // If we're running in parallel, then we want to gather some
        // state statistics across the nodes.  We start by writing the
        // original state file to the local disk.
        // TODO: Do we really want to do this?  Maybe just echo the stats to the log?
        // WriteState.write_to_file(P.WriteStateDirectory, NodeString);

// #define MPI_REDUCE_IN_PLACE(vec,len,type,op) MPI_Reduce(MPI_rank!=0?(vec):MPI_IN_PLACE, vec, len, type, op, 0, comm_global)
        // Now we need to do MPI reductions for stats
        // These stats are all in double precision (or int)
        // Maximize MaxAcceleration
        MPI_REDUCE_TO_ZERO(&WriteState.MaxAcceleration, 1, MPI_DOUBLE, MPI_MAX);
        // Maximize MaxVelocity
        MPI_REDUCE_TO_ZERO(&WriteState.MaxVelocity, 1, MPI_DOUBLE, MPI_MAX);
        // Minimize MinVrmsOnAmax
        MPI_REDUCE_TO_ZERO(&WriteState.MinVrmsOnAmax, 1, MPI_DOUBLE, MPI_MIN);
        // Minimize MinCellSize
        MPI_REDUCE_TO_ZERO(&WriteState.MinCellSize, 1, MPI_INT, MPI_MIN);
        // Maximize MaxCellSize
        MPI_REDUCE_TO_ZERO(&WriteState.MaxCellSize, 1, MPI_INT, MPI_MAX);
        // sqrt(Sum(SQR of RMS_Velocity))
        MPI_REDUCE_TO_ZERO(&WriteState.RMS_Velocity, 1, MPI_DOUBLE, MPI_SUM);
        // sqrt(Sum(SQR of StdDevCellSize))
        MPI_REDUCE_TO_ZERO(&WriteState.StdDevCellSize, 1, MPI_DOUBLE, MPI_SUM);
        // Maximize MaxGroupDiameter
        MPI_REDUCE_TO_ZERO(&WriteState.MaxGroupDiameter, 1, MPI_INT, MPI_MAX);
        // Maximize MaxL0GroupSize
        MPI_REDUCE_TO_ZERO(&WriteState.MaxL0GroupSize, 1, MPI_INT, MPI_MAX);
        // Sum WriteState.DirectsPerParticle
        MPI_REDUCE_TO_ZERO(&WriteState.DirectsPerParticle, 1, MPI_DOUBLE, MPI_SUM);
        // sum subsample counts
        MPI_REDUCE_TO_ZERO(&WriteState.np_subA_state, 1, MPI_INT, MPI_SUM);
        MPI_REDUCE_TO_ZERO(&WriteState.np_subB_state, 1, MPI_INT, MPI_SUM);
        MPI_REDUCE_TO_ZERO(&WriteState.np_lightcone, 1, MPI_INT64_T, MPI_SUM);
// #undef MPI_REDUCE_IN_PLACE

        // Note that we're not summing up any timing or group finding reporting;
        // these just go in the logs
    #endif

    WriteState.StdDevCellSize = sqrt(WriteState.StdDevCellSize);
        // This is the standard deviation of the fractional overdensity in cells.
        // But for the parallel code: this has been divided by CPD^3, not the number of cells on the node

    STDLOG(0,"MinCellSize = %d, MaxCellSize = %d, RMS Fractional Overdensity = %g\n", 
        WriteState.MinCellSize, WriteState.MaxCellSize, WriteState.StdDevCellSize);

    double code_to_kms = WriteState.VelZSpace_to_kms/WriteState.VelZSpace_to_Canonical;
    WriteState.RMS_Velocity = sqrt(WriteState.RMS_Velocity/P.np);

    STDLOG(0,"Rms |v| in simulation is %f km/s.\n", WriteState.RMS_Velocity * code_to_kms);
    STDLOG(0,"Maximum v_j in simulation is %f km/s.\n", WriteState.MaxVelocity * code_to_kms);
    STDLOG(0,"Maximum a_j in simulation is %f code units.\n", WriteState.MaxAcceleration);
    STDLOG(0,"Minimum cell Vrms/Amax in simulation is %f code units.\n", WriteState.MinVrmsOnAmax);
    STDLOG(0,"Maximum group diameter in simulation is %d.\n", WriteState.MaxGroupDiameter); 
    STDLOG(0,"Maximum L0 group size in simulation is %d.\n", WriteState.MaxL0GroupSize); 
    STDLOG(0,"Mean Directs per particle in simulation is %d.\n", WriteState.DirectsPerParticle); 

    //NAM TODO put an intelligent assertf here. 
    // the things we'd really like to check are: 
    //    1. is the # of subsampled particles in A/B constant for each snapshot z? 
    //    2. is the # of subsampled particles what we expect from the subsample fractions? 


    // if (ReadState.DoSubsampleOutput){ 
    //     assertf(WriteState.np_subA_state == (int) ( P.ParticleSubsampleA * P.np), "Subsample A contains %d particles, expected %d.\n", WriteState.np_subA_state, (int) (P.ParticleSubsampleA * P.np) ); 
    //     assertf(WriteState.np_subB_state == (int) ( P.ParticleSubsampleB * P.np), "Subsample A contains %d particles, expected %d.\n", WriteState.np_subB_state, (int) (P.ParticleSubsampleB * P.np) ); 
    // }


    // If we're writing lightcones, we only want the header to be written once
    // But a-priori there's no good way to know which nodes/slabs will have LC particles,
    // Now that we've done the reduction, we know if any LC particles were written, thus rank 0 can write the header
    if(WriteState.np_lightcone && MPI_rank == 0){
        char dir[32];
        sprintf(dir, "Step%04d", ReadState.FullStepNumber);
        
        std::string headerfn = "";
        headerfn = headerfn + P.LightConeDirectory + "/" + dir + "/header";
        LightCone::WriteHeaderFile(headerfn.c_str());
    }

    return;
}
