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
#include <sys/time.h>
#include <sys/resource.h> 

#define COMPUTE_FOF_DENSITY
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
#include "PTimer.cc"

STimer FinishPartition;
STimer FinishSort;
STimer FinishCellIndex;
STimer FinishMerge;
PTimer FinishFreeSlabs;
STimer ComputeMultipoles;
STimer WriteMergeSlab;
STimer WriteMultipoleSlab;

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

STimer DriftMove, DriftRebin, DriftInsert;

STimer prologue;
STimer epilogue;
STimer WallClockDirect;
STimer SingleStepSetup;
STimer SingleStepTearDown;

uint64 naive_directinteractions = 0;    
//********************************************************************************

#include "file.h"
#include "grid.cpp"

#include "particlestruct.cpp"

#include "Parameters.cpp"
#include "statestructure.cpp"
State ReadState, WriteState;

// #include "ParticleCellInfoStructure.cpp"
// #include "maxcellsize.cpp"
#include "IC_classes.h"
#include "slabtypes.cpp"
SlabBuffer *LBW;

// Two quick functions so that the I/O routines don't need to know 
// about the LBW object. TODO: Move these to an io specific file
void IO_SetIOCompleted(int arena) { LBW->SetIOCompleted(arena); }
void IO_DeleteArena(int arena)    { LBW->DeAllocateArena(arena); }


#include "threadaffinity.h"

#ifdef IOTHREADED
#include "io_thread.cpp"
#else
#include "io_dio.cpp"
//#include "io_fopen.cpp"
#endif

#include "slabsize.cpp"
SlabSize *Slab;

#include "particles.cpp"
Particles *PP;

#include "dependency.cpp"

// Need this for both insert.cpp and timestep.cpp.
int FINISH_WAIT_RADIUS = 1;

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

#include "multipole_taylor.cpp"

#include "check.cpp"

// Forward-declare GFC
class GroupFindingControl;
GroupFindingControl *GFC;

#include "Cosmology.cpp"
Cosmology *cosm;
#include "lpt.cpp"
#include "output_timeslice.cpp"
#include "LightCones.cpp"

#include "loadIC.cpp"

#include "binning.cpp"
FLOAT * density; //!< Array to accumulate gridded densities in for low resolution inline power-spectra.

#include "groupfinding.cpp"
#include "microstep.cpp"

#include "timestep.cpp"
#include "reporting.cpp"

#include <fenv.h>

void load_slabsize(Parameters &P){
    char filename[1024];
    sprintf(filename,"%s/slabsize",P.ReadStateDirectory);
    Slab->read(filename);
    STDLOG(1,"Reading SlabSize file from %s\n", filename);
}


/*! \brief Initializes global objects
 *
 */
void Prologue(Parameters &P, bool ic) {
    omp_set_nested(true);

    STDLOG(1,"Entering Prologue()\n");
    prologue.Clear();
    prologue.Start();
    
    int cpd = P.cpd;
    int order = P.order;
    long long int np = P.np;
    assert(np>0);

    LBW = new SlabBuffer(cpd, order, cpd*MAXIDS, P.MAXRAMMB*1024*1024);
    PP = new Particles(cpd, LBW);
    STDLOG(1,"Initializing Multipoles()\n");
    MF  = new SlabMultipoles(order, cpd);

    STDLOG(1,"Setting up insert list\n");
    uint64 maxILsize = P.np+1;
    if (ic) {
        if (P.NumSlabsInsertListIC>0) maxILsize =(maxILsize* P.NumSlabsInsertListIC)/P.cpd+1;
    } else {
        if (P.NumSlabsInsertList>0) maxILsize   =(maxILsize* P.NumSlabsInsertList)/P.cpd+1;
    }
    IL = new InsertList(cpd, maxILsize);

    STDLOG(1,"Setting up IO\n");

    char logfn[1050];
    sprintf(logfn,"%s/lastrun.iolog", P.LogDirectory);
    io_ramdisk_global = P.RamDisk;
    STDLOG(0,"Setting RamDisk == %d\n", P.RamDisk);
    IO_Initialize(logfn);

    P.DensityKernelRad2 = 0.0;   // Don't compute densities

    if(!ic) {
            // ReadMaxCellSize(P);
        load_slabsize(P);
        TY  = new SlabTaylor(order,cpd);
            RL = new Redlack(cpd);

            SlabForceTime = new STimer[cpd];
            SlabForceLatency = new STimer[cpd];
            SlabFarForceTime = new STimer[cpd];

            RL->ReadInAuxiallaryVariables(P.ReadStateDirectory);
        
		// ForceOutputDebug outputs accelerations as soon as we compute them
		// i.e. before GroupFinding has a chance to rearrange them
        if(P.AllowGroupFinding && !P.ForceOutputDebug){
            GFC = new GroupFindingControl(P.FoFLinkingLength[0]/pow(P.np,1./3),
                                          P.FoFLinkingLength[1]/pow(P.np,1./3),
                                          P.FoFLinkingLength[2]/pow(P.np,1./3),
                                          P.cpd, PP->invcpd, P.GroupRadius, P.MinL1HaloNP, P.np);
	    #ifdef COMPUTE_FOF_DENSITY
		P.DensityKernelRad2 = GFC.linking_length;
		P.DensityKernelRad2 *= P.DensityKernelRad2; 
		// We use square radii
	    #endif
	}
    } else {
            TY = NULL;
            RL = NULL;
            JJ = NULL;
    }

    prologue.Stop();
    STDLOG(1,"Leaving Prologue()\n");
}

/*! \brief Tears down global objects
 *
 */
void Epilogue(Parameters &P, bool ic) {
    STDLOG(1,"Entering Epilogue()\n");
    epilogue.Clear();
    epilogue.Start();

    char timingfn[1050];
    sprintf(timingfn,"%s/lastrun.time", P.LogDirectory);
    FILE * timingfile = fopen(timingfn,"w");
    assertf(timingfile != NULL, "Couldn't open timing file \"%s\"\n", timingfile);
    ReportTimings(timingfile);
    fclose(timingfile);
    STDLOG(0,"Wrote Timing File to %s\n",timingfn);

    // IO_Terminate();

    if(IL->length!=0) { IL->DumpParticles(); assert(IL->length==0); }
    
    if(Slab != NULL){
        char filename[1024];
        sprintf(filename,"%s/slabsize",P.WriteStateDirectory);
        Slab->write(filename);
        STDLOG(1,"Writing SlabSize file to %s\n", filename);
    }

    // Some pipelines, like standalone_fof, don't use multipoles
    if(MF != NULL){
        MF->ComputeRedlack();  // NB when we terminate SlabMultipoles we write out these
        MF->WriteOutAuxiallaryVariables(P.WriteStateDirectory);
        delete MF;
    }

    if(ReadState.DoBinning){
            STDLOG(1,"Outputting Binned Density\n");
            char denfn[2048];
            sprintf(denfn,"%s/density",P.ReadStateDirectory);
            FILE * densout = fopen(denfn,"wb");
            fwrite(density,sizeof(FLOAT),P.PowerSpectrumN1d*P.PowerSpectrumN1d*P.PowerSpectrumN1d,densout);
            fclose(densout);
            delete density; density = 0;
    }
    
    if(WriteState.Do2LPTVelocityRereading)
        finish_2lpt_rereading();

    LBW->report();
    delete LBW;
    delete PP;
    delete IL;
    delete Slab;


    if(!ic) {
        if(0 and P.ForceOutputDebug){
            #ifdef CUDADIRECT
            STDLOG(1,"Direct Interactions: CPU (%llu) and GPU (%llu)\n",
                        JJ->DirectInteractions_CPU,JJ->TotalDirectInteractions_GPU);
            if(!(JJ->DirectInteractions_CPU == JJ->TotalDirectInteractions_GPU)){
                printf("Error:\n\tDirect Interactions differ between CPU (%llu) and GPU (%llu)\n",
                        JJ->DirectInteractions_CPU,JJ->TotalDirectInteractions_GPU);
                //assert(JJ->DirectInteractions_CPU == JJ->TotalDirectInteractions_GPU);
            }
            #endif
        }
        
            delete TY;
            delete RL;
            delete[] SlabForceLatency;
            delete[] SlabForceTime;
            delete[] SlabFarForceTime;
            delete JJ;
            delete GFC;
    }

    STDLOG(0,"MinCellSize = %d, MaxCellSize = %d\n", 
        WriteState.MinCellSize, WriteState.MaxCellSize);
    WriteState.RMS_Velocity = sqrt(WriteState.RMS_Velocity/P.np);
    STDLOG(0,"Rms |v| in simulation is %f.\n", WriteState.RMS_Velocity);
    STDLOG(0,"Maximum v_j in simulation is %f.\n", WriteState.MaxVelocity);
    STDLOG(0,"Maximum a_j in simulation is %f.\n", WriteState.MaxAcceleration);
    STDLOG(0,"Minimum cell Vrms/Amax in simulation is %f.\n", WriteState.MinVrmsOnAmax);
    
    // Report peak memory usage
    struct rusage rusage;
    assert(getrusage(RUSAGE_SELF, &rusage) == 0);
    STDLOG(0, "Peak resident memory usage was %.3g GB\n", (double) rusage.ru_maxrss / 1024 / 1024);
    
    epilogue.Stop();
    STDLOG(1,"Leaving Epilogue()\n");
}

std::vector<std::vector<int>> free_cores;  // list of free cores for each socket
void init_openmp(){
    // Tell singlestep to use the desired number of threads
    int max_threads = omp_get_max_threads();
    int ncores = omp_get_num_procs();
    int nthreads = P.OMP_NUM_THREADS > 0 ? P.OMP_NUM_THREADS : max_threads + P.OMP_NUM_THREADS;
    
    assertf(nthreads <= max_threads, "Trying to use more OMP threads (%d) than omp_get_max_threads() (%d)!  This will cause global objects that have already used omp_get_max_threads() to allocate thread workspace (like PTimer) to fail.\n",
        nthreads, max_threads);
    assertf(nthreads <= ncores, "Trying to use more threads (%d) than cores (%d).  This will probably be very slow.\n", nthreads, ncores);
    
    omp_set_num_threads(nthreads);
    STDLOG(1, "Initializing OpenMP with %d threads (system max is %d; P.OMP_NUM_THREADS is %d)\n", nthreads, max_threads, P.OMP_NUM_THREADS);
    
    // If threads are bound to cores via OMP_PROC_BIND,
    // then identify free cores for use by GPU and IO threads
    if(omp_get_proc_bind() == omp_proc_bind_false){
        //free_cores = NULL;  // signal that cores are not bound to threads
        STDLOG(1, "OMP_PROC_BIND = false; threads will not be bound to cores\n");
    }
    else{
        int core_assignments[nthreads];
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
        int main_thread_core = 0;
        set_core_affinity(main_thread_core);
        STDLOG(1, "Assigning main singlestep thread to core %d\n", main_thread_core);
    }
}

void setup_log(){
    std::setvbuf(stdout,(char *)_IONBF,0,0);
    std::setvbuf(stderr,(char *)_IONBF,0,0);

    stdlog_threshold_global = P.LogVerbosity;
    char logfn[1050];
    sprintf(logfn,"%s/lastrun.log", P.LogDirectory);
    stdlog.open(logfn);
    STDLOG_TIMESTAMP;
}

void check_read_state(int AllowIC, bool &MakeIC, double &da){
    // Check if ReadStateDirectory is accessible, or if we should 
    // build a new state from the IC file
    char rstatefn[1050];
    sprintf(rstatefn,"%s/state",P.ReadStateDirectory);

    if(access(rstatefn,0) ==-1){
        STDLOG(0,"Can't find ReadStateDirectory %s\n", P.ReadStateDirectory);
        if(AllowIC != 1){
            QUIT("Read State Directory ( %s ) is inaccessible and initial state creation is prohibited. Terminating.\n",P.ReadStateDirectory);

        } else{
            STDLOG(0,"Generating initial State from initial conditions\n");
        // We have to fill in a few items, just to bootstrap the rest of the code.
            ReadState.ScaleFactor = 1.0/(1+P.InitialRedshift);
            ReadState.FullStepNumber = -1;  
        // So that this number is the number of times forces have been computed.
        // The IC construction will yield a WriteState that is number 0,
        // so our first time computing forces will read from 0 and write to 1.
            da = 0;
            MakeIC = true;
        }
    } else {
    // We're doing a normal step
        CheckDirectoryExists(P.ReadStateDirectory);
        STDLOG(0,"Reading ReadState from %s\n",P.ReadStateDirectory);
        ReadState.read_from_file(P.ReadStateDirectory);
        ReadState.AssertStateLegal(P);
        MakeIC = false;
    // Handle some special cases
        if (P.ForceOutputDebug==1) {
            STDLOG(0,"ForceOutputDebug option invoked; setting time step to 0.\n");
            da = 0;
        }
    }
}
