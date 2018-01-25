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

#include "header.cpp"
#include "threevector.hh"
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
#include "kick.cpp"

#include "direct.h"
#include "direct.cpp"
#include "directdriver.cpp"

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

#include "Cosmology.cpp"
Cosmology *cosm;
#include "output_timeslice.cpp"
#include "LightCones.cpp"

#include "loadIC.cpp"
#include "lpt.cpp"

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

    // IO_Terminate();

    if(IL->length!=0) { IL->DumpParticles(); assert(IL->length==0); }
    
    char filename[1024];
    sprintf(filename,"%s/slabsize",P.WriteStateDirectory);
    Slab->write(filename);
    STDLOG(1,"Writing SlabSize file to %s\n", filename);

    MF->ComputeRedlack();  // NB when we terminate SlabMultipoles we write out these
    MF->WriteOutAuxiallaryVariables(P.WriteStateDirectory);

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

    delete MF;
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
    
    epilogue.Stop();
    STDLOG(1,"Leaving Epilogue()\n");
}
