#include "header.cpp"
#include "threevector.hh"

#include "stdlog.cc"

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


STimer TaylorCompute;

STimer DriftMoveRebin, DriftInsert;
STimer LPTDriftICReRead;
PTimer DriftMove, DriftRebin;

STimer prologue;
STimer epilogue;
STimer WallClockDirect;
STimer SingleStepSetup;
STimer SingleStepTearDown;

STimer TotalRead;
STimer TotalWrite;
uint64 total_disk_read = 0;
uint64 total_disk_write = 0;



uint64 naive_directinteractions = 0;    

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
// about the LBW object.
void IO_SetIOCompleted(int arena) { LBW->SetIOCompleted(arena); }
void IO_DeleteArena(int arena)    { LBW->DeAllocateArena(arena); }


#include "threadaffinity.h"

#ifdef IOTHREADED
#include "io_thread.cpp"
#else
//#include "io_dio.cpp"
#include "io_fopen.cpp"
#endif

#include "particles.cpp"
Particles *PP;

#include "slabsize.cpp"
SlabSize *Slab;

#include "dependency.cpp"


#include "direct.h"
#include "direct.cpp"
#include "directdriver.cpp"
NearFieldDriver *JJ;

#include "insert.cpp"
#include "drift.cpp"
#include "merge.cpp"
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

#include"Cosmology.cpp"
Cosmology *cosm;
#include "output_timeslice.cpp"
#include "LightCones.cpp"

// Bookkeeping for 2LPT velocity re-reading
typedef struct {
    uint64 n_part;
    //uint64 n_read;
} VelIC;
VelIC* vel_ics;  // Array of VelIC structs

//FIXME:These will be slow for any large problem and should be refactored

#include "loadIC.cpp"
#include "lpt.cpp"

#include "binning.cpp"
FLOAT * density;

#include "timestep.cpp"
#include "reporting.cpp"

#include <fenv.h>



void Prologue(Parameters &P, bool ic) {
    omp_set_nested(true);    

    STDLOG(1,"Entering Prologue()\n");
    prologue.Clear();
    prologue.Start();
    
    int cpd = P.cpd;
    int order = P.order;
    long long int np = P.np;
    assert(np>0);


    // TODO: Need a better number to set the maximum size of the linear buffer
    LBW = new SlabBuffer(cpd, order, cpd*MAXIDS, 64e9);
    PP = new Particles(cpd, LBW);
    Slab = new SlabSize(cpd);
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
    	char filename[1024];
    	sprintf(filename,"%s/slabsize",P.ReadStateDirectory);
    	Slab->read(filename);
    	STDLOG(1,"Reading SlabSize file from %s\n", filename);
    	TY  = new SlabTaylor(order,cpd);
    	RL = new Redlack(cpd);

    	SlabForceTime = new STimer[cpd];
    	SlabForceLatency = new STimer[cpd];
    	SlabFarForceTime = new STimer[cpd];

    	JJ = new NearFieldDriver();
    	RL->ReadInAuxiallaryVariables(P.ReadStateDirectory);
    } else {
    	TY = NULL;
    	RL = NULL;
    	JJ = NULL;
    }

    prologue.Stop();
    STDLOG(1,"Leaving Prologue()\n");
}

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

    delete MF;
    delete LBW;
    delete PP;
    delete IL;
    delete Slab;


    if(!ic) {
        if(P.ForceOutputDebug){
            #ifdef CUDADIRECT
            STDLOG(1,"Direct Interactions: CPU (%llu) and GPU (%llu)\n",
                        JJ->DirectInteractions_CPU,JJ->DirectInteractions_GPU());
            if(!(JJ->DirectInteractions_CPU == JJ->DirectInteractions_GPU())){
                printf("Error:\n\tDirect Interactions differ between CPU (%llu) and GPU (%llu)\n",
                        JJ->DirectInteractions_CPU,JJ->DirectInteractions_GPU());
                assert(JJ->DirectInteractions_CPU == JJ->DirectInteractions_GPU());
            }
            #endif
        }
    	delete TY;
    	delete RL;
    	delete[] SlabForceLatency;
    	delete[] SlabForceTime;
    	delete[] SlabFarForceTime;
    	delete JJ;
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
