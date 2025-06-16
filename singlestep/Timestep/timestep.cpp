/* timestep.cpp

This is the main routine to evolve the particles.  This defines the
slab-based pipeline that we will use by setting up a set of
preconditions and actions.  It also handles most of the allocation
and deallocation of arenas.

We should endeavor to make the basic outline of the pipeline very
clear in this source file.

`timestep_ic.cpp` contains a similar pipeline, suitable for creating
the initial state.

*/

int FORCE_RADIUS = -1;
int GROUP_RADIUS = -1;

#define FETCHAHEAD (std::max(2*GROUP_RADIUS + FORCE_RADIUS, FINISH_WAIT_RADIUS) + 2)

// Recall that all of these Dependencies have a built-in STimer
// to measure the amount of time spent on Actions.
SlabDependency *FetchSlabs;
SlabDependency *TransposePos;
SlabDependency *NearForce;
SlabDependency *TaylorTranspose;  // 2D
SlabDependency *TaylorForce;
SlabDependency *Kick;
SlabDependency *MakeCellGroups;
SlabDependency *FindCellGroupLinks;
SlabDependency *DoGlobalGroups;
SlabDependency *Output;
SlabDependency *FinishGroups;
SlabDependency *Drift;
SlabDependency *NeighborSend;  // 2D
SlabDependency *FinishParticles;
SlabDependency *FinishMultipoles;
SlabDependency *CheckForMultipoles; // 1D

EventDependency *DoReceiveManifest;
EventDependency *DoSendManifest;
EventDependency *DoNeighborRecv;
EventDependency *Check2DMultipoleMPI;
EventDependency *Check2DTaylorMPI;

// TODO: should we consider de-coupling PARALLEL from the concept of a merged convolve/singlestep?
#ifdef PARALLEL
#include "ConvolutionParametersStatistics.cpp"
#include "InCoreConvolution.cpp"
#include "ParallelConvolution.cpp"
STimer ConvolutionWallClock;
STimer BarrierWallClock;

#include "manifest.cpp"
#endif

// The wall-clock time minus all of the above Timers might be a measure
// of the spin-locked time in the timestep() loop.
STimer TimeStepWallClock;

#include "pipeline.cpp"

void InitializeForceRadius(int NoForces){
    FORCE_RADIUS = NoForces ? 0 : P.NearFieldRadius;
    assertf(FORCE_RADIUS >= 0, "Illegal FORCE_RADIUS: {:d}\n", FORCE_RADIUS);
    STDLOG(0,"Adopting FORCE_RADIUS = {:d}\n", FORCE_RADIUS);

#ifdef PARALLEL
    // The IL will use this to size itself
    if(MPI_size_z > 1) NeighborRecvEvent::receive_ahead = 3;
    else NeighborRecvEvent::receive_ahead = 0;
#endif
}


// This happens much later, after outputs and group finding are planned
void InitializePipelineWidths(int MakeIC){
    GROUP_RADIUS = GFC != NULL ? P.GroupRadius : 0;
    assertf(GROUP_RADIUS >= 0, "Illegal GROUP_RADIUS: {:d}\n", GROUP_RADIUS);
    STDLOG(0,"Adopting GROUP_RADIUS = {:d}\n", GROUP_RADIUS);

    // The 2LPT pipeline is short (no group finding). We can afford to wait an extra slab to allow for large IC displacements
    FINISH_WAIT_RADIUS = (MakeIC || LPTStepNumber()) > 0 ? 2 : 1;

#ifdef PARALLEL
    /* In the parallel code, we're about to send all of the info up to
    slab-1 to the neighbor.  This can cause a problem if the pipeline
    is thin (e.g., no group finding), because the PosXYZSlabs are needed
    over a domain of +-FORCE_RADIUS.

    For the first slab to finish, we have to assure that PosXYZSlab[slab]
    is not needed to Kick any slabs on the neighbor.  That means we must
    have done Kick[slab-FORCE_RADIUS] on this node.

    Further, we have to assure that PosXYZSlab[slab-1] is not still needed
    as a source to any slabs on this node.  Need Kick[slab-1+FORCE_RADIUS]
    to be done to avoid this.

    We fix this by forcing FINISH_WAIT_RADIUS to be big enough.  */
    if (FINISH_WAIT_RADIUS+2*GROUP_RADIUS<FORCE_RADIUS)
        FINISH_WAIT_RADIUS = FORCE_RADIUS-2*GROUP_RADIUS;

    // TODO: Not sure inflating FINISH_WAIT_RADIUS is the best way to deal with this
    int PAD = 0;
    assertf(total_slabs_on_node >= (2*GROUP_RADIUS + 1) + 2*FORCE_RADIUS + 1 + PAD,
        "Not enough slabs on node to close first group!\n");
    assertf(total_slabs_on_node >= 2*GROUP_RADIUS + 2*FORCE_RADIUS + 1 + 2 * FINISH_WAIT_RADIUS + PAD,
        "Not enough slabs on node to finish any slabs!\n");
#endif

    STDLOG(0,"Adopting FINISH_WAIT_RADIUS = {:d}\n", FINISH_WAIT_RADIUS);
}

void timestep(int NoForces) {

#ifdef PARALLEL
    ParallelConvolveDriver = new ParallelConvolution(P.cpd, P.order, P.MultipoleDirectory);
#endif

    TimeStepWallClock.Start();
    STDLOG(1,"Initiating timestep()\n");

    int nslabs = P.cpd;
    int first = first_slab_on_node;  // First slab to load
    STDLOG(1,"First slab to load will be {:d}\n", first);

    #ifdef ONE_SIDED_GROUP_FINDING
        int first_outputslab = FORCE_RADIUS + 2*GROUP_RADIUS + (int)(GROUP_RADIUS > 0);
    #else
        int first_outputslab = FORCE_RADIUS + 2*GROUP_RADIUS;
    #endif
    if(LPTStepNumber() == 2){
        first_outputslab = std::max(FINISH_WAIT_RADIUS,first_outputslab);
    }

    STDLOG(1, "first_outputslab = {:d}, first_outputslab + FINISH_WAIT_RADIUS = {:d}\n", first_outputslab, first_outputslab + FINISH_WAIT_RADIUS);


    FetchSlabs          = new         FetchSlabsDep(nslabs, first);
    Output              = new             OutputDep(nslabs, first + first_outputslab);
    Drift               = new              DriftDep(nslabs, first + first_outputslab);
    FinishParticles     = new    FinishParticlesDep(nslabs, first + first_outputslab + FINISH_WAIT_RADIUS);
    FinishMultipoles    = new   FinishMultipolesDep(nslabs, first + first_outputslab + FINISH_WAIT_RADIUS);

    if(!NoForces){
        TransposePos        = new       TransposePosDep(nslabs, first + 0);
        NearForce           = new          NearForceDep(nslabs, first + FORCE_RADIUS);
        TaylorForce         = new        TaylorForceDep(nslabs, first + FORCE_RADIUS);
        Kick                = new               KickDep(nslabs, first + FORCE_RADIUS);
    } else {
        TransposePos = new NoopDep(nslabs);
        NearForce    = new NoopDep(nslabs);
        TaylorForce  = new NoopDep(nslabs);
        Kick         = new NoopDep(nslabs);
    }

#ifdef PARALLEL
    if(!NoForces){
        TaylorTranspose     = new    TaylorTransposeDep(nslabs, first + FORCE_RADIUS);  // 2D
        Check2DTaylorMPI    = new Taylor2DMPIEvent();
    } else {
        TaylorTranspose     = new NoopDep(nslabs);
        Check2DTaylorMPI    = new NoopEvent();
    }
    
    CheckForMultipoles  = new CheckForMultipolesDep(nslabs, first + first_outputslab + FINISH_WAIT_RADIUS);
    
    DoReceiveManifest   = new ReceiveManifestEvent();
    DoSendManifest      = new SendManifestEvent();
    NeighborSend        = new NeighborSendDep(nslabs, first + first_outputslab + FINISH_WAIT_RADIUS);  // 2D
    DoNeighborRecv      = new NeighborRecvEvent();
    Check2DMultipoleMPI = new Multipole2DMPIEvent();
#else
    TaylorTranspose    = new NoopDep(nslabs);
    NeighborSend       = new NoopDep(nslabs);
    CheckForMultipoles = new NoopDep(nslabs);

    DoReceiveManifest   = new NoopEvent();
    DoSendManifest      = new NoopEvent();
    DoNeighborRecv      = new NoopEvent();
    Check2DMultipoleMPI = new NoopEvent();
    Check2DTaylorMPI    = new NoopEvent();
#endif

    // If group finding is disabled, we can make the dependencies no-ops so they don't hold up the pipeline
    #ifdef ONE_SIDED_GROUP_FINDING
        int first_groupslab = FORCE_RADIUS+1;
    #else
        int first_groupslab = FORCE_RADIUS+2*GROUP_RADIUS;
    #endif
    if(GFC != NULL){
        MakeCellGroups      = new     MakeCellGroupsDep(nslabs, first + FORCE_RADIUS);
        FindCellGroupLinks  = new FindCellGroupLinksDep(nslabs, first + FORCE_RADIUS + 1);
        DoGlobalGroups      = new     DoGlobalGroupsDep(nslabs, first + first_groupslab);
        FinishGroups        = new       FinishGroupsDep(nslabs, first + first_groupslab);
    } else {
        MakeCellGroups      = new NoopDep(nslabs);
        FindCellGroupLinks  = new NoopDep(nslabs);
        DoGlobalGroups      = new NoopDep(nslabs);
        FinishGroups        = new NoopDep(nslabs);
    }

    // Let FetchSlabs start early, in case we want to overlap convolve and IO
    //while(FetchSlabs->Attempt()){}

    #ifdef PARALLEL
    TimeStepWallClock.Stop();

    if(!NoForces){
        ConvolutionWallClock.Start();

        ParallelConvolveDriver->Convolve();

        for (int slab = first + FORCE_RADIUS; slab < first + FORCE_RADIUS + total_slabs_on_node; slab++){
            SB->AllocateArena(TaylorSlab, slab, RAMDISK_NO);
            ParallelConvolveDriver->RecvTaylorSlab(slab);
        }

        ParallelConvolveDriver->SendTaylors(FORCE_RADIUS);

        ConvolutionWallClock.Stop();
        ParallelConvolveDriver->CS.ConvolveWallClock = ConvolutionWallClock.Elapsed();
    }

    TimeStepWallClock.Start();

    // Lightweight setup of z-dimension exchanges
    SetupNeighborExchange(first + first_outputslab + FINISH_WAIT_RADIUS, total_slabs_on_node);
    #endif

	Dependency::SpinTimer.Start();
	int timestep_loop_complete = 0; 
	while (!timestep_loop_complete){

               FetchSlabs->Attempt();
             TransposePos->Attempt();
                NearForce->Attempt();
    
          TaylorTranspose->Attempt();  // 2D
              TaylorForce->Attempt();
                     Kick->Attempt();
    
           MakeCellGroups->Attempt();
       FindCellGroupLinks->Attempt();
           DoGlobalGroups->Attempt();
    
                   Output->Attempt();
             FinishGroups->Attempt();
    
                    Drift->Attempt();
             NeighborSend->Attempt();  // 2D
    
          FinishParticles->Attempt();
         FinishMultipoles->Attempt();
       CheckForMultipoles->Attempt();  // 1D
           
        DoReceiveManifest->Attempt();  // 1D
           DoSendManifest->Attempt();  // 1D
    
           DoNeighborRecv->Attempt();  // 2D
      Check2DMultipoleMPI->Attempt();  // 2D
         Check2DTaylorMPI->Attempt();  // 2D

#ifdef PARALLEL
		timestep_loop_complete = CheckForMultipoles->alldone(total_slabs_on_node);
#else
		timestep_loop_complete = FinishMultipoles->alldone(total_slabs_on_node);
#endif
    }

    Dependency::SpinTimer.Stop();

    if(IL->length!=0)
        IL->DumpParticles();

    assertf(IL->length==0,
        "Insert List not empty ({:d}) at the end of timestep().  Time step too big?\n", IL->length);

    STDLOG(1,"Finished timestep SlabDependency loop!\n");

    if (GFC != NULL) assertf(GFC->GLL->length==0,
	"GroupLinkList not empty ({:d}) at the end of timestep.  Global group finding didn't run properly.\n", GFC->GLL->length);

    uint64 total_n_output = ((OutputDep *) Output)->n_output;  // TODO: how should we do field access while preserving the ability to NoOp?
    if(GFC != NULL)
        total_n_output += GFC->n_L0_output;
	
	uint64 merged_primaries = ((FinishParticlesDep *)FinishParticles)->merged_primaries;
    #ifdef PARALLEL	
        BarrierWallClock.Start();
        // These reductions force some synchronization, at least!
    	MPI_REDUCE_TO_ZERO(&total_n_output,   1, MPI_UINT64_T, MPI_SUM);
        MPI_REDUCE_TO_ZERO(&merged_primaries, 1, MPI_UINT64_T, MPI_SUM);

        MPI_Barrier(comm_global);
		BarrierWallClock.Stop();

        STDLOG(1,"Ready to proceed to the remaining work\n");

        // This MPI call also forces a synchronization over the MPI processes,
        // so things like Reseting GPUs could fire multiple times on one node.
        SendManifest->FreeAfterSend();
        // Run this again, just in case the SlabDependency loop on this node finished
        // before the neighbor received the non-blocking MPI transfer.

	    TimeStepWallClock.Stop(); ConvolutionWallClock.Start(); 
        convtimebuffer = (char*) malloc(CONVTIMEBUFSIZE);   // Need to allocate space for the timings
    	delete ParallelConvolveDriver;
	    ConvolutionWallClock.Stop(); TimeStepWallClock.Start();

    #endif
    if (MPI_rank==0)
        assertf( merged_primaries == P.np, "Merged slabs contain {:d} particles instead of {:d}!\n", merged_primaries, P.np);


    if(ReadState.DoTimeSliceOutput && MPI_rank==0){
        assertf(total_n_output == P.np, "TimeSlice output contains {:d} particles instead of {:d}!\n", total_n_output, P.np);
	}

    STDLOG(1,"Completing timestep()\n");
    TimeStepWallClock.Stop();
}


void free_dependencies(){
    delete FetchSlabs;
    delete TransposePos;
    delete NearForce;
    delete TaylorTranspose;
    delete TaylorForce;
    delete Kick;
    delete MakeCellGroups;
    delete FindCellGroupLinks;
    delete DoGlobalGroups;
    delete Output;
    delete FinishGroups;
    delete Drift;
    delete NeighborSend;
    delete FinishParticles;
    delete FinishMultipoles;
    delete CheckForMultipoles;

    delete DoReceiveManifest;
    delete DoSendManifest;
    delete DoNeighborRecv;
    delete Check2DMultipoleMPI;
    delete Check2DTaylorMPI;
}
