/* timestep.cpp

This is the main routine to evolve the particles.  This defines the
slab-based pipeline that we will use by setting up a set of 
preconditions and actions.  It also handles most of the allocation
and deallocation of arenas.  

We should endeavor to make the basic outline of the pipeline very
clear in this source file.

The bottom of this file contains a redacted version of the pipeline
for the creation of the initial state.  This just loads particles
to the insert list and then calls finish.

We also provide another simplified pipeline to recover multipoles
from position slabs.  This is invoked via the `recover_multipoles`
executable.

*/

int FORCE_RADIUS = -1;
int GROUP_RADIUS = -1;

// I think in most cases we would prefer to read ahead until memory limited
//#define FETCHAHEAD (2*FORCE_RADIUS)
//#define FETCHAHEAD 1000
#define FETCHAHEAD FORCE_RADIUS + 3
#define FETCHPERSTEP 1
// Recall that all of these Dependencies have a built-in STimer
// to measure the amount of time spent on Actions.
Dependency FetchSlabs;
Dependency TransposePos;
Dependency NearForce;
Dependency TaylorForce;
Dependency Kick;
Dependency MakeCellGroups;
Dependency FindCellGroupLinks;
Dependency DoGlobalGroups;
Dependency Output;
Dependency Microstep;
Dependency FinishGroups;
Dependency Drift;
Dependency Finish;

Dependency LPTVelocityReRead;

// The wall-clock time minus all of the above Timers might be a measure
// of the spin-locked time in the timestep() loop.
STimer TimeStepWallClock;

#include "manifest.cpp"

// -----------------------------------------------------------------
/*
 * The precondition for loading new slabs into memory
 * We limit the additional slabs read to FETCHAHEAD
 */
int FetchSlabsPrecondition(int slab) {
    if(slab > Kick.last_slab_executed + FETCHAHEAD)
        // This was +1, but for non-blocking reads 
        // I think we want to work one more ahead
        return 0;
    
    #ifdef PARALLEL
    if (FetchSlabs.raw_number_executed>=total_slabs_on_node) return 0;
    	// This prevents FetchSlabAction from reading beyond the 
	// range of slabs on the node.  In the PARALLEL code, these
	// data will arrive from the Manifest.  We have to implement
	// this on the raw_number because the reported number is adjusted
	// by the Manifest, which leads to a race condition when running
	// the PARALLEL code on a single node test.
	// In the non-PARALLEL code, we have intentionally re-queued some slabs,
	// so don't apply this test.
    #endif
    return 1;
}

/*
 * Loads a set of slabs at a common x slice into memory
 * All "normal" slabtypes should be loaded here. Note that loads may be async.
 */
void FetchSlabsAction(int slab) {
    STDLOG(0,"Fetching slab %d with %d particles\n", slab, SS->size(slab));
    // Load all of the particle files together
    SB->LoadArenaNonBlocking(CellInfoSlab,slab);
    SB->LoadArenaNonBlocking(PosSlab,slab);
    assertf(SS->size(slab)*sizeof(posstruct)<=
        fsize(SB->ReadSlabPath(PosSlab,slab).c_str()),
        "PosSlab size doesn't match prediction\n");

    // Don't bother to load the vel/aux/taylors for slabs that won't be kicked until the wrap
    #ifndef PARALLEL
    if(FetchSlabs.number_of_slabs_executed < FORCE_RADIUS)
        return;
    #endif

    SB->LoadArenaNonBlocking(VelSlab, slab);
    assertf(SS->size(slab)*sizeof(velstruct)<=
        fsize(SB->ReadSlabPath(VelSlab,slab).c_str()),
        "VelSlab size doesn't match prediction\n");
    SB->LoadArenaNonBlocking(AuxSlab, slab);
    assertf(SS->size(slab)*sizeof(auxstruct)<=
        fsize(SB->ReadSlabPath(AuxSlab, slab).c_str()),
        "AuxSlab size doesn't match prediction\n");
    SB->LoadArenaNonBlocking(TaylorSlab,slab);
}

// -----------------------------------------------------------------

int TransposePosPrecondition(int slab){
    if(   !SB->IsIOCompleted(PosSlab,      slab) ||
          !SB->IsIOCompleted(CellInfoSlab, slab)   ) {
        Dependency::NotifySpinning(WAITING_FOR_IO);
        return 0;
    }
        return 1;
}

void TransposePosAction(int slab){
    STDLOG(0,"Transposing position slab %d with %d particles\n", slab, SS->size(slab));
    
    SB->AllocateArena(PosXYZSlab, slab);
    int cpd = P.cpd;
    
    // Could do this over skewers; should make a skewersize(slab, y) function somewhere
    #pragma omp parallel for schedule(static)
    for(int y = 0; y < cpd; y++){
        for(int z = 0; z < cpd; z++){
            posstruct *pos = CP->PosCell(slab, y, z);
            List3<FLOAT> posxyz = CP->PosXYZCell(slab, y, z);
            int count = CP->NumberParticle(slab,y,z);
            
            #pragma ivdep
            for(int i = 0; i < count; i++){
                posxyz.X[i] = pos[i].x;
                posxyz.Y[i] = pos[i].y;
                posxyz.Z[i] = pos[i].z;
            }
        }
    }
    
    #ifndef PARALLEL
    // If this is a "ghost" slab, we only need its transpose
    if(TransposePos.number_of_slabs_executed < FORCE_RADIUS)
        SB->DeAllocate(PosSlab, slab);
    #endif
}


// -----------------------------------------------------------------

int NearForcePrecondition(int slab) {
    for(int i=-FORCE_RADIUS;i<=FORCE_RADIUS;i++){
        // Technically, I think we only need the CellInfo to construct pencils
        // But it's convenient to have pos so the GPU can immediately execute any pencil
        if(TransposePos.notdone(slab+i))
            return 0;
        if( !SB->IsIOCompleted( CellInfoSlab, slab+i ) ){
            Dependency::NotifySpinning(WAITING_FOR_IO);
            return 0;
        }
    }
    
    return 1;
}

void NearForceAction(int slab) {
    // Do some data checks
    assertf(are_cellinfo_legal(slab, SS->size(slab)),
            "Cell info of slab %d contain out of bounds data\n", slab);
    // Could also check that the sum of the cell counts add up to SS->size(slab);

    STDLOG(1,"Computing near-field force for slab %d\n", slab);
    SlabForceTime[slab].Start();
        
    NFD->ExecuteSlab(slab, P.ForceOutputDebug);
    //NFD->ExecuteSlab(slab, 1);  // Use this line instead to force blocking GPU work

    SlabForceLatency[slab].Start();
    if (P.ForceOutputDebug) {
        // We want to output the AccSlab to the NearAcc file.
        // This must be a blocking write.
        NFD->Finalize(slab);

#ifdef DIRECTSINGLESPLINE
        // Single spline requires a prefactor multiplication, which we defer to the kick for efficiency
        // But analysis routines that use ForceOutputDebug, like Ewald, expect this prefactor to already be applied
        // So apply it here, storing the original in a temporary copy
        uint64 npslab = SS->size(slab);
        accstruct *nearacctmp = new accstruct[npslab];
        accstruct *nearacc = (accstruct *) SB->GetSlabPtr(AccSlab, slab);
        memcpy(nearacctmp, nearacc, npslab*sizeof(accstruct));
        FLOAT inv_eps3 = 1./(NFD->SofteningLengthInternal*NFD->SofteningLengthInternal*NFD->SofteningLengthInternal);
        for(uint64 i = 0; i < npslab; i++)
            nearacc[i] *= inv_eps3; 
#endif
        SB->WriteArena(AccSlab, slab, IO_KEEP, IO_BLOCKING,
            SB->WriteSlabPath(NearAccSlab,slab).c_str());

#ifdef DIRECTSINGLESPLINE
        // restore the original
        memcpy(nearacc, nearacctmp, npslab*sizeof(accstruct));
        delete[] nearacctmp;
#endif

    }

    // Busy-wait for all GPU work for this slab to finish
    // while(!NFD->SlabDone(slab)) ;
}

// -----------------------------------------------------------------

int TaylorForcePrecondition(int slab) {
    if( !SB->IsIOCompleted( CellInfoSlab,  slab ) ||
        !SB->IsIOCompleted( PosSlab,       slab ) ||
        !SB->IsIOCompleted( TaylorSlab,    slab ) ){
        Dependency::NotifySpinning(WAITING_FOR_IO);
        return 0;
    }
    return 1;
}

void TaylorForceAction(int slab) {
     // We finished reading this TaylorSlab, so we can delete it to save space
    if (WriteState.OverwriteConvState){
        STDLOG(1, "Deleting TaylorSlab %d since we have finished reading it\n",slab);
        assertf(remove(SB->ReadSlabPath(TaylorSlab,slab).c_str()) == 0, "Could not remove TaylorSlab %d\n",slab);
    }

    STDLOG(1,"Computing far-field force for slab %d\n", slab);
    SlabFarForceTime[slab].Start();
    SB->AllocateArena(FarAccSlab, slab);
    
    TaylorCompute.Start();
    ComputeTaylorForce(slab);
    TaylorCompute.Stop();

    if(P.ForceOutputDebug){
        // We want to output the FarAccSlab to the FarAcc file.
        // This must be a blocking write.
        SB->WriteArena(FarAccSlab, slab, IO_KEEP, IO_BLOCKING);
    }
    SB->DeAllocate(TaylorSlab,slab);
    SlabFarForceTime[slab].Stop();
}


// -----------------------------------------------------------------

int KickPrecondition(int slab) {
    if( !SB->IsIOCompleted( VelSlab, slab ) ) {
        Dependency::NotifySpinning(WAITING_FOR_IO);
        return 0;
    }
    
    // must have far forces
    if(TaylorForce.notdone(slab)){
        return 0;
    }
    
    //must have near forces
    if (NearForce.notdone(slab)) 
        return 0;
    else {
        // If pencil construction (NearForce) is finished, but not NFD, then we're waiting for the GPU result
        if(!NFD->SlabDone(slab)){
#ifdef CUDADIRECT
            Dependency::NotifySpinning(WAITING_FOR_GPU);
#endif
            return 0;
        }
    }

    return 1;
}

void KickAction(int slab) {
    SlabForceTime[slab].Stop();
    SlabForceLatency[slab].Stop();

    // Release the trailing slab if it won't be needed at the wrap
    // Technically we could release it anyway and re-do the transpose from PosSlab,
    // but if we're not doing group finding we may have already written and released PosSlab
    if(Kick.number_of_slabs_executed >= 2*FORCE_RADIUS)
        SB->DeAllocate(PosXYZSlab, slab - FORCE_RADIUS);

    // Special case: if this is the last slab, free all +/- FORCE_RADIUS
    STDLOG(1,"%d slabs have been Kicked so far\n", Kick.number_of_slabs_executed);
    if(Kick.number_of_slabs_executed == CP->cpd-1)
        for(int j = slab - FORCE_RADIUS+1; j <= slab + FORCE_RADIUS; j++)
            SB->DeAllocate(PosXYZSlab, j);

    // Queue up slabs near the wrap to be loaded again later
    // This way, we don't have idle slabs taking up memory while waiting for the pipeline to wrap around
    #ifndef PARALLEL
    if(Kick.number_of_slabs_executed < FORCE_RADIUS){
        STDLOG(2,"Marking slab %d for repeat\n", slab - FORCE_RADIUS);
        TransposePos.mark_to_repeat(slab - FORCE_RADIUS);
	// BUG FIXED: This DeAllocation was missing
        SB->DeAllocate(PosXYZSlab, slab - FORCE_RADIUS);
        // The first two won't need PosSlab until the second time around
        //SB->DeAllocate(PosSlab, slab - FORCE_RADIUS);
        SB->DeAllocate(CellInfoSlab, slab - FORCE_RADIUS);
        FetchSlabs.mark_to_repeat(slab - FORCE_RADIUS);
    }
    #endif

    //If we are doing blocking forces, the finalization happens in NearForceAction
    if(!P.ForceOutputDebug && !P.ForceCPU)
        NFD->Finalize(slab);
    AddAccel.Start();
    RescaleAndCoAddAcceleration(slab);
    SB->DeAllocate(FarAccSlab,slab);
    AddAccel.Stop();
    int step = LPTStepNumber();
    KickCellTimer.Start();
    if (step) {
        // We have LPT IC work to do
        if (step==1) {
            STDLOG(1,"Kicking slab %d as LPT step 1\n", slab);
            KickSlab(slab, 0, 0, KickCell_2LPT_1);
        } else if (step==2) {
            STDLOG(1,"Kicking slab %d as LPT step 2\n", slab);
            KickSlab(slab, 0, 0, KickCell_2LPT_2);
        } else if (step==3) {
            STDLOG(1,"Kicking slab %d as LPT step 3\n", slab);
            KickSlab(slab, 0, 0, KickCell_2LPT_3);
        } else QUIT("LPT Kick %d not implemented\n", step);
    } else {
        // This is just a standard step
        FLOAT kickfactor1 =  ReadState.LastHalfEtaKick;
        FLOAT kickfactor2 =  WriteState.FirstHalfEtaKick;
        STDLOG(1,"Kicking slab %d by %f + %f\n", slab, kickfactor1, kickfactor2);
        KickSlab(slab, kickfactor1, kickfactor2, KickCell);
    }
    KickCellTimer.Stop();
}

// -----------------------------------------------------------------

int MakeCellGroupsPrecondition(int slab) {
    // Only PosXYZ is used for sources, so we're free to rearrange PosSlab
    // in group finding after the transpose
    if( TransposePos.notdone(slab) ) return 0;
    
    // Need the new velocities because we're going to rearrange particles
    if( Kick.notdone(slab) ) return 0;
    
    // Also need the auxs, because we're going to re-order
    if( !SB->IsIOCompleted( AuxSlab, slab ) ) {
        Dependency::NotifySpinning(WAITING_FOR_IO);
        return 0;
    }
    return 1;
}

void MakeCellGroupsAction(int slab) {
	STDLOG(1,"Making Cell Groups in slab %d\n", slab);
	GFC->ConstructCellGroups(slab);
}

// -----------------------------------------------------------------

int FindCellGroupLinksPrecondition(int slab) {
    // We want to find all links between this slab and the one just behind
    for (int j=-1; j<=0; j++) 
        if (MakeCellGroups.notdone(slab+j)) return 0;
    return 1;
}

void FindCellGroupLinksAction(int slab) {
    // Find links between slab and slab-1
	STDLOG(1,"Finding Group Links between slab %d and %d\n", slab, slab-1);
	FindGroupLinks(slab);
}

// -----------------------------------------------------------------

int DoGlobalGroupsPrecondition(int slab) {
    // We're going to close all CellGroups in this slab.
    // GlobalGroups can span 2*GroupRadius+1.
    // But even though we usually encounter a CellGroup in its minimum slab,
    // we could be anywhere in the first instance.  So we have to query a big range.
    // That said, if the nearby slab has already closed global groups, then
    // we can proceed.  This particularly matters in the parallel version, where
    // we may already have closed groups in higher numbered slabs.
    if (Kick.notdone(slab)) return 0;
    // Look behind; can stop as soon as one finds a closed slab
    // The lower bound has a +1 (> not >=) because FindLinks(n) connects n and n-1
    for (int j=0; j>-2*GROUP_RADIUS; j--) {
        if (DoGlobalGroups.done(slab+j-1)) break;
        if (FindCellGroupLinks.notdone(slab+j)) return 0;
    }
    // Look ahead; can stop as soon as one finds a closed slab
    for (int j=1; j<=2*GROUP_RADIUS; j++){
        if (DoGlobalGroups.done(slab+j)) break;
        if (FindCellGroupLinks.notdone(slab+j)) return 0;
    }
    return 1;
}

void DoGlobalGroupsAction(int slab) {
    STDLOG(0,"Finding Global Groups in slab %d\n", slab);
    FindAndProcessGlobalGroups(slab);
}

// -----------------------------------------------------------------
/*
 * Checks if we are ready to do all outputs for this step
 * Anything that modifies the particles at the current time should happen before here
 */
int OutputPrecondition(int slab) {
    if (DoGlobalGroups.notdone(slab)) return 0;  // Must have found groups to be able to output light cones
    // note that group outputs were already done

    // Note the following conditions only have any effect if group finding is turned off
    
    if (Kick.notdone(slab)) return 0;  // Must have accelerations

    // Also obviously need the aux!
    if( !SB->IsIOCompleted( AuxSlab, slab ) ) {
        Dependency::NotifySpinning(WAITING_FOR_IO);
        return 0;
    }
    
    return 1;
}

uint64 n_output = 0;
void OutputAction(int slab) {
    STDLOG(1,"Output slab %d\n", slab);

    int step = WriteState.FullStepNumber;
    if (LPTStepNumber()>0) return;
    // Some output might want to be skipped during an IC step,
    // e.g., no light cones

    OutputTimeSlice.Start();

    if (ReadState.DoTimeSliceOutput) {
        // We've already done a K(1) and thus need a K(-1/2)
        FLOAT unkickfactor = WriteState.FirstHalfEtaKick;
        STDLOG(1,"Outputting slab %d with unkick factor %f\n",slab, unkickfactor);
        n_output += Output_TimeSlice(slab, unkickfactor);
    }
    OutputTimeSlice.Stop();

    OutputLightCone.Start();
    if (ReadState.OutputIsAllowed) {
        // TODO: LightCones may need a half un-kick if GFC == NULL
        // but we can probably handle that in the interpolation
        for(int i = 0; i < P.NLightCones; i++){
            STDLOG(1,"Outputting LightCone %d (origin (%f,%f,%f)) for slab %d\n",i,LCOrigin[i].x,LCOrigin[i].y,LCOrigin[i].z,slab);
            makeLightCone(slab,i);
        }
    }
    OutputLightCone.Stop();

    OutputBin.Start();
    if(ReadState.DoBinning){
        int zstride = CP->cpd /omp_get_max_threads();
        int ystride = CP->cpd /omp_get_max_threads();
        int minstride = 12;
        if (ystride < minstride) ystride = minstride;
        if (zstride < minstride) zstride = minstride;
        int cpd = CP->cpd;
        STDLOG(1,"Binning particles for slab %d\n",slab);
        #pragma omp parallel for schedule(dynamic,ystride)
        for (int y=0;y<cpd;y++) {
            for (int z=0;z<cpd;z++) {
                Cell c = CP->GetCell(slab, y, z);
                tsc(c.pos,CP->CellCenter(slab,y,z),density,c.count(),P.PowerSpectrumN1d,1.0);
            }
        }
    }
    OutputBin.Stop();

}

// -----------------------------------------------------------------

int MicrostepPrecondition(int slab){
    // We are going to second-half kick this slab
    if (Output.notdone(slab))
        return 0;
    return 1;
}

void MicrostepAction(int slab){
    STDLOG(1,"Starting microsteps for slab %d\n", slab);

    // All kicks (and half-unkicks) for output are done; discard accels.
    // We de-allocate in Drift if we aren't doing group finding
    SB->DeAllocate(AccSlab,slab);

    return;
    MicrostepCPU.Start();
    // Do microstepping here
    if(MicrostepEpochs != NULL){
        STDLOG(1,"Beginning microsteps for slab %d\n", slab);
        MicrostepControl *MC = new MicrostepControl;
        MC->setup(GFC->globalslabs[slab], *MicrostepEpochs, P.MicrostepTimeStep, NFD->eps);
        //MC->LaunchGroupsGPU();
        MC->ComputeGroupsCPU();

        GFC->microstepcontrol[slab] = MC;
    }
    MicrostepCPU.Stop();
}

// -----------------------------------------------------------------

int FinishGroupsPrecondition(int slab){
    // Is the asychronous GPU microstepping done?
    //if (!GFC->microstepcontrol[slab]->GPUGroupsDone()) return 0

    // We are going to release these groups.
    if (Microstep.notdone(slab)) return 0;
    
    return 1;
}

void FinishGroupsAction(int slab){
    // Scatter pos,vel updates to slabs, and release GGS
    STDLOG(0, "Finishing groups in slab %d\n", slab);
    delete GFC->microstepcontrol[slab];
    GFC->microstepcontrol[slab] = NULL;
    FinishGlobalGroups(slab);
    GFC->DestroyCellGroups(slab);
}

// -----------------------------------------------------------------
/*
 * Checks if we are ready to load the LPT velocities during an IC step.
 * Should not happen in normal execution
 */
int FetchLPTVelPrecondition(int slab){
    // Don't read too far ahead
    if(LPTVelocityReRead.number_of_slabs_executed > 
            Drift.number_of_slabs_executed + 2*FINISH_WAIT_RADIUS + 1) {
        return 0;
    }

    return 1;
}

void FetchLPTVelAction(int slab){
    // This is blocking because it uses the LoadIC module, not SB
    load_ic_vel_slab(slab);
}

// -----------------------------------------------------------------
/* Checks if we are ready to apply a full-timestep drift to the particles
 * Any operation that relies on the positions and velocities being synced
 * should be checked for completion this year.
 */
int DriftPrecondition(int slab) {
    // We must have finished scattering into this slab
    if (FinishGroups.notdone(slab)) return 0;
    
    // We will move the particles, so we must have done outputs
    if (Output.notdone(slab)) return 0;
    
    // We can't move particles until they've been used as gravity sources
    // However, we only use PosXYZSlab as sources, so we're free to move PosSlab
    
    // We also must have the 2LPT velocities
    // The finish radius is a good guess of how ordered the ICs are
    for(int i=-FINISH_WAIT_RADIUS;i<=FINISH_WAIT_RADIUS;i++) 
        if (LPTVelocityReRead.notdone(slab+i)) return 0;
        
    return 1;
}

void DriftAction(int slab) {
    int step = LPTStepNumber();
    if (step) {
        // We have LPT IC work to do
        if (step==1) {
            STDLOG(1,"Drifting slab %d as LPT step 1\n", slab);
            DriftAndCopy2InsertList(slab, 0, DriftCell_2LPT_1);
        } else if (step==2) {
            STDLOG(1,"Drifting slab %d as LPT step 2\n", slab);
            DriftAndCopy2InsertList(slab, 0, DriftCell_2LPT_2);
        } else if (step==3) {
            STDLOG(1,"Drifting slab %d as LPT step 3\n", slab);
            DriftAndCopy2InsertList(slab, 0, DriftCell_2LPT_3);
        } else QUIT("LPT Drift %d not implemented\n", step);
    } else {
        //         This is just a normal drift
        FLOAT driftfactor = WriteState.DeltaEtaDrift;
        // WriteState.etaD-ReadState.etaD; 
        STDLOG(1,"Drifting slab %d by %f\n", slab, driftfactor);
        //DriftAndCopy2InsertList(slab, driftfactor, DriftCell);
        DriftSlabAndCopy2InsertList(slab, driftfactor, DriftSlab);
    }
    
    // We freed AccSlab in Microstep to save space
    if (GFC == NULL){
	    // We kept the accelerations until here because of third-order LPT
	    if (P.StoreForces && !P.ForceOutputDebug) {
	        STDLOG(1,"Storing Forces in slab %d\n", slab);
	        SB->StoreArenaBlocking(AccSlab,slab);
	    }
	    else{
	        SB->DeAllocate(AccSlab,slab);
	    }
	}
}

// -----------------------------------------------------------------

int FinishPrecondition(int slab) {
    for(int j=-FINISH_WAIT_RADIUS;j<=FINISH_WAIT_RADIUS;j++) {
        if( Drift.notdone(slab+j) ) return 0;
    }
    
    return 1;
}

uint64 merged_particles = 0;
void FinishAction(int slab) {
    STDLOG(1,"Finishing slab %d\n", slab);
    
    if (WriteState.Do2LPTVelocityRereading)
        SB->DeAllocate(VelLPTSlab, slab);
    
    // Gather particles from the insert list and make the merge slabs
    uint64 n_merge = FillMergeSlab(slab);
    merged_particles += n_merge;
    
    // This may be the last time be need any of the CellInfo slabs that we just used
    // We can't immediately free CellInfo before NearForce might need it until we're FORCE_RADIUS away
    // An alternative to this would be to just wait for FORCE_RADIUS before finishing
    for(int j = -2*FORCE_RADIUS, consec = 0; j <= 2*FORCE_RADIUS; j++){
        if(Finish.done(slab + j) || j == 0)
            consec++;
        else
            consec = 0;
        if (consec >= 2*FORCE_RADIUS + 1)
            SB->DeAllocate(CellInfoSlab, slab + j - FORCE_RADIUS);
    }
    
    // Now delete the original particles
    SB->DeAllocate(PosSlab,slab);
    SB->DeAllocate(VelSlab,slab);
    SB->DeAllocate(AuxSlab,slab);
    
    // Make the multipoles
    SB->AllocateArena(MultipoleSlab,slab);
    ComputeMultipoleSlab(slab);
    
    // Write out the particles and multipoles and delete
    WriteMergeSlab.Start();
    SB->StoreArenaNonBlocking(MergePosSlab,slab);
    SB->StoreArenaNonBlocking(MergeVelSlab,slab);
    SB->StoreArenaNonBlocking(MergeAuxSlab,slab);
    SB->StoreArenaNonBlocking(MergeCellInfoSlab,slab);
    WriteMergeSlab.Stop();

    WriteMultipoleSlab.Start();
    SB->StoreArenaNonBlocking(MultipoleSlab,slab);
    WriteMultipoleSlab.Stop();

    int pwidth = FetchSlabs.number_of_slabs_executed - Finish.number_of_slabs_executed;
    STDLOG(1, "Current pipeline width (N_fetch - N_finish) is %d\n", pwidth);

    #ifdef PARALLEL
    if (Finish.number_of_slabs_executed==0) SendManifest.QueueToSend(slab);
    #endif
    // TODO: is there a different place in the code where we would rather report this?
    ReportMemoryAllocatorStats();
}

// -----------------------------------------------------------------
// A no-op precondition that always passes
int NoopPrecondition(int slab){
    return 1;
}

// A no-op action that does nothing
void NoopAction(int slab){
    return;
}

// ===================================================================

/*
 * Registers all of the dependencies and their associated actions.
 * The Dependency module is responsible for running the registered steps.
 */
void timestep(void) {
    TimeStepWallClock.Clear();
    TimeStepWallClock.Start();
    STDLOG(1,"Initiating timestep()\n");

    FORCE_RADIUS = P.NearFieldRadius;
    GROUP_RADIUS = GFC != NULL ? P.GroupRadius : 0;
    // The 2LPT pipeline is short (no group finding). We can afford to wait an extra slab to allow for large IC displacements
    FINISH_WAIT_RADIUS = LPTStepNumber() > 0 ? 2 : 1;
    assertf(FORCE_RADIUS >= 0, "Illegal FORCE_RADIUS: %d\n", FORCE_RADIUS);
    assertf(GROUP_RADIUS >= 0, "Illegal GROUP_RADIUS: %d\n", GROUP_RADIUS); 
    STDLOG(0,"Adopting FORCE_RADIUS = %d\n", FORCE_RADIUS);
    STDLOG(0,"Adopting GROUP_RADIUS = %d\n", GROUP_RADIUS);
    STDLOG(0,"Adopting FINISH_WAIT_RADIUS = %d\n", FINISH_WAIT_RADIUS);

    int nslabs = P.cpd;
    int first = first_slab_on_node;  // First slab to load
    STDLOG(1,"First slab to load will be %d\n", first);

        FetchSlabs.instantiate(nslabs, first, &FetchSlabsPrecondition,          &FetchSlabsAction         );
      TransposePos.instantiate(nslabs, first, &TransposePosPrecondition,       &TransposePosAction      );
         NearForce.instantiate(nslabs, first + FORCE_RADIUS, &NearForcePrecondition,          &NearForceAction         );
       TaylorForce.instantiate(nslabs, first + FORCE_RADIUS, &TaylorForcePrecondition,        &TaylorForceAction       );
              Kick.instantiate(nslabs, first + FORCE_RADIUS, &KickPrecondition,               &KickAction              );
            Output.instantiate(nslabs, first + FORCE_RADIUS + 2*GROUP_RADIUS, &OutputPrecondition,             &OutputAction            );
             Drift.instantiate(nslabs, first + FORCE_RADIUS + 2*GROUP_RADIUS, &DriftPrecondition,              &DriftAction             );
            Finish.instantiate(nslabs, first + FORCE_RADIUS + 2*GROUP_RADIUS + FINISH_WAIT_RADIUS, &FinishPrecondition,             &FinishAction            );
            
    // If group finding is disabled, we can make the dependencies no-ops so they don't hold up the pipeline
    if(GFC != NULL){
        MakeCellGroups.instantiate(nslabs, first + FORCE_RADIUS, &MakeCellGroupsPrecondition,     &MakeCellGroupsAction    );
    FindCellGroupLinks.instantiate(nslabs, first + FORCE_RADIUS + 1, &FindCellGroupLinksPrecondition, &FindCellGroupLinksAction);
        DoGlobalGroups.instantiate(nslabs, first + FORCE_RADIUS + 2*GROUP_RADIUS, &DoGlobalGroupsPrecondition,     &DoGlobalGroupsAction    );
             Microstep.instantiate(nslabs, first + FORCE_RADIUS + 2*GROUP_RADIUS, &MicrostepPrecondition,          &MicrostepAction         );
          FinishGroups.instantiate(nslabs, first + FORCE_RADIUS + 2*GROUP_RADIUS, &FinishGroupsPrecondition,       &FinishGroupsAction      );
    } else {
        MakeCellGroups.instantiate(nslabs, first, &NoopPrecondition, &NoopAction );
    FindCellGroupLinks.instantiate(nslabs, first, &NoopPrecondition, &NoopAction );
        DoGlobalGroups.instantiate(nslabs, first, &NoopPrecondition, &NoopAction );
             Microstep.instantiate(nslabs, first, &NoopPrecondition, &NoopAction );
          FinishGroups.instantiate(nslabs, first, &NoopPrecondition, &NoopAction );
    }
           
    if(WriteState.Do2LPTVelocityRereading)
        LPTVelocityReRead.instantiate(nslabs, first + FORCE_RADIUS + 2*GROUP_RADIUS - FINISH_WAIT_RADIUS,
                                          &FetchLPTVelPrecondition,   &FetchLPTVelAction   );
    else
        LPTVelocityReRead.instantiate(nslabs, first, &NoopPrecondition, &NoopAction );

    while( !Finish.alldone(total_slabs_on_node) ) {
           for(int i =0; i < FETCHPERSTEP; i++) FetchSlabs.Attempt();
         TransposePos.Attempt();
            NearForce.Attempt();
          TaylorForce.Attempt();
                 Kick.Attempt();
       MakeCellGroups.Attempt();
   FindCellGroupLinks.Attempt();
       DoGlobalGroups.Attempt();
               Output.Attempt();
            Microstep.Attempt();
         FinishGroups.Attempt();
    LPTVelocityReRead.Attempt();
                Drift.Attempt();
               Finish.Attempt();
	    // TODO: The following line will be omitted once the MPI monitoring thread is in place.
           SendManifest.FreeAfterSend();
	    ReceiveManifest.Check();   // This checks if Send is ready; no-op in non-blocking mode
	
	    // If the manifest has been received, install it.
	    if (ReceiveManifest.is_ready()) ReceiveManifest.ImportData();
    }

    if(IL->length!=0)
        IL->DumpParticles();
    
    assertf(IL->length==0, 
        "Insert List not empty (%d) at the end of timestep().  Time step too big?\n", IL->length);
    
    #ifdef PARALLEL
        usleep(1e6);
        STDLOG(1,"Finished timestep loop!!\n");
        unsigned int tmp = merged_particles;
        MPI_Allreduce(MPI_IN_PLACE, &tmp, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
        merged_particles = tmp;
    #endif 
    assertf(merged_particles == P.np, "Merged slabs contain %d particles instead of %d!\n", merged_particles, P.np);

    if (GFC != NULL) assertf(GFC->GLL->length==0,
	"GroupLinkList not empty (%d) at the end of timestep.  Global group finding didn't run properly.\n", GFC->GLL->length);

    uint64 total_n_output = n_output;
    if(GFC != NULL)
        total_n_output += GFC->n_L0_output;
    
    if(ReadState.DoTimeSliceOutput)
        assertf(total_n_output == P.np, "TimeSlice output contains %d particles instead of %d!\n", total_n_output, P.np);

    STDLOG(1,"Completing timestep()\n");
    TimeStepWallClock.Stop();
}


// ===================================================================
// Other timesteps that re-use dependencies above

#include "timestep_ic.cpp"

#include "timestep_recover_multipoles.cpp"

#include "timestep_benchmark_io.cpp"

#include "timestep_standalone_fof.cpp"
