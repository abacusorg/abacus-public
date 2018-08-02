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
from position slabs.  This is invoked via the `make_multipoles`
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

// -----------------------------------------------------------------
/*
 * The precondition for loading new slabs into memory
 * We limit the additional slabs read to FETCHAHEAD
 */
int FetchSlabPrecondition(int slab) {
    if(slab > Kick.last_slab_executed + FETCHAHEAD)
        // This was +1, but for non-blocking reads 
        // I think we want to work one more ahead
        return 0;
    
    return 1;
}

/*
 * Loads a set of slabs at a common x slice into memory
 * All "normal" slabtypes should be loaded here. Note that loads may be async.
 */
void FetchSlabAction(int slab) {
    STDLOG(0,"Fetching slab %d with %d particles\n", slab, Slab->size(slab));
    // Load all of the particle files together
    LBW->LoadArenaNonBlocking(CellInfoSlab,slab);
    LBW->LoadArenaNonBlocking(PosSlab,slab);
    assertf(Slab->size(slab)*sizeof(posstruct)<=
        fsize(LBW->ReadSlabDescriptorName(PosSlab,slab).c_str()),
        "PosSlab size doesn't match prediction\n");

    // Don't bother to load the vel/aux/taylors for slabs that won't be kicked until the wrap
    if(FetchSlabs.number_of_slabs_executed < FORCE_RADIUS)
        return;

    LBW->LoadArenaNonBlocking(VelSlab, slab);
    assertf(Slab->size(slab)*sizeof(velstruct)<=
        fsize(LBW->ReadSlabDescriptorName(VelSlab,slab).c_str()),
        "VelSlab size doesn't match prediction\n");
    LBW->LoadArenaNonBlocking(AuxSlab, slab);
    assertf(Slab->size(slab)*sizeof(auxstruct)<=
        fsize(LBW->ReadSlabDescriptorName(AuxSlab, slab).c_str()),
        "AuxSlab size doesn't match prediction\n");
    LBW->LoadArenaNonBlocking(TaylorSlab,slab);
}

// -----------------------------------------------------------------

int TransposePosPrecondition(int slab){
    if(   !LBW->IOCompleted(PosSlab,      slab) ||
          !LBW->IOCompleted(CellInfoSlab, slab)   ) {
        Dependency::NotifySpinning(WAITING_FOR_IO);
        return 0;
    }
        return 1;
}

void TransposePosAction(int slab){
    STDLOG(0,"Transposing position slab %d with %d particles\n", slab, Slab->size(slab));
    
    LBW->AllocateArena(PosXYZSlab, slab);
    int cpd = P.cpd;
    
    // Could do this over skewers; should make a skewersize(slab, y) function somewhere
    #pragma omp parallel for schedule(static)
    for(int y = 0; y < cpd; y++){
        for(int z = 0; z < cpd; z++){
            posstruct *pos = PP->PosCell(slab, y, z);
            List3<FLOAT> posxyz = PP->PosXYZCell(slab, y, z);
            int count = PP->NumberParticle(slab,y,z);
            
            #pragma ivdep
            for(int i = 0; i < count; i++){
                posxyz.X[i] = pos[i].x;
                posxyz.Y[i] = pos[i].y;
                posxyz.Z[i] = pos[i].z;
            }
        }
    }
    
    // If this is a "ghost" slab, we only need its transpose
    if(TransposePos.number_of_slabs_executed < FORCE_RADIUS)
        LBW->DeAllocate(PosSlab, slab);
}


// -----------------------------------------------------------------

int NearForcePrecondition(int slab) {
    for(int i=-FORCE_RADIUS;i<=FORCE_RADIUS;i++){
        if(TransposePos.notdone(slab+i))
            return 0;
        if( !LBW->IOCompleted( CellInfoSlab, slab+i ) ){
            Dependency::NotifySpinning(WAITING_FOR_IO);
            return 0;
        }
    }
    
    return 1;
}

void NearForceAction(int slab) {
    // Do some data checks
    assertf(are_cellinfo_legal(slab, Slab->size(slab)),
            "Cell info of slab %d contain out of bounds data\n", slab);
    // Could also check that the sum of the cell counts add up to Slab.size

    STDLOG(1,"Computing near-field force for slab %d\n", slab);
    SlabForceTime[slab].Start();
        
    JJ->ExecuteSlab(slab, P.ForceOutputDebug);
    //JJ->ExecuteSlab(slab, 1);  // Uncomment to force blocking

    SlabForceLatency[slab].Start();
    if (P.ForceOutputDebug) {
        // We want to output the AccSlab to the NearAcc file.
        // This must be a blocking write.
        JJ->Finalize(slab);

#ifdef DIRECTSINGLESPLINE
        // Single spline requires a prefactor multiplication, which we defer to the kick for efficiency
        // But analysis routines that use ForceOutputDebug, like Ewald, expect this prefactor to already be applied
        // So apply it here, storing the original in a temporary copy
        uint64 npslab = Slab->size(slab);
        accstruct *nearacctmp = new accstruct[npslab];
        accstruct *nearacc = (accstruct *) LBW->ReturnIDPtr(AccSlab, slab);
        memcpy(nearacctmp, nearacc, npslab*sizeof(accstruct));
        FLOAT inv_eps3 = 1./(JJ->SofteningLengthInternal*JJ->SofteningLengthInternal*JJ->SofteningLengthInternal);
        for(int i = 0; i < npslab; i++)
            nearacc[i] *= inv_eps3;
#endif
        LBW->WriteArena(AccSlab, slab, IO_KEEP, IO_BLOCKING,
        LBW->WriteSlabDescriptorName(NearAccSlab,slab).c_str());

#ifdef DIRECTSINGLESPLINE
        // restore the original
        memcpy(nearacc, nearacctmp, npslab*sizeof(accstruct));
        delete[] nearacctmp;
#endif

    }
}

// -----------------------------------------------------------------

int TaylorForcePrecondition(int slab) {
    if( !LBW->IOCompleted( CellInfoSlab,  slab ) ||
        !LBW->IOCompleted( PosSlab,       slab ) ||
        !LBW->IOCompleted( TaylorSlab,    slab ) ){
        Dependency::NotifySpinning(WAITING_FOR_IO);
        return 0;
    }
    return 1;
}

void TaylorForceAction(int slab) {
    // We finished reading this TaylorSlab, so we can delete it to save space
    if (P.OverwriteState){
        STDLOG(1, "Deleting TaylorSlab %d since we have finished reading it\n",slab);
        assertf(remove(LBW->ReadSlabDescriptorName(TaylorSlab,slab).c_str()) == 0, "Could not remove TaylorSlab %d\n",slab);
    }
    
    STDLOG(1,"Computing far-field force for slab %d\n", slab);
    SlabFarForceTime[slab].Start();
    LBW->AllocateArena(FarAccSlab,slab);
    
    TaylorCompute.Start();
    ComputeTaylorForce(slab);
    TaylorCompute.Stop();

    if(P.ForceOutputDebug){
        // We want to output the FarAccSlab to the FarAcc file.
        // This must be a blocking write.
        LBW->WriteArena(FarAccSlab, slab, IO_KEEP, IO_BLOCKING,
                LBW->WriteSlabDescriptorName(FarAccSlab,slab).c_str());
    }
    LBW->DeAllocate(TaylorSlab,slab);
    SlabFarForceTime[slab].Stop();
}


// -----------------------------------------------------------------

int KickPrecondition(int slab) {
    if( !LBW->IOCompleted( VelSlab, slab ) ) {
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
        // If pencil construction (NearForce) is finished, but not JJ, then we're waiting for the GPU result
        if(!JJ->SlabDone(slab)){
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
        LBW->DeAllocate(PosXYZSlab, slab - FORCE_RADIUS);

    // Special case: if this is the last slab, free all +/- FORCE_RADIUS
    if(Kick.number_of_slabs_executed == PP->cpd-1)
        for(int j = slab - FORCE_RADIUS+1; j <= slab + FORCE_RADIUS; j++)
            LBW->DeAllocate(PosXYZSlab, j);

    // Queue up slabs near the wrap to be loaded again later
    // This way, we don't have idle slabs taking up memory while waiting for the pipeline to wrap around
    // TODO: Could this be 2*FORCE_RADIUS for PosXYZ?
    if(Kick.number_of_slabs_executed < FORCE_RADIUS){
        STDLOG(2,"Marking slab %d for repeat\n", slab - FORCE_RADIUS);
        TransposePos.mark_to_repeat(slab - FORCE_RADIUS);
	// BUG FIXED: This DeAllocation was missing
        LBW->DeAllocate(PosXYZSlab, slab - FORCE_RADIUS);
        // The first two won't need PosSlab until the second time around
        //LBW->DeAllocate(PosSlab, slab - FORCE_RADIUS);
        LBW->DeAllocate(CellInfoSlab, slab - FORCE_RADIUS);
        FetchSlabs.mark_to_repeat(slab - FORCE_RADIUS);
    }

    //If we are doing blocking forces, the finalization happens in NearForceAction
    if(!P.ForceOutputDebug && !P.ForceCPU)
        JJ->Finalize(slab);
    AddAccel.Start();
    RescaleAndCoAddAcceleration(slab);
    LBW->DeAllocate(FarAccSlab,slab);
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
    if( !LBW->IOCompleted( AuxSlab, slab ) ) {
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
    // we can proceed.
    // The lower bound has a +1 because FindLinks looks one slab back
    for (int j=-2*GROUP_RADIUS+1; j<=2*GROUP_RADIUS; j++){
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
    
    if (Kick.notdone(slab)) return 0;  // Must have accelerations
    // note that this condition only has any effect if group finding is turned off
    
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
        int zstride = PP->cpd /omp_get_max_threads();
        int ystride = PP->cpd /omp_get_max_threads();
        int minstride = 12;
        if (ystride < minstride) ystride = minstride;
        if (zstride < minstride) zstride = minstride;
        int cpd = PP->cpd;
        STDLOG(1,"Binning particles for slab %d\n",slab);
        #pragma omp parallel for schedule(dynamic,ystride)
        for (int y=0;y<cpd;y++) {
            for (int z=0;z<cpd;z++) {
                Cell c = PP->GetCell(slab, y, z);
                tsc(c.pos,PP->CellCenter(slab,y,z),density,c.count(),P.PowerSpectrumN1d,1.0);
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
    LBW->DeAllocate(AccSlab,slab);

    return;
    MicrostepCPU.Start();
    // Do microstepping here
    if(MicrostepEpochs != NULL){
        STDLOG(1,"Beginning microsteps for slab %d\n", slab);
        MicrostepControl *MC = new MicrostepControl;
        MC->setup(GFC->globalslabs[slab], *MicrostepEpochs, P.MicrostepTimeStep, JJ->eps);
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
    // This is blocking because it uses the LoadIC module, not LBW
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
        DriftAndCopy2InsertList(slab, driftfactor, DriftCell);
    }
    
    // We freed AccSlab in Microstep to save space
    if (GFC == NULL){
	    // We kept the accelerations until here because of third-order LPT
	    if (P.StoreForces && !P.ForceOutputDebug) {
	        STDLOG(1,"Storing Forces in slab %d\n", slab);
	        LBW->StoreArenaBlocking(AccSlab,slab);
	    }
	    else{
	        LBW->DeAllocate(AccSlab,slab);
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
        LBW->DeAllocate(VelLPTSlab, slab);
    
    // Gather particles from the insert list and make the merge slabs
    uint64 n_merge = FillMergeSlab(slab);
    merged_particles += n_merge;

    FinishFreeSlabs.Start();
    
    // This may be the last time be need any of the CellInfo slabs that we just used
    // We can't immediately free CellInfo before NearForce might need it until we're FORCE_RADIUS away
    // An alternative to this would be to just wait for FORCE_RADIUS before finishing
    for(int j = -2*FORCE_RADIUS, consec = 0; j <= 2*FORCE_RADIUS; j++){
        if(Finish.done(slab + j) || j == 0)
            consec++;
        else
            consec = 0;
        if (consec >= 2*FORCE_RADIUS + 1)
            LBW->DeAllocate(CellInfoSlab, slab + j - FORCE_RADIUS);
    }
    
    // Now delete the original particles
    LBW->DeAllocate(PosSlab,slab);
    LBW->DeAllocate(VelSlab,slab);
    LBW->DeAllocate(AuxSlab,slab);

    FinishFreeSlabs.Stop();
    
    // Make the multipoles
    LBW->AllocateArena(MultipoleSlab,slab);
    ComputeMultipoleSlab(slab);
    
    // Write out the particles and multipoles and delete
    WriteMergeSlab.Start();
    LBW->StoreArenaNonBlocking(MergePosSlab,slab);
    LBW->StoreArenaNonBlocking(MergeVelSlab,slab);
    LBW->StoreArenaNonBlocking(MergeAuxSlab,slab);
    LBW->StoreArenaNonBlocking(MergeCellInfoSlab,slab);
    WriteMergeSlab.Stop();

    WriteMultipoleSlab.Start();
    LBW->StoreArenaNonBlocking(MultipoleSlab,slab);
    WriteMultipoleSlab.Stop();

    int pwidth = FetchSlabs.number_of_slabs_executed - Finish.number_of_slabs_executed;
    STDLOG(1, "Current pipeline width (N_fetch - N_finish) is %d\n", pwidth);
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
    assertf(FORCE_RADIUS >= 0, "Illegal FORCE_RADIUS: %d\n", FORCE_RADIUS);
    assertf(GROUP_RADIUS >= 0, "Illegal GROUP_RADIUS: %d\n", GROUP_RADIUS); 
    STDLOG(0,"Adopting FORCE_RADIUS = %d\n", FORCE_RADIUS);
    STDLOG(0,"Adopting GROUP_RADIUS = %d\n", GROUP_RADIUS);

    int cpd = P.cpd;
    int first = 0;  // First slab to load
    STDLOG(1,"First slab to load will be %d\n", first);

        FetchSlabs.instantiate(cpd, first, &FetchSlabPrecondition,          &FetchSlabAction         );
      TransposePos.instantiate(cpd, first, &TransposePosPrecondition,       &TransposePosAction      );
         NearForce.instantiate(cpd, first + FORCE_RADIUS, &NearForcePrecondition,          &NearForceAction         );
       TaylorForce.instantiate(cpd, first + FORCE_RADIUS, &TaylorForcePrecondition,        &TaylorForceAction       );
              Kick.instantiate(cpd, first + FORCE_RADIUS, &KickPrecondition,               &KickAction              );
            Output.instantiate(cpd, first + FORCE_RADIUS + 2*GROUP_RADIUS, &OutputPrecondition,             &OutputAction            );
             Drift.instantiate(cpd, first + FORCE_RADIUS + 2*GROUP_RADIUS, &DriftPrecondition,              &DriftAction             );
            Finish.instantiate(cpd, first + FORCE_RADIUS + 2*GROUP_RADIUS + FINISH_WAIT_RADIUS, &FinishPrecondition,             &FinishAction            );
            
    // If group finding is disabled, we can make the dependencies no-ops so they don't hold up the pipeline
    if(GFC != NULL){
        MakeCellGroups.instantiate(cpd, first + FORCE_RADIUS, &MakeCellGroupsPrecondition,     &MakeCellGroupsAction    );
    FindCellGroupLinks.instantiate(cpd, first + FORCE_RADIUS + 1, &FindCellGroupLinksPrecondition, &FindCellGroupLinksAction);
        DoGlobalGroups.instantiate(cpd, first + FORCE_RADIUS + 2*GROUP_RADIUS, &DoGlobalGroupsPrecondition,     &DoGlobalGroupsAction    );
             Microstep.instantiate(cpd, first + FORCE_RADIUS + 2*GROUP_RADIUS, &MicrostepPrecondition,          &MicrostepAction         );
          FinishGroups.instantiate(cpd, first + FORCE_RADIUS + 2*GROUP_RADIUS, &FinishGroupsPrecondition,       &FinishGroupsAction      );
    } else {
        MakeCellGroups.instantiate(cpd, first, &NoopPrecondition, &NoopAction );
    FindCellGroupLinks.instantiate(cpd, first, &NoopPrecondition, &NoopAction );
        DoGlobalGroups.instantiate(cpd, first, &NoopPrecondition, &NoopAction );
             Microstep.instantiate(cpd, first, &NoopPrecondition, &NoopAction );
          FinishGroups.instantiate(cpd, first, &NoopPrecondition, &NoopAction );
    }
           
    if(WriteState.Do2LPTVelocityRereading)
        LPTVelocityReRead.instantiate(cpd, first + FORCE_RADIUS + 2*GROUP_RADIUS - FINISH_WAIT_RADIUS,
                                          &FetchLPTVelPrecondition,   &FetchLPTVelAction   );
    else
        LPTVelocityReRead.instantiate(cpd, first, &NoopPrecondition, &NoopAction );

    while( !Finish.alldone() ) {
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
    }

    if(IL->length!=0)
        IL->DumpParticles();
    
    assertf(IL->length==0, 
        "Insert List not empty (%d) at the end of timestep().  Time step too big?\n", IL->length);
    
    assertf(merged_particles == P.np, "Merged slabs contain %d particles instead of %d!\n", merged_particles, P.np);

    uint64 total_n_output = n_output;
    if(GFC != NULL)
        total_n_output += GFC->n_L0_output;
    
    if(ReadState.DoTimeSliceOutput)
        assertf(total_n_output == P.np, "TimeSlice output contains %d particles instead of %d!\n", total_n_output, P.np);

    STDLOG(1,"Completing timestep()\n");
    TimeStepWallClock.Stop();
}


// ===================================================================

uint64 NP_from_IC = 0;

int FetchICPrecondition(int slab) {
    // We always do this.
    return 1;
}
void FetchICAction(int slab) {
    STDLOG(1,"Fetching slab %d\n", slab);
    // Get a slab of particles and put them on the InsertList
    NP_from_IC += LoadSlab2IL(slab);
    
    // We also need to create a null slab
    LBW->AllocateSpecificSize(PosSlab,slab, 0);
    LBW->AllocateSpecificSize(VelSlab,slab, 0);
    LBW->AllocateSpecificSize(AuxSlab,slab, 0);
    LBW->AllocateArena(CellInfoSlab,slab);
    int cpd = PP->cpd;
    for (int y=0; y<cpd; y++)
        for (int z=0; z<cpd; z++) {
            PP->CellInfo(slab,y,z)->makenull();
        }
    return;
}
/*
 * Registers the preconditions and actions for an IC step
 * When doing IC loading, we require that the neighboring slabs be loaded
 * just to be sure that no particles have crossed the boundary.  This is trivial
 * if we overload the Drift dependency with the FetchIC condition/actions.
 */

void timestepIC(void) {
    STDLOG(0,"Initiating timestepIC()\n");
    TimeStepWallClock.Clear();
    TimeStepWallClock.Start();
    
    FORCE_RADIUS = 0;  // so we know when we can free CellInfo in Finish
    GROUP_RADIUS = 0;

    int cpd = P.cpd; int first = 0;
    Drift.instantiate(cpd, first, &FetchICPrecondition, &FetchICAction );
    Finish.instantiate(cpd, first + FINISH_WAIT_RADIUS,  &FinishPrecondition,  &FinishAction );

    while( !Finish.alldone() ) {
        Drift.Attempt();
       Finish.Attempt();
    }

    assertf(NP_from_IC == P.np, "Expected to read a total of %llu particles from IC files, but only read %llu.\n", P.np, NP_from_IC);
    assertf(merged_particles == P.np, "Merged slabs contain %d particles instead of %d!\n", merged_particles, P.np);
    
    char filename[1024];
    sprintf(filename,"%s/slabsize",P.WriteStateDirectory);
    Slab->write(filename);
    assertf(IL->length==0,
        "Insert List not empty (%d) at the end of timestep().  Particles in IC files not sufficiently sorted?\n", IL->length);

    STDLOG(1,"Completing timestepIC()\n");
    TimeStepWallClock.Stop();
}

// ===================================================================
// Multipole recovery mode

int FetchPosSlabPrecondition(int slab) {
    if(LBW->total_allocation > .5*P.MAXRAMMB*1024LLU*1024LLU){
        // Are we spinning because we need more RAM?
        Dependency::NotifySpinning(NOT_ENOUGH_RAM);
        return 0;
    }
    return 1;
}

void FetchPosSlabAction(int slab) {
    STDLOG(0,"Fetching slab %d with %d particles\n", slab, Slab->size(slab));
    // Load all of the particle files together
    LBW->LoadArenaNonBlocking(MergeCellInfoSlab,slab);
    LBW->LoadArenaNonBlocking(MergePosSlab,slab);  // Load directly into the merge slabs
    assertf(Slab->size(slab)*sizeof(posstruct)<=
        fsize(LBW->ReadSlabDescriptorName(MergePosSlab,slab).c_str()),
        "PosSlab size doesn't match prediction\n");
    /*LBW->LoadArenaNonBlocking(AuxSlab,slab);
    assertf(Slab->size(slab)*sizeof(auxstruct)<=
        fsize(LBW->ReadSlabDescriptorName(AuxSlab,slab).c_str()),
        "AuxSlab size doesn't match prediction\n");*/
}

int FinishMultipolesPrecondition(int slab) {
    if( !LBW->IOCompleted( MergePosSlab,      slab )
        || !LBW->IOCompleted( MergeCellInfoSlab, slab )
        //|| !LBW->IOCompleted( AuxSlab,      slab )
        ) return 0;
    return 1;
}

void FinishMultipolesAction(int slab) {
    STDLOG(1,"Finishing multipole slab %d\n", slab);
        
    // Make the multipoles
    LBW->AllocateArena(MultipoleSlab,slab);
    ComputeMultipoleSlab(slab);
    
    WriteMultipoleSlab.Start();
    LBW->StoreArenaNonBlocking(MultipoleSlab,slab);
    WriteMultipoleSlab.Stop();
    
    LBW->DeAllocate(MergePosSlab,slab);
    LBW->DeAllocate(MergeCellInfoSlab,slab);
}


void timestepMultipoles(void) {
    STDLOG(0,"Initiating timestepMultipoles()\n");
    TimeStepWallClock.Clear();
    TimeStepWallClock.Start();
    
    FORCE_RADIUS = 0;  // so we know when we can free CellInfo in Finish
    GROUP_RADIUS = 0;

    int cpd = P.cpd; int first = 0;
    FetchSlabs.instantiate(cpd, first, &FetchPosSlabPrecondition, &FetchPosSlabAction );
    Finish.instantiate(cpd, first,  &FinishMultipolesPrecondition,  &FinishMultipolesAction );

    while( !Finish.alldone() ) {
        FetchSlabs.Attempt();
            Finish.Attempt();
    }

    STDLOG(1,"Completing timestepMultipoles()\n");
    TimeStepWallClock.Stop();
}

// =========================================================
// IO Benchmark mode

int FinishBenchmarkIOPrecondition(int slab) {
    // Wait for everything to be read
    if( !LBW->IOCompleted( CellInfoSlab,      slab )
        || !LBW->IOCompleted( PosSlab, slab )
        || !LBW->IOCompleted( VelSlab, slab )
        || !LBW->IOCompleted( AuxSlab,      slab )
        || !LBW->IOCompleted( TaylorSlab,      slab )
        ) return 0;
    return 1;
}

void FinishBenchmarkIOAction(int slab) {
    STDLOG(1,"Finishing benchmark IO slab %d\n", slab);
        
    /*// Make the multipoles
    LBW->AllocateArena(MultipoleSlab,slab);
    ComputeMultipoleSlab(slab);
    
    WriteMultipoleSlab.Start();
    LBW->StoreArenaNonBlocking(MultipoleSlab,slab);
    WriteMultipoleSlab.Stop();
    
    LBW->DeAllocate(MergePosSlab,slab);
    LBW->DeAllocate(MergeCellInfoSlab,slab);*/

    LBW->WriteArena(CellInfoSlab, slab, IO_DELETE, IO_NONBLOCKING, 
                    LBW->WriteSlabDescriptorName(MergeCellInfoSlab,slab).c_str());
    LBW->WriteArena(PosSlab, slab, IO_DELETE, IO_NONBLOCKING, 
                    LBW->WriteSlabDescriptorName(MergePosSlab,slab).c_str());
    LBW->WriteArena(VelSlab, slab, IO_DELETE, IO_NONBLOCKING, 
                    LBW->WriteSlabDescriptorName(MergeVelSlab,slab).c_str());
    LBW->WriteArena(AuxSlab, slab, IO_DELETE, IO_NONBLOCKING, 
                    LBW->WriteSlabDescriptorName(MergeAuxSlab,slab).c_str());
    LBW->WriteArena(TaylorSlab, slab, IO_DELETE, IO_NONBLOCKING, 
                    LBW->WriteSlabDescriptorName(MultipoleSlab,slab).c_str());
}

void timestepBenchmarkIO(int nslabs) {
    // We want to read slabs from the read directory and write them to the write directory, probably without modification
    // We can probably reuse the main FetchSlabs depdendency and write a new Finish dependency
    // One may not want to have to read and write all the slabs for a large box, so `nslabs` can be specified to use fewer

    STDLOG(0,"Initiating timestepBenchmarkIO()\n");
    TimeStepWallClock.Clear();
    TimeStepWallClock.Start();
    
    FORCE_RADIUS = 0;
    GROUP_RADIUS = 0;

    int cpd = P.cpd; int first = 0;
    assertf(nslabs <= cpd, "nslabs (%d) cannot be larger than cpd (%d)\n", nslabs, cpd);
    if (nslabs <= 0)
        nslabs = cpd;

    // Use the Kick as finish because FetchSlabs fetches FETCHAHEAD past the kick
    FetchSlabs.instantiate(nslabs, first, &FetchSlabPrecondition, &FetchSlabAction );
    Kick.instantiate(nslabs, first,  &FinishBenchmarkIOPrecondition,  &FinishBenchmarkIOAction );

    while( !Kick.alldone() ) {
        FetchSlabs.Attempt();
              Kick.Attempt();
    }

    STDLOG(1,"Completing timestepBenchmarkIO()\n");
    TimeStepWallClock.Stop();
}

// =============================================================================================== //

#include "read_pack14.cpp"

const char* StandaloneFOF_slice_dir;
int StandaloneFOFLoadSlabPrecondition(int slab) {
    if(LBW->total_allocation > .5*P.MAXRAMMB*1024LLU*1024LLU){
        Dependency::NotifySpinning(NOT_ENOUGH_RAM);
        return 0;
    }
    return 1;
}

void StandaloneFOFLoadSlabAction(int slab) {
    char fname[1024];
    // TODO: Add support for L0 slabs?
    sprintf(fname, "%s/%s.z%5.3f.slab%04d.dat", StandaloneFOF_slice_dir, P.SimName, ReadState.Redshift, slab);
    STDLOG(1,"Load Slab %d from \"%s\"\n", slab, fname);

    size_t s = fsize(fname);
    LBW->AllocateSpecificSize(TimeSlice, slab, s);
    // We will read the raw pack14 asynchronously with LBW
    // then unpack it in a separate dependency
    // TODO: support states as well as time slices
    LBW->ReadArena(TimeSlice, slab, IO_NONBLOCKING, fname);
}

int StandaloneFOFUnpackSlabPrecondition(int slab) {
    if (! LBW->IOCompleted(TimeSlice, slab)) return 0;
    return 1;
}

void StandaloneFOFUnpackSlabAction(int slab) {
    printf("Unpacking slab %d\n", slab);
    STDLOG(1, "Unpacking slab %d\n", slab);
    int nump = unpack_slab_pack14(slab, P.HaloTaggableFraction);
    STDLOG(1,"Found %d particles in slab %d\n", nump, slab);

    LBW->DeAllocate(TimeSlice, slab);
}

int StandaloneFOFMakeCellGroupsPrecondition(int slab) {
    if (TransposePos.notdone(slab)) return 0;
    return 1;
}

int StandaloneFOFFinishPrecondition(int slab) {
    if (DoGlobalGroups.notdone(slab)) return 0;
    return 1;
}

void StandaloneFOFFinishAction(int slab) {
    STDLOG(1,"Deleting slab %d\n", slab);

    // Release the group-local copies of the particles
    GlobalGroupSlab *GGS = GFC->globalslabs[slab];
    delete GGS;
    GFC->globalslabs[slab] = NULL;

    LBW->DeAllocate(PosSlab, slab);
    LBW->DeAllocate(VelSlab, slab);
    LBW->DeAllocate(AuxSlab, slab);
    LBW->DeAllocate(CellInfoSlab, slab);
}


void timestepStandaloneFOF(const char* slice_dir) {
    STDLOG(0,"Initiating timestepStandaloneFOF()\n");
    TimeStepWallClock.Clear();
    TimeStepWallClock.Start();

    int cpd = GFC->cpd;

    StandaloneFOF_slice_dir = slice_dir;

    FORCE_RADIUS = 0;
    GROUP_RADIUS = P.GroupRadius;
    assertf(GROUP_RADIUS >= 0, "Illegal GROUP_RADIUS: %d\n", GROUP_RADIUS); 
    STDLOG(0,"Adopting GROUP_RADIUS = %d\n", GROUP_RADIUS);

    int first = 0;
            FetchSlabs.instantiate(cpd, first, &StandaloneFOFLoadSlabPrecondition, &StandaloneFOFLoadSlabAction);
          TransposePos.instantiate(cpd, first, &StandaloneFOFUnpackSlabPrecondition, &StandaloneFOFUnpackSlabAction);
        MakeCellGroups.instantiate(cpd, first, &StandaloneFOFMakeCellGroupsPrecondition, &MakeCellGroupsAction);
    FindCellGroupLinks.instantiate(cpd, first + 1, &FindCellGroupLinksPrecondition, &FindCellGroupLinksAction);
        DoGlobalGroups.instantiate(cpd, first + 2*GFC->GroupRadius, &DoGlobalGroupsPrecondition, &DoGlobalGroupsAction);
                Finish.instantiate(cpd, first + 2*GFC->GroupRadius, &StandaloneFOFFinishPrecondition, &StandaloneFOFFinishAction);

    while (!Finish.alldone()) {
        FetchSlabs.Attempt();
        TransposePos.Attempt();
        MakeCellGroups.Attempt();
        FindCellGroupLinks.Attempt();
        DoGlobalGroups.Attempt();
        Finish.Attempt();
    }

    TimeStepWallClock.Stop();
}
