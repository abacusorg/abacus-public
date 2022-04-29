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

//#define WAIT_MULTIPOLES 1


int FORCE_RADIUS = -1;
int GROUP_RADIUS = -1;

//#define FETCHAHEAD (2*GROUP_RADIUS + FORCE_RADIUS + FINISH_WAIT_RADIUS + 2)
//#define FETCHAHEAD (2*GROUP_RADIUS + FORCE_RADIUS + 2)
#define FETCHAHEAD (max(2*GROUP_RADIUS + FORCE_RADIUS, FINISH_WAIT_RADIUS) + 2)
#define FETCHPERSTEP 1
// Recall that all of these Dependencies have a built-in STimer
// to measure the amount of time spent on Actions.
Dependency FetchSlabs;
Dependency TransposePos;
Dependency NearForce;
Dependency TaylorTranspose;  // 2D
Dependency TaylorForce;
Dependency Kick;
Dependency MakeCellGroups;
Dependency FindCellGroupLinks;
Dependency DoGlobalGroups;
Dependency Output;
Dependency Microstep;
Dependency FinishGroups;
Dependency Drift;
Dependency NeighborSend;  // 2D
Dependency FinishParticles;
Dependency FinishMultipoles;
Dependency CheckForMultipoles; //only for parallel case. otherwise NOOP. 
Dependency UnpackLPTVelocity;

// TODO: should we consider de-coupling PARALLEL from the concept of a merged convolve/singlestep?
#ifdef PARALLEL
#include "ConvolutionParametersStatistics.cpp"
#include "InCoreConvolution.cpp"
#include "ParallelConvolution.cpp"
STimer ConvolutionWallClock;
STimer BarrierWallClock;
STimer MultipoleTransferCheck; 
STimer TaylorTransferCheck;
#endif

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
    // We want to read ahead enough that we are reading while computing forces
    // i.e. we would like to read a few slabs ahead of the Kick.
    // But if the steps after Kick are slow enough that we aren't Finishing promptly,
    // we can get a buildup of slabs and run out of memory. So we tie to the Drift instead.
    if(FetchSlabs.wrap(slab - Drift.last_slab_executed) > FETCHAHEAD
        && FetchSlabs.done(Drift.last_slab_executed))
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
    STDLOG(1,"Fetching slab %d with %d particles (%d ghost)\n",
        slab, SS->size(slab), SS->size_with_ghost(slab));

    // Load all of the particle files together
    SB->LoadArenaNonBlocking(CellInfoSlab,slab);
    SB->LoadArenaNonBlocking(PosSlab,slab);

    if(WriteState.Do2LPTVelocityRereading && UnpackLPTVelocity.notdone(slab)){
        ICFile::FromFormat(P.ICFormat, slab)->read_vel_nonblocking();
    }

    // Don't bother to load the vel/aux/taylors for slabs that won't be kicked until the wrap
    //#ifndef PARALLEL
    // LHG: only need this when we're memory starved
    #if 0
    if(FetchSlabs.number_of_slabs_executed < FORCE_RADIUS)
        return;
    #endif

#ifdef PARALLEL
    // SB->AllocateArena(TaylorSlab, slab + FORCE_RADIUS, RAMDISK_NO);
// 	ParallelConvolveDriver->RecvTaylorSlab(slab + FORCE_RADIUS);
// 	STDLOG(2, "Received Taylor slab via MPI%d\n", slab + FORCE_RADIUS);

#else
    SB->LoadArenaNonBlocking(TaylorSlab,slab);
#endif

    SB->LoadArenaNonBlocking(VelSlab, slab);
    SB->LoadArenaNonBlocking(AuxSlab, slab);
}

// -----------------------------------------------------------------

int TransposePosPrecondition(int slab){
    // TODO: technically need separate WAITING_FOR_IO flags for each slab, otherwise "reasons" can cross-talk
    if(!SB->IsIOCompleted(PosSlab, slab)){
        if(SB->IsSlabPresent(PosSlab, slab))
            Dependency::NotifySpinning(WAITING_FOR_IO);
        return 0;
    }

    if(!SB->IsIOCompleted(CellInfoSlab, slab)){
        if(SB->IsSlabPresent(CellInfoSlab, slab))
            Dependency::NotifySpinning(WAITING_FOR_IO);
        return 0;
    }

        return 1;
}

void TransposePosAction(int slab){
    SB->AllocateArena(PosXYZSlab, slab);
    int cpd = P.cpd;

    NUMA_FOR(y,0,cpd)
        // Execute the transpose in pencils
        posstruct *pos = CP->PosCell(slab, y, node_z_start_ghost);
        List3<FLOAT> posxyz = CP->PosXYZCell(slab, y, node_z_start_ghost);
        uint64 Npencil = CP->PencilLenWithGhost(slab,y);

        for(uint64 i = 0; i < Npencil; i++){
            posxyz.X[i] = pos[i].x;
            posxyz.Y[i] = pos[i].y;
            posxyz.Z[i] = pos[i].z;
        }
    }
}


// -----------------------------------------------------------------

int NearForcePrecondition(int slab) {
    for(int i=-FORCE_RADIUS;i<=FORCE_RADIUS;i++){
        // Technically, I think we only need the CellInfo to construct pencils
        // But it's convenient to have pos so the GPU can immediately execute any pencil
        if(TransposePos.notdone(slab+i))
            return 0;
        if( !SB->IsIOCompleted( CellInfoSlab, slab+i ) ){
            if(SB->IsSlabPresent(CellInfoSlab, slab+i))
                Dependency::NotifySpinning(WAITING_FOR_IO);
            return 0;
        }
    }

    return 1;
}

void NearForceAction(int slab) {
    // Do some data checks
    assertf(are_cellinfo_legal(slab, SS->size(slab), SS->size_with_ghost(slab)),
            "Cell info of slab %d contain out of bounds data\n", slab);
    // Could also check that the sum of the cell counts add up to SS->size(slab);

    SlabForceTime[slab].Start();

    int blocking = 0;
    NFD->ExecuteSlab(slab, blocking);

    SlabForceLatency[slab].Start();
    
    // Busy-wait for all GPU work for this slab to finish
    // while(!NFD->SlabDone(slab)) ;
}

// -----------------------------------------------------------------

int TaylorTransposePrecondition(int slab){
    // TODO: install in manifest
#ifdef PARALLEL
    TaylorTransferCheck.Start(); 
    int ready = ParallelConvolveDriver->CheckTaylorSlabReady(slab);
    TaylorTransferCheck.Stop();
    if(!ready) {
        if(SB->IsSlabPresent(TaylorSlab, slab)) Dependency::NotifySpinning(WAITING_FOR_MPI);
        return 0;
    }
#endif

    return 1;
}

void TaylorTranposeAction(int slab){
    if(MPI_size_z > 1){
        // FFT and launch MPI All-to-all
        TY->ComputeIFFTZAndMPI(slab);
    }

    // for the 1D parallel code, this is a no-op
}

// -----------------------------------------------------------------

int TaylorForcePrecondition(int slab) {
    if(TaylorTranspose.notdone(slab) ||  // 1D & 2D
        !TY->IsMPIDone(slab)){  // 2D
        return 0;
    }

    if( !SB->IsIOCompleted( CellInfoSlab, slab ) ){
        if(SB->IsSlabPresent(CellInfoSlab, slab))
            Dependency::NotifySpinning(WAITING_FOR_IO);
        return 0;
    }
    if( !SB->IsIOCompleted( PosSlab, slab ) ){
        if(SB->IsSlabPresent(PosSlab, slab))
            Dependency::NotifySpinning(WAITING_FOR_IO);
        return 0;
    }

#ifndef PARALLEL
    if( !SB->IsIOCompleted( TaylorSlab, slab ) ){
        if(SB->IsSlabPresent(TaylorSlab, slab))
            Dependency::NotifySpinning(WAITING_FOR_IO);
        return 0;
    }
#endif

    return 1;
}

void TaylorForceAction(int slab) {
	MTCOMPLEX *t = (MTCOMPLEX *) SB->GetSlabPtr(TaylorSlab, slab);

    SlabFarForceTime[slab].Start();
    SB->AllocateArena(FarAccSlab, slab);

    TaylorCompute.Start();
    ComputeTaylorForce(slab);
    TaylorCompute.Stop();

    if(P.StoreForces == 2){
        // We want to output the FarAccSlab to the FarAcc file.
        // This must be a blocking write.
        SB->WriteArena(FarAccSlab, slab, IO_KEEP, IO_BLOCKING);
    }

#ifdef PARALLEL
    // If parallel, the Taylors came via MPI and don't have a corresponding slab file
    int delete_taylors_file = 0;
#else
    // Deallocate and delete the underlying file if we're overwriting
    int delete_taylors_file = WriteState.OverwriteConvState;
#endif

    SB->DeAllocate(TaylorSlab, slab, delete_taylors_file);
    SlabFarForceTime[slab].Stop();
}


// -----------------------------------------------------------------

int KickPrecondition(int slab) {
    if( !SB->IsIOCompleted( VelSlab, slab ) ){
        if(SB->IsSlabPresent(VelSlab, slab))
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
    
    // computing CPU forces must be done before the pos releases below
    if(!P.ForceCPU && P.ForceOutputDebug)
        NFD->CheckGPUCPU(slab);

    // Release the trailing slab if it won't be needed at the wrap
    // Technically we could release it anyway and re-do the transpose from PosSlab,
    // but if we're not doing group finding we may have already written and released PosSlab
    if(Kick.raw_number_executed >= 2*FORCE_RADIUS)
        SB->DeAllocate(PosXYZSlab, slab - FORCE_RADIUS);

    // Special case: if this is the last slab, free all +/- FORCE_RADIUS
    // Not worrying about this in the PARALLEL case; we have other non-destructions
    STDLOG(1,"%d slabs have been Kicked so far\n", Kick.raw_number_executed);
    if(Kick.raw_number_executed == total_slabs_on_node-1)
        for(int j = slab - FORCE_RADIUS+1; j <= slab + FORCE_RADIUS; j++)
            SB->DeAllocate(PosXYZSlab, j);

    //#ifndef PARALLEL
    #if 0
    // Queue up slabs near the wrap to be loaded again later
    // This way, we don't have idle slabs taking up memory while waiting for the pipeline to wrap around
    if(Kick.number_of_slabs_executed < FORCE_RADIUS){
        STDLOG(3,"Marking slab %d for repeat\n", slab - FORCE_RADIUS);
        TransposePos.mark_to_repeat(slab - FORCE_RADIUS);
        SB->DeAllocate(PosXYZSlab, slab - FORCE_RADIUS);
        // The first two won't need PosSlab until the second time around
        //SB->DeAllocate(PosSlab, slab - FORCE_RADIUS);
        SB->DeAllocate(CellInfoSlab, slab - FORCE_RADIUS);
        FetchSlabs.mark_to_repeat(slab - FORCE_RADIUS);
    }
    #endif
    
    // Accumulate stats from the SICs and release them
    NFD->Finalize(slab);
    
    if (P.StoreForces == 2) {
        // We want to output the AccSlab to the NearAcc file.
        // This must be a blocking write.

#ifdef DIRECTSINGLESPLINE
        // Single spline requires a prefactor multiplication, which we defer to the kick for efficiency
        // But analysis routines that use ForceOutputDebug, like Ewald, expect this prefactor to already be applied
        // So apply it here, storing the original in a temporary copy
        uint64 npslab = SS->size(slab);
        accstruct *nearacctmp = new accstruct[npslab];
        accstruct *nearacc = (accstruct *) SB->GetSlabPtr(AccSlab, slab);
        memcpy(nearacctmp, nearacc, npslab*sizeof(accstruct));
        FLOAT inv_eps3 = 1./(NFD->SofteningLengthInternal*NFD->SofteningLengthInternal*NFD->SofteningLengthInternal);
        #pragma omp parallel for schedule(static)
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

    ReleaseFreeMemoryToKernel();
}

// -----------------------------------------------------------------

int MakeCellGroupsPrecondition(int slab) {
    // Only PosXYZ is used for sources, so we're free to rearrange PosSlab
    // in group finding after the transpose
    if( TransposePos.notdone(slab) ) return 0;

    // Need the new velocities because we're going to rearrange particles
    if( Kick.notdone(slab) ) return 0;

    // Also need the auxs, because we're going to re-order
    if( !SB->IsIOCompleted( AuxSlab, slab ) ){
        if(SB->IsSlabPresent(AuxSlab, slab))
            Dependency::NotifySpinning(WAITING_FOR_IO);
        return 0;
    }
    return 1;
}

void MakeCellGroupsAction(int slab) {
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
	FindGroupLinks(slab);
}

// -----------------------------------------------------------------

int DoGlobalGroupsPrecondition(int slab) {
#ifdef ONE_SIDED_GROUP_FINDING
    /* We're going to search for GlobalGroups that include cells from
    slabs [slab,slab+2*GroupRadius].  However, this also includes the idea
    that we will *not* find groups if they include anything in slab-1.
    So we must have GroupLinks between [slab-1,slab] as well.  Moreover,
    we must go all the way up to [slab+2*GR-1,slab+2*GR], as there may
    be groups that are only now becoming eligible.
    */

    if (Kick.notdone(slab)) return 0;
    for (int j=0; j<=2*GROUP_RADIUS; j++){
        if (FindCellGroupLinks.notdone(slab+j)) return 0;
    }

#else
    // We're going to close all CellGroups in this slab.
    // GlobalGroups can span 2*GroupRadius+1.
    // But even though we usually encounter a CellGroup in its minimum slab,
    // we could be anywhere in the first instance.  So we have to query a big range.
    // That said, if the nearby slab has already closed global groups, then
    // we can proceed.  This particularly matters in the parallel version, where
    // we may already have closed groups in higher numbered slabs.

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
#endif
    return 1;
}

void DoGlobalGroupsAction(int slab) {
    FindAndProcessGlobalGroups(slab);

    // The first 2*GroupRadius times we get here, we can attempt to free
    // info from slab.  The Manifest code sends everything <S, so we need S=slab+1
    #ifdef ONE_SIDED_GROUP_FINDING
        if (DoGlobalGroups.raw_number_executed<2*GROUP_RADIUS) {
            SendManifest->QueueToSend(slab+1);
            SendManifest++;
        }
    #endif
    STDLOG(0,"Exiting Find Global Groups action in slab %d\n", slab);

}

// -----------------------------------------------------------------

int MicrostepPrecondition(int slab){
    // We are going to second-half kick this slab
    if (DoGlobalGroups.notdone(slab))
        return 0;
    return 1;
}

void MicrostepAction(int slab){

    // TODO: This is now not the place to do this.
    // All kicks (and half-unkicks) for output are done; discard accels.
    // We de-allocate in Drift if we aren't doing group finding
    // SB->DeAllocate(AccSlab,slab);

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

	//Dependency do_action() assumes that each dependency processes all particles in a given slab.
	//Microstepping is an exception; it only does the group particles! Correct the bookkeeping here.
	Microstep.num_particles += GFC->globalslabs[slab]->np - SS->size(slab);

    MicrostepCPU.Stop();
}

// -----------------------------------------------------------------

int FinishGroupsPrecondition(int slab){
    // Is the asychronous GPU microstepping done?
    //if (!GFC->microstepcontrol[slab]->GPUGroupsDone()) return 0

    // We are going to release these groups.
    // TODO: If Microstep changes, this may change.  At present, we really need Output to be done
    // because this step triggers the writeback of the new Pos/Vel.
    if (Microstep.notdone(slab)) return 0;

    return 1;
}

void FinishGroupsAction(int slab){
    // Scatter pos,vel updates to slabs, and release GGS
    FinishGlobalGroups(slab);   // This will Scatter Pos/Vel
    delete GFC->microstepcontrol[slab];
    GFC->microstepcontrol[slab] = NULL;

    // TODO: When we're ready to send Group-based Manifests, it would go here.
    // Would pass slab+1 to the manifest code as the faux finished slab.
}

// -----------------------------------------------------------------
/*
 * Checks if we are ready to do all outputs for this step
 * Anything that modifies the particles at the current time should happen before here
 * Importantly: the OutputAction is only for non-L0 particles.
 * L0 particle outputs need to happen in the DoGlobalGroupsAction(),
 * before Microstepping.
 */
int OutputPrecondition(int slab) {

    #ifdef ONE_SIDED_GROUP_FINDING
    /* This must wait until all groups including the slab have been found.
    It used to be that closing groups on the current slab S would do this,
    but now groups are only found looking upwards from S, so we need to have
    closed groups on all slabs from S to S-2*GroupRadius, inclusive.
    */
    for (int s=0; s<=2*GROUP_RADIUS; s++)
        if (FinishGroups.notdone(slab-s)) return 0;
    // Must have found groups to be able to output light cones
    // note that group outputs were already done
    #else
        if (FinishGroups.notdone(slab)) return 0;
    #endif

    // Note the following conditions only have any effect if group finding is turned off

    if (Kick.notdone(slab)) return 0;  // Must have accelerations

    // Also obviously need the aux!  This was checked in CellGroups,
    // but that may have been skipped if there's no group finding.
    if( !SB->IsIOCompleted( AuxSlab, slab ) ){
        if(SB->IsSlabPresent(AuxSlab, slab))
            Dependency::NotifySpinning(WAITING_FOR_IO);
        return 0;
    }

    return 1;
}

uint64 n_output = 0;
void OutputAction(int slab) {

    // We are finally done with the groups for this slab.
    // Delete the Cell Groups.
    if (GFC!=NULL) GFC->DestroyCellGroups(slab);

    int step = WriteState.FullStepNumber;
    if (LPTStepNumber()>0) return;
    // Some output might want to be skipped during an IC step,
    // e.g., no light cones

    OutputTimeSlice.Start();


    // Having found all groups, we should output the Non-L0 (i.e., field) Taggable subsample.
    if(ReadState.DoSubsampleOutput) {
        assertf(ReadState.DoGroupFindingOutput == 1, "Subsample output should turn on group finding!\n"); // Currently Subsample Output requires GroupFinding Output.
        OutputNonL0Taggable(slab);
    }

    if (ReadState.DoTimeSliceOutput) {
        // If we are doing group finding, then we are doing group finding output and subsample output
        assertf(GFC == NULL || (ReadState.DoSubsampleOutput == 1 && ReadState.DoGroupFindingOutput == 1),
            "Preparing for timeslice output, expected either no group finding, or group finding and subsampling output!\n");

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
            // Start timing
            STimer lightConeTimer;
            lightConeTimer.Start();
            makeLightCone(slab,i);
            lightConeTimer.Stop();
            STDLOG(1, "LightCone %d for slab %d creation took %f seconds\n", i, slab, lightConeTimer.Elapsed());
        }
    }
    OutputLightCone.Stop();

    OutputBin.Start();
    if(ReadState.DoBinning){
        int ystride = CP->cpd /omp_get_max_threads();
        int minstride = 12;
        if (ystride < minstride) ystride = minstride;
        int cpd = CP->cpd;
        STDLOG(1,"Binning particles for slab %d\n",slab);
        #pragma omp parallel for schedule(dynamic,ystride)
        for (int y=0;y<cpd;y++) {
            for (int z = node_z_start; z < node_z_start + node_z_size; z++) {
                Cell c = CP->GetCell(slab, y, z);
                tsc(c.pos,CP->CellCenter(slab,y,z),density,c.count(),P.PowerSpectrumN1d,1.0);
            }
        }
    }
    OutputBin.Stop();
}

// -----------------------------------------------------------------
/*
 * Checks if we are ready to load the LPT velocities during an IC step.
 * Should not happen in normal execution
 */
int UnpackLPTVelocityPrecondition(int slab){
    if(FetchSlabs.notdone(slab))
        return 0;
    
    if(!ICFile::FromFormat(P.ICFormat, slab)->check_vel_read_done())
        return 0;

    return 1;
}

void UnpackLPTVelocityAction(int slab){
    // TODO: decide how to handle 2LPT velocity re-reading in the 2D code.
    // Do we need a neighbor exchange step, or is each node going to read +/- z?
    // SlabBuffer nominally only allows for one z, so we'd have to extend that
    unpack_ic_vel_slab(slab);
}

// -----------------------------------------------------------------
/* Checks if we are ready to apply a full-timestep drift to the particles
 * Any operation that relies on the positions and velocities being synced
 * should be checked for completion this year.
 */
int DriftPrecondition(int slab) {
    // We must have finished scattering into this slab
    // if (FinishGroups.notdone(slab)) return 0;

    // We will move the particles, so we must have done outputs
    if (Output.notdone(slab)) return 0;

    // We can't move particles until they've been used as gravity sources
    // However, we only use PosXYZSlab as sources, so we're free to move PosSlab

    // We also must have the 2LPT velocities
    // The finish radius is a good guess of how ordered the ICs are
    if(WriteState.Do2LPTVelocityRereading)
        for(int i=-FINISH_WAIT_RADIUS;i<=FINISH_WAIT_RADIUS;i++) 
            if (UnpackLPTVelocity.notdone(slab+i)) {
                return 0;
            }

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
        DriftPencilsAndCopy2InsertList(slab, driftfactor, DriftPencil);
    }

    // We freed AccSlab in Microstep to save space
    // if (GFC == NULL){
    // TODO: Remove that condition; we always can do this here.
	    // We kept the accelerations until here because of third-order LPT
	    if (P.StoreForces == 1) {
	        STDLOG(1,"Storing Forces in slab %d\n", slab);
	        SB->StoreArenaBlocking(AccSlab,slab);
	    }
	    else{
	        SB->DeAllocate(AccSlab,slab);
            ReleaseFreeMemoryToKernel();
	    }
	// }
}

// -----------------------------------------------------------------

int NeighborSendPrecondition(int slab){
    // Must have all particles drifted into this slab
    for(int j=-FINISH_WAIT_RADIUS;j<=FINISH_WAIT_RADIUS;j++) {
        if( Drift.notdone(slab+j) )
            return 0;
    }

    if (NeighborSend.raw_number_executed == total_slabs_on_node) return 0;

    return 1;
}

void NeighborSendAction(int slab){
    // Send the full list of neighbor particles for this slab
    DoNeighborSend(slab);  // in parallel.cpp
}

// -----------------------------------------------------------------

int FinishParticlesPrecondition(int slab) {
    for(int j=-FINISH_WAIT_RADIUS;j<=FINISH_WAIT_RADIUS;j++) {
        if( Drift.notdone(slab+j) )
            return 0;
    }

	if (FinishParticles.alldone(total_slabs_on_node)) return 0;

    if( !IsNeighborReceiveDone(slab) ){
        // This is an effective dependency on NeighborSend, because we won't
        // receive before sending.
        // We only need to receive 1 slab, not FWR, because the remote node
        // waits for FWR before sending.
        return 0;
    }

    return 1;
}


uint64 merged_particles = 0;
void FinishParticlesAction(int slab) {
	FinishPreamble.Start();

    SB->report_current();
    SB->report_peak();

    if (WriteState.Do2LPTVelocityRereading)
        SB->DeAllocate(VelLPTSlab, slab);
	FinishPreamble.Stop();

    // Gather particles from the insert list and make the merge slabs
    uint64 n_merge = FillMergeSlab(slab);
    merged_particles += n_merge;

	FinishPreamble.Start();

    // This may be the last time be need any of the CellInfo slabs that we just used
    // We can't immediately free CellInfo before NearForce might need it until we're FORCE_RADIUS away
    // An alternative to this would be to just wait for FORCE_RADIUS before finishing
    for(int j = -2*FORCE_RADIUS, consec = 0; j <= 2*FORCE_RADIUS; j++){
        if(FinishParticles.done(slab + j) || j == 0)
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

	STDLOG(2,"Done deallocing pos, vel, aux for slab %d\n", slab);
	FinishPreamble.Stop();

    // Make the multipoles
    ComputeMultipoleSlab(slab);

    // Write and free the merge particles
    WriteMergeSlab.Start();
    SB->StoreArenaNonBlocking(MergePosSlab,slab);
    SB->StoreArenaNonBlocking(MergeVelSlab,slab);
    SB->StoreArenaNonBlocking(MergeAuxSlab,slab);
    SB->StoreArenaNonBlocking(MergeCellInfoSlab,slab);
    WriteMergeSlab.Stop();

    #ifdef PARALLEL
    if (FinishParticles.raw_number_executed==0) SendManifest->QueueToSend(slab);
    #endif

}

// -----------------------------------------------------------------

int FinishMultipolesPrecondition(int slab){
    return FinishParticles.done(slab) && MF->IsMPIDone(slab);
}

void FinishMultipolesAction(int slab){
    // In the 2D code, the MPI transpose is now complete.
    // In the 1D code, this is just a continuation of FinishParticles.

    if(MPI_size_z > 1){
        STDLOG(1, "Executing multipoles z-FFT for slab %d\n", slab);
        MTCOMPLEX *slabptr = (MTCOMPLEX *) SB->AllocateArena(MultipoleSlab, slab, ramdisk_multipole_flag);
        MF->ComputeFFTZ(slab, slabptr);
    }

#ifdef PARALLEL
	QueueMultipoleMPI.Start();
    STDLOG(2, "Attempting to SendMultipoleSlab %d\n", slab);
    // distribute z's to appropriate nodes for this node's x domain.
 	ParallelConvolveDriver->SendMultipoleSlab(slab);
    // if we are finishing the first slab, set up receive MPI calls for incoming multipoles.
	if (FinishMultipoles.raw_number_executed==0){
		STDLOG(2, "Attempting to RecvMultipoleSlab %d\n", slab);
        // receive z's from other nodes for all x's.
		ParallelConvolveDriver->RecvMultipoleSlab(slab);
	}
	QueueMultipoleMPI.Stop();
#else
    // Write and free the multipoles
    WriteMultipoleSlab.Start();
    SB->StoreArenaNonBlocking(MultipoleSlab,slab);
    WriteMultipoleSlab.Stop();
#endif

    int pwidth = FetchSlabs.raw_number_executed - FinishMultipoles.raw_number_executed;
    STDLOG(1, "Current pipeline width (N_fetch - N_finish) is %d\n", pwidth);
    // release is cheap but not totally free, so might run every few Finishes
    //if (FinishMultipoles.raw_number_executed % 3 == 0)
        ReleaseFreeMemoryToKernel();
    ReportMemoryAllocatorStats();
}

// -----------------------------------------------------------------


#ifdef PARALLEL
int CheckForMultipolesPrecondition(int slab) {

    if( FinishMultipoles.notdone(slab) ) return 0;
	
	MultipoleTransferCheck.Start();
	int multipole_transfer_complete = ParallelConvolveDriver->CheckForMultipoleTransferComplete(slab);
	MultipoleTransferCheck.Stop();
	if (multipole_transfer_complete) return 1;
    else {
		if(SB->IsSlabPresent(MultipoleSlab, slab))
				Dependency::NotifySpinning(WAITING_FOR_MPI);
		return 0;
	}
}

void CheckForMultipolesAction(int slab) {
	SB->DeAllocate(MultipoleSlab, slab);
}

#endif
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

void AttemptReceiveManifest(){
    // If the manifest has been received, install it.
    if (ReceiveManifest->is_ready()) {
        ReceiveManifest->ImportData();
        ReceiveManifest++;
        STDLOG(1, "Readying the next Manifest, number %d\n", ReceiveManifest-_ReceiveManifest);
        ReceiveManifest->SetupToReceive();
    }
    return;
}

#define INSTANTIATE(dependency, first_relative) do { dependency.instantiate(nslabs, first + first_relative, &dependency##Precondition, &dependency##Action, #dependency); } while(0)
#define INSTANTIATE_NOOP(dependency, first_relative) do { dependency.instantiate(nslabs, first + first_relative, &NoopPrecondition, &NoopAction, ""); } while(0)


void InitializePipelineWidths(int MakeIC){
    FORCE_RADIUS = MakeIC ? 0 : P.NearFieldRadius;
    // The 2LPT pipeline is short (no group finding). We can afford to wait an extra slab to allow for large IC displacements
    FINISH_WAIT_RADIUS = (MakeIC || LPTStepNumber()) > 0 ? 2 : 1;
    assertf(FORCE_RADIUS >= 0, "Illegal FORCE_RADIUS: %d\n", FORCE_RADIUS);

    STDLOG(0,"Adopting FORCE_RADIUS = %d\n", FORCE_RADIUS);
    STDLOG(0,"Adopting FINISH_WAIT_RADIUS = %d\n", FINISH_WAIT_RADIUS);
}

// This happens much later, after outputs and group finding are planned
void InitializeGroupRadius(){
    GROUP_RADIUS = GFC != NULL ? P.GroupRadius : 0;
    assertf(GROUP_RADIUS >= 0, "Illegal GROUP_RADIUS: %d\n", GROUP_RADIUS);
    STDLOG(0,"Adopting GROUP_RADIUS = %d\n", GROUP_RADIUS);
}

void timestep(void) {

#ifdef PARALLEL
    ParallelConvolveDriver = new ParallelConvolution(P.cpd, P.order, P.MultipoleDirectory);
#endif

    TimeStepWallClock.Clear();  TimeStepWallClock.Start();
    STDLOG(1,"Initiating timestep()\n");

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

    // TODO: I'm not sure inflating FINISH_WAIT_RADIUS is the best way to deal with this
    // TODO: Also not sure this is the minimum number of slabs, even in that case
    // assertf(total_slabs_on_node >= 2*FINISH_WAIT_RADIUS + 1 + 2*FORCE_RADIUS + 4*GROUP_RADIUS, "Not enough slabs on node to finish any slabs!\n");
    int PAD = 0;
    assertf(total_slabs_on_node >= (2*GROUP_RADIUS + 1) + 2*FORCE_RADIUS + 1 + PAD, "Not enough slabs on node to close first group!\n");
    assertf(total_slabs_on_node >= 2*GROUP_RADIUS + FORCE_RADIUS + 2 * FINISH_WAIT_RADIUS + 1 + PAD, "Not enough slabs on node to finish any slabs!\n");
#endif

    int nslabs = P.cpd;
    int first = first_slab_on_node;  // First slab to load
    STDLOG(1,"First slab to load will be %d\n", first);

    #ifdef ONE_SIDED_GROUP_FINDING
        int first_outputslab = FORCE_RADIUS + 2*GROUP_RADIUS + (int)(GROUP_RADIUS > 0);
    #else
        int first_outputslab = FORCE_RADIUS + 2*GROUP_RADIUS;
    #endif
    if(LPTStepNumber() == 2){
        first_outputslab = max(FINISH_WAIT_RADIUS,first_outputslab);
    }


    INSTANTIATE(                  FetchSlabs, 0);
    INSTANTIATE(                TransposePos, 0);
    INSTANTIATE(                   NearForce, FORCE_RADIUS);
    INSTANTIATE(                 TaylorForce, FORCE_RADIUS);
    INSTANTIATE(                        Kick, FORCE_RADIUS);
    INSTANTIATE(                      Output, first_outputslab);
    INSTANTIATE(                       Drift, first_outputslab);
    INSTANTIATE(             FinishParticles, first_outputslab + FINISH_WAIT_RADIUS);
    INSTANTIATE(            FinishMultipoles, first_outputslab + FINISH_WAIT_RADIUS);
#ifdef PARALLEL
    INSTANTIATE(             TaylorTranspose, FORCE_RADIUS);  // 2D
    INSTANTIATE(                NeighborSend, first_outputslab + FINISH_WAIT_RADIUS);  // 2D
    INSTANTIATE(          CheckForMultipoles, first_outputslab + FINISH_WAIT_RADIUS);
#else
    INSTANTIATE_NOOP(        TaylorTranspose, FORCE_RADIUS);
    INSTANTIATE_NOOP(           NeighborSend, first_outputslab + FINISH_WAIT_RADIUS);
    INSTANTIATE_NOOP(     CheckForMultipoles, first_outputslab + FINISH_WAIT_RADIUS);
#endif

    // If group finding is disabled, we can make the dependencies no-ops so they don't hold up the pipeline
    #ifdef ONE_SIDED_GROUP_FINDING
        int first_groupslab = FORCE_RADIUS+1;
    #else
        int first_groupslab = FORCE_RADIUS+2*GROUP_RADIUS;
    #endif
    if(GFC != NULL){
        INSTANTIATE(          MakeCellGroups, FORCE_RADIUS);
        INSTANTIATE(      FindCellGroupLinks, FORCE_RADIUS + 1);
        INSTANTIATE(          DoGlobalGroups, first_groupslab);
        INSTANTIATE(               Microstep, first_groupslab);
        INSTANTIATE(            FinishGroups, first_groupslab);
    } else {
        INSTANTIATE_NOOP(     MakeCellGroups, FORCE_RADIUS);
        INSTANTIATE_NOOP( FindCellGroupLinks, FORCE_RADIUS + 1);
        INSTANTIATE_NOOP(     DoGlobalGroups, first_groupslab);
        INSTANTIATE_NOOP(          Microstep, first_groupslab);
        INSTANTIATE_NOOP(       FinishGroups, first_groupslab);
    }

    if(WriteState.Do2LPTVelocityRereading)
        INSTANTIATE(       UnpackLPTVelocity, 0);
    else
        INSTANTIATE_NOOP(  UnpackLPTVelocity, 0);


    // Let FetchSlabs start early, in case we want to overlap convolve and IO
    //while(FetchSlabs.Attempt()){}

    #ifdef PARALLEL
    TimeStepWallClock.Stop();
    ConvolutionWallClock.Clear(); ConvolutionWallClock.Start();

    ParallelConvolveDriver->Convolve();
    ParallelConvolveDriver->SendTaylors(FORCE_RADIUS);

    ConvolutionWallClock.Stop();
    ParallelConvolveDriver->CS.ConvolveWallClock = ConvolutionWallClock.Elapsed();
    MultipoleTransferCheck.Clear(); TaylorTransferCheck.Clear();

    TimeStepWallClock.Start();

    for (int slab = first + FORCE_RADIUS; slab < first + FORCE_RADIUS + total_slabs_on_node; slab ++ ){
        SB->AllocateArena(TaylorSlab, slab, RAMDISK_NO);
        ParallelConvolveDriver->RecvTaylorSlab(slab);
        STDLOG(2, "Set up to receive Taylor slab %d via MPI\n", slab);
    }

    // Lightweight setup of z-dimension exchanges
    SetupNeighborExchange(first + first_outputslab + FINISH_WAIT_RADIUS, total_slabs_on_node);
    #endif

	
	int timestep_loop_complete = 0; 
	while (!timestep_loop_complete){

        for(int i =0; i < FETCHPERSTEP; i++) FetchSlabs.Attempt();
         TransposePos.Attempt();
            NearForce.Attempt();

       ReceiveManifest->Check();  // This checks if Send is ready; no-op in non-blocking mode

          TaylorForce.Attempt();
                 Kick.Attempt();

       AttemptReceiveManifest();
       AttemptNeighborReceive(0,P.cpd);  // 2D
          MF->CheckAnyMPIDone();  // 2D

       MakeCellGroups.Attempt();
   FindCellGroupLinks.Attempt();
       DoGlobalGroups.Attempt();

       ReceiveManifest->Check();  // This checks if Send is ready; no-op in non-blocking mode

               Output.Attempt();
            Microstep.Attempt();
         FinishGroups.Attempt();

       AttemptReceiveManifest();
       AttemptNeighborReceive(0,P.cpd);  // 2D
          MF->CheckAnyMPIDone();  // 2D

    UnpackLPTVelocity.Attempt();
                Drift.Attempt();
         NeighborSend.Attempt();

       ReceiveManifest->Check();  // This checks if Send is ready; no-op in non-blocking mode

      FinishParticles.Attempt();
     FinishMultipoles.Attempt();

            CheckSendManifest();  // We look at each Send Manifest to see if there's material to free.
                        //   SendManifest->FreeAfterSend();

       AttemptReceiveManifest();
       AttemptNeighborReceive(0,P.cpd);  // 2D
          MF->CheckAnyMPIDone();  // 2D

   CheckForMultipoles.Attempt();

#ifdef PARALLEL
		timestep_loop_complete = CheckForMultipoles.alldone(total_slabs_on_node);
#else
		timestep_loop_complete = FinishMultipoles.alldone(total_slabs_on_node);
#endif
    }

    if(IL->length!=0)
        IL->DumpParticles();

    assertf(IL->length==0,
        "Insert List not empty (%d) at the end of timestep().  Time step too big?\n", IL->length);

    STDLOG(1,"Finished timestep dependency loop!\n");

    if (GFC != NULL) assertf(GFC->GLL->length==0,
	"GroupLinkList not empty (%d) at the end of timestep.  Global group finding didn't run properly.\n", GFC->GLL->length);

    uint64 total_n_output = n_output;
    if(GFC != NULL)
        total_n_output += GFC->n_L0_output;
	
	
    #ifdef PARALLEL	
        BarrierWallClock.Clear(); BarrierWallClock.Start();
        // These reductions force some synchronization, at least!
    	MPI_REDUCE_TO_ZERO(&total_n_output,   1, MPI_UINT64_T, MPI_SUM);		
        MPI_REDUCE_TO_ZERO(&merged_particles, 1, MPI_UINT64_T, MPI_SUM);		
		
        STDLOG(2,"Ready to proceed to the remaining work\n");

        MPI_Barrier(comm_global);
		BarrierWallClock.Stop();

        // This MPI call also forces a synchronization over the MPI processes,
        // so things like Reseting GPUs could fire multiple times on one node.
        SendManifest->FreeAfterSend();
        // Run this again, just in case the dependency loop on this node finished
        // before the neighbor received the non-blocking MPI transfer.

	    TimeStepWallClock.Stop(); ConvolutionWallClock.Start(); 
        convtimebuffer = (char*) malloc(CONVTIMEBUFSIZE);   // Need to allocate space for the timings
    	delete ParallelConvolveDriver;
	    ConvolutionWallClock.Stop(); TimeStepWallClock.Start();

    #endif
    if (MPI_rank==0)
        assertf(merged_particles == P.np, "Merged slabs contain %d particles instead of %d!\n", merged_particles, P.np);


    if(ReadState.DoTimeSliceOutput && MPI_rank==0){
        assertf(total_n_output == P.np, "TimeSlice output contains %d particles instead of %d!\n", total_n_output, P.np);
	}



    STDLOG(1,"Completing timestep()\n");
    TimeStepWallClock.Stop();
}
