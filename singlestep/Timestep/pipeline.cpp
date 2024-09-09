/* pipeline.cpp
 *
 * pipeline.cpp defines the dependencies used in timestep.cpp. Each is a
 * SlabDependency object that defines an action() and precondition().
 *
 */

 // -----------------------------------------------------------------

class FetchSlabsDep : public SlabDependency {
public:

    FetchSlabsDep(int cpd, int initialslab)
        : SlabDependency("FetchSlabs", cpd, initialslab) { };

    int precondition(int slab){
        /*
        * The precondition for loading new slabs into memory
        * We limit the additional slabs read to FETCHAHEAD
        */
        // We want to read ahead enough that we are reading while computing forces
        // i.e. we would like to read a few slabs ahead of the Kick
        // But if the steps after Kick are slow enough that we aren't Finishing promptly,
        // we can get a buildup of slabs and run out of memory. So we tie to the Drift instead.
        if(wrap(slab - Drift->last_slab_executed) > FETCHAHEAD
            && done(Drift->last_slab_executed))
            return 0;

        #ifdef PARALLEL
        if (raw_number_executed>=total_slabs_on_node) return 0;
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
    void action(int slab) {
        STDLOG(1,"Fetching slab {:d} with {:d} particles ({:d} w/ ghost)\n",
            slab, SS->size(slab), SS->size_with_ghost(slab));

        // Load all of the particle files together
        SB->LoadArenaNonBlocking(CellInfoSlab,slab);
        SB->LoadArenaNonBlocking(PosSlab,slab);

        // Don't bother to load the vel/aux/taylors for slabs that won't be kicked until the wrap
        //#ifndef PARALLEL
        // LHG: only need this when we're memory starved
        #if 0
        if(FetchSlabs->number_of_slabs_executed < FORCE_RADIUS)
            return;
        #endif

    #ifdef PARALLEL
        // SB->AllocateArena(TaylorSlab, slab + FORCE_RADIUS, RAMDISK_NO);
    // 	ParallelConvolveDriver->RecvTaylorSlab(slab + FORCE_RADIUS);
    // 	STDLOG(2, "Received Taylor slab via MPI{:d}\n", slab + FORCE_RADIUS);

    #else
        SB->LoadArenaNonBlocking(TaylorSlab,slab);
    #endif

        SB->LoadArenaNonBlocking(VelSlab, slab);
        SB->LoadArenaNonBlocking(AuxSlab, slab);
    }
};


// -----------------------------------------------------------------


class TransposePosDep : public SlabDependency {
public:

    TransposePosDep(int cpd, int initialslab)
        : SlabDependency("TransposePos", cpd, initialslab){ }

    int precondition(int slab){
        // TODO: technically need separate WAITING_FOR_IO flags for each slab, otherwise "reasons" can cross-talk
        if(!SB->IsIOCompleted(PosSlab, slab)){
            if(SB->IsSlabPresent(PosSlab, slab))
                NotifySpinning(WAITING_FOR_IO);
            return 0;
        }

        if(!SB->IsIOCompleted(CellInfoSlab, slab)){
            if(SB->IsSlabPresent(CellInfoSlab, slab))
                NotifySpinning(WAITING_FOR_IO);
            return 0;
        }

        return 1;
    }

    void action(int slab){
        SB->AllocateArena(PosXYZSlab, slab);
        int cpd = P.cpd;

        NUMA_FOR(y,0,cpd, NO_CLAUSE, FALLBACK_DYNAMIC){
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
        NUMA_FOR_END;
    }
};


// -----------------------------------------------------------------
class NearForceDep : public SlabDependency {
public:

    NearForceDep(int cpd, int initialslab)
        : SlabDependency("NearForce", cpd, initialslab){ }

    int precondition(int slab) {
        for(int i=-FORCE_RADIUS;i<=FORCE_RADIUS;i++){
            // Technically, I think we only need the CellInfo to construct pencils
            // But it's convenient to have pos so the GPU can immediately execute any pencil
            if(TransposePos->notdone(slab+i))
                return 0;
            if( !SB->IsIOCompleted( CellInfoSlab, slab+i ) ){
                if(SB->IsSlabPresent(CellInfoSlab, slab+i))
                    SlabDependency::NotifySpinning(WAITING_FOR_IO);
                return 0;
            }
        }

        return 1;
    }

    void action(int slab) {
        // Do some data checks
        assertf(are_cellinfo_legal(slab, SS->size(slab), SS->size_with_ghost(slab)),
                "Cell info of slab {:d} contain out of bounds data\n", slab);
        // Could also check that the sum of the cell counts add up to SS->size(slab);

        SlabForceTime[slab].Start();

        int blocking = 0;
        NFD->ExecuteSlab(slab, blocking);

        SlabForceLatency[slab].Start();
        
        // Busy-wait for all GPU work for this slab to finish
        // while(!NFD->SlabDone(slab)) ;
    }
};

// -----------------------------------------------------------------

class TaylorTransposeDep : public SlabDependency {
public:
    TaylorTransposeDep(int cpd, int initialslab)
        : SlabDependency("TaylorTranspose", cpd, initialslab){ }

    int precondition(int slab){
    #ifdef PARALLEL
        if(!SB->IsSlabPresent(TaylorSlab, slab))
            return 0;

        int ready = ParallelConvolveDriver->CheckTaylorSlabReady(slab);
        if(!ready) {
            SlabDependency::NotifySpinning(WAITING_FOR_MPI);
            return 0;
        }
    #endif

        return 1;
    }

    void action(int slab){
        if(MPI_size_z > 1){
            // FFT and launch MPI All-to-all
            MTCOMPLEX *t = (MTCOMPLEX *) SB->GetSlabPtr(TaylorSlab, slab);
            TY->ComputeIFFTZAndMPI(slab, t);
        }

        // for the 1D parallel code, this is a no-op
    }
};

// -----------------------------------------------------------------

class TaylorForceDep : public SlabDependency {
public:

    TaylorForceDep(int cpd, int initialslab)
        : SlabDependency("TaylorForce", cpd, initialslab) { }

    int precondition(int slab) {
        if(TaylorTranspose->notdone(slab) ||  // 1D & 2D
            !TY->IsMPIDone(slab)){  // 2D
            return 0;
        }

        if( !SB->IsIOCompleted( CellInfoSlab, slab ) ){
            if(SB->IsSlabPresent(CellInfoSlab, slab))
                SlabDependency::NotifySpinning(WAITING_FOR_IO);
            return 0;
        }
        if( !SB->IsIOCompleted( PosSlab, slab ) ){
            if(SB->IsSlabPresent(PosSlab, slab))
                SlabDependency::NotifySpinning(WAITING_FOR_IO);
            return 0;
        }

    #ifndef PARALLEL
        if( !SB->IsIOCompleted( TaylorSlab, slab ) ){
            if(SB->IsSlabPresent(TaylorSlab, slab))
                SlabDependency::NotifySpinning(WAITING_FOR_IO);
            return 0;
        }
    #endif

        return 1;
    }

    void action(int slab) {
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
};

// -----------------------------------------------------------------

class KickDep : public SlabDependency {
    int set_aux_dens;

public:
    KickDep(int cpd, int initialslab)
        : SlabDependency("Kick", cpd, initialslab){

        set_aux_dens = ReadState.SetAuxDensity;

        if(set_aux_dens){
            STDLOG(0, "Will store densities in aux\n");
        }
    }

    int precondition(int slab) {
        if( !SB->IsIOCompleted( VelSlab, slab ) ){
            if(SB->IsSlabPresent(VelSlab, slab))
                SlabDependency::NotifySpinning(WAITING_FOR_IO);
            return 0;
        }

        // must have far forces
        if(TaylorForce->notdone(slab)){
            return 0;
        }

        //must have near forces
        if (NearForce->notdone(slab))
            return 0;
        else {
            // If pencil construction (NearForce) is finished, but not NFD, then we're waiting for the GPU result
            if(!NFD->SlabDone(slab)){
    #ifdef CUDADIRECT
                SlabDependency::NotifySpinning(WAITING_FOR_GPU);
    #endif
                return 0;
            }
        }

        if(set_aux_dens || LPTStepNumber() == 1){
            // Will write aux density for 2D group finding
            // or packing vel in aux for 2LPT
            if( !SB->IsIOCompleted( AuxSlab, slab ) ){
                if(SB->IsSlabPresent(AuxSlab, slab))
                    SlabDependency::NotifySpinning(WAITING_FOR_IO);
                return 0;
            }
        }

        if(LPTStepNumber() == 1){
            // Need linear velocities for 2LPT
            if( !SB->IsIOCompleted( PosSlab, slab ) ){
                if(SB->IsSlabPresent(PosSlab, slab))
                    SlabDependency::NotifySpinning(WAITING_FOR_IO);
                return 0;
            }
        }

        return 1;
    }

    void action(int slab) {
        SlabForceTime[slab].Stop();
        SlabForceLatency[slab].Stop();
        
        // computing CPU forces must be done before the pos releases below
        if(!P.ForceCPU && P.ForceOutputDebug)
            NFD->CheckGPUCPU(slab);

        // Release the trailing slab if it won't be needed at the wrap
        // Technically we could release it anyway and re-do the transpose from PosSlab,
        // but if we're not doing group finding we may have already written and released PosSlab
        KickDealloc.Start();
        if(Kick->raw_number_executed >= 2*FORCE_RADIUS)
            SB->DeAllocate(PosXYZSlab, slab - FORCE_RADIUS);

        // Special case: if this is the last slab, free all +/- FORCE_RADIUS
        // Not worrying about this in the PARALLEL case; we have other non-destructions
        STDLOG(1,"{:d} slabs have been Kicked so far\n", Kick->raw_number_executed);
        if(Kick->raw_number_executed == total_slabs_on_node-1)
            for(int j = slab - FORCE_RADIUS+1; j <= slab + FORCE_RADIUS; j++)
                SB->DeAllocate(PosXYZSlab, j);

        //#ifndef PARALLEL
        #if 0
        // Queue up slabs near the wrap to be loaded again later
        // This way, we don't have idle slabs taking up memory while waiting for the pipeline to wrap around
        if(Kick->number_of_slabs_executed < FORCE_RADIUS){
            STDLOG(3,"Marking slab {:d} for repeat\n", slab - FORCE_RADIUS);
            TransposePos->mark_to_repeat(slab - FORCE_RADIUS);
            SB->DeAllocate(PosXYZSlab, slab - FORCE_RADIUS);
            // The first two won't need PosSlab until the second time around
            //SB->DeAllocate(PosSlab, slab - FORCE_RADIUS);
            SB->DeAllocate(CellInfoSlab, slab - FORCE_RADIUS);
            FetchSlabs->mark_to_repeat(slab - FORCE_RADIUS);
        }
        #endif

        KickDealloc.Stop();
        
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
                SB->WriteSlabPath(NearAccSlab,slab));

    #ifdef DIRECTSINGLESPLINE
            // restore the original
            memcpy(nearacc, nearacctmp, npslab*sizeof(accstruct));
            delete[] nearacctmp;
    #endif
        }

        AddAccel.Start();
        RescaleAndCoAddAcceleration(slab);
        AddAccel.Stop();
        
        KickDealloc.Start();
        SB->DeAllocate(FarAccSlab,slab);
        KickDealloc.Stop();

        int step = LPTStepNumber();
        KickCellTimer.Start();
        if (step) {
            // We have LPT IC work to do
            if (step==1) {
                STDLOG(1,"Kicking slab {:d} as LPT step 1\n", slab);
                KickSlab(slab, 0, 0, 0, KickCell_2LPT_1);
            } else if (step==2) {
                STDLOG(1,"Kicking slab {:d} as LPT step 2\n", slab);
                KickSlab(slab, 0, 0, 0, KickCell_2LPT_2);
            } else if (step==3) {
                STDLOG(1,"Kicking slab {:d} as LPT step 3\n", slab);
                KickSlab(slab, 0, 0, 0, KickCell_2LPT_3);
            } else QUIT("LPT Kick {:d} not implemented\n", step);
        } else {
            // This is just a standard step
            FLOAT kickfactor1 =  ReadState.LastHalfEtaKick;
            FLOAT kickfactor2 =  WriteState.FirstHalfEtaKick;
            STDLOG(1,"Kicking slab {:d} by {:f} + {:f}\n", slab, kickfactor1, kickfactor2);
            KickSlab(slab, kickfactor1, kickfactor2, set_aux_dens, KickCell);
        }
        KickCellTimer.Stop();

        ReleaseFreeMemoryToKernel();
    }
};

// -----------------------------------------------------------------

class MakeCellGroupsDep : public SlabDependency {
public:
    MakeCellGroupsDep(int cpd, int initialslab)
        : SlabDependency("MakeCellGroups", cpd, initialslab){ }

    int precondition(int slab) {
        // Only PosXYZ is used for sources, so we're free to rearrange PosSlab
        // in group finding after the transpose
        if( TransposePos->notdone(slab) ) return 0;

        // Need the new velocities because we're going to rearrange particles
        if( Kick->notdone(slab) ) return 0;

        // Also need the auxs, because we're going to re-order
        if( !SB->IsIOCompleted( AuxSlab, slab ) ){
            if(SB->IsSlabPresent(AuxSlab, slab))
                SlabDependency::NotifySpinning(WAITING_FOR_IO);
            return 0;
        }
        return 1;
    }

    void action(int slab) {
        GFC->ConstructCellGroups(slab);
    }
};

// -----------------------------------------------------------------

class FindCellGroupLinksDep : public SlabDependency {
public:
    FindCellGroupLinksDep(int cpd, int initialslab)
        : SlabDependency("FindCellGroupLinks", cpd, initialslab){ }

    int precondition(int slab) {
        // We want to find all links between this slab and the one just behind
        for (int j=-1; j<=0; j++)
            if (MakeCellGroups->notdone(slab+j)) return 0;
        return 1;
    }

    void action(int slab) {
        // Find links between slab and slab-1
        FindGroupLinks(slab);
    }
};

// -----------------------------------------------------------------

class DoGlobalGroupsDep : public SlabDependency {
public:
    DoGlobalGroupsDep(int cpd, int initialslab)
        : SlabDependency("DoGlobalGroups", cpd, initialslab){ }

    int precondition(int slab) {
    #ifdef ONE_SIDED_GROUP_FINDING
        /* We're going to search for GlobalGroups that include cells from
        slabs [slab,slab+2*GroupRadius].  However, this also includes the idea
        that we will *not* find groups if they include anything in slab-1.
        So we must have GroupLinks between [slab-1,slab] as well.  Moreover,
        we must go all the way up to [slab+2*GR-1,slab+2*GR], as there may
        be groups that are only now becoming eligible.
        */

        if (Kick->notdone(slab)) return 0;
        for (int j=0; j<=2*GROUP_RADIUS; j++){
            if (FindCellGroupLinks->notdone(slab+j)) return 0;
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
            if (DoGlobalGroups->done(slab+j-1)) break;
            if (FindCellGroupLinks->notdone(slab+j)) return 0;
        }
        // Look ahead; can stop as soon as one finds a closed slab
        for (int j=1; j<=2*GROUP_RADIUS; j++){
            if (DoGlobalGroups->done(slab+j)) break;
            if (FindCellGroupLinks->notdone(slab+j)) return 0;
        }
    #endif
        return 1;
    }

    void action(int slab) {
        FindAndProcessGlobalGroups(slab);

        // The first 2*GroupRadius times we get here, we can attempt to free
        // info from slab.  The Manifest code sends everything <S, so we need S=slab+1
        #ifdef ONE_SIDED_GROUP_FINDING
            if (DoGlobalGroups->raw_number_executed<2*GROUP_RADIUS) {
                SendManifest->QueueToSend(slab+1);
                SendManifest++;
            }
        #endif
        STDLOG(0,"Exiting Find Global Groups action in slab {:d}\n", slab);

    }
};

// -----------------------------------------------------------------

class MicrostepDep : public SlabDependency {
public:
    MicrostepDep(int cpd, int initialslab)
        : SlabDependency("Microstep", cpd, initialslab){ }

    int precondition(int slab){
        // We are going to second-half kick this slab
        if (DoGlobalGroups->notdone(slab))
            return 0;
        return 1;
    }

    void action(int slab){

        // TODO: This is now not the place to do this.
        // All kicks (and half-unkicks) for output are done; discard accels.
        // We de-allocate in Drift if we aren't doing group finding
        // SB->DeAllocate(AccSlab,slab);

        return;
        MicrostepCPU.Start();
        // Do microstepping here
        if(MicrostepEpochs != NULL){
            STDLOG(1,"Beginning microsteps for slab {:d}\n", slab);
            MicrostepControl *MC = new MicrostepControl;
            MC->setup(GFC->globalslabs[slab], *MicrostepEpochs, P.MicrostepTimeStep, NFD->eps);
            //MC->LaunchGroupsGPU();
            MC->ComputeGroupsCPU();

            GFC->microstepcontrol[slab] = MC;
        }

        //SlabDependency do_action() assumes that each SlabDependency processes all particles in a given slab.
        //Microstepping is an exception; it only does the group particles! Correct the bookkeeping here.
        Microstep->num_particles += GFC->globalslabs[slab]->np - SS->size(slab);

        MicrostepCPU.Stop();
    }
};

// -----------------------------------------------------------------

class FinishGroupsDep : public SlabDependency {
public:
    FinishGroupsDep(int cpd, int initialslab)
        : SlabDependency("FinishGroups", cpd, initialslab){ }

    int precondition(int slab){
        // Is the asychronous GPU microstepping done?
        //if (!GFC->microstepcontrol[slab]->GPUGroupsDone()) return 0

        // We are going to release these groups.
        // TODO: If Microstep changes, this may change.  At present, we really need Output to be done
        // because this step triggers the writeback of the new Pos/Vel.
        if (Microstep->notdone(slab)) return 0;

        return 1;
    }

    void action(int slab){
        // Scatter pos,vel updates to slabs, and release GGS
        FinishGlobalGroups(slab);   // This will Scatter Pos/Vel
        delete GFC->microstepcontrol[slab];
        GFC->microstepcontrol[slab] = NULL;

        // TODO: When we're ready to send Group-based Manifests, it would go here.
        // Would pass slab+1 to the manifest code as the faux finished slab.
    }
};

// -----------------------------------------------------------------
/*
 * Checks if we are ready to do all outputs for this step
 * Anything that modifies the particles at the current time should happen before here
 * Importantly: the OutputAction is only for non-L0 particles.
 * L0 particle outputs need to happen in the DoGlobalGroupsAction(),
 * before Microstepping.
 */
class OutputDep : public SlabDependency {
public:
    uint64 n_output = 0;
    uint64 np_lightcone = 0;

    OutputDep(int cpd, int initialslab)
        : SlabDependency("Output", cpd, initialslab){ }

    int precondition(int slab) {

        #ifdef ONE_SIDED_GROUP_FINDING
        /* This must wait until all groups including the slab have been found.
        It used to be that closing groups on the current slab S would do this,
        but now groups are only found looking upwards from S, so we need to have
        closed groups on all slabs from S to S-2*GroupRadius, inclusive.
        */
        for (int s=0; s<=2*GROUP_RADIUS; s++)
            if (FinishGroups->notdone(slab-s)) return 0;
        // Must have found groups to be able to output light cones
        // note that group outputs were already done
        #else
            if (FinishGroups->notdone(slab)) return 0;
        #endif

        // Note the following conditions only have any effect if group finding is turned off

        if (Kick->notdone(slab)) return 0;  // Must have accelerations

        // Also obviously need the aux!  This was checked in CellGroups,
        // but that may have been skipped if there's no group finding.
        if( !SB->IsIOCompleted( AuxSlab, slab ) ){
            if(SB->IsSlabPresent(AuxSlab, slab))
                SlabDependency::NotifySpinning(WAITING_FOR_IO);
            return 0;
        }

        return 1;
    }

    void action(int slab) {

        // We are finally done with the groups for this slab.
        // Delete the Cell Groups.
        if (GFC!=NULL) GFC->DestroyCellGroups(slab);

        if (LPTStepNumber()>0) return;
        // Some output might want to be skipped during an IC step,
        // e.g., no light cones

        OutputTimeSlice.Start();


        // Having found all groups, we should output the Non-L0 (i.e., field) Taggable subsample.
        if(ReadState.DoSubsampleOutput) {
            assertf(ReadState.DoGroupFindingOutput, "Subsample output should turn on group finding!\n"); // Currently Subsample Output requires GroupFinding Output
            OutputNonL0Taggable(slab);
        }

        if (ReadState.DoTimeSliceOutput) {
            // If we are doing group finding, then we are doing group finding output and subsample output
            assertf(GFC == NULL || (ReadState.DoSubsampleOutput && ReadState.DoGroupFindingOutput),
                "Preparing for timeslice output, expected either no group finding, or group finding and subsampling output!\n");

            // We've already done a K(1) and thus need a K(-1/2)
            FLOAT unkickfactor = WriteState.FirstHalfEtaKick;
            STDLOG(1,"Outputting slab {:d} with unkick factor {:f}\n",slab, unkickfactor);
            n_output += Output_TimeSlice(slab, unkickfactor);
        }
        OutputTimeSlice.Stop();

        OutputLightCone.Start();
        if (ReadState.OutputIsAllowed) {
            // TODO: LightCones may need a half un-kick if GFC == NULL
            // but we can probably handle that in the interpolation
            for(size_t i = 0; i < LCOrigin.size(); i++){
                // STDLOG(2,"Outputting LightCone {:d} (origin ({:f},{:f},{:f})) for slab {:d}\n",i,LCOrigin[i].x,LCOrigin[i].y,LCOrigin[i].z,slab);
                // Start timing
                // STimer lightConeTimer;
                // lightConeTimer.Start();
                this->np_lightcone += makeLightCone(slab,i);
                // lightConeTimer.Stop();
                // STDLOG(2, "LightCone {:d} for slab {:d} creation took {:f} seconds\n", i, slab, lightConeTimer.Elapsed());
            }
        }
        OutputLightCone.Stop();

        OutputBin.Start();
        if(ReadState.DoBinning){
            int ystride = CP->cpd /omp_get_max_threads();
            int minstride = 12;
            if (ystride < minstride) ystride = minstride;
            int cpd = CP->cpd;
            STDLOG(1,"Binning particles for slab {:d}\n",slab);
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
};

// -----------------------------------------------------------------
/* Checks if we are ready to apply a full-timestep drift to the particles
 * Any operation that relies on the positions and velocities being synced
 * should be checked for completion this year.
 */
class DriftDep : public SlabDependency {
public:
    static int drift_ahead;  // how far ahead of Finish to work

    DriftDep(int cpd, int initialslab)
        : SlabDependency("Drift", cpd, initialslab){ }

    int precondition(int slab) {
        // We must have finished scattering into this slab
        // if (FinishGroups->notdone(slab)) return 0;

        // We will move the particles, so we must have done outputs
        if (Output->notdone(slab)) return 0;

        // We can't move particles until they've been used as gravity sources
        // However, we only use PosXYZSlab as sources, so we're free to move PosSlab

        if(this->raw_number_executed > FinishParticles->raw_number_executed + drift_ahead){
            // Don't flood the IL
            return 0;
        }

        return 1;
    }

    void action(int slab) {
        int step = LPTStepNumber();
        if (step) {
            // We have LPT IC work to do
            if (step==1) {
                STDLOG(1,"Drifting slab {:d} as LPT step 1\n", slab);
                DriftAndCopy2InsertList(slab, 0, DriftCell_2LPT_1);
            } else if (step==2) {
                STDLOG(1,"Drifting slab {:d} as LPT step 2\n", slab);
                DriftAndCopy2InsertList(slab, 0, DriftCell_2LPT_2);
            } else if (step==3) {
                STDLOG(1,"Drifting slab {:d} as LPT step 3\n", slab);
                DriftAndCopy2InsertList(slab, 0, DriftCell_2LPT_3);
            } else QUIT("LPT Drift {:d} not implemented\n", step);
        } else {
            //         This is just a normal drift
            FLOAT driftfactor = WriteState.DeltaEtaDrift;
            STDLOG(1,"Drifting slab {:d} by {:f}\n", slab, driftfactor);
            //DriftAndCopy2InsertList(slab, driftfactor, DriftCell);
            DriftPencilsAndCopy2InsertList(slab, driftfactor, DriftPencil);
        }

        if(NFD){
            // We kept the accelerations until here because of third-order LPT
            if (P.StoreForces == 1 || (P.StoreForces == 3 && ReadState.DoTimeSliceOutput)) {
                STDLOG(1,"Storing Forces in slab {:d}\n", slab);
                if(this->raw_number_executed == 0){
                    fs::path dir = P.OutputDirectory / "acc" / fmt::format("Step{:04d}", ReadState.FullStepNumber);
                    fs::create_directories(dir);

                }
                SB->StoreArenaBlocking(AccSlab,slab);
            }
            else{
                SB->DeAllocate(AccSlab,slab);
                ReleaseFreeMemoryToKernel();
            }
        }
    }
};

int DriftDep::drift_ahead = 2*FINISH_WAIT_RADIUS + 1 + 3;  // +3 for slosh

// -----------------------------------------------------------------

class NeighborSendDep : public SlabDependency {
public:
    NeighborSendDep(int cpd, int initialslab)
        : SlabDependency("NeighborSend", cpd, initialslab){ }

    int precondition(int slab){
        // Must have all particles drifted into this slab
        for(int j=-FINISH_WAIT_RADIUS;j<=FINISH_WAIT_RADIUS;j++) {
            if( Drift->notdone(slab+j) )
                return 0;
        }

        if (NeighborSend->raw_number_executed == total_slabs_on_node) return 0;

        return 1;
    }

    void action(int slab){
        // Send the full list of neighbor particles for this slab
        DoNeighborSend(slab);  // in neighbor_exchange.cpp
    }
};

// -----------------------------------------------------------------

class NeighborRecvEvent : public EventDependency {
public:
    static int receive_ahead;
    NeighborRecvEvent()
        : EventDependency("NeighborReceive") { }

    int action(){
        return AttemptNeighborReceive(FinishParticles->last_slab_executed + 1, receive_ahead);
    }
};

int NeighborRecvEvent::receive_ahead = 3;  // slosh

// -----------------------------------------------------------------

class FinishParticlesDep : public SlabDependency {
public:
    uint64 merged_primaries = 0;

    FinishParticlesDep(int cpd, int initialslab)
        : SlabDependency("FinishParticles", cpd, initialslab){ }

    int precondition(int slab) {
        for(int j=-FINISH_WAIT_RADIUS;j<=FINISH_WAIT_RADIUS;j++) {
            if( Drift->notdone(slab+j) )
                return 0;
        }

        if (FinishParticles->alldone(total_slabs_on_node)) return 0;

        if( !IsNeighborReceiveDone(slab) ){
            // This is an effective dependency on NeighborSend, because we won't
            // receive before sending.
            // We only need to receive 1 slab, not FWR, because the remote node
            // waits for FWR before sending.
            return 0;
        }

        return 1;
    }

    void action(int slab) {
        // Usually actions are log level 1, but let's put this at level 0
        STDLOG(0, "Executing FinishParticles on slab {:d}\n", slab);
        
        FinishFreeSlabs.Start();

        SB->report_current();
        SB->report_peak();

        FinishFreeSlabs.Stop();

        // Gather particles from the insert list and make the merge slabs
        uint64 n_merge, n_merge_with_ghost;
        FillMergeSlab(slab, &n_merge, &n_merge_with_ghost);
        merged_primaries += n_merge;

        // manually adjust Mpart/s rates to reflect new size
        FinishParticles->num_particles += n_merge - SS->size(slab);
        FinishParticles->num_particles_with_ghost += n_merge_with_ghost - SS->size_with_ghost(slab);

        FinishFreeSlabs.Start();

        // This may be the last time be need any of the CellInfo slabs that we just used
        // We can't immediately free CellInfo before NearForce might need it until we're FORCE_RADIUS away
        // An alternative to this would be to just wait for FORCE_RADIUS before finishing
        for(int j = -2*FORCE_RADIUS, consec = 0; j <= 2*FORCE_RADIUS; j++){
            if(FinishParticles->done(slab + j) || j == 0)
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

        STDLOG(2,"Done deallocing pos, vel, aux for slab {:d}\n", slab);
        FinishFreeSlabs.Stop();

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
        if (FinishParticles->raw_number_executed==0) SendManifest->QueueToSend(slab);
        #endif

    }
};

// -----------------------------------------------------------------

class FinishMultipolesDep : public SlabDependency {
public:
    FinishMultipolesDep(int cpd, int initialslab)
        : SlabDependency("FinishMultipoles", cpd, initialslab){ }

    int precondition(int slab){
        return FinishParticles->done(slab) && MF->IsMPIDone(slab);
    }

    void action(int slab){
        // In the 2D code, the MPI transpose is now complete.
        // In the 1D code, this is just a continuation of FinishParticles

        if(MPI_size_z > 1){
            STDLOG(1, "Executing multipoles z-FFT for slab {:d}\n", slab);
            MTCOMPLEX *slabptr = (MTCOMPLEX *) SB->AllocateArena(MultipoleSlab, slab, ramdisk_multipole_flag);
            MF->ComputeFFTZ(slab, slabptr);
        }

    #ifdef PARALLEL
        QueueMultipoleMPI.Start();
        STDLOG(2, "Attempting to SendMultipoleSlab {:d}\n", slab);
        // distribute z's to appropriate nodes for this node's x domain.
        ParallelConvolveDriver->SendMultipoleSlab(slab);
        // if we are finishing the first slab, set up receive MPI calls for incoming multipoles.
        if (FinishMultipoles->raw_number_executed==0){
            STDLOG(2, "Attempting to RecvMultipoleSlab {:d}\n", slab);
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

        int pwidth = FetchSlabs->raw_number_executed - FinishMultipoles->raw_number_executed;
        STDLOG(1, "Current pipeline width (N_fetch - N_finish) is {:d}\n", pwidth);
        // release is cheap but not totally free, so might run every few Finishes
        //if (FinishMultipoles->raw_number_executed % 3 == 0)
            ReleaseFreeMemoryToKernel();
        ReportMemoryAllocatorStats();
    }
};

// -----------------------------------------------------------------

#ifdef PARALLEL
class CheckForMultipolesDep : public SlabDependency {
public:
    CheckForMultipolesDep(int cpd, int initialslab)
        : SlabDependency("CheckForMultipoles", cpd, initialslab){ }

    int precondition(int slab) {
        if( FinishMultipoles->notdone(slab) ) return 0;
        
        if ( ParallelConvolveDriver->CheckForMultipoleTransferComplete(slab) ) return 1;
        else {
            if(SB->IsSlabPresent(MultipoleSlab, slab)) SlabDependency::NotifySpinning(WAITING_FOR_MPI);
            return 0;
        }
    }

    void action(int slab) {
        SB->DeAllocate(MultipoleSlab, slab);
    }
};

#endif

// -----------------------------------------------------------------

class NoopDep : public SlabDependency {
public:
    NoopDep(int cpd)
        : SlabDependency("NOOP", cpd, 0){
        
        number_of_slabs_executed = cpd;
        raw_number_executed = cpd;
    }

    int precondition(int slab){ return 1; }
    void action(int slab){ return; }

    int Attempt(void) { return 0; }
    int done(int s) { return 1; }
};

class NoopEvent : public EventDependency {
public:
    NoopEvent() : EventDependency("NOOP"){ }

    int precondition(){ return 1; }
    int action(){ return 0; }

    int Attempt(){ return 0; }
};

// ===================================================================

class ReceiveManifestEvent : public EventDependency {
public:
    ReceiveManifestEvent()
        : EventDependency("ReceiveManifest") { }

    int action(){
        int didsomething = ReceiveManifest->Check();

        if(ReceiveManifest->is_ready()){
            // If the manifest has been received, install it.
            ReceiveManifest->ImportData();
            ReceiveManifest++;
            STDLOG(1, "Readying the next Manifest, number {:d}\n", ReceiveManifest-_ReceiveManifest);
            ReceiveManifest->SetupToReceive();
            didsomething = 1;
        }
        return didsomething;
    }
};

// -----------------------------------------------------------------

class SendManifestEvent : public EventDependency {
public:
    SendManifestEvent()
        : EventDependency("SendManifest") { }

    int action(){
        // CheckSendManifest()
        // We look at each Send Manifest to see if there's material to free.
        int ret = 0;
        for (int j=0; j<nManifest; j++) ret |= _SendManifest[j].FreeAfterSend();
        return ret;
    }
};

// -----------------------------------------------------------------

class Multipole2DMPIEvent : public EventDependency {
public:
    Multipole2DMPIEvent()
        : EventDependency("Check Multipoles 2D MPI") { }

    int action(){
        return MF->CheckAnyMPIDone();
    }
};

// -----------------------------------------------------------------

class Taylor2DMPIEvent : public EventDependency {
public:
    Taylor2DMPIEvent()
        : EventDependency("Check Taylors 2D MPI") { }

    int action(){
        return TY->CheckAnyMPIDone();
    }
};

// -----------------------------------------------------------------
