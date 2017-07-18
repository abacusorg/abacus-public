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

#define FETCHAHEAD (3*FORCE_RADIUS +3)
#define FETCHPERSTEP 1
// Recall that all of these Dependencies have a built-in STimer
// to measure the amount of time spent on Actions.
Dependency FetchSlabs;
Dependency TransposePos;
Dependency NearForce;
Dependency TaylorForce;
Dependency Kick;
Dependency Drift;
Dependency Finish;
Dependency Output;
Dependency Group;
Dependency LPTVelocityReRead;


// The wall-clock time minus all of the above Timers might be a measure
// of the spin-locked time in the timestep() loop.
STimer TimeStepWallClock;

// -----------------------------------------------------------------

int FetchSlabPrecondition(int slab) {
    if(FetchSlabs.number_of_slabs_executed > 
            Kick.number_of_slabs_executed + FETCHAHEAD ) {
        // This was +1, but for non-blocking reads 
        // I think we want to work one more ahead
        return 0;
    }

    if(LBW->total_allocation > .5*P.MAXRAMMB*1024LLU*1024LLU){
        Dependency::NotifySpinning(NOT_ENOUGH_RAM);
        STDLOG(0,"Warning: unable to load more slabs due to RAM limits. Currently using %.2f GB, and the limit is 0.5*MAXRAMMB = %.2f GB.\n",
                LBW->total_allocation/1024./1024./1024., .5*P.MAXRAMMB/1024);
        return 0;
    }
    
    return 1;
}

void FetchSlabAction(int slab) {
    STDLOG(0,"Fetching slab %d with %d particles\n", slab, Slab->size(slab));
    // Load all of the particle files together
    LBW->LoadArenaNonBlocking(CellInfoSlab,slab);
    LBW->LoadArenaNonBlocking(PosSlab,slab);
    assertf(Slab->size(slab)*sizeof(posstruct)<=
        fsize(LBW->ReadSlabDescriptorName(PosSlab,slab).c_str()),
        "PosSlab size doesn't match prediction\n");
    LBW->LoadArenaNonBlocking(VelSlab,slab);
    assertf(Slab->size(slab)*sizeof(velstruct)<=
        fsize(LBW->ReadSlabDescriptorName(VelSlab,slab).c_str()),
        "VelSlab size doesn't match prediction\n");
    LBW->LoadArenaNonBlocking(AuxSlab,slab);
    assertf(Slab->size(slab)*sizeof(auxstruct)<=
        fsize(LBW->ReadSlabDescriptorName(AuxSlab,slab).c_str()),
        "AuxSlab size doesn't match prediction\n");
    LBW->LoadArenaNonBlocking(TaylorSlab,slab+FORCE_RADIUS);
}

// -----------------------------------------------------------------

int TransposePosPrecondition(int slab){
    if(   !LBW->IOCompleted(PosSlab, slab)
       || !LBW->IOCompleted(CellInfoSlab, slab)) {
        Dependency::NotifySpinning(WAITING_FOR_IO);
        return 0;
    }
        return 1;
}

void TransposePosAction(int slab){
#ifdef CUDADIRECT
    STDLOG(0,"Transposing position slab %d with %d particles\n", slab, Slab->size(slab));
    
    LBW->AllocateArena(PosXYZSlab, slab);
    int cpd = P.cpd;
    
    // Could do this over skewers; should make a skewersize(slab, y) function somewhere
    #pragma omp parallel for schedule(static)
    for(int y = 0; y < cpd; y++){
        for(int z = 0; z < cpd; z++){
            posstruct *pos = PP->PosCell(slab, y, z);
            List3<FLOAT> posxyz = PP->PosXYZCell(slab, y, z);
            int count = PP->CellInfo(slab,y,z)->count;
            
            #pragma ivdep
            for(int i = 0; i < count; i++){
                posxyz.X[i] = pos[i].x;
                posxyz.Y[i] = pos[i].y;
                posxyz.Z[i] = pos[i].z;
            }
        }
    }
#endif
}


// -----------------------------------------------------------------

int NearForcePrecondition(int slab) {
    for(int i=-FORCE_RADIUS;i<=FORCE_RADIUS;i++){
        if(TransposePos.notdone(slab+i))
            return 0;
        if( !LBW->IOCompleted( CellInfoSlab, slab+i ) ||
            !LBW->IOCompleted( PosSlab,      slab+i ) ||
            !LBW->IOCompleted( VelSlab,      slab+i ) ||  // are vel and aux actually near force preconditions?
            !LBW->IOCompleted( AuxSlab,      slab+i ) ) {
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

    SlabForceLatency[slab].Start();
    if (P.ForceOutputDebug) {
        // We want to output the NearAccSlab to the NearAcc file.
        // This must be a blocking write.
        JJ->Finalize(slab);
        LBW->WriteArena(NearAccSlab, slab, IO_KEEP, IO_BLOCKING,
        LBW->WriteSlabDescriptorName(NearAccSlab,slab).c_str());
    }
    
    //
}

// -----------------------------------------------------------------

int TaylorForcePrecondition(int slab) {
    if( NearForce.notdone(slab) ) return 0;  // Is this actually a Taylor precondition?
    if( !LBW->IOCompleted(TaylorSlab,slab) ){
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
    LBW->AllocateArena(AccSlab,slab);
    //ZeroAcceleration(slab,AccSlab);  // do we need this?  might be slow
    
    TaylorCompute.Start();
    ComputeTaylorForce(slab);
    TaylorCompute.Stop();

    if(P.ForceOutputDebug){
        // We want to output the AccSlab to the FarAcc file.
        // This must be a blocking write.
        LBW->WriteArena(AccSlab, slab, IO_KEEP, IO_BLOCKING,
                LBW->WriteSlabDescriptorName(FarAccSlab,slab).c_str());
    }
    LBW->DeAllocate(TaylorSlab,slab);
    SlabFarForceTime[slab].Stop();
}


// -----------------------------------------------------------------

int KickPrecondition(int slab) {
    // must have far forces
    if(TaylorForce.notdone(slab)){
        return 0;
    }
    
    //must have near forces
    if (NearForce.notdone(slab) || !JJ->SlabDone(slab)) {
#ifdef CUDADIRECT
        // Start the timer if we've gone one full loop without executing anything
        Dependency::NotifySpinning(WAITING_FOR_GPU);
#endif
        return 0;
    }

    return 1;
}

void KickAction(int slab) {
    SlabForceTime[slab].Stop();
    SlabForceLatency[slab].Stop();
    //If we are doing blocking forces, the finalization happens in NearForceAction
    if(!P.ForceOutputDebug && !P.ForceCPU)
        JJ->Finalize(slab);
    AddAccel.Start();
    RescaleAndCoAddAcceleration(slab);
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

int GroupPrecondition(int slab) {
    // We need to have kicked everything that matters
    for(int j=-GROUP_RADIUS;j<=GROUP_RADIUS;j++)  {
        if( Kick.notdone(slab+j) ) return 0;
    }
    return 1;
}

void GroupAction(int slab) {
    if (LPTStepNumber()) return;
    // We can't be doing group finding during an IC step

    STDLOG(1,"Finding groups for slab %d\n", slab);
    // No action yet, but this is where one would find groups and coevolution sets
    // One could also use the Accelerations to set individual particle microstep levels.
    // (one could imagine doing this in Force, but perhaps one wants the group context?)
}

// -----------------------------------------------------------------

int OutputPrecondition(int slab) {
    if (Kick.notdone(slab)) return 0;  // Must have kicked
    if (Group.notdone(slab)) return 0;  // Must have found groups
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
        STDLOG(1,"Outputting slab %d\n",slab);
        n_output += Output_TimeSlice(slab);
    }
    OutputTimeSlice.Stop();

    OutputLightCone.Start();
    if (ReadState.OutputIsAllowed) {
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

int FetchLPTVelPrecondition(int slab){
    // Don't read too far ahead
    if(LPTVelocityReRead.number_of_slabs_executed > 
            Drift.number_of_slabs_executed + 2*FINISH_WAIT_RADIUS + 1) {
        return 0;
    }
    
    // Allow .75 instead of .5 here because FetchSlabs might
    // eat up all the RAM.  The right way to do this is really to
    // make load_ic_vel_slab non-blocking through the IO module,
    // which first requires routing loadIC through the IO module.
    // Then, we can move load_ic_vel_slab to FetchSlabs.
    if(LBW->total_allocation > .75*P.MAXRAMMB*1024LLU*1024LLU){
        // Are we spinning because we need more RAM?
        Dependency::NotifySpinning(NOT_ENOUGH_RAM);
        return 0;
    }

    return 1;
}

void FetchLPTVelAction(int slab){
    // This is blocking because it uses the LoadIC module, not LBW
    load_ic_vel_slab(slab);
}

// -----------------------------------------------------------------

int DriftPrecondition(int slab) {
    // We cannot move these particles until they have been fully used as gravity sources
    for(int j=-FORCE_RADIUS;j<=FORCE_RADIUS;j++)  {
        if( TaylorForce.notdone(slab+j) ) return 0;
        if( NearForce.notdone(slab+j) ) return 0;
    }
    // We also must have Output this slab 
    if (Output.notdone(slab)) return 0;
    
    // We also must have the 2LPT velocities
    // The finish radius is a good guess of how ordered the ICs are
    if (WriteState.Do2LPTVelocityRereading)
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
    
    // We kept the accelerations until here because of third-order LPT
    if (P.StoreForces && !P.ForceOutputDebug) {
        STDLOG(1,"Storing Forces in slab %d\n", slab);
        LBW->StoreArenaBlocking(AccSlab,slab);
    }
    else{
        LBW->DeAllocate(AccSlab,slab);
    }
    LBW->DeAllocate(NearAccSlab,slab);
    LBW->DeAllocate(PosXYZSlab,slab);
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
    
    // Now delete the original particles
    LBW->DeAllocate(CellInfoSlab,slab);
    LBW->DeAllocate(PosSlab,slab);
    LBW->DeAllocate(VelSlab,slab);
    LBW->DeAllocate(AuxSlab,slab);

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
}


// ===================================================================


void timestep(void) {
    TimeStepWallClock.Clear();
    TimeStepWallClock.Start();
    STDLOG(1,"Initiating timestep()\n");

    FORCE_RADIUS = P.NearFieldRadius;
    GROUP_RADIUS = 0;
    assertf(FORCE_RADIUS != -1, "Illegal FORCE_RADIUS: %d\n", FORCE_RADIUS);
    assertf(GROUP_RADIUS != -1, "Illegal GROUP_RADIUS: %d\n", GROUP_RADIUS); 
    STDLOG(0,"Adopting FORCE_RADIUS = %d\n", FORCE_RADIUS);
    STDLOG(0,"Adopting GROUP_RADIUS = %d\n", GROUP_RADIUS);

    int cpd = P.cpd;
    int first_slab_to_load = 0;   // Might eventually be different
    int first = first_slab_to_load; 
    STDLOG(1,"First slab to load will be %d\n", first);

       FetchSlabs.instantiate(cpd, first, &FetchSlabPrecondition,     &FetchSlabAction     );
     TransposePos.instantiate(cpd, first, &TransposePosPrecondition,  &TransposePosAction  );
        NearForce.instantiate(cpd, first, &NearForcePrecondition,     &NearForceAction     );
      TaylorForce.instantiate(cpd, first, &TaylorForcePrecondition,   &TaylorForceAction   );
             Kick.instantiate(cpd, first, &KickPrecondition,          &KickAction          );
            Group.instantiate(cpd, first, &GroupPrecondition,         &GroupAction         );
           Output.instantiate(cpd, first, &OutputPrecondition,        &OutputAction        );
            Drift.instantiate(cpd, first, &DriftPrecondition,         &DriftAction         );
           Finish.instantiate(cpd, first, &FinishPrecondition,        &FinishAction        );
           
if(WriteState.Do2LPTVelocityRereading)
LPTVelocityReRead.instantiate(cpd, first + 2*FORCE_RADIUS + 1 - FINISH_WAIT_RADIUS,
                                          &FetchLPTVelPrecondition,   &FetchLPTVelAction   );

    while( !Finish.alldone() ) {
           for(int i =0; i < FETCHPERSTEP; i++) FetchSlabs.Attempt();
         TransposePos.Attempt();
            NearForce.Attempt();
          TaylorForce.Attempt();
                 Kick.Attempt();
                Group.Attempt();
               Output.Attempt();
  if(WriteState.Do2LPTVelocityRereading)
    LPTVelocityReRead.Attempt();
                Drift.Attempt();
               Finish.Attempt();
    }

    if(IL->length!=0)
        IL->DumpParticles();
    
    assertf(IL->length==0, 
        "Insert List not empty (%d) at the end of timestep().  Time step too big?\n", IL->length);
    
    assertf(merged_particles == P.np, "Merged slabs contain %"PRIu64" particles instead of %"PRIu64"!\n", merged_particles, P.np);
    
    if(ReadState.DoTimeSliceOutput)
        assertf(n_output == P.np, "TimeSlice output contains %"PRIu64" particles instead of %"PRIu64"!\n", n_output, P.np);

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

// When doing IC loading, we require that the neighboring slabs be loaded
// just to be sure that no particles have crossed the boundary.  This is trivial
// if we overload the Drift dependency with the FetchIC condition/actions.

void timestepIC(void) {
    STDLOG(0,"Initiating timestepIC()\n");
    TimeStepWallClock.Clear();
    TimeStepWallClock.Start();

    int cpd = P.cpd; int first = 0;
    Drift.instantiate(cpd, first, &FetchICPrecondition, &FetchICAction );
    Finish.instantiate(cpd, first,  &FinishPrecondition,  &FinishAction );

    while( !Finish.alldone() ) {
        Drift.Attempt();
       Finish.Attempt();
    }

    assertf(NP_from_IC == P.np, "Expected to read a total of %llu particles from IC files, but only read %llu.\n", P.np, NP_from_IC);
    assertf(merged_particles == P.np, "Merged slabs contain %"PRIu64" particles instead of %"PRIu64"!\n", merged_particles, P.np);
    
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
