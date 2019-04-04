/* timestep_recover_multipoles.cpp

This is a minimal slab pipeline that will load position slabs
and compute multipoles and write them out.  This allows us to
recover from a state where multipoles were lost, e.g. corrupted
or the disk died.  Also, when we backup a state, we usually don't
bother to back up the multipoles because this pipeline is quick
to run.

This pipeline will be automatically invoked when missing multipoles
are detected by the abacus.py wrapper.  make_multipoles.cpp is the
top-level entry point that will invoke this pipeline.

This file is directly #include'd in timestep.cpp since it borrows a 
lot of infrastructure from there.

*/

int FetchPosSlabPrecondition(int slab) {
    // read 10 slabs ahead
    if(slab > Finish.last_slab_executed + 10){
        return 0;
    }
    return 1;
}

void FetchPosSlabAction(int slab) {
    STDLOG(0,"Fetching slab %d with %d particles\n", slab, SS->size(slab));

    // Load directly into the merge slabs
    SB->LoadArenaNonBlocking(MergeCellInfoSlab, slab);
    SB->LoadArenaNonBlocking(MergePosSlab, slab);

    assertf(SS->size(slab)*sizeof(posstruct) <= fsize(SB->ReadSlabPath(PosSlab,slab).c_str()),
        "PosSlab size doesn't match prediction\n");
}

int FinishMultipolesPrecondition(int slab) {
    if( !SB->IsIOCompleted( MergePosSlab,      slab )
        || !SB->IsIOCompleted( MergeCellInfoSlab, slab )
        ) return 0;
    return 1;
}

void FinishMultipolesAction(int slab) {
    STDLOG(1,"Finishing multipole slab %d\n", slab);
        
    // Make the multipoles
    SB->AllocateArena(MultipoleSlab,slab);
    ComputeMultipoleSlab(slab);
    
    WriteMultipoleSlab.Start();
    SB->StoreArenaNonBlocking(MultipoleSlab,slab);
    WriteMultipoleSlab.Stop();
    
    SB->DeAllocate(MergePosSlab,slab);
    SB->DeAllocate(MergeCellInfoSlab,slab);
}


void timestepMultipoles(void) {
    STDLOG(0,"Initiating timestepMultipoles()\n");
    TimeStepWallClock.Clear();
    TimeStepWallClock.Start();
    
    FORCE_RADIUS = 0;  // so we know when we can free CellInfo in Finish
    GROUP_RADIUS = 0;

    int cpd = P.cpd; int first = first_slab_on_node;
    FetchSlabs.instantiate(cpd, first, &FetchPosSlabPrecondition, &FetchPosSlabAction );
    Finish.instantiate(cpd, first,  &FinishMultipolesPrecondition,  &FinishMultipolesAction );

    while( !Finish.alldone(total_slabs_on_node) ) {
        FetchSlabs.Attempt();
            Finish.Attempt();
           SendManifest->FreeAfterSend();
        ReceiveManifest->Check();   // This checks if Send is ready; no-op in non-blocking mode
        // If the manifest has been received, install it.
        if (ReceiveManifest->is_ready()) ReceiveManifest->ImportData();
    }

    STDLOG(1,"Completing timestepMultipoles()\n");
    TimeStepWallClock.Stop();
}
