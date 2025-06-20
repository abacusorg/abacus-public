// Copyright 2012-2025 The Abacus Developers
// SPDX-License-Identifier: GPL-3.0-or-later

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
    if(FetchSlabs.raw_number_executed >= total_slabs_on_node){
        return 0;
    }
    return 1;
}

void FetchPosSlabAction(int slab) {
    STDLOG(0,"Fetching slab {:d} with {:d} particles\n", slab, SS->size(slab));

    // Load directly into the merge slabs
    SB->LoadArenaNonBlocking(MergeCellInfoSlab, slab);
    SB->LoadArenaNonBlocking(MergePosSlab, slab);

    assertf(SS->size(slab)*sizeof(posstruct) <= fs::file_size(SB->ReadSlabPath(PosSlab,slab)),
        "PosSlab size doesn't match prediction\n");
}

int FinishMultipolesPrecondition(int slab) {
    if( !SB->IsIOCompleted( MergePosSlab,      slab )
        || !SB->IsIOCompleted( MergeCellInfoSlab, slab )
        ) return 0;
    return 1;
}

void FinishMultipolesAction(int slab) {
    STDLOG(1,"Finishing multipole slab {:d}\n", slab);
        
    // Make the multipoles
	int ramdisk_multipole_flag; 
	#ifdef PARALLEL
		ramdisk_multipole_flag = RAMDISK_NO;
	#else
		ramdisk_multipole_flag = RAMDISK_AUTO;
	#endif
		
    SB->AllocateArena(MultipoleSlab,slab, ramdisk_multipole_flag);
    ComputeMultipoleSlab(slab);
    
    WriteMultipoleSlab.Start();
#ifndef PARALLEL
    SB->StoreArenaNonBlocking(MultipoleSlab,slab);
#endif
    WriteMultipoleSlab.Stop();

#ifdef PARALLEL	
	QueueMultipoleMPI.Start();
	 STDLOG(2, "Attempting to SendMultipoleSlab {:d}\n", slab);
	 	ParallelConvolveDriver->SendMultipoleSlab(slab); //distribute z's to appropriate nodes for this node's x domain.
	if (Finish.raw_number_executed==0){ //if we are finishing the first slab, set up receive MPI calls for incoming multipoles.
		STDLOG(2, "Attempting to RecvMultipoleSlab {:d}\n", slab);
		ParallelConvolveDriver->RecvMultipoleSlab(slab); //receive z's from other nodes for all x's.
	}

	QueueMultipoleMPI.Stop();
#endif
	    
    SB->DeAllocate(MergePosSlab,slab);
    SB->DeAllocate(MergeCellInfoSlab,slab);
}


void timestepMultipoles(void) {	
    STDLOG(0,"Initiating multipole recovery timestep.()\n");
    TimeStepWallClock.Clear();
    TimeStepWallClock.Start();
    
    FORCE_RADIUS = 0;  // so we know when we can free CellInfo in Finish
    GROUP_RADIUS = 0;

#ifdef PARALLEL
    int create_MT_file = 1;  // that's... why we're here
    ParallelConvolveDriver = new ParallelConvolution(P.cpd, P.order, P.MultipoleDirectory, create_MT_file);
#endif

    int nslabs = P.cpd;
    int first = first_slab_on_node;

    FetchSlabs.instantiate(         nslabs,  first,  &FetchPosSlabPrecondition,           &FetchPosSlabAction,           "FetchPosSlab");
    Finish.instantiate(             nslabs,  first,  &FinishMultipolesPrecondition,       &FinishMultipolesAction,       "FinishMultipoles");
#ifdef PARALLEL
	CheckForMultipoles.instantiate( nslabs,  first,  &CheckForMultipolesPrecondition,  &CheckForMultipolesAction,  "CheckMultipoles"); 
#else
	CheckForMultipoles.instantiate( nslabs,  first,  &NoopPrecondition,  &NoopAction,  "CheckMultipoles"); 
#endif
	
	int timestep_loop_complete = 0; 
	while (!timestep_loop_complete){
        FetchSlabs.Attempt();
            Finish.Attempt();
           SendManifest->FreeAfterSend();
        ReceiveManifest->Check();   // This checks if Send is ready; no-op in non-blocking mode
        // If the manifest has been received, install it.
        if (ReceiveManifest->is_ready()) ReceiveManifest->ImportData();
		
		CheckForMultipoles.Attempt();	

	#ifdef PARALLEL
		timestep_loop_complete = CheckForMultipoles.alldone(total_slabs_on_node);
	#else
		timestep_loop_complete = Finish.alldone(total_slabs_on_node);
	#endif		
    }

    STDLOG(1,"Completing timestepMultipoles()\n");
    TimeStepWallClock.Stop();
}
