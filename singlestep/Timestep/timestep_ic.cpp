/* timestep_ic.cpp

This file contains a minimal slab pipeline that is run during an IC step.
singlestep.cpp will call timestepIC() (instead of the usual timestep())
to invoke this pipeline.

This file is directly #include'd in timestep.cpp since it borrows a 
lot of infrastructure from there.

*/

uint64 NP_from_IC = 0;

int FetchICPrecondition(int slab) {
    // We always do this.
    #ifdef PARALLEL
    if (Drift.raw_number_executed>=total_slabs_on_node) return 0;
    // This prevents FetchSlabAction from reading beyond the 
    // range of slabs on the node.  In the PARALLEL code, these
    // data will arrive from the Manifest.  We have to implement
    // this on the raw_number because the reported number is adjusted
    // by the Manifest, which leads to a race condition when running
    // the PARALLEL code on a single node test.
    #endif
    return 1;
}
void FetchICAction(int slab) {
    STDLOG(1,"Fetching slab %d\n", slab);
    // Get a slab of particles and put them on the InsertList
    NP_from_IC += LoadSlab2IL(slab);
    
    // We also need to create a null slab
    // These slabs will never come off ramdisk because this is the first timestep
    SB->AllocateSpecificSize(PosSlab,slab, 0, RAMDISK_NO);
    SB->AllocateSpecificSize(VelSlab,slab, 0, RAMDISK_NO);
    SB->AllocateSpecificSize(AuxSlab,slab, 0, RAMDISK_NO);
    SB->AllocateArena(CellInfoSlab,slab, RAMDISK_NO);
    int cpd = CP->cpd;
    for (int y=0; y<cpd; y++)
        for (int z=0; z<cpd; z++) {
            CP->CellInfo(slab,y,z)->makenull();
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
    FINISH_WAIT_RADIUS = 2;  // The IC pipeline is very short; we have plenty of RAM to allow for large IC displacements

    int cpd = P.cpd; int first = first_slab_on_node;
    Drift.instantiate(cpd, first, &FetchICPrecondition, &FetchICAction );
    Finish.instantiate(cpd, first + FINISH_WAIT_RADIUS,  &FinishPrecondition,  &FinishAction );

    while( !Finish.alldone(total_slabs_on_node) ) {
        Drift.Attempt();
       Finish.Attempt();
       SendManifest->FreeAfterSend();
    ReceiveManifest->Check();   // This checks if Send is ready; no-op in non-blocking mode
    // If the manifest has been received, install it.
    if (ReceiveManifest->is_ready()) ReceiveManifest->ImportData();

    }

    STDLOG(1, "Read %d particles from IC files\n", NP_from_IC);
    #ifdef PARALLEL
        MPI_REDUCE_TO_ZERO(&merged_particles, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM);
        MPI_REDUCE_TO_ZERO(&NP_from_IC, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM);
        STDLOG(1,"Ready to proceed to the remaining work\n");
        
        MPI_Barrier(MPI_COMM_WORLD);
        // This MPI call also forces a syncrhonization over the MPI processes, 
        // so things like Reseting GPUs could fire multiple times on one node.
        SendManifest->FreeAfterSend();
        // Run this again, just in case the dependency loop on this node finished
       // before the neighbor received the non-blocking MPI transfer.
    #endif
    STDLOG(1, "Particles remaining on insert list: %d\n", IL->length);
    if (MPI_rank==0) {
        STDLOG(1, "Merged %d particles\n", merged_particles);
        assertf(merged_particles == P.np, "Merged slabs contain %d particles instead of %d!\n", merged_particles, P.np);
    }

    if(IL->length!=0)
        IL->DumpParticles();

    if(MPI_rank == 0)
        assertf(NP_from_IC == P.np, "Expected to read a total of %u particles from IC files, but only read %u.\n", P.np, NP_from_IC);

    assertf(IL->length==0,
        "Insert List not empty (%d) at the end of timestep().  Particles in IC files not sufficiently sorted?\n", IL->length);
    
    char filename[1024];
    sprintf(filename,"%s/slabsize",P.WriteStateDirectory);
    SS->write(filename);

    STDLOG(1,"Completing timestepIC()\n");
    TimeStepWallClock.Stop();
}
