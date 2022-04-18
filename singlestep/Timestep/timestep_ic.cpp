/* timestep_ic.cpp

This file contains a minimal slab pipeline that is run during an IC step.
singlestep.cpp will call timestepIC() (instead of the usual timestep())
to invoke this pipeline.

This file is directly #include'd in timestep.cpp since it borrows a 
lot of infrastructure from there.

*/

Dependency ReadIC;

uint64 NP_from_IC = 0;

int ReadICPrecondition(int slab) {
    // We always do this.
    #ifdef PARALLEL
    if (ReadIC.raw_number_executed>=total_slabs_on_node) return 0;
    // This prevents ReadICAction from reading beyond the 
    // range of slabs on the node.  In the PARALLEL code, these
    // data will arrive from the Manifest.  We have to implement
    // this on the raw_number because the reported number is adjusted
    // by the Manifest, which leads to a race condition when running
    // the PARALLEL code on a single node test.
    #endif
    return 1;
}

void ReadICAction(int slab) {
    // Read an IC slab file into an arena
    // It will be unpacked into particles in a later dependency
    unique_ptr<ICFile> ic = ICFile::FromFormat(P.ICFormat, slab);
    ic->read_nonblocking();
}

int UnpackICPrecondition(int slab){
    if(ReadIC.notdone(slab))
        return 0;
    
    unique_ptr<ICFile> ic = ICFile::FromFormat(P.ICFormat, slab);
    if(!ic->check_read_done())
        return 0;

    return 1;
}

void UnpackICAction(int slab){
    uint64 NP_thisslab = UnpackICtoIL(slab);
    NP_from_IC += NP_thisslab;

    // Record the number of particles read in SlabSize so the timing log uses the right particle count
    SS->setold(slab, NP_thisslab, NP_thisslab);

    // We also need to create a null slab
    // These slabs will never come off ramdisk because this is the first timestep
    SB->AllocateSpecificSize(PosSlab,slab, 0, RAMDISK_NO);
    SB->AllocateSpecificSize(VelSlab,slab, 0, RAMDISK_NO);
    SB->AllocateSpecificSize(AuxSlab,slab, 0, RAMDISK_NO);
    SB->AllocateArena(CellInfoSlab,slab, RAMDISK_NO);
    
    int cpd = CP->cpd;

    #pragma omp parallel for schedule(static)
    for (int y=0; y<cpd; y++)
        for (int z=node_z_start_ghost; z<node_z_start_ghost + node_z_size_with_ghost; z++)
            CP->CellInfo(slab,y,z)->makenull();

    DoNeighborSend(slab);
}

/*
 * Registers the preconditions and actions for an IC step
 * When doing IC loading, we require that the neighboring slabs be loaded
 * just to be sure that no particles have crossed the boundary.  This is trivial
 * if we overload the Drift dependency with the FetchIC condition/actions.
 */

void timestepIC(void) {
#ifdef PARALLEL
	
    // The IC step doesn't convolve, but it does need to distribute multipoles
	int create_MT_file = 1; 
	ConvolutionWallClock.Clear(); ConvolutionWallClock.Start();
	ParallelConvolveDriver = new ParallelConvolution(P.cpd, P.order, P.MultipoleDirectory, create_MT_file);
	ConvolutionWallClock.Stop(); 
	ParallelConvolveDriver->CS.ConvolveWallClock = ConvolutionWallClock.Elapsed(); 

#endif
	
    STDLOG(0,"Initiating timestepIC()\n");
    STDLOG(0,"Adopting FINISH_WAIT_RADIUS = %d\n", FINISH_WAIT_RADIUS);

    TimeStepWallClock.Clear();
    TimeStepWallClock.Start();
    
    FORCE_RADIUS = 0;  // so we know when we can free CellInfo in Finish
    GROUP_RADIUS = 0;
    FINISH_WAIT_RADIUS = 2;  // The IC pipeline is very short; we have plenty of RAM to allow for large IC displacements

    int nslabs = P.cpd;
    int first = first_slab_on_node;

    // Lightweight setup of z-dimension exchanges
    SetupNeighborExchange(first, total_slabs_on_node);

    INSTANTIATE(ReadIC, 0);
    Drift.instantiate(nslabs, first, &UnpackICPrecondition, &UnpackICAction, "UnpackIC");
    INSTANTIATE(FinishParticles, FINISH_WAIT_RADIUS);
    INSTANTIATE(FinishMultipoles, FINISH_WAIT_RADIUS);

#ifdef PARALLEL
	INSTANTIATE(CheckForMultipoles, FINISH_WAIT_RADIUS);
#else
	INSTANTIATE_NOOP(CheckForMultipoles, FINISH_WAIT_RADIUS);
#endif
	
	
	int timestep_loop_complete = 0; 
	while (!timestep_loop_complete){
        ReadIC.Attempt();
        Drift.Attempt();
       	FinishParticles.Attempt();
        FinishMultipoles.Attempt();
       	SendManifest->FreeAfterSend();
    	ReceiveManifest->Check();   // This checks if Send is ready; no-op in non-blocking mode
    // If the manifest has been received, install it.
    	if (ReceiveManifest->is_ready()) ReceiveManifest->ImportData();
        AttemptNeighborReceive(0,P.cpd);  // 2D
        MF->CheckAnyMPIDone();
   		CheckForMultipoles.Attempt();
		
#ifdef PARALLEL
		timestep_loop_complete = CheckForMultipoles.alldone(total_slabs_on_node);
#else
		timestep_loop_complete = FinishMultipoles.alldone(total_slabs_on_node);
#endif
    }


    STDLOG(1, "Read %d particles from IC files\n", NP_from_IC);
    #ifdef PARALLEL
        BarrierWallClock.Start();
        MPI_REDUCE_TO_ZERO(&merged_particles, 1, MPI_UINT64_T, MPI_SUM);
        MPI_REDUCE_TO_ZERO(&NP_from_IC, 1, MPI_UINT64_T, MPI_SUM);
        STDLOG(1,"Ready to proceed to the remaining work\n");
        
        MPI_Barrier(comm_global);
        BarrierWallClock.Stop(); 
        
        // This MPI call also forces a synchronization over the MPI processes, 
        // so things like Reseting GPUs could fire multiple times on one node.
        SendManifest->FreeAfterSend();
        // Run this again, just in case the dependency loop on this node finished
       // before the neighbor received the non-blocking MPI transfer.
	    TimeStepWallClock.Stop(); ConvolutionWallClock.Start(); 
    	delete ParallelConvolveDriver;
	    ConvolutionWallClock.Stop(); TimeStepWallClock.Start(); 
		
    #endif
    STDLOG(1, "Particles remaining on insert list: %d\n", IL->length);
    if (MPI_rank==0) {
        STDLOG(1, "Merged %d particles\n", merged_particles);
        assertf(merged_particles == P.np, "Merged slabs contain %d particles instead of %d!  Read %d IC particles.\n",
            merged_particles, P.np, NP_from_IC);
    }

    if(IL->length!=0)
        IL->DumpParticles();

    if(MPI_rank == 0)
        assertf(NP_from_IC == P.np, "Expected to read a total of %u particles from IC files, but only read %u.\n", P.np, NP_from_IC);

    assertf(IL->length==0,
        "Insert List not empty (%d) at the end of timestep().  Particles in IC files not sufficiently sorted?\n", IL->length);
    
    char filename[1024];
    int ret = snprintf(filename,1024,"%s/slabsize",P.WriteStateDirectory);
    assert(ret >= 0 && ret < 1024);
    SS->write(filename);

    STDLOG(1,"Completing timestepIC()\n");
    TimeStepWallClock.Stop();
}
