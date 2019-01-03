/* timestep_benchmark_io.cpp

This file contains a minimal slab pipeline that just reads
slabs from the read state and writes them to the write state.
It is useful for performing realistic IO benchmarks (although
it may be unrealistic in the sense that contention from other
parts of the code will not be present, to the extend that matters).

This pipeline is called from benchmark_io.cpp.  I don't think
we have a top-level wrapper for that right now; we just invoke
the executable manually when we want to benchmark the IO.

This file is directly #include'd in timestep.cpp since it borrows a 
lot of infrastructure from there.

*/

int FinishBenchmarkIOPrecondition(int slab) {
    // Wait for everything to be read
    if( !SB->IsIOCompleted( CellInfoSlab,      slab )
        || !SB->IsIOCompleted( PosSlab, slab )
        || !SB->IsIOCompleted( VelSlab, slab )
        || !SB->IsIOCompleted( AuxSlab,      slab )
        || !SB->IsIOCompleted( TaylorSlab,      slab )
        ) return 0;
    return 1;
}

void FinishBenchmarkIOAction(int slab) {
    STDLOG(1,"Finishing benchmark IO slab %d\n", slab);
        
    /*// Make the multipoles
    SB->AllocateArena(MultipoleSlab,slab);
    ComputeMultipoleSlab(slab);
    
    WriteMultipoleSlab.Start();
    SB->StoreArenaNonBlocking(MultipoleSlab,slab);
    WriteMultipoleSlab.Stop();
    
    SB->DeAllocate(MergePosSlab,slab);
    SB->DeAllocate(MergeCellInfoSlab,slab);*/

    // Ignore details about ramdisk and just write, I guess
    // TODO: the more canonical way to do this would be something like multipole recovery
    SB->WriteArena(CellInfoSlab, slab, IO_DELETE, IO_NONBLOCKING, 
                    SB->WriteSlabPath(MergeCellInfoSlab,slab).c_str());
    SB->WriteArena(PosSlab, slab, IO_DELETE, IO_NONBLOCKING, 
                    SB->WriteSlabPath(MergePosSlab,slab).c_str());
    SB->WriteArena(VelSlab, slab, IO_DELETE, IO_NONBLOCKING, 
                    SB->WriteSlabPath(MergeVelSlab,slab).c_str());
    SB->WriteArena(AuxSlab, slab, IO_DELETE, IO_NONBLOCKING, 
                    SB->WriteSlabPath(MergeAuxSlab,slab).c_str());
    SB->WriteArena(TaylorSlab, slab, IO_DELETE, IO_NONBLOCKING, 
                    SB->WriteSlabPath(MultipoleSlab,slab).c_str());
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

    int cpd = P.cpd; int first = first_slab_on_node;
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
