/* timestep_standalone_fof.cpp

This file contains a minimal slab pipeline that runs on-the-fly
group finding on a pack14 output.  It is primarily invoked
from standalone_fof.cpp via the standalone_fof.py interface.

This file is directly #include'd in timestep.cpp since it borrows a 
lot of infrastructure from there.

*/

#include "read_pack14.cpp"

const char* StandaloneFOF_slice_dir;
int StandaloneFOFLoadSlabPrecondition(int slab) {
    if(slab > Finish.last_slab_executed + 2*GFC->GroupRadius + FETCHAHEAD)
        return 0;
    return 1;
}

void StandaloneFOFLoadSlabAction(int slab) {
    char fname[1024];
    // TODO: Add support for L0 slabs?
    sprintf(fname, "%s/%s.z%5.3f.slab%04d.dat", StandaloneFOF_slice_dir, P.SimName, ReadState.Redshift, slab);
    STDLOG(1,"Load Slab %d from \"%s\"\n", slab, fname);

    size_t s = fsize(fname);
    SB->AllocateSpecificSize(TimeSlice, slab, s);
    // We will read the raw pack14 asynchronously with SB
    // then unpack it in a separate dependency
    // TODO: support states as well as time slices
    // TODO: ignores ramdisk
    SB->ReadArena(TimeSlice, slab, IO_NONBLOCKING, fname);
}

int StandaloneFOFUnpackSlabPrecondition(int slab) {
    if (! SB->IsIOCompleted(TimeSlice, slab)) return 0;
    return 1;
}

void StandaloneFOFUnpackSlabAction(int slab) {
    printf("Unpacking slab %d\n", slab);
    STDLOG(1, "Unpacking slab %d\n", slab);
    int nump = unpack_slab_pack14(slab, P.HaloTaggableFraction);
    STDLOG(1,"Found %d particles in slab %d\n", nump, slab);

    SB->DeAllocate(TimeSlice, slab);
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

    SB->DeAllocate(PosSlab, slab);
    SB->DeAllocate(VelSlab, slab);
    SB->DeAllocate(AuxSlab, slab);
    SB->DeAllocate(CellInfoSlab, slab);
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