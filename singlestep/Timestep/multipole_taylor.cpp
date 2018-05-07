// These routines provide the interface to the multipole and taylor
// evaluation modules.

// TODO: The SlabTaylor and SlabMultipole lists have inconsistent arguments.
// Taylors has long long, Multipoles has ints.  Need to beware of spilling 32-bit
// integers on the particle lists.

void ComputeTaylorForce(int slab) {
    TY->ConstructOffsets.Start();
    int cpd = PP->cpd;
    int *count = new  int[cpd*cpd];
    int *offset = new int[cpd*cpd];
    FLOAT3 *cc;
    assert(posix_memalign((void **) &cc, 4096, cpd*cpd*sizeof(FLOAT3)) == 0);

    #pragma omp parallel for schedule(static)
    for(int y=0;y<cpd;y++)
        for(int z=0;z<cpd;z++) {
            cellinfo *ciptr = PP->CellInfo(slab,y,z);
            int i = y*cpd + z;
            offset[i] = ciptr->startindex;
            count[i] =  ciptr->count;
            // For box-centered, this will be the cell center.
            // For local cell-centered positions, this will return 0!
            cc[i] = PP->LocalCellCenter(slab,y,z);
        }
    TY->ConstructOffsets.Stop();
    
    // Recall that our accelerations/positions may be in single-precision.
    posstruct *pos = (posstruct *) LBW->ReturnIDPtr(PosSlab,slab);
    accstruct *acc = (accstruct *) LBW->ReturnIDPtr(AccSlab,slab);
    uint64 slabsize;

    MTCOMPLEX *TaylorCoefficients = (MTCOMPLEX *) LBW->ReturnIDPtr(TaylorSlab,slab);
    TY->EvaluateSlabTaylor( slab, acc, pos, count, offset, cc, TaylorCoefficients);
    RL->ApplyRedlack(slab, acc, pos, count, offset, cc,P.np);

    delete[] count;
    delete[] offset;
    free(cc);
}

void ComputeMultipoleSlab(int slab) {
    // This routine must use the Merged slabs!
    STDLOG(1,"Beginning to compute multipoles for slab %d\n", slab);
    ComputeMultipoles.Start();
    MF->ConstructOffsets.Start();
    int cpd = PP->cpd;
    int *count = new int[cpd*cpd];
    int *offset = new int[cpd*cpd];
    FLOAT3 *cc;
    assert(posix_memalign((void **) &cc, 4096, cpd*cpd*sizeof(FLOAT3)) == 0);

    #pragma omp parallel for schedule(static)
    for(int y=0;y<cpd;y++)
        for(int z=0;z<cpd;z++) {
	        cellinfo *ciptr = PP->MergeCellInfo(slab,y,z);
            int i = y*cpd + z;
            offset[i] = ciptr->startindex;
            count[i] =  ciptr->count;
            // For box-centered, this will be the cell center.
            // For local cell-centered positions, this will return 0!
            cc[i] = PP->LocalCellCenter(slab,y,z);
        }
    MF->ConstructOffsets.Stop();

    // Recall that our accelerations/positions may be in single-precision.
    posstruct *pos = (posstruct *) LBW->ReturnIDPtr(MergePosSlab,slab);
    uint64 slabsize;
    STDLOG(1,"Calling multipole module.\n");
    MTCOMPLEX *slabmultipoles = (MTCOMPLEX *) LBW->ReturnIDPtr(MultipoleSlab,slab);
    MF->ComputeMultipoleFFTYZ( slab, pos, count, offset, cc, slabmultipoles);

    delete[] count;
    delete[] offset;
    free(cc);
    ComputeMultipoles.Stop();
    STDLOG(1,"Done with multipoles.\n");
}

// Thin wrappers to convert from arenas to pointers
void DoCellMultipoles(int slab, int y, int z){ return;};
void DoZRowMultipoles(int slab, int y){ return;};
void DoYRowMultipoles(int slab){ return;};
