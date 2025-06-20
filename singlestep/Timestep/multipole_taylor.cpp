// Copyright 2012-2025 The Abacus Developers
// SPDX-License-Identifier: GPL-3.0-or-later

/// These routines provide the interface to the multipole and taylor
/// evaluation modules.

// TODO: The SlabTaylor and SlabMultipole lists have inconsistent arguments.
// Taylors has long long, Multipoles has ints.  Need to beware of spilling 32-bit
// integers on the particle lists.

#ifdef PARALLEL
int ramdisk_multipole_flag = RAMDISK_NO;
#else
int ramdisk_multipole_flag = RAMDISK_AUTO;
#endif

/// The driver routine to compute the Taylors on one slab
void ComputeTaylorForce(int slab) {
    TY->ConstructOffsets.Start();
    int cpd = CP->cpd;

    // TODO: the concept of a per-pencil ghost offset is fairly general. Where to store?
    int *ghost_offsets = new int[cpd];

    #pragma omp parallel for schedule(static)
    for(int y=0;y<cpd;y++){
        
        cellinfo *startci = CP->CellInfo(slab,y,node_z_start);
        ghost_offsets[y] = startci->startindex_with_ghost - startci->startindex;
        
        for(int z = 0; z < node_z_size; z++) {
            cellinfo *ciptr = CP->CellInfo(slab,y,node_z_start + z);
            int i = y*node_z_size + z;
            // startindex indexes accs;
            // startindex + ghost_offsets[y] indexes pos
            TY->offset[i] = ciptr->startindex;
            TY->count[i] =  ciptr->count;
            // For box-centered, this will be the cell center.
            // For local cell-centered positions, this will return 0!
            TY->cc[i] = CP->LocalCellCenter(slab,y,node_z_start + z);
        }
    }
    TY->ConstructOffsets.Stop();
    
    // Recall that our accelerations/positions may be in single-precision.
    posstruct *pos = (posstruct *) SB->GetSlabPtr(PosSlab,slab);
    acc3struct *acc = (acc3struct *) SB->GetSlabPtr(FarAccSlab,slab);

    MTCOMPLEX *TaylorCoefficients = (MTCOMPLEX *) SB->GetSlabPtr(TaylorSlab,slab);
    TY->EvaluateSlabTaylor( slab, acc, pos, TY->count, TY->offset, ghost_offsets,
                            TY->cc, TaylorCoefficients);
    RL->ApplyRedlack(slab, acc, pos, TY->count, TY->offset, ghost_offsets,
                            TY->cc,P.np);

    delete[] ghost_offsets;
}

/// The driver routine to compute the Multipoles on one slab
void ComputeMultipoleSlab(int slab) {
    // This routine must use the Merged slabs!
    STDLOG(1,"Computing multipoles for slab {:d}\n", slab);
    ComputeMultipoles.Start();
    MF->ConstructOffsets.Start();
    int cpd = CP->cpd;

    #pragma omp parallel for schedule(static)
    for(int y=0;y<cpd;y++)
        for(int z = 0; z < node_z_size; z++) {
	        cellinfo *ciptr = CP->MergeCellInfo(slab, y, node_z_start + z);
            int i = y*node_z_size + z;
            MF->offset[i] = ciptr->startindex_with_ghost;
            MF->count[i] =  ciptr->count;
            // For box-centered, this will be the cell center.
            // For local cell-centered positions, this will return 0!
            MF->cc[i] = CP->LocalCellCenter(slab, y, node_z_start + z);
        }
    MF->ConstructOffsets.Stop();
    
    MTCOMPLEX *slabmultipoles = NULL;
    if(MPI_size_z == 1){
        // Only allocate if it's going to be used
        slabmultipoles = (MTCOMPLEX *) SB->AllocateArena(MultipoleSlab, slab, ramdisk_multipole_flag);
    }
    
    // Recall that our accelerations/positions may be in single-precision.
    posstruct *pos = (posstruct *) SB->GetSlabPtr(MergePosSlab,slab);
    STDLOG(3,"Calling multipole module.\n");

    // This will do 1D or 2D FFT, depending on whether we've split z
    MF->ComputeMultipoleFFT( slab, pos, MF->count, MF->offset, MF->cc, slabmultipoles);

    ComputeMultipoles.Stop();
    STDLOG(1,"Done with multipoles.\n");
}
