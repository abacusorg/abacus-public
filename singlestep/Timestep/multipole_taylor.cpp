// These routines provide the interface to the multipole and taylor
// evaluation modules.

// TODO: The SlabTaylor and SlabMultipole lists have inconsistent arguments.
// Taylors has long long, Multipoles has ints.  Need to beware of spilling 32-bit
// integers on the particle lists.

void ComputeTaylorForce(int slab) {
    int cpd = PP->cpd;
    int *count = new  int[cpd*cpd];
    int *offset = new int[cpd*cpd];
    FLOAT3 *cc = (FLOAT3 *) malloc(cpd*cpd*sizeof(FLOAT3));

    #pragma omp parallel for schedule(dynamic,1)
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

    // Recall that our accelerations/positions may be in single-precision.
    posstruct *pos = (posstruct *) LBW->ReturnIDPtr(PosSlab,slab);
    accstruct *acc = (accstruct *) LBW->ReturnIDPtr(AccSlab,slab);
    uint64 slabsize;

    Complex *TaylorCoefficients = (Complex *) LBW->ReturnIDPtr(TaylorSlab,slab);
    TY->EvaluateSlabTaylor( slab, acc, pos, count, offset, cc, TaylorCoefficients);//,Accslab_mutex_lock);
    RL->ApplyRedlack(slab, acc, pos, count, offset, cc,P.np);//,Accslab_mutex_lock);

    //wake up the near force
    //FIXME: May not be necessary any more since we block on outputting the accelerations 
    //JJ->wakeup();

    delete[] count;
    delete[] offset;
    free(cc);
}

void ComputeMultipoleSlab(int slab) {
    // This routine must use the Merged slabs!
    STDLOG(1,"Beginning to compute multipoles for slab %d\n", slab);
    ComputeMultipoles.Start();
    int cpd = PP->cpd;
    int *count = new int[cpd*cpd];
    int *offset = new int[cpd*cpd];
    FLOAT3 *cc = (FLOAT3 *) malloc(cpd*cpd*sizeof(FLOAT3));

    #pragma omp parallel for schedule(dynamic,1)
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

    // Recall that our accelerations/positions may be in single-precision.
    posstruct *pos = (posstruct *) LBW->ReturnIDPtr(MergePosSlab,slab);
    uint64 slabsize;
    STDLOG(1,"Calling multipole module.\n");
    Complex *slabmultipoles = (Complex *) LBW->ReturnIDPtr(MultipoleSlab,slab);
    MF->ComputeMultipoleFFTYZ( slab, pos, count, offset, cc, slabmultipoles);


    delete[] count;
    delete[] offset;
    free(cc);
    ComputeMultipoles.Stop();
    STDLOG(1,"Done with multipoles.\n");
}
