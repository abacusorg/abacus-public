// ================  Routines for the non-L0 taggable output =============

/** Gather all of the taggable particles that aren't in L0 groups into
two vectors, converting to global positions.

Space must be allocated beforehand.  Returns the number of elements used.
Warning: This must be called after ScatterGlobalGroupsAux() and before
ScatterGlobalGroups(), because the latter may contain microstepped 
updates of the positions and velocities.
*/

// These use data models from halostats.hh.
// TODO: Consider refactoring that?

void GatherTaggableFieldParticles(int slab, SlabAccum<RVfloat> ** rv, SlabAccum<TaggedPID> ** pid, FLOAT unkickfactor, uint64 * nfield) {
    slab = GFC->WrapSlab(slab);
    double vel_convert_units = ReadState.VelZSpace_to_kms/ReadState.VelZSpace_to_Canonical; 
    int nfield_A = nfield_B = 0; 

    #pragma omp parallel for schedule(dynamic,1) reduction(+: nfield_A, nfield_B)
    for (int j=0; j<GFC->cpd; j++){

        if ( rv[0] != NULL) PencilAccum<RVfloat>   *pRVA   =  rv[0].StartPencil(j);
        if ( rv[1] != NULL) PencilAccum<RVfloat>   *pRVB   =  rv[1].StartPencil(j);
        if (pid[0] != NULL) PencilAccum<TaggedPID> *pPIDsA = pid[0].StartPencil(j);
        if (pid[1] != NULL) PencilAccum<TaggedPID> *pPIDsB = pid[1].StartPencil(j);

        for (int k=0; k<GFC->cpd; k++) {
            // Loop over cells
            posstruct offset = CP->CellCenter(slab, j, k);
            Cell c = CP->GetCell(slab, j, k);
            for (int p=0; p<c.count(); p++){
                int tag = c.aux[p].is_taggable(); //0 non taggable, TAGGABLE_SUB_A or TAGGABLE_SUB_B for the subsamples. 
                if (tag > 0 && !c.aux[p].is_L0()) {
                    // We found a taggable field particle
                    posstruct r = c.pos[p] + offset;
                    velstruct v = c.vel[p];
                    if(c.acc != NULL){ v -= unkickfactor*TOFLOAT3(c.acc[p]);}

                    v *= vel_convert_units;  
                      
                    if (tag == TAGGABLE_SUB_A){
                        assertf(rv[0] != NULL and pid[0] != NULL, 
                            "This particle is taggable in subset A, but OutputNonL0Taggable thinks P.ParticleSubsampleA should be zero.\n");
                        pRVA->append(RVfloat(r.x, r.y, r.z, v.x, v.y, v.z));
                        pPIDsA->append(TaggedPID(c.aux[p]);)
                        nfield_A++;
                    }
                    else if (tag == TAGGABLE_SUB_B){
                        assertf(rv[1] != NULL and pid[1] != NULL, 
                            "This particle is taggable in subset B, but OutputNonL0Taggable thinks P.ParticleSubsampleB should be zero.\n");
                        pRVB->append(RVfloat(r.x, r.y, r.z, v.x, v.y, v.z));
                        pPIDsB->append(TaggedPID(c.aux[p]));
                        nfield_B++;
                    }
                }
            }
            pRVA->FinishCell();
            pRVB->FinishCell();
            pPIDsA->FinishCell();
            pPIDsB->FinishCell();
        }

        pRVA->FinishPencil();
        pRVB->FinishPencil();
        pPIDsA->FinishPencil();
        pPIDsB->FinishPencil();
    }

    nfield[0] = nfield_A; nfield[1] = nfield_B; 
}

void OutputNonL0Taggable(int slab) {
    // If subsampling output is requested,
    // Write out the taggable particles not in 01 halos.
    // This has to get called after all GlobalGroups in this slab
    // have been found.

    double subsample_fracs[NUM_SUBSAMPLES] = {P.ParticleSubsampleA, P.ParticleSubsampleB}; 
    int slab_type[2*NUM_SUBSAMPLES] = {FieldRVSlabA, FieldRVSlabB, FieldPIDSlabA, FieldPIDSlabB}; 

    //RVfloat   *  rvSlabs[NUM_SUBSAMPLES];
    //TaggedPID * pidSlabs[NUM_SUBSAMPLES];

    SlabAccum<RVfloat>   rvA;    
    SlabAccum<RVfloat>   rvB; 
    SlabAccum<TaggedPID> pidA;  
    SlabAccum<TaggedPID> pidB; 

    if (P.ParticleSubsampleA > 0) {
        rvA.setup( CP->cpd, P.np/P.cpd*P.ParticleSubsampleA);   
        pidA.setup(CP->cpd, P.np/P.cpd*P.ParticleSubsampleA);   
    }
    else { rvA = NULL; pidA = NULL;}

    if (P.ParticleSubsampleB > 0) {
        rvB.setup( CP->cpd, P.np/P.cpd*P.ParticleSubsampleB); //TODO: what's a good guess here? 
        pidB.setup(CP->cpd, P.np/P.cpd*P.ParticleSubsampleB); 
    }
    else { rvB = NULL; pidB = NULL;}

    // for (int i = 0; i < NUM_SUBSAMPLES; i++){
    //         // TODO: better heuristic? what will happen in very small sims?  
    //     if (subsample_fracs[i] > 0) {
    //         // uint64 maxsize;
    //         // if (subsample_fracs[i] > 0) maxsize = SS->size(slab) * subsample_fracs[i]; 
    //         // else maxsize = SS->size(slab) * 0.1; //dummy size. 
    //         // maxsize += 6*sqrt(maxsize); // 6-sigma buffer

    //         // SB->AllocateSpecificSize(slab_type[i],                slab, maxsize*sizeof(RVfloat));
    //         // SB->AllocateSpecificSize(slab_type[i+NUM_SUBSAMPLES], slab, maxsize*sizeof(TaggedPID));

    //         // rvSlabs[i]   = (RVfloat*)   SB->GetSlabPtr(slab_type[i]  , slab);
    //         // pidSlabs[i] = (TaggedPID*) SB->GetSlabPtr(slab_type[i+2], slab);
    //     }
    //     else {
    //         STDLOG(4, "Setting rvSlabs and pidSlabs to NULL for %d\n", i);
    //         rvSlabs[i] = NULL;
    //         pidSlabs[i] = NULL; 
    //     }
    // }

    uint64 nfield[NUM_SUBSAMPLES] = {0, 0}; 

    SlabAccum<RVfloat> *  rv[NUM_SUBSAMPLES];    
    SlabAccum<RVfloat> * pid[NUM_SUBSAMPLES]; 
    GatherTaggableFieldParticles(slab, rv, pid, WriteState.FirstHalfEtaKick, nfield);

    WriteState.np_subA_state += nfield[0]; 
    WriteState.np_subB_state += nfield[1]; 

    //NAM stopped here. 
    for (int i = 0; i < NUM_SUBSAMPLES; i++){
        if(subsample_fracs[i] > 0){
            if(nfield[i] > 0){
                // only write the uniform subsample files if they will have non-zero size
                SB->ResizeSlab(slab_type[i],                slab, nfield[i]*sizeof(RVfloat));
                SB->ResizeSlab(slab_type[i+NUM_SUBSAMPLES], slab, nfield[i]*sizeof(TaggedPID));
                SB->StoreArenaNonBlocking(slab_type[i], slab);
                SB->StoreArenaNonBlocking(slab_type[i+NUM_SUBSAMPLES], slab);
            } else {
                SB->DeAllocate(slab_type[i], slab);
                SB->DeAllocate(slab_type[i+NUM_SUBSAMPLES], slab);
            }
        }
        STDLOG(1,"Writing %d non-L0 Taggable particles in slab %d in subsample %d of %d.\n", nfield[i], slab, i+1, NUM_SUBSAMPLES);
    }

     rvA.destroy();
     rvB.destroy();
    pidA.destroy();
    pidB.destroy();

}
