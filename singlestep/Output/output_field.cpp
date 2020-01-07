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

void GatherTaggableFieldParticles(int slab, SlabAccum<RVfloat> * rv, SlabAccum<TaggedPID> * pid, FLOAT unkickfactor, uint64 * nfield) {
    slab = GFC->WrapSlab(slab);
    double vel_convert_units = ReadState.VelZSpace_to_kms/ReadState.VelZSpace_to_Canonical; 
    int nfield_A = 0; 
    int nfield_B = 0;

    #pragma omp parallel for schedule(dynamic,1) reduction(+: nfield_A, nfield_B)
    for (int j=0; j<GFC->cpd; j++){

        PencilAccum<RVfloat>   *pRVA  ;
        PencilAccum<RVfloat>   *pRVB  ;
        PencilAccum<TaggedPID> *pPIDsA;
        PencilAccum<TaggedPID> *pPIDsB;

        pRVA   =  rv[0].StartPencil(j);
        pRVB   =  rv[1].StartPencil(j);
        pPIDsA = pid[0].StartPencil(j);
        pPIDsB = pid[1].StartPencil(j);

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
                        assertf(P.ParticleSubsampleA != 0, 
                            "This particle is taggable in subset A, but P.ParticleSubsampleA is zero.\n");
                        pRVA->append(RVfloat(r.x, r.y, r.z, v.x, v.y, v.z));
                        pPIDsA->append(TaggedPID(c.aux[p]));
                        nfield_A++;
                    }
                    else if (tag == TAGGABLE_SUB_B){
                        assertf(P.ParticleSubsampleB != 0, 
                            "This particle is taggable in subset B, but P.ParticleSubsampleA is zero.\n");
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

    SlabAccum<RVfloat>   rvA;    
    SlabAccum<RVfloat>   rvB; 
    SlabAccum<TaggedPID> pidA;  
    SlabAccum<TaggedPID> pidB; 

     rvA.setup(CP->cpd, P.np/P.cpd*P.ParticleSubsampleA);   
    pidA.setup(CP->cpd, P.np/P.cpd*P.ParticleSubsampleA);   
     rvB.setup(CP->cpd, P.np/P.cpd*P.ParticleSubsampleB); 
    pidB.setup(CP->cpd, P.np/P.cpd*P.ParticleSubsampleB); 

    uint64 nfield[NUM_SUBSAMPLES] = {0, 0}; 
    SlabAccum<RVfloat>    rv[NUM_SUBSAMPLES] = { rvA,  rvB};    
    SlabAccum<TaggedPID> pid[NUM_SUBSAMPLES] = {pidA, pidB}; 
    GatherTaggableFieldParticles(slab, rv, pid, WriteState.FirstHalfEtaKick, nfield);

    WriteState.np_subA_state += nfield[0]; 
    WriteState.np_subB_state += nfield[1]; 

    if (P.ParticleSubsampleA > 0){
        SB->AllocateSpecificSize(FieldRVSlabA, slab, rvA.get_slab_bytes());
        rvA.copy_to_ptr((RVfloat *)SB->GetSlabPtr(FieldRVSlabA, slab));
        SB->StoreArenaBlocking(FieldRVSlabA, slab); //NAM TODO: temporarily turned off nonblocking writes. 

        SB->AllocateSpecificSize(FieldPIDSlabA, slab, pidA.get_slab_bytes());
        pidA.copy_to_ptr((TaggedPID *)SB->GetSlabPtr(FieldPIDSlabA, slab));
        SB->StoreArenaBlocking(FieldPIDSlabA, slab);
    }
    if (P.ParticleSubsampleB > 0) {
        SB->AllocateSpecificSize(FieldRVSlabB, slab, rvB.get_slab_bytes());
        rvB.copy_to_ptr((RVfloat *)SB->GetSlabPtr(FieldRVSlabB, slab));
        SB->StoreArenaBlocking(FieldRVSlabB, slab);

        SB->AllocateSpecificSize(FieldPIDSlabB, slab, pidB.get_slab_bytes());
        pidB.copy_to_ptr((TaggedPID *)SB->GetSlabPtr(FieldPIDSlabB, slab));
        SB->StoreArenaBlocking(FieldPIDSlabB, slab);
    }

    STDLOG(1,"Wrote %d, %d non-L0 Taggable particles in subsamples A, B for slab %d.\n", nfield[0], nfield[1], slab);
}



