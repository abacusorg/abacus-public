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

void GatherTaggableFieldParticles(int slab, RVfloat ** pv, TaggedPID ** pid, FLOAT unkickfactor, uint64 * nfield) {
    slab = GFC->WrapSlab(slab);
    for (int j=0; j<GFC->cpd; j++)
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

                    if (tag == TAGGABLE_SUB_A){
                        pv[0][nfield[0]]  = RVfloat(r.x, r.y, r.z, v.x, v.y, v.z);
                        pid[0][nfield[0]] = c.aux[p].pid();
                        nfield[0]++;
                    }
                    else if (tag == TAGGABLE_SUB_B){
                        pv[1][nfield[1]]  = RVfloat(r.x, r.y, r.z, v.x, v.y, v.z);
                        pid[1][nfield[1]] = c.aux[p].pid();
                        nfield[1]++;
                    }
                }
            }
        }
}

void OutputNonL0Taggable(int slab) {
    // If subsampling output is requested,
    // Write out the taggable particles not in 01 halos.
    // This has to get called after all GlobalGroups in this slab
    // have been found.

    //NAM not happy with this... 

    double subsample_fracs[NUM_SUBSAMPLES] = {P.ParticleSubsampleA, P.ParticleSubsampleB}; 
    int slab_type[2*NUM_SUBSAMPLES] = {FieldRVSlabA, FieldRVSlabB, FieldPIDSlabA, FieldPIDSlabB}; 

    for (int i = 0; i < NUM_SUBSAMPLES; i++){
            // TODO: better heuristic? what will happen in very small sims?  
        uint64 maxsize = SS->size(slab)* subsample_fracs[i];
        maxsize += 6*sqrt(maxsize); // 6-sigma buffer

        SB->AllocateSpecificSize(slab_type[i],    slab, maxsize*sizeof(RVfloat));
        SB->AllocateSpecificSize(slab_type[i+NUM_SUBSAMPLES],  slab, maxsize*sizeof(TaggedPID));
    }

    uint64 nfield[NUM_SUBSAMPLES] = {0, 0}; 
    RVfloat   *  rvSlabs[NUM_SUBSAMPLES] = {(RVfloat*)   SB->GetSlabPtr(slab_type[0], slab), (RVfloat*)   SB->GetSlabPtr(slab_type[1], slab)};
    TaggedPID * pidSlabs[NUM_SUBSAMPLES] = {(TaggedPID*) SB->GetSlabPtr(slab_type[2], slab), (TaggedPID*) SB->GetSlabPtr(slab_type[3], slab)};

    GatherTaggableFieldParticles(slab, rvSlabs, pidSlabs, WriteState.FirstHalfEtaKick, nfield);

    for (int i = 0; i < NUM_SUBSAMPLES; i++){
        if(nfield[i] > 0){
                // only write the uniform subsample files if they will have non-zero size
            SB->ResizeSlab(slab_type[i], slab, nfield[i]*sizeof(RVfloat));
            SB->ResizeSlab(slab_type[i+NUM_SUBSAMPLES], slab, nfield[i]*sizeof(TaggedPID));
            SB->StoreArenaNonBlocking(slab_type[i], slab);
            SB->StoreArenaNonBlocking(slab_type[i+NUM_SUBSAMPLES], slab);
        } else {
            SB->DeAllocate(slab_type[i], slab);
            SB->DeAllocate(slab_type[i+NUM_SUBSAMPLES], slab);
        }
        STDLOG(1,"Writing %d non-L0 Taggable particles in slab %d in subsample %d of %d.\n", nfield[i], slab, i, NUM_SUBSAMPLES);
    }
    

}
