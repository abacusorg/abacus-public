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

uint64 GatherTaggableFieldParticles(int slab, RVfloat *pv, TaggedPID *pid, FLOAT unkickfactor) {
    slab = GFC->WrapSlab(slab);
    uint64 nfield = 0;
    for (int j=0; j<GFC->cpd; j++)
        for (int k=0; k<GFC->cpd; k++) {
            // Loop over cells
            posstruct offset = CP->CellCenter(slab, j, k);
            Cell c = CP->GetCell(slab, j, k);
            for (int p=0; p<c.count(); p++)
                if (c.aux[p].is_taggable() && !c.aux[p].is_L0()) {
                    // We found a taggable field particle
                    posstruct r = c.pos[p] + offset;
                    velstruct v = c.vel[p];
            if(c.acc != NULL)
                v -= unkickfactor*TOFLOAT3(c.acc[p]);
                    pv[nfield] = RVfloat(r.x, r.y, r.z, v.x, v.y, v.z);
                    pid[nfield] = c.aux[p].pid();
                    nfield++;
                }
        }
    return nfield;
}

void OutputNonL0Taggable(int slab) {
    // Write out the taggable particles not in 01 halos.
    // This has to get called after all GlobalGroups in this slab
    // have been found.

    // TODO: better heuristic? what will happen in very small sims?  
    // Also technically HaloTaggableFraction is only used in the IC step
    uint64 maxsize = SS->size(slab)*P.HaloTaggableFraction;
            maxsize += 6*sqrt(maxsize);  // 6-sigma buffer

    SB->AllocateSpecificSize(TaggableFieldSlab, slab, maxsize*sizeof(RVfloat));
    SB->AllocateSpecificSize(TaggableFieldPIDSlab, slab, maxsize*sizeof(TaggedPID));

    uint64 nfield = GatherTaggableFieldParticles(slab,
            (RVfloat *) SB->GetSlabPtr(TaggableFieldSlab, slab),
            (TaggedPID *) SB->GetSlabPtr(TaggableFieldPIDSlab, slab),
            WriteState.FirstHalfEtaKick);
    if(nfield > 0){
        // only write the uniform subsample files if they will have non-zero size
        SB->ResizeSlab(TaggableFieldSlab, slab, nfield*sizeof(RVfloat));
        SB->ResizeSlab(TaggableFieldPIDSlab, slab, nfield*sizeof(TaggedPID));
        SB->StoreArenaNonBlocking(TaggableFieldSlab, slab);
        SB->StoreArenaNonBlocking(TaggableFieldPIDSlab, slab);
    } else {
        SB->DeAllocate(TaggableFieldSlab, slab);
        SB->DeAllocate(TaggableFieldPIDSlab, slab);
    }
    STDLOG(1,"Writing %d non-L0 Taggable particles in slab %d\n", nfield, slab);
}
