/*
 * This class creates an AppendArena object of the desired output type
 * and repeatedly calls "addparticle()" on it.  The actual implementation
 * of the writing functions is in appendarena.cpp.
 *
*/

#include "appendarena.cpp"

AppendArena *get_AA_by_format(const char* format){
    AppendArena *AA;

    if (strcmp(format,"RVdouble")==0) {
        STDLOG(1,"Using Output Format RVdouble\n");
        AA = new OutputRVdouble();
    } else if (strcmp(format,"Packed")==0) {
        STDLOG(1,"Using Output Format Packed\n");
        AA = new OutputPacked();
    } else if (strcmp(format,"Heitmann")==0) {
        STDLOG(1,"Using Output Format Heitmann\n");
        AA = new OutputHeitmann();
    } else if (strcmp(format,"RVdoubleTag")==0) {
        STDLOG(1,"Using Output Format RVdoubleTag\n");
        AA = new OutputRVdoubleTag();
    } else if (strcmp(format,"RVZel")==0) {
        STDLOG(1,"Using Output Format RVZel\n");
        AA = new OutputRVZel();
    }
    else {
        QUIT("Unrecognized case: OutputFormat = %s\n", format);
    }
    
    return AA;
}

void WriteHeaderFile(const char* fn){
	std::ofstream headerfile;
	headerfile.open(fn);
	headerfile << P.header();
	headerfile << ReadState.header();
	headerfile << "\nOutputType = \"TimeSlice\"\n";
	headerfile.close();
}

uint64 Output_TimeSlice(int slab, FLOAT unkickfactor) {
    AppendArena *AA;
    FLOAT vscale;

    AA = get_AA_by_format(P.OutputFormat);

    // Setup the Arena
    int headersize = 1024*1024;
    SB->AllocateSpecificSize(TimeSlice, slab, 
    	   SS->size(slab)*(AA->sizeof_particle())
	+ CP->cpd*(CP->cpd)*(AA->sizeof_cell()) + headersize);
    AA->initialize(TimeSlice, slab, CP->cpd, ReadState.VelZSpace_to_Canonical);

    // Write the header to its own file
    if(slab == 0){
        char filename[1024];
        sprintf(filename, "%s/slice%5.3f/header",  
            P.OutputDirectory, 
            ReadState.Redshift);
        WriteHeaderFile(filename);
    }
        
    // and also add the header to the slab file
    if (!P.OmitOutputHeader) {
        AA->addheader((const char *) P.header());
        AA->addheader((const char *) ReadState.header());
        char head[1024];
        sprintf(head, "\nOutputType = \"TimeSlice\"\n"); 
        AA->addheader((const char *) head);
        sprintf(head, "SlabNumber = %d\n", slab);
        AA->addheader((const char *) head);
        // For sanity, be careful that the previous lines end with a \n!
        AA->finalize_header();
    }

    // Now scan through the cells
    velstruct vel;
    integer3 ijk(slab,0,0);
    uint64 n_added = 0;
    for (ijk.y=0; ijk.y<CP->cpd; ijk.y++) 
        for (ijk.z=0;ijk.z<CP->cpd;ijk.z++) {
            Cell c = CP->GetCell(ijk);

            // We sometimes use the maximum velocity to scale.
            // But we do not yet have the global velocity (slab max will be set in Finish,
            // while the global max has to wait for all slabs to be done).
            // What is available after the kick is the max_component_velocity in each cell.
            vscale = c.ci->max_component_velocity/ReadState.VelZSpace_to_Canonical;	
            // The maximum velocity of this cell, converted to ZSpace unit-box units.

            // Start the cell
            AA->addcell(ijk, vscale);
            
            // Now pack the particles
            accstruct *acc = CP->AccCell(ijk);
            for (int p=0;p<c.count();p++) {
                vel = (c.vel[p] - TOFLOAT3(acc[p])*unkickfactor);    // We supply in code units
                // Detail: we write particles with their L0 bits intact.  So if we want to run a non-group-finding step
                // after a group-finding step (e.g. for debugging), we need to know that we can ignore the L0 bit
                if(GFC == NULL || !c.aux[p].is_L0()){
                    AA->addparticle(c.pos[p], vel, c.aux[p]);
                    n_added++;
                }
            }
            AA->endcell();
        }

    SB->ResizeSlab(TimeSlice, slab, AA->bytes_written());

    // Write out this time slice
    SB->StoreArenaNonBlocking(TimeSlice, slab);
    delete AA;
    
    return n_added;
}


// ================  Routines for the non-L1 taggable output =============

/** Gather all of the taggable particles that aren't in L1 groups into
two vectors, converting to global positions.

Space must be allocated beforehand.  Returns the number of elements used.
Warning: This must be called after ScatterGlobalGroupsAux() and before
ScatterGlobalGroups(), because the latter may contain microstepped 
updates of the positions and velocities.
*/

uint64 GatherTaggableFieldParticles(int slab, RVfloat *pv, TaggedPID *pid, FLOAT unkickfactor) {
    slab = GFC->WrapSlab(slab);
    uint64 nfield = 0;
    for (int j=0; j<GFC->cpd; j++)
        for (int k=0; k<GFC->cpd; k++) {
            // Loop over cells
            posstruct offset = CP->CellCenter(slab, j, k);
            Cell c = CP->GetCell(slab, j, k);
            for (int p=0; p<c.count(); p++)
                if (c.aux[p].is_taggable() && !c.aux[p].is_L1()) {
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

void OutputNonL1Taggable(int slab) {
    // Write out the taggable particles not in L1 halos.
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

}
