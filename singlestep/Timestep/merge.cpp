// Copyright 2012-2025 The Abacus Developers
// SPDX-License-Identifier: GPL-3.0-or-later

/* merge.cpp
 *
 *This routine combines the particles surviving in a slab (post-drift)
 *with the new particles arriving from the insert list.  We do this
 *by inspecting the insert list to build a tabulation of the number of
 *new particles per cell, then allocating a new arena and copying 
 *both old and new into the appropriate locations.
 *
 *This is one of the few places that arenas are 
 *created or destroyed outside of the timestep.cpp file.  Otherwise
 *we try to keep that behavior explicit in timestep.
 *
 *We also aggregate the cell-based statistics into global statistics
 *in this function, placing the results in WriteState so that it can
 *be output (recall that the statistics are from the ReadState time,
 *however!).
 *
 */

// Here are some shorthand functions.  
// Beware that these will evaluate 'newvalue' twice.
#define TRACK_MAX(maxsofar,newvalue) if ((newvalue)>maxsofar) maxsofar=newvalue
#define TRACK_MIN(minsofar,newvalue) if ((newvalue)<minsofar) minsofar=newvalue


class SkewerIndex {
public:
    // We need to collect some items for one skewer's worth of indexing,
    // so that we can multi-thread easily.
    // In the 2D code, a "skewer" spans the local z domain (with ghosts),
    // rather than [0,cpd).
    uint64 ilskewerlength;    // How many items in this skewer
    uint64 ilskewerstart;     // Where the indexing in the full list starts

    uint64 ilprimarystart;    // Where the indexing of the primary particles starts
    uint64 ilprimarylength;   // How many items in this skewer, without ghosts

    uint64 activelength;    // How many active original particles in this skewer

    uint64 mergelength;     // How many particles in the Merged skewer
    uint64 mergestart;      // Where the Merged skewer will start
    uint64 mergestart_no_ghost;  // Indexing skipping ghosts

    // The stats to accumulate
    int maxcellsize, mincellsize;
    double stddev_cellsize;
    double sum_square_velocity;
    FLOAT max_velocity;
    FLOAT max_acceleration;
    FLOAT min_vrms_on_amax;
    
    // A helper variable
    double mean_cellsize;

    SkewerIndex() {
        // Initialize the per-skewer statistics
        maxcellsize = 0;
        mincellsize = 1e9;
        mean_cellsize = P.np*CP->invcpd3;
        stddev_cellsize = 0.0;
        sum_square_velocity = 0.0;
        max_velocity = 0;
        max_acceleration = 0;
        min_vrms_on_amax = 1e30;
    }
    ~SkewerIndex() { }

    void stats(cellinfo *ci, cellinfo *mci) {
        // Track some statistics
        TRACK_MAX(maxcellsize, mci->count);
        TRACK_MIN(mincellsize, mci->count);
        TRACK_MAX(max_acceleration, ci->max_component_acceleration);
        TRACK_MAX(max_velocity, ci->max_component_velocity);

        if (ci->max_component_acceleration>0) {
            FLOAT vrms_on_amax = sqrt(ci->mean_square_velocity)/
                                    ci->max_component_acceleration;
            TRACK_MIN(min_vrms_on_amax, vrms_on_amax);
        }
        double res = mci->count/mean_cellsize-1;
        stddev_cellsize += res*res;                // Fractional overdensity
        sum_square_velocity += (double) ci->mean_square_velocity*(ci->count);
    }        

    void addstats(SkewerIndex &s) {
        // Given another skewer, compare/coadd stats to this one
        TRACK_MAX(maxcellsize, s.maxcellsize);
        TRACK_MIN(mincellsize, s.mincellsize);
        TRACK_MAX(max_acceleration, s.max_acceleration);
        TRACK_MAX(max_velocity, s.max_velocity);
        TRACK_MIN(min_vrms_on_amax, s.min_vrms_on_amax);
        stddev_cellsize += s.stddev_cellsize;
        sum_square_velocity += s.sum_square_velocity;
    }

    uint64 _search(ilstruct *list, uint64 len, int ygoal) {
        if(len == 0)
            return 0;
        // Set ilskewerstart so that list[ilskewerstart] is the first 
        // element with y >= ygoal.
        // This is done with bisection.
        uint64 low = 0, high = len-1, mid;
        if (list[low].celly()>=ygoal) return low;
        if (list[high].celly()<ygoal) return high+1;
        // Otherwise, we search.
        // list[low] always has y<ygoal.
        // list[high] always has y>=ygoal.
        while (high-low>1) {
            mid = (low+high)/2;
            if (list[mid].celly()>=ygoal) high=mid; else low=mid;
        }
        // We end when the two are only 1 apart.
        assertf(list[low].celly()<list[high].celly()
                && list[low].celly()<ygoal && list[high].celly()>=ygoal, 
             "Error in IL index search; insert list was not sorted!");
        return high;
    }


    void search(ilstruct *list, uint64 len, int ygoal) {
        ilskewerstart = _search(list, len, ygoal);
    }

    uint64 search_forward(ilstruct *list, int zgoal) {
        uint64 i;
        for(i = ilskewerstart; i < ilskewerstart + ilskewerlength; i++){
            if (list[i].local_cellz() >= zgoal){
                break;
            }
        }
        return i;
    }

    uint64 search_backward(ilstruct *list, int zgoal) {
        int64_t i;  // using a signed i
        for(i = (int64_t) ilskewerstart + (int64_t) ilskewerlength - 1; i >= (int64_t) ilskewerstart; i--){
            if (list[i].local_cellz() < zgoal){
                break;
            }
        }
        return i + 1;  // one past the end, so end-start is size
    }
};


void FillMergeSlab(int slab, uint64 *nmerge, uint64 *nmerge_with_ghost) {
    // This routine allocates the MergePos, MergeVel, MergeAux, 
    // and MergeCellInfo slabs.  They will be destroyed when written.

    /* With regard to ghosts, all of the slabs built by this routine
     * will hold MERGE_GHOST_RADIUS ghost columns.  This does not require
     * special treatment, except that the cellinfos will record an extra count,
     * so we know the particle offsets in slabs without ghosts
     * (probably just for AccSlab).
     */

    int cpd = P.cpd;

    // Sort the insert list
    uint64 ilslablength;

    STDLOG(3,"Insert list contains a total of {:d} particles.\n", IL->length);
    ilstruct *ILnew = IL->PartitionAndSort(slab,&ilslablength);
    STDLOG(1,"Insert list contains {:d} new particles for slab {:d}; {:d} remaining\n", ilslablength, slab, IL->length);

    FinishCellIndex.Start();

    SB->AllocateArena(MergeCellInfoSlab, slab);
    SB->AllocateArena(InsertCellInfoSlab, slab);   // Will delete at bottom

    // Make the SkewerIndex objects
    SkewerIndex *skewer = new SkewerIndex[cpd];
    ilstruct *ilhead = ILnew;     // The start of the slab's IL 
    uint64 inslab = ilslablength;    // Total particles to merge in the slab, primary + ghost
    uint64 inslab_no_ghost = 0;  // Just primary

    // Search the insert list to find the start of each skewer
    skewer[0].ilskewerstart = 0;    // First one is obvious!
    // Use static scheduling to get some cache re-use
    #pragma omp parallel for schedule(static)
    for (int y=1; y < cpd; y++) 
        skewer[y].search(ilhead, ilslablength, y);

    // Now find the bounds of each skewer's primary region.
    // The bounds will be near the skewer ends, so we'll search from the front and back.
    #pragma omp parallel for schedule(static)
    for (int y=0; y < cpd; y++) {
        // first establish the skewer length
        skewer[y].ilskewerlength = (y<cpd-1) ?
                (skewer[y+1].ilskewerstart - skewer[y].ilskewerstart)
                : (ilslablength - skewer[y].ilskewerstart);

        // now search from the front and back
        // zgoal is relative to node_z_start - MERGE_GHOST_RADIUS
        skewer[y].ilprimarystart = skewer[y].search_forward(ilhead, MERGE_GHOST_RADIUS);
        uint64 ilprimaryend = skewer[y].search_backward(ilhead, MERGE_GHOST_RADIUS + node_z_size);
        skewer[y].ilprimarylength = ilprimaryend - skewer[y].ilprimarystart;
        
        assertf(skewer[y].ilprimarylength <= skewer[y].ilskewerlength,
            "Primary length ({:d}) must be <= primary + ghost length ({:d})\n",
            skewer[y].ilprimarylength, skewer[y].ilskewerlength);
    }

    // Count the particles from CellInfo in each skewer, loading activelength
    // Use static scheduling because this is fixed work
    // TODO: In principle, this could have been computed in Drift,
    //       thereby avoiding an extra touch of this moderately sized array.
    #pragma omp parallel for schedule(static)
    for (int y=0; y < cpd; y++) {
        skewer[y].activelength = 0;
        // We only consider primary particles from CellInfo; ghosts always come fresh from the IL.
        // And more than that, we know we completely emptied the MERGE_GHOST_RADIUS boundary cells.
        for (int z=node_z_start + MERGE_GHOST_RADIUS; z < node_z_start + node_z_size - MERGE_GHOST_RADIUS; z++) {
            cellinfo *ci =  CP->CellInfo(slab,y,z);
            skewer[y].activelength += ci->active;
        }
    }

    // Cumulate the SkewerIndex's, filling in start and length variables
    // Serial loop
    for (int y=0; y < cpd; y++) {
        // inslab already has all IL primary + ghost from ilslablength.
        // inslab_no_ghost must sum up the primary count skewer-by-skewer
        inslab += skewer[y].activelength;
        inslab_no_ghost += skewer[y].activelength + skewer[y].ilprimarylength;
        
        skewer[y].mergelength = 
                skewer[y].activelength + skewer[y].ilskewerlength;
        
        if(y == 0){
            skewer[y].mergestart = 0;
            skewer[y].mergestart_no_ghost = 0;
        } else {
            skewer[y].mergestart = 
                skewer[y-1].mergestart + skewer[y-1].mergelength;

            uint64 last_mergelength_no_ghost =
                skewer[y-1].activelength + skewer[y-1].ilprimarylength;
            skewer[y].mergestart_no_ghost =
                skewer[y-1].mergestart_no_ghost + last_mergelength_no_ghost;
        }
    }

    // Build the InsertCellInfo and MergeCellInfo indexing
    // Use NUMA_FOR, not because we expect load imbalancing, but to get the merge slabs on the right NUMA nodes
    NUMA_FOR(y,0,cpd, NO_CLAUSE, FALLBACK_STATIC){
        // For this skewer, set up pointers and counters
        ilstruct *ilread = ilhead + skewer[y].ilskewerstart;
            // This pointer will walk along the insert list
        ilstruct *ilend = ilread + skewer[y].ilskewerlength;
            // The end of the insert list for this skewer
        uint64 il_index = skewer[y].ilskewerstart;
        uint64 mci_index = skewer[y].mergestart;
        uint64 mci_index_noghost = skewer[y].mergestart_no_ghost;

        if(ilread < ilend){
            assertf(ilread->local_cellz() < node_z_size + 2*MERGE_GHOST_RADIUS,
                "First particle in insert list skewer falls outside z domain. Invalid ghost?\n"
                );
            assertf((ilend-1)->local_cellz() < node_z_size + 2*MERGE_GHOST_RADIUS,
                "Last particle in insert list skewer falls outside z domain (global z {:d}). Invalid ghost?\n",
                (ilend-1)->global_cellz()
                );
        }

        for(int z = node_z_start - MERGE_GHOST_RADIUS; z < node_z_start + node_z_size + MERGE_GHOST_RADIUS; z++) {
            cellinfo *ici = CP->InsertCellInfo(slab,y,z);
            cellinfo *mci = CP->MergeCellInfo(slab,y,z);

            // The IL will bring in ghosts, but we'll only use the old CellInfos for primaries
            int in_primary = (z >= node_z_start) && (z < node_z_start + node_z_size);

            cellinfo *ci = NULL;
            if(in_primary){
                // Use the real cellinfo
                ci = CP->CellInfo(slab,y,z);
            }

            // Now count the particles on the insert list
            // Note that this requires the insert list sorted to be in y,z order
            int zlocal = z - (node_z_start - MERGE_GHOST_RADIUS);
            ici->startindex_with_ghost = il_index;
            ici->count=0;
            while( ilread < ilend &&
                   ilread->celly() == y && 
                   ilread->local_cellz() == zlocal ) {
                ici->count++;
                il_index++; ilread++;    // Move along the insert list
            }

            // Setup the MergeCellInfo index
            mci->startindex = mci_index_noghost;
            mci->startindex_with_ghost = mci_index;
            if(in_primary)
                assertf(mci->startindex <= mci->startindex_with_ghost,
                    "mci->startindex ({:d}) must be <= mci->startindex_with_ghost ({:d}); cell ({:d},{:d},{:d})\n",
                    mci->startindex, mci->startindex_with_ghost, slab, y, z
                    );
            mci->count = ici->count + (in_primary ? ci->active : 0);
            mci->active = mci->count;
            mci_index += mci->count;
            mci_index_noghost += in_primary ? mci->count : 0;

            // Copy the stats, just to have
            if(in_primary){
                mci->mean_square_velocity = ci->mean_square_velocity;
                mci->max_component_velocity  = ci->max_component_velocity ;
                mci->max_component_acceleration = ci->max_component_acceleration;
                skewer[y].stats(ci, mci);    // Gather the skewer stats
            }
        }
    }
    NUMA_FOR_END;

    STDLOG(2,"Merge slab {:d} will have {:d} old primary particles and {:d} new primaries and {:d} new ghosts\n",
        slab, inslab - ilslablength, inslab_no_ghost - (inslab - ilslablength), inslab - inslab_no_ghost);

    // Accumulate the stats for the full slab
    // Scalar loop
    for (int y=1; y < cpd; y++) skewer[0].addstats(skewer[y]);
    // Now skewer[0] has the slab-based statistics.
    // Can refer to these as skewer->variable

    // Write out the stats
    STDLOG(2,"Cells in slab {:d} range from {:d} to {:d} particles\n", slab, skewer->mincellsize, skewer->maxcellsize);
    TRACK_MAX(WriteState.MaxCellSize, skewer->maxcellsize);
    TRACK_MIN(WriteState.MinCellSize, skewer->mincellsize);

    WriteState.StdDevCellSize += skewer->stddev_cellsize*CP->invcpd3;

    TRACK_MAX(WriteState.MaxVelocity, skewer->max_velocity);
    TRACK_MAX(WriteState.MaxAcceleration, skewer->max_acceleration);
    TRACK_MIN(WriteState.MinVrmsOnAmax, skewer->min_vrms_on_amax);

    STDLOG(2, "Slab {:d}: Max v_j {:f}, Max a_j {:f}, Min <|v|>/amax {:f}\n", 
    	slab, skewer->max_velocity, skewer->max_acceleration, skewer->min_vrms_on_amax);

    WriteState.RMS_Velocity += skewer->sum_square_velocity;

    delete[] skewer;
    FinishCellIndex.Stop();


    // All done with the Cell indexing, now we're ready to copy the particles!
    FinishMerge.Start();

    // Allocate the Merging Slabs to hold 'inslab' number of particles
    SS->set(slab, inslab_no_ghost, inslab);
    SB->AllocateSpecificSize(MergePosSlab, slab, inslab*sizeof(posstruct));
    SB->AllocateSpecificSize(MergeVelSlab, slab, inslab*sizeof(velstruct));
    SB->AllocateSpecificSize(MergeAuxSlab, slab, inslab*sizeof(auxstruct));
    STDLOG(2,"Allocating Merge Slabs to contain {:d} particles\n", inslab);

    size_t nwritten = 0;
    NUMA_FOR(y,0,cpd, reduction(+:nwritten), FALLBACK_DYNAMIC){
        for(int z = node_z_start - MERGE_GHOST_RADIUS; z < node_z_start + node_z_size + MERGE_GHOST_RADIUS; z++) {
            Cell mc = CP->GetMergeCell(slab, y, z);
            cellinfo *ici = CP->InsertCellInfo(slab,y,z);
            int insert_count = ici->count;
            
            int cell_written = 0;
            if(z >= node_z_start + MERGE_GHOST_RADIUS && z < node_z_start + node_z_size - MERGE_GHOST_RADIUS) {
                // Active particles only come from the primary zone; ghosts always come fresh from the IL.
                // And more than that, we emptied MERGE_GHOST_RADIUS boundary cells in the drift.

                Cell c = CP->GetCell(slab, y, z);

                // Copy the active particles from the old cell
                if(c.active() > 0){  // memcpy NULL is undefined, even with count 0
                    memcpy(mc.pos, c.pos, c.active()*sizeof(posstruct));
                    memcpy(mc.vel, c.vel, c.active()*sizeof(velstruct));
                    memcpy(mc.aux, c.aux, c.active()*sizeof(auxstruct));
                }
                cell_written += c.active();
            }

            // Now copy the particles from the insert list
            ilstruct *ilpart = ilhead+ici->startindex_with_ghost;
            for (int j = 0; j < insert_count; j++) {
                    mc.pos[cell_written+j] = ilpart[j].pos;
                    mc.vel[cell_written+j] = ilpart[j].vel;
                    mc.aux[cell_written+j] = ilpart[j].aux;
            }
            cell_written += insert_count;

            assertf(cell_written == mc.count(),
                "Predicted merge cell size doesn't match insert list plus old particles\n");
            nwritten += cell_written;

            // We should reset the L0/L1 bits, so they don't affect
            // future outputs on steps that don't involve group finding.
            if (GFC!=NULL) {  // Group finding happened
                for (int j = 0; j < mc.count(); j++) 
                    mc.aux[j].reset_L01_bits();
            }
			
            // In the vast majority of cases, the PIDs are numbered [0,NP).
            // One could check this condition here; failure could indicate
            // PID corruption
            if(P.MaxPID >= 0){
    			for (int j = 0; j < mc.count(); j++) {
    				assertf(mc.aux[j].unpackpid() < static_cast<size_t>(P.MaxPID),
                        "PID {:d} bigger than MaxPID {:d}\n",
                        mc.aux[j].unpackpid(),
                        P.MaxPID);
    			}
            }
        }
    }
    NUMA_FOR_END;

    assertf(nwritten == inslab, "Wrote {:d} merge particles, expected {:d}. Indexing failure/invalid ghosts?\n",
        nwritten, inslab);

    // Delete the particles from the insert list.
    free(ILnew);   // Need to free this space!
    // IL->ShrinkIL(IL->length - ilslablength);
    STDLOG(2,"After merge, insert list contains a total of {:d} particles.\n", IL->length);
    SB->DeAllocate(InsertCellInfoSlab, slab);
    FinishMerge.Stop();
    
    *nmerge = inslab_no_ghost;
    *nmerge_with_ghost = inslab;
}
