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
    // so that we can multi-thread easily
    uint64 ilskewerlength;    // How many items in this skewer
    uint64 ilskewerstart;     // Where the indexing in the full list starts

    uint64 activelength;    // How many active original particles in this skewer
    uint64 activestart;     // Where the skewer starts

    uint64 mergelength;     // How many particles in the Merged skewer
    uint64 mergestart;      // Where the Merged skewer will start

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
        if (list[low].xyz.y>=ygoal) return low;
        if (list[high].xyz.y<ygoal) return high+1;
        // Otherwise, we search.
        // list[low] always has y<ygoal.
        // list[high] always has y>=ygoal.
        while (high-low>1) {
            mid = (low+high)/2;
            if (list[mid].xyz.y>=ygoal) high=mid; else low=mid;
        }
        // We end when the two are only 1 apart.
        assertf(list[low].xyz.y<list[high].xyz.y
                && list[low].xyz.y<ygoal && list[high].xyz.y>=ygoal, 
             "Error in IL index search; insert list was not sorted!");
        return high;
    }
    void search(ilstruct *list, uint64 len, int ygoal) {
        ilskewerstart = _search(list, len, ygoal);
    }
};


uint64 FillMergeSlab(int slab) {
    // This routine allocates the MergePos, MergeVel, MergeAux, 
    // and MergeCellInfo slabs.  They will be destroyed when written.
    int cpd = P.cpd;

    // Sort the insert list
    uint64 ilslablength;

    STDLOG(3,"Insert list contains a total of %d particles.\n", IL->length);
    ilstruct *ILnew = IL->PartitionAndSort(slab,&ilslablength);
    STDLOG(1,"Insert list contains %d new particles for slab %d; %d remaining\n", ilslablength, slab, IL->length);

    FinishCellIndex.Start();

    SB->AllocateArena(MergeCellInfoSlab, slab);
    SB->AllocateArena(InsertCellInfoSlab, slab);   // Will delete at bottom

    // Make the SkewerIndex objects
    SkewerIndex *skewer = new SkewerIndex[cpd];
    ilstruct *ilhead = ILnew;     // The start of the slab's IL 
    uint64 inslab = 0;    // Total particles in the slab

    // Search the insert list to find the start of each skewer
    skewer[0].ilskewerstart = 0;    // First one is obvious!
    // Use static scheduling to get some cache re-use
    #pragma omp parallel for schedule(static)
    for (int y=1; y < cpd; y++) 
        skewer[y].search(ilhead, ilslablength, y);
    

    // Count the particles from CellInfo in each skewer, loading activelength
    // Use static scheduling because this is fixed work
    // TODO: In principle, this could have been computed in Drift,
    //       thereby avoiding an extra touch of this moderately sized array.
    #pragma omp parallel for schedule(static)
    for (int y=0; y < cpd; y++) {
        skewer[y].activelength = 0;
        for (int z=0; z < cpd; z++) {
            cellinfo *ci =  CP->CellInfo(slab,y,z);
            skewer[y].activelength += ci->active;
        }
    }

    // Cumulate the SkewerIndex's, filling in start and length variables
    // Scalar loop
    for (int y=0; y < cpd; y++) {
        inslab += skewer[y].activelength;   // Total active particles
        skewer[y].activestart = (y>0)?
                (skewer[y-1].activestart+skewer[y-1].activelength )
                :0;
        skewer[y].ilskewerlength = (y<cpd-1)?
                (skewer[y+1].ilskewerstart-skewer[y].ilskewerstart)
                : (ilslablength-skewer[y].ilskewerstart);
        skewer[y].mergelength = 
                skewer[y].activelength + skewer[y].ilskewerlength;
        skewer[y].mergestart = (y>0)?
                (skewer[y-1].mergestart+skewer[y-1].mergelength )
                :0;
    }

    // Build the InsertCellInfo and MergeCellInfo indexing
    #pragma omp parallel for schedule(static)
    for(int y=0;y<cpd;y++) {
        // For this skewer, set up pointers and counters
        ilstruct *ilread = ilhead + skewer[y].ilskewerstart;
            // This pointer will walk along the insert list
        ilstruct *ilend = ilread + skewer[y].ilskewerlength;
            // The end of the insert list for this skewer
        uint64 il_index = skewer[y].ilskewerstart;                
        uint64 mci_index = skewer[y].mergestart;

        for(int z=0;z<cpd;z++) {
            cellinfo *ici = CP->InsertCellInfo(slab,y,z);
            cellinfo *mci = CP->MergeCellInfo(slab,y,z);
            cellinfo *ci =  CP->CellInfo(slab,y,z);

            // Now count the particles on the insert list
            // Note that this requires the insert list sorted to be in y,z order
            ici->startindex = il_index;
            ici->count=0;
            // TODO: The x test is unneeded, since we 
            // partitioned the slab to the end of the insert list.
            while( ilread<ilend && ((ilread->xyz.x == slab) && 
                    (ilread->xyz.y == y   ) && 
                    (ilread->xyz.z == z   )) ) {
                ici->count++;
                il_index++; ilread++;    // Move along the insert list
            }

            // Setup the MergeCellInfo index
            mci->startindex = mci_index;
            mci->count = ci->active + ici->count;
            mci->active = mci->count;
            mci_index += mci->count;

            // Copy the stats, just to have
            mci->mean_square_velocity = ci->mean_square_velocity;
            mci->max_component_velocity  = ci->max_component_velocity ;
            mci->max_component_acceleration = ci->max_component_acceleration;
            skewer[y].stats(ci, mci);    // Gather the skewer stats
        }
    }
    STDLOG(2,"Slab %d contains %d old particles and %d new particles\n", slab, inslab, ilslablength);

    // Accumulate the stats for the full slab
    // Scalar loop
    for (int y=1; y < cpd; y++) skewer[0].addstats(skewer[y]);
    // Now skewer[0] has the slab-based statistics.
    // Can refer to these as skewer->variable

    // Write out the stats
    STDLOG(2,"Cells in slab %d range from %d to %d particles\n", slab, skewer->mincellsize, skewer->maxcellsize);
    TRACK_MAX(WriteState.MaxCellSize, skewer->maxcellsize);
    TRACK_MIN(WriteState.MinCellSize, skewer->mincellsize);

    WriteState.StdDevCellSize += skewer->stddev_cellsize*CP->invcpd3;

    TRACK_MAX(WriteState.MaxVelocity, skewer->max_velocity);
    TRACK_MAX(WriteState.MaxAcceleration, skewer->max_acceleration);
    TRACK_MIN(WriteState.MinVrmsOnAmax, skewer->min_vrms_on_amax);

    STDLOG(2, "Slab %d: Max v_j %f, Max a_j %f, Min <|v|>/amax %f\n", 
    	slab, skewer->max_velocity, skewer->max_acceleration, skewer->min_vrms_on_amax);

    WriteState.RMS_Velocity += skewer->sum_square_velocity;

    delete[] skewer;
    FinishCellIndex.Stop();


    // All done with the Cell indexing, now we're ready to copy the particles!
    FinishMerge.Start();
    inslab += ilslablength;   
            // This should be the total number of merged particles

    // Allocate the Merging Slabs to hold 'inslab' number of particles
    SS->set(slab, inslab);
    SB->AllocateSpecificSize(MergePosSlab, slab, inslab*sizeof(posstruct));
    SB->AllocateSpecificSize(MergeVelSlab, slab, inslab*sizeof(velstruct));
    SB->AllocateSpecificSize(MergeAuxSlab, slab, inslab*sizeof(auxstruct));
    STDLOG(2,"Allocating Merge Slabs to contain %d particles\n", inslab);

    #pragma omp parallel for schedule(static)
    for(int y=0;y<cpd;y++){
        for(int z=0;z<cpd;z++) {
            Cell c;
            c = CP->GetCell(slab, y, z);
            Cell mc = CP->GetMergeCell(slab, y, z);
            cellinfo *ici = CP->InsertCellInfo(slab,y,z);
            int insert_count = ici->count;
            
            assertf(insert_count + c.active() == mc.count(),
                    "Predicted merge cell size doesn't match insert list plus old particles\n");


            // Copy the active particles from the old cell
            memcpy(mc.pos, c.pos, c.active()*sizeof(posstruct));
            memcpy(mc.vel, c.vel, c.active()*sizeof(velstruct));
            memcpy(mc.aux, c.aux, c.active()*sizeof(auxstruct));
            mc.pos += c.active();
            mc.vel += c.active();
            mc.aux += c.active();

            // Now copy the particles from the insert list
            ilstruct *ilpart = ilhead+ici->startindex;
            #pragma ivdep
            for (int j = 0; j < insert_count; j++) {
                    mc.pos[j] = ilpart[j].pos;
                    mc.vel[j] = ilpart[j].vel;
                    mc.aux[j] = ilpart[j].aux;
            }

            // We should reset the L0/L1 bits, so they don't affect
            // future outputs on steps that don't involve group finding.
            if (GFC!=NULL) {  // Group finding happened
                auxstruct *origaux = mc.aux-c.active();
                for (int j = 0; j < mc.count(); j++) 
                    origaux[j].reset_L01_bits();
            }
			
            // In the vast majority of cases, the PIDs are numbered [0,NP).
            // One could check this condition here; failure could indicate
            // PID corruption
            if(P.MaxPID >= 0){
    			mc = CP->GetMergeCell(slab, y, z);
    			for (int j = 0; j < mc.count(); j++) {
    				assertf(mc.aux[j].pid() < P.MaxPID, "PID %d bigger than MaxPID %d\n", mc.aux[j].pid(), P.MaxPID);
    			}
            }
            
            //DoCellMultipoles(slab, y, z);
        }
        //DoZRowMultipoles(slab, y);
    }

    /*#pragma omp parallel for schedule(static)
    for(int z=0;z<cpd;z++)
        DoYRowMultipoles(slab, z);

    // Do the final thread reductions
    FinishSlabMultipoles(slab);*/

    // Delete the particles from the insert list.
    free(ILnew);   // Need to free this space!
    // IL->ShrinkIL(IL->length - ilslablength);
    STDLOG(2,"After merge, insert list contains a total of %d particles.\n", IL->length);
    SB->DeAllocate(InsertCellInfoSlab, slab);
    FinishMerge.Stop();
    return inslab;
}
