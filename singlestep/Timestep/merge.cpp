/* merge.cpp

This routine combines the particles surviving in a slab (post-drift)
with the new particles arriving from the insert list.  We do this
by inspecting the insert list to build a tabulation of the number of
new particles per cell, then allocating a new arena and copying 
both old and new into the appropriate locations.

This is one of the few places that linear buffer arenas are 
created or destroyed outside of the timestep.cpp file.  Otherwise
we try to keep that behavior explicit in timestep.

We also aggregate the cell-based statistics into global statistics
in this function, placing the results in WriteState so that it can
be output (recall that the statistics are from the ReadState time,
however!).

*/

// Here are some shorthand functions.  
// Beware that these will evaluate 'newvalue' twice.
#define TRACK_MAX(maxsofar,newvalue) if ((newvalue)>maxsofar) maxsofar=newvalue
#define TRACK_MIN(minsofar,newvalue) if ((newvalue)<minsofar) minsofar=newvalue

uint64 FillMergeSlab(int slab) {
    // This routine allocates the MergePos, MergeVel, MergeAux, 
    // and MergeCellInfo slabs.  They will be destroyed when written.

    // Sort the insert list
    uint64 head, ilslablength;
    IL->PartitionAndSort(slab,&head,&ilslablength);
    STDLOG(1,"Insert list contains a total of %d particles.\n", IL->length);
    STDLOG(1,"Insert list contains %d new particles for slab %d\n", ilslablength, slab);

    FinishCellIndex.Start();

    LBW->AllocateArena(MergeCellInfoSlab, slab);
    LBW->AllocateArena(InsertCellInfoSlab, slab);   // Will delete at bottom
    
    #pragma omp parallel for schedule(static)
    for (int y=0; y < P.cpd; y++)
        for (int z=0; z < P.cpd; z++){
            PP->MergeCellInfo(slab,y,z)->makenull();
            PP->InsertCellInfo(slab,y,z)->makenull();
        }

    // Build the InsertCellInfo and MergeCellInfo indexing
    ilstruct *ilread = IL->il + head;
        // This pointer will walk along the insert list
    ilstruct *ilend = ilread + ilslablength;
	// The end of the insert list
    uint64 il_index = 0;		// We will index from IL->il+head, not IL->il
    uint64 mci_index = 0;
    uint64 inslab = 0;

    int maxcellsize = 0, mincellsize = 1e9;
    double mean_cellsize = P.np*PP->invcpd3;
    double stddev_cellsize = 0.0;
    double sum_square_velocity = 0.0;
    FLOAT max_velocity = 0;
    FLOAT max_acceleration = 0;
    FLOAT min_vrms_on_amax = 1e30;

    int cpd = P.cpd;
    for(int y=0;y<cpd;y++)
        for(int z=0;z<cpd;z++) {
            cellinfo *ici = PP->InsertCellInfo(slab,y,z);
            cellinfo *mci = PP->MergeCellInfo(slab,y,z);
            cellinfo *ci =  PP->CellInfo(slab,y,z);

            // Now count the particles on the insert list
            // Note that this requires the insert list sorted to be in y,z order
            ici->startindex = il_index;
            ici->count=0;
            while( ilread<ilend && ((ilread->xyz.x == slab) && 
                    (ilread->xyz.y == y   ) && 
                    (ilread->xyz.z == z   )) ) {
                ici->count++;
                il_index++; ilread++;    // Move along the insert list
            }
            // STDLOG(1,"Counted %d particles in cell %d, %d, %d\n", 
            // ici->count, slab, y, z);

            // Now count the particles in the original cell
            inslab += ci->active;

            // Setup the MergeCellInfo index
            mci->startindex = mci_index;
            mci->count = ci->active + ici->count;
            mci->active = mci->count;
            mci_index += mci->count;

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
            stddev_cellsize += res*res;		// Fractional overdensity
            sum_square_velocity += (double) ci->mean_square_velocity*(ci->count);
    }
    STDLOG(0,"Slab %d contains %d old particles and %d new particles\n", slab, inslab, il_index);

    STDLOG(1,"Cell in slab %d range from %d to %d particles\n", slab, mincellsize, maxcellsize);
    TRACK_MAX(WriteState.MaxCellSize, maxcellsize);
    TRACK_MIN(WriteState.MinCellSize, mincellsize);

    WriteState.StdDevCellSize += stddev_cellsize*PP->invcpd3;

    STDLOG(1,"Maximum v_j in slab %d is %f.\n", slab, max_velocity);
    TRACK_MAX(WriteState.MaxVelocity, max_velocity);

    STDLOG(1,"Maximum a_j in slab %d is %f.\n", slab, max_acceleration);
    TRACK_MAX(WriteState.MaxAcceleration, max_acceleration);

    STDLOG(1,"Minimum <|v|>/amax in slab %d is %f.\n", slab, min_vrms_on_amax);
    TRACK_MIN(WriteState.MinVrmsOnAmax, min_vrms_on_amax);

    WriteState.RMS_Velocity += sum_square_velocity;

    FinishCellIndex.Stop();
    FinishMerge.Start();

    assertf(il_index == ilslablength, 
    	"Failed to count all particles on insert list\n");   
    assertf(mci_index == inslab+ilslablength,
    	"Didn't account for all particles in the merged list\n");
    	// Did we get the right total for the merged list?
    inslab = mci_index;   // This should be the total number of merged particles


    // Allocate the Merging Slabs to hold 'inslab' number of particles
    Slab->set(slab, inslab);
    LBW->AllocateSpecificSize(MergePosSlab, slab, inslab*sizeof(posstruct));
    LBW->AllocateSpecificSize(MergeVelSlab, slab, inslab*sizeof(velstruct));
    LBW->AllocateSpecificSize(MergeAuxSlab, slab, inslab*sizeof(auxstruct));
    STDLOG(1,"Allocating Merge Slabs to contain %d particles\n", inslab);

    ilstruct *ilhead = IL->il+head;

    #pragma omp parallel for schedule(static)
    for(int y=0;y<cpd;y++)
        for(int z=0;z<cpd;z++) {
            Cell c  = PP->GetCell(slab, y, z);
            Cell mc = PP->GetMergeCell(slab, y, z);
            cellinfo *ici = PP->InsertCellInfo(slab,y,z);
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
            for (int j = 0; j<insert_count; j++) {
                mc.pos[j] = ilpart[j].pos;
                mc.vel[j] = ilpart[j].vel;
                mc.aux[j] = ilpart[j].aux;
            }

            // Copy the stats, just to have
            mc.ci->mean_square_velocity = c.ci->mean_square_velocity;
            mc.ci->max_component_velocity  = c.ci->max_component_velocity ;
            mc.ci->max_component_acceleration = c.ci->max_component_acceleration;
        }

    // Delete the particles from the insert list.
    IL->ResetILlength(IL->length - ilslablength);
    STDLOG(1,"After merge, insert list contains a total of %d particles.\n", IL->length);
    LBW->DeAllocate(InsertCellInfoSlab, slab);
    FinishMerge.Stop();
    return inslab;
}
