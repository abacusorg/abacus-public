// This is the top-level code for group-finding.
// We include this file in the program.


// TODO: Are the SlabAccum maxsize set sensibly?

#include "fof_sublist.cpp"
	// Code to do FOF on a set (e.g., a cell)

#include "slab_accum.cpp"
	// Code to establish templated slab-based storage of flexible size 
	// that is cell indexed and multi-threaded by pencil

#include "grouplink.cpp"
	// Code to store the list of links between CellGroups.
	// This is GCC->GLL

#include "cellgroup.cpp"
	// Code to define and create CellGroups.

#include "multiplicity_stats.cpp"



class GroupFindingControl {
    // This is the control class for all group finding.
  public:
    int cpd;		
    FOFloat linking_length;    // In code units
    FOFloat boundary;	// The distance from the origin that indicates that something
    		// is within linking_length of the edge
    int GroupRadius;    // Cell radius for group finding
    Float invcpd;
    uint64 np; // the number of particles in the simulation; used to compute buffers

    // Parameters for finding subgroups
    FOFloat linking_length_level1;
    FOFloat linking_length_level2;
    int minhalosize;	// The minimum size for a level 1 halo to be outputted

    uint64 pPtot, fPtot, fGtot, CGtot, GGtot, Ltot;
    int largest_GG;

    MultiplicityStats L0stats, L1stats, L2stats;

    SlabAccum<CellGroup> *cellgroups;   // Allocated [0,cpd), one for each slab
    int *cellgroups_status;     // Allocated [0,cpd) for each slab. 
    	// 0 means CellGroups not found.
    	// 1 means CellGroups found but slab not closed.
    	// 2 means slab closed (so CellGroups may be deleted).
    GroupLinkList *GLL;

    STimer CellGroupTime, CreateFaceTime, FindLinkTime, 
    	SortLinks, IndexLinks, FindGlobalGroupTime, IndexGroups, 
	GatherGroups, ScatterGroups, ProcessLevel1;
    PTimer L1FOF, L2FOF, L1Tot;

    GroupFindingControl(FOFloat _linking_length, 
    	FOFloat _linking_length_level1, FOFloat _linking_length_level2,
    int _cpd, FOFloat _invcpd, int _GroupRadius, int _minhalosize, uint64 _np) {
	cpd = _cpd; 
	linking_length = _linking_length;
	linking_length_level1 = _linking_length_level1;
	linking_length_level2 = _linking_length_level2;
	minhalosize = _minhalosize;
	invcpd = _invcpd;
	boundary = (invcpd/2.0-linking_length);
	np = _np;
	GroupRadius = _GroupRadius;
	cellgroups = new SlabAccum<CellGroup>[cpd];
	cellgroups_status = new int[cpd];
	for (int j=0;j<cpd;j++) cellgroups_status[j] = 0;
	GLL = new GroupLinkList(cpd, np/cpd*linking_length/_invcpd*3*15);    
	// This is a MultiAppendList, so the buffer cannot grow. 
	// It will be storing links between cells.  Most particles do 
	// not generate such a link; on average 3*linking_length/cellsize do. 
	// But we will have up to 10 slabs worth in here.
	// So we need to be a bit generous.
        pPtot = fPtot = fGtot = CGtot = GGtot = Ltot = 0;
	largest_GG = 0;
	return;
    }

    ~GroupFindingControl() { 
        delete[] cellgroups;
        delete[] cellgroups_status;
	delete GLL;
    }
    int WrapSlab(int s) {
	while (s<0) s+=cpd;
	while (s>=cpd) s-=cpd;
	return s;
    }
    integer3 WrapCell(int i, int j, int k) {
        return integer3(WrapSlab(i), WrapSlab(j), WrapSlab(k));
    }
    
    void ConstructCellGroups(int slab);
    void DestroyCellGroups(int slab);


    void report() {
	 // TODO: Eventually these should be STDLOG
	 printf("Found %lld cell groups (including boundary singlets)\n", CGtot);
	 printf("Used %lld pseudoParticles, %lld faceParticles, %lld faceGroups\n",
	     pPtot, fPtot, fGtot);
	 printf("Found %lld links between groups.\n", Ltot);
	 printf("Found %lld global groups\n", GGtot);
	 printf("Longest GroupLink list was %lld, compared to %lld allocation\n", GLL->longest, GLL->maxlist);
	 printf("Largest Global Group has %d particles\n", largest_GG);

	 printf("\nL1 group multiplicity distribution:\n");
	 L1stats.report();

	 float total_time = CellGroupTime.Elapsed()+
			CreateFaceTime.Elapsed()+
			FindLinkTime.Elapsed()+
			SortLinks.Elapsed()+
			IndexLinks.Elapsed()+
			FindGlobalGroupTime.Elapsed()+
			IndexGroups.Elapsed()+
			GatherGroups.Elapsed()+
			ProcessLevel1.Elapsed()+
			ScatterGroups.Elapsed();
	 printf("\nTimings: \n");
	 #define RFORMAT(a) a.Elapsed(), a.Elapsed()/total_time*100.0
	 printf("Finding Cell Groups:     %8.4f sec (%5.2f%%)\n",
			RFORMAT(CellGroupTime));
	 printf("Creating Faces:          %8.4f sec (%5.2f%%)\n",
			RFORMAT(CreateFaceTime));
	 printf("Finding Group Links:     %8.4f sec (%5.2f%%)\n",
			RFORMAT(FindLinkTime));
	 printf("Sort Links:              %8.4f sec (%5.2f%%)\n",
			RFORMAT(SortLinks));
	 printf("Index Links:             %8.4f sec (%5.2f%%)\n",
			RFORMAT(IndexLinks));
	 printf("Find Global Groups:      %8.4f sec (%5.2f%%)\n",
			RFORMAT(FindGlobalGroupTime));
	 printf("Index Global Groups:     %8.4f sec (%5.2f%%)\n",
			RFORMAT(IndexGroups));
	 printf("Gather Group Particles:  %8.4f sec (%5.2f%%)\n",
			RFORMAT(GatherGroups));
	 printf("Level 1 & 2 Processing:  %8.4f sec (%5.2f%%)\n",
			RFORMAT(ProcessLevel1));
	 printf("Level 1 FOF:                  %8.4f sec (%5.2f%%)\n",
			RFORMAT(L1FOF));
	 printf("Level 2 FOF:                  %8.4f sec (%5.2f%%)\n",
			RFORMAT(L2FOF));
	 printf("Level 1 Total:                %8.4f sec (%5.2f%%)\n",
			RFORMAT(L1Tot));
	 printf("Scatter Group Particles: %8.4f sec (%5.2f%%)\n",
			RFORMAT(ScatterGroups));
	 printf("Total Booked Time:       %8.4f sec (%5.2f Mp/sec)\n", total_time, np/total_time*1e-6);
	 #undef RFORMAT
    }
};

// We keep the code for Constructing CellGroups in here, because
// it is where the cellgroups[] are created and destroyed.

void GroupFindingControl::ConstructCellGroups(int slab) {
    // Construct the Cell Groups for this slab
    CellGroupTime.Start();
    slab = WrapSlab(slab);
    cellgroups[slab].setup(cpd, 1024);     // TODO: Pick a better maxsize!
    FOFcell doFOF[omp_get_max_threads()];
    #pragma omp parallel for schedule(static)
    for (int g=0; g<omp_get_max_threads(); g++) 
    	doFOF[g].setup(linking_length, boundary);

    #pragma omp parallel for schedule(dynamic,1)
    for (int j=0; j<cpd; j++) {
	float *aligned;
	int ret = posix_memalign((void **)&(aligned), 64, 8*sizeof(float)); assert(ret==0);
	int g = omp_get_thread_num();
        PencilAccum<CellGroup> *cg = cellgroups[slab].StartPencil(j);
        for (int k=0; k<cpd; k++) {
	    // Find CellGroups in (slab,j,k).  Append results to cg.
	    Cell c = PP->GetCell(slab, j, k);
	    doFOF[g].findgroups(c.pos, c.vel, c.aux, c.acc, c.count());
	    // printf("Cell %d %d %d: Found %d cell groups with %d particles, plus %d boundary singlets\n", slab,j,k,doFOF[g].ngroups, doFOF[g].nmultiplets, doFOF[g].nsinglet_boundary-doFOF[g].nmultiplets);
	    // We need to clear the L0 & L1 bits for this timestep
	    for (int p=0; p<c.count(); p++) c.aux[p].reset_L01_bits();
	    for (int gr=0; gr<doFOF[g].ngroups; gr++) {
		CellGroup tmp(doFOF[g].groups[gr], boundary, aligned);
		// printf("Group %d: %d %d, %f %f %f to %f %f %f, %f %02x\n",
		// 	gr, tmp.start, tmp.size(), 
		// 	aligned[0], aligned[1], aligned[2],
		// 	aligned[4], aligned[5], aligned[6], boundary, tmp.n >> 25);
		cg->append(tmp);
	    }
	    // Also need to look at the singlets!
	    for (int p=doFOF[g].nmultiplets; p<doFOF[g].nsinglet_boundary; p++) {
		CellGroup tmp(p, c.pos[p], boundary);
		// printf("Singl %d: %f %f %f %lld %02x\n",
		   // p, c.pos[p].x, c.pos[p].y, c.pos[p].z, c.aux[p].val, tmp.n >>25);
	        cg->append(tmp);
	    }
	    // for (int p=0; p<doFOF[g].np; p++) {
	        // if (c.aux[p].val == 1516) {
		    // printf("Singl %d: %f %f %f %lld     %d\n",
		       // p, c.pos[p].x, c.pos[p].y, c.pos[p].z, c.aux[p].val,
		       // doFOF[g].nsinglet_boundary);
		// }
	    // }
	    cg->FinishCell();
	}
	cg->FinishPencil();
	free(aligned);
    }
    #pragma omp parallel for schedule(static)
    for (int g=0; g<omp_get_max_threads(); g++) 
    	doFOF[g].destroy();
    uint64 tot = cellgroups[slab].get_slab_size();
    CGtot += tot;
    cellgroups_status[slab] = 1;
    STDLOG(1,"Found %lld cell groups in slab %d\n", tot, slab);
    CellGroupTime.Stop();
    return;
}

void GroupFindingControl::DestroyCellGroups(int slab) {
    slab = WrapSlab(slab);
    cellgroups[slab].destroy();
    cellgroups_status[slab] = 2;
    return;
}

GroupFindingControl *GFC;	
	// We have one global instance of this, initialized in proepi


#include "findgrouplinks.cpp"
	// Code to search between pairs of cells and find the linked groups,
	// which get added to GLL.

#include "globalgroup.cpp"
	// Code to traverse the links and find the GlobalGroups as 
	// sets of CellGroups (stored by their LinkIDs).

#include "halostat.cpp"
	// Code to compute L1 halo properties

// #include "particle_subsample.cpp"

// ===================== Global Groups: Find & Process ===============


class GlobalGroupSlab {
  public:
    SlabAccum<GlobalGroup> globalgroups;
    	// The global group information
    SlabAccum<LinkID> globalgrouplist;
    	// The cell group decompositions of these global groups

    // The following accumulate possible output 
    SlabAccum<HaloStat> L1halos;	// Stats about each L1 halo
    SlabAccum<TaggedPID> TaggedPIDs;	// The tagged PIDs in each L1 halo
    SlabAccum<RVfloat> L1Particles;     // The taggable subset in each L1 halo, pos/vel
    SlabAccum<TaggedPID> L1PIDs;	// The taggable subset in each L1 halo, PID
    
    int slab;    // This slab number
    posstruct *pos;  // Big vectors of all of the pos/vel/aux for these global groups
    velstruct *vel;
    auxstruct *aux;
    accstruct *acc;
    uint64 np;

    int largest_group;

    GlobalGroupSlab(int _slab) {
        pos = NULL; vel = NULL; aux = NULL; acc = NULL; np = 0;
	slab = _slab; largest_group = 0;
    }
    void destroy() {
	// TODO: If we use arenas, change this
        if (pos!=NULL) free(pos); pos = NULL;
        if (vel!=NULL) free(vel); vel = NULL;
        if (aux!=NULL) free(aux); aux = NULL;
        if (acc!=NULL) free(acc); acc = NULL;
	np = 0;
    }
    ~GlobalGroupSlab() { destroy(); }
    void setup(uint64 _np) {
	destroy();
	np = _np;
	// TODO: May eventually prefer these to be arenas.
	int ret;
        ret = posix_memalign((void **)&pos, 4096, sizeof(posstruct)*np); assert(ret==0);
        ret = posix_memalign((void **)&vel, 4096, sizeof(velstruct)*np); assert(ret==0);
        ret = posix_memalign((void **)&aux, 4096, sizeof(auxstruct)*np); assert(ret==0);
        ret = posix_memalign((void **)&acc, 4096, sizeof(accstruct)*np); assert(ret==0);
    }

    void GatherGlobalGroups() {
	// Copy all of the particles into a single list
	GFC->IndexGroups.Start();
	int diam = 2*GFC->GroupRadius+1;
	assert(GFC->cpd>=4*GFC->GroupRadius+1);
	// TODO: This registers the periodic wrap using the cells.
	// However, this will misbehave if CPD is smaller than the group diameter,
	// because the LinkIDs have already been wrapped.
	uint64 *this_pencil = new uint64[GFC->cpd];
	int *largest = new int[GFC->cpd];

	#pragma omp parallel for schedule(static) 
	for (int j=0; j<GFC->cpd; j++) {
	    uint64 local_this_pencil = 0;
	    int local_largest = 0;
	    for (int k=0; k<GFC->cpd; k++)
		for (int n=0; n<globalgroups[j][k].size(); n++) {
		    int size = globalgroups[j][k][n].np;
		    local_this_pencil += size;
		    local_largest = std::max(local_largest, size);
		}
	    this_pencil[j] = local_this_pencil;
	    largest[j] = local_largest;
	}

	// Now for the serial accumulation over the pencils
	uint64 *pstart = new uint64[GFC->cpd];
	uint64 total_particles = 0;
	largest_group = 0;
	for (int j=0; j<GFC->cpd; j++) {
	    pstart[j] = total_particles;    // So we have the starting indices
	    total_particles += this_pencil[j];
	    largest_group = std::max(largest[j], largest_group);
	}
	delete[] largest;
	delete[] this_pencil;

	// Now we can allocate these buffers
	setup(total_particles);
	GFC->IndexGroups.Stop();

	GFC->GatherGroups.Start();
	// Now copy the particles into these structures
	#pragma omp parallel for schedule(static)
	for (int j=0; j<GFC->cpd; j++)
	    for (int k=0; k<GFC->cpd; k++)
		for (int n=0; n<globalgroups[j][k].size(); n++) {
		    // Process globalgroups[j][k][n]
		    // Compute where we'll put the particles, and update this starting point
		    globalgroups[j][k].ptr(n)->start += pstart[j];
		    uint64 start = globalgroups[j][k][n].start;

		    LinkID *cglink = globalgrouplist.pencils[j].data
		    			+globalgroups[j][k][n].cellgroupstart;
		    	// This is where we'll find the CG LinkIDs for this GG
		    integer3 firstcell(slab,j,k);

		    for (int c=0; c<globalgroups[j][k][n].ncellgroups; c++, cglink++) {
		        // Loop over CellGroups
			integer3 cellijk = cglink->cell();
			CellGroup *cg = LinkToCellGroup(*cglink);
			// printf("%d %d %d %d in %d %d %d n=%d\n", 
			    // j, k, n, c, cellijk.x, cellijk.y, cellijk.z, cg->size());
			// CellGroup cg = GFC->cellgroups[cellijk.x][cellijk.y][cellijk.z][cglink->cellgroup()];
			Cell cell = PP->GetCell(cellijk);
			// Copy the particles, including changing positions to 
			// the cell-centered coord of the first cell.  
			// Note periodic wrap.
			memcpy(vel+start, cell.vel+cg->start, sizeof(velstruct)*cg->size());
			memcpy(aux+start, cell.aux+cg->start, sizeof(auxstruct)*cg->size());
			for (int p=0; p<cg->size(); p++) aux[p].set_L0();
				// This particle is in L0
			memcpy(acc+start, cell.acc+cg->start, sizeof(accstruct)*cg->size());
			cellijk -= firstcell;
			if (cellijk.x> diam) cellijk.x-=GFC->cpd;
			if (cellijk.x<-diam) cellijk.x+=GFC->cpd;
			if (cellijk.y> diam) cellijk.y-=GFC->cpd;
			if (cellijk.y<-diam) cellijk.y+=GFC->cpd;
			if (cellijk.z> diam) cellijk.z-=GFC->cpd;
			if (cellijk.z<-diam) cellijk.z+=GFC->cpd;
			posstruct offset = GFC->invcpd*(cellijk);
			// printf("Using offset %f %f %f\n", offset.x, offset.y, offset.z);
			for (int p=0; p<cg->size(); p++) pos[start+p] = offset+cell.pos[cg->start+p];
			start += cg->size();
		    } // End loop over cellgroups in this global group

		    // TODO: Might need to compute the COM for light cone output
		} // End loop over globalgroups in a cell
	// End loop over cells
	delete[] pstart;
	GFC->GatherGroups.Stop();
	return;
    }

    void ScatterGlobalGroupsAux() {
        // Write the information from pos,vel,aux back into the original Slabs
	int diam = 2*GFC->GroupRadius+1;
	GFC->ScatterGroups.Start();

	#pragma omp parallel for schedule(static)
	for (int j=0; j<GFC->cpd; j++)
	    for (int k=0; k<GFC->cpd; k++)
		for (int n=0; n<globalgroups[j][k].size(); n++) {
		    // Process globalgroups[j][k][n]
		    // Recall where the particles start
		    uint64 start = globalgroups[j][k][n].start;

		    LinkID *cglink = globalgrouplist.pencils[j].data
		    			+globalgroups[j][k][n].cellgroupstart;
		    	// This is where we'll find the CG LinkIDs for this GG
		    integer3 firstcell(slab,j,k);
		    for (int c=0; c<globalgroups[j][k][n].ncellgroups; c++, cglink++) {
		        // Loop over CellGroups
			integer3 cellijk = cglink->cell();
			CellGroup *cg = LinkToCellGroup(*cglink);
			Cell cell = PP->GetCell(cellijk);
			// Copy the aux back
			memcpy(cell.aux+cg->start, aux+start, sizeof(auxstruct)*cg->size());
			start += cg->size();
		    } // End loop over cellgroups in this global group
		} // End loop over globalgroups in a cell
	// End loop over cells
	GFC->ScatterGroups.Stop();
	return;
    }

    void ScatterGlobalGroups() {
        // Write the information from pos,vel,acc back into the original Slabs
	// Note that aux is handled separately!
	int diam = 2*GFC->GroupRadius+1;
	GFC->ScatterGroups.Start();

	#pragma omp parallel for schedule(static)
	for (int j=0; j<GFC->cpd; j++)
	    for (int k=0; k<GFC->cpd; k++)
		for (int n=0; n<globalgroups[j][k].size(); n++) {
		    // Process globalgroups[j][k][n]
		    // Recall where the particles start
		    uint64 start = globalgroups[j][k][n].start;

		    LinkID *cglink = globalgrouplist.pencils[j].data
		    			+globalgroups[j][k][n].cellgroupstart;
		    	// This is where we'll find the CG LinkIDs for this GG
		    integer3 firstcell(slab,j,k);
		    for (int c=0; c<globalgroups[j][k][n].ncellgroups; c++, cglink++) {
		        // Loop over CellGroups
			integer3 cellijk = cglink->cell();
			CellGroup *cg = LinkToCellGroup(*cglink);
			Cell cell = PP->GetCell(cellijk);
			// Copy the particles, including changing positions to 
			// the cell-centered coord of the first cell.  
			// Note periodic wrap.
			memcpy(cell.vel+cg->start, vel+start, sizeof(velstruct)*cg->size());
			memcpy(cell.acc+cg->start, acc+start, sizeof(accstruct)*cg->size());
			cellijk -= firstcell;
			if (cellijk.x> diam) cellijk.x-=GFC->cpd;
			if (cellijk.x<-diam) cellijk.x+=GFC->cpd;
			if (cellijk.y> diam) cellijk.y-=GFC->cpd;
			if (cellijk.y<-diam) cellijk.y+=GFC->cpd;
			if (cellijk.z> diam) cellijk.z-=GFC->cpd;
			if (cellijk.z<-diam) cellijk.z+=GFC->cpd;
			posstruct offset = GFC->invcpd*(cellijk);
			// printf("Using offset %f %f %f\n", offset.x, offset.y, offset.z);
			for (int p=0; p<cg->size(); p++) 
			     cell.pos[cg->start+p] = pos[start+p] - offset;
			start += cg->size();
		    } // End loop over cellgroups in this global group
		} // End loop over globalgroups in a cell
	// End loop over cells
	GFC->ScatterGroups.Stop();
	return;
    }


    void FindSubGroups() {
	// Process each group, looking for L1 and L2 subgroups.
	GFC->ProcessLevel1.Start();
	L1halos.setup(GFC->cpd, 1024);    // TODO: Need better start value
	TaggedPIDs.setup(GFC->cpd, 1024);    // TODO: Need better start value
	L1Particles.setup(GFC->cpd, 1024);    // TODO: Need better start value
	L1PIDs.setup(GFC->cpd, 1024);    // TODO: Need better start value
	FOFcell FOFlevel1[omp_get_max_threads()], FOFlevel2[omp_get_max_threads()];
	posstruct **L1pos = new posstruct *[omp_get_max_threads()];
	velstruct **L1vel = new velstruct *[omp_get_max_threads()];
	auxstruct **L1aux = new auxstruct *[omp_get_max_threads()];
	accstruct **L1acc = new accstruct *[omp_get_max_threads()];
	MultiplicityStats L1stats[omp_get_max_threads()];

	#pragma omp parallel for schedule(static)
	for (int g=0; g<omp_get_max_threads(); g++) {
	    FOFlevel1[g].setup(GFC->linking_length_level1, 1e10);
	    FOFlevel2[g].setup(GFC->linking_length_level2, 1e10);
	    L1pos[g] = (posstruct *)malloc(sizeof(posstruct)*largest_group);
	    L1vel[g] = (velstruct *)malloc(sizeof(velstruct)*largest_group);
	    L1aux[g] = (auxstruct *)malloc(sizeof(auxstruct)*largest_group);
	    L1acc[g] = (accstruct *)malloc(sizeof(accstruct)*largest_group);
	}
	
	#pragma omp parallel for schedule(dynamic,1)
	for (int j=0; j<GFC->cpd; j++) {
	    GFC->L1Tot.Start();
	    int g = omp_get_thread_num();
	    PencilAccum<HaloStat> *pL1halos = L1halos.StartPencil(j);
	    PencilAccum<TaggedPID> *pTaggedPIDs = TaggedPIDs.StartPencil(j);
	    PencilAccum<RVfloat> *pL1Particles = L1Particles.StartPencil(j);
	    PencilAccum<TaggedPID> *pL1PIDs = L1PIDs.StartPencil(j);
	    for (int k=0; k<GFC->cpd; k++) {
		// uint64 groupid = ((slab*GFC->cpd+j)*GFC->cpd+k)*4096;
		uint64 groupid = (((uint64)slab*10000+(uint64)j)*10000+(uint64)k)*1000;
			// A basic label for this group
		for (int n=0; n<globalgroups[j][k].size(); n++) {
		    if (globalgroups[j][k][n].np<GFC->minhalosize) continue;
		    // We have a large-enough global group to process
		    posstruct *grouppos = pos+globalgroups[j][k][n].start;
		    velstruct *groupvel = vel+globalgroups[j][k][n].start;
		    auxstruct *groupaux = aux+globalgroups[j][k][n].start;
		    accstruct *groupacc = acc+globalgroups[j][k][n].start;
		    int groupn = globalgroups[j][k][n].np;
		    GFC->L1FOF.Start();
		    FOFlevel1[g].findgroups(grouppos, NULL, NULL, NULL, groupn);
		    GFC->L1FOF.Stop();
		    // Now we've found the L1 groups
		    for (int a=0; a<FOFlevel1[g].ngroups; a++) {
			int size = FOFlevel1[g].groups[a].n;
		        if (size<GFC->minhalosize) continue;
			L1stats[g].push(size);
			// The group is big enough.
			FOFparticle *start = FOFlevel1[g].p+FOFlevel1[g].groups[a].start;
			// Particle indices are in start[0,size).index()
			// which deref the grouppos, groupvel, groupaux list.

			// We now need to find the L2 subgroups.
			// Make a temporary list of the particles in the L1 group
			for (int b=0; b<size; b++) {
			    L1pos[g][b] = grouppos[start[b].index()];
			    L1vel[g][b] = groupvel[start[b].index()];
			    groupaux[start[b].index()].set_L1();
			    L1aux[g][b] = groupaux[start[b].index()];
			    L1acc[g][b] = groupacc[start[b].index()];
			}
			GFC->L2FOF.Start();
			FOFlevel2[g].findgroups(L1pos[g], NULL, NULL, NULL, size);
			std::sort(FOFlevel2[g].groups, FOFlevel2[g].groups+FOFlevel2[g].ngroups);
			    // Groups now in descending order of multiplicity
			GFC->L2FOF.Stop();

			// Merger trees require tagging the taggable particles 
			// of the biggest L2 group in the original aux array.  
			// This can be done:
			// The L2 index() gives the position in the L1 array,
			// and that index() gets back to aux.
			FOFparticle *L2start = FOFlevel2[g].p + FOFlevel2[g].groups[0].start;
			for (int p=0; p<FOFlevel2[g].groups[0].n; p++) {
			    if (groupaux[start[L2start[p].index()].index()].is_taggable())
				groupaux[start[L2start[p].index()].index()].set_tagged();
			}

			uint64 taggedstart = pTaggedPIDs->get_pencil_size();
			uint64 npstart = pL1Particles->get_pencil_size();

			// Output the Tagged PIDs
			for (int b=0; b<size; b++)
			    if (groupaux[start[b].index()].is_tagged()) 
			    	pTaggedPIDs->append(TaggedPID(groupaux[start[b].index()].pid()));

			// Output the Taggable Particles
			posstruct offset = PP->CellCenter(slab, j, k);
			for (int b=0; b<size; b++)
			    if (groupaux[start[b].index()].is_taggable()) {
				posstruct r = WrapPosition(grouppos[start[b].index()]+offset);
				velstruct v = groupvel[start[b].index()];
			    	pL1Particles->append(RVfloat(r.x, r.y, r.z, v.x, v.y, v.z));
			    	pL1PIDs->append(TaggedPID(groupaux[start[b].index()].pid()));
			    }

			HaloStat h = ComputeStats(size, L1pos[g], L1vel[g], L1aux[g], FOFlevel2[g], offset);
			h.id = groupid+n*50+a;
			h.L0_N = groupn;
			h.taggedstart = taggedstart;
			h.ntagged = pTaggedPIDs->get_pencil_size()-taggedstart;
			h.npstart = npstart;
			h.npout = pL1Particles->get_pencil_size()-npstart;
			pL1halos->append(h);
		    } // Done with this L1 halo
		} // Done with this group
		pL1halos->FinishCell();
		pTaggedPIDs->FinishCell();
		pL1Particles->FinishCell();
		pL1PIDs->FinishCell();
	    }
	    pL1halos->FinishPencil();
	    pTaggedPIDs->FinishPencil();
	    pL1Particles->FinishPencil();
	    pL1PIDs->FinishPencil();
	    GFC->L1Tot.Stop();
	}

	// Need to update the pL1halos.npstart values for their pencil starts!
	TaggedPIDs.build_pstart();
	L1Particles.build_pstart();
	for (int j=0; j<GFC->cpd; j++) 
	    for (int k=0; k<GFC->cpd; k++) 
		for (int n=0; n<L1halos[j][k].size(); n++) {
		    HaloStat *h = L1halos[j][k].ptr(n);
		    h->npstart += L1Particles.pstart[j];
		    h->taggedstart += TaggedPIDs.pstart[j];
		}

	// Coadd the stats
	for (int g=0; g<omp_get_max_threads(); g++) {
	    GFC->L1stats.add(L1stats[g]);
	}

	// Now delete all of the temporary storage!
	for (int g=0; g<omp_get_max_threads(); g++) {
	    FOFlevel1[g].destroy();
	    FOFlevel2[g].destroy();
	    free(L1pos[g]);
	    free(L1vel[g]);
	    free(L1aux[g]);
	    free(L1acc[g]);
	}
	delete[] L1pos;
	delete[] L1vel;
	delete[] L1aux;
	delete[] L1acc;
	GFC->ProcessLevel1.Stop();
    }

};


void SimpleOutput(int slab, GlobalGroupSlab &GGS) {
    // This just writes two sets of ASCII files to see the outputs,
    // but it also helps to document how the global group information is stored.
    char fname[200];
    sprintf(fname, "/tmp/out.group.%03d", slab);
    FILE *fp = fopen(fname,"w");
    for (int j=0; j<GFC->cpd; j++)
	for (int k=0; k<GFC->cpd; k++)
	    for (int n=0; n<GGS.globalgroups[j][k].size(); n++)
		fprintf(fp, "%d %d %d %d %d\n", (int)GGS.globalgroups[j][k][n].start,
				GGS.globalgroups[j][k][n].np, GGS.slab, j, k);
    fclose(fp);

    sprintf(fname, "/tmp/out.halo.%03d", slab);
    fp = fopen(fname,"w");
    for (int j=0; j<GFC->cpd; j++)
	for (int k=0; k<GFC->cpd; k++)
	    for (int n=0; n<GGS.L1halos[j][k].size(); n++) {
		HaloStat h = GGS.L1halos[j][k][n];
		fprintf(fp, "%4d %7.4f %7.4f %7.4f %f %4d %3d %3d %7.4f %7.4f %7.4f %f %lu %lu\n", 
		    h.N, h.x[0], h.x[1], h.x[2], h.r50,
		    h.subhalo_N[0], h.subhalo_N[1], h.subhalo_N[2], 
		    h.subhalo_x[0], h.subhalo_x[1], h.subhalo_x[2], 
		    h.subhalo_r50, h.id, h.npout);
	    }
    fclose(fp);

/*
    fp = fopen(fname,"wb");
    for (int j=0; j<GFC->cpd; j++) {
	PencilAccum<TaggedPID> pids = GGS.TaggedPIDs[j];
	fwrite((void *)pids.data, sizeof(TaggedPID), pids.get_pencil_size(), fp);
    }
    fclose(fp);
*/

    // Remember that these positions are relative to the first-cell position,
    // which is why we include that cell ijk in the first outputs
    sprintf(fname, "/tmp/out.pos.%03d", slab);
    fp = fopen(fname,"w");
    for (int p=0; p<GGS.np; p++)
	fprintf(fp, "%f %f %f %d\n", GGS.pos[p].x, GGS.pos[p].y, GGS.pos[p].z, (int)GGS.aux[p].pid());
    fclose(fp);
}

uint64 GatherTaggableFieldParticles(int slab, RVfloat *pv, TaggedPID *pid) {
    // Gather all of the taggable particles that aren't in L1 groups into a file,
    // converting to global positions.
    // Warning: This must be called after ScatterGlobalGroupsAux()
    // and before ScatterGlobalGroups()
    slab = GFC->WrapSlab(slab);
    uint64 nfield = 0;
    for (int j=0; j<GFC->cpd; j++)
	for (int k=0; k<GFC->cpd; k++) {
	    // Loop over cells
	    posstruct offset = PP->CellCenter(slab, j, k);
	    Cell c = PP->GetCell(slab, j, k);
	    for (int p=0; p<c.count(); p++)
		if (c.aux[p].is_taggable() && !c.aux[p].is_L1()) {
		    // We found a taggable field particle
		    posstruct r = c.pos[p] + offset;
		    velstruct v = c.vel[p];
		    pv[nfield] = RVfloat(r.x, r.y, r.z, v.x, v.y, v.z);
		    pid[nfield] = c.aux[p].pid();
		    nfield++;
		}

	}
    return nfield;
}

#ifdef STANDALONE_FOF
void HaloOutput(int slab, GlobalGroupSlab &GGS) {
    char fname[200];
    sprintf(fname, "/tmp/out.binhalo.%03d", slab);
    GGS.L1halos.dump_to_file(fname);

    sprintf(fname, "/tmp/out.tagged.%03d", slab);
    GGS.TaggedPIDs.dump_to_file(fname);

    sprintf(fname, "/tmp/out.L1rv.%03d", slab);
    GGS.L1Particles.dump_to_file(fname);

    sprintf(fname, "/tmp/out.L1pid.%03d", slab);
    GGS.L1PIDs.dump_to_file(fname);
    // TODO: When we're ready to send this to arenas, we can use copy_to_ptr()
}

#else   // !STANDALONE_FOF
void HaloOutput(int slab, GlobalGroupSlab &GGS) {
    slab = LBW->WrapSlab(slab);

    // TODO: Need to include the control logic regarding whether 
    // a given file should be written.  
    // TODO: Need to add these arena names to slabtypes (or change them here).
    // DJE suspects there are duplications in function.

    // Write out the taggable particles not in L1 halos
    uint64 maxsize = PP->np*taggable_fraction*1.05;
    uint64 nfield = GatherTaggableFieldParticles(slab,
	(RVfloat *) LBW->AllocateSpecificSize(TaggableFieldSlab, slab, sizeof(RVfloat)
	(TaggedPID *) LBW->AllocateSpecificSize(TaggableFieldPIDSlab, slab, sizeof(TaggedPID));
    LBW->ResizeSlab(TaggableFieldSlab, slab, nfield*sizeof(RVfloat));
    LBW->ResizeSlab(TaggableFieldPIDSlab, slab, nfield*sizeof(TaggedPID));
    StoreArenaNonBlocking(TaggableFieldSlab, slab);
    StoreArenaNonBlocking(TaggableFieldPIDSlab, slab);

    if (GGS.L1halos.cpd==0) return;
    	// If this is true, then FindSubgroups() wasn't run and
	// nothing about L1 groups is even defined.  No point making
	// empty files!

    // Write out the stats on the L1 halos
    GGS.L1halos.copy_to_ptr((HaloStat *)LBW->AllocateSpecificSize(L1halosSlab, 
    		slab, GGS.L1halos.get_slab_bytes());
    StoreArenaNonBlocking(L1halosSlab, slab);

    // Write out tagged PIDs from the L1 halos
    GGS.TaggedPIDs.copy_to_ptr((TaggedPID *)LBW->AllocateSpecificSize(TaggedPIDsSlab, 
    		slab, GGS.TaggedPIDs.get_slab_bytes());
    StoreArenaNonBlocking(TaggedPIDsSlab, slab);

    // Write out the pos/vel of the taggable particles in L1 halos
    GGS.L1Particles.copy_to_ptr((RVfloat *)LBW->AllocateSpecificSize(L1ParticlesSlab, 
    		slab, GGS.L1Particles.get_slab_bytes());
    StoreArenaNonBlocking(L1ParticlesSlab, slab);

    // Write out the PIDs of the taggable particles in the L1 halos
    GGS.L1PIDs.copy_to_ptr((TaggedPID *)LBW->AllocateSpecificSize(L1PIDsSlab, 
    		slab, GGS.L1PIDs.get_slab_bytes());
    StoreArenaNonBlocking(L1PIDsSlab, slab);

    return;
}
#endif

void FindAndProcessGlobalGroups(int slab) {
    slab = GFC->WrapSlab(slab);
    GlobalGroupSlab GGS(slab);
    GGS.globalgroups.setup(GFC->cpd, 1024);   // TODO: Correct size to start?
    GGS.globalgrouplist.setup(GFC->cpd, 1024);   // TODO: Correct size to start?
    CreateGlobalGroups(GGS.slab, GGS.globalgroups, GGS.globalgrouplist);
    STDLOG(1,"Closed global groups in slab %d, finding %lld groups involving %lld cell groups\n", slab, GGS.globalgroups.size(), GGS.globalgrouplist.size());
    GFC->GGtot += GGS.globalgroups.get_slab_size();

    // Now process and output each one....
    // Start by gathering all of the particles into a contiguous set.
    GGS.GatherGlobalGroups();
    STDLOG(1,"Gathered %lld particles from global groups in slab %d\n", GGS.np, slab);
    GFC->largest_GG = std::max(GFC->largest_GG, GGS.largest_group);
    // TODO: Write to STDLOG
    // The GGS.globalgroups[j][k][n] now reference these as [start,start+np)

    // TODO: Lots to do here
    // For now, maybe we should just output the group multiplicity and the PIDs,
    // as a way of demonstrating that we have something.
    GGS.FindSubGroups();
    GGS.ScatterGlobalGroupsAux();

    #ifdef ASCII_TEST_OUTPUT
    SimpleOutput(slab, GGS);
    #endif
    HaloOutput(slab, GGS);

    GGS.ScatterGlobalGroups();
    return;
}


/* Having found the global groups, we should:

Gather a list of the particles, including the cell-offset positions.

Compute sub-group finding on that list, perhaps just returning a list of 
particle indices.

Tag the tagable particles in the L2 groups.
Compute statistics on the L1 sub-groups (do we need to make a permuted list
for this, or just dereference the indices?).
Output the L1 group information.
Output the tagable particles in the groups.

Execute microstepping on the groups.

Copy those pos/vel back into the original slabs, removing cell-offset positions.

*/







#ifdef DONT_COMPILE


// ======================= Code for the top-level driver ===============
// TODO: Move this to timestep.cpp and add to the setup and looping.
// TODO: Need to add timing and logging.

Dependency MakeCellGroups;
int MakeCellGroupsPrecondition(int slab) {
    // Need the new velocities because we're going to rearrange particles
    if( Kick.notdone(slab) ) return 0;
    
    // Also need the auxs, because we're going to re-order
    if( !LBW->IOCompleted( AuxSlab, slab ) ) {
        Dependency::NotifySpinning(WAITING_FOR_IO);
        return 0;
    }
    return 1;
    // TODO: Is it ok that the AccSlab is not being reordered now?
}

void MakeCellGroupAction(int slab) {
    GFC->CreateCellGroups(slab);
    return;
}

Dependency FindGroupLinks;
int FindGroupLinksPrecondition(int slab) {
    // We want to find all links between this slab and the one just behind
    for (int j=-1; j<=0; j++) 
	if (MakeCellGroups.notdone(slab+j)) return 0;
    return 1;
}

void FindGroupLinksAction(int slab) {
    // Find links between slab and slab-1
    FindGroupLinks(slab);
    return;
}
 
Dependency DoGlobalGroups;
int DoGlobalGroupsPrecondition(int slab) {
    // We're going to close all CellGroups in this slab.
    // GlobalGroups can span 2*GroupRadius+1.
    // But even though we usually encounter a CellGroup in its minimum slab,
    // we could be anywhere in the first instance.  So we have to query a big range.
    for (int j=-2*GFC->GroupRadius; j<=2*GFC->GroupRadius; j++) 
	if (FindGroupLinks.notdone(slab+j)) return 0;
    return 1;
}
void DoGlobalGroupsAction(int slab) {
    FindAndProcessGlobalGroups(slab);
    // TODO: LightCone Outputs might have to go here, since we need the GlobalGroups

    // At this point, all CellGroups in this slab have been closed and used
    // so we can delete them.
    GFC->DestroyCellGroups(slab);
    return;
}

#endif 
