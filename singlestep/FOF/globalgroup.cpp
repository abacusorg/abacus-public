/* This is the code that traverses the graph of links to make the GlobalGroups.
A GlobalGroup is just a list of CellGroups, stored by their LinkID.

This code creates two SlabAccums, one of the GlobalGroups and one of the
associated LinkID, with indexing between them.

This code also sorts the GLL and eventually removes all elements of that 
list that have been traversed.  It also marks all CellGroups it touches as
closed.

*/

// This class will accumulate an estimate of the amount of L1 work in each pencil,
// in the hopes that ordering the threads by this will yield better load balancing.
class PencilStats {
  public:
    float work;   // The estimated amount of work this pencil will require
    int pnum; 	// The pencil number
    int pad[14];   // To avoid cache line contention.
    
    PencilStats() { pnum = -1; work = 0.0; }
    // We provide a sort operator that will yield a decreasing list
    bool operator< (const PencilStats& c) const { return (work>c.work); }

    inline void add(float newwork) { work += newwork; }
    void reset(int _pnum) { pnum = _pnum; work = 0.0; }
};

class GlobalGroup {
    public:
    int ncellgroups;    // The number of CellGroups in this GlobalGroup
    int cellgroupstart;	// We will be storing a contiguous list of pointers to 
    	// CellGroups for all GlobalGroups in each Pencil.  This index is the 
	// offset from the start of the Pencil of that list.
	// Note that the CellGroups themselves are not necessarily in the Pencil!
	// The Pencil Number will be that of the CellGroup that initiated the
	// GlobalGroup.
    int np;    // The number of particles in this group
    uint64 start;  // The particle starting point.  
    	// Offset from the pencil beginning when constructed.
	// Updated to be offset from the slab start, when we gather the particles

    GlobalGroup(int _ncellgroups, int _cellgroupstart, int _np, int _start) {
	// Our constructor just gets the CellGroups stored.
	// We'll do any processing later.
    	ncellgroups = _ncellgroups; cellgroupstart = _cellgroupstart; 
	np = _np; start = _start;
    }
	
	bool operator>(const GlobalGroup &other) const {
		return np > other.np;
	}
};

class LinkIndex {     // These are the links for a given cell
  public:
    int start, n;     // Indexing from the start of the pencil
};

class LinkPencil {
  public:
    GroupLink *data;   // Where this pencil starts
    LinkIndex *cells; 	// [0,cpd)
    // uint64 pad[6];	// To fill up 64 bytes

    LinkPencil() { data = NULL; cells = NULL; }
    ~LinkPencil() {}

    inline void IndexPencilOld(LinkIndex *_cells, GroupLink *start, GroupLink *end, int slab, int j) {
	// GroupLinks for this Pencil begin at *start, end at *end
	// cells is where we put results, external array [0,cpd)
        cells = _cells;  // Copy the external pointer
	data = start;    // Where we'll index this pencil from
	for (int k=0; k<GFC->cpd; k++) {
	    GroupLink ref; 
	    cells[k].start = start-data;
	    // Now search for the end of this cell
	    ref.a = LinkID(slab, j, k+1, 0);    // This is the starting point we're seeking
	    while (start<end && *start<ref) start++;   // Advance to the end
	    cells[k].n = (start-data)-cells[k].start;
	}
	return;
    }

    inline void IndexPencil(LinkIndex *_cells, GroupLink *start, GroupLink *end, int slab, int j) {
	// GroupLinks for this Pencil begin at *start, end at *end
	// cells is where we put results, external array [0,cpd)
        cells = _cells;  // Copy the external pointer
	data = start;    // Where we'll index this pencil from
	uint64 iptr = 0, iend = end-start;
	for (int k=0; k<GFC->cpd; k++) {
	    GroupLink ref; 
	    uint64 ibegin = iptr;
	    cells[k].start = ibegin;
	    // Now search for the end of this cell
	    ref.a = LinkID(slab, j, k+1, 0);    // This is the starting point we're seeking
	    while (iptr<iend && start[iptr]<ref) iptr++;   // Advance to the end
	    cells[k].n = iptr-ibegin;
	}
	return;
    }

};

inline CellGroup *LinkToCellGroup(LinkID link) {
    // For this LinkID, return a pointer to the matching CellGroup
    integer3 c = link.cell();
    return GFC->cellgroups[c.x][c.y][c.z].ptr(link.cellgroup());
}



/* =================  GlobalGroupSlab ======================== */


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

    // We're going to accumulate an estimate of work in each pencil
    PencilStats *pstat;    // Will be [0,cpd)

    int largest_group;


    void setup(int _slab) {
	assertf(pstat==NULL, "GlobalGroupSlab setup is not legal");    // Otherwise, we're re-allocating something!
        pos = NULL; vel = NULL; aux = NULL; acc = NULL; np = 0;
	slab = GFC->WrapSlab(_slab); largest_group = 0;
	int ret = posix_memalign((void **) &pstat, 64, sizeof(PencilStats)*GFC->cpd);
	assert(ret==0);
	for (int j=0; j<GFC->cpd; j++) pstat[j].reset(j);
	globalgroups.setup(GFC->cpd, GFC->particles_per_pencil);   
	globalgrouplist.setup(GFC->cpd, GFC->particles_per_pencil);   
    }
    void destroy() {
        if (pos!=NULL) free(pos); pos = NULL;
        if (vel!=NULL) free(vel); vel = NULL;
        if (aux!=NULL) free(aux); aux = NULL;
        if (acc!=NULL) free(acc); acc = NULL;
	if (pstat!=NULL) free(pstat); pstat = NULL;
	np = 0;
    }

    void allocate(uint64 _np) {
	np = _np;
	// LONG-TERM: May eventually prefer these to be arenas.
	int ret;
        ret = posix_memalign((void **)&pos, 4096, sizeof(posstruct)*np); assert(ret==0);
        ret = posix_memalign((void **)&vel, 4096, sizeof(velstruct)*np); assert(ret==0);
        ret = posix_memalign((void **)&aux, 4096, sizeof(auxstruct)*np); assert(ret==0);
        ret = posix_memalign((void **)&acc, 4096, sizeof(accstruct)*np); assert(ret==0);
    }

    GlobalGroupSlab() { pstat = NULL; }    // We have a null constructor
    ~GlobalGroupSlab() { destroy(); }

    void CreateGlobalGroups();
    void GatherGlobalGroups();
    void ScatterGlobalGroupsAux();
    void ScatterGlobalGroups();
    void FindSubGroups();
    void SimpleOutput();
    void HaloOutput();

};




void GlobalGroupSlab::CreateGlobalGroups() {
    // For this slab, we want to traverse the graph of links to find all GlobalGroups.
    // Given pstat[0,cpd)
    assertf(pstat!=NULL,"setup() not yet called\n");   // Setup not yet called!
    GFC->SortLinks.Start();

    // The groups can span 2*R+1 slabs, but we might be on either end,
    // so we have to consider 4*R+1.
    int rad = GFC->GroupRadius;	
    int diam = 4*rad+1;
    int cpd = GFC->cpd;
    // int cpdpad = (cpd/8+1)*8;   // A LinkIndex is 8 bytes, so let's get each pencil onto a different cacheline
    int cpdpad = cpd;
    
    // We're going to sort the GLL and then figure out the starting index 
    // of every cell in this slab.
    GFC->GLL->Sort();		// Sorts by LinkID
    GFC->SortLinks.Stop();
    GFC->IndexLinks.Start();
    // GFC->GLL->AsciiPrint();
    LinkPencil **links;		// For several slabs and all pencils in each
    links = (LinkPencil **)malloc(sizeof(LinkPencil *)*diam);
    links[0] = (LinkPencil *)malloc(sizeof(LinkPencil)*diam*cpd);

    {int ret = posix_memalign((void **)&(links[0]), 64, sizeof(LinkPencil)*diam*cpd); assert(ret==0);}
    for (int j=1; j<diam; j++) links[j] = links[j-1]+cpd;
    LinkIndex *cells;
    {int ret = posix_memalign((void **)&cells, 64, sizeof(LinkIndex)*diam*cpd*cpdpad); assert(ret==0);}
    // cells = (LinkIndex *)malloc(sizeof(LinkIndex)*diam*cpd*cpdpad);
    for (int s=0; s<diam; s++) {
	int thisslab = GFC->WrapSlab(slab+s-diam/2);
	assertf(GFC->cellgroups_status[thisslab]>0, "Cellgroup slab not present.  Something is wrong in dependencies!");
	    // Just to check that the CellGroups are present or already closed.
	#pragma omp parallel for schedule(dynamic,1)
	for (int j=0; j<cpd; j++) {
	    // Now find the starting point for this Pencil
	    GroupLink *start = GFC->GLL->Search(thisslab, j);
	    // printf("Pencil %d %d starts at %d\n", thisslab, j, (int)(start-GFC->GLL->list));
	    // Find the start for each Cell in the pencil
	    // This step does dominate the time
	    GFC->IndexLinksIndex.Start();
	    links[s][j].IndexPencil(cells+(s*cpd+j)*cpdpad, start, GFC->GLL->list+GFC->GLL->length, thisslab, j);
	    GFC->IndexLinksIndex.Stop();

	    // for (int k=0; k<cpd; k++) 
		// printf("%d %d %d starts at %d %d\n", thisslab, j, k, links[s][j].cells[k].start, links[s][j].cells[k].n);

	}
    }
    GFC->IndexLinks.Stop();
    GFC->FindGlobalGroupTime.Start();

    // Then need to go cell by cell to take unclosed groups and try to close them.
    // That step will need to be run every (2R+1) pencil to avoid from contention.
    // 4R+1 might be better, just to avoid any cache overlap.
    // Have to guard against contention across the periodic wrap, as well.
    // Pick a division no less than diam that divides into CPD so that 
    // we can skip this issue.
    int split = diam;
    MultiplicityStats L0stats[omp_get_max_threads()];
    while (split<cpd && cpd%split!=0) split++;
    for (int w=0; w<split; w++) {
	#pragma omp parallel for schedule(dynamic,1) 
	for (int j=w; j<cpd; j+=split) {
	    // Do this pencil
	    int cumulative_np = 0;    // We count all particles in this pencil
	    PencilAccum<GlobalGroup> *gg_pencil = globalgroups.StartPencil(j);
	    PencilAccum<LinkID> *gg_list = globalgrouplist.StartPencil(j);
	    std::vector<LinkID> cglist;    // List of CellGroups in this GlobalGroup
	    cglist.reserve(64);
	    for (int k=0; k<cpd; k++) {
		// Do this cell
		CellPtr<CellGroup> c = GFC->cellgroups[slab][j][k];
		for (int g=0; g<c.size(); g++) if (c[g].is_open()) {
		    // Loop over groups, skipping closed ones
		    int ggsize = 0; 
		    cglist.resize(0);   // Reset to null list
		    cglist.push_back(LinkID(slab, j, k, g));   // Prime the queue
		    int searching = 0;
		    // printf("GG %d %d %d %d\n", slab, j, k, g);
		    while (searching<cglist.size()) {
		        // We have an unresolved group to search
			integer3 thiscell = cglist[searching].cell();  
			CellGroup *thiscg = LinkToCellGroup(cglist[searching]);
			ggsize += thiscg->size();
			assertf(thiscg->is_open(), "Cellgroup found to be closed.  Likely need to increase GroupRadius!\n");
				// If this fails, probably a group has spanned
				// beyond 2*R+1 cells and something got closed
				// prematurely.
			thiscg->close_group();
			// if (searching>0) {
			    // printf("Link to %d %d %d %d\n",
			    	// thiscell.x, thiscell.y, thiscell.z, cglist[searching].cellgroup());
			// }
			// Now get these links
			int s = GFC->WrapSlab(thiscell.x-slab+diam/2);  // Map to [0,diam)
			LinkPencil *lp = links[s]+thiscell.y;
			GroupLink *g; int t;
			for (g = lp->data + lp->cells[thiscell.z].start, t = 0;
				t<lp->cells[thiscell.z].n; t++, g++) {
			    // Consider all the links from this cell
			    if (g->a.id == cglist[searching].id) {
				// This link originates from the group we're using
				// But now we need to check if its destination
				// is already on the list.
				int seek = 0;
				while (seek<cglist.size()) 
				    if (g->b.id == cglist[seek++].id) goto FoundIt;
				cglist.push_back(g->b);  // Didn't find it; add to queue
				FoundIt:
				// Either way, this link has been used its one time.
				g->mark_for_deletion();
			    }
			} // Done with the links from this cell
			searching++;  // Ready for the next unresolved cellgroup
		    } // End search over graph of cell groups
		    // printf("Closed with size %d\n", ggsize);

		    // Time to build the GlobalGroup
		    // We only track the group if it has more than one particle
		    if (ggsize>1) {
			int start = gg_list->buffer->get_pencil_size();
			for (int t = 0; t<cglist.size(); t++) {
			    integer3 tmp = cglist[t].cell();
			    // printf("GGlist: %d %d %d %d\n",
			    	// tmp.x, tmp.y, tmp.z, cglist[t].cellgroup());
			    gg_list->append(cglist[t]);
			}
			// printf("    GGpencil: %d %d %d %d\n", 
				// (int)cglist.size(), start, ggsize, cumulative_np);
			gg_pencil->append(GlobalGroup(cglist.size(), start, ggsize, cumulative_np));
			cumulative_np += ggsize;
			L0stats[omp_get_thread_num()].push(ggsize);
			pstat[j].add(ggsize*ggsize);
			    // Track the later work, assuming scaling as N^2
		    }
		} // End this group
		gg_pencil->FinishCell();
		gg_list->FinishCell();
		// printf("Done with cell %d %d\n", j, k);
	    } // End this cell
	    gg_pencil->FinishPencil();
	    gg_list->FinishPencil();
	} // End this pencil
    } // End loop over non-contending sub-slabs

    // Get rid of the links that have been used.
    GFC->GLL->PartitionAndDiscard();

    // Cumulate the L0 statistics
    for (int j=0; j<omp_get_max_threads(); j++) GFC->L0stats.add(L0stats[j]);

    // Free the indexing space
    free(cells);
    free(links[0]);
    free(links);
    GFC->FindGlobalGroupTime.Stop();
}




void GlobalGroupSlab::GatherGlobalGroups() {
    // Copy all of the particles into a single list
    GFC->IndexGroups.Start();
    int diam = 2*GFC->GroupRadius+1;
    assertf(GFC->cpd>=4*GFC->GroupRadius+1, "CPD is too small compared to GroupRadius\n");
    // This registers the periodic wrap using the cells.
    // However, this will misbehave if CPD is smaller than the group diameter,
    // because the LinkIDs have already been wrapped.
    // One just can't use a small CPD
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
    allocate(total_particles);
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




void GlobalGroupSlab::ScatterGlobalGroupsAux() {
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




void GlobalGroupSlab::ScatterGlobalGroups() {
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





void GlobalGroupSlab::FindSubGroups() {
    // Process each group, looking for L1 and L2 subgroups.
    GFC->ProcessLevel1.Start();
    L1halos.setup(GFC->cpd, GFC->particles_per_pencil);    // TUNING: This start value is a guess
    TaggedPIDs.setup(GFC->cpd, GFC->particles_per_pencil);    
    L1Particles.setup(GFC->cpd, GFC->particles_per_pencil);   
    L1PIDs.setup(GFC->cpd, GFC->particles_per_pencil);    
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

    // It seems that the work between pencils is so heterogeneous that even the
    // dynamic scheduling can't smooth it out.  So we're going to try ordering the
    // pencils by the work estimate (largest first)
    std::sort(pstat, pstat+GFC->cpd);
    
    // for (int j=0; j<GFC->cpd; j++) 
    #pragma omp parallel for schedule(dynamic,1)
    for (int jj=0; jj<GFC->cpd; jj++) {
	int j = pstat[jj].pnum;    // Get the pencil number from the list
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
		    GFC->L2FOF.Stop();

		    // Merger trees require tagging the taggable particles 
		    // of the biggest L2 group in the original aux array.  
		    // This can be done:
		    // The L2 index() gives the position in the L1 array,
		    // and that index() gets back to aux.
			if(FOFlevel2[g].ngroups > 0){
				std::sort(FOFlevel2[g].groups, FOFlevel2[g].groups+FOFlevel2[g].ngroups);
				// Groups now in descending order of multiplicity
			
				FOFparticle *L2start = FOFlevel2[g].p + FOFlevel2[g].groups[0].start;
				for (int p=0; p<FOFlevel2[g].groups[0].n; p++) {
				if (groupaux[start[L2start[p].index()].index()].is_taggable())
					groupaux[start[L2start[p].index()].index()].set_tagged();
				}
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





void GlobalGroupSlab::SimpleOutput() {
    // This just writes two sets of ASCII files to see the outputs,
    // but it also helps to document how the global group information is stored.
    char fname[200];
    sprintf(fname, "/tmp/out.group.%03d", slab);
    FILE *fp = fopen(fname,"w");
    for (int j=0; j<GFC->cpd; j++)
	for (int k=0; k<GFC->cpd; k++)
	    for (int n=0; n<globalgroups[j][k].size(); n++)
		fprintf(fp, "%d %d %d %d %d\n", (int)globalgroups[j][k][n].start,
				globalgroups[j][k][n].np, slab, j, k);
    fclose(fp);

    sprintf(fname, "/tmp/out.halo.%03d", slab);
    fp = fopen(fname,"w");
    for (int j=0; j<GFC->cpd; j++)
	for (int k=0; k<GFC->cpd; k++)
	    for (int n=0; n<L1halos[j][k].size(); n++) {
		HaloStat h = L1halos[j][k][n];
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
	PencilAccum<TaggedPID> pids = TaggedPIDs[j];
	fwrite((void *)pids.data, sizeof(TaggedPID), pids.get_pencil_size(), fp);
    }
    fclose(fp);
*/

    // Remember that these positions are relative to the first-cell position,
    // which is why we include that cell ijk in the first outputs
    sprintf(fname, "/tmp/out.pos.%03d", slab);
    fp = fopen(fname,"w");
    for (int p=0; p<np; p++)
	fprintf(fp, "%f %f %f %d\n", pos[p].x, pos[p].y, pos[p].z, (int)aux[p].pid());
    fclose(fp);
}

#ifdef STANDALONE_FOF
void GlobalGroupSlab::HaloOutput() {
    char fname[200];
    sprintf(fname, "/tmp/out.binhalo.%03d", slab);
    L1halos.dump_to_file(fname);

    sprintf(fname, "/tmp/out.tagged.%03d", slab);
    TaggedPIDs.dump_to_file(fname);

    sprintf(fname, "/tmp/out.L1rv.%03d", slab);
    L1Particles.dump_to_file(fname);

    sprintf(fname, "/tmp/out.L1pid.%03d", slab);
    L1PIDs.dump_to_file(fname);
}

#else   // !STANDALONE_FOF
void GlobalGroupSlab::HaloOutput() {
	/* Outputs a uniform subsample of the particles ("taggable particles")
	 * and L1 halo stats and tagged particles. Currently we always output
	 * taggable particles if we are also outputting halos. Halo outputs will
	 * be skipped if no L1 halos were found (but taggable particles will still be written).
	 */
	 
    // TODO: Need to include the control logic regarding whether 
    // a given file should be written.
	STDLOG(0,"Beginning halo output for slab %d\n", slab);
	
	// Ensure the output directory for this step exists
	char dir[16];
	sprintf(dir, "Step%04d", ReadState.FullStepNumber);
	CreateSubDirectory(P.GroupDirectory, dir);

    // Write out the taggable particles not in L1 halos
	// TODO: better heuristic? what will happen in very small sims?  Also technically HaloTaggableFraction is only used in the IC step
	uint64 maxsize = P.np*P.HaloTaggableFraction*1.05;
	LBW->AllocateSpecificSize(TaggableFieldSlab, slab, maxsize*sizeof(RVfloat));
	LBW->AllocateSpecificSize(TaggableFieldPIDSlab, slab, maxsize*sizeof(TaggedPID));
	
    uint64 nfield = GatherTaggableFieldParticles(slab,
        (RVfloat *) LBW->ReturnIDPtr(TaggableFieldSlab, slab),
        (TaggedPID *) LBW->ReturnIDPtr(TaggableFieldPIDSlab, slab));
	if(nfield > 0){
		// only write the uniform subsample files if they will have non-zero size
		LBW->ResizeSlab(TaggableFieldSlab, slab, nfield*sizeof(RVfloat));
		LBW->ResizeSlab(TaggableFieldPIDSlab, slab, nfield*sizeof(TaggedPID));
		LBW->StoreArenaNonBlocking(TaggableFieldSlab, slab);
		LBW->StoreArenaNonBlocking(TaggableFieldPIDSlab, slab);
	} else {
		LBW->DeAllocate(TaggableFieldSlab, slab);
		LBW->DeAllocate(TaggableFieldPIDSlab, slab);
	}
		
    if (L1halos.pencils == NULL || L1halos.get_slab_size() == 0) return;
	// If pencils is NULL, then FindSubgroups() wasn't run and
	// nothing about L1 groups is even defined.
	// And if get_slab_size() is 0, then we found ran FOF but found nothing.
	// No point making empty files!

    // Write out the stats on the L1 halos
	LBW->AllocateSpecificSize(L1halosSlab, slab, L1halos.get_slab_bytes());
    L1halos.copy_to_ptr((HaloStat *)LBW->ReturnIDPtr(L1halosSlab, slab));
    LBW->StoreArenaNonBlocking(L1halosSlab, slab);

    // Write out tagged PIDs from the L1 halos
	LBW->AllocateSpecificSize(TaggedPIDsSlab, slab, TaggedPIDs.get_slab_bytes());
    TaggedPIDs.copy_to_ptr((TaggedPID *)LBW->ReturnIDPtr(TaggedPIDsSlab, slab));
    LBW->StoreArenaNonBlocking(TaggedPIDsSlab, slab);

    // Write out the pos/vel of the taggable particles in L1 halos
	LBW->AllocateSpecificSize(L1ParticlesSlab, slab, L1Particles.get_slab_bytes());
    L1Particles.copy_to_ptr((RVfloat *)LBW->ReturnIDPtr(L1ParticlesSlab, slab));
    LBW->StoreArenaNonBlocking(L1ParticlesSlab, slab);

    // Write out the PIDs of the taggable particles in the L1 halos
	LBW->AllocateSpecificSize(L1PIDsSlab, slab, L1PIDs.get_slab_bytes());
    L1PIDs.copy_to_ptr((TaggedPID *)LBW->ReturnIDPtr(L1PIDsSlab, slab));
    LBW->StoreArenaNonBlocking(L1PIDsSlab, slab);

    return;
}
#endif





