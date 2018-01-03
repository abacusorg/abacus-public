/* This is the code that traverses the graph of links to make the GlobalGroups.
A GlobalGroup is just a list of CellGroups, stored by their LinkID.

This code creates two SlabAccums, one of the GlobalGroups and one of the
associated LinkID, with indexing between them.

This code also sorts the GLL and eventually removes all elements of that 
list that have been traversed.  It also marks all CellGroups it touches as
closed.

*/

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
};

class LinkIndex {     // These are the links for a given cell
  public:
    int start, n;     // Indexing from the start of the pencil
};

class LinkPencil {
  public:
    GroupLink *data;   // Where this pencil starts
    LinkIndex *cells; 	// [0,cpd)

    LinkPencil() { data = NULL; cells = NULL; }
    ~LinkPencil() {}

    void IndexPencil(LinkIndex *_cells, GroupLink *start, GroupLink *end, int slab, int j) {
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
};

inline CellGroup *LinkToCellGroup(LinkID link) {
    // For this LinkID, return a pointer to the matching CellGroup
    integer3 c = link.cell();
    return GFC->cellgroups[c.x][c.y][c.z].ptr(link.cellgroup());
}

void CreateGlobalGroups(int slab,
		SlabAccum<GlobalGroup> &globalgroups,
		SlabAccum<LinkID> &globalgrouplist) {
    // For this slab, we want to traverse the graph of links to find all GlobalGroups
    GFC->SortLinks.Start();
    slab = GFC->WrapSlab(slab);

    // The groups can span 2*R+1 slabs, but we might be on either end,
    // so we have to consider 4*R+1.
    int rad = GFC->GroupRadius;	
    int diam = 4*rad+1;
    int cpd = GFC->cpd;
    
    // We're going to sort the GLL and then figure out the starting index 
    // of every cell in this slab.
    GFC->GLL->Sort();		// Sorts by LinkID
    GFC->SortLinks.Stop();
    GFC->IndexLinks.Start();
    // GFC->GLL->AsciiPrint();
    LinkPencil **links;		// For several slabs and all pencils in each
    links = (LinkPencil **)malloc(sizeof(LinkPencil *)*diam);
    links[0] = (LinkPencil *)malloc(sizeof(LinkPencil)*diam*cpd);
    for (int j=1; j<diam; j++) links[j] = links[j-1]+cpd;
    LinkIndex *cells;
    cells = (LinkIndex *)malloc(sizeof(LinkIndex)*diam*cpd*cpd);
    for (int s=0; s<diam; s++) {
	int thisslab = GFC->WrapSlab(slab+s-diam/2);
	assert(GFC->cellgroups_status[thisslab]>0);
	    // Just to check that the CellGroups are present or already closed.
	#pragma omp parallel for schedule(dynamic,1) 
	for (int j=0; j<cpd; j++) {
	    // Now find the starting point for this Pencil
	    GroupLink *start = GFC->GLL->Search(thisslab, j);
	    // printf("Pencil %d %d starts at %d\n", thisslab, j, (int)(start-GFC->GLL->list));
	    // Find the start for each Cell in the pencil
	    links[s][j].IndexPencil(cells+(s*cpd+j)*cpd, start, GFC->GLL->list+GFC->GLL->length, thisslab, j);

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
			assert(thiscg->is_open());
				// If this fails, probably a group has spanned
				// beyond 2*R+1 cells and something got closed
				// prematurely.
			thiscg->close_group();
			if (searching>0) {
			    // printf("Link to %d %d %d %d\n",
			    	// thiscell.x, thiscell.y, thiscell.z, cglist[searching].cellgroup());
			}
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

    // Free the indexing space
    free(cells);
    free(links[0]);
    free(links);
    GFC->FindGlobalGroupTime.Stop();
}
