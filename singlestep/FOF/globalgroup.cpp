/** \file This is the code that traverses the graph of links to make the GlobalGroups.
A GlobalGroup is just a list of CellGroups, stored by their LinkID.

This code creates two SlabAccums, one of the GlobalGroups and one of the
associated LinkID, with indexing between them.

This code also sorts the GLL and eventually removes all elements of that 
list that have been traversed.  It also marks all CellGroups it touches as
closed.

*/

/// This class will accumulate an estimate of the amount of L1 work in each pencil,
/// in the hopes that ordering the threads by this will yield better load balancing.
class PencilStats {
  public:
    float work;   // The estimated amount of work this pencil will require
    int pnum;         // The pencil number
    uint8_t pad[CACHE_LINE_SIZE-8];   // To avoid cache line contention.
    
    PencilStats() { pnum = -1; work = 0.0; }
    // We provide a sort operator that will yield a decreasing list
    bool operator< (const PencilStats& c) const { return (work>c.work); }

    inline void add(float newwork) { work += newwork; }
    void reset(int _pnum) { pnum = _pnum; work = 0.0; }
};

/** A GlobalGroup contains the information about one group,
composed of multiple CellGroups.
*/

class GlobalGroup {
    public:
    int ncellgroups;    ///< The number of CellGroups in this GlobalGroup
    int cellgroupstart;        ///< We will be storing a contiguous list of pointers to 
            ///< CellGroups for all GlobalGroups in each Pencil.  This index is the 
        ///< offset from the start of the Pencil of that list.
        ///< Note that the CellGroups themselves are not necessarily in the Pencil!
        ///< The Pencil Number will be that of the CellGroup that initiated the
        ///< GlobalGroup.
    int np;    ///< The number of particles in this group
    uint64 start;  ///< The particle starting point.  
            ///< Offset from the pencil beginning when constructed.
        ///< Updated to be offset from the slab start, when we gather the particles

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

/// This is the lookup table into the GroupLink list for a given
/// pencil.
class LinkPencil {
  public:
    GroupLink *data;   ///< Where this pencil starts
    LinkIndex *cells;         ///< [0,zwidth)
    // uint8_t pad[CACHE_LINE_SIZE];        // To fill up 64 bytes

    LinkPencil() { data = NULL; cells = NULL; }
    ~LinkPencil() {}

    inline void IndexPencil(LinkIndex *_cells, GroupLink *start, GroupLink *end, int slab, int j) {
        // GroupLinks for this Pencil begin at *start, end at *end
        // cells is where we put results, external array [0,zwidth)
        cells = _cells;  // Copy the external pointer
        data = start;    // Where we'll index this pencil from
        uint64 iptr = 0, iend = end-start;
        for (int k=0; k<GFC->zwidth; k++) {  // local k
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
    integer3 c = link.localcell();
    assertf(GFC->cellgroups_status[c.x] == 1, "Failed to find cellgroup in slab %d with status %d\n", c.x, GFC->cellgroups_status[c.x]);
    CellGroup *ret = GFC->cellgroups[c.x][c.y][c.z].ptr(link.cellgroup());
    assertf(ret != NULL, "Bad LinkToCellGroup?\n");
    return ret;
}



/* =================  GlobalGroupSlab ======================== */


/** This contains all of the GlobalGroups in a given slab.
It contains the methods to construct those Global Groups from the 
GroupLinkList, as well as the methods to gather/scatter the particles
and then analyze the groups.

Belonging to a slab means that it contains a CellGroup that hadn't
been previously included into a GlobalGroup in a different slab.

In past versions, we searched all remaining CellGroups in the supplied
slab, meaning that when we were done with this search, all cellgroups in 
this slab were completed, i.e., belonged to GlobalGroups either in
this slab or a previous one.

In the new version, we are only allowing the GlobalGroup search to
run in the positive direction in slab number, searching out to
S+GroupDiam.  If a search would take one to slab-1, then we mark
those cellgroups as deferred for the present search and will return
to them in a future slab's attempt.  The CellGroups in slab S will
all have been used when CreateGlobalGroups() has been called on
slab=S through S-GroupDiam.  It remains the case that each GlobalGroup
is found only once, and each CellGroup is used only once.
*/

class GlobalGroupSlab {
  public:
    SlabAccum<GlobalGroup> globalgroups;
            ///< The global group information
    SlabAccum<LinkID> globalgrouplist;
            ///< The cell group decompositions of these global groups

    // The following accumulate possible output 
    SlabAccum<HaloStat> L1halos;              ///< Stats about each L1 halo
    SlabAccum<RVfloat> HaloRVA;       ///< The taggable subset in each L1 halo, pos/vel, subsample A
    SlabAccum<RVfloat> HaloRVB;       ///< The taggable subset in each L1 halo, pos/vel, subsample B
    SlabAccum<TaggedPID> HaloPIDsA;        ///< The taggable subset in each halo, PID, subsample A. If this is an L1OutputRedshift, we only output L1 halo pids. If this is a subsample redshift, we include the PIDs of the L0 halos.
    SlabAccum<TaggedPID> HaloPIDsB;        ///< The taggable subset in each halo, PID, subsample B. If this is an L1OutputRedshift, we only output L1 halo pids. If this is a subsample redshift, we include the PIDs of the L0 halos.
    
    int slab;    ///< This slab number
    int rad; 
    int diam, slabbias;
    int split;    ///< The pencil splitting
    int cpd;
    int zstart;  ///< includes ghosts
    int zwidth;
    
    int ghost_radius;  ///< for determining 2D ownership
    int dedup;  ///< whether to de-duplicate 2D groups

    posstruct *pos;  ///< Big vectors of all of the pos/vel/aux for these global groups
    velstruct *vel;
    auxstruct *aux;
    accstruct *acc;
    uint64 np;
    uint64 ncellgroups;

    LinkPencil **links;                // For several slabs and all pencils in each
    LinkIndex *cells;               // For several slabs and all cells in each

    // We're going to accumulate an estimate of work in each pencil
    PencilStats *pstat;    ///< Will be [0,cpd)

    int largest_group;

    uint64 *pstart;
    uint64 *pstart_cg;


    /// We have a boring constructor; call setup() to actually get started.
    void setup(int _slab) {
        assertf(pstat==NULL, "GlobalGroupSlab setup is not legal");    // Otherwise, we're re-allocating something!
        pos = NULL; vel = NULL; aux = NULL; acc = NULL; np = 0;
        slab = GFC->WrapSlab(_slab); largest_group = 0;

        cpd = GFC->cpd;
        zstart = GFC->zstart;
        zwidth = GFC->zwidth;
        ghost_radius = GHOST_RADIUS;
        dedup = MPI_size_z > 1;
        rad = GFC->GroupRadius;
        // We will split the work by pencils, and this means that
        // we have to reserve 2*r+1 in either direction.
        split = 4*rad+1;
        while (split<cpd && (cpd%split)!=0) split++;

        #ifdef ONE_SIDED_GROUP_FINDING
            // The groups can span 2*R+1 slabs.  
            // We will only look in the upward direction.
            diam = 2*rad+1; 
            slabbias = 0;
        #else
            // The groups can span 2*R+1 slabs in either direction,
            // so 4*R+1
            diam = 4*rad+1; 
            slabbias = diam/2;
        #endif

        int ret = posix_memalign((void **) &pstat, CACHE_LINE_SIZE, sizeof(PencilStats)*cpd);
        assert(ret==0);
        for (int j=0; j<cpd; j++) pstat[j].reset(j);

        // These are counters for where the pencils start in the lists
        pstart = new uint64[cpd];
        pstart_cg = new uint64[cpd];  

        // How many global groups do we expect?  Let's guess 10% of the particles
        // The cellgroup list is modestly bigger, but not enough to care.
        int maxsize = GFC->particles_per_slab/10;
        globalgroups.setup(cpd, zwidth, maxsize);
        globalgrouplist.setup(cpd, zwidth, maxsize);
    }
    void destroy() {
        if (pos!=NULL) free(pos); pos = NULL;
        if (vel!=NULL) free(vel); vel = NULL;
        if (aux!=NULL) free(aux); aux = NULL;
        if (acc!=NULL) free(acc); acc = NULL;
        if (pstat!=NULL) free(pstat); pstat = NULL;
        np = 0;
        delete[] pstart;
        delete[] pstart_cg;
        
        globalgroups.destroy();
        globalgrouplist.destroy();
        L1halos.destroy();
        HaloRVA.destroy();
        HaloRVB.destroy();
        HaloPIDsA.destroy();
        HaloPIDsB.destroy();
    }

    void allocate(uint64 _np) {
        np = _np;
        // LONG-TERM: May eventually prefer these to be arenas.
        int ret;
        ret = posix_memalign((void **)&pos, PAGE_SIZE, sizeof(posstruct)*np); assert(ret==0);
        ret = posix_memalign((void **)&vel, PAGE_SIZE, sizeof(velstruct)*np); assert(ret==0);
        ret = posix_memalign((void **)&aux, PAGE_SIZE, sizeof(auxstruct)*np); assert(ret==0);
        if(SB->IsSlabPresent(AccSlab,slab))  // leave NULL if no acc
            ret = posix_memalign((void **)&acc, PAGE_SIZE, sizeof(accstruct)*np); assert(ret==0);
    }

    GlobalGroupSlab() { pstat = NULL; }    // We have a null constructor
    ~GlobalGroupSlab() { destroy(); }

    void IndexLinks();
    void DeferGlobalGroups();
    void ClearDeferrals();
    void AppendParticleToPencil(PencilAccum<RVfloat> ** pHaloRVs, PencilAccum<TaggedPID> ** pHaloPIDs, 
                                            posstruct * grouppos, velstruct * groupvel, accstruct * groupacc, auxstruct * groupaux, 
                                            int index, posstruct offset, int & nA, int & nB);
    void CreateGlobalGroups();
    void GatherGlobalGroups();
    void ScatterGlobalGroupsAux();
    void ScatterGlobalGroups();
    void FindSubGroups();
    void SimpleOutput();
    void HaloOutput();
    void WriteGroupHeaderFile(const char* fn);
    uint64 L0TimeSliceOutput(FLOAT unkick_factor);

};




/** For this slab, we want to traverse the graph of links to find all GlobalGroups.

This involves going to every unassigned CellGroup and then looking
for its neighbors in the GroupLinkList.  By doing a FOF-like repeating
search, we build up all of the CGs that are associated.  These become
a GlobalGroup.  

As links are used, they are marked for deletion from the GLL.
As CGs are included, they are marked as assigned, so they can't
start another GG.

GlobalGroups cannot span more than 2*GroupRadius+1 cells, and we
use this fact to parallelize the work so that threads are working
on distinct portions of the slab.

We also make heavy use of SlabAccum's to order the work by Pencil;
this will aid in the later bookkeeping.
*/

void GlobalGroupSlab::IndexLinks() {
    // We're going to sort the GLL and then figure out the starting index 
    // of every cell in this slab.
    // We only care about links that initiate from slabs between S and S+diam.
    GFC->SortLinks.Start();
    GFC->GLL->Sort();                // Sorts by LinkID
    GFC->SortLinks.Stop();

    GFC->IndexLinks.Start();
    // GFC->GLL->AsciiPrint();
    links = new LinkPencil*[diam];

    // int cpdpad = (cpd/8+1)*8;   // A LinkIndex is 8 bytes, so let's get each pencil onto a different cacheline
    int zpad = zwidth;
    // Allocate space for the links[slab][y] LinkPencils.
    {int ret = posix_memalign((void **)&(links[0]), CACHE_LINE_SIZE, sizeof(LinkPencil)*diam*cpd); assert(ret==0);}
    for (int j=1; j<diam; j++) links[j] = links[j-1]+cpd;
    // Allocate space for the LinkIndex array, which is on a per-cell basis.
    {int ret = posix_memalign((void **)&cells, CACHE_LINE_SIZE, sizeof(LinkIndex)*diam*cpd*zpad); assert(ret==0);}

    // Now loop over slabs to construct the lookup tables
    for (int s=0; s<diam; s++) {
        int thisslab = GFC->WrapSlab(slab+s-slabbias);
        assertf(GFC->cellgroups_status[thisslab]>0, "Cellgroup slab %d not present (value %d).  Something is wrong in dependencies!\n", thisslab, GFC->cellgroups_status[thisslab]);
            // Just to check that the CellGroups are present or already closed.
        NUMA_FOR(j,0,cpd)
            // Now find the starting point for this Pencil
            GroupLink *start = GFC->GLL->Search(thisslab, j);
            // printf("Pencil %d %d starts at %d\n", thisslab, j, (int)(start-GFC->GLL->list));
            // Find the start for each Cell in the pencil
            // This step does dominate the time
            GFC->IndexLinksIndex.Start();
            links[s][j].IndexPencil(cells+(s*cpd+j)*zpad, start, GFC->GLL->list+GFC->GLL->length, thisslab, j);
            GFC->IndexLinksIndex.Stop();

            // for (int k=0; k<zwidth; k++) 
                // printf("%d %d %d starts at %d %d\n", thisslab, j, k, links[s][j].cells[k].start, links[s][j].cells[k].n);

        }
    }
    GFC->IndexLinks.Stop();
}



void GlobalGroupSlab::DeferGlobalGroups() {
    /* We want to find all of the global groups that cross to S-1 and mark those
    cell groups as deferred.  To do this, we start by finding links that cross
    from S to S-1.  We take the CellGroups marked by the S end of the link and 
    follow the links for that group, considering only links to >=S.  Note that 
    some GlobalGroups could be split into pieces when we ignore the connections
    through <=S-1, but we don't care because we will still find all of the pieces
    in the domain >=S and mark them as deferred.
        Once we have the CellGroups deferred, then we can ignore those CellGroups
    as starting points for valid GlobalGroups in this slab.
        At the end, we must clear the deferrals.
        // TODO: Or we could do it before we start?
    */
    STDLOG(1, "Starting Deferral search on slab %d\n", slab);
    GFC->DeferGroups.Start();
    int lastslab = GFC->WrapSlab(slab-slabbias-1);

    for (int w=0; w<split; w++) {
        #pragma omp parallel for schedule(dynamic,1)
        for (int j=w; j<cpd; j+=split) {
            // Do this pencil within the boundary slab.  
            LinkPencil *slablp = links[0]+j;
            std::vector<LinkID> cglist;    // List of CellGroups in this GlobalGroup
            cglist.reserve(64);
            for (int k=0; k<zwidth; k++) {   // Loop over cells, local k
                GroupLink *gg; int t;
                // Loop over links from this cell
                for (gg = slablp->data + slablp->cells[k].start, t = 0;
                    t<slablp->cells[k].n; t++, gg++) {
                    // If the link doesn't point to S-1, skip to next
                    integer3 linkcell = gg->b.localcell();
                    if (linkcell.x != lastslab) continue;
                    // If this CellGroup is not open, we could skip the search
                    
                    // We've found a CellGroup that connects to slab S-1,
                    // so we need to track it.
                    cglist.resize(0);   // Reset to null list
                    cglist.push_back(gg->a);   // Prime the queue
                    uint64 searching=0;

                    while (searching<cglist.size()) {
                        // We have an unresolved group to search
                        integer3 thiscell = cglist[searching].localcell();  
                        CellGroup *thiscg = LinkToCellGroup(cglist[searching]);
                        // There's a chance we've already been here; if yes, skip
                        if (!(thiscg->is_open())) { searching++; continue; }

                        int s = GFC->WrapSlab(thiscell.x-slab+slabbias);  // Map to [0,diam)
                        assertf(s < diam, "GroupLink points to slab offset %d which violates GroupDiameter %d.  Likely need to increase GroupRadius! (slab,j,k,t = %d,%d,%d,%d)\n",
                            s, diam, slab, j, k, t);
                        LinkPencil *lp = links[s]+thiscell.y;

                        // Loop over the links that leave from this cell
                        GroupLink *g; int t;
                        for (g = lp->data + lp->cells[thiscell.z].start, t = 0;
                                t<lp->cells[thiscell.z].n; t++, g++) {
                            if (g->a.id == cglist[searching].id) {
                                // The link originates from the group we're considering
                                // But now we need to check its destination
                                uint64 seek = 0;
                                integer3 thatcell;
                                CellGroup *thatcg;
                                while (seek<cglist.size())
                                    if (g->b.id == cglist[seek++].id) goto FoundIt;
                                // We didn't find the link in the list, but maybe we don't care.
                                thatcell = g->b.localcell();
                                if (thatcell.x==lastslab) goto IgnoreIt;
                                    // Cell is in the no-go slab
                                thatcg = LinkToCellGroup(g->b);
                                if (!(thatcg->is_open())) goto IgnoreIt;
                                    // CellGroup is already marked as deferred
                                    // We've been there before; don't return!.
                                cglist.push_back(g->b);  // Yes, we'll need to consider that CellGroup
                                FoundIt:
                                IgnoreIt:
                                continue;  // Just effectively a no-op
                            }
                        }  // Done with links from this CellGroup
                        thiscg->defer_group();    // Mark it as deferred
                        searching++;   // Ready for next CellGroup in the queue
                    } // Done processing this item in the CellGroup queue
                } // Done with this seed link
            } // Done with this cell 
        }  // Done with this pencil
    }  // Done with this split
    GFC->DeferGroups.Stop();
    STDLOG(1, "Done looking for Deferrals in slab %d\n", slab);
    return;
}

void GlobalGroupSlab::ClearDeferrals() {
    // We need to loop over the relevant slabs to clear the deferral flags from all CellGroups
    GFC->ClearDefer.Start();
    STDLOG(1,"Clearing Deferrals for slab %d\n", slab);
    for (int s=0; s<diam; s++) {
        int thisslab = GFC->WrapSlab(slab+s-slabbias);
        NUMA_FOR(j,0,cpd)
            PencilAccum<CellGroup> pencil = GFC->cellgroups[thisslab][j];
            for (int g=0; g<pencil._size; g++) pencil.data[g].clear_deferral();

            /*
            for (int k=0; k<zwidth; k++) { // Loop over cells
                CellPtr<CellGroup> cg = GFC->cellgroups[thisslab][j][k];
                for (int g=0; cg.size(); g++) cg[g].clear_deferral();
            }
            */
        }
    }
    STDLOG(1,"Done Deferrals for slab %d\n", slab);
    GFC->ClearDefer.Stop();
    return;
}

void GlobalGroupSlab::CreateGlobalGroups() {
    // Given pstat[0,cpd)
    assertf(pstat!=NULL,"setup() not yet called\n");   // Setup not yet called!
    IndexLinks();
    #ifdef ONE_SIDED_GROUP_FINDING
    DeferGlobalGroups();
    #endif
    int max_group_diameter = 0;

    GFC->FindGlobalGroupTime.Start();
    // Then need to go cell by cell to take unclosed groups and try to close them.
    // That step will need to be run every (2R+1) pencil to avoid from contention.
    // 4R+1 might be better, just to avoid any cache overlap.
    // Have to guard against contention across the periodic wrap, as well.
    // Pick a division no less than diam that divides into CPD so that 
    // we can skip this issue.
    MultiplicityStats L0stats[omp_get_max_threads()];

    for (int w=0; w<split; w++) {
        #pragma omp parallel for schedule(dynamic,1) reduction(max:max_group_diameter)
        for (int j=w; j<cpd; j+=split) {
            // Do this pencil
            int cumulative_np = 0;    // We count all particles in this pencil
            PencilAccum<GlobalGroup> *gg_pencil = globalgroups.StartPencil(j);
            PencilAccum<LinkID> *gg_list = globalgrouplist.StartPencil(j);
            std::vector<LinkID> cglist;    // List of CellGroups in this GlobalGroup
            cglist.reserve(64);
            for (int k=0; k<zwidth; k++) {  // local k
                // Do this cell
                CellPtr<CellGroup> c = GFC->cellgroups[slab][j][k];
                for (int g=0; g<c.size(); g++) if (c[g].is_open()) {
                    // Loop over groups, skipping closed and deferred ones
                    int ggsize = 0;
                    int minzslab = k;  // to determine 2D ownership

                    cglist.resize(0);   // Reset to null list
                    cglist.push_back(LinkID(slab, j, k, g));   // Prime the queue
                    uint64 searching = 0;
                    // printf("GG %d %d %d %d\n", slab, j, k, g);
                    while (searching<cglist.size()) {
                        // We have an unresolved group to search
                        integer3 thiscell = cglist[searching].localcell();  
                        CellGroup *thiscg = LinkToCellGroup(cglist[searching]);
                        ggsize += thiscg->size();
                        minzslab = std::min(minzslab, thiscell.z);
                        assertf(thiscg->is_open(), "Cellgroup in slab %d found to be closed while making groups in slab %d.  Likely need to increase GroupRadius!\n", thiscell.x, slab);
                                // If this fails, probably a group has spanned
                                // beyond 2*R+1 cells and something got closed
                                // prematurely.
                        thiscg->close_group();
                        // if (searching>0) {
                            // printf("Link to %d %d %d %d\n",
                                    // thiscell.x, thiscell.y, thiscell.z, cglist[searching].cellgroup());
                        // }
                        // Now get these links
                        int s = GFC->WrapSlab(thiscell.x-slab+slabbias);  // Map to [0,diam)
                        assertf(s < diam, "GroupLink points to slab offset %d which violates GroupDiameter %d.  Likely need to increase GroupRadius! (slab,j,k,g = %d,%d,%d,%d)\n",
                            s, diam, slab, j, k, g);
                        LinkPencil *lp = links[s]+thiscell.y;
                        // Keep track of the maximum slab accessed
                        max_group_diameter = std::max(max_group_diameter, s);
                        GroupLink *gl; int t;
                        for (gl = lp->data + lp->cells[thiscell.z].start, t = 0;
                                t<lp->cells[thiscell.z].n; t++, gl++) {
                            // Consider all the links from this cell
                            if (gl->a.id == cglist[searching].id) {
                                // This link originates from the group we're using
                                // But now we need to check if its destination
                                // is already on the list.
                                uint64 seek = 0;
                                while (seek<cglist.size()) 
                                    if (gl->b.id == cglist[seek++].id) goto FoundIt;
                                cglist.push_back(gl->b);  // Didn't find it; add to queue
                                FoundIt:
                                // Either way, this link has been used its one time.
                                gl->mark_for_deletion();
                            }
                        } // Done with the links from this cell
                        searching++;  // Ready for the next unresolved cellgroup
                    } // End search over graph of cell groups
                    // printf("Closed with size %d\n", ggsize);

                    // recall that we are using local z
                    int minz_in_primary = minzslab >= ghost_radius && minzslab < zwidth - ghost_radius;

                    // Time to build the GlobalGroup
                    // We only track the group if it has more than one particle
                    if (ggsize>1 && !(dedup && !minz_in_primary)) {
                        int start = gg_list->buffer->get_pencil_size();
                        // We're going to sort the cellgroup list, so that
                        // multiple groups within one cell are contiguous
                        // But there's no point in doing this unless there are 3+ CG.
                        if (cglist.size()>2) {
                            std::sort(cglist.data(), cglist.data()+cglist.size());
                            /*
                            for (uint64 t = 1; t<cglist.size(); t++) 
                                assertf(cglist[t-1].id <= cglist[t].id,
                                    "Failed to sort propertly: %lld > %lld\n", 
                                        cglist[t-1].id,
                                        cglist[t].id);
                            */
                        }

                        for (uint64 t = 0; t<cglist.size(); t++) {
                            //integer3 tmp = cglist[t].localcell();
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
    GFC->max_group_diameter = std::max(GFC->max_group_diameter, max_group_diameter);

    // Free the indexing space
    free(cells);
    free(links[0]);
    delete[] links;
    GFC->FindGlobalGroupTime.Stop();

    #ifdef ONE_SIDED_GROUP_FINDING
    ClearDeferrals();
    #endif
    return;
}



/** Copy all of the particles from the disparate CellGroups into a single list,
so that each GlobalGroup is contiguous.  This also establishes an index
from each GlobalGroup into this slab-scale list.
*/

void GlobalGroupSlab::GatherGlobalGroups() {
    GFC->IndexGroups.Start();
    assertf(cpd>=4*GFC->GroupRadius+1, "CPD is too small compared to GroupRadius\n");
    assertf(zwidth>=4*GFC->GroupRadius+1, "zwidth=%d is too small compared to GroupRadius\n", zwidth);
    // This registers the periodic wrap using the cells.
    // However, this will misbehave if CPD is smaller than the group diameter,
    // because the LinkIDs have already been wrapped.
    // One just can't use a small CPD
    uint64 *this_pencil = new uint64[cpd];
    uint64 *this_pencil_cg = new uint64[cpd];
    int *largest = new int[cpd];

    int nthread = omp_get_max_threads();
    pint64 local_ncellgroups[nthread];
    for(int i = 0; i < nthread; i++)
        local_ncellgroups[i] = 0;

    NUMA_FOR(j,0,cpd)
        int g = omp_get_thread_num();
        uint64 local_this_pencil = 0;
        uint64 _local_ncellgroups = 0;
        int local_largest = 0;
        for (int k=0; k<zwidth; k++){  // local k
            for (int n=0; n<globalgroups[j][k].size(); n++) {
                int size = globalgroups[j][k][n].np;
                _local_ncellgroups += globalgroups[j][k][n].ncellgroups;
                local_this_pencil += size;
                local_largest = std::max(local_largest, size);
            }
        }
        this_pencil[j] = local_this_pencil;
        this_pencil_cg[j] = _local_ncellgroups;
        local_ncellgroups[g] += _local_ncellgroups;
        largest[j] = local_largest;
    }

    ncellgroups = 0;
    for(int i = 0; i < nthread; i++)
        ncellgroups += local_ncellgroups[i];

    // Now for the serial accumulation over the pencils
    uint64 total_particles = 0;
    uint64 total_cg = 0;
    largest_group = 0;
    for (int j=0; j<cpd; j++) {
        pstart[j] = total_particles;    // So we have the starting indices
        pstart_cg[j] = total_cg;        // So we have the starting indices
        total_particles += this_pencil[j];
        total_cg += this_pencil_cg[j];
        largest_group = std::max(largest[j], largest_group);
    }
    delete[] largest;
    delete[] this_pencil;
    delete[] this_pencil_cg;

    // Now we can allocate these buffers
    allocate(total_particles);
    GFC->IndexGroups.Stop();

    GFC->GatherGroups.Start();
    // Now copy the particles into these structures
    NUMA_FOR(j,0,cpd)
        for (int k=0; k<zwidth; k++){  // local k
            for (int n=0; n<globalgroups[j][k].size(); n++) {
                // Process globalgroups[j][k][n]
                // Compute where we'll put the particles, and update this starting point
                globalgroups[j][k].ptr(n)->start += pstart[j];
                uint64 start = globalgroups[j][k][n].start;

                LinkID *cglink = globalgrouplist.pencils[j].data
                                    +globalgroups[j][k][n].cellgroupstart;
                    // This is where we'll find the CG LinkIDs for this GG
                integer3 firstcell(slab,j,k + zstart);

                for (int c=0; c<globalgroups[j][k][n].ncellgroups; c++, cglink++) {
                    // Loop over CellGroups
                    integer3 cellijk = cglink->localcell();
                    cellijk.z += zstart;
                    CellGroup *cg = LinkToCellGroup(*cglink);
                    // printf("%d %d %d %d in %d %d %d n=%d\n", 
                        // j, k, n, c, cellijk.x, cellijk.y, cellijk.z, cg->size());
                    // CellGroup cg = GFC->cellgroups[cellijk.x][cellijk.y][cellijk.z][cglink->cellgroup()];
                    Cell cell = CP->GetCell(cellijk);
                    // Copy the particles, including changing positions to 
                    // the cell-centered coord of the first cell.  
                    // Note periodic wrap.
                    memcpy(vel+start, cell.vel+cg->start, sizeof(velstruct)*cg->size());
                    memcpy(aux+start, cell.aux+cg->start, sizeof(auxstruct)*cg->size());
                    for (int p=0; p<cg->size(); p++) aux[p+start].set_L0();
                            // This particle is in L0
                    if(cell.acc != NULL)
                        memcpy(acc+start, cell.acc+cg->start, sizeof(accstruct)*cg->size());
                    cellijk -= firstcell;
                    if (cellijk.x> diam) cellijk.x-=cpd;
                    if (cellijk.x<-diam) cellijk.x+=cpd;
                    if (cellijk.y> diam) cellijk.y-=cpd;
                    if (cellijk.y<-diam) cellijk.y+=cpd;
                    if (cellijk.z> diam) cellijk.z-=cpd;
                    if (cellijk.z<-diam) cellijk.z+=cpd;
                    posstruct offset = GFC->invcpd*(cellijk);
                        // Ok to use single precision because this is only a few cells 
                        // (and L1 group finding is in limited precision too)
                    // printf("Using offset %f %f %f\n", offset.x, offset.y, offset.z);
                    for (int p=0; p<cg->size(); p++) pos[start+p] = offset+cell.pos[cg->start+p]; 
                    start += cg->size();
                } // End loop over cellgroups in this global group

                // TODO: Might need to compute the COM for light cone output
            } // End loop over globalgroups in a cell
        }
    }
    // End loop over cells
    GFC->GatherGroups.Stop();
    return;
}




/** This writes the information from the GlobalGroup slab-scale list
of particles back into the original cell-based lists.  Here, we write
back only the AuxSlab information.
*/

void GlobalGroupSlab::ScatterGlobalGroupsAux() {
    GFC->ScatterAux.Start();

    NUMA_FOR(j,0,cpd)
        for (int k=0; k<zwidth; k++){  // local k
            for (int n=0; n<globalgroups[j][k].size(); n++) {
                // Process globalgroups[j][k][n]
                // Recall where the particles start
                uint64 start = globalgroups[j][k][n].start;

                LinkID *cglink = globalgrouplist.pencils[j].data
                                    +globalgroups[j][k][n].cellgroupstart;
                    // This is where we'll find the CG LinkIDs for this GG
                for (int c=0; c<globalgroups[j][k][n].ncellgroups; c++, cglink++) {
                    // Loop over CellGroups
                    integer3 cellijk = cglink->localcell();
                    cellijk.z += zstart;
                    CellGroup *cg = LinkToCellGroup(*cglink);
                    Cell cell = CP->GetCell(cellijk);
                    // Copy the aux back
                    memcpy(cell.aux+cg->start, aux+start, sizeof(auxstruct)*cg->size());
                    start += cg->size();
                } // End loop over cellgroups in this global group
            } // End loop over globalgroups in a cell
        }
    }
    // End loop over cells
    GFC->ScatterAux.Stop();
    return;
}




/** This writes the information from the GlobalGroup slab-scale list
of particles back into the original cell-based lists.  Here, we write
back only the PosSlab and VelSlab information.  AuxSlab is separate,
as some uses don't need Pos/Vel to be returned.
*/

void GlobalGroupSlab::ScatterGlobalGroups() {
    GFC->ScatterGroups.Start();

    STDLOG(1,"Scattering global group pos/vel from slab %d\n", slab);
    NUMA_FOR(j,0,cpd)
        for (int k=0; k<zwidth; k++){  // local k
            for (int n=0; n<globalgroups[j][k].size(); n++) {
                // Process globalgroups[j][k][n]
                // Recall where the particles start
                uint64 start = globalgroups[j][k][n].start;

                LinkID *cglink = globalgrouplist.pencils[j].data
                                    +globalgroups[j][k][n].cellgroupstart;
                    // This is where we'll find the CG LinkIDs for this GG
                integer3 firstcell(slab,j,k + zstart);
                for (int c=0; c<globalgroups[j][k][n].ncellgroups; c++, cglink++) {
                    // Loop over CellGroups
                    integer3 cellijk = cglink->localcell();
                    cellijk.z += zstart;
                    CellGroup *cg = LinkToCellGroup(*cglink);
                    Cell cell = CP->GetCell(cellijk);
                    // Copy the particles, including changing positions to 
                    // the cell-centered coord of the first cell.  
                    // Note periodic wrap.
                    memcpy(cell.vel+cg->start, vel+start, sizeof(velstruct)*cg->size());
                    //memcpy(cell.acc+cg->start, acc+start, sizeof(accstruct)*cg->size());
                    cellijk -= firstcell;
                    if (cellijk.x> diam) cellijk.x-=cpd;
                    if (cellijk.x<-diam) cellijk.x+=cpd;
                    if (cellijk.y> diam) cellijk.y-=cpd;
                    if (cellijk.y<-diam) cellijk.y+=cpd;
                    if (cellijk.z> diam) cellijk.z-=cpd;
                    if (cellijk.z<-diam) cellijk.z+=cpd;
                    posstruct offset = GFC->invcpd*(cellijk);
                        // Ok to use single precision, because this is only a few cells
                    // printf("Using offset %f %f %f\n", offset.x, offset.y, offset.z);
                    for (int p=0; p<cg->size(); p++) 
                         cell.pos[cg->start+p] = pos[start+p] - offset;
                    start += cg->size();
                } // End loop over cellgroups in this global group
            } // End loop over globalgroups in a cell
        }
    }
    // End loop over cells
    GFC->ScatterGroups.Stop();
    STDLOG(1,"Done scattering global group pos/vel from slab %d\n", slab);
    return;
}


/** This considers each GlobalGroup (L0) and performs addiitonal
work to find the L1 and L2 halos and generate output statistics
about L1.  This can either use FOF or SO.
*/

void GlobalGroupSlab::AppendParticleToPencil(PencilAccum<RVfloat> ** pHaloRVs, PencilAccum<TaggedPID> ** pHaloPIDs, 
                                            posstruct * grouppos, velstruct * groupvel, accstruct * groupacc, auxstruct * groupaux, int index, posstruct offset, int & nA, int & nB) {
    int taggable = groupaux[index].is_taggable();
    int npA = 0; int npB = 0;

    if (taggable != 0 || P.OutputAllHaloParticles) {
        posstruct r = WrapPosition(grouppos[index]+offset);
        velstruct v = groupvel[index];
        // Velocities were full kicked; half-unkick before halostats
        if (groupacc != NULL)
            v -= TOFLOAT3(groupacc[index])*WriteState.FirstHalfEtaKick;
        v *= ReadState.VelZSpace_to_kms/ReadState.VelZSpace_to_Canonical; 
        
        if ((taggable == TAGGABLE_SUB_A) || P.OutputAllHaloParticles){ // if we request all halo particles, output them in subsample 'A' files.
            pHaloPIDs[0]->append(TaggedPID(groupaux[index]));
            pHaloRVs[0] ->append(RVfloat(r.x, r.y, r.z, v.x, v.y, v.z));
            npA ++; 

        }
        else if (taggable == TAGGABLE_SUB_B){
            pHaloPIDs[1]->append(TaggedPID(groupaux[index]));
            pHaloRVs[1] ->append(RVfloat(r.x, r.y, r.z, v.x, v.y, v.z));
            npB ++; 
        }
    }

    nA += npA; nB += npB;
}


void GlobalGroupSlab::FindSubGroups() {
    GFC->ProcessLevel1.Start();
    int maxthreads = omp_get_max_threads();

    // Guessing L1 halo count is 3% of particles, given minsize
    L1halos.setup(cpd, zwidth, GFC->particles_per_slab/30);    
    // Guessing L1 Particles are 3% of total (taggable)
    HaloRVA.setup(cpd, zwidth, GFC->particles_per_slab/30);   
    HaloRVB.setup(cpd, zwidth, GFC->particles_per_slab/30);  
     
    HaloPIDsA.setup(cpd, zwidth, GFC->particles_per_slab/30);    
    HaloPIDsB.setup(cpd, zwidth, GFC->particles_per_slab/30);    
    posstruct **L1pos = new posstruct *[maxthreads];
    velstruct **L1vel = new velstruct *[maxthreads];
    auxstruct **L1aux = new auxstruct *[maxthreads];
    accstruct **L1acc = new accstruct *[maxthreads];
    MultiplicityStats L1stats[maxthreads];

    // BOT
    MultiplicityStats L2stats[maxthreads];

    #ifdef SPHERICAL_OVERDENSITY
	SOcell FOFlevel1[maxthreads], FOFlevel2[maxthreads];
	#pragma omp parallel for schedule(static,1)
	for (int g=0; g<maxthreads; g++) {
	    FOFlevel1[g].setup(GFC->SOdensity1, P.SO_NPForMinDensity);
	    FOFlevel2[g].setup(GFC->SOdensity2, P.SO_NPForMinDensity*GFC->SOdensity2/GFC->SOdensity1);
	}
	STDLOG(1,"Seeking SO halos, L1 = %f, L2 = %f, with min_central = %f and %f\n", 
		FOFlevel1[0].threshold, FOFlevel2[0].threshold,
        FOFlevel1[0].min_central/FOFlevel1[0].FOFunitdensity, 
        FOFlevel2[0].min_central/FOFlevel2[0].FOFunitdensity);
    #else
	FOFcell FOFlevel1[maxthreads], FOFlevel2[maxthreads];
	#pragma omp parallel for schedule(static,1)
	for (int g=0; g<maxthreads; g++) {
	    FOFlevel1[g].setup(GFC->linking_length_level1, 1e10);
	    FOFlevel2[g].setup(GFC->linking_length_level2, 1e10);
	}
	STDLOG(1,"Seeking FOF halos, L1 = %f, L2 = %f (comoving unit-box units)\n", 
		FOFlevel1[0].linking_length, FOFlevel2[0].linking_length);
    #endif

    #pragma omp parallel for schedule(static,1)
    for (int g=0; g<maxthreads; g++) {
        L1pos[g] = new posstruct[largest_group];
        L1vel[g] = new velstruct[largest_group];
        L1aux[g] = new auxstruct[largest_group];
        L1acc[g] = new accstruct[largest_group];
    }

    // It seems that the work between pencils is so heterogeneous that even the
    // dynamic scheduling can't smooth it out.  So we're going to try ordering the
    // pencils by the work estimate (largest first)
    std::sort(pstat, pstat+cpd);
    
    int nthread = omp_get_max_threads();
    padded<int> local_np_subA[nthread];
    padded<int> local_np_subB[nthread];
    for(int i = 0; i < nthread; i++){
        local_np_subA[i] = 0;
        local_np_subB[i] = 0;
    }


    NUMA_FOR(jj,0,cpd)
        int j = pstat[jj].pnum;    // Get the pencil number from the list
        GFC->L1Tot.Start();
        const int g = omp_get_thread_num();
        PencilAccum<HaloStat>  *pL1halos   =   L1halos.StartPencil(j);
        PencilAccum<RVfloat>   *pHaloRVA   =   HaloRVA.StartPencil(j);
        PencilAccum<RVfloat>   *pHaloRVB   =   HaloRVB.StartPencil(j);
        PencilAccum<TaggedPID> *pHaloPIDsA = HaloPIDsA.StartPencil(j);
        PencilAccum<TaggedPID> *pHaloPIDsB = HaloPIDsB.StartPencil(j);

        PencilAccum<TaggedPID> * pHaloPIDs[NUM_SUBSAMPLES] = {pHaloPIDsA, pHaloPIDsB}; 
        PencilAccum<RVfloat>   *  pHaloRVs[NUM_SUBSAMPLES] = {pHaloRVA,   pHaloRVB  }; 


        for (int k=0; k<zwidth; k++){  // local k
            // uint64 groupid = ((slab*cpd+j)*cpd+k)*4096;
            uint64 groupidstart = (((uint64)slab*10000+(uint64)j)*10000+(uint64) (k + zstart))*1000000;
            uint64 nextid; 
            posstruct offset = CP->CellCenter(slab, j, k + zstart);
                    // A basic label for this group
            for (int n=0; n<globalgroups[j][k].size(); n++) { //loop over L0 groups. 
                nextid = groupidstart + n*1000;
                posstruct *grouppos = pos+globalgroups[j][k][n].start;
                velstruct *groupvel = vel+globalgroups[j][k][n].start;
                auxstruct *groupaux = aux+globalgroups[j][k][n].start;
                accstruct *groupacc = NULL;
                if(acc != NULL) groupacc = acc+globalgroups[j][k][n].start;
                int groupn = globalgroups[j][k][n].np;

                if (groupn >= GFC->minhalosize){ //This L0 group is too small. Need to output its particles if we're doing subsampling, though. 

                    GFC->L1FOF.Start();
                    #ifndef SPHERICAL_OVERDENSITY
                    if (GFC->linking_length_level1==GFC->linking_length)
                        FOFlevel1[g].assign_to_one_group(grouppos, NULL, NULL, groupacc, groupn);
                    else    // This grabs the findgroups statement!
                    #endif
                    FOFlevel1[g].findgroups(grouppos, NULL, groupaux, groupacc, groupn);
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
                            if (groupacc != NULL){
                                L1acc[g][b] = groupacc[start[b].index()];
                                // Velocities were full kicked; half-unkick before halostats
                                L1vel[g][b] -= TOFLOAT3(L1acc[g][b])*WriteState.FirstHalfEtaKick;
                            }
                        }

                        GFC->L2FOF.Start();
                        FOFlevel2[g].findgroups(L1pos[g], NULL, L1aux[g], L1acc[g], size);
                        GFC->L2FOF.Stop();

                        // BOT
                        for (int b=0; b<FOFlevel2[g].ngroups; b++) {
                            int size_L2 = FOFlevel2[g].groups[b].n;
                            if (size_L2<GFC->minhalosize) continue;
                            L2stats[g].push(size_L2);
                        }
                    
                        // Merger trees require tagging the taggable particles 
                        // of the biggest L2 group in the original aux array.  
                        // This can be done:
                        // The L2 index() gives the position in the L1 array,
                        // and that index() gets back to aux.
                        uint32_t ntaggedA = 0;
                        uint32_t ntaggedB = 0;
                        if(FOFlevel2[g].ngroups > 0){
                            std::sort(FOFlevel2[g].groups, FOFlevel2[g].groups+FOFlevel2[g].ngroups);
                            // Groups now in descending order of multiplicity
                    
                            FOFparticle *L2start = FOFlevel2[g].p + FOFlevel2[g].groups[0].start;
                            for (int p=0; p<FOFlevel2[g].groups[0].n; p++) {
                                if (groupaux[start[L2start[p].index()].index()].is_taggable()) //1 or 2 = taggable. 0 = not taggable. 
                                    groupaux[start[L2start[p].index()].index()].set_tagged();                
                            }
                        }

                        uint64 npstartA    = pHaloRVA->get_pencil_size();
                        uint64 npstartB    = pHaloRVB->get_pencil_size();

                        // Output the Tagged Particles and Taggable Particles. 
                        // If P.OutputAllHaloParticles is set, then we output 
                        //      all L1 particles
                        for (int b=0; b<size; b++){
                            int index = start[b].index(); 

                            if (groupaux[index].is_tagged()) {
                                int taggable = groupaux[index].is_taggable();
                                assertf(taggable != 0, "Uh oh, this particle is tagged but not taggable\n");
                                if      (taggable == TAGGABLE_SUB_A) ntaggedA++;
                                else if (taggable == TAGGABLE_SUB_B) ntaggedB++;
                            }

                            AppendParticleToPencil(pHaloRVs, pHaloPIDs, grouppos, groupvel, groupacc, groupaux, index, offset, local_np_subA[g].i, local_np_subB[g].i);
                        }

                        HaloStat h = ComputeStats(size, L1pos[g], L1vel[g], L1aux[g], FOFlevel2[g], offset);
                        h.id = nextid++;
                        h.L0_N = groupn;
                        h.npstartA = npstartA;
                        h.npstartB = npstartB;
                        h.ntaggedA = ntaggedA;  
                        h.ntaggedB = ntaggedB;


                        #ifdef SPHERICAL_OVERDENSITY
                        //fetch SO stats for this L1 halo.
                        int idx = FOFlevel1[g].groups[a].center_particle;  // Index of the central particle in the L0 set
                        posstruct central_particle = WrapPosition(grouppos[idx] + offset);
                        h.SO_central_particle[0] = central_particle.x; 
                        h.SO_central_particle[1] = central_particle.y; 
                        h.SO_central_particle[2] = central_particle.z; 
                        h.SO_central_density  = (GFC->use_aux_dens ? groupaux[idx].get_density() : groupacc[idx].w) / FOFlevel1[g].FOFunitdensity; 
                        h.SO_radius = sqrt(FOFlevel1[g].groups[a].halo_radius2) / FOF_RESCALE; 

                        //now repeat for the largest L2 halo.
                        idx = FOFlevel2[g].groups[0].center_particle;      // Index of the central particle in the L1 set
                        posstruct L2max_central_particle = WrapPosition(L1pos[g][idx] + offset);
                        h.SO_L2max_central_particle[0] = L2max_central_particle.x;
                        h.SO_L2max_central_particle[1] = L2max_central_particle.y;
                        h.SO_L2max_central_particle[2] = L2max_central_particle.z;
                        h.SO_L2max_central_density  =    (GFC->use_aux_dens ? L1aux[g][idx].get_density() : L1acc[g][idx].w) / FOFlevel2[g].FOFunitdensity; 
                        h.SO_L2max_radius = sqrt(FOFlevel2[g].groups[0].halo_radius2) / FOF_RESCALE; 

                        #else
                        //set SO stats to zero.
                        for (int i = 0; i < 3; i ++){
                            h.SO_central_particle[i]       = 0.0; 
                            h.SO_L2max_central_particle[i] = 0.0;
                        }
                        h.SO_central_density = 0.0; h.SO_L2max_central_density = 0.0; 
                        h.SO_radius          = 0.0; h.SO_L2max_radius          = 0.0; 
                        #endif 

                        h.npoutA = pHaloRVA->get_pencil_size()-npstartA;
                        h.npoutB = pHaloRVB->get_pencil_size()-npstartB;
                   
                        pL1halos->append(h);
                    } // Done with this L1 halo
                    // If this is a subsampling redshift, 
                    // we now output the non-L1 particles in this L0 group.
                    // This must happen here because the microstepping might 
                    // alter the L0 pos/vel before we output the field particles.
                }

                if (ReadState.DoSubsampleOutput){ //Regardless of whether this L0 group is big enough to do L1 group finding, 
                                                    //if we're outputing the particle subsample, output all of its L0 particles. 
                    for (int b=0; b<groupn; b++) {
                        if (groupaux[b].is_L1()) continue;  // Already in the L1 set
                        AppendParticleToPencil(pHaloRVs, pHaloPIDs, grouppos, groupvel, groupacc, groupaux, b, offset, local_np_subA[g].i, local_np_subB[g].i); 
                    }
                }
                
                // If we want to do anything else with the L0 particles,
                // we could do so here.
            } // Done with this L0 group
            pL1halos->FinishCell();
            pHaloRVA->FinishCell();
            pHaloRVB->FinishCell();
            pHaloPIDsA->FinishCell();
            pHaloPIDsB->FinishCell();

        }
        pL1halos->FinishPencil();
        pHaloRVA->FinishPencil();
        pHaloRVB->FinishPencil();
        pHaloPIDsA->FinishPencil();
        pHaloPIDsB->FinishPencil();
        GFC->L1Tot.Stop();
    }

    int np_subA = 0, np_subB = 0;
    for(int i = 0; i < nthread; i++){
        np_subA += local_np_subA[i];
        np_subB += local_np_subB[i];
    }

    // Need to update the pL1halos.npstart values for their pencil starts!
    HaloRVA.build_pstart();
    HaloRVB.build_pstart();
    for (int j=0; j<cpd; j++) 
        for (int k=0; k<zwidth; k++)  // local k
            for (int n=0; n<L1halos[j][k].size(); n++) {
                HaloStat *h = L1halos[j][k].ptr(n);
                h->npstartA += HaloRVA.pstart[j];
                h->npstartB += HaloRVB.pstart[j];
            }

    // Coadd the stats
    if (ReadState.DoSubsampleOutput){
        WriteState.np_subA_state += np_subA; 
        WriteState.np_subB_state += np_subB; 
    }

    uint64 previous = GFC->L1stats.ngroups;
    // BOT
    uint64 previous_L2 = GFC->L2stats.ngroups;
    for (int g=0; g<omp_get_max_threads(); g++) {
        GFC->L1stats.add(L1stats[g]);

        // BOT
        GFC->L2stats.add(L2stats[g]);
    	GFC->numdists1 += FOFlevel1[g].numdists;
    	GFC->numdists2 += FOFlevel2[g].numdists;
    	GFC->numsorts1 += FOFlevel1[g].numsorts;
    	GFC->numsorts2 += FOFlevel2[g].numsorts;
    	GFC->numcenters1 += FOFlevel1[g].numcenters;
    	GFC->numcenters2 += FOFlevel2[g].numcenters;
    	#ifdef SPHERICAL_OVERDENSITY
    	GFC->numcg1 += FOFlevel1[g].numcg;
    	GFC->numcg2 += FOFlevel2[g].numcg;
    	GFC->numgroups1 += FOFlevel1[g].numgroups;
    	GFC->numgroups2 += FOFlevel2[g].numgroups;
    	if (g>0) {
    	    FOFlevel1[0].coadd_timers(FOFlevel1[g]);
    	    FOFlevel2[0].coadd_timers(FOFlevel2[g]);
        }
        #endif
    }
    previous = GFC->L1stats.ngroups-previous;
    STDLOG(1,"Found %l L1 halos\n", previous);
    // BOT
    previous_L2 = GFC->L2stats.ngroups-previous_L2;
    STDLOG(1,"Found %l L2 halos\n", previous_L2);
    
    #ifdef SPHERICAL_OVERDENSITY
    if (FOFlevel1[0].Total.Elapsed()>0.0) {
      // B.H. trying to do timing
      //STDLOG(3,"L1 Timing: %f = %f %f %f %f\n",
      STDLOG(1,"L1 Timing: %f = %f %f %f %f %f\n",
	    FOFlevel1[0].Total.Elapsed(),
	    FOFlevel1[0].Copy.Elapsed(),
	    FOFlevel1[0].Sweep.Elapsed(),
	    FOFlevel1[0].Distance.Elapsed(),
	    FOFlevel1[0].Search.Elapsed(),
	    FOFlevel1[0].Sort.Elapsed());
	STDLOG(1,"L2 Timing: %f = %f %f %f %f %f\n",
	    FOFlevel2[0].Total.Elapsed(),
	    FOFlevel2[0].Copy.Elapsed(),
	    FOFlevel2[0].Sweep.Elapsed(),
	    FOFlevel2[0].Distance.Elapsed(),
	    FOFlevel2[0].Search.Elapsed(),
	    FOFlevel2[0].Sort.Elapsed());
	// These show that the timing is 50% dominated by the Search step,
	// with most of the rest split between Copy and Sweep.
    }
    #endif

    // Now delete all of the temporary storage!
    // Want free's to be on the same threads as the original
    #pragma omp parallel for schedule(static,1)
    for (int g=0; g<omp_get_max_threads(); g++) {
        FOFlevel1[g].destroy();
        FOFlevel2[g].destroy();
        delete[] L1pos[g];
        delete[] L1vel[g];
        delete[] L1aux[g];
        delete[] L1acc[g];
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
    for (int j=0; j<cpd; j++)
        for (int k=0; k<zwidth; k++)  // local k
            for (int n=0; n<globalgroups[j][k].size(); n++)
                fprintf(fp, "%d %d %d %d %d\n", (int)globalgroups[j][k][n].start,
                                globalgroups[j][k][n].np, slab, j, k);
    fclose(fp);

    sprintf(fname, "/tmp/out.halo.%03d", slab);
    fp = fopen(fname,"w");
    for (int j=0; j<cpd; j++)
        for (int k=0; k<zwidth; k++)  // local k
            for (int n=0; n<L1halos[j][k].size(); n++) {
                HaloStat h = L1halos[j][k][n];
                fprintf(fp, "%4d %7.4f %7.4f %7.4f %d %4d %3d %3d %7.4f %7.4f %7.4f %d %lu %u %u\n", 
                    h.N, h.x_com[0], h.x_com[1], h.x_com[2], h.r50_com,
                    h.L2_N[0], h.L2_N[1], h.L2_N[2], 
                    h.x_L2com[0], h.x_L2com[1], h.x_L2com[2], 
                    h.r50_L2com, h.id, h.npoutA, h.npoutB);
            }
    fclose(fp);

    // Remember that these positions are relative to the first-cell position,
    // which is why we include that cell ijk in the first outputs
    sprintf(fname, "/tmp/out.pos.%03d", slab);
    fp = fopen(fname,"w");
    for (uint64 p=0; p<np; p++)
        fprintf(fp, "%f %f %f %d\n", pos[p].x, pos[p].y, pos[p].z, (int)aux[p].pid());
    fclose(fp);
}

#if 0
void GlobalGroupSlab::HaloOutput() {
    char fname[200];
    sprintf(fname, "/tmp/out.binhalo.%03d", slab);
    L1halos.dump_to_file(fname);

    sprintf(fname, "/tmp/out.L1rv.%03d", slab);
    HaloRV.dump_to_file(fname);

    sprintf(fname, "/tmp/out.L1pid.%03d", slab);
    HaloPIDs.dump_to_file(fname);
}

#else   // !STANDALONE_FOF

/** Outputs a uniform subsample of the particles ("taggable particles")
and L1 halo stats and tagged particles. Currently we always output
taggable particles if we are also outputting halos. Halo outputs will
be skipped if no L1 halos were found (but taggable particles will still be written).
 */


void GlobalGroupSlab::HaloOutput() {
    GFC->OutputLevel1.Start();
    STDLOG(0,"Beginning halo output for slab %d\n", slab);

    // This will create the directory if it doesn't exist (and is parallel safe)
    char dir[32];
    sprintf(dir, "Step%04d_z%5.3f", ReadState.FullStepNumber, ReadState.Redshift);
    CreateSubDirectory(P.GroupDirectory, dir);

    std::string headerfn = "";
    headerfn = headerfn + P.GroupDirectory + "/" + dir + "/header";
    WriteGroupHeaderFile(headerfn.c_str());

    if (ReadState.DoSubsampleOutput){
         // Write out the pos/vel of the taggable particles in L1 halos
        if (P.ParticleSubsampleA > 0){
            SB->AllocateSpecificSize(HaloRVSlabA, slab, HaloRVA.get_slab_bytes());
            HaloRVA.copy_to_ptr((RVfloat *)SB->GetSlabPtr(HaloRVSlabA, slab));
            SB->StoreArenaNonBlocking(HaloRVSlabA, slab);
        }
        if (P.ParticleSubsampleB > 0) {
            SB->AllocateSpecificSize(HaloRVSlabB, slab, HaloRVB.get_slab_bytes());
            HaloRVB.copy_to_ptr((RVfloat *)SB->GetSlabPtr(HaloRVSlabB, slab));
            SB->StoreArenaNonBlocking(HaloRVSlabB, slab);
        }
    }

    // Write out the PIDs of the taggable particles in the halos. If DoSubsampleOutput, store L0 and L1. If DoGrpFindingOutput only, do L1 only. 
    if (P.ParticleSubsampleA > 0 and HaloPIDsA.get_slab_size() != 0){
        SB->AllocateSpecificSize(HaloPIDsSlabA, slab, HaloPIDsA.get_slab_bytes());
        HaloPIDsA.copy_to_ptr((TaggedPID *)SB->GetSlabPtr(HaloPIDsSlabA, slab));
        SB->StoreArenaNonBlocking(HaloPIDsSlabA, slab);
    }

    if (P.ParticleSubsampleB > 0 and HaloPIDsB.get_slab_size() != 0) {
        SB->AllocateSpecificSize(HaloPIDsSlabB, slab, HaloPIDsB.get_slab_bytes());
        HaloPIDsB.copy_to_ptr((TaggedPID *)SB->GetSlabPtr(HaloPIDsSlabB, slab));
        SB->StoreArenaNonBlocking(HaloPIDsSlabB, slab);
    }

    // If we have catalogues to output, do so. 
    if (not (L1halos.pencils == NULL || L1halos.get_slab_size() == 0)){
        // Write out the stats on the L1 halos
        SB->AllocateSpecificSize(L1halosSlab, slab, L1halos.get_slab_bytes());
        L1halos.copy_to_ptr((HaloStat *)SB->GetSlabPtr(L1halosSlab, slab));
        SB->StoreArenaNonBlocking(L1halosSlab, slab);
    }

    GFC->OutputLevel1.Stop();
    return;
}

void GlobalGroupSlab::WriteGroupHeaderFile(const char* fn){
    std::ofstream headerfile;
    headerfile.open(fn);
    headerfile << P.header();
    headerfile << ReadState.header();
    headerfile << "\nOutputType = \"GroupOutput\"\n";
    headerfile.close();
}

#endif






#ifndef STANDALONE_FOF
/** When we do group finding, we write all the non-L0 particles into normal
time-slice slabs and the L0 group particles into their own slabs.
This routine does the latter.  We are passed full-kicked particles,
so we need to half-unkick.
The particles are in group-centered units, which we will convert to
global units.  For pack14, each cell group will get its own header.
*/

uint64 GlobalGroupSlab::L0TimeSliceOutput(FLOAT unkick_factor){

    AppendArena *AA;
    uint64 n_added = 0;

    AA = get_AA_by_format(P.OutputFormat);

    // Setup the Arena
    int headersize = 1024*1024;
    SB->AllocateSpecificSize(L0TimeSlice, slab, np*AA->sizeof_particle()
                                + ncellgroups*AA->sizeof_cell() + headersize);
    AA->initialize(L0TimeSlice, slab, CP->cpd, ReadState.VelZSpace_to_Canonical);

    // add the ParseHeader
    AA->addheader((const char *) P.header());
    AA->addheader((const char *) ReadState.header());
    char head[1024];
    sprintf(head, "\nOutputType = \"L0TimeSlice\"\n"); 
    AA->addheader((const char *) head);
    sprintf(head, "SlabNumber = %d\n", slab);
    AA->addheader((const char *) head);
    // For sanity, be careful that the previous lines end with a \n!
    AA->finalize_header();


    SlabAccum<TaggedPID> TimeSlicePIDs;
    TimeSlicePIDs.setup(cpd, zwidth, GFC->particles_per_slab);

    // Now scan through cell groups
    #pragma omp parallel for schedule(static) reduction(+:n_added)
    for (int j=0; j<cpd; j++){

        PencilAccum<TaggedPID> *pTimeSlicePIDs = TimeSlicePIDs.StartPencil(j);
        AA->start_pencil(j, pstart[j]*AA->sizeof_particle() + pstart_cg[j]*AA->sizeof_cell());

        for (int k = 0; k < zwidth; k++){  // local k
            // At this point, the only global groups are those that belong uniquely to this node,
            // so it's okay to scan the ghost cells (a group may start there)
            integer3 firstcell(slab,j,k + zstart);
            for (int n=0; n<globalgroups[j][k].size(); n++) {
                // Process globalgroups[j][k][n]
                // Recall where the particles start
                uint64 start = globalgroups[j][k][n].start;

                LinkID *cglink = globalgrouplist.pencils[j].data
                                    +globalgroups[j][k][n].cellgroupstart;
                    // This is where we'll find the CG LinkIDs for this GG

                for (int c=0; c<globalgroups[j][k][n].ncellgroups; c++, cglink++) {
                    // Loop over CellGroups
                    integer3 cellijk = cglink->localcell();
                    cellijk.z += zstart;
                    CellGroup *cg = LinkToCellGroup(*cglink);
                    Cell _cell = CP->GetCell(cellijk);
                    // Convert pos back to global. Note periodic wrap.
                    cellijk -= firstcell;
                    /// TODO: Is this logic still correct??
                    if (cellijk.x> diam) cellijk.x-=cpd;
                    if (cellijk.x<-diam) cellijk.x+=cpd;
                    if (cellijk.y> diam) cellijk.y-=cpd;
                    if (cellijk.y<-diam) cellijk.y+=cpd;
                    if (cellijk.z> diam) cellijk.z-=cpd;
                    if (cellijk.z<-diam) cellijk.z+=cpd;
                    posstruct offset = GFC->invcpd*cellijk;

                    // The max vel tracking happens in the kick, which operates on all particles, so it should be valid for L0 groups
                    FLOAT vscale = _cell.ci->max_component_velocity/ReadState.VelZSpace_to_Canonical;
                    // Start the cell group
                    AA->addcell(j, _cell.ijk, vscale);

                    for(int pi = start; pi < start + cg->size(); pi++){
                        posstruct _p = pos[pi] - offset;
                        velstruct _v = vel[pi] - TOFLOAT3(acc[pi])*unkick_factor;
                        AA->addparticle(j, _p, _v, aux[pi]);
                        if (strcmp(P.OutputFormat,"Pack9")==0) pTimeSlicePIDs->append(TaggedPID(aux[pi]));
                        n_added++;
                    }
                    AA->endcell(j);
                    start += cg->size();
                }
            }
            pTimeSlicePIDs->FinishCell();

        }

        pTimeSlicePIDs->FinishPencil();

    }

    if (strcmp(P.OutputFormat,"Pack9")==0){
        SB->AllocateSpecificSize(L0TimeSlicePIDs, slab, TimeSlicePIDs.get_slab_bytes());
        TimeSlicePIDs.copy_to_ptr((TaggedPID *)SB->GetSlabPtr(L0TimeSlicePIDs, slab));
        SB->StoreArenaNonBlocking(L0TimeSlicePIDs, slab);
    }

    SB->ResizeSlab(L0TimeSlice, slab, AA->finalize_arena());
    // Write out this time slice
    SB->StoreArenaNonBlocking(L0TimeSlice, slab);
    delete AA;

    return n_added;
}
#endif
