/* Initial code to support data transfer between multiple nodes */


// ================== Manifest helpers =====================================

class ManifestArena {
  // The information we're passing for a single arena
  public:
    uint64 size;    // Arena size in bytes of usable space, not including guards
    char *ptr;	// The pointer to the space
    int type;    // The SlabType
    int slab; 	// slab number
    ManifestArena() { size = 0; type = slab = 0; ptr = NULL; }
};

class DependencyRecord {
  // The information we're passing for a single dependency
  public:
    int begin;
    int end;
    // We're passing the information that [begin,end) should be set by force
    // on the receiving node.
    // We require begin<=end; the wrapping is done within the Dependency class.
    DependencyRecord() { begin = end = 0; }

    void Load(Dependency &d, int finished_slab) {
	// Load the status of this dependency 
	end = finished_slab;
	for (begin=end-1; begin>end-d.cpd; begin--) {
	    if (d.notdone(begin)) break;
	    d.mark_to_repeat(begin);
	}
	// We've found the first notdone slab
	begin++;   // Want to pass the first done one
	STDLOG(1, "Dependency [%d,%d)\n", begin, end);
	return;
    }
    void Set(Dependency &d) {
	for (int s=begin; s<end; s++) d.force_done(s);
	return;
    }

    // We also are going to use this same logic to hold the cellgroup_status[]
    // information from GFC.  In particular, we need to transmit all cases
    // where this is ==1.  We don't care about ==2.
    void LoadCG(int finished_slab) {
	if (GFC==NULL) { begin=end = finished_slab-1; return;}
	STDLOG(1, "LoadCG start\n");
	end = finished_slab-1;
	while (GFC->cellgroups_status[PP->WrapSlab(end)]==2) end--;
	end++;   // Now this marks the first ==2.
	for (begin=end-1; begin>end-GFC->cpd; begin--) {
	    int s = PP->WrapSlab(begin);
	    if (GFC->cellgroups_status[s]==0) break;
	    // Need to load this over to the Arenas
	    STDLOG(1,"Packing CellGroupArena for slab %d\n", s);
	    GFC->cellgroups[s].pack(CellGroupArena,s);
	    // Having moved it, we can delete the original
	    GFC->DestroyCellGroups(s);
	}
	begin++;   // Now this marks the first ==1.
	STDLOG(1, "CG Dependency [%d,%d)\n", begin, end);
	return;
    }
    void SetCG() {
	// Set the cellgroups_status to 1 for the indicated slabs
	for (int s=begin; s<end; s++) {
	    GFC->cellgroups_status[PP->WrapSlab(s)]=1;
	    // And move the information from the Arenas back into the SlabArray
	    GFC->cellgroups[PP->WrapSlab(s)].unpack(CellGroupArena,s);
	    // Now we can delete the CellGroupArena
	    LBW->DeAllocate(CellGroupArena,s);
	}
	return;
    }
};

// Need some new Partition functions, to get ranges of slabs
// Will select things between [x,slab), where x=slab-CPD/2, handling the wrap
inline bool is_below_slab(ilstruct *particle, int slab) {
    int x = particle->xyz.x-slab;  // We want x<0, but not too much.
    x = PP->WrapSlab(x);
    return x>PP->cpd/2;
}

// When we partition the GroupLink list, GlobalGroups have already
// been closed in slab, so there are no GroupLinks that involve the 
// slab.  By definition, the two ends of the link can differ by at 
// most one slab.  So it is enough to test one end.
inline bool link_below_slab(GroupLink *link, int slab) {
    int sa = link->a.slab();
    sa = PP->WrapSlab(sa-slab);
    return (sa>PP->cpd/2);
}

// ================== Manifest Class =====================================

// These are oversized, but still tiny compared to what will be transferred.
#define MAXMANIFEST 1024
#define MAXDEPENDENCY 64

class Manifest {
  // This is all of the information that we'll be passing between nodes
  public:
    ManifestArena arenas[MAXMANIFEST];
    int numarenas;   // The number of arenas being sent
    ilstruct *il;
    int numil;
    GroupLink *links;
    int numlinks;
    DependencyRecord dep[MAXDEPENDENCY];
    int numdep;
    int completed;	// ==1 if the Manifest is ready for use, ==0 if not yet, ==2 if already done

    Manifest() {
    	numarenas = numil = numlinks = numdep = 0;
	completed = 0;
	il = NULL;
	links = NULL;
	return;
    }
    ~Manifest() { }

    inline int is_ready() { if (completed==1) return 1; else return 0;}
	// Call this to see if the Manifest is ready to retrieve
    void done() { completed=2; }
    	// Call this to mark that we're done with the Manifest

    void LoadArena(int type, int s) {
	ManifestArena *a = arenas+numarenas;
	a->type = type;
	a->slab = s;
	a->size = LBW->IDSizeBytes(type, s);
	a->ptr =  LBW->ReturnIDPtr(type, s);
	STDLOG(1, "Queuing slab %d of type %d, size %l\n", s, type, a->size);
	numarenas++;
	assertf(numarenas<MAXMANIFEST, "numarenas has overflowed; increase MAXMANIFEST.");
	return;
    }

    // Here's the prototypes for the main routines
    void QueueToSend(int finished_slab);
    void Send();
    void Receive();
    void ImportData();
};    

void Manifest::QueueToSend(int finished_slab) {
    // We've just finished the given slab.  
    // Load all of the information from lower slabs into the Manifest.
    // Then call the non-blocking communication and return.
    int cpd = P.cpd;

    // TODO: Consider some safety measures here.
    // E.g., we really should not have any slabs that have open GlobalGroupSlabs
    // in our wake.  Might check, since one could screw up the Dependencies.

    STDLOG(1,"Queueing the SendManifest at slab=%d\n", finished_slab);

    // Load the information from the Dependencies
    numdep = 0;
    dep[numdep++].Load(FetchSlabs, finished_slab);
    dep[numdep++].Load(TransposePos, finished_slab);
    dep[numdep++].Load(NearForce, finished_slab);
    dep[numdep++].Load(TaylorForce, finished_slab);
    dep[numdep++].Load(Kick, finished_slab);
    dep[numdep++].Load(MakeCellGroups, finished_slab);
    dep[numdep++].Load(FindCellGroupLinks, finished_slab);
    dep[numdep++].Load(DoGlobalGroups, finished_slab);
    dep[numdep++].Load(Output, finished_slab);
    dep[numdep++].Load(Microstep, finished_slab);
    dep[numdep++].Load(FinishGroups, finished_slab);
    dep[numdep++].Load(Drift, finished_slab);
    dep[numdep++].Load(Finish, finished_slab);
    dep[numdep++].Load(LPTVelocityReRead, finished_slab);
    dep[numdep++].LoadCG(finished_slab);
    	// LoadCG() includes moving info into the CellGroupArenas
    assertf(numdep<MAXDEPENDENCY, "numdep has overflowed its MAX value");

    STDLOG(1,"Loading Arenas into the SendManifest\n");
    // Now load all of the arenas into the Manifest
    for (int type=0; type<MAXIDS; type++) {
	// Loop over all SlabTypes
	for (int s=finished_slab-1; s>finished_slab-cpd; s--) {
	    // Check each trailing slab; if present, load it up
	    if (LBW->IDPresent(type,s)) LoadArena(type,s);
		else break;
	}
    }

    STDLOG(1,"Loading Insert List into the SendManifest\n");
    // Partition the Insert List, malloc *il, and save it off
    /// for (int j=0;j<IL->length;j++) assertf(IL->list[j].xyz.x>=0&&IL->list[j].xyz.x<P.cpd, "Bad value at IL[%d] %d %d %d\n", j, IL->list[j].xyz.x, IL->list[j].xyz.y, IL->list[j].xyz.z);

    uint64 mid = ParallelPartition(IL->list, IL->length, finished_slab, is_below_slab);
    /// for (int j=0;j<IL->length;j++) assertf(IL->list[j].xyz.x>=0&&IL->list[j].xyz.x<P.cpd, "Bad value after PP at IL[%d] %d %d %d\n", j, IL->list[j].xyz.x, IL->list[j].xyz.y, IL->list[j].xyz.z);

    numil = IL->length-mid;
    int ret = posix_memalign((void **)&il, 4096, sizeof(ilstruct)*numil);
    memcpy(il, IL->list+mid, sizeof(ilstruct)*numil);
	// Possible TODO: Consider whether this copy should be multi-threaded
    STDLOG(1, "Insert list had size %l, now size %l\n", IL->length, mid);
    IL->ShrinkMAL(mid);
    /// for (int j=0;j<numil;j++) assertf(il[j].xyz.x>=0&&il[j].xyz.x<P.cpd, "Bad value in queued il[%d] %d %d %d\n", j, il[j].xyz.x, il[j].xyz.y, il[j].xyz.z);

    STDLOG(1,"Loading GroupLink List into the SendManifest\n");
    // Partition the GroupLink List, malloc *links, and save it off
    // TODO: Do these group finding variables always exist?
    if (GFC!=NULL) {
	mid = ParallelPartition(GFC->GLL->list, GFC->GLL->length, finished_slab, link_below_slab);
	ret = posix_memalign((void **)&links, 4096, sizeof(GroupLink)*(GFC->GLL->length-mid));
	numlinks = GFC->GLL->length-mid;
	memcpy(links, GFC->GLL->list+mid, sizeof(GroupLink)*numlinks);
	    // Possible TODO: Consider whether this copy should be multi-threaded
	STDLOG(1, "Grouplink list had size %l, now size %l\n", GFC->GLL->length, mid);
	GFC->GLL->ShrinkMAL(mid);
    }

    // TODO: Fork the communication thread and have it invoke this->Send()
    this->Send();
    return;
}



void Manifest::Send() {
    // This is the routine called by the communication thread.
    // It should send to its partner using blocking communication.
    // It can delete arenas as it finishes.

    // TODO: Send the Manifest.  Maybe wait for handshake?
    char fname[1024];
    size_t retval;
    sprintf(fname, "%s/manifest", P.WriteStateDirectory);
    FILE *fp = fopen(fname, "wb");
    retval = fwrite(this, sizeof(Manifest), 1, fp);
    STDLOG(1,"Sending the Manifest to %s\n", fname);

    for (int n=0; n<numarenas; n++) {
	// TODO: Send arenas[n].size bytes from arenas[n].ptr
	retval = fwrite(arenas[n].ptr, 1, arenas[n].size, fp);
	STDLOG(1,"Writing %l bytes of arenas to file\n", retval);
	// Now we can delete this arena
	LBW->DeAllocate(arenas[n].type, arenas[n].slab);
    }
    // TODO: Send the insert list: numil*sizeof(ilstruct) bytes from *il
    retval = fwrite(il, sizeof(ilstruct), numil, fp);
    STDLOG(1,"Writing %l objects of insert list to file\n", retval);
    STDLOG(1,"IL particle on slab %d\n", il[0].xyz.x);
    free(il);
    // TODO: Send the list: numlinks*sizeof(GroupLink) bytes from *links
    if (GFC!=NULL) {
	retval = fwrite(links, sizeof(GroupLink), numlinks, fp);
	STDLOG(1,"Writing %l objects of group links to file\n", retval);
	free(links);
    }
    completed = 1;
    // TODO: Can terminate the communication thread after this.
    fclose(fp);
}

void Manifest::Receive() {
    // This is the routine invoked by the communication thread.
    // It should store the results.
    // It allocates the needed space.

    // TODO: Receive the Manifest and overload the variables in *this.
    char fname[1024];
    size_t retval;
    sprintf(fname, "%s/manifest", P.WriteStateDirectory);
    FILE *fp = fopen(fname, "rb");
    retval = fread(this, sizeof(Manifest), 1, fp);
    STDLOG(1,"Reading Manifest from %s\n", fname);

    for (int n=0; n<numarenas; n++) {
	LBW->AllocateSpecificSize(arenas[n].type, arenas[n].slab, arenas[n].size);
	arenas[n].ptr = LBW->ReturnIDPtr(arenas[n].type, arenas[n].slab);
	// TODO: Receive arenas[n].size bytes into arenas[n].ptr
	retval = fread(arenas[n].ptr, 1, arenas[n].size, fp);
	STDLOG(1,"Reading %l bytes of arenas to file\n", retval);
    }

    int ret = posix_memalign((void **)&il, 64, numil*sizeof(ilstruct));
    assert(il!=NULL);
    // TODO: Receive numil*sizeof(ilstruct) bytes into *il
    retval = fread(il, sizeof(ilstruct), numil, fp);
    STDLOG(1,"Reading %l objects of insert list from file\n", retval);
    STDLOG(1,"IL particle on slab %d\n", il[0].xyz.x);

    ret = posix_memalign((void **)&links, 64, numlinks*sizeof(GroupLink));
    assert(links!=NULL);
    // TODO: Receive numlinks*sizeof(GroupLink) bytes into *links
    if (GFC!=NULL) {
	retval = fread(links, sizeof(GroupLink), numlinks, fp);
	STDLOG(1,"Reading %l objects of links from file\n", retval);
    }
    // TODO: Can terminate the communication thread after this
    completed = 1;
    fclose(fp);
}

void Manifest::ImportData() {
    // We invoke this within the timestep loop in a blocking manner.
    // This should get everything ready to use.

    // The arenas are already in place from the non-blocking thread,
    // but we want to keep them inactive until we're sure all the info
    // has been transmitted.  This can be done because the first event
    // in a slab requires the IOComplete flag on the particle data.
    // The Arena allocation did not set this flag; we do so here, thereby
    // opening this region for use.
    // TODO: Consider whether there are any other race conditions; DJE thinks not
    assertf(completed==1, "ImportData has been called when completed==%d\n", completed);

    STDLOG(1,"Importing ReceiveManifest into the flow\n");
    for (int n=0; n<numarenas; n++) {
	LBW->SetIOCompleted(arenas[n].type, arenas[n].slab);
    }

    // Set the dependencies.   Be careful that this keeps the Dependencies 
    // matched to the dep[] order set in QueueToSend().
    int n = 0;
    dep[n++].Set(FetchSlabs);
    dep[n++].Set(TransposePos);
    dep[n++].Set(NearForce);
    dep[n++].Set(TaylorForce);
    dep[n++].Set(Kick);
    dep[n++].Set(MakeCellGroups);
    dep[n++].Set(FindCellGroupLinks);
    dep[n++].Set(DoGlobalGroups);
    dep[n++].Set(Output);
    dep[n++].Set(Microstep);
    dep[n++].Set(FinishGroups);
    dep[n++].Set(Drift);
    dep[n++].Set(Finish);
    dep[n++].Set(LPTVelocityReRead);
    dep[n++].SetCG();
    	// This will copy data back to GFC from CellGroupArenas
    assert(n==numdep);

    // Add *il to the insert list
    /// for (int j=0;j<numil;j++) assertf(il[j].xyz.x>=0&&il[j].xyz.x<P.cpd, "Bad value at il[%d] %d %d %d\n", j, il[j].xyz.x, il[j].xyz.y, il[j].xyz.z);
    uint64 len = IL->length;
    STDLOG(1, "Growing IL list from %d by %d\n", len, numil);
    IL->GrowMAL(len+numil);
    memcpy(IL->list+len, il, numil*sizeof(ilstruct));
	// Possible TODO: Should this copy be multi-threaded?
    free(il);
    /// for (int j=0;j<IL->length;j++) assertf(IL->list[j].xyz.x>=0&&IL->list[j].xyz.x<P.cpd, "Bad value at IL[%d] %d %d %d\n", j, IL->list[j].xyz.x, IL->list[j].xyz.y, IL->list[j].xyz.z);

    /* 
    for (int j=len; j<IL->length; j++) {
    	if (IL->list[j].xyz.x>0&&IL->list[j].xyz.x<10)
	    fprintf(stderr,"Illegal IL particle: %d %d %d\n", 
	    	IL->list[j].xyz.x, IL->list[j].xyz.y, IL->list[j].xyz.z);
    	if (IL->list[j].xyz.x<0 || IL->list[j].xyz.x>=P.cpd
    	   || IL->list[j].xyz.x<0 || IL->list[j].xyz.x>=P.cpd
    	   || IL->list[j].xyz.x<0 || IL->list[j].xyz.x>=P.cpd)
	    fprintf(stderr,"Illegal IL particle: %d %d %d\n", 
	    	IL->list[j].xyz.x, IL->list[j].xyz.y, IL->list[j].xyz.z);
    }
    */

    // Add *links to the GroupLink list
    if (GFC!=NULL) {
	len = GFC->GLL->length;
	GFC->GLL->GrowMAL(len+numlinks);
	memcpy(GFC->GLL->list+len, links, numlinks*sizeof(GroupLink));
	// Possible TODO: Should this copy be multi-threaded?
	free(links);
	STDLOG(1, "Growing GroupLink list from %d by %d\n", len, numil);
    }
    
    // We're done with this Manifest!
    done();
    // TODO: Could erase the file, but we won't for now
}

/* ================= Putting this in the code ======================

Probably easiest if we early on query the disk to get the range of
slabs in the ReadState and hence the number of slabs.  Or maybe the
range of slabs is a property of the state, on the idea that in the
parallel code, each node has its own state directory.

Need to adjust ForceSlabPrecondition to not attempt to read
beyond the slabs on disk.

Need to adjust the overall timestep.cpp loop to complete the
correct number of slabs, which I think should be the number that
are on disk.  E.g., !Finish.alldone() needs to change to
while(Finish.number_of_slabs_executed<number_of_local_slabs).  

Probably create two global instances of Manifest, one to SendManifest and one
to ReceiveManifest.  Need to launch the Receive thread at the beginning, 
running Receive().  The Send thread could in principle wait until 
QueueToSend() is executed.

Add some logic to end of the Finish action, e.g.,
    if (Finish.number_of_slabs_executed==0) SendManifest.QueueToSend(slab);

And add a test to inside the timestep loop, e.g.,
    if (ReceiveManifest.is_ready()) ReceiveManifest.ImportData();



For group finding, there is information buried in the GFC that we need 
to pass along.  When we finish the first slab, we have necessarily 
closed GlobablGroups in the neighboring slabs.  I think there is a 
rare condition in which GlobalGroups have been formed in a slab, but
not yet finished/deleted, because the Kick.notdone has blocked the
Output.  That's a pain, because it means that we have to save the 
GlobalGroupSlab GFC->globalslabs[].

CellGroups are also destroyed when GlobalGroups are Finished.
At that time, we mark cellgroups_status to 2.  

The destructor for GroupFindingControl requires all cellgroups_status
to be 2, but maybe we want this to simply be !=1.  We never actually 
test >0 in our dependencies.

If we changed the dependencies so that we never had an open GlobalGroupSlab
at the time when QueueToSend is being called, then we would only need to
send CellGroups and cellgroups_status.  Simpler!

*/


