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
	}
	// We've found the first notdone slab
	begin++;   // Want to pass the first done one
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
	end = finished_slab-1;
	while (GFC->cellgroups_status[PP->WrapSlab(end)]==2) end--;
	end++;   // Now this marks the first ==2.
	for (begin=end-1; begin>end-GFC->cpd; begin--) {
	    int s = PP->WrapSlab(begin);
	    if (GFC->cellgroups_status[s]==0) break;
	    // Need to load this over to the Arenas
	    GFC->cellgroups[s].pack(CellGroupArena,s);
	    // Having moved it, we can delete the original
	    GFC->DestroyCellGroups(s);
	}
	begin++;   // Now this marks the first ==1.
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
    int completed;	// ==1 if the Manifest is completed, ==0 if not.

    Manifest() {
    	numarenas = numil = numlinks = numdep = 0;
	completed = 0;
	il = NULL;
	links = NULL;
	return;
    }
    ~Manifest() { }

    inline int is_ready() { return completed; }
	// Call this to see if the Manifest is ready to retrieve

    void LoadArena(int type, int s) {
	ManifestArena *a = arenas+numarenas;
	a->type = type;
	a->slab = s;
	a->size = LBW->IDSizeBytes(type, s);
	a->ptr =  LBW->ReturnIDPtr(type, s);
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

    // Now load all of the arenas into the Manifest
    for (int type=0; type<MAXIDS; type++) {
	// Loop over all SlabTypes
	for (int s=finished_slab-1; s>finished_slab-cpd; s--) {
	    // Check each trailing slab; if present, load it up
	    if (LBW->IDPresent(type,s)) LoadArena(type,s);
		else break;
	}
    }

    // Partition the Insert List, malloc *il, and save it off
    uint64 mid = ParallelPartition(IL->list, IL->length, finished_slab, is_below_slab);
    int ret = posix_memalign((void **)&il, 4096, sizeof(ilstruct)*(IL->length-mid));
    memcpy(il, IL->list+mid, IL->length-mid);
	// Possible TODO: Consider whether this copy should be multi-threaded
    IL->ShrinkMAL(mid);

    // Partition the GroupLink List, malloc *links, and save it off
    // TODO: Do these group finding variables always exist?
    mid = ParallelPartition(GFC->GLL->list, GFC->GLL->length, finished_slab, link_below_slab);
    ret = posix_memalign((void **)&links, 4096, sizeof(GroupLink)*(GFC->GLL->length-mid));
    memcpy(links, GFC->GLL->list+mid, GFC->GLL->length-mid);
	// Possible TODO: Consider whether this copy should be multi-threaded
    GFC->GLL->ShrinkMAL(mid);

    // TODO: Fork the communication thread and have it invoke this->Send()
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

    for (int n=0; n<numarenas; n++) {
	// TODO: Send arenas[n].size bytes from arenas[n].ptr
	retval = fwrite(arenas[n].ptr, 1, arenas[n].size, fp);
	// Now we can delete this arena
	LBW->DeAllocate(arenas[n].type, arenas[n].slab);
    }
    // TODO: Send the insert list: numil*sizeof(ilstruct) bytes from *il
    retval = fwrite(il, sizeof(ilstruct), numil, fp);
    free(il);
    // TODO: Send the list: numlinks*sizeof(GroupLink) bytes from *links
    retval = fwrite(links, sizeof(GroupLink), numlinks, fp);
    free(links);
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

    for (int n=0; n<numarenas; n++) {
	LBW->AllocateSpecificSize(arenas[n].type, arenas[n].slab, arenas[n].size);
	arenas[n].ptr = LBW->ReturnIDPtr(arenas[n].type, arenas[n].slab);
	// TODO: Receive arenas[n].size bytes into arenas[n].ptr
	retval = fread(arenas[n].ptr, 1, arenas[n].size, fp);
    }

    int ret = posix_memalign((void **)&il, 64, numil*sizeof(ilstruct));
    assert(il!=NULL);
    // TODO: Receive numil*sizeof(ilstruct) bytes into *il
    retval = fread(il, sizeof(ilstruct), numil, fp);

    ret = posix_memalign((void **)&links, 64, numlinks*sizeof(GroupLink));
    assert(links!=NULL);
    // TODO: Receive numlinks*sizeof(GroupLink) bytes into *links
    retval = fread(links, sizeof(GroupLink), numlinks, fp);
    completed = 1;
    // TODO: Can terminate the communication thread after this
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
    uint64 len = IL->length;
    IL->GrowMAL(len+numil);
    memcpy(IL->list+len, il, numil*sizeof(ilstruct));
	// Possible TODO: Should this copy be multi-threaded?
    free(il);

    // Add *links to the GroupLink list
    len = GFC->GLL->length;
    GFC->GLL->GrowMAL(len+numlinks);
    memcpy(GFC->GLL->list+len, links, numlinks*sizeof(GroupLink));
	// Possible TODO: Should this copy be multi-threaded?
    free(links);
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


