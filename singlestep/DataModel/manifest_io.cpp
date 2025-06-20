// Copyright 2012-2025 The Abacus Developers
// SPDX-License-Identifier: GPL-3.0-or-later

/** \file We can implement a simple parallelism into abacus by
doing a slab-oriented decomposition, in which different nodes are
responsible for different slabs.

A given node will start with a particular range of slabs on disk.
It starts by loading the first one in the sequence and continues
upward, doing work, until it has finished the first slab it 
can (say, N).  This means that all information relevant to slab N
is completed, which implies that the information relevant to N-1 and
to N+1 are now rigorously disconnected.  At that point, we package
up all of the slabs <N and send them to the lower neighbor.  We
similarly send the portion of the insert list and GroupLink list
that involve <N.  That node will then have all of the info needed
to finish slab N-1.

Notice that this means that the domain of a given node shifts 
upward from one epoch to the next, which is fine.  We simply need
to detect the range of slabs on disk, and act appropriately.

This file contains the code to do this transfer.  Most of the work
is done in a non-blocking manner by two threads per node, one sending,
one receiving.  We execute a blocking routine to load up the 
Manifest of data to be sent and then a blocking routine to store
the data when it is received.  

*/

// ================== Manifest helpers =====================================

/// The information we're passing for a single arena
class ManifestArena {
  public:
    uint64 size;    ///< Arena size in bytes of usable space, not including guards
    char *ptr;	///< The pointer to the space
    int type;    ///< The SlabType
    int slab; 	///< slab number
    ManifestArena() { size = 0; type = slab = 0; ptr = NULL; }
};

/// The information we're passing for a single dependency
class DependencyRecord {
  public:
    int begin;
    int end;
    ///< We're passing the information that [begin,end) should be set by force
    ///< on the receiving node.
    ///< We require begin<=end; the wrapping is done within the Dependency class.
    DependencyRecord() { begin = end = 0; }

    /// Load the status of this dependency 
    void Load(Dependency &d, int finished_slab) {
        end = finished_slab;
        for (begin=end-1; begin>end-d.cpd; begin--) {
            if (d.notdone(begin)) break;
            //// d.mark_to_repeat(begin);  // We're not unmarking anymore
        }
        // We've found the first notdone slab
        begin++;   // Want to pass the first done one
        // Now look for the last done slab
        for (;end<finished_slab+d.cpd;end++) {
            if (d.notdone(end)) break;
        } end--;
        STDLOG(2, "Load Dependency [{:d},{:d})\n", begin, end);
        return;
    }
    void Set(Dependency &d) {
        for (int s=begin; s<end; s++) d.force_done(s);
        STDLOG(2, "Set Dependency [{:d},{:d})\n", begin, end);
        return;
    }

    /// We also are going to use this same logic to hold the cellgroup_status[]
    /// information from GFC.  In particular, we need to transmit all cases
    /// where this is ==1.  We don't care about ==2.
    void LoadCG(int finished_slab) {
        if (GFC==NULL) { begin=end = finished_slab-1; return;}
        STDLOG(2, "LoadCG start\n");
        end = finished_slab-1;
        while (GFC->cellgroups_status[CP->WrapSlab(end)]==2) end--;
        end++;   // Now this marks the first ==2.
        for (begin=end-1; begin>end-GFC->cpd; begin--) {
            int s = CP->WrapSlab(begin);
            if (GFC->cellgroups_status[s]==0) break;
            // Need to load this over to the Arenas
            STDLOG(2,"Packing CellGroupArena for slab {:d}\n", s);
            GFC->cellgroups[s].pack(CellGroupArena,s);
            // Having moved it, we can delete the original
            GFC->DestroyCellGroups(s);
        }
        begin++;   // Now this marks the first ==1.
        STDLOG(2, "CG Dependency [{:d},{:d})\n", begin, end);
        return;
    }

    /// Set the cellgroups_status to 1 for the indicated slabs
    void SetCG(int finished_slab) {
        if (GFC==NULL) return;
        for (int s=begin; s<end; s++) {
            GFC->cellgroups_status[CP->WrapSlab(s)]=1;
            // And move the information from the Arenas back into the SlabArray
            GFC->cellgroups[CP->WrapSlab(s)].unpack(CellGroupArena,s);
            // Now we can delete the CellGroupArena
            SB->DeAllocate(CellGroupArena,s);
        }
        STDLOG(2, "Marking cellgroups_status = 1 for [{:d},{:d})\n", begin, end);
        // We need to set cellgroups_status=2 for [end,end+2*GroupRadius]
        // This is so that CreateGlobalGroups() doesn't crash
        for (int s=end; s<=end+2*GFC->GroupRadius; s++)
            GFC->cellgroups_status[CP->WrapSlab(s)]=2;
        STDLOG(2, "Marking cellgroups_status = 2 for [{:d},{:d}]\n", end, end+2*GFC->GroupRadius);
        return;
    }
};

int global_minslab_search;    // Used to pass an extra variable into partition

// Need some new Partition functions, to get ranges of slabs
// Will select things between [x,slab), where x=slab-CPD/2, handling the wrap
inline bool is_below_slab(ilstruct *particle, int slab) {
    int x = particle->xyz.x-slab;  // We want x<0, but not too much.
    x = CP->WrapSlab(x);
    return x>=global_minslab_search;
}

/// When we partition the GroupLink list, GlobalGroups have already
/// been closed in slab, so there are no GroupLinks that involve the 
/// slab.  By definition, the two ends of the link can differ by at 
/// most one slab.  So it is enough to test one end.
inline bool link_below_slab(GroupLink *link, int slab) {
    int sa = link->a.slab();
    sa = CP->WrapSlab(sa-slab);
    return sa>=global_minslab_search;
}

// prototypes
void *ManifestSendThread(void *p);
void *ManifestReceiveThread(void *p);

// ================== Manifest Class =====================================

// These are oversized, but still tiny compared to what will be transferred.
#define MAXMANIFEST 1024
#define MAXDEPENDENCY 64

/// This is all of the information that we'll be passing between nodes
struct ManifestCore {
    ManifestArena arenas[MAXMANIFEST];   ///< All of the arena info
    int numarenas;   ///< The number of arenas being sent
    int numil;	///< The number of IL objects
    int numlinks;	///< The number of GroupLink objects
    DependencyRecord dep[MAXDEPENDENCY];   ///< The dependency info
    int numdep;  ///< And the number of dependencies, just to check.
    int remote_first_slab_finished;
};

/// This is the class for sending information between the nodes.
class Manifest {
  public:
    ManifestCore m;   ///< The info we're sending over

    // Temporary storage buffers
    ilstruct *il;	///< Storage for the IL
    GroupLink *links;   ///< Storage for the GLL
    
    // Control logic
    int completed;	///< ==2 if the Manifest is ready for use, ==0 if not yet , ==3 if we've already imported it, ==1 if transfer in progress.
    STimer Load;        ///< The timing for the Queue & Import blocking routines
    STimer Transmit;	///< The timing for the Send & Receive routines, usually non-blocking
    STimer CheckCompletion;	///< The timing for the blocking checks
    size_t bytes;       ///< The number of bytes received

    int blocking;	///< =1 if blocking, =0 if non-blocking
    pthread_t thread;	///< The thread object
    int launched;	///< =1 if the thread was launced

    std::string RecNodeString;     // The node string we're reading from

    Manifest() {
    	m.numarenas = m.numil = m.numlinks = m.numdep = 0;
        completed = 0;
        bytes = 0;
        blocking = 1;
        launched = 0;
        il = NULL;
        links = NULL;
        return;
    }
    ~Manifest() { 
        void *p;
        // We want to do a blocking pthread_join here, because
        // if one node has finished early, before the Send has been
        // completed, then we need to hang on until it has.
        if (launched==1) int retval = pthread_join(thread, NULL); 
    }

    void set_blocking() { blocking = 1; }
    void set_nonblocking() { blocking = 0; }

    inline int is_ready() { if (completed==2) return 1; return 0;}
	// Call this to see if the Manifest is ready to retrieve
    void done() { completed = 3; }

    void LoadArena(int type, int s) {
        ManifestArena *a = m.arenas+m.numarenas;
        a->type = type;
        a->slab = s;
        a->size = SB->SlabSizeBytes(type, s);
        a->ptr =  SB->GetSlabPtr(type, s);
        STDLOG(2, "Queuing slab {:d} of type {:d}, size {:d}\n", s, type, a->size);
        m.numarenas++;
        assertf(m.numarenas<MAXMANIFEST, "numarenas has overflowed; increase MAXMANIFEST.");
        return;
    }

    // Here's the prototypes for the main routines
    void QueueToSend(int finished_slab);
    void Send();
    int FreeAfterSend();
    int Check();
    void SetupToReceive();
    void Receive();
    void ImportData();

    /// Create the non-blocking Send thread
    void LaunchSendThread() {
        assertf(blocking==0, "Asked to start Send Thread, but blocking = {:d}\n", blocking);
        int retval = pthread_create(&thread, NULL, ManifestSendThread, (void *)this);
        assertf(retval==0, "Failed to create SendManifest thread! {:d}", retval);
        launched = 1;
        STDLOG(2, "Launched SendManifest Thread\n");
    }

    /// Create the non-blocking Receive thread
    void LaunchReceiveThread() {
        assertf(blocking==0, "Asked to start Receive Thread, but blocking = {:d}\n", blocking);
        int retval = pthread_create(&thread, NULL, ManifestReceiveThread, (void *)this);
        assertf(retval==0, "Failed to create ReceiveManifest thread! {:d}", retval);
        launched = 1;
        STDLOG(2, "Launched ReceiveManifest Thread\n");
    }
};    


/// Here are our outgoing and incoming Manifest instances
Manifest *SendManifest, *ReceiveManifest; 

/// Call this routine to turn on Non-blocking Manifest
void SetupManifest() {
    SendManifest = new Manifest;
    ReceiveManifest = new Manifest;
    #ifdef PARALLEL
    int rank = MPI_rank+1;
    if (rank>=MPI_size) rank-=MPI_size;
    RecNodeString = fmt::format(".{:04d}", rank);
    #ifdef IOTHREADED
        // The following routines turn on the non-blocking I/O
        // For now we tie this to IO Threading, but it's not the same threads
        STDLOG(2,"Turning on non-blocking Manifest I/O Code\n");
        SendManifest->set_nonblocking();
        ReceiveManifest->set_nonblocking();
        ReceiveManifest->LaunchReceiveThread();
    #else
        STDLOG(2,"Turning on blocking Manifest I/O Code\n");
    #endif
    #endif
}

void FreeManifest() {
    delete SendManifest;
    delete ReceiveManifest;
}

void *ManifestSendThread(void *p) {
    Manifest *m = (Manifest *)p;
    m->Send();
    return NULL;
}

void *ManifestReceiveThread(void *p) {
    Manifest *m = (Manifest *)p;
    const fs::path fname = P.WriteStateDirectory / fmt::format("manifest_done{:s}", m->RecNodeString);

    while (fs::exists(fname)) usleep(10000);
    	// Wait until the file exists
    m->Receive();
    return NULL;
}




// ================  Routine to define the outgoing information =======

/** The blocking routine to prepare information to Send

We've just finished the given slab.  
Load all of the information from lower slabs into the Manifest.
Then call the non-blocking communication and return.
*/

void Manifest::QueueToSend(int finished_slab) {
    int cpd = P.cpd;

    // TODO: Consider some safety measures here.
    // E.g., we really should not have any slabs that have open GlobalGroupSlabs
    // in our wake.  Might check, since one could screw up the Dependencies.

    first_slab_finished = finished_slab;   // A global record of this
    m.remote_first_slab_finished = finished_slab; 
    STDLOG(1,"Queueing the SendManifest at slab={:d}\n", finished_slab);
    Load.Start();

    // Load the information from the Dependencies
    m.numdep = 0;
    m.dep[m.numdep++].Load(FetchSlabs, finished_slab);
    m.dep[m.numdep++].Load(TransposePos, finished_slab);
    m.dep[m.numdep++].Load(NearForce, finished_slab);
    m.dep[m.numdep++].Load(TaylorForce, finished_slab);
    m.dep[m.numdep++].Load(Kick, finished_slab);
    m.dep[m.numdep++].Load(MakeCellGroups, finished_slab);
    m.dep[m.numdep++].Load(FindCellGroupLinks, finished_slab);
    int min_links_slab = m.dep[m.numdep-1].begin-1;
    	// We just determined that FindCellGroupLinks has executed on begin, so
	// the GLL might contain links including begin-1.
    m.dep[m.numdep++].Load(DoGlobalGroups, finished_slab);
    m.dep[m.numdep++].Load(Output, finished_slab);
    m.dep[m.numdep++].Load(FinishGroups, finished_slab);
    m.dep[m.numdep++].Load(Drift, finished_slab);
    int min_il_slab = m.dep[m.numdep-1].begin-1;
    	// We just determined that Drift has executed on begin, so
	// the rebinning might have taken particles to begin-1.
    m.dep[m.numdep++].Load(Finish, finished_slab);
    m.dep[m.numdep++].LoadCG(finished_slab);
    	// LoadCG() includes moving info into the CellGroupArenas
    assertf(m.numdep<MAXDEPENDENCY, "m.numdep has overflowed its MAX value");

    STDLOG(2,"Queuing Arenas into the SendManifest\n");
    // Now load all of the arenas into the Manifest
    int min_slab = finished_slab;
    for (int type=0; type<SB->NumTypes; type++) {
	    // Loop over all SlabTypes
	    for (int s=finished_slab-1; s>finished_slab-cpd; s--) {
	        // Check each trailing slab; if present, load it up
	        if (SB->IsSlabPresent(type,s)) {
	    	    LoadArena(type,s);
		        min_slab = std::min(min_slab, s);
                if (type==PosSlab) {
                    // Zero out the SlabSize for all particles being sent
                    SS->setold(s,0);
                }
	        }
	        else if (s<min_slab) break;
	        // Some slabs have been already been deallocated by the finish slab,            // but we need the ones that were waiting for the periodic wrap.
	        // This relies on the fact that the first SlabType, PosSlab, 
	        // stretches back to the beginning.
	    }
    }
    STDLOG(2,"Done Queuing Arenas, spanning [{:d},{:d})\n", min_slab, finished_slab);

    STDLOG(2,"Queuing Insert List into the SendManifest, extracting [{:d},{:d})\n",
    	min_il_slab, finished_slab);
    // Partition the Insert List, malloc *il, and save it off
    IL->CollectGaps();     // Assure there are no MALgaps
    global_minslab_search = CP->WrapSlab(min_il_slab-finished_slab);
    uint64 mid = ParallelPartition(IL->list, IL->length, finished_slab, is_below_slab);

    m.numil = IL->length-mid;
    int ret = posix_memalign((void **)&il, 4096, sizeof(ilstruct)*m.numil);
    memcpy(il, IL->list+mid, sizeof(ilstruct)*m.numil);
	// Possible TODO: Consider whether this copy should be multi-threaded
    STDLOG(2, "Insert list had size {:d}, now size {:d}; sending {:d}\n", IL->length, mid, m.numil);
    IL->ShrinkMAL(mid);

    // Partition the GroupLink List, malloc *links, and save it off
    // TODO: Do these group finding variables always exist?
    if (GFC!=NULL) {
	STDLOG(2,"Queuing GroupLink List into the SendManifest, extracting [{:d},{:d})\n", min_links_slab, finished_slab);
	global_minslab_search = CP->WrapSlab(min_links_slab-finished_slab);
	mid = ParallelPartition(GFC->GLL->list, GFC->GLL->length, finished_slab, link_below_slab);
	ret = posix_memalign((void **)&links, 4096, sizeof(GroupLink)*(GFC->GLL->length-mid));
	m.numlinks = GFC->GLL->length-mid;
	memcpy(links, GFC->GLL->list+mid, sizeof(GroupLink)*m.numlinks);
	    // Possible TODO: Consider whether this copy should be multi-threaded
	STDLOG(2, "Grouplink list had size {:d}, now size {:d}; sending {:d}\n", GFC->GLL->length, mid, m.numlinks);
	GFC->GLL->ShrinkMAL(mid);
    }
    Load.Stop();

    // TODO: Fork the communication thread and have it invoke this->Send()
    if (blocking) {
	STDLOG(2, "Preparing to Send the Manifest by blocking method.\n");
    	this->Send(); 
    } else {
    	STDLOG(2, "Forking the SendManifest->Send() thread.\n");
    	LaunchSendThread();
    }
    return;
}

// =============== Routine to actually transmit the outgoing info ====

/// This is the non-blocking routine called by the communication thread.
/// It should send to its partner using blocking communication.
/// It can delete arenas as it finishes.

void Manifest::Send() {
    Transmit.Start();
    // TODO: Send the ManifestCore.  Maybe wait for handshake?
    size_t retval;
    FILE *fp = fopen( (P.WriteStateDirectory / ("manifest" + NodeString)).c_str(), "wb");
    bytes = fwrite(&m, sizeof(ManifestCore), 1, fp)*sizeof(ManifestCore);
    STDLOG(1,"Sending the Manifest Core to {}\n", fname);

    for (int n=0; n<m.numarenas; n++) {
        // TODO: Send arenas[n].size bytes from arenas[n].ptr
        retval = fwrite(m.arenas[n].ptr, 1, m.arenas[n].size, fp);
        bytes += retval;
        if (blocking) STDLOG(2,"Writing {:d} bytes of arenas to file\n", retval);
        // Now we can delete this arena
        SB->DeAllocate(m.arenas[n].type, m.arenas[n].slab);
    }
    // TODO: Send the insert list: m.numil*sizeof(ilstruct) bytes from *il
    retval = fwrite(il, sizeof(ilstruct), m.numil, fp);
    bytes += retval*sizeof(ilstruct);
    if (blocking) STDLOG(2,"Writing {:d} objects of insert list to file\n", retval);
    if (blocking) STDLOG(2,"IL particle on slab {:d}\n", il[0].xyz.x);
    free(il);
    // TODO: Send the list: m.numlinks*sizeof(GroupLink) bytes from *links
    if (GFC!=NULL) {
        retval = fwrite(links, sizeof(GroupLink), m.numlinks, fp);
        bytes += retval*sizeof(GroupLink);
        if (blocking) STDLOG(2,"Writing {:d} objects of group links to file\n", retval);
        free(links);
    }
    fclose(fp);
    // usleep(5e6);   // Just force a wait here, to see what the code does.
    completed = 2;
    // TODO: Can terminate the communication thread after this.
    // Touch a file to indicate that we're done
    fp=(fopen(P.WriteStateDirectory / ("manifest_done" + NodeString)).c_str(),"w");
    fclose(fp);
    Transmit.Stop();
}

/// This is the routine to call frequently to try to clean up 
/// space after Send's have happened.
inline int Manifest::FreeAfterSend() {
    return 0;
}


// =============== Routine to actually receive the incoming info ====

/// We have to issue the request to receive the ManifestCore
void Manifest::SetupToReceive() {
}

/// This can be called in timestep.cpp to manually check (blocking)
/// whether Receive is ready to run.
// TODO: This may have a rather different meaning in MPI
inline int Manifest::Check() {
    #ifndef PARALLEL
    return;	// If we're not doing PARALLEL, let this optimize to a no-op
    #endif

    if (blocking==0) return; 	// We're doing this by a thread
    if (completed>0) return;	// We've already been here once
    fs::path fname = P.WriteStateDirectory / fmt::format("manifest_done{:s}", RecNodeString);
    if (fs::exists(fname)) return;
    	// The signal file doesn't yet exist, so try again later
    // Otherwise, we're ready to go
    STDLOG(2,"Check indicates we are ready to Receive the Manifest\n");
    Receive();
}

/// This is the routine invoked by the communication thread.
/// It should store the results.
/// It allocates the needed space.
void Manifest::Receive() {
    Transmit.Start();
    // TODO: Receive the Manifest and overload the variables in *this.
    size_t retval;
    fs::path fname = P.WriteStateDirectory / fmt::format("manifest{:s}", RecNodeString);
    FILE *fp = fopen(fname.c_str(), "rb");
    bytes += fread(&m, sizeof(ManifestCore), 1, fp)*sizeof(ManifestCore);
    STDLOG(1,"Reading Manifest Core from {}\n", fname);

    for (int n=0; n<m.numarenas; n++) {
        SB->AllocateSpecificSize(m.arenas[n].type, m.arenas[n].slab, m.arenas[n].size);
        // TODO: Need to force these to be in memory with RamDisk=NO
        m.arenas[n].ptr = SB->GetSlabPtr(m.arenas[n].type, m.arenas[n].slab);
        // TODO: Receive arenas[n].size bytes into arenas[n].ptr
        retval = fread(m.arenas[n].ptr, 1, m.arenas[n].size, fp);
        bytes += retval;
        if (blocking) STDLOG(2,"Reading {:d} bytes of arenas to file\n", retval);
    }

    int ret = posix_memalign((void **)&il, 64, m.numil*sizeof(ilstruct));
    assert(il!=NULL);
    // TODO: Receive m.numil*sizeof(ilstruct) bytes into *il
    retval = fread(il, sizeof(ilstruct), m.numil, fp);
    bytes += retval*sizeof(ilstruct);
    if (blocking) STDLOG(2,"Reading {:d} objects of insert list from file\n", retval);
    if (blocking) STDLOG(2,"IL particle on slab {:d}\n", il[0].xyz.x);

    ret = posix_memalign((void **)&links, 64, m.numlinks*sizeof(GroupLink));
    assert(links!=NULL);
    // TODO: Receive m.numlinks*sizeof(GroupLink) bytes into *links
    if (GFC!=NULL) {
        retval = fread(links, sizeof(GroupLink), m.numlinks, fp);
        bytes += retval*sizeof(GroupLink);
        if (blocking) STDLOG(2,"Reading {:d} objects of links from file\n", retval);
    }
    // TODO: Can terminate the communication thread after this
    completed = 2;
    fclose(fp);
    Transmit.Stop();
}



// ============== Routine to import the incoming info =============

/** We invoke this within the timestep loop in a blocking manner
after the ReceiveManifest is marked as completed.
This should get everything ready to use.

The arenas are already in place from the non-blocking thread,
but we want to keep them inactive until we're sure all the info
has been transmitted.  This can be done because the first event
in a slab requires the IOComplete flag on the particle data.
The Arena allocation did not set this flag; we do so here, thereby
opening this region for use.
*/

void Manifest::ImportData() {
    assertf(completed==2, "ImportData has been called when completed=={:d}\n", completed);

    STDLOG(1,"Importing ReceiveManifest of {:d} bytes into the flow\n", bytes);
    Load.Start();
    for (int n=0; n<m.numarenas; n++) {
	    SB->SetIOCompleted(m.arenas[n].type, m.arenas[n].slab);
	    STDLOG(2,"Completing Import of arena slab {:d} of type {:d} and size {:d}\n", 
	    	m.arenas[n].slab, m.arenas[n].type, m.arenas[n].size);
        if (m.arenas[n].type==PosSlab) {
            // Set the SlabSize based on the newly arrived PosSlab
            SS->setold(m.arenas[n].slab, m.arenas[n].size/sizeof(posstruct));
        }
    }

    // Set the dependencies.   Be careful that this keeps the Dependencies 
    // matched to the dep[] order set in QueueToSend().
    int n = 0;
    m.dep[n++].Set(FetchSlabs);
    m.dep[n++].Set(TransposePos);
    m.dep[n++].Set(NearForce);
    m.dep[n++].Set(TaylorForce);
    m.dep[n++].Set(Kick);
    m.dep[n++].Set(MakeCellGroups);
    m.dep[n++].Set(FindCellGroupLinks);
    m.dep[n++].Set(DoGlobalGroups);
    m.dep[n++].Set(Output);
    m.dep[n++].Set(FinishGroups);
    m.dep[n].end++;
        // We need to also mark the first_finished_slab on the other node
        // as having been completed.  This is because Finish requires this
        // as a precondition, and the particles incoming to Finish are on the
        // insert list.
    m.dep[n++].Set(Drift);
    m.dep[n++].Set(Finish);
    m.dep[n++].SetCG(m.remote_first_slab_finished);
    	// This will copy data back to GFC from CellGroupArenas
    assert(n==m.numdep);

    // Add *il to the insert list
    IL->CollectGaps();     // Assure there are no MALgaps
    uint64 len = IL->length;
    IL->GrowMAL(len+m.numil);
    STDLOG(2, "Growing IL list from {:d} by {:d} = {:d}\n", len, m.numil, IL->length);
    memcpy(IL->list+len, il, m.numil*sizeof(ilstruct));
	// Possible TODO: Should this copy be multi-threaded?
    free(il);

    // Add *links to the GroupLink list
    if (GFC!=NULL) {
	len = GFC->GLL->length;
	GFC->GLL->GrowMAL(len+m.numlinks);
	memcpy(GFC->GLL->list+len, links, m.numlinks*sizeof(GroupLink));
	// Possible TODO: Should this copy be multi-threaded?
	free(links);
	STDLOG(2, "Growing GroupLink list from {:d} by {:d} = {:d}\n", len, m.numil, GFC->GLL->length);
    }
    
    // We're done with this Manifest!
    // TODO: Could erase the file, but we won't for now
    done();
    Load.Stop();
}
