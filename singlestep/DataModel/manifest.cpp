
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

#include "mpi_limiter.h"

// #define NO_MPI    // Just keep these routines blocked out for now
#ifdef NO_MPI
#include "manifest_io.cpp"
#else
// It was too confusing to keep the MPI and I/O based codes in the same
// file.  The I/O code is just a touch-stone single-node version.

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
    void Load(Dependency &d, int finished_slab, const char label[]) {
        end = finished_slab;
		
        // Don't look past the slabs this node owns
        for (begin=end-1; CP->WrapSlab(begin - first_slab_on_node) < total_slabs_on_node; begin--) {
            if (!d.instantiated || d.notdone(begin)) break;
            //// d.mark_to_repeat(begin);   // We're not going to unmark
        }
		
        // We've found the first notdone slab
        begin++;   // Want to pass the first done one
		
        // Now look for the last done slab
        // Don't look past the slabs this node owns
        for (; CP->WrapSlab(end - first_slab_on_node) < total_slabs_on_node; end++) {
            if (!d.instantiated || d.notdone(end)) break;
        } end--;
        /* The manifest code is called in two different use cases:
           1) when the GlobalGroups have happened and we want to pass 
           along the earlier info so that Groups can proceed on the neighbor.
           2) when Finish has happened and we want to get rid of everything
           relevant to the earlier slabs.
           The issue is that case (1) needs to indicate the existence of 
           data that hasn't yet been shipped.  In particular, that FindLinks has happened on forward slabs.
        */
        if (Drift.notdone(finished_slab)) {
            // We're in the first case; don't pass along anything ahead.
            STDLOG(3, "Changing end from %d to %d\n", end, finished_slab);
            end = finished_slab;
        }
        STDLOG(2, "Load Dependency %s [%d,%d)\n", label, begin, end);
        return;
    }
    void Set(Dependency &d, const char label[]) {
        for (int s=begin; s<end; s++) d.force_done(s);
        STDLOG(2, "Set Dependency %s [%d,%d)\n", label, begin, end);
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
            if (GFC->cellgroups_status[s]!=1) break;
            // Need to load this over to the Arenas
            STDLOG(2,"Packing CellGroupArena for slab %d\n", s);
            GFC->cellgroups[s].pack(CellGroupArena,s);
            // Having moved it, we can delete the original
            GFC->DestroyCellGroups(s);  // This sets cellgroups_status to 2
        }
        begin++;   // Now this marks the first ==1.
        STDLOG(2, "CG Dependency [%d,%d)\n", begin, end);
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
        STDLOG(2, "Marking cellgroups_status = 1 for [%d,%d)\n", begin, end);
        #ifndef ONE_SIDED_GROUP_FINDING
            // In the two-sided case, there is only one Manifest.
            // We need to set cellgroups_status=2 for [end,end+2*GroupRadius]
            // This is so that CreateGlobalGroups() doesn't crash.
            // This isn't needed in the one-sided case; the dependencies 
            // take care of the need.
            for (int s=end; s<=end+2*GFC->GroupRadius; s++)
                GFC->cellgroups_status[CP->WrapSlab(s)]=2;
            STDLOG(2, "Marking cellgroups_status = 2 for [%d,%d]\n", end, end+2*GFC->GroupRadius);
        #endif
        return;
    }
};

int global_minslab_search;    // Used to pass an extra variable into partition

// Need some new Partition functions, to get ranges of slabs
// Will select things between [x,slab), where x=slab-CPD/2, handling the wrap
inline bool is_below_slab(ilstruct *particle, int slab) {
    int x = particle->newslab-slab;  // We want x<0, but not too much.
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
    AbacusMPILimiter mpi_limiter;

  public:
    ManifestCore m;   ///< The info we're sending over

    // Temporary storage buffers
    ilstruct *il;	///< Storage for the IL
    GroupLink *links;   ///< Storage for the GLL
    
    // Control logic
    enum ManifestStatus {
        MANIFEST_NOTREADY,
        MANIFEST_TRANSFERRING,
        MANIFEST_READY,
        MANIFEST_ALREADYIMPORTED,
    };

    enum ManifestStatus completed;	///< The current status of the manifest
    STimer Load;        ///< The timing for the Queue & Import blocking routines
    STimer Transmit;	///< The timing for the Send & Receive routines, usually non-blocking
    STimer CheckCompletion;        ///< The timing to check for completion
    STimer Communication;  ///< The timing of the main MPI communication time.
    size_t bytes;       ///< The number of bytes received

    MPI_Request *requests;    ///< A list of the non-blocking requests issued, each used once
    int numpending;        ///< The number of requests pending
    int maxpending;	   ///< The maximum number we're using

    int tag_offset;    ///< We'll add this to our MPI tags, just to keep them separate

//    void free_requests() {
//        assertf(numpending<=0, "We've been asked to free the MPI listing before all is completed, %d.\n", numpending);
//        if (requests!=NULL) delete[] requests;
//        numpending = 0;
//	max_pending = 0;
//    }


    #define MAX_REQUESTS 1024    // This is the most MPI work we can handle; hopefully very generous
    #define SIZE_MPI (1<<30)     // The maximum size of each MPI Isend.

    Manifest() : mpi_limiter(P.MPICallRateLimit_ms) {
    	m.numarenas = m.numil = m.numlinks = m.numdep = 0;
        #ifdef PARALLEL
            completed = MANIFEST_NOTREADY;
        #else
            completed = MANIFEST_ALREADYIMPORTED;   // We're not doing anything
        #endif
        bytes = 0;
        il = NULL;
        links = NULL;
        requests = new MPI_Request[MAX_REQUESTS];
        for (int j=0; j<MAX_REQUESTS; j++) requests[j]=MPI_REQUEST_NULL;
        numpending = 0;
        maxpending = 0;
        return;
    }
    ~Manifest() { 
        void *p;
        delete[] requests;
        // free_requests();
    }

    void set_tag(int j) {
        tag_offset = j*MAX_REQUESTS;
    }

    /// Allocate N things to track
//    void set_pending(int n) {
//        requests = new MPI_Request[n];
//        for (int j=0; j<n; j++) requests[j]=MPI_REQUEST_NULL;
//        numpending = n;
//        // STDLOG(1,"Establishing Manifest MPI listing of %d activities\n", numpending);
//    }

    inline void mark_as_done(int j) {
        assertf(requests[j]==MPI_REQUEST_NULL,"MPI_Request %d wasn't NULL\n", j);
        numpending--;
        // STDLOG(1,"Marked Manifest activity %d as done.  %d remain\n", j, numpending);
    }

    /// Return 1 if newly found to be done, 0 otherwise.
    /// This will only return 1 once; it's ok to call again, but will return 0
    inline int check_if_done(int j) {
        if (requests[j]!=MPI_REQUEST_NULL) {
            int err = 0, sent=0;
            #ifdef PARALLEL
            err = MPI_Test(requests+j,&sent,MPI_STATUS_IGNORE);
            #endif
            if (sent) { mark_as_done(j); return 1; }
        }
        return 0;
    }

    /// Return 1 if there is no more work
    inline int check_all_done() {
        if (numpending==0) return 1; 
        else return 0; 
    }

    /// This launches a send of N bytes, perhaps split into multiple Isend calls
    /// size must be in bytes; we only deal with bytes
    void do_MPI_Isend(void *ptr, uint64 size, int rank) {
        #ifdef PARALLEL
        bytes += size;
        while (size>0) {
            assertf(maxpending<MAX_REQUESTS, "Too many MPI requests %d\n", maxpending);
            int thissize = std::min(size, (uint64) SIZE_MPI);
            MPI_Isend(ptr, thissize, MPI_BYTE, rank, tag_offset+maxpending, comm_manifest,requests+maxpending);
            numpending++; maxpending++; size -= thissize; ptr = (char *)ptr+thissize;
        }
        #endif
    }

    /// This launches a receive of N bytes, perhaps split into multiple Irecv calls
    /// size must be in bytes; we only deal with bytes
    void do_MPI_Irecv(void *ptr, uint64 size, int rank) {
        #ifdef PARALLEL
        bytes += size;
        while (size>0) {
            assertf(maxpending<MAX_REQUESTS, "Too many MPI requests %d\n", maxpending);
            int thissize = std::min(size, (uint64) SIZE_MPI);
            MPI_Irecv(ptr, thissize, MPI_BYTE, rank, tag_offset+maxpending, comm_manifest,requests+maxpending);
            numpending++; maxpending++; size -= thissize; ptr = (char *)ptr+thissize;
        }
        #endif
    }


    inline int is_ready() { if (completed==MANIFEST_READY) return 1; return 0;}
	// Call this to see if the Manifest is ready to retrieve
    void done() { completed = MANIFEST_ALREADYIMPORTED; }

    void LoadArena(int type, int s) {
        ManifestArena *a = m.arenas+m.numarenas;
        a->type = type;
        a->slab = s;
        a->size = SB->SlabSizeBytes(type, s);
        a->ptr =  SB->GetSlabPtr(type, s);
        SB->MarkSlabUnavailable(type,s);   // Place this arena off limits
        STDLOG(3, "Queuing slab %d of type %d, size %l\n", s, type, a->size);
        m.numarenas++;
        assertf(m.numarenas<MAXMANIFEST, "numarenas has overflowed; increase MAXMANIFEST.");
        return;
    }

    // Here's the prototypes for the main routines
    void QueueToSend(int finished_slab);
    void Send();
    void FreeAfterSend();
    void Check();
    void SetupToReceive();
    void Receive();
    void ImportData();

};    


/// Here are our outgoing and incoming Manifest instances
/// We're going to have a sequence of these, which operate disjointly in time.  
/// We can just increment the pointer when one is done.
/// The original pointers are in the _vars; we keep these for deleting.
Manifest *SendManifest, *ReceiveManifest; 
Manifest *_SendManifest, *_ReceiveManifest; 
int nManifest;

/// Call this routine at the beginning of the timestep
void SetupManifest(int _nManifest) {
    nManifest = _nManifest+1;
        // We put on an extra one, just to avoid accidental overrunning.
    SendManifest = _SendManifest = new Manifest[nManifest];
    ReceiveManifest = _ReceiveManifest = new Manifest[nManifest];
    for (int j=0;j<nManifest;j++) {
        SendManifest[j].set_tag(j);
        ReceiveManifest[j].set_tag(j);
    }
    #ifdef PARALLEL
        assertf(MPI_size>1, "Can't run MPI-based manifest code with only 1 process.\n"); 
        // TODO: I don't see a way around this.  One ends up with the destination and source arenas being the same.
    #endif
    ReceiveManifest->SetupToReceive();
}

void FreeManifest() {
    STDLOG(2,"Freeing SendManifest\n");
    delete[] _SendManifest;
    STDLOG(2,"Freeing ReceiveManifest\n");
    delete[] _ReceiveManifest;
}

void CheckSendManifest() {
    for (int j=0; j<nManifest; j++) _SendManifest[j].FreeAfterSend();
    return;
}


// ================  Routine to define the outgoing information =======

/** The blocking routine to prepare information to Send

We've just finished the given slab.  
Load all of the information from lower slabs into the Manifest.
Then call the non-blocking communication and return.
*/

void Manifest::QueueToSend(int finished_slab) {
    #ifdef PARALLEL
    Load.Start();
    int cpd = P.cpd;

    // TODO: Consider some safety measures here.
    // E.g., we really should not have any slabs that have open GlobalGroupSlabs
    // in our wake.  Might check, since one could screw up the Dependencies.

    first_slab_finished = finished_slab;   // A global record of this
    m.remote_first_slab_finished = finished_slab;
    STDLOG(1,"Queueing the SendManifest at slab=%d\n", finished_slab);

    // Load the information from the Dependencies
    m.numdep = 0;
    m.dep[m.numdep++].Load(FetchSlabs, finished_slab, "FetchSlabs");
    m.dep[m.numdep++].Load(TransposePos, finished_slab, "TransposePos");
    m.dep[m.numdep++].Load(NearForce, finished_slab, "NearForce");
    m.dep[m.numdep++].Load(TaylorForce, finished_slab, "TaylorForce");
    m.dep[m.numdep++].Load(Kick, finished_slab, "Kick");
    m.dep[m.numdep++].Load(MakeCellGroups, finished_slab, "MakeCellGroups");
    m.dep[m.numdep++].Load(FindCellGroupLinks, finished_slab, "FindCellGroupLinks");
    int min_links_slab = m.dep[m.numdep-1].begin-1;
    	// We just determined that FindCellGroupLinks has executed on begin, so
	// the GLL might contain links including begin-1.
    m.dep[m.numdep++].Load(DoGlobalGroups, finished_slab, "DoGlobalGroups");
    m.dep[m.numdep++].Load(Output, finished_slab, "Output");
    m.dep[m.numdep++].Load(Microstep, finished_slab, "Microstep");
    m.dep[m.numdep++].Load(FinishGroups, finished_slab, "FinishGroups");
    m.dep[m.numdep++].Load(Drift, finished_slab, "Drift");
    int min_il_slab = m.dep[m.numdep-1].begin-FINISH_WAIT_RADIUS;
    	// We just determined that Drift has executed on begin, so
	// the rebinning might have taken particles to begin-FINISH_WAIT_RADIUS.
    m.dep[m.numdep++].Load(Finish, finished_slab, "Finish");
    m.dep[m.numdep++].Load(UnpackLPTVelocity, finished_slab, "UnpackLPTVelocity");
    m.dep[m.numdep++].LoadCG(finished_slab);
    	// LoadCG() includes moving info into the CellGroupArenas
    assertf(m.numdep<MAXDEPENDENCY, "m.numdep has overflowed its MAX value");

    STDLOG(2,"Queuing Arenas into the SendManifest\n");
    // Now load all of the arenas into the Manifest
    int min_slab = finished_slab;
    for (int type=0; type<NUMTYPES; type++) {
	    // Loop over all SlabTypes
	    for (int s=finished_slab-1; s>finished_slab-cpd; s--) {
	        // Check each trailing slab; if present, load it up
	        if (SB->IsSlabPresent(type,s) and not SB->IsOutputSlab(type) ) {
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
    STDLOG(2,"Done Queuing Arenas, spanning [%d,%d)\n", min_slab, finished_slab);

    STDLOG(2,"Queuing Insert List into the SendManifest, extracting [%d,%d)\n",
    	min_il_slab, finished_slab);
    // Partition the Insert List, malloc *il, and save it off
    IL->CollectGaps();   // Want to assure no MALgaps
    global_minslab_search = CP->WrapSlab(min_il_slab-finished_slab);
    uint64 mid = ParallelPartition(IL->list, IL->length, finished_slab, is_below_slab);

    m.numil = IL->length-mid;
    int ret = posix_memalign((void **)&il, PAGE_SIZE, sizeof(ilstruct)*m.numil);
    memcpy(il, IL->list+mid, sizeof(ilstruct)*m.numil);
	// Possible TODO: Consider whether this copy should be multi-threaded
    STDLOG(2, "Insert list had size %l, now size %l; sending %l\n", IL->length, mid, m.numil);
    IL->ShrinkMAL(mid);

    // Partition the GroupLink List, malloc *links, and save it off
    // TODO: Do these group finding variables always exist?
    if (GFC!=NULL) {
        STDLOG(2,"Queuing GroupLink List into the SendManifest, extracting [%d,%d)\n", min_links_slab, finished_slab);
        global_minslab_search = CP->WrapSlab(min_links_slab-finished_slab);
        mid = ParallelPartition(GFC->GLL->list, GFC->GLL->length, finished_slab, link_below_slab);
        ret = posix_memalign((void **)&links, PAGE_SIZE, sizeof(GroupLink)*(GFC->GLL->length-mid));
        m.numlinks = GFC->GLL->length-mid;
        memcpy(links, GFC->GLL->list+mid, sizeof(GroupLink)*m.numlinks);
            // Possible TODO: Consider whether this copy should be multi-threaded
        STDLOG(2, "Grouplink list had size %l, now size %l; sending %l\n", GFC->GLL->length, mid, m.numlinks);
        GFC->GLL->ShrinkMAL(mid);
    }
    Load.Stop();

    this->Send(); 
    // usleep(2e6);   // TODO: Don't forget to remove this
    #endif
    return;
}

// =============== Routine to actually transmit the outgoing info ====

void Manifest::Send() {
    #ifdef PARALLEL
    Transmit.Start();
    /// set_pending(m.numarenas+3);
    int rank = MPI_rank;
    rank--; if (rank<0) rank+=MPI_size;   // Now rank is the destination node
    STDLOG(1,"Will send the SendManifest to node rank %d\n", rank);
    // Send the ManifestCore
    STDLOG(2,"Isend Manifest Core\n");
    do_MPI_Isend(&m, sizeof(ManifestCore), rank);
    // Send all the Arenas
    for (int n=0; n<m.numarenas; n++) {
        STDLOG(2,"Isend Manifest Arena %d (slab %d of type %d) of size %d\n", 
            n, m.arenas[n].slab, m.arenas[n].type, m.arenas[n].size);
        do_MPI_Isend(m.arenas[n].ptr, m.arenas[n].size, rank);
    }
    // Send all the Insert List fragment
    STDLOG(2,"Isend Manifest Insert List of length %d\n", m.numil);
    do_MPI_Isend(il, sizeof(ilstruct)*m.numil, rank);
    if (GFC!=NULL) {
	// Send all the GroupLink List fragment
        STDLOG(2,"Isend Manifest GroupLink List of length %d\n", m.numlinks);
        do_MPI_Isend(links, sizeof(GroupLink)*m.numlinks, rank);
    } 
    // Victory!
    STDLOG(1,"Done queuing the SendManifest, %d total MPI parts, %l bytes\n", maxpending, bytes);
    completed = MANIFEST_TRANSFERRING;
    Transmit.Stop();
    Communication.Start();
    #endif
    return;
}

/// This is the routine to call frequently to try to clean up 
/// space after Send's have happened.
inline void Manifest::FreeAfterSend() {
    #ifdef PARALLEL
    if (completed != MANIFEST_TRANSFERRING) return;   // No active Isend's yet

    // Has it been long enough since our last time querying MPI?
    if(!mpi_limiter.Try())
        return;

    CheckCompletion.Start();

    /* OLDCODE
    check_if_done(0);    // Manifest Core
    if (check_if_done(1)) {
        free(il);    // Insert List
        STDLOG(1,"Freeing the Send Manifest Insert List, %d left\n", numpending);
    }
    if (check_if_done(2)) {
        free(links);    // GroupLink List; won't get called if already marked_as_done
        STDLOG(1,"Freeing the Send Manifest GroupLink List, %d left\n", numpending);
    }
    for (int n=0; n<m.numarenas; n++) 
        if (check_if_done(n+3)) {  // Arenas
            SB->DeAllocate(m.arenas[n].type, m.arenas[n].slab, 1);  // Deallocate, deleting any underlying file
            STDLOG(1,"Freeing the Send Manifest Arena, slab %d of type %d, %d left\n",
                m.arenas[n].slab, m.arenas[n].type, numpending);
        }
    END OLDCODE */

    for (int n=0; n<maxpending; n++) 
	if (check_if_done(n)) {
	    STDLOG(2,"Sent Manifest part %d, %d left\n", n, numpending);
	}
    if (check_all_done()) {
        // At present, we don't know which MPI send fragment maps to which arena, 
        // so we have to wait until all are sent to delete.
        Communication.Stop();
        STDLOG(2,"Send Manifest appears to be completely sent\n");
        free(il);    // Insert List
        if (GFC!=NULL) free(links);    // GroupLink List
        for (int n=0; n<m.numarenas; n++) {
            SB->DeAllocate(m.arenas[n].type, m.arenas[n].slab, 1);  // Deallocate, deleting any underlying file
            STDLOG(2,"Freeing the Send Manifest Arena, slab %d of type %d\n",
                m.arenas[n].slab, m.arenas[n].type);
        }
        completed=MANIFEST_READY;
        STDLOG(1,"Marking the Send Manifest as completely sent: %6.3f sec for %l bytes, %6.3f GB/sec\n", Communication.Elapsed(),
            bytes, bytes/(Communication.Elapsed()+1e-15)/(1024.0*1024.0*1024.0));
    }
    CheckCompletion.Stop();
    #endif
    return;
}


// =============== Routine to actually receive the incoming info ====

/// We have to issue the request to receive the ManifestCore
void Manifest::SetupToReceive() {
    #ifdef PARALLEL
    Transmit.Start();
    /// set_pending(1);
    int rank = MPI_rank;
    rank++; if (rank>=MPI_size) rank-=MPI_size;   // Now rank is the source node
    STDLOG(2,"Ireceive the Manifest Core from node rank %d\n", rank);
    do_MPI_Irecv(&m, sizeof(ManifestCore), rank);
    Transmit.Stop();
    #endif
    return;
}

/// This can be called in timestep.cpp to manually check 
/// whether Receive is ready to run.
inline void Manifest::Check() {
    #ifdef PARALLEL
    if (completed >= MANIFEST_READY) return;   // Nothing's active now

    // Has it been long enough since our last time querying MPI?
    if(!mpi_limiter.Try())
        return;
    
    CheckCompletion.Start();
    for (int n=0; n<maxpending; n++) 
	if (check_if_done(n)) {
	    STDLOG(2,"Received Manifest part %d, %d left\n", n, numpending);
	}
    CheckCompletion.Stop();
    if (check_all_done()) {
        if (completed == MANIFEST_TRANSFERRING) {
            STDLOG(2,"Marking the Receive Manifest as fully received\n");
            completed = MANIFEST_READY;
        }
        if (completed == MANIFEST_NOTREADY) {
            STDLOG(2,"Received the Manifest Core\n");
            ReleaseFreeMemoryToKernel();
            Receive();
            completed = MANIFEST_TRANSFERRING;
        } 
    }
    #endif
    return;
}

/* OLDCODE
    if (completed==0) {
        CheckCompletion.Start();
        if (check_if_done(0)) { 
            // The ManifestCore has arrived!  
            // We need to allocate the space and trigger Irecv's for the rest
            STDLOG(1,"Received the Manifest Core\n");
	    if (!check_all_done())
	        STDLOG(1,"Error: Receive Manifest should be done, but claims it isn't\n");
            CheckCompletion.Stop();
            Receive();
        } else CheckCompletion.Stop();
    }
    if (completed==1) {
        CheckCompletion.Start();
        for (int n=0; n<maxpending; n++) 
	    if (check_if_done(n)) {
		STDLOG(1,"Received Manifest part %d, %d left\n", n, numpending);
	    }
        if (check_all_done()) {
            completed=2;
            STDLOG(1,"Marking the Receive Manifest as fully received\n");
        }
        CheckCompletion.Stop();
    }
END OLDCODE */

/// This is the routine invoked by the communication thread.
/// It should store the results.
/// It allocates the needed space.
void Manifest::Receive() {
    #ifdef PARALLEL
    Transmit.Start();
    /// set_pending(m.numarenas+2);
    int rank = MPI_rank;
    MPI_Status retval;
    rank++; if (rank>=MPI_size) rank-=MPI_size;   // Now rank is the source node
    STDLOG(2,"Will receive the ReceiveManifest from node rank %d\n", rank);
    // Receive all the Arenas
    for (int n=0; n<m.numarenas; n++) {
        // TODO: is it true that none of the manifest slabs should land on ramdisk?
        SB->AllocateSpecificSize(m.arenas[n].type, m.arenas[n].slab, m.arenas[n].size, RAMDISK_NO);
        m.arenas[n].ptr = SB->GetSlabPtr(m.arenas[n].type, m.arenas[n].slab);
        //memset(m.arenas[n].ptr, 0, m.arenas[n].size);   // TODO remove
        STDLOG(2,"Ireceive Manifest Arena %d (slab %d of type %d) of size %d\n", 
            n, m.arenas[n].slab, m.arenas[n].type, m.arenas[n].size);
        do_MPI_Irecv(m.arenas[n].ptr, m.arenas[n].size, rank);
    }
    // Receive all the Insert List fragment
    int ret = posix_memalign((void **)&il, 64, m.numil*sizeof(ilstruct));
    assertf(ret == 0, "Failed to allocate %d bytes for receive manifest insert list\n", m.numil*sizeof(ilstruct));
    //memset(il, 0, sizeof(ilstruct)*m.numil);   // TODO remove
    STDLOG(2,"Ireceive Manifest Insert List of length %d\n", m.numil);
    do_MPI_Irecv(il, sizeof(ilstruct)*m.numil, rank);
    // Receive all the GroupLink List fragment
    if (GFC!=NULL) {
        ret = posix_memalign((void **)&links, 64, m.numlinks*sizeof(GroupLink));
        assertf(ret == 0, "Failed to allocate %d bytes for receive manifest group link list\n", m.numlinks*sizeof(GroupLink));
        //memset(links, 0, sizeof(GroupLink)*m.numlinks);   // TODO remove
        STDLOG(2,"Ireceive Manifest GroupLink List of length %d\n", m.numlinks);
        do_MPI_Irecv(links, sizeof(GroupLink)*m.numlinks, rank);
    } 
    // Victory!
    STDLOG(1,"Done issuing the ReceiveManifest, %d MPI parts, %l bytes\n", maxpending, bytes);
    Transmit.Stop();
    Communication.Start();
    #endif
    return;
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
    #ifdef PARALLEL
    assertf(completed==MANIFEST_READY, "ImportData has been called when completed==%d\n", completed);

    Communication.Stop();
    STDLOG(1,"Importing ReceiveManifest of %l bytes into the flow; took %6.3f sec, %6.3f GB/sec\n", bytes,
        Communication.Elapsed(), bytes/(Communication.Elapsed()+1e-15)/(1024.0*1024*1024));
    Load.Start();
    for (int n=0; n<m.numarenas; n++) {
	    SB->SetIOCompleted(m.arenas[n].type, m.arenas[n].slab);
	    STDLOG(3,"Completing Import of arena slab %d of type %d and size %l\n", 
	    	m.arenas[n].slab, m.arenas[n].type, m.arenas[n].size);
        if (m.arenas[n].type==PosSlab) {
            // Set the SlabSize based on the newly arrived PosSlab
            SS->setold(m.arenas[n].slab, m.arenas[n].size/sizeof(posstruct));
        }
    }

    // Set the dependencies.   Be careful that this keeps the Dependencies 
    // matched to the dep[] order set in QueueToSend().
    int n = 0;
    m.dep[n++].Set(FetchSlabs, "FetchSlabs");
    m.dep[n++].Set(TransposePos, "TransposePos");
    m.dep[n++].Set(NearForce, "NearForce");
    m.dep[n++].Set(TaylorForce, "TaylorForce");
    m.dep[n++].Set(Kick, "Kick");
    m.dep[n++].Set(MakeCellGroups, "MakeCellGroups");
    m.dep[n++].Set(FindCellGroupLinks, "FindCellGroupLinks");
    m.dep[n++].Set(DoGlobalGroups, "DoGlobalGroups");
    m.dep[n++].Set(Output, "Output");
    m.dep[n++].Set(Microstep, "Microstep");
    m.dep[n++].Set(FinishGroups, "FinishGroups");
    
    // We only want to lie about an extra Drift being completed on this node
    // if this is a Finish manifest and not a Group Finding manifest.
    // We don't distinguish explicitly between these two,
    // but the latter won't have any Drifts done, so we can check that
    if(m.dep[n].begin < m.dep[n].end)
        m.dep[n].end++;
        // We need to also mark the first_finished_slab on the other node
        // as having been completed.  This is because Finish requires this
        // as a precondition, and the particles incoming to Finish are on the
        // insert list.
    m.dep[n++].Set(Drift, "Drift");
    m.dep[n++].Set(Finish, "Finish");
    m.dep[n++].Set(UnpackLPTVelocity, "UnpackLPTVelocity");
    m.dep[n++].SetCG(m.remote_first_slab_finished);
    	// This will copy data back to GFC from CellGroupArenas
    assert(n==m.numdep);

    // Add *il to the insert list
    IL->CollectGaps();    // Want to assure we have no MALgaps
    uint64 len = IL->length;
    IL->GrowMAL(len+m.numil);
    STDLOG(2, "Growing IL list from %d by %d = %l\n", len, m.numil, IL->length);
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
        STDLOG(2, "Growing GroupLink list from %d by %d = %l\n", len, m.numlinks, GFC->GLL->length);
    }
    
    // We're done with this Manifest!
    done();
    Load.Stop();
    #endif
    return;
}

#endif    // NO_MPI
