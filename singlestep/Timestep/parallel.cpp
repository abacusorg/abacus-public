/* file: parallel.cpp
 *
 * This file is #included in proepi.cpp
 * and is responsible for initializing the vars in mpi_header.cpp.
 *
 * It also defines the NeighborExchange routines.
 */

// A thin function to call MPI_Init() as early as possible,
// as recommended by the standard.
// MPI_rank, etc, will *not* be available until InitializeParallel()
void StartMPI() {
    // No STDLOG yet!

    #ifdef PARALLEL
        int init = 1;
        MPI_Initialized(&init);
        assertf(!init, "MPI was already initialized!\n");

        int ret = -1;
        MPI_Init_thread(NULL, NULL, MPI_THREAD_FUNNELED, &ret);
        assertf(ret>=MPI_THREAD_FUNNELED, "MPI_Init_thread() reports it supports level %d, not MPI_THREAD_FUNNELED.\n", ret);
    #endif
}


// Initialize the parallel topology: communicators, ranks, etc; but not the z domain.
// Requires the Parameters, but not the state.
void InitializeParallelTopology() {
    #ifdef PARALLEL
        MPI_Comm_rank(MPI_COMM_WORLD, &_world_rank);  // purely informational, we will not use _world_rank
        MPI_Comm_size(MPI_COMM_WORLD, &MPI_size);

        MPI_size_z = P.NumZRanks;
        MPI_size_x = MPI_size/MPI_size_z;
        assertf(MPI_size_x * MPI_size_z == MPI_size,
            "NumZRanks=%d does not divide MPI_size=%d evenly!", MPI_size_z, MPI_size);
        assertf(MPI_size_x > 1,
            "Must have at least two X ranks to run the parallel code!");

        // Establish the 2D topology
        int ndims = 2;
        int dims[] = {MPI_size_x, MPI_size_z};
        int periodic[] = {1, 1};
        int reorder = 1;
        MPI_Cart_create(MPI_COMM_WORLD, ndims, dims, periodic, reorder, &comm_2d);
        assertf(comm_2d != MPI_COMM_NULL, "Topology doesn't divide MPI size evenly?");
        
        // We let MPI_Cart_create() assign a rank that may be different from the COMM_WORLD rank!
        // But we'll only ever use the comm_2d rank
        MPI_Comm_rank(comm_2d, &MPI_rank);
        sprintf(NodeString,".%04d",MPI_rank);

        // Now establish this node's 2D rank
        int coords[ndims];
        MPI_Cart_coords(comm_2d, MPI_rank, ndims, coords);
        MPI_rank_x = coords[0];
        MPI_rank_z = coords[1];

        // Make a communicator across x
        int remain_dims_x[] = {1, 0};
        MPI_Cart_sub(comm_2d, remain_dims_x, &comm_1d_x);
        
        // Check that the x comm rank is the same as the x rank in the 2D topology
        int _comm_1d_x_rank = -1;
        MPI_Comm_rank(comm_1d_x, &_comm_1d_x_rank);
        assertf(_comm_1d_x_rank == MPI_rank_x,
            "Rank in comm_1d_x = %d, 2D x rank is %d\n", _comm_1d_x_rank, MPI_rank_x);

        // Make a communicator across z
        // TODO: unclear if one 2D communicator or multiple 1D communicators better
        int remain_dims_z[] = {0, 1};
        MPI_Cart_sub(comm_2d, remain_dims_z, &comm_1d_z);
        
        // Check that z comm rank is the same as the z rank in the 2D topology
        int _comm_1d_z_rank = -1;
        MPI_Comm_rank(comm_1d_z, &_comm_1d_z_rank);
        assertf(_comm_1d_z_rank == MPI_rank_z,
            "Rank in comm_1d_z = %d, 2D z rank is %d\n", _comm_1d_z_rank, MPI_rank_z);

        #ifdef MULTIPLE_MPI_COMM
            MPI_Comm_dup(comm_1d_x, &comm_taylors);
            MPI_Comm_dup(comm_1d_x, &comm_multipoles);
            MPI_Comm_dup(comm_1d_x, &comm_manifest);
            comm_global = comm_2d;
        #else
            comm_taylors = comm_1d_x;
            comm_multipoles = comm_1d_x;
            comm_manifest = comm_1d_x;
            comm_global = comm_2d;
        #endif

        MPI_Comm_dup(comm_1d_z, &comm_multipoles_z);
        MPI_Comm_dup(comm_1d_z, &comm_taylors_z);

        // Any custom MPI types
        MPI_Type_contiguous(sizeof(ilstruct), MPI_BYTE, &MPI_ilstruct);  // for 2D transfers off the IL
        MPI_Type_commit(&MPI_ilstruct);
    #endif
}

// Initialize `node_z_start` and friends, as well as ghost info
// Needs ReadState
void InitializeParallelDomain(){
    all_node_z_start = new int[MPI_size_z+1];
    all_node_z_size = new int[MPI_size_z];
    all_node_ky_start = new int[MPI_size_z+1];
    all_node_ky_size = new int[MPI_size_z];

    #ifdef PARALLEL
        // will always reflect the memory layout of the input slabs
        GHOST_RADIUS = ReadState.GhostRadius;

        if(MPI_size_z > 1){
            assertf(GHOST_RADIUS >= FORCE_RADIUS,
                "GHOST_RADIUS=%d not big enough for FORCE_RADIUS=%d\n",
                GHOST_RADIUS, FORCE_RADIUS);
                // InitGroupFinding() will check GHOST_RADIUS against GROUP_RADIUS
        }

        // The primary, "corporeal" domain boundaries will be registered to y=0
        // This is arbitrary, but it makes it so the multipoles never span the wrap
        for(int i = 0; i < MPI_size_z; i++){
            all_node_z_start[i] = P.cpd * i / MPI_size_z;
            all_node_z_size[i] = P.cpd * (i + 1) / MPI_size_z - all_node_z_start[i];
        }
        all_node_z_start[MPI_size_z] = P.cpd;
        node_z_start = all_node_z_start[MPI_rank_z];
        node_z_size = all_node_z_size[MPI_rank_z];
        
        // and compute the domain with ghosts
        node_z_start_ghost = node_z_start - GHOST_RADIUS;
        if (node_z_start_ghost < 0) node_z_start_ghost += P.cpd;
        node_z_size_with_ghost = node_z_size + 2*GHOST_RADIUS;

        assertf(node_z_start + node_z_size <= P.cpd, "Bad z split calculation?\n");  // A node won't span the wrap
    
    #else
    
        node_z_start = 0;
        node_z_size = P.cpd;
        node_z_start_ghost = 0;
        node_z_size_with_ghost = P.cpd;

        all_node_z_start[0] = 0;
        all_node_z_start[1] = P.cpd;
        all_node_z_size[0] = P.cpd;

    #endif

    if (MPI_size_z > 1){
        // After the multipoles y-FFT, each node gets a sub-domain of ky for all z
        all_node_ky_start[0] = 0;
        for(int i = 1; i < MPI_size_z + 1; i++){
            all_node_ky_start[i] = i*(P.cpd+1)/2/MPI_size_z;
            all_node_ky_size[i-1] = all_node_ky_start[i] - all_node_ky_start[i-1];
        }
        node_cpdp1half = all_node_ky_size[MPI_rank_z];
    } else {
        all_node_ky_start[0] = 0;
        all_node_ky_start[1] = P.cpd;
        all_node_ky_size[0] = P.cpd;
        node_cpdp1half = (P.cpd+1)/2;
    }
}

// Initialize the ghost radius for the merge slabs
// Has to be done after deciding that the next step will do group finding
void InitializeParallelMergeDomain(){
    #ifdef PARALLEL
        if(MPI_size_z > 1){
            MERGE_GHOST_RADIUS = P.NearFieldRadius;
            if(WriteState.DoGroupFindingOutput && WriteState.VelIsSynchronous){
                MERGE_GHOST_RADIUS = std::max(MERGE_GHOST_RADIUS, 2*P.GroupRadius+1);
            }
        } else {
            MERGE_GHOST_RADIUS = 0;
        }

        // we could track a set of indices like "node_z_merge_start",
        // but if those are *only* used by the merge, we can keep them
        // local to the merge

    int zleft = MPI_rank_z - 1; if (zleft < 0) zleft += MPI_size_z;
    int zright = MPI_rank_z + 1; if (zright >= MPI_size_z) zright -= MPI_size_z;
    assertf(MERGE_GHOST_RADIUS <= all_node_z_size[zleft] && 
            MERGE_GHOST_RADIUS <= all_node_z_size[zright],
            "MERGE_GHOST_RADIUS (%d) bigger than neighbor primary regions!\n"
            );

    WriteState.GhostRadius = MERGE_GHOST_RADIUS;
    #endif
}


void LogParallelTopology(){
    // once we have STDLOG

    #ifdef PARALLEL
        STDLOG(0,"Initialized parallel topology.\n");
        STDLOG(0,"1D node rank %d of %d total\n", MPI_rank, MPI_size);
        STDLOG(0,"2D (x,z) node rank is (%d,%d) of (%d,%d) total\n",
            MPI_rank_x, MPI_rank_z, MPI_size_x, MPI_size_z);
        STDLOG(2,"Original 1D rank was %d, remapped to %d\n", _world_rank, MPI_rank);

        STDLOG(0,"This node's primary domain is z=[%d,%d) (%d columns)\n",
            node_z_start, node_z_start + node_z_size, node_z_size);
        int ghostend = node_z_start_ghost + node_z_size_with_ghost;
        if (ghostend >= P.cpd) ghostend -= P.cpd;
        STDLOG(0,"Including ghost: z=[%d,%d) (%d columns)\n",
            node_z_start_ghost, ghostend, node_z_size_with_ghost);

        STDLOG(0,"The ReadState has GHOST_RADIUS=%d ghost columns\n", GHOST_RADIUS);
        STDLOG(0,"The WriteState will have MERGE_GHOST_RADIUS=%d\n", MERGE_GHOST_RADIUS);
    #endif

    char hostname[1024];
    gethostname(hostname,1024);
    STDLOG(0,"Host machine name is %s\n", hostname);
}


void FinalizeParallel() {
    #ifdef PARALLEL

        STDLOG(1, "Tearing down MPI communicators\n");

        #ifdef MULTIPLE_MPI_COMM

        MPI_Comm_free(&comm_taylors);
        MPI_Comm_free(&comm_multipoles);
        MPI_Comm_free(&comm_manifest);

        #endif

        MPI_Comm_free(&comm_2d);
        MPI_Comm_free(&comm_1d_x);
        MPI_Comm_free(&comm_1d_z);

        // Finalize MPI
        STDLOG(0,"Calling MPI_Finalize()\n");
        MPI_Finalize();
    #else
    #endif

    delete[] all_node_z_start;
    delete[] all_node_z_size;
    delete[] all_node_ky_start;
    delete[] all_node_ky_size;
}

/* Begin Neighbor Exchange routines */

#ifdef PARALLEL
// Partition functions
inline bool is_right_ghost(ilstruct *particle, int slab, int first_ghost){
    
    if(particle->newslab != slab)
        return 0;

    int cpdm1half = (P.cpd - 1)/2;

    int cellz = particle->global_cellz();

    // With just two splits, one node will be bigger than half the domain!
    if(MPI_size_z == 2 && cellz == node_z_start)
        return 0;

    int dist = cellz - first_ghost;
    if (dist > cpdm1half){
        dist -= P.cpd;
    } else if (dist < -cpdm1half){
        dist += P.cpd;
    }

    return dist >= 0;
}

inline bool is_left_ghost(ilstruct *particle, int slab, int first_ghost){
    if(particle->newslab != slab)
        return 0;

    int cpdm1half = (P.cpd - 1)/2;

    int cellz = particle->global_cellz();

    // With just two splits, one node will be bigger than half the domain!
    if(MPI_size_z == 2 && cellz == (node_z_start + node_z_size - 1))
        return 0;

    int dist = cellz - first_ghost;
    if (dist > cpdm1half){
        dist -= P.cpd;
    } else if (dist < -cpdm1half){
        dist += P.cpd;
    }

    return dist <= 0;
}

/* NeighborExchange happens when a slab and both adjacent slabs have Drifted,
   and does the following for each neighbor rank:
 * 1. Partition the IL to find particles the neighbor needs (for its primary or ghost)
 * 2. Copy them to a buffer
 * 3. Initiate an MPI_Isend
 * 4. Partition the tail of the IL, this time for particles that are *only* in the neighbor's primary
 * 5. "Delete" them
 * 
 * And then later...
 * 6. Check if we have an incoming MPI message
 * 7. Allocate space to receive the message
 * 8. Start the MPI_Irecv
 * 
 * And then later...
 * 9. Check if the receive is done
 * 10. Push the buffer onto the insert list
 * 11. Free the receive buffer and release the MPI handle
 *
 * And sometime after the send...
 * 12. Check if the send is complete
 * 13. Free the send buffer and release the MPI handle
 */

class NeighborExchanger {
/* Represents a z exchange with one neighbor. */
private:
    int slab;
    int zneigh;  // neighbor z rank
    int right;   // sending the left or right side of our domain

    int send_status = 0;  // 0 = not sent; 1 = sending; 2 = done
    ilstruct *sendbuf = NULL;  // staging the MPI_Isend
    size_t sendelem = 0;    // number of ilstructs

    int recv_status = 0;    // 0 = waiting; 1 = receiving; 2 = done
    ilstruct *recvbuf = NULL;  // staging the MPI_Irecv
    size_t recvelem = 0;

    MPI_Request send_handle = MPI_REQUEST_NULL;
    MPI_Request recv_handle = MPI_REQUEST_NULL;

    // If MPI_size_z=2, then the right and left neighbor will be the same.
    // So we'll need separate tag "namespaces".
    static const int RIGHT_TAG_OFFSET = 100000;
    static const int LEFT_TAG_OFFSET = 200000;

public:
    NeighborExchanger(int slab, int zneigh, int right) :
        slab(slab), zneigh(zneigh), right(right) { }

    
    void send(){
        // 1. partition
        ilstruct *start = partition(&sendelem);

        // 2. copy
        copy_to_sendbuf(start, sendelem);

        // 3. MPI_Isend
        launch_mpi_send();

        // 4. tail partition
        size_t removed = tail_partition(start, sendelem);

        // 5. delete
        IL->ShrinkMAL(IL->length - removed);
    }

    int try_free_sendbuf(){
        int ret = 0;
        if(send_status == 1){
            int done = 0;
            // 12. done send?
            MPI_Test(&send_handle, &done, MPI_STATUS_IGNORE);
            if (done) {
                STDLOG(1, "Done ilstruct send to zrank %d for slab %d\n", zneigh, slab);

                // 13. free
                free(sendbuf);
                sendbuf = NULL;
                send_status = 2;
                ret = 1;
            }
        }
        return ret;
    }

    int done_send(){
        return send_status == 2;
    }

private:
    ilstruct *partition(size_t *nhigh){
        IL->CollectGaps();
        
        uint64 mid;
        if(right){
            // [0..mid) are not in slab, [mid..length) are in slab
            mid = ParallelPartition(IL->list, IL->length,
                is_right_ghost,
                slab,
                node_z_start + node_z_size - MERGE_GHOST_RADIUS);
        }
        else {
            mid = ParallelPartition(IL->list, IL->length,
                is_left_ghost,
                slab,
                node_z_start + MERGE_GHOST_RADIUS - 1);
        }
        
        *nhigh = (size_t) (IL->length - mid);
        return IL->list + mid;
    }


    size_t tail_partition(ilstruct *start, size_t nelem){
        // No need to collect gaps; still contiguous from partition()
        
        /* Normal steps will not need to execute this.
           This only does anything when the IL has particles that are so
           far away from our domain that we will not need them as ghosts.
           Since we can only drift 1 cell under normal conditions, and
           MERGE_GHOST_RADIUS is usually 2, this usually cannot happen.
           Common exceptions are the ICs (which have some disorder),
           and MERGE_GHOST_RADIUS=0 (but we probably never do that).
        */

        if(WriteState.FullStepNumber != 0 &&
            WriteState.LPTStepNumber == 0 &&
            MERGE_GHOST_RADIUS > 0)
                return 0;

        uint64 mid;
        if(right){
            mid = ParallelPartition(start, nelem,
                is_right_ghost,
                slab,
                node_z_start + node_z_size + MERGE_GHOST_RADIUS);
                // only delete particles that are neighbor's primary
                // and we do not need as ghosts
        }
        else {
            mid = ParallelPartition(start, nelem,
                is_left_ghost,
                slab,
                node_z_start - MERGE_GHOST_RADIUS - 1);
        }

        size_t nhigh = (size_t) nelem - mid;
        STDLOG(2,"Removing %d particles that are strictly neighbor's primary\n", nhigh);
        return nhigh;
    }

    void copy_to_sendbuf(ilstruct *start, size_t nelem){
        size_t sz = sizeof(ilstruct)*nelem;
        posix_memalign((void **) &(this->sendbuf), PAGE_SIZE, sz);
        assertf(this->sendbuf != NULL,
            "Failed neighbor exchange sendbuf alloc of %d bytes\n", sz);

        #pragma omp parallel for schedule(static)
        for(size_t i = 0; i < nelem; i++){
            this->sendbuf[i] = start[i];
        }
    }

    void launch_mpi_send(){
        assertf(send_status == 0,
            "send() called but exchange already sent to zrank %d, slab %d?\n", zneigh, slab);

        assertf(sendelem < INT32_MAX,
            "2D neighbor exchanger trying to send %d particles, which overflows 32-bit int\n",
            sendelem);

        STDLOG(1, "Sending %d ilstructs (%.4g MB) to zrank %d for slab %d\n", sendelem, sendelem*sizeof(ilstruct)/1e6, zneigh, slab);
        // each pair of ranks has only one send and one recv per slab
        int tag = right ? (slab + RIGHT_TAG_OFFSET) : (slab + LEFT_TAG_OFFSET);
        MPI_Isend((const void *) sendbuf, (int) sendelem, MPI_ilstruct, zneigh, tag, comm_1d_z, &send_handle);

        send_status = 1;  // sending
    }

    void push_buffer_to_IL(){
        // Bypass IL->Push, which builds new ilstructs.
        // We just need to do a raw copy!
        IL->CollectGaps();
        uint64 oldlen = IL->length;
        IL->GrowMAL(oldlen + recvelem);

        #pragma omp parallel for schedule(static)
        for(uint64 i = 0; i < recvelem; i++){
            // Adjust sorting key to be relative to the node-local ghost_z_start.
            // TODO: could probably emplace, copying all fields except zlocal
            ilstruct p = recvbuf[i];
            p.setzlocal( CP->WrapSlab(p.local_cellz() - (node_z_start - all_node_z_start[zneigh])) );
            IL->list[oldlen + i] = p;
        }
    }

public:
    int try_receive(){
        int ret = 0;
        if (send_status > 0      // started send
            && recv_status == 1  // receiving
            ){
            // Not allowed to import recvbuf until we send! Otherwise we will send things we received

            int done = 0;
            // 9. done receive?
            MPI_Test(&recv_handle, &done, MPI_STATUS_IGNORE);
            if (done) {
                STDLOG(1, "Done ilstruct receive from zrank %d for slab %d\n", zneigh, slab);
                // 10. push to IL
                push_buffer_to_IL();

                // 11. free
                free(recvbuf);
                recvbuf = NULL;
                recv_status = 2;
                ret = 1;
            }
            return ret;
        }

        if (recv_status == 2)
            return ret;  // all done!

        if (recv_status == 0){
            // 6. check for incoming
            MPI_Status status;
            // if we sent to the right, receive from the left
            int tag = right ? (slab + LEFT_TAG_OFFSET) : (slab + RIGHT_TAG_OFFSET);
            int ready = 0;
            MPI_Iprobe(zneigh, tag, comm_1d_z, &ready, &status);
            if(ready){
                // 7. alloc recvbuf
                int _recvelem32;
                MPI_Get_count(&status, MPI_ilstruct, &_recvelem32);
                recvelem = (size_t) _recvelem32;
                size_t sz = sizeof(ilstruct)*recvelem;
                posix_memalign((void **) &recvbuf, PAGE_SIZE, sz);
                assertf(recvbuf != NULL, "Failed neighbor exchange recvbuf alloc of %d bytes\n", sz);
                
                // 8. start recv
                STDLOG(1, "Receiving %d ilstructs (%.4g MB) from zrank %d for slab %d\n",
                    recvelem, recvelem*sizeof(ilstruct)/1e6, zneigh, slab);
                MPI_Irecv(recvbuf, (int) recvelem, MPI_ilstruct, zneigh, tag, comm_1d_z, &recv_handle);
                recv_status = 1;
                ret = 1;
            }
            return ret;
        }

        assertf(recv_status == 1, "Unknown recv_status %d\n", recv_status);
        return ret;
    }

    int done_receive(){
        return recv_status == 2;
    }

    void force_done(){
        assertf(recv_status == 0, "Tried to force done an active Neighbor Exchange?");
        assertf(send_status == 0, "Tried to force done an active Neighbor Exchange?");

        recv_status = 2;
        send_status = 2;
    }
};


NeighborExchanger **left_exchanger;
NeighborExchanger **right_exchanger;
int neighbor_exchange_is_noop;

// Lightweight setup of data structures
void SetupNeighborExchange(int first, int nslab){
    if(MPI_size_z == 1){
        neighbor_exchange_is_noop = 1;
        STDLOG(1,"Skipping Neighbor Exchange setup\n");
        return;
    }

    STDLOG(1,"Setting up Neighbor Exchange to do exchanges on slabs [%d,%d)\n",
        first, CP->WrapSlab(first+nslab));

    left_exchanger = new NeighborExchanger*[P.cpd]();
    right_exchanger = new NeighborExchanger*[P.cpd]();

    int leftz = MPI_rank_z - 1;
    if(leftz < 0) leftz += MPI_size_z;
    
    int rightz = MPI_rank_z + 1;
    if(rightz >= MPI_size_z) rightz -= MPI_size_z;

    for(int i = first; i < first+nslab; i++){
        int iw = CP->WrapSlab(i);
        left_exchanger[iw] = new NeighborExchanger(iw, leftz, 0);
        right_exchanger[iw] = new NeighborExchanger(iw, rightz, 1);
    }
    // slabs not on this node have NULL exchangers
}


// Called immediately after drift
void DoNeighborSend(int slab){
    if(neighbor_exchange_is_noop)
        return;

    STDLOG(1,"Triggering Neighbor Exchange sends on slab %d\n", slab);

    assertf(left_exchanger[slab] != NULL,
        "Slab %d not set up to do neighbor send?\n", slab);
    assertf(right_exchanger[slab] != NULL,
        "Slab %d not set up to do neighbor send?\n", slab);

    left_exchanger[slab]->send();
    right_exchanger[slab]->send();
}

// Called repeatedly in the timestep loop
int AttemptNeighborReceive(int first, int receive_ahead){
    // This will typically be called with
    // first = (Finish.last_slab_executed - FINISH_WAIT_RADIUS), receive_ahead ~ 3
    // It's also safe to call with (0,cpd) if one always wants to receive everything
    int ret = 0;

    if(neighbor_exchange_is_noop){
        return ret;  // safely no-op
    }

    // Try to receive data from one or more slabs
    for(int i = first; i < first + receive_ahead; i++){
        int iw = CP->WrapSlab(i);
        // try_receive() will run twice for every exchange: once to launch the receive,
        // and again to install the data
        if(left_exchanger[iw] != NULL) ret |= left_exchanger[iw]->try_receive();
        if(right_exchanger[iw] != NULL) ret |= right_exchanger[iw]->try_receive();
    }

    // Do a lightweight release of the send buffer and MPI handle
    for(int i = 0; i < P.cpd; i++){
        if(left_exchanger[i] != NULL){
            ret |= left_exchanger[i]->try_free_sendbuf();
        }
        if(right_exchanger[i] != NULL){
            ret |= right_exchanger[i]->try_free_sendbuf();
        }
    }

    return ret;
}

// Checked in the Finish precondition
int IsNeighborReceiveDone(int slab){
    if(neighbor_exchange_is_noop){
        return 1;  // safely no-op
    }

    int sw = CP->WrapSlab(slab);

    assertf(left_exchanger[sw] != NULL,
        "Checking neighbor receive status on slab %d not set up to receive?\n", sw);
    assertf(right_exchanger[sw] != NULL,
        "Checking neighbor receive status on slab %d not set up to receive?\n", sw);

    return left_exchanger[sw]->done_receive() && right_exchanger[sw]->done_receive();
}


// Called in the Epilogue
void TeardownNeighborExchange(){
    if(neighbor_exchange_is_noop){
        return;  // safely no-op
    }

    STDLOG(3,"Tearing down Neighbor Exchange\n");
    
    for(int i = 0; i < P.cpd; i++){
        if(left_exchanger[i] != NULL){
            assertf(left_exchanger[i]->done_send(), "left exchanger for slab %d not done send?\n", i);
            assertf(left_exchanger[i]->done_receive(), "left exchanger for slab %d not done receive?\n", i);
            delete left_exchanger[i];
        }
        if(right_exchanger[i] != NULL){
            assertf(right_exchanger[i]->done_send(), "right exchanger for slab %d not done send?\n", i);
            assertf(right_exchanger[i]->done_receive(), "right exchanger for slab %d not done receive?\n", i);
            delete right_exchanger[i];
        }
    }
    delete[] left_exchanger;
    delete[] right_exchanger;
}

#else // PARALLEL

void SetupNeighborExchange(int first, int nslab) { }
int AttemptNeighborReceive(int first, int receive_ahead){ return 0; }
void DoNeighborSend(int slab){ }
int IsNeighborReceiveDone(int slab){ return 1; }
void TeardownNeighborExchange(){ }

#endif
