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

        STDLOG(0,"Initialized MPI.\n");
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
            MPI_Comm_dup(comm_2d, &comm_taylors);
            MPI_Comm_dup(comm_2d, &comm_multipoles);
            MPI_Comm_dup(comm_1d_x, &comm_manifest);
            MPI_Comm_dup(comm_2d, &comm_global);
        #else
            comm_taylors = comm_2d;
            comm_multipoles = comm_2d;
            comm_manifest = comm_1d_x;
            comm_global = comm_2d;
        #endif

        // Any custom MPI types
        MPI_Type_contiguous(sizeof(ilstruct), MPI_BYTE, &MPI_ilstruct);  // for 2D transfers off the IL
        MPI_Type_commit(&MPI_ilstruct);
    #endif
}

// Initialize `node_z_start` and friends, as well as ghost info
// Needs ReadState
void InitializeParallelDomain(){
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
        node_z_start = P.cpd * MPI_rank_z / MPI_size_z;
        node_z_size = P.cpd * (MPI_rank_z + 1) / MPI_size_z - node_z_start;
        
        node_z_start_ghost = node_z_start - GHOST_RADIUS;
        if (node_z_start_ghost < 0) node_z_start_ghost += P.cpd;
        node_z_size_with_ghost = node_z_size + 2*GHOST_RADIUS;

        assertf(node_z_start + node_z_size <= P.cpd, "Bad z split calculation?");  // A node won't span the wrap

    #else
    
        node_z_start = 0;
        node_z_size = P.cpd;
        node_z_start_ghost = 0;
        node_z_size_with_ghost = P.cpd;
    
    #endif
}

// Initialize the ghost radius for the merge slabs
// Has to be done after deciding that the next step will do group finding
void InitializeParallelMergeDomain(){
    #ifdef PARALLEL
        if(MPI_size_z > 1){
            MERGE_GHOST_RADIUS = P.NearFieldRadius;  // TODO: this will be std::max(P.NearFieldRadius, GROUP_RADIUS_NEXT_STEP)
        } else {
            MERGE_GHOST_RADIUS = 0;
        }

        // we could track a set of indices like "node_z_merge_start",
        // but if those are *only* used by the merge, we can keep them
        // local to the merge
    #endif

    WriteState.GhostRadius = MERGE_GHOST_RADIUS;
}


void LogParallelTopology(){
    // once we have STDLOG

    #ifdef PARALLEL
        STDLOG(0,"Initialized parallel topology.\n");
        STDLOG(0,"1D node rank %d of %d total\n", MPI_rank, MPI_size);
        STDLOG(0,"2D (x,z) node rank is (%d,%d) of (%d,%d) total\n",
            MPI_rank_x, MPI_rank_z, MPI_size_x, MPI_size_z);
        STDLOG(2,"Original 1D rank was %d, remapped to %d\n", _world_rank, MPI_rank);

        STDLOG(0,"This node will do %d columns in z=[%d,%d)\n",
            node_z_size, node_z_start, node_z_start + node_z_size);
        STDLOG(0,"z domain including ghost is %d columns in z=[%d,%d)\n",
            node_z_size_with_ghost, node_z_start_ghost, node_z_start_ghost + node_z_size_with_ghost);

        STDLOG(0,"The ReadState has GHOST_RADIUS=%d\n ghost columns\n", GHOST_RADIUS);
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
        MPI_Comm_free(&comm_global);

        #endif

        MPI_Comm_free(&comm_2d);
        MPI_Comm_free(&comm_1d_x);
        MPI_Comm_free(&comm_1d_z);

        // Finalize MPI
        STDLOG(0,"Calling MPI_Finalize()\n");
        MPI_Finalize();
    #else
    #endif
}

/* Begin Neighbor Exchange routines */

// Partition functions
inline bool is_right_ghost(ilstruct *particle, int first_ghost){
    // TODO: this modulus could be slow. First let's check if it is;
    // if so, we could store z in the ilstruct.
    // Or fix CPD at compile time and use libdivide!

    // k = y*cpd + z, so k % cpd gives z
    return (particle->k % P.cpd) >= first_ghost;
}

inline bool is_left_ghost(ilstruct *particle, int first_ghost){
    return (particle->k % P.cpd) <= first_ghost;
}

/* NeighborExchange happens immediately after the Drift, and does the following for each neighbor:
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
 */

/* Represents a z exchange with one neighbor.
 */
class NeighborExchanger {
private:
    int slab;
    int zneigh;  // neighbor z rank
    int right;   // sending the left or right side of our domain

    ilstruct *sendbuf;  // staging the MPI_Isend
    size_t sendelem;    // number of ilstructs

    int recv_status;    // 0 = waiting; 1 = receiving; 2 = done
    ilstruct *recvbuf;  // staging the MPI_Irecv
    size_t recvelem;

    MPI_Request send_handle;
    MPI_Request recv_handle;

public:
    NeighborExchanger(int slab, int zneigh, int right){
        this->slab = slab;
        this->zneigh = zneigh;
        this->right = right;

        this->recv_status = 0;
        this->recvbuf = NULL;
        this->recvelem = 0;
        this->recv_handle = MPI_REQUEST_NULL;

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

private:
    ilstruct *partition(size_t *nhigh){
        IL->CollectGaps();
        
        uint64 mid;
        if(right){
            // [0..mid) are not in slab, [mid..length) are in slab
            mid = ParallelPartition(IL->list, IL->length,
                node_z_start + node_z_size - MERGE_GHOST_RADIUS,
                is_right_ghost);
        }
        else {
            mid = ParallelPartition(IL->list, IL->length,
                node_z_start + MERGE_GHOST_RADIUS - 1,
                is_left_ghost);
        }
        
        *nhigh = (size_t) (IL->length - mid);
        return IL->list + mid;
    }

    size_t tail_partition(ilstruct *start, size_t nelem){
        // No need to collect gaps!
        
        uint64 mid;
        if(right){
            mid = ParallelPartition(start, nelem,
                node_z_start + node_z_size,
                is_right_ghost);
        }
        else {
            mid = ParallelPartition(start, nelem,
                node_z_start - 1,
                is_left_ghost);
        }
        
        size_t nhigh = (size_t) nelem - mid;
        return nhigh;
    }

    void copy_to_sendbuf(ilstruct *start, size_t nelem){
        #pragma omp parallel for schedule(static)
        for(size_t i = 0; i < nelem; i++){
            this->sendbuf[i] = start[i];
        }
    }

    void launch_mpi_send(){
        assertf(sendelem < INT32_MAX,
            "2D neighbor exchanger trying to send %d particles, which overflows 32-bit int\n",
            sendelem);

        STDLOG(1, "Sending %d ilstructs (%.3g GB) to rank %d for slab %d\n", sendelem, sendelem*sizeof(ilstruct), zneigh, slab);
        int tag = slab;  // each pair of ranks has only one send and one recv per slab
        MPI_Isend((const void *) sendbuf, (int) sendelem, MPI_ilstruct, zneigh, tag, comm_1d_z, &send_handle);
    }

    void push_buffer_to_IL(){
        // Bypass IL->Push, which builds new ilstructs.
        // We just need to do a raw copy!
        IL->CollectGaps();
        uint64 oldlen = IL->length;
        IL->GrowMAL(oldlen + recvelem);

        #pragma omp parallel for schedule(static)
        for(uint64 i = 0; i < recvelem; i++){
            IL->list[oldlen + i] = recvbuf[i];
        }
    }

public:
    int try_mpi_recv(){
        if (recv_status == 1){
            int done = 0;
            // 9. done receive?
            MPI_Test(&recv_handle, &done, MPI_STATUS_IGNORE);
            if (done) {
                STDLOG(1, "Done ilstruct receive from rank %d for slab %d\n", recvelem, recvelem*sizeof(ilstruct), zneigh, slab);
                // 10. push to IL
                push_buffer_to_IL();

                // 11. free
                free(recvbuf);
                recvbuf = NULL;
                recv_status = 2;
            }
        }
        if (recv_status == 2)
            return 1;  // all done!

        assertf(recv_status == 0, "Unknown recv_status %d\n", recv_status);

        // 6. check for incoming
        MPI_Status status;
        int tag = slab;
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
            STDLOG(1, "Receiving %d ilstructs (%.3g GB) from rank %d for slab %d\n", recvelem, recvelem*sizeof(ilstruct), zneigh, slab);
            MPI_Irecv(recvbuf, (int) recvelem, MPI_ilstruct, zneigh, tag, comm_1d_z, &recv_handle);
            recv_status = 1;
        }
    }
};

NeighborExchanger **left_exchanger;
NeighborExchanger **right_exchanger;


void StartNeighborExchange(int slab){
    if(left_exchanger == NULL){
        left_exchanger = new NeighborExchanger*[P.cpd]();
    }
    if(right_exchanger == NULL){
        right_exchanger = new NeighborExchanger*[P.cpd]();
    }

    // The construction of these objects starts the MPI work

    int leftz = MPI_rank_z - 1;
    if(leftz < 0) leftz += MPI_size_z;
    left_exchanger[slab] = new NeighborExchanger(slab, leftz, 0);

    int rightz = MPI_rank_z + 1;
    if(rightz >= MPI_size_z) rightz -= MPI_size_z;
    right_exchanger[slab] = new NeighborExchanger(slab, rightz, 1);
}


void TeardownNeighborExchange(){
    for(int i = 0; i < P.cpd; i++){
        assertf(left_exchanger[i] == NULL, "Left exchanger %d still present?\n", i);
        assertf(right_exchanger[i] == NULL, "Right exchanger %d still present?\n", i);
    }
    delete[] left_exchanger;
    delete[] right_exchanger;
}
