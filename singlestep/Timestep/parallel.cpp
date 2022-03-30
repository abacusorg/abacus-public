/* file: parallel.cpp
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
