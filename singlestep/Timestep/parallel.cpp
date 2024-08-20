/* file: parallel.cpp
 *
 * This file is #included in proepi.cpp
 * and is responsible for initializing the vars in mpi_header.cpp.
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

        // The primary, "corporeal" domain boundaries will be registered to z=0
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

    // A stress test: go up, then down
    /*MERGE_GHOST_RADIUS = WriteState.FullStepNumber + 2;
    if(2*MERGE_GHOST_RADIUS > P.cpd/MPI_size_z){
        MERGE_GHOST_RADIUS = P.cpd/MPI_size_z - MERGE_GHOST_RADIUS;
    }*/

    // We will make the simplifying assumption of non-overlapping ghosts.
    // One could likely make this work as long as one thought about the MPI_size_z = 2 case
    // and avoided sending/deleting duplicates
    assertf(2*MERGE_GHOST_RADIUS <= node_z_size,
        "MERGE_GHOST_RADIUS (%d) big enough to overlap on node with primary size %d!",
        MERGE_GHOST_RADIUS, node_z_size
    );

    int zleft = MPI_rank_z - 1; if (zleft < 0) zleft += MPI_size_z;
    int zright = MPI_rank_z + 1; if (zright >= MPI_size_z) zright -= MPI_size_z;
    assertf(MERGE_GHOST_RADIUS <= all_node_z_size[zleft] && 
            MERGE_GHOST_RADIUS <= all_node_z_size[zright],
            "MERGE_GHOST_RADIUS (%d) bigger than neighbor primary regions!\n",
            MERGE_GHOST_RADIUS
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
