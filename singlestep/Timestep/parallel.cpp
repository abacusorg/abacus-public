// file: parallel.cpp
// This file is #included in proepi.cpp

// These are the logical (i.e. global) indices that track the z-range belonging to this node

// TODO: vars might go to parallel.h
int node_z_start;
int node_z_size;
int node_z_start_ghost;
int node_z_size_with_ghost;

char NodeString[8] = "";     // Set to "" for serial, ".NNNN" for MPI
int MPI_size = 1, MPI_rank = 0;     // We'll set these globally, so that we don't have to keep fetching them

// The 2D rank and size
int MPI_rank_x = 0, MPI_size_x = 1;
int MPI_rank_z = 0, MPI_size_z = 1;

int GHOST_RADIUS = 0;        // in the read state
int MERGE_GHOST_RADIUS = 0;  // in the write state

// A thin function to call MPI_Init as early as possible,
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


// Initialize the parallel state: topology, communicators, ranks, etc.
// This is separated from StartMPI() because it requires that we have read
// the parameter file and done some other work best done after MPI_Init()
void InitializeParallel() {
    #ifdef PARALLEL
        int _world_rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &_world_rank);  // purely informational
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
        MPI_Cart_coords(comm_2d, MPI_rank, ndims, &coords);
        MPI_rank_x = coords[0];
        MPI_rank_z = coords[1];

        #ifdef MULTIPLE_MPI_COMM
            MPI_Comm_dup(comm_2d, &comm_taylors);
            MPI_Comm_dup(comm_2d, &comm_multipoles);
            MPI_Comm_dup(comm_2d, &comm_manifest);  // TODO: does this really need to be an all-to-all comm?
            MPI_Comm_dup(comm_2d, &comm_global);
        #else
            comm_taylors = comm_2d;
            comm_multipoles = comm_2d;
            comm_manifest = comm_2d;
            comm_global = comm_2d;
        #endif

        if(MPI_size_z > 1){
            GHOST_RADIUS = ReadState.GhostRadius;
            assertf(GHOST_RADIUS >= FORCE_RADIUS,
                "GHOST_RADIUS=%d not big enough for FORCE_RADIUS=%d\n",
                GHOST_RADIUS, FORCE_RADIUS);
            assertf(GHOST_RADIUS >= GROUP_RADIUS,
                "GHOST_RADIUS=%d not big enough for GROUP_RADIUS=%d\n",
                GHOST_RADIUS, GROUP_RADIUS);
            MERGE_GHOST_RADIUS = FORCE_RADIUS;  // TODO: this will be std::max(FORCE_RADIUS, GROUP_RADIUS_NEXT_STEP)
        } else {
            GHOST_RADIUS = 0;
            MERGE_GHOST_RADIUS = 0;
        }

        // The primary, "corporeal" domain boundaries will be registered to y=0
        // This is arbitrary, but it makes it so the multipoles never span the wrap
        node_z_start = P.cpd * MPI_rank_z / MPI_size_z;
        node_z_size = P.cpd * (MPI_rank_z + 1) / MPI_size_z - node_z_start;
        
        node_z_start_ghost = node_z_start - GHOST_RADIUS;
        if (node_z_start_ghost < 0) node_z_start_ghost += P.cpd;
        node_z_size_with_ghost = node_z_size + 2*GHOST_RADIUS;

        assertf(node_z_start + node_z_size <= P.cpd, "Bad z split calculation?");  // A node can't span the wrap

        STDLOG(0,"Initialized MPI.\n");
        STDLOG(0,"1D node rank %d of %d total\n", MPI_rank, MPI_size);
        STDLOG(0,"2D (x,z) node rank is (%d,%d) of (%d,%d) total\n", MPI_rank_x, MPI_rank_z, MPI_size_x, MPI_size_z);
        STDLOG(2,"Original 1D rank was %d, remapped to %d\n", _world_rank, MPI_rank);

    #else
    
        int node_z_start = 0;
        int node_z_size = P.cpd;
        int node_z_start_ghost = 0;
        int node_z_size_with_ghost = P.cpd;
    
    #endif
    
    char hostname[1024];
    gethostname(hostname,1024);
    STDLOG(0,"Host machine name is %s\n", hostname);
}

void FinalizeParallel() {
    #ifdef PARALLEL

        #ifdef MULTIPLE_MPI_COMM

        STDLOG(1, "MULTIPLE_MPI_COMM was used, tearing down communicators\n");

        MPI_Comm_free(&comm_taylors);
        MPI_Comm_free(&comm_multipoles);
        MPI_Comm_free(&comm_manifest);
        MPI_Comm_free(&comm_global);

        #endif

        // Finalize MPI
        STDLOG(0,"Calling MPI_Finalize()\n");
        MPI_Finalize();
    #else
    #endif
}