#include "config.h"

// These are the logical (i.e. global) indices that track the z-range belonging to this node
// Initialized in parallel.cpp
int node_z_start;
int node_z_size;
int node_z_start_ghost;
int node_z_size_with_ghost;

// The above names are pretty verbose, which is unfortunate. And we often end up writing
// "node_z_start_ghost + node_z_size_with_ghost" to mark the end of iteration, which is long.
// But a new var like "node_z_end_ghost" seems ambiguous to LHG: is it the last cell, or one past that?
// LHG is open to new names, but for now we opt for descriptive ones.

int MPI_size = 1, MPI_rank = 0;     // We'll set these globally, so that we don't have to keep fetching them

// The 2D rank and size
int MPI_rank_x = 0, MPI_size_x = 1;
int MPI_rank_z = 0, MPI_size_z = 1;

int GHOST_RADIUS = 0;        // in the read state
int MERGE_GHOST_RADIUS = 0;  // in the write state

char NodeString[8] = "";     // Set to "" for serial, ".NNNN" for MPI

int _world_rank;  // purely informational, do not use

int *all_node_z_start;  // z domains of all nodes
int *all_node_z_size;

int node_cpdp1half;  // (CPD+1)/2/MPI_size_z
int *all_node_ky_start;  // 2D: start of the ky multipole splits
int *all_node_ky_size;  // 2D: sizes of the multipole splits

#ifdef PARALLEL
    #include "mpi.h"

	// These are the communicators we will use instead of MPI_COMM_WORLD
	// They will be aliased to MPI_COMM_WORLD if ./configure --disable-multiple-comm
	MPI_Comm comm_taylors;
	MPI_Comm comm_multipoles;
	MPI_Comm comm_manifest;
	MPI_Comm comm_global; // e.g. for state reductions
    
    MPI_Comm comm_2d;  // a communicator for the nodes splitting a slab
    MPI_Comm comm_1d_x;  // the nodes for all x, of the same z split
    MPI_Comm comm_1d_z;  // the nodes for all z, of the same x
    
    MPI_Comm comm_multipoles_z;  // 2D
    MPI_Comm comm_taylors_z;  // 2D

    MPI_Datatype MPI_ilstruct;  // for the 2D code

    // This does an in-place reduction to rank 0
    #define MPI_REDUCE_TO_ZERO(vec,len,type,op) MPI_Reduce(MPI_rank!=0?(vec):MPI_IN_PLACE, vec, len, type, op, 0, comm_global)
    // #define MPI_REDUCE_TO_ZERO(vec,len,type,op) while(0)
    //#define MPI_REDUCE_TO_ZERO(vec,len,type,op) MPI_Allreduce(MPI_IN_PLACE, vec, len, type, op, MPI_COMM_WORLD)
#else
    // Just put in some stubs to help compilation
    typedef void * MPI_Request;
    #define MPI_REQUEST_NULL NULL
#endif

void SetupNeighborExchange(int first, int nslab, int make_noop);
void DoNeighborSend(int slab);
int AttemptNeighborReceive(int first, int receive_ahead);
int IsNeighborReceiveDone(int slab);
void TeardownNeighborExchange();
