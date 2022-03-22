#include "config.h"

#ifdef PARALLEL
    #include "mpi.h"

	// These are the communicators we will use instead of MPI_COMM_WORLD
	// They will be aliased to MPI_COMM_WORLD if ./configure --disable-multiple-comm
	MPI_Comm comm_taylors;
	MPI_Comm comm_multipoles;
	MPI_Comm comm_manifest;
	MPI_Comm comm_global; // e.g. for state reductions
    
    MPI_Comm comm_2d;  // a communicator for the nodes splitting a slab

    // This does an in-place reduction to rank 0
    #define MPI_REDUCE_TO_ZERO(vec,len,type,op) MPI_Reduce(MPI_rank!=0?(vec):MPI_IN_PLACE, vec, len, type, op, 0, comm_global)
    // #define MPI_REDUCE_TO_ZERO(vec,len,type,op) while(0)
    //#define MPI_REDUCE_TO_ZERO(vec,len,type,op) MPI_Allreduce(MPI_IN_PLACE, vec, len, type, op, MPI_COMM_WORLD)
#else
    // Just put in some stubs to help compilation
    typedef void * MPI_Request;
    #define MPI_REQUEST_NULL NULL
#endif

// For some unresolved issue, this has to be loaded before header.cpp
