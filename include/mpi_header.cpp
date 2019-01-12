#include "config.h"

#ifdef PARALLEL
    #include "mpi.h"
    // This does an in-place reduction to rank 0
    #define MPI_REDUCE_TO_ZERO(vec,len,type,op) MPI_Reduce(MPI_rank!=0?(vec):MPI_IN_PLACE, vec, len, type, op, 0, MPI_COMM_WORLD)
    //#define MPI_REDUCE_TO_ZERO(vec,len,type,op) MPI_Allreduce(MPI_IN_PLACE, vec, len, type, op, MPI_COMM_WORLD)
#else
    // Just put in some stubs to help compilation
    typedef void * MPI_Request;
    #define MPI_REQUEST_NULL NULL
#endif

// For some unresolved issue, this has to be loaded before header.cpp
