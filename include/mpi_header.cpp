#include "config.h"

#ifdef PARALLEL
#include "mpi.h"
#else
    // Just put in some stubs to help compilation
    typedef void * MPI_Request;
    #define MPI_REQUEST_NULL NULL
#endif

// For some unresolved issue, this has to be loaded before header.cpp
