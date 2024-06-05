//cuda error handling wrappers

// This will output the proper CUDA error strings in the event that a CUDA host call returns an error
#define checkCudaErrors(err)           __checkCudaErrors (err, __FILE__, __LINE__)

inline void __checkCudaErrors( cudaError err, const char *file, const int line )
{
    if( cudaSuccess != err) {
        int dev=0; cudaGetDevice(&dev);
        fprintf(stderr, "%s:%i : CUDA Runtime API error %d on device %d: %s.\n",
                file, line, (int)err,dev,cudaGetErrorString( err ) );
        assert(cudaSuccess == err);//exit(-1);
    }
}

// This will output the proper error string when calling cudaGetLastError
#define checkLastCuda(msg)      __checkLastCuda (__FILE__, __LINE__)

inline void __checkLastCuda(const char *file, const int line )
{       
    cudaError_t err = cudaGetLastError(); 
    if( cudaSuccess != err) { 
        fprintf(stderr, "%s:%i : checkLastCuda() CUDA error code %d (%s).\n",
                file, line, (int)err, cudaGetErrorString( err ) );
        assert(cudaSuccess == err);
    }
}
