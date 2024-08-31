// This will output the proper CUDA error strings in the event that a CUDA host call returns an error
#define checkCudaErrors(err)           __checkCudaErrors (err, __FILE__, __LINE__)

inline void __checkCudaErrors( cuda::cudaError err, const fs::path &file, const int line )
{
    if( cuda::cudaSuccess != err) {
        int dev=0; cuda::cudaGetDevice(&dev);
        fprintf(stderr, "%s(%i) : CUDA Runtime API error %d on device %d: %s.\n",
                file, line, (int)err,dev,cuda::cudaGetErrorString( err ) );
        assert(cuda::cudaSuccess == err);//exit(-1);
    }
}

// This will output the proper error string when calling cudaGetLastError
#define getLastCudaError(msg)      __getLastCudaError (msg, __FILE__, __LINE__)

inline void __getLastCudaError( const std::string &errorMessage, const fs::path &file, const int line )
{
    cuda::cudaError_t err = cuda::cudaGetLastError();
    if( cuda::cudaSuccess != err) {
        fprintf(stderr, "%s(%i) : getLastCudaError() CUDA error : %s : (%d) %s.\n",
                file, line, errorMessage, (int)err, cuda::cudaGetErrorString( err ) );
        assert(cuda::cudaSuccess == err);
    }
}


#define checkCufftErrors(err)           __checkCudaErrors (err, __FILE__, __LINE__)

inline void __checkCudaErrors( cuda::cufftResult_t err, const fs::path &file, const int line )
{
    if( cuda::CUFFT_SUCCESS != err) {
        int dev=0; cuda::cudaGetDevice(&dev);
        fprintf(stderr, "%s(%i) : CUDA Runtime API error %d on device %d\n",
                file, line, (int)err,dev);
        assert(cuda::CUFFT_SUCCESS == err);//exit(-1);
    }
}

