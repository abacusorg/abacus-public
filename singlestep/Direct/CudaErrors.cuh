// Copyright 2012-2025 The Abacus Developers
// SPDX-License-Identifier: GPL-3.0-or-later

//cuda error handling wrappers

// This will output the proper CUDA error strings in the event that a CUDA host call returns an error
#define checkCudaErrors(err)           __checkCudaErrors (err, __FILE__, __LINE__)

inline void __checkCudaErrors( cudaError err, const fs::path &file, const int line )
{
    if( cudaSuccess != err) {
        int dev=0; cudaGetDevice(&dev);
        fmt::print(stderr, "{}:{:d} : CUDA Runtime API error {:d} on device {:d}: {:s}.\n",
                file, line, (int)err,dev,cudaGetErrorString( err ) );
        assert(cudaSuccess == err);//exit(-1);
    }
}

// This will output the proper error string when calling cudaGetLastError
#define checkLastCuda(msg)      __checkLastCuda (__FILE__, __LINE__)

inline void __checkLastCuda(const fs::path &file, const int line )
{       
    cudaError_t err = cudaGetLastError(); 
    if( cudaSuccess != err) { 
        fmt::print(stderr, "{}:{:d} : checkLastCuda() CUDA error code {:d} ({:s}).\n",
                file, line, (int)err, cudaGetErrorString( err ) );
        assert(cudaSuccess == err);
    }
}
