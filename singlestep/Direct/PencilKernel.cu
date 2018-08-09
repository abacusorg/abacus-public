/// \file These are the GPU-side routines that implement the Block-on-Block computation.

/// We specialize these Block-based kernels to the choice of 
/// COMPUTE_FOF_DENSITY, since we don't expect to re-use them in
/// the microstepping.
///
/// If the Source Block is full, then we can unroll some of the loop
#ifndef PTXDIRECT
__device__ inline void FullDirectTile(
        FLOAT * source_cache_x, FLOAT * source_cache_y, FLOAT * source_cache_z,
        FLOAT   *sinkx,     FLOAT * sinky,    FLOAT * sinkz,
        FLOAT *ax, FLOAT *ay, FLOAT *az, 
	#ifdef COMPUTE_FOF_DENSITY
	FLOAT *aw,
	#endif
        FLOAT *eps, FLOAT *b2){

    #pragma unroll 16
    for(int n = 0; n <NFBlockSize;n++){
        FLOAT sourcex = source_cache_x[n];
        FLOAT sourcey = source_cache_y[n];
        FLOAT sourcez = source_cache_z[n];
	#ifdef COMPUTE_FOF_DENSITY
        direct_density<0>(*sinkx,*sinky,*sinkz,
                sourcex,sourcey,sourcez,
                *ax,*ay,*az,*aw,*eps,*b2);
	#else
        direct<0>(*sinkx,*sinky,*sinkz,
                sourcex,sourcey,sourcez,
                *ax,*ay,*az,*eps,*b2);
	#endif
    }
}
#else
#include "ASMDirectTile.cu"
#endif

/// As for FullDirectTile(), but if one doesn't have all of the sources, 
/// then one has to actually check the loop bound.
__device__ inline void PartialDirectTile(
        FLOAT * source_cache_x, FLOAT * source_cache_y, FLOAT * source_cache_z,
        FLOAT  * sinkx,     FLOAT *sinky,    FLOAT *sinkz,
        FLOAT *ax, FLOAT *ay, FLOAT *az, 
	#ifdef COMPUTE_FOF_DENSITY
	FLOAT *aw,
	#endif
        FLOAT *eps, FLOAT *b2, int &i){

    for(int n = 0; n <i;n++){
        FLOAT sourcex = source_cache_x[n];
        FLOAT sourcey = source_cache_y[n];
        FLOAT sourcez = source_cache_z[n];
	#ifdef COMPUTE_FOF_DENSITY
        direct_density<0>(*sinkx,*sinky,*sinkz,
                sourcex,sourcey,sourcez,
                *ax,*ay,*az,*aw,*eps,*b2);
	#else
        direct<0>(*sinkx,*sinky,*sinkz,
                sourcex,sourcey,sourcez,
                *ax,*ay,*az,*eps,*b2);
	#endif
    }
}



// ======================  ComputeDirects ======================

/// The routine that accepts a DeviceData instance and unpacks
/// for its thread number what action should be taken, then invokes that.
///
/// eps may be 1/eps, eps^2, eps^3, etc, depending on the type of softening
__global__ void ComputeDirects(DeviceData d, FLOAT eps){

    __shared__ FLOAT SourceCacheX[NFBlockSize];
    __shared__ FLOAT SourceCacheY[NFBlockSize];
    __shared__ FLOAT SourceCacheZ[NFBlockSize];
    

    int id = blockDim.x*blockIdx.x + threadIdx.x;
    int myDI = 0;

    FLOAT sinkX, sinkY, sinkZ;
    int sinkIdx = d.SinkBlockParentPencil[blockIdx.x];
    if(id < d.SinkSetIdMax[sinkIdx]){
        sinkX = d.SinkSetPositions.X[id];
        sinkY = d.SinkSetPositions.Y[id];
        sinkZ = d.SinkSetPositions.Z[id];
    }else{
         sinkX =0;
         sinkY =0;
         sinkZ =0;
    }


#ifdef COMPUTE_FOF_DENSITY
    accstruct a = {(FLOAT)0.0,(FLOAT)0.0,(FLOAT)0.0,(FLOAT)0.0};
#else
    accstruct a = {(FLOAT)0.0,(FLOAT)0.0,(FLOAT)0.0};
#endif
    
    int InteractionStart = sinkIdx * WIDTH;
    int InteractionMax =  InteractionStart + WIDTH;

    #pragma unroll
    for(int c = InteractionStart; c < InteractionMax; c++){
        int sourceIdx = d.SinkSourceInteractionList[c];
        FLOAT yOffset = d.SinkSourceYOffset[c];
        int sourceStart = d.SourceSetStart[sourceIdx];
        int sourceCount = d.SourceSetCount[sourceIdx];
        int nB = sourceCount/NFBlockSize;

        for(int b = 0; b < nB; b+=1){
            int idx = sourceStart + b*NFBlockSize + threadIdx.x;
            SourceCacheX[threadIdx.x] = d.SourceSetPositions.X[idx];
            SourceCacheY[threadIdx.x] = d.SourceSetPositions.Y[idx]+yOffset;
            SourceCacheZ[threadIdx.x] = d.SourceSetPositions.Z[idx];
            __syncthreads();
            
            myDI += NFBlockSize;
            FullDirectTile( SourceCacheX, SourceCacheY, SourceCacheZ,
                    &sinkX, &sinkY, &sinkZ,
                    &(a.x),&(a.y),&(a.z),
		    #ifdef COMPUTE_FOF_DENSITY
		    &(a.w),
		    #endif
                    &eps,&d.b2);  // try non-pointer?
            __syncthreads();

        }

        int remaining = sourceCount%NFBlockSize;

        if(threadIdx.x < remaining){
            int idx = sourceStart + nB*NFBlockSize + threadIdx.x;
            SourceCacheX[threadIdx.x] = d.SourceSetPositions.X[idx];
            SourceCacheY[threadIdx.x] = d.SourceSetPositions.Y[idx]+yOffset;
            SourceCacheZ[threadIdx.x] = d.SourceSetPositions.Z[idx];
        }
        __syncthreads();
        
        myDI += remaining;
        PartialDirectTile(SourceCacheX, SourceCacheY, SourceCacheZ,
                &sinkX, &sinkY, &sinkZ,
                &(a.x),&(a.y),&(a.z),
		#ifdef COMPUTE_FOF_DENSITY
		&(a.w),
		#endif
                &eps, &d.b2, remaining);
        __syncthreads();
    }

    if(id < d.SinkSetIdMax[sinkIdx]){
        // TODO: haven't seen these trip in a long, long time. Remove them?
        /*assert(::isfinite(a.x));
        assert(::isfinite(a.y));
        assert(::isfinite(a.z));
        assert(::isfinite(a.w));*/
        d.SinkSetAccelerations[id] = a;
        atomicAdd(&DI, myDI);
    }

}
