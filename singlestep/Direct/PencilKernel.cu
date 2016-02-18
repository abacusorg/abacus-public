__global__ void ComputeDirects(DeviceData d, FLOAT eps2){

    __shared__ FLOAT SourceCacheX[NFBlockSize];
    __shared__ FLOAT SourceCacheY[NFBlockSize];
    __shared__ FLOAT SourceCacheZ[NFBlockSize];
    

    int id = blockDim.x*blockIdx.x + threadIdx.x;
    int myDI = 0;

    FLOAT sinkX, sinkY, sinkZ;

    sinkX = d.SinkSetPositions.X[id];
    sinkY = d.SinkSetPositions.Y[id];
    sinkZ = d.SinkSetPositions.Z[id];

    int sinkIdx = d.SinkBlockParentPencil[blockIdx.x];

    FLOAT3 a = {(FLOAT) 0.0,(FLOAT) 0.0,(FLOAT) 0.0};

    #ifdef DIRECTSPLINE
    eps2 = RSQRT(eps2);  // Direct spline uses 1/eps instead of eps^2
    #elif defined DIRECTCUBIC
    eps2 = eps2*eps2*RSQRT(eps2); // Direct cubic uses eps^3 instead of eps^2
    #endif
    

    #pragma unroll
    for(int c = d.SinkSetInteractionListStart[sinkIdx]; c < d.SinkSetInteractionListStart[sinkIdx] + WIDTH; c++){
        int sourceIdx = d.SinkSetInteractionListStart[c];
        int sourceStart = d.SourceSetStart[sourceIdx];
        int sourceCount = d.SourceSetCount[sourceIdx];
        int nB = sourceCount/NFBlockSize;

        for(int b = 0; b < nB; b+=1){
            SourceCacheX[threadIdx.x] = d.SourceSetPositions.X[sourceStart + b*NFBlockSize + threadIdx.x];
            SourceCacheY[threadIdx.x] = d.SourceSetPositions.Y[sourceStart + b*NFBlockSize + threadIdx.x];
            SourceCacheZ[threadIdx.x] = d.SourceSetPositions.Z[sourceStart + b*NFBlockSize + threadIdx.x];
            __syncthreads();
            
            myDI += NFBlockSize;
            FullDirectTile( SourceCacheX, SourceCacheY, SourceCacheZ,
                    &(a.x),&(a.y),&(a.z),
                    &sinkX, &sinkY, &sinkZ,
                    &eps2);
            __syncthreads();
        }

        int remaining = sourceCount%NFBlockSize;

        if(threadIdx.x < remaining){
            SourceCacheX[threadIdx.x] = d.SourceSetPositions.X[sourceStart + nB*NFBlockSize + threadIdx.x];
            SourceCacheY[threadIdx.x] = d.SourceSetPositions.Y[sourceStart + nB*NFBlockSize + threadIdx.x];
            SourceCacheZ[threadIdx.x] = d.SourceSetPositions.Z[sourceStart + nB*NFBlockSize + threadIdx.x];
            __syncthreads();
        }
        myDI += remaining;
        PartialDirectTile(SourceCacheX, SourceCacheY, SourceCacheZ,
                &(a.x),&(a.y),&(a.z),
                &sinkX, &sinkY, &sinkZ,
                &eps2, remaining);
        __syncthreads();
        
    }
    if(id < d.SinkSetStart[sinkIdx] + d.SinkSetCount[sinkIdx]){
        d.SinkSetAccelerations[id] = a;
        atomicAdd(&DI, myDI);
    }
}
