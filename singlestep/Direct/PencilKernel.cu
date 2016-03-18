__global__ void ComputeDirects(DeviceData d, FLOAT eps2){

    __shared__ FLOAT SourceCacheX[NFBlockSize];
    __shared__ FLOAT SourceCacheY[NFBlockSize];
    __shared__ FLOAT SourceCacheZ[NFBlockSize];
    

    int id = blockDim.x*blockIdx.x + threadIdx.x;
    int myDI = 0;

    FLOAT sinkX, sinkY, sinkZ;
    int sinkIdx = d.SinkBlockParentPencil[blockIdx.x];
    if(id < d.SinkSetStart[sinkIdx] + d.SinkSetCount[sinkIdx]){
        sinkX = d.SinkSetPositions.X[id];
        sinkY = d.SinkSetPositions.Y[id];
        sinkZ = d.SinkSetPositions.Z[id];
    }else{
         sinkX =0;
         sinkY =0;
         sinkZ =0;
    }


    FLOAT3 a = {(FLOAT) 0.0,(FLOAT) 0.0,(FLOAT) 0.0};

    #ifdef DIRECTSPLINE
    eps2 = RSQRT(eps2);  // Direct spline uses 1/eps instead of eps^2
    #elif defined DIRECTCUBIC
    eps2 = eps2*eps2*RSQRT(eps2); // Direct cubic uses eps^3 instead of eps^2
    #endif
    
    int InteractionStart = sinkIdx * WIDTH;
    int InteractionMax =  InteractionStart + WIDTH;

    #pragma unroll
    for(int c = InteractionStart; c < InteractionMax; c++){
        int sourceIdx = d.SinkSourceInteractionList[c];
        int sourceStart = d.SourceSetStart[sourceIdx];
        int sourceCount = d.SourceSetCount[sourceIdx];
        int nB = sourceCount/NFBlockSize;

        for(int b = 0; b < nB; b+=1){
            int idx = sourceStart + b*NFBlockSize + threadIdx.x;
            SourceCacheX[threadIdx.x] = d.SourceSetPositions.X[idx];
            SourceCacheY[threadIdx.x] = d.SourceSetPositions.Y[idx];
            SourceCacheZ[threadIdx.x] = d.SourceSetPositions.Z[idx];
            __syncthreads();
            
            myDI += NFBlockSize;
            FullDirectTile( SourceCacheX, SourceCacheY, SourceCacheZ,
                    &sinkX, &sinkY, &sinkZ,
                    &(a.x),&(a.y),&(a.z),
                    &eps2);
            __syncthreads();

        }

        int remaining = sourceCount%NFBlockSize;

        if(threadIdx.x < remaining){
            int idx = sourceStart + nB*NFBlockSize + threadIdx.x;
            SourceCacheX[threadIdx.x] = d.SourceSetPositions.X[idx];
            SourceCacheY[threadIdx.x] = d.SourceSetPositions.Y[idx];
            SourceCacheZ[threadIdx.x] = d.SourceSetPositions.Z[idx];
        }
        __syncthreads();
        
        myDI += remaining;
        PartialDirectTile(SourceCacheX, SourceCacheY, SourceCacheZ,
                &sinkX, &sinkY, &sinkZ,
                &(a.x),&(a.y),&(a.z),
                &eps2, remaining);
        __syncthreads();
    }

    if(id < d.SinkSetStart[sinkIdx] + d.SinkSetCount[sinkIdx]){
        assert(isfinite(a.x));
        assert(isfinite(a.y));
        assert(isfinite(a.z));
        //atomicAdd(&(d.SinkSetAccelerations[id].x),a.x);
        //atomicAdd(&(d.SinkSetAccelerations[id].y),a.y);
        //atomicAdd(&(d.SinkSetAccelerations[id].z),a.z);
        d.SinkSetAccelerations[id] = a;
        atomicAdd(&DI, myDI);
    }

}
