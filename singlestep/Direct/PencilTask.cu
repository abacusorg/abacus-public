/// This is code for the Global Pencil-on-Pencil work.
///
/// In particular, the function GPUPencilTask receives a pointer
/// to a SetInteractionCollection and an assignment of a GPUBuffer
/// and Stream.
///
/// It must load the host-side data into Pinned Memory and then the GPU,
/// then execute the direct kernel, then copy the data back into the 
/// SIC-supplied memory.

// ======================== DeviceData =====================

/// This is the structure that is actually passed to the GPU,
/// containing all of the information from a SIC instance.

struct DeviceData{
    List3<FLOAT>    SinkSetPositions;
    List3<FLOAT>    SourceSetPositions;
    
    int *           SinkSetIdMax;
    accstruct *        SinkSetAccelerations;
    FLOAT	    b2;
    int             NSinkBlocks;
    int *           SinkBlockParentPencil;

    int *           SourceSetStart;
    int *           SourceSetCount;
    int             NSourceBlocks;

    int             NSinkSets;
    int             NSourceSets;
    int             InteractionCount;

    int *           SinkSourceInteractionList;
    FLOAT *         SinkSourceYOffset;
};

// Here's the GPU code that knows what to do with a DeviceData structure
#include "PencilKernel.cu"

// ================ Configure Buffer as DeviceData ==============

// Just a little thing to avoid compiler warnings
void SetPointer(char **p, char *val) { *p = val; }

// Some macros to help readability
#define   CudaConfig(ptr, size) SetPointer((char **)&ptr, buf.device+used_gpu); used_gpu+=(size/4096+1)*4096;

#define     WCConfig(ptr, size) SetPointer((char **)&ptr, buf.hostWC+used_hostWC); used_hostWC+=(size/4096+1)*4096;

#define PinnedConfig(ptr, size) SetPointer((char **)&ptr, buf.host+used_host); used_host+=(size/4096+1)*4096;


/// Given a Buffer, we want to set pointers in the DeviceData's that
/// divide it up to our needed vectors.
void ConfigureBufferAsDeviceData(GPUBuffer &buf, 
	DeviceData &gpu, DeviceData &pinned) {
    uint64 used_gpu = 0;
    uint64 used_host = 0;
    uint64 used_hostWC = 0;

    // Allocate GPU-side memory
    CudaConfig(gpu.SinkSetIdMax,              sizeof(int) * MaxNSink);
    CudaConfig(gpu.SourceSetStart,            sizeof(int) * MaxNSource);
    CudaConfig(gpu.SourceSetCount,            sizeof(int) * MaxNSource);
    CudaConfig(gpu.SinkSourceInteractionList, sizeof(int) * MaxNSink * WIDTH);
    CudaConfig(gpu.SinkSourceYOffset,         sizeof(FLOAT) * MaxNSink*WIDTH);
    CudaConfig(gpu.SinkBlockParentPencil,     sizeof(int) * MaxSinkBlocks);
    CudaConfig(gpu.SinkSetPositions.X,        sizeof(FLOAT) * MaxSinkSize);
    CudaConfig(gpu.SinkSetPositions.Y,        sizeof(FLOAT) * MaxSinkSize);
    CudaConfig(gpu.SinkSetPositions.Z,        sizeof(FLOAT) * MaxSinkSize);
    CudaConfig(gpu.SinkSetAccelerations,      sizeof(accstruct) * MaxSinkSize);
    CudaConfig(gpu.SourceSetPositions.X,      sizeof(FLOAT) * MaxSourceSize);
    CudaConfig(gpu.SourceSetPositions.Y,      sizeof(FLOAT) * MaxSourceSize);
    CudaConfig(gpu.SourceSetPositions.Z,      sizeof(FLOAT) * MaxSourceSize);
    assertf(used_gpu<buf.size, "Configuration of Buffer requesting %ld bytes on device, but only %ld available\n", used_gpu, buf.size);   // Check that we didn't overflow

    // Allocate host-side buffers
    WCConfig(pinned.SinkSetIdMax,              sizeof(int) * MaxNSink);
    WCConfig(pinned.SourceSetStart,            sizeof(int) * MaxNSource);
    WCConfig(pinned.SourceSetCount,            sizeof(int) * MaxNSource);
    WCConfig(pinned.SinkSourceInteractionList, sizeof(int) * MaxNSink * WIDTH);
    WCConfig(pinned.SinkSourceYOffset,         sizeof(FLOAT) * MaxNSink * WIDTH);
    WCConfig(pinned.SinkBlockParentPencil,     sizeof(int) * MaxSinkBlocks);
    WCConfig(pinned.SinkSetPositions.X,        sizeof(FLOAT) * MaxSinkSize);
    WCConfig(pinned.SinkSetPositions.Y,        sizeof(FLOAT) * MaxSinkSize);
    WCConfig(pinned.SinkSetPositions.Z,        sizeof(FLOAT) * MaxSinkSize);
    PinnedConfig(pinned.SinkSetAccelerations,  sizeof(accstruct) * MaxSinkSize);
    WCConfig(pinned.SourceSetPositions.X,      sizeof(FLOAT) * MaxSourceSize);
    WCConfig(pinned.SourceSetPositions.Y,      sizeof(FLOAT) * MaxSourceSize);
    WCConfig(pinned.SourceSetPositions.Z,      sizeof(FLOAT) * MaxSourceSize);
    assert(used_host<buf.sizeDef);
    assert(used_hostWC<buf.sizeWC);
    assertf(used_host<buf.sizeDef, "Configuration of Buffer requesting %ld bytes on host Def, but only %ld available\n", used_host, buf.sizeDef);   // Check that we didn't overflow
    assertf(used_hostWC<buf.sizeWC, "Configuration of Buffer requesting %ld bytes on host WC, but only %ld available\n", used_host, buf.sizeWC);   // Check that we didn't overflow
    return;
}

#undef CudaConfig
#undef WCConfig
#undef PinnedConfig

// ============= GPUPencilTask: executing one SetInteractionCollection =======

// Space-saving macros.  We rely on the idea that the name of a field
// is the same in the SetInteractionCollection and the DeviceData structures.

// This copies from the SIC task to the Pinned host-side DeviceData,
// and then to the GPU.
#define CopyToGPU(name,size)\
    thissize = size;\
    memcpy(PinnedBuffer.name, task->name, thissize);\
    checkCudaErrors(cudaMemcpyAsync(StreamData.name, PinnedBuffer.name, thissize, cudaMemcpyHostToDevice, DeviceStreams[g]));\
    size_to_gpu += thissize

// For the particle lists, we've loaded the pinned memory separately,
// so we only need to copy to the GPU
#define CopyListToGPU(name,d,size)\
    thissize = size;\
    checkCudaErrors(cudaMemcpyAsync(StreamData.name.d, PinnedBuffer.name.d, thissize, cudaMemcpyHostToDevice, DeviceStreams[g]));\
    size_to_gpu += thissize

// This copies back from the GPU to Pinned memory
#define CopyFromGPU(name,size)\
    thissize = size;\
    checkCudaErrors(cudaMemcpyAsync(PinnedBuffer.name, StreamData.name, thissize, cudaMemcpyDeviceToHost, DeviceStreams[g]));\
    size_from_gpu += thissize



/// This routine is invoked on a single SetInteractionCollection.
///
/// It must copy the relevant contents to a DeviceData instance in
/// the pinned memory, then copy that to the matching DeviceData 
/// instance on the GPU.  Then it invokes the GPU kernel, with
/// the GPU-side DeviceData as the argument.  Then it copies the
/// partial accelerations back from the GPU to pinned memory, and
/// then to space that was allocated in the SIC.  Then marks the SIC
/// as completed.
///
/// Importantly, this task has substantial host-side work too, because
/// the SIC does not contain the actual particle data, but only the 
/// instructions on where to find it.  So we must invoke the SinkPencilPlan's
/// and SourcePencilPlan's to actually copy the data from the slab locations
/// into the Pinned memory.

void GPUPencilTask(void *item, int g){
#ifdef CUDADIRECT
    // Given a buffer/stream number g, and a pointer to our task, do it.
    SetInteractionCollection *task = (SetInteractionCollection *)item;
    task->DeviceThreadTimer.Start();
    task->LaunchDeviceKernels.Start();
    task->AssignedDevice = g;
    checkCudaErrors(cudaSetDevice(g%NGPU));  // only needed for blocking directs

    DeviceData StreamData, PinnedBuffer;

    ConfigureBufferAsDeviceData(Buffers[g], StreamData, PinnedBuffer);

    //Schedule data copy to GPU
    size_t thissize = 0;
    uint64 size_to_gpu = 0, size_from_gpu = 0;
    StartThroughputTimer(DeviceStreams[g], cudaSuccess, (void *) task);
    // Zero the device-side accelerations
    checkCudaErrors(cudaMemsetAsync(StreamData.SinkSetAccelerations,
                0 , sizeof(accstruct) * NFBlockSize * task->NSinkBlocks, DeviceStreams[g]));
    StreamData.InteractionCount = task->InteractionCount;
    StreamData.b2 = task->b2;

    // Need to load the particles to the PinnedBuffer.
    // Copy the sinks into position
    task->LaunchDeviceKernels.Stop();
    task->FillSinks.Start();
    for (int j=0; j<task->NSinkSets; j++) {
        task->SinkPlan[j].copy_into_pinned_memory(PinnedBuffer.SinkSetPositions, task->SinkSetStart[j], task->SinkSetCount[j]);
    }
    task->FillSinks.Stop();
    task->LaunchDeviceKernels.Start();
    
    // Now copy these to the GPU
    CopyListToGPU(SinkSetPositions, X, sizeof(FLOAT) * task->NSinkBlocks * NFBlockSize);
    CopyListToGPU(SinkSetPositions, Y, sizeof(FLOAT) * task->NSinkBlocks * NFBlockSize);
    CopyListToGPU(SinkSetPositions, Z, sizeof(FLOAT) * task->NSinkBlocks * NFBlockSize);

    task->LaunchDeviceKernels.Stop();

    // Repeat this with the sources
    task->FillSources.Start();
    for (int j=0; j<task->NSourceSets; j++) {
        task->SourcePlan[j].copy_into_pinned_memory(PinnedBuffer.SourceSetPositions, task->SourceSetStart[j], task->SourceSetCount[j]);
    }
    task->FillSources.Stop();
    task->LaunchDeviceKernels.Start();
    
    CopyListToGPU(SourceSetPositions, X, sizeof(FLOAT) * task->NSourceBlocks * NFBlockSize);
    CopyListToGPU(SourceSetPositions, Y, sizeof(FLOAT) * task->NSourceBlocks * NFBlockSize);
    CopyListToGPU(SourceSetPositions, Z, sizeof(FLOAT) * task->NSourceBlocks * NFBlockSize);

    // Now copy other information from the SIC to the GPU
    CopyToGPU(SinkSetIdMax,                 sizeof(int)*task->NSinkSets);
    CopyToGPU(SinkBlockParentPencil,        sizeof(int)*task->NSinkBlocks);
    CopyToGPU(SourceSetStart,               sizeof(int)*task->NSourceSets);
    CopyToGPU(SourceSetCount,               sizeof(int)*task->NSourceSets);
    CopyToGPU(SinkSourceInteractionList,    sizeof(int)*task->InteractionCount);
    CopyToGPU(SinkSourceYOffset,            sizeof(FLOAT)*task->InteractionCount);

    // Run the GPU kernel for this Interaction collection!
    dim3 dimGrid(task->NSinkBlocks);
    dim3 dimBlock(NFBlockSize);
    ComputeDirects<<<dimGrid,dimBlock,0,DeviceStreams[g]>>>(StreamData,task->eps);
    // Control won't return until it's done
    
    // Copy back results from GPU
    // If the memory is unpinned, this is blocking
    CopyFromGPU(SinkSetAccelerations, sizeof(accstruct) * NFBlockSize * task->NSinkBlocks);
    
    task->bytes_to_device = size_to_gpu;
    task->bytes_from_device = size_from_gpu;
    task->LaunchDeviceKernels.Stop();

    // It's safe for the next work unit to wipe the host PinnedBuffer once everything's been transferred
    // We could allow sink/source copying while waiting for accels, but we'd need a signalling mechanism from the callback
    task->WaitForResult.Start();
    checkCudaErrors(cudaStreamSynchronize(DeviceStreams[g]));
    task->WaitForResult.Stop();

    // Now copy the data from Pinned back to the SIC buffer
    task->CopyAccelFromPinned.Start();
    // SetInteractionCollection * task = (SetInteractionCollection *) data;
    memcpy(task->SinkSetAccelerations, PinnedBuffer.SinkSetAccelerations, sizeof(accstruct) * NFBlockSize * task->NSinkBlocks);

    // Declare victory!
    task->SetCompleted();
    MarkCompleted(DeviceStreams[g], cudaSuccess, (void *) task);
    task->CopyAccelFromPinned.Stop();
    task->DeviceThreadTimer.Stop();
#endif
}

#undef CopyToGPU
#undef CopyListToGPU
#undef CopyFromGPU
