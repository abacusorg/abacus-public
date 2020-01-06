/** \file The top-level CUDA code to control the GPUs, particularly in 
the non-blocking mode.

In particular, the preferred non-blocking mode operates several 
independent threads, typically 3 per CUDA Device.  A task operates
on one thread end-to-end, including communication to & from the GPU
and the computation itself.  Each thread controls one CUDA stream.
Because there are multiple streams on each GPU, the GPU can be 
busy computing on one task while communication happens on others.

For the global near-field forces, the work is organized into Sink
Pencils and Source Pencils.  A given cell is part of multiple sink
pencils and hence generates multiple partial accelerations that must
be coadded later.  A task would involve a portion of the sink pencils
as well as the associated source pencils, subdivided so that the 
memory on the GPU is not exhausted. 

For the microstepping, the task involves executing multiple kick-drift-kick
tree-based steps on the GPU, returning the positions and velocities.
The task includes many groups.

Each task therefore has to be prepared on the CPU, assembling a large 
amount of data and the inter-relations needed for the computation.
The GPU task then copies this data to a pinned memory buffer on the
host side and then onto the GPU.  Then the computation kernel is called.
When that is done, the result is copied back to the pinned memory
and then to the destination.

This file contains the functions to set up the GPU threads, 
allocate the pinned memory, and then to run the GPU thread by
receiving tasks and dispensing them to the correct function.

The data model and specifics of the two task use-cases are in 
other files.

Pinned memory is slow to allocate and so we set up a large buffer
for each thread at the beginning of the time slice.  Each task
then partitions up its thread's buffer to match its needs.

This file is compiled separately (hence the headers) and then linked
to the rest of the code, because it needs CUDA.

*/

#include "config.h"
#include "header.cpp"


#ifdef __INTEL_COMPILER
#include <x86intrin.h>  // needed for _rdtsc(). only for nvcc using icc with certain gcc versions?
#endif

#include <cstring>
#include <cstdio>
#include <cassert>
#include <pthread.h>
#include "omp.h"

#ifdef HAVE_LIBNUMA
#include <numaif.h>
#endif

#include "CudaErrors.cuh"

#ifndef NFRADIUS
#define NFRADIUS 2
#endif

#define WIDTH (2*NFRADIUS +1)

// TODO: these are NOT using the Abacus-defined float3/double3
// CUDA unfortunately defines the same names in the global namespace; not sure if we can containerize that
// For now, the types remain binary-compatible though
#ifdef DOUBLEPRECISION
    #define FLOAT double
    #define FLOAT4 double4
    #define FLOAT3 double3
    #define RSQRT rsqrt
    #define MIN fmin
    #define ItemSize 2
#else
    #define FLOAT float
    #define FLOAT4 float4
    #define FLOAT3 float3
    #define RSQRT rsqrtf
    #define MIN fminf
    #define ItemSize 1
#endif

#define posstruct FLOAT3
#ifdef COMPUTE_FOF_DENSITY
    #define accstruct FLOAT4
#else
    #define accstruct FLOAT3
#endif

#include "STimer.h"
#include "DeviceFunctions.h"
#include "IncludeGPUKernels.cuh"

// Provide a mechanism to call singlestep's STDLOG from this compilation unit
// TODO: consider writing to a new log file, as we do with the IO
void stdlog_hook(int verbosity, const char* str);
#define STDLOG(_verbosity,...) do { \
    char logstr[1024]; \
    sprintf(logstr, __VA_ARGS__);\
    stdlog_hook(_verbosity, logstr);\
} while(0)

#define assertf(_mytest,...) do { \
    if (!(_mytest)) { \
        STDLOG(0,"Failed Assertion: %s\n", #_mytest);\
        fprintf(stderr,"Failed Assertion: %s\n", #_mytest);\
        STDLOG(0,__VA_ARGS__);\
        assert(0==99); \
    }} while(0)


// Forward declarations that will get linked later
int set_core_affinity(int core_id);

#include "SetInteractionCollection.hh"

int NGPU = 0;                // Number of actual CUDA devices (not threads)
int NGPUQueues = 0;  // Number of actual work queues

// ===================  The Queue of Tasks =====================

/// We communicate to the queues with a simple FIFO queue. 
/// Tasks go in, and are read out by the QueueWatchers.
/// This class is the item that goes into that queue for
/// the inter-thread communication.
struct GPUQueueTask {
    void *task; ///< The pointer to the task data
    int type;   ///< == TASKTYPE_KILL signals to kill the Watcher
                ///< == TASKTYPE_SIC is for SetInteractionCollections
};

enum GPUQueueTaskType { TASKTYPE_KILL,
                        TASKTYPE_SIC, };

#include "tbb/concurrent_queue.h"
tbb::concurrent_bounded_queue<GPUQueueTask> *work_queues;



// =============== GPU Buffers =============================

/// The memory in each GPU is divided into a handful of equal sized
/// buffers, with matching pinned space on the host.  
/// Each buffer is paired with one thread.  Tasks will occupy one
/// buffer.  The memory allocations persist until the thread ends.

struct GPUBuffer {
    uint64 size;        ///< The size in bytes
    uint64 sizeWC;        ///< The size in bytes of the WC part
    uint64 sizeDef;        ///< The size in bytes of the non-WC part
    char * device;        ///< The device-side memory
    char * host;        ///< The host-side memory allocated as default
    char * hostWC;        ///< The host-side memory allocated as write combined
    volatile int ready;            ///< Have we finished initializing the CUDA stream for this buffer?
};

// GPUSetup() spins up one thread per GPU buffer
// Each will be passed this struct which contains information
// about which GPU to attach to, whether to use pinned memory, etc.
struct ThreadInfo {
    int thread_num;  // thread/buffer number
    int core;
    int queue;  // which queue are we watching?
    int UsePinnedGPUMemory;
    pthread_barrier_t *barrier; // this is a pointer because threads share barriers
    pthread_t thread;  // the thread itself
};

cudaStream_t * DeviceStreams;        ///< The actual CUDA streams, one per thread
GPUBuffer * Buffers;                 ///< The pointers to the allocated memory
ThreadInfo *DeviceThreads;        ///< The threads!


// =============== Code to invoke execution on a SIC ============

// int CurrentGPU = 0;

void GPUPencilTask(void *, int);

/// This submits a single SetInteractionCollection to the GPU queue.
void SetInteractionCollection::GPUExecute(int blocking){
    // This method is invoked on a SIC in order to schedule it for execution.
    // At this point, we are directing work to individual GPUs.

    // Push to the NUMA-appropriate queue
    int q = ((j_low + j_high) / 2. / cpd) * NGPUQueues;
    if(j_low == j_high && j_high == cpd)
        q--;
    assert(q < NGPUQueues);

    Blocking = blocking;
    
    if(blocking){
        // We're running the task from the main thread;
        // make sure the GPU threads have finished the necessary CUDA setup
        while(!Buffers[0].ready){}
        
        // Send to buffer 0. No need to load balance when blocking!
        GPUPencilTask(this, 0);
    }
    else {
        GPUQueueTask t;
        t.type = TASKTYPE_SIC;
        t.task = (void *)this;
        work_queues[q].push(t);
    }
}


// ====================== GPU configuration =====================

/// Function to get the number of GPUs
extern "C"
// Use ./configure --with-max-gpus=1 to force 1 GPU
int GetNGPU(){
    int ngpu;
    checkCudaErrors(cudaGetDeviceCount(&ngpu));
    if(ngpu > MAX_GPUS)
        ngpu = MAX_GPUS;
    return ngpu;
}

/// Function to get the amount of memory on the GPUs
extern "C" double GetDeviceMemory(){
    int ngpu;
    checkCudaErrors(cudaGetDeviceCount(&ngpu));
    double mem = 1e99;
    for(int g = 0; g < ngpu; g++){
        cudaDeviceProp p;
        checkCudaErrors(cudaGetDeviceProperties(&p,g));
        double m = p.totalGlobalMem/1e9;
        mem = min(.8*m,mem);
    }
    return mem;
}


// =============== Setting up the GPUs =========================

// These are the sizes we allocate on the GPU/pinned memory
int MaxSinkBlocks, MaxSourceBlocks;
int MaxNSink, MaxNSource;
size_t MaxSinkSize, MaxSourceSize;

int init = 0;
int BPD;                 // Buffers per device

static volatile uint64 host_alloc_bytes;

void *QueueWatcher(void *);

// Here is the routine that is called to configure the GPU threads
// and then actually initiate the threads.

// We given the number of GPUs, buffers per GPU, and bytes per buffer.
// Also CPD and the number of particles.
// We use this to compute and return maximum task sizes.

/*
For the Pencil task, the average memory usage:

Memory per Block:
Sink pos and accel:  
        FLOAT3*NFBlockSize bytes of WC memory
        + accstruct*NFBlockSize bytes of Default memory
Source Pos: FLOAT3*NFBLockSize, but we expect 5x more SourceBlocks
SinkBlockParentPencil: An extra 4 bytes per block.

Sink Pencils: 4+(4+FLOAT)*WIDTH
Source Pencils: 4+4
    This is very small, typically <1%
    We can be conservative and assume this to be per block rather than per cell.
    This slightly neglects the fact that the sources need 4 boundary
    rows, but we just built in a little extra.

In summary:
Per SourceBlock: FLOAT3*NFBlockSize+8/WIDTH of WC memory
Per SinkBlock:   FLOAT3*NFBlockSize+4+FLOAT+4/WIDTH of WC memory
                 accstruct*NFBlockSize of default memory

Assuming that the number of SourceBlocks is WIDTHx SinkBlocks,
this sets the ratio of WC to default memory.  And we can then
return the maximum number of Source and Sink Blocks, so that
the upstream routine can size the SIC to match.

We ignore the fact that the Sink/Source Pencil count may be 
limited by the fraction of CPD**2

*/



/// This launches all of the GPU threads and associated control.
///
/// We will be running multiple GPU threads, usually a few per device.
/// Each thread can only do one task (e.g., a SIC) at a time, and 
/// it controls a single CUDA stream that guarantees sequential 
/// execution of the load, compute, and store cycle.
/// We have multiple threads per GPU so that load/store can be
/// overlapped with compute.
/// 
/// This routine is supplied with the size (in bytes) of each
/// GPU buffer, and it returns estimates on that basis of how
/// big the tasks can be.  This is used in the planning of tasks.
///
/// Note that we size this to NFRADIUS, which is allowed to be bigger
/// than P.NearFieldRadius.

extern "C" void GPUSetup(int cpd, uint64 MaxBufferSize, 
    int numberGPUs, int bufferperdevice, 
    int *ThreadCoreStart, int NThreadCores,
    int *GPUQueueAssignments,
    int *maxsinkblocks, int *maxsourceblocks,
    int UsePinnedGPUMemory) {

    if (init != 0) return;
    BPD = bufferperdevice;
    NGPU = numberGPUs;
    int NBuf = BPD*NGPU;

    // Determine the number of queues by the max assigned queue
    for (int i = 0; i < NGPU; i++){
        NGPUQueues = max(NGPUQueues, GPUQueueAssignments[i]+1);
    }
    STDLOG(1, "Using %d GPU work queues\n", NGPUQueues);

    // The following estimates are documented above, and assume sizeof(int)=4
    float BytesPerSourceBlockWC = sizeof(FLOAT)*3*NFBlockSize+8.0;
    float BytesPerSinkBlockWC = sizeof(FLOAT)*3*NFBlockSize
                    +(4.0+sizeof(FLOAT))*WIDTH+4.0;
    float BytesPerSinkBlockDef = sizeof(accstruct)*NFBlockSize;

    // Set the levels assuming WIDTH Sources per Sink
    float TotalBytesPerSinkBlock = BytesPerSinkBlockDef+
            BytesPerSinkBlockWC+WIDTH*BytesPerSourceBlockWC;
    STDLOG(2, "Bytes per Block = %5.1f+%5.1f+%d*%5.1f = %5.1f\n",
            BytesPerSinkBlockDef, BytesPerSinkBlockWC, WIDTH, BytesPerSourceBlockWC,
            TotalBytesPerSinkBlock);

    // We're splitting the Buffer into WC and Def memory
    float RatioDeftoAll = BytesPerSinkBlockDef/TotalBytesPerSinkBlock;
    assert(RatioDeftoAll<=1.0);  // Guard against screwup

    // This is how many blocks we'll allocate
    MaxSinkBlocks = (MaxBufferSize-1e5)/TotalBytesPerSinkBlock;
            // Remove 100KB for alignment factors
    MaxSourceBlocks = WIDTH*MaxSinkBlocks;
    *maxsinkblocks = MaxSinkBlocks;
    *maxsourceblocks = MaxSourceBlocks;
    // The number of particles corresponding to these
    MaxSinkSize     = NFBlockSize * MaxSinkBlocks;
    MaxSourceSize   = NFBlockSize * MaxSourceBlocks;

    STDLOG(2, "Planning for %d sink and %d source blocks, each %d particles\n",
            MaxSinkBlocks, MaxSourceBlocks, NFBlockSize);

    // And then we have storage for the Pencils.
    // Here we pessimistically assume one Pencil per SinkBlock.
    MaxNSink = MaxSinkBlocks;
    MaxNSource = MaxNSink;
    assertf(MaxNSink>2*cpd, "MaxNSink = %d is too small\n", MaxNSink);
    assertf(MaxNSource>(2+WIDTH)*cpd, "MaxNSource = %d is too small\n", MaxNSource);
    // Formally, this one could be slightly larger because of the boundary
    // rows.  However, we neglected this in the Memory estimate, so we'll
    // overflow our buffer if we put it in here.  Fortunately, 1 Pencil
    // per SinkBlock assumption is really very conservative.  Only likely
    // counter-case is very low PPC sims, where one would typically have
    // a small problem that fits very easily on the GPU (i.e., MaxBufferSize
    // will have led to a MaxNSink that is > CPD**2.


    DeviceStreams = new cudaStream_t[NBuf];
    Buffers = new GPUBuffer[NBuf];
    DeviceThreads = new ThreadInfo[NBuf];
    work_queues = new tbb::concurrent_bounded_queue<GPUQueueTask>[NGPUQueues];
    
    // Start one thread per buffer
    
    // Assign threads to cores
    // TODO: make this automatic using libnuma
    //int n_socket = numa_available() == -1 ? 1 : (numa_max_node() + 1);
    //STDLOG(1, "Detected %d sockets/NUMA nodes\n", n_socket);
    //assertf(n_socket > 0, "n_socket %d less than 1\n", n_socket);

    if(UsePinnedGPUMemory < 0)
        UsePinnedGPUMemory = MaxSinkBlocks >= 10000;  // Pinning is slow, so for very small problems it's faster to use unpinned memory

    if(UsePinnedGPUMemory)
        STDLOG(1, "Allocating pinned memory\n");
    else
        STDLOG(1, "Allocating host-side memory, but not pinning because this is a small problem\n");

    pthread_barrier_t *thread_startup_barriers[NGPU];
    for(int g = 0; g < NBuf; g++){
        Buffers[g].size = MaxBufferSize;
        Buffers[g].sizeWC = Buffers[g].size*(1.0-RatioDeftoAll);
        Buffers[g].sizeDef = Buffers[g].size*RatioDeftoAll;
        Buffers[g].ready = 0;

        int core_start = ThreadCoreStart[g % NGPU];
        int core = -1;
        // If either the core start or the core count are invalid, do not bind this thread to a core
        if(core_start >= 0 && NThreadCores > 0)
            // cycle through cores, but keep each GPU within P.NGPUThreadCores of its starting core
            core = core_start + ((g/NGPU) % NThreadCores);

        if(g < NGPU){
            thread_startup_barriers[g] = new pthread_barrier_t;
            int p_ret = pthread_barrier_init(thread_startup_barriers[g], NULL, BPD);
            assertf(p_ret == 0, "pthread_barrier_init failed with value %d\n", p_ret);
        }
        
        ThreadInfo *info = DeviceThreads + g;
        info->thread_num = g;
        info->core = core;
        info->UsePinnedGPUMemory = UsePinnedGPUMemory;
        info->barrier = thread_startup_barriers[g%NGPU];
        info->queue = GPUQueueAssignments[g%NGPU];

        if(core >= 0)
            STDLOG(0, "GPU buffer thread %d (GPU %d) assigned to core %d, watching queue %d\n", g, g % NGPU, core, info->queue);
        else
            STDLOG(0, "GPU buffer thread %d (GPU %d) not bound to a core, watching queue %d\n", g, g % NGPU, info->queue);
        
        // Start one thread per buffer
        int p_retval = pthread_create(&(DeviceThreads[g].thread), NULL, QueueWatcher, info);
        assertf(p_retval == 0, "pthread_create failed with value %d\n", p_retval);

        host_alloc_bytes += Buffers[g].sizeWC + Buffers[g].sizeDef;
    }

    STDLOG(1, "Allocated %f MB host-side memory\n", host_alloc_bytes/1024./1024);

    init = 1;
}


/// This is the routine to turn off the GPU threads,
/// reset the device, and delete the host-side memory.

void GPUReset(){
    GPUQueueTask t;
    t.type = TASKTYPE_KILL;
    t.task = NULL;
    for(int i = 0; i < BPD*NGPU; i++)
        work_queues[DeviceThreads[i].queue].push(t);
    for(int i = 0; i < BPD*NGPU; i++)
        assert(pthread_join(DeviceThreads[i].thread, NULL) == 0);
    for(int i = 0; i < NGPUQueues; i++){
        assert(work_queues[i].empty());
    }

    cudaDeviceReset();

    delete[] Buffers;
    delete[] DeviceStreams;
    delete[] DeviceThreads;
    delete[] work_queues;
}


// =============== Main code for the GPU thread =================

#define CudaAllocate(ptr,size) checkCudaErrors(cudaMalloc((void **)&(ptr), size))

#define PinnedAllocate(ptr,size) if(UsePinnedGPUMemory) checkCudaErrors(cudaHostAlloc((void **)&(ptr), size, cudaHostAllocDefault)); \
                                    else assert(posix_memalign((void **)&(ptr), CACHE_LINE_SIZE, size) == 0);

#define WCAllocate(ptr,size) if(UsePinnedGPUMemory) checkCudaErrors(cudaHostAlloc((void **)&(ptr), size, cudaHostAllocWriteCombined)); \
                                    else assert(posix_memalign((void **)&(ptr), CACHE_LINE_SIZE, size) == 0);

#include "cuda_profiler_api.h"

/// Here is the top-level function for each GPU thread.
/// It allocates the space, then monitors a queue to get its work tasks
/// and then invoke the Task function.  When it receives the queue
/// signal to quit, it deletes its buffer memory.
///
/// A wrapper to GPUPencilTask that executes work units from the queue
/// Intended to be executed by several threads/buffers in parallel
/// Our original intention was to let the CUDA stream be our "queue"
/// but the copyback operation is blocking unless we pin memory
/// So this serves as a non-blocking queue.
///
/// By handling buffer memory here, we get a favorable NUMA situation.
/// However, it is also useful to have the thread running on the 
/// correct socket for its PCI slot, so that is passed in as well
/// to run set_core_affinity().

void *QueueWatcher(void *arg){
    ThreadInfo* info = (ThreadInfo *) arg;
    int assigned_device = info->thread_num;
    int n = assigned_device;                 // The buffer number
    int gpu = assigned_device % NGPU;        // The GPU device number
    int queue = info->queue;
    int assigned_core = info->core;
    if (assigned_core >= 0)
        set_core_affinity(assigned_core);
    int UsePinnedGPUMemory = info->UsePinnedGPUMemory;

    STDLOG(1,"Running GPU thread %d, core %d\n", n, assigned_core);

    checkCudaErrors(cudaSetDevice(gpu));
    if (assigned_device < NGPU) {
        // Only run these the first time a GPU is seen
        checkCudaErrors(cudaSetDeviceFlags(cudaDeviceScheduleBlockingSync));
        checkCudaErrors(cudaDeviceSetCacheConfig(cudaFuncCachePreferShared));
    }

    // Wait on the barrier: the "head" thread for each GPU must complete setting device flags before any CUDA operations can occur
    pthread_barrier_wait(info->barrier);
    STDLOG(1,"Barrier passed on stream %d\n", n);

    // Initiate the stream
    checkCudaErrors(cudaStreamCreateWithFlags(&DeviceStreams[n], cudaStreamNonBlocking));
    STDLOG(1,"GPU stream %d initiated\n", n);

    // Allocate CUDA memory
    CudaAllocate(Buffers[n].device,     Buffers[n].size);
    WCAllocate(Buffers[n].hostWC,       Buffers[n].sizeWC);
    PinnedAllocate(Buffers[n].host,     Buffers[n].sizeDef);
    STDLOG(1,"GPU thread %d memory allocated\n", n);
    // We make 2/3 of the host pinned memory WriteCombined, which is
    // good for sending data to the GPU.  The other 1/3 is normal, 
    // better for returning the data from the GPU.
    
#ifdef HAVE_LIBNUMA
    // Query the current NUMA node of the allocated buffers
    // We are using the move_pages function purely to query NUMA state, not move anything
    int page = -1, ret = 0;
    ret = move_pages(0, 1, (void **) &(Buffers[n].host), NULL, &page, 0);
    if(ret == 0)
        STDLOG(1, "Host buffer for GPU %d allocated on NUMA node %d on core %d\n", gpu, page, assigned_core);
    else
        STDLOG(1, "NUMA page query failed for GPU %d on core %d\n", gpu, assigned_core);
    ret = move_pages(0, 1, (void **) &(Buffers[n].hostWC), NULL, &page, 0);
    if(ret == 0)
        STDLOG(1, "Host write-combined buffer for GPU %d allocated on NUMA node %d on core %d\n", gpu, page, assigned_core);
    else
        STDLOG(1, "NUMA page query failed on host write-combined buffer for GPU %d on core %d\n", gpu, assigned_core);
#endif

    Buffers[n].ready = 1;

    // Main work loop: watch the queue
    while(true){
        GPUQueueTask item;
        work_queues[queue].pop(item);
        if (item.type == TASKTYPE_KILL)
            break;
        // Each thread (not work unit) is bound to a device

        if (item.type == TASKTYPE_SIC)
            GPUPencilTask(item.task, assigned_device);
    }

    STDLOG(1, "Received item signaling termination in GPU thread %d\n", n);
    
    // All done; make sure profiling info is sync'd
    checkCudaErrors(cudaStreamSynchronize(DeviceStreams[assigned_device]));
    checkCudaErrors(cudaProfilerStop());

    // Free our memory
    /* Leaking these saves some time, maybe a few seconds
    checkCudaErrors(cudaFree(Buffers[n].device));
    if (UsePinnedGPUMemory) {
        // The following cudaFreeHost lines cause this error on summitdev
        // unless jsrun is used with "--smpiargs off"
        // "CUDA Hook Library: Failed to find symbol mem_find_dreg_entries, /dev/shm/lgarrison/abacus/singlestep/singlestep: undefined symbol: __PAMI_Invalidate_region"
        // Worst case, it's probably safe to comment out these lines and leak these, since we're about to exit anyway
        checkCudaErrors(cudaFreeHost(Buffers[n].hostWC));
        checkCudaErrors(cudaFreeHost(Buffers[n].host));
    } else {
        free(Buffers[n].hostWC);
        free(Buffers[n].host);
    }*/
    
    STDLOG(1, "Terminated GPU thread %d\n", n);

    if(assigned_device < NGPU)
        delete info->barrier;  // this was allocated in GPUSetup

    return NULL;
}

#undef CudaAllocate
#undef PinnedAllocate
#undef WCAllocate


// ================= Timing the GPUs ==================

// TODO: This may need some adjustment when there's more than one task type.
// Or maybe not?  


/// Timing the GPU is tricky because there are multiple threads 
/// issuing work to it.  So we keep a global counter of all of 
/// active threads, and run a timer when that count is non-zero.
/// That gives a global throughput number.

void CUDART_CB StartThroughputTimer(cudaStream_t stream, cudaError_t status, void *data){
    assert(pthread_mutex_lock(&SetInteractionCollection::GPUTimerMutex) == 0);
    SetInteractionCollection::ActiveThreads++;
    if (SetInteractionCollection::ActiveThreads == 1)
        SetInteractionCollection::GPUThroughputTimer.Start();
    assert(pthread_mutex_unlock(&SetInteractionCollection::GPUTimerMutex) == 0);
}

void CUDART_CB MarkCompleted( cudaStream_t stream, cudaError_t status, void *data){
#ifdef CUDADIRECT
    
    assert(pthread_mutex_lock(&SetInteractionCollection::GPUTimerMutex) == 0);
    SetInteractionCollection::ActiveThreads--;
    if (SetInteractionCollection::ActiveThreads == 0)
        SetInteractionCollection::GPUThroughputTimer.Stop();
    assert(pthread_mutex_unlock(&SetInteractionCollection::GPUTimerMutex) == 0);
#endif
}


// ==================  The specific tasks =======================

/// Counting the number of direct interactions performed
__device__ unsigned long long int DI;

__global__ void getDI(unsigned long long int * h_di){
    *h_di = DI;
}

// Here is the code that knows how to execute a single Pencil task.
#include "PencilTask.cu"
