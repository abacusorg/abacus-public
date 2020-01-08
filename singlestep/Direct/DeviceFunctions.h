#ifdef CUDADIRECT

//from direct_gpu_driver.cu
/*extern "C"
void DeviceAcceleration(accstruct* acc,int slab, int slabNP,int cpd,FLOAT eps2, int *whendone,int blocking, int * pred, STimer* timers);

extern "C"
void SetupGPU(int* slabID, posstruct ** slabs,int * slabnp, int ** cs, int ** cnp,int cpd, STimer* timers);

extern "C" void DeviceSlabCleanUp(accstruct * acc_slab);

extern "C" void PinSlab(FLOAT3 * slab, int id, int slabNP,int cpd);

extern "C" void UnpinSlab(FLOAT3 * slab, int id);

extern "C" void UnpinSlabs(void** slabs_to_unpin, int cpd);

extern "C" int * GetCellStart(int cpd,int s);

extern "C" int * GetCellNP(int cpd,int s);
*/
extern "C" double GetDeviceMemory();

extern "C" int GetNGPU();
extern "C" void GPUSetup(int cpd, uint64 MaxBufferSize,
	int numberGPUs, int bufferperdevice,
        int *ThreadCoreStart, int NThreadCores,
	int *maxsinkblocks, int *maxsourceblocks, int UsePinnedGPUMemory);
extern "C" void GPUReset();

extern "C" void print_gpu_mem_usage();

extern "C" void get_cuda_timers(void*);

// These are the sizes we allocate on the GPU/pinned memory
extern int MaxSinkBlocks, MaxSourceBlocks;
extern int MaxNSink, MaxNSource;
extern size_t MaxSinkSize, MaxSourceSize;

#else
// No GPU memory limits; could be RAM limits if we wanted
int MaxSinkBlocks = INT32_MAX, MaxSourceBlocks = INT32_MAX;
int MaxNSink = INT32_MAX, MaxNSource = INT32_MAX;
size_t MaxSinkSize = INT32_MAX, MaxSourceSize = INT32_MAX;
#endif
