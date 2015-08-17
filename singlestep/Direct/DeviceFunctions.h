//from direct_gpu_driver.cu
extern "C"
void DeviceAcceleration(accstruct* acc,int slab, int slabNP,int cpd,FLOAT eps2, int *whendone,int blocking, int * pred, STimer* timers);

extern "C"
void SetupGPU(int* slabID, posstruct ** slabs,int * slabnp, int ** cs, int ** cnp,int cpd, STimer* timers);

extern "C" void DeviceSlabCleanUp(accstruct * acc_slab);

extern "C" void PinSlab(FLOAT3 * slab, int id, int slabNP,int cpd);

extern "C" void UnpinSlab(FLOAT3 * slab, int id);

extern "C" void UnpinSlabs(void** slabs_to_unpin, int cpd);

extern "C" int * GetCellStart(int cpd,int s);

extern "C" int * GetCellNP(int cpd,int s);

extern "C" unsigned long long int DeviceGetDI();
extern "C" unsigned long long int DeviceGetOneDI(int g);

extern "C" int GetNGPU();
