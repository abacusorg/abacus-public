//from direct_gpu_driver.cu
extern "C"
void DeviceAcceleration(accstruct* acc,int slab, int slabNP,int cpd,FLOAT eps2, int *whendone,int blocking, int * pred);

extern "C"
void SetupGPU(int* slabID, posstruct ** slabs,int * slabnp, int ** cs, int ** cnp,int cpd);

extern "C" void DeviceSlabCleanUp(accstruct * acc_slab);

extern "C" unsigned long long int DeviceGetDI();
