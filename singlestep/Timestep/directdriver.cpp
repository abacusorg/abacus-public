//driver for directs

#ifdef CUDADIRECT
#include "DeviceFunctions.h"
#endif

class NearFieldDriver{
    public:
        NearFieldDriver();
        ~NearFieldDriver();

        void ExecuteSlab(int slabID, int blocking);
        int SlabDone(int slabID);
        int UnpinStatus(int slab);
        void** FindUnpinnableSlabs();
        void CleanupSlab(int slabID);
        uint64 DirectInteractions_CPU;
        uint64 DirectInteractions_GPU();
        uint64 DirectInteractions_GPU(int g);
        uint64 NSink_CPU;
        uint64 NSink_GPU();
        uint64 NSink_GPU(int g);

    private:
        int *slabcomplete;
        int WIDTH;
        int RADIUS;
        Direct *DD;
        int NGPU;
        uint64 *NSink_GPU_final;

        void ExecuteSlabGPU(int slabID, int blocking);
        void ExecuteSlabCPU(int slabID,int * predicate);
        void ExecuteSlabCPU(int slabID);
        void CheckGPUCPU(int slabID);
        
};



NearFieldDriver::NearFieldDriver(){
    int nthread = omp_get_max_threads();
    STDLOG(1,"Initializing NearFieldDriver with %d OMP threads (OMP is aware of %d procs).\n",nthread, omp_get_num_procs());
    DD = new Direct[nthread];
    DirectInteractions_CPU = 0;
    NSink_CPU = 0;

    assert(NFRADIUS==P.NearFieldRadius);
    slabcomplete = new int[P.cpd];
    for(int i = 0; i < P.cpd;i++) slabcomplete[i]=0;
    WIDTH = P.NearFieldRadius*2+1;
    RADIUS = P.NearFieldRadius;
    
#ifdef CUDADIRECT
    NGPU = GetNGPU();
    GPUTimers = new STimer[NGPU+1];
    SetupGPUTimers = new STimer[6];
    NSink_GPU_final = (uint64*) malloc(NGPU*sizeof(uint64));
    for(int g = 0; g < NGPU; g++)
        NSink_GPU_final[g] = 0;
#endif
}

NearFieldDriver::~NearFieldDriver()
{
    delete[] DD;
    delete[] slabcomplete;
#ifdef CUDADIRECT
    delete[] GPUTimers;
    delete[] SetupGPUTimers;
#endif
}

void** NearFieldDriver::FindUnpinnableSlabs(){
    // Return an array of pointers to slabs that can be unpinned
    void** slabs_to_unpin = (void**) malloc((2*P.cpd)*sizeof(void*));
    for(int s = 0; s < P.cpd; s++){
        slabs_to_unpin[2*s] = NULL;  // PosSlab
        slabs_to_unpin[2*s+1] = NULL;  // NearAccSlab
        
        // slabcomplete == 0: nothing ready to unpin
        // slabcomplete == 1: ready, but nothing has been unpinned
        // slabcomplete == 2: NearAcc has been unpinned
        // slabcomplete == 3: NearAcc and PosSlab have been unpinned; all done
        if(slabcomplete[s] == 1){
            slabs_to_unpin[2*s+1] = LBW->ReturnIDPtr(NearAccSlab,s);
            slabcomplete[s] = 2;
        }
        
        if(slabcomplete[s] == 1 || slabcomplete[s] == 2){
            slabs_to_unpin[2*s] = LBW->ReturnIDPtr(PosSlab,s);
            slabcomplete[s] = 3;
            
            // We can only unpin the PosSlab if it has been fully used as a source
            for(int j = -P.NearFieldRadius; j <= P.NearFieldRadius; j++)  {
                if( !SlabDone(s+j) ){
                    slabs_to_unpin[2*s] = NULL;
                    slabcomplete[s] = 2;
                }
            }
        }
    }
    
    return slabs_to_unpin;
}

void NearFieldDriver::ExecuteSlabGPU(int slabID, int blocking){
    #ifdef CUDADIRECT
    CellBookkeeping.Start();
    int * CellStart[WIDTH];
    int * CellNP[WIDTH];
    #ifdef AVXDIRECT
    int pred[P.cpd*P.cpd];
    #else
    int * pred = NULL;
    #endif
    int SlabNP[WIDTH];
    int SlabIDs[WIDTH];
    posstruct *SlabPos[WIDTH];

    for(int i = 0; i < WIDTH; i++){
        int sid = PP->WrapSlab(slabID+i-RADIUS);
        SlabIDs[i] = sid;
        cellinfo * SlabCI = (cellinfo *) LBW->ReturnIDPtr(CellInfoSlab,sid);
        SlabNP[i] = Slab->size(sid);
        SlabPos[i] = (posstruct *) LBW->ReturnIDPtr(PosSlab,sid);
        PinSlab(SlabPos[i],sid,SlabNP[i],P.cpd);

        CellStart[i] = GetCellStart(P.cpd,i);
        CellNP[i] = GetCellNP(P.cpd,i);

        #pragma omp parallel for schedule(dynamic,1)
        for(int y = 0; y <P.cpd*P.cpd;y++){
            (CellStart[i])[y] = SlabCI[y].startindex;
            (CellNP[i])[y] = SlabCI[y].count;
            #ifdef AVXDIRECT
            pred[y] = (SlabCI[y].count < P.GPUMinCellSinks);
            #endif

        }
    }
    CellBookkeeping.Stop();
    
    CountSinks.Start();
    // Count sink particles per GPU
    cellinfo * ThisSlabCI = (cellinfo *) LBW->ReturnIDPtr(CellInfoSlab,PP->WrapSlab(slabID));
    for(int g = 0; g < NGPU; g++){
        int ystart = (P.cpd/NGPU)*g;
        int ystart_next = (P.cpd/NGPU)*(g+1);
        
        // Last GPU gets the leftover cells
        if(g == NGPU -1)
            ystart_next = P.cpd;
        
        uint64 NSink_GPU_cell = 0;
        #pragma omp parallel for schedule(dynamic,1) reduction(+:NSink_GPU_cell)
        for(int idx = ystart*P.cpd; idx < ystart_next*P.cpd; idx++){
            if(pred != NULL && pred[idx])
                continue;
            NSink_GPU_cell += ThisSlabCI[idx].count;
        }
        
        NSink_GPU_final[g] += NSink_GPU_cell;
    }
    CountSinks.Stop();
    
    FindUnpin.Start();
    // Tell SetupGPU to unpin any Pos or NearAcc slabs that are finished
    void** slabs_to_unpin = FindUnpinnableSlabs();
    FindUnpin.Stop();
        
    SetupGPUTimers[0].Start();
    SetupGPU(SlabIDs,SlabPos,SlabNP,CellStart,CellNP,P.cpd, SetupGPUTimers);
    SetupGPUTimers[0].Stop();
    
    // Now that we've just synchronized (in SetupGPU), unpin any slabs
    GPUUnpinTimer.Start();
    UnpinSlabs(slabs_to_unpin, P.cpd);
    GPUUnpinTimer.Stop();
    
    // The "total throughput" timer
    // Only for the GPU
    GPUTimers[NGPU].Start();
    
    GPUDirectsLaunchTimer.Start();
    DeviceAcceleration((accstruct *)LBW->ReturnIDPtr(NearAccSlab,slabID), slabID, Slab->size(slabID),
            P.cpd, P.SofteningLength*P.SofteningLength, slabcomplete+slabID, blocking,pred, GPUTimers);
    GPUDirectsLaunchTimer.Stop();
    
    if(pred != NULL){
        CPUFallbackTimer.Start();
        ExecuteSlabCPU(slabID,pred);
        CPUFallbackTimer.Stop();
    }
    #endif
}

void NearFieldDriver::CheckGPUCPU(int slabID){
// Computes the CPU result to compare to the GPU result
// but does not overwrite the GPU forces.
    size_t len = Slab->size(slabID) *sizeof(accstruct);
    accstruct * a_cpu = (accstruct *)malloc(len);
    accstruct * a_tmp = (accstruct *)malloc(len);
    accstruct * a_gpu = (accstruct *) LBW->ReturnIDPtr(NearAccSlab,slabID);
    auxstruct * aux = (auxstruct *)LBW->ReturnIDPtr(AuxSlab,slabID);
    memcpy(a_tmp,a_gpu,len);  // Save the GPU result in tmp
    memset(a_gpu,0,len);
    ExecuteSlabCPU(slabID);  // Compute the CPU result
    memcpy(a_cpu,a_gpu,len);  // Save the CPU result
    memcpy(a_gpu,a_tmp,len); // Load the GPU result
    #ifdef DOUBLEPRECISION
    FLOAT target = 1e-10;
    #else
    FLOAT target = 1e-1;
    #endif

    for(int i = 0; i < Slab->size(slabID);i++){
        accstruct ai_g = a_gpu[i];
        accstruct ai_c = a_cpu[i];
        FLOAT delta =2* (ai_g-ai_c).norm()/(ai_g.norm() + ai_c.norm());
        if(!(delta < target)){
            printf("Error in slab %d:\n\ta_gpu[%d]: (%5.4f,%5.4f,%5.4f)\n\ta_cpu[%d]: (%5.4f,%5.4f,%5.4f)\n\tdelta:%f\n",
                    slabID,i,ai_g.x,ai_g.y,ai_g.z,i,ai_c.x,ai_c.y,ai_c.z,delta);
            //assert(delta < target);
        }
        assert(isfinite(ai_g.x));
        assert(isfinite(ai_g.y));
        assert(isfinite(ai_g.z));

        assert(isfinite(ai_c.x));
        assert(isfinite(ai_c.y));
        assert(isfinite(ai_c.z));
    }
    free(a_cpu);
    free(a_tmp);
}

void NearFieldDriver::ExecuteSlabCPU(int slabID){
    ExecuteSlabCPU(slabID,(int *)NULL);
}


void NearFieldDriver::ExecuteSlabCPU(int slabID, int * predicate){
    uint64 DI_slab = 0;
    uint64 NSink_CPU_slab = 0;
    #pragma omp parallel for schedule(dynamic,1) reduction(+:DI_slab,NSink_CPU_slab)
    for(int y = 0; y < P.cpd; y++){
        int g = omp_get_thread_num();
        //STDLOG(1,"Executing directs on pencil y=%d in slab %d, in OMP thread %d on CPU %d (nprocs: %d)\n", y, slabID, g, sched_getcpu(), omp_get_num_procs());
        for(int z = 0; z < P.cpd; z++){
            if(predicate != NULL && !predicate[y*P.cpd +z]) continue;
            posstruct * sink_pos = PP->PosCell(slabID,y,z);
            accstruct * sink_acc = PP->NearAccCell(slabID,y,z);
            uint64 np_sink = PP->NumberParticle(slabID,y,z);
            NSink_CPU_slab += np_sink;
            if (np_sink == 0) continue;
            for(int i = slabID - RADIUS; i <= slabID + RADIUS; i++){
                for(int j = y - RADIUS; j <= y + RADIUS; j++){
                    for(int k = z - RADIUS; k <= z + RADIUS; k++){
                        posstruct *source_pos = PP->PosCell(i,j,k);
                        FLOAT3 delta = PP->CellCenter(slabID,y,z)-PP->CellCenter(i,j,k);
                        uint64 np_source = PP->NumberParticle(i,j,k);
                        if(np_source >0) DD[g].AVXExecute(sink_pos,source_pos,np_sink,np_source,
                                delta,P.SofteningLength*P.SofteningLength,sink_acc);
                        DI_slab+=np_sink*np_source;
                    }
                }
            }
        }
    }
    DirectInteractions_CPU += DI_slab;
    NSink_CPU += NSink_CPU_slab;
}

void NearFieldDriver::ExecuteSlab(int slabID, int blocking){
    #ifdef CUDADIRECT
    if(!P.ForceCPU){
        ExecuteSlabGPU(slabID,blocking);
        if(P.ForceOutputDebug) CheckGPUCPU(slabID);
    }
    else{
    #else
        {
    #endif
        ExecuteSlabCPU(slabID);
        slabcomplete[slabID] = 1;
        }
}

void NearFieldDriver::CleanupSlab(int slab){
#ifdef CUDADIRECT
    if(!P.ForceCPU){
        UnpinSlab((accstruct *)LBW->ReturnIDPtr(PosSlab,slab),slab);
    }
#endif
}


int NearFieldDriver::SlabDone(int slab){
    slab = PP->WrapSlab(slab);
    return slabcomplete[slab] != 0;
}

int NearFieldDriver::UnpinStatus(int slab){
    slab = PP->WrapSlab(slab);
    return slabcomplete[slab];
}

uint64 NearFieldDriver::DirectInteractions_GPU(){
    #ifdef CUDADIRECT
    return DeviceGetDI();
    #endif
    return 0;
}

uint64 NearFieldDriver::DirectInteractions_GPU(int g){
    #ifdef CUDADIRECT
    assert(g < NGPU);
    return DeviceGetOneDI(g);
    #endif
    return 0;
}

uint64 NearFieldDriver::NSink_GPU(){
    #ifdef CUDADIRECT
    uint64 n = 0;
    for(int g = 0; g < NGPU; g++)
        n += NSink_GPU_final[g];
    return n;
    #endif
    return 0;
}

uint64 NearFieldDriver::NSink_GPU(int g){
    #ifdef CUDADIRECT
    assert(g < NGPU);
    return NSink_GPU_final[g];
    #endif
    return 0;
}
