//driver for directs


class NearFieldDriver{
    public:
        NearFieldDriver();
        ~NearFieldDriver();

        void ExecuteSlab(int slabID, int blocking);
        int SlabDone(int slabID);
        void CleanupSlab(int slabID);
        uint64 DirectInteractions_CPU;
        uint64 DirectInteractions_GPU();

    private:
        int *slabcomplete;
        int WIDTH;
        int RADIUS;
        Direct *DD;

        void ExecuteSlabGPU(int slabID, int blocking);
        void ExecuteSlabCPU(int slabID,int * predicate);
        void ExecuteSlabCPU(int slabID);
        void CheckGPUCPU(int slabID);
        
};



NearFieldDriver::NearFieldDriver(){
    DD = new Direct[omp_get_max_threads()];
    DirectInteractions_CPU = 0;

    assert(NFRADIUS==P.NearFieldRadius);
    slabcomplete = new int[P.cpd];
    for(int i = 0; i < P.cpd;i++) slabcomplete[i]=0;
    WIDTH = P.NearFieldRadius*2+1;
    RADIUS = P.NearFieldRadius;
}

NearFieldDriver::~NearFieldDriver()
{
    delete[] DD;
    delete[] slabcomplete;
}



#ifdef CUDADIRECT
#include "DeviceFunctions.h"
#endif

void NearFieldDriver::ExecuteSlabGPU(int slabID, int blocking){
    #ifdef CUDADIRECT
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
    SetupGPU(SlabIDs,SlabPos,SlabNP,CellStart,CellNP,P.cpd);

    DeviceAcceleration((accstruct *)LBW->ReturnIDPtr(NearAccSlab,slabID), slabID, Slab->size(slabID),
            P.cpd, P.SofteningLength*P.SofteningLength, slabcomplete+slabID, blocking,pred);
    if(pred != NULL) ExecuteSlabCPU(slabID,pred);
    #endif
}

void NearFieldDriver::CheckGPUCPU(int slabID){
// Computes the CPU result to compare to the GPU result
// but does not overwrite the GPU forces.
#ifndef DIRECTSPLINE // Spline softening isn't implemented for the CPU yet, so it doesn't make sense to compare the two
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
    }
    free(a_cpu);
    free(a_tmp);
#endif
}

void NearFieldDriver::ExecuteSlabCPU(int slabID){
    ExecuteSlabCPU(slabID,(int *)NULL);
}


void NearFieldDriver::ExecuteSlabCPU(int slabID, int * predicate){
    uint64 DI_slab = 0;
    #pragma omp parallel for schedule(dynamic,1) reduction(+:DI_slab)
    for(int y = 0; y < P.cpd; y++){
        int g = omp_get_thread_num();
        for(int z = 0; z < P.cpd; z++){
            if(predicate != NULL && !predicate[y*P.cpd +z]) continue;
            posstruct * sink_pos = PP->PosCell(slabID,y,z);
            accstruct * sink_acc = PP->NearAccCell(slabID,y,z);
            uint64 np_sink = PP->NumberParticle(slabID,y,z);
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
    DirectInteractions_CPU+=DI_slab;
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
    if(!P.ForceCPU) UnpinSlab((accstruct *)LBW->ReturnIDPtr(PosSlab,slab),slab);
#endif
}


int NearFieldDriver::SlabDone(int slab){
    slab = PP->WrapSlab(slab);
    if (slabcomplete[slab] == 1) {
        #ifdef CUDADIRECT
        if(!P.ForceCPU) {
            DeviceSlabCleanUp((accstruct *)LBW->ReturnIDPtr(NearAccSlab,slab));
        }

        #endif
        slabcomplete[slab] = 2;
    }
    return slabcomplete[slab] != 0;
}

uint64 NearFieldDriver::DirectInteractions_GPU(){
    #ifdef CUDADIRECT
    return DeviceGetDI();
    #endif
    return 0;
}
