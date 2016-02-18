//driver for directs

#include "PencilOnPencil.cc"

#ifdef CUDADIRECT
#include "DeviceFunctions.h"
#endif

class NearFieldDriver{
    public:
        NearFieldDriver();
        ~NearFieldDriver();

        void ExecuteSlab(int slabID, int blocking);
        int SlabDone(int slabID);
        void Finalize(int slabID);
        uint64 DirectInteractions_CPU;
        uint64 DirectInteractions_GPU();
        uint64 DirectInteractions_GPU(int g);
        uint64 NSink_CPU;
        uint64 NSink_GPU();
        uint64 NSink_GPU(int g);

    private:
        int *slabcomplete;
        int * SlabNSplit;

        SetInteractionCollection *** SlabInteractionCollections; 
        int WIDTH;
        int RADIUS;
        Direct *DD;
        int NGPU;
        double GPUMemoryGB;
        uint64 *NSink_GPU_final;

        void ExecuteSlabGPU(int slabID, int blocking);
        void ExecuteSlabCPU(int slabID,int * predicate);
        void ExecuteSlabCPU(int slabID);
        void CheckGPUCPU(int slabID);
        
};



NearFieldDriver::NearFieldDriver(){
    int nthread = omp_get_max_threads();
    STDLOG(1,
            "Initializing NearFieldDriver with %d OMP threads (OMP is aware of %d procs).\n",
            nthread, omp_get_num_procs());
    DD = new Direct[nthread];
    DirectInteractions_CPU = 0;
    NSink_CPU = 0;

    assert(NFRADIUS==P.NearFieldRadius);
    slabcomplete = new int [P.cpd];
    SlabNSplit = new int[P.cpd];
    SlabInteractionCollections = new SetInteractionCollection **[P.cpd];
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

void NearFieldDriver::ExecuteSlabGPU(int slabID, int blocking){
    //calculate the required subdivisions of the slab to fit in GPU memory
    
    //Add up the sources
    uint64 NSource = 0;
    for(int i = slabID-RADIUS; i <= slabID + RADIUS; i++) NSource+=Slab->size(i);

    uint64 NSink = Slab->size(slabID);

    //Pencils are aligned to GPUBlocksize boundaries
    //In the best case we will have ~NSource/GPUBlocksize blocks
    //Worst case:
    //      We have all but one pencil GPUBlocksize + 1 long with last pencil holding the rest
    int N1 = 2*( P.cpd*(P.cpd + WIDTH) -1 );
    int N2 =std::ceil( (NSource -( NFBlockSize + 1.0) * ( P.cpd*(P.cpd + WIDTH) -1 ))/(1.0*NFBlockSize) + 1);
    if (N2 < 0) N2 = 0; //if we don't have at least GPUBlocksize particles per pencil, this calculation fails
    
    int MaxSourceBlocks = N1 + N2;
    
    N1 = 2*(P.cpd*P.cpd -1);
    N2 =std::ceil( (NSink -( NFBlockSize + 1.0) * ( P.cpd*(P.cpd + WIDTH) -1 ))/(1.0*NFBlockSize) + 1);
    if (N2 < 0) N2 = 0;
    int MaxSinkBlocks = N1 + N2; 

    double MaxMemGB = sizeof(FLOAT)* 4 * ((MaxSourceBlocks + MaxSinkBlocks) * 1e-9)*NFBlockSize;
    
    int NSplit = std::ceil(MaxMemGB/(GPUMemoryGB));
    NSplit = std::max(NSplit,NGPU);
    STDLOG(2,"Splitting slab %d into %d blocks for directs.\n",slabID,NSplit);
    SlabNSplit[slabID] = NSplit;
    
    uint64 NPTarget = NSink/NSplit;
    SlabInteractionCollections[slabID] = new SetInteractionCollection *[WIDTH*NSplit];
    int SplitPoint[NSplit];

    int yl = 0;
    for(int i =0; i < NSplit -1; i++){
        uint64 NPBlock = 0;
        int yh = yl;
        while(NPBlock < NPTarget && yh < P.cpd){
            for(int k =0; k < P.cpd; k++) NPBlock+=PP->NumberParticle(slabID,yh,k);
            yh++;
        }
        SplitPoint[i] = yh;
    }
    SplitPoint[NSplit-1] = P.cpd;

    STDLOG(2,"Direct splits for slab %d are at: ",slabID);
    for(int i = 0; i < NSplit -1; i++) STDLOG(2,"%d, ",SplitPoint[i]);
    STDLOG(2,".\n");

    for(int w = 0; w < WIDTH; w++){
        int kl =0;
        for(int n = 0; n < NSplit; n++){
            int kh = SplitPoint[n];
            SlabInteractionCollections[slabID][w*NSplit + n] = 
                new SetInteractionCollection(slabID,w,kl,kh);
            SlabInteractionCollections[slabID][w*NSplit + n]->GPUExecute(blocking);
            kl = kh;
        }
    }

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


int NearFieldDriver::SlabDone(int slab){
    slab = PP->WrapSlab(slab);

    if (slabcomplete[slab] == 0){

        int NSplit = SlabNSplit[slab];

        int complete = 1;
        for(int i = 0; i < WIDTH*NSplit; i++)
            complete = complete && SlabInteractionCollections[slab][i]->CheckCompletion();

        slabcomplete[slab] = complete;
        return complete;
    }

    return 1;
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

void NearFieldDriver::Finalize(int slab){
    slab = PP->WrapSlab(slab);

    assertf(slabcomplete[slab] != 0,
            "Finalize called for slab %d but it is not complete\n",slab);

    SetInteractionCollection ** Slices = SlabInteractionCollections[slab];
    int NSplit = SlabNSplit[slab];

    int cpd = P.cpd;
    int nfr = P.NearFieldRadius;


    for(int n = 0; n < NSplit; n++){
        for(int w = 0; w < WIDTH; w++){
            int sliceIdx = n*WIDTH + w;
            SetInteractionCollection *Slice = Slices[sliceIdx];
            int kl = Slice->K_low;
            int kh = Slice->K_high;

            #pragma omp parallel for schedule(dynamic,1)
            for(int k = kl; k < kh; k++){
                for(int j = w; j < cpd; j+=WIDTH){
                     int zmid = PP->WrapSlab(j + nfr);
                     int sinkIdx = k * cpd + zmid;
                     int SinkCount = Slice->SinkSetCount[sinkIdx];
                     int Start = Slice->SinkSetStart[sinkIdx];
                     for(int z=zmid-nfr; j < zmid+nfr;z++){
                         int CellNP = PP->NumberParticle(slab,k,z);
                         accstruct *a = PP->AccCell(slab,k,z);
                         
                         #pragma simd
                         for(int p = 0; p <CellNP; p++)
                             a[p] += Slice->SinkSetAccelerations[Start +p];

                         Start+=CellNP;
                         SinkCount-=CellNP;
                     }
                     assert(SinkCount == 0);
                }
            }
            delete Slice;
        }
    }
    delete[] Slices;
}
