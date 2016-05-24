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
        uint64 *DirectInteractions_GPU;
        uint64 TotalDirectInteractions_GPU;
        uint64 NSink_CPU;

        // Total timings from building and running the SetInteractionCollections
        // These are gathered in Finalize(slab)
        double  Construction;
        double      FillSinkLists;
        double          CountSinks;
        double          CalcSinkBlocks;
        double          FillSinks;
        double      FillSourceLists;
        double          CountSources;
        double          CalcSourceBlocks;
        double          FillSources;
        double      FillInteractionList;
        
        double LaunchDeviceKernels;
        double *DeviceCopyTimes;
        double *DeviceExecutionTimes;
        double *DeviceCopybackTimes;
        double *DeviceTotalTimes;
        
        STimer CalcSplitDirects;
        STimer SICExecute;
        STimer CPUFallbackTimer;
        STimer FinalizeTimer;

    private:
        int *slabcomplete;
        int * SlabNSplit;

        SetInteractionCollection *** SlabInteractionCollections; 
        int WIDTH;
        int RADIUS;
        Direct *DD;
        int NGPU;
        double GPUMemoryGB;
        int MaxSourceBlocks;
        int MaxSinkBlocks;
        uint64 *NSink_GPU_final;
        pthread_mutex_t * CellLock;

        void ExecuteSlabGPU(int slabID, int blocking);
        void ExecuteSlabCPU(int slabID,int * predicate);
        void ExecuteSlabCPU(int slabID);
        void CheckGPUCPU(int slabID);
        void CheckInteractionList(int slabID);
        
};



NearFieldDriver::NearFieldDriver() :         
        Construction{0},
              FillSinkLists{0},
                  CountSinks{0},
                  CalcSinkBlocks{0},
                  FillSinks{0},
              FillSourceLists{0},
                  CountSources{0},
                  CalcSourceBlocks{0},
                  FillSources{0},
        FillInteractionList{0},
        LaunchDeviceKernels{0},
        TotalDirectInteractions_GPU{0}
{
    int nthread = omp_get_max_threads();
    STDLOG(1,
            "Initializing NearFieldDriver with %d OMP threads (OMP is aware of %d procs).\n",
            nthread, omp_get_num_procs());
    DD = new Direct[nthread];
    DirectInteractions_CPU = 0;
    DirectInteractions_GPU = new uint64[NGPU*DirectBPD]();

    assert(NFRADIUS==P.NearFieldRadius);
    slabcomplete = new int [P.cpd];
    SlabNSplit = new int[P.cpd];
    SlabInteractionCollections = new SetInteractionCollection **[P.cpd];
    for(int i = 0; i < P.cpd;i++) slabcomplete[i]=0;
    WIDTH = P.NearFieldRadius*2+1;
    RADIUS = P.NearFieldRadius;
    CellLock = (pthread_mutex_t *) malloc(P.cpd*P.cpd*sizeof(pthread_mutex_t));
    for(int i =0; i < P.cpd*P.cpd; i ++) pthread_mutex_init(&CellLock[i],NULL);
    
#ifdef CUDADIRECT
    NGPU = GetNGPU();
    GPUMemoryGB = GetDeviceMemory();
    GPUMemoryGB = std::min(5.0e-9*P.np*sizeof(FLOAT3),GPUMemoryGB/(DirectBPD));
    size_t BlockSizeBytes = sizeof(FLOAT) *3 * NFBlockSize;
    MaxSinkBlocks = floor(1e9 * GPUMemoryGB/(6*BlockSizeBytes));
    MaxSourceBlocks = 5 * MaxSinkBlocks;
    STDLOG(1,"Initializing GPU with %f x10^6 sink blocks and %f x10^6 source blocks\n",
            MaxSinkBlocks/1e6,MaxSourceBlocks/1e6);
    GPUSetup(P.cpd,P.cpd,MaxSinkBlocks,MaxSourceBlocks,DirectBPD);
    SICExecute.Clear();

    NSink_GPU_final = (uint64*) malloc(NGPU*sizeof(uint64));
    for(int g = 0; g < NGPU; g++)
        NSink_GPU_final[g] = 0;
    
    DeviceCopyTimes = new double[NGPU*DirectBPD]();
    DeviceExecutionTimes = new double[NGPU*DirectBPD]();
    DeviceCopybackTimes = new double[NGPU*DirectBPD]();
    DeviceTotalTimes = new double[NGPU*DirectBPD]();
#endif
}

NearFieldDriver::~NearFieldDriver()
{
    delete[] DD;
    delete[] slabcomplete;
#ifdef CUDADIRECT
    delete[] DirectInteractions_GPU;
    delete[] DeviceCopyTimes;
    delete[] DeviceExecutionTimes;
    delete[] DeviceCopybackTimes;
    delete[] DeviceTotalTimes;
#endif
    for(int i =0; i < P.cpd*P.cpd; i++)
        pthread_mutex_destroy(&CellLock[i]);
    free(CellLock);
}

void NearFieldDriver::ExecuteSlabGPU(int slabID, int blocking){
    //calculate the required subdivisions of the slab to fit in GPU memory
    CalcSplitDirects.Start();
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
    
    int SourceBlocks = N1 + N2;
    
    N1 = 2*(P.cpd*P.cpd -1);
    N2 =std::ceil( (NSink -( NFBlockSize + 1.0) * ( P.cpd*(P.cpd + WIDTH) -1 ))/(1.0*NFBlockSize) + 1);
    if (N2 < 0) N2 = 0;
    int SinkBlocks = N1 + N2;
    
    int NSplitSource = std::ceil((1.0*SourceBlocks)/MaxSourceBlocks);
    int NSplitSink = std::ceil((1.0*SinkBlocks)/MaxSinkBlocks);
    int NSplit = std::max(NSplitSource,NSplitSink);
    STDLOG(1,"Splitting slab %d into %d blocks for directs.\n",slabID,NSplit);
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

    STDLOG(1,"Direct splits for slab %d are at: \n",slabID);
    for(int i = 0; i < NSplit -1; i++) STDLOG(2,"\t\t\t\t %d: %d \n",i,SplitPoint[i]);
    
    CalcSplitDirects.Stop();

    
    for(int w = 0; w < WIDTH; w++){
        int kl =0;
        for(int n = 0; n < NSplit; n++){
            int kh = SplitPoint[n];
            // The construction and execution are timed internally in each SIC then reduced in Finalize(slab)
            SlabInteractionCollections[slabID][w*NSplit + n] = 
                new SetInteractionCollection(slabID,w,kl,kh);
            STDLOG(1,"Executing directs for slab %d w = %d k: %d - %d\n",slabID,w,kl,kh);
            SICExecute.Start();
            SlabInteractionCollections[slabID][w*NSplit + n]->GPUExecute(blocking);
            SICExecute.Stop();
            kl = kh;
        }
    }

}

void NearFieldDriver::CheckGPUCPU(int slabID){
    return;
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
                                delta,WriteState.SofteningLength*WriteState.SofteningLength,sink_acc);
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

void NearFieldDriver::Finalize(int slab){
    FinalizeTimer.Start();
    slab = PP->WrapSlab(slab);

    #ifdef CUDADIRECT

    assertf(SlabDone(slab) != 0,
            "Finalize called for slab %d but it is not complete\n",slab);

    if(P.ForceOutputDebug)
        CheckInteractionList(slab);

    SetInteractionCollection ** Slices = SlabInteractionCollections[slab];
    int NSplit = SlabNSplit[slab];

    int cpd = P.cpd;
    int nfr = P.NearFieldRadius;
    int width = 2*nfr +1;


    for(int n = 0; n < NSplit; n++){
        for(int w = 0; w < WIDTH; w++){
            int sliceIdx = w*NSplit + n;
            SetInteractionCollection *Slice = Slices[sliceIdx];
            int kl = Slice->K_low;
            int kh = Slice->K_high;
            int Nj = (cpd - w)/width;
            int jr = (cpd - w)%width;
            if(Nj * width + w < cpd)
                Nj++;
            //Slice->PrintInteractions();
            
            // Gather timings from each SlabInteractionCollection
            Construction +=Slice->Construction.Elapsed();
                FillSinkLists+=Slice->FillSinkLists.Elapsed();
                    CountSinks+=Slice->CountSinks.Elapsed();
                    CalcSinkBlocks+=Slice->CalcSinkBlocks.Elapsed();
                    FillSinks+=Slice->FillSinks.Elapsed();          
                FillSourceLists+=Slice->FillSourceLists.Elapsed();
                    CountSources+=Slice->CountSources.Elapsed();
                    CalcSourceBlocks+=Slice->CalcSinkBlocks.Elapsed();
                    FillSources+=Slice->FillSources.Elapsed();
                FillInteractionList+=Slice->FillInteractionList.Elapsed();
            LaunchDeviceKernels += Slice->LaunchDeviceKernels.Elapsed();
            
            int g = Slice->AssignedDevice;
            DeviceCopyTimes[g] += Slice->CopyTime;
            DeviceExecutionTimes[g] += Slice->ExecutionTime;
            DeviceCopybackTimes[g] += Slice->CopybackTime;
            DeviceTotalTimes[g] += Slice->TotalTime;
            
            DirectInteractions_GPU[g] += Slice->DirectTotal;
            TotalDirectInteractions_GPU += Slice->DirectTotal;

            int NThread = omp_get_max_threads();

            #pragma omp parallel for schedule(dynamic,1)
            for(int sinkIdx = 0; sinkIdx < Slice->NSinkList; sinkIdx++){
                int SinkCount = Slice->SinkSetCount[sinkIdx];
                int j = sinkIdx%Nj;
                int k = sinkIdx/Nj + Slice->K_low;
                int jj = w + j * width;
                int zmid = PP->WrapSlab(jj + P.NearFieldRadius);
                int Start = Slice->SinkSetStart[sinkIdx];


                for(int z=zmid-nfr; z <= zmid+nfr;z++){
                    int cid = PP->CellID(k,z);
                    pthread_mutex_lock(&CellLock[cid]);
                    int CellNP = PP->NumberParticle(slab,k,z);
                    accstruct *a = PP->NearAccCell(slab,k,z);

                    if(P.ForceOutputDebug == 1){
                        for(int p =0; p < CellNP; p++){
                            assert(isfinite(Slice->SinkSetAccelerations[Start +p].x));
                            assert(isfinite(Slice->SinkSetAccelerations[Start +p].y));
                            assert(isfinite(Slice->SinkSetAccelerations[Start +p].z));
                        }
                    }

                    #pragma simd
                    for(int p = 0; p <CellNP; p++)
                        a[p] += Slice->SinkSetAccelerations[Start +p];
                    pthread_mutex_unlock(&CellLock[cid]);

                    Start+=CellNP;
                    SinkCount-=CellNP;
                }
                assert(SinkCount == 0);
            }

            delete Slice;
        }
    }
    delete[] Slices;
    if(P.ForceOutputDebug) CheckGPUCPU(slab);
    #endif
    FinalizeTimer.Stop();
}

#include <vector>
#include <algorithm>
void NearFieldDriver::CheckInteractionList(int slab){

    vector<int> ** il = new std::vector<int> *[P.cpd*P.cpd];

    for(int i = 0; i < P.cpd*P.cpd; i++){
        il[i] = new vector<int>();
        il[i]->reserve(WIDTH*WIDTH*WIDTH);
    }

    SetInteractionCollection ** Slices = SlabInteractionCollections[slab];
    int NSplit = SlabNSplit[slab];

    for(int s = 0; s < WIDTH*NSplit; s++)
        Slices[s]->AddInteractionList(il);


    //check number of interactions
    for(int i = 0; i < P.cpd*P.cpd; i++){
        if (il[i]->size() != WIDTH*WIDTH*WIDTH){
            printf("Error: cell %d has %d =/= %d interactions\n",i,il[i]->size(),WIDTH*WIDTH*WIDTH);
            assert(il[i]->size() == WIDTH*WIDTH*WIDTH);
        }
    }

    //Now we do a full check of the interaction list
    
    //build a correct interaction list for comparison
    
    vector<int> ** il_test = new std::vector<int> *[P.cpd*P.cpd];
    for(int i = 0; i < P.cpd*P.cpd; i++){
        il_test[i] = new vector<int>();
        il_test[i]->reserve(WIDTH*WIDTH*WIDTH);
    }
    for(int y = 0; y < P.cpd; y++){
        for(int z = 0; z < P.cpd; z++){
            for(int i = slab - RADIUS; i <= slab + RADIUS; i++){
                for(int j = y - RADIUS; j <= y + RADIUS; j++){
                    for(int k = z - RADIUS; k <= z + RADIUS; k++){
                        int sinkCellId = y*P.cpd + z;
                        int sourceCellId = P.cpd*P.cpd * PP->WrapSlab(i) +
                            P.cpd*PP->WrapSlab(j) +
                            PP->WrapSlab(k);
                        il_test[sinkCellId]->push_back(sourceCellId);
                    }
                }
            }
        }
    }
    
    
    //sort the lists and check them
    for(int i=0; i <P.cpd*P.cpd; i++){
        std::sort(il[i]->begin(),il[i]->end());
        std::sort(il_test[i]->begin(),il_test[i]->end());

        assert(il[i]->size() == il_test[i]->size());

        for(int j = 0; j < il[i]->size(); j ++)
                assertf(il[i]->at(j) == il_test[i]->at(j), "Interaction list %d at %d failed in slab %d\n", i, j, slab);
        delete il[i];
        delete il_test[i];
    }

    delete[] il;
    delete[] il_test;
    
    //STDLOG(1,"Checking the interaction list for slab %d passed.\n", slab);
    //Checking the interaction list is time-consuming
    //it should never be done in an actual run
}
