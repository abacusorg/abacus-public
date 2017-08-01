//driver for directs

#include "StructureOfLists.cc"
#include "SetInteractionCollection.hh"

#include "DeviceFunctions.h"

class NearFieldDriver{
    public:
        NearFieldDriver();
        ~NearFieldDriver();
    
        double SofteningLength;  // Effective Plummer length, used for timestepping.  Unit-box units.
        double SofteningLengthInternal;  // The equivalent length for the current softening technique.  Unit-box units.
        FLOAT eps; // Some power of SofteningLengthInternal, like ^2 or ^3, precomputed for the softening kernel

        void ExecuteSlab(int slabID, int blocking);
        int SlabDone(int slabID);
        void Finalize(int slabID);
        uint64 DirectInteractions_CPU;
        uint64 *DirectInteractions_GPU;
        uint64 TotalDirectInteractions_GPU = 0;
        uint64 NSink_CPU;

        // Total timings from building and running the SetInteractionCollections
        // These are gathered in Finalize(slab)
        double  Construction = 0;
        double      FillSinkLists = 0;
        double          CountSinks = 0;
        double          CalcSinkBlocks = 0;
        double          AllocAccels = 0;
        double      FillSourceLists = 0;
        double          CountSources = 0;
        double          CalcSourceBlocks = 0;
        double      FillInteractionList = 0;
        
        double DeviceThreadTimer = 0;
        double     LaunchDeviceKernels = 0;
        double     FillSinks = 0;
        double     FillSources = 0;
        double     WaitForResult = 0;
        double     CopyAccelFromPinned = 0;

    
        double *GB_to_device, *GB_from_device;
        uint64 *DeviceSinks;
        uint64 *DeviceSources;
        
        STimer CalcSplitDirects;
        STimer SICExecute;
        STimer CPUFallbackTimer;
        STimer FinalizeTimer;
        STimer FinalizeBookkeeping;
        STimer CopyPencilToSlab;
        STimer ZeroAccel;

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
        int MinSplits;
        uint64 *NSink_GPU_final;
        pthread_mutex_t * CellLock;

        void ExecuteSlabGPU(int slabID, int blocking);
        void ExecuteSlabCPU(int slabID,int * predicate);
        void ExecuteSlabCPU(int slabID);
        void CheckGPUCPU(int slabID);
        void CheckInteractionList(int slabID);
        int GetNSplit(uint64 NSource, uint64 NSink);
        
};



NearFieldDriver::NearFieldDriver() :
        SofteningLength{WriteState.SofteningLength/P.BoxSize},
        SofteningLengthInternal{WriteState.SofteningLengthInternal/P.BoxSize},
            
        #ifdef DIRECTCUBICSPLINE
        eps{(FLOAT) (1./SofteningLengthInternal)}
        #elif defined DIRECTCUBICPLUMMER
        eps{(FLOAT) (SofteningLengthInternal*SofteningLengthInternal*SofteningLengthInternal)}
        #elif defined DIRECTSINGLESPLINE
        eps{(FLOAT) (1./(SofteningLengthInternal*SofteningLengthInternal))}
        #else
        eps{(FLOAT) (SofteningLengthInternal*SofteningLengthInternal)}
        #endif
{
    assertf(isfinite(eps), "Infinite eps!  Softening length too small for this precision?\n");
    
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
    CellLock = (pthread_mutex_t *) malloc(P.cpd*P.cpd*sizeof(pthread_mutex_t));
    for(int i =0; i < P.cpd*P.cpd; i ++) pthread_mutex_init(&CellLock[i],NULL);
    
#ifdef CUDADIRECT
    NGPU = GetNGPU();
    GPUMemoryGB = GetDeviceMemory();
    GPUMemoryGB = std::min(5.0e-9*P.np*sizeof(FLOAT3),GPUMemoryGB/(DirectBPD));
    STDLOG(1, "Using %f GB of GPU memory (per device)\n", GPUMemoryGB);
    size_t BlockSizeBytes = sizeof(FLOAT) *3 * NFBlockSize;
    MaxSinkBlocks = floor(1e9 * GPUMemoryGB/(6*BlockSizeBytes));
    MaxSourceBlocks = WIDTH * MaxSinkBlocks;
    STDLOG(1,"Initializing GPU with %f x10^6 sink blocks and %f x10^6 source blocks\n",
            MaxSinkBlocks/1e6,MaxSourceBlocks/1e6);
    
    // We want to compute maxkwidth for GPUSetup
    // maxkwidth will occur when we have the fewest splits
    // The fewest splits will occur when we are operating on the smallest slabs
    MinSplits = GetNSplit(WIDTH*Slab->min, Slab->min);
    MinSplits = ceil(.8*MinSplits);  // fudge factor to account for uneven slabs
    // This may not account for unequal splits, though.  Unless we really need to save GPU memory, just use maxkwidth=cpd
    //MinSplits = 1;
    STDLOG(1,"MinSplits = %d\n", MinSplits);
    if(MinSplits*WIDTH < omp_get_num_threads())
        STDLOG(1, "*** WARNING: MinSplits*WIDTH is less than the number of threads! Finalize() might be inefficient\n");
    GPUSetup(P.cpd, std::ceil(1.*P.cpd/MinSplits), MaxSinkBlocks, MaxSourceBlocks, DirectBPD, P.GPUThreadCoreStart, P.NGPUThreadCores);
    SICExecute.Clear();

    NSink_GPU_final = (uint64*) malloc(NGPU*sizeof(uint64));
    for(int g = 0; g < NGPU; g++)
        NSink_GPU_final[g] = 0;
    
    DirectInteractions_GPU = new uint64[NGPU*DirectBPD]();
    GB_to_device = new double[NGPU*DirectBPD]();
    GB_from_device = new double[NGPU*DirectBPD]();
    DeviceSinks = new uint64[NGPU*DirectBPD]();
    DeviceSources = new uint64[NGPU*DirectBPD]();
#endif
}

NearFieldDriver::~NearFieldDriver()
{
    delete[] DD;
    delete[] slabcomplete;
    delete[] SlabNSplit;
    delete[] SlabInteractionCollections;
#ifdef CUDADIRECT
    free(NSink_GPU_final);
    delete[] DirectInteractions_GPU;
    delete[] GB_to_device;
    delete[] GB_from_device;
    delete[] DeviceSinks;
    delete[] DeviceSources;
    
    GPUReset();
#endif
    for(int i =0; i < P.cpd*P.cpd; i++)
        pthread_mutex_destroy(&CellLock[i]);
    free(CellLock);
}

// Compute the number of splits for a given number of sources and sinks
int NearFieldDriver::GetNSplit(uint64 NSource, uint64 NSink){
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
    const int NSplit = std::max(NSplitSource,NSplitSink);
    
    return NSplit;
}

void NearFieldDriver::ExecuteSlabGPU(int slabID, int blocking){
    //calculate the required subdivisions of the slab to fit in GPU memory
    CalcSplitDirects.Start();
    //Add up the sources
    uint64 NSource = 0;
    for(int i = slabID-RADIUS; i <= slabID + RADIUS; i++) NSource+=Slab->size(i);

    uint64 NSink = Slab->size(slabID);

    const int NSplit = GetNSplit(NSource, NSink);
    assertf(NSplit >= MinSplits, "NSplit (%d) is less than MinSplits (%d)\n", NSplit, MinSplits);
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
        yl = yh;
    }
    SplitPoint[NSplit-1] = P.cpd;

    if (NSplit > 1) {
        STDLOG(1,"Direct splits for slab %d are at: \n",slabID);
        for(int i = 0; i < NSplit -1; i++)
            STDLOG(1,"\t\t\t\t %d: %d \n",i,SplitPoint[i]);
    }
    
    CalcSplitDirects.Stop();

    // If we thread over y-splits, would that help with NUMA locality?
    for(int w = 0; w < WIDTH; w++){
        int kl =0;
        for(int n = 0; n < NSplit; n++){
            int kh = SplitPoint[n];
            // The construction and execution are timed internally in each SIC then reduced in Finalize(slab)
            SlabInteractionCollections[slabID][w*NSplit + n] = 
                new SetInteractionCollection(slabID,w,kl,kh);
            
            // Check that we have enough blocks
            int NSinkBlocks = SlabInteractionCollections[slabID][w*NSplit + n]->NSinkBlocks;
            int NSourceBlocks = SlabInteractionCollections[slabID][w*NSplit + n]->NSourceBlocks;
            assertf(NSinkBlocks <= MaxSinkBlocks,
                    "NSinkBlocks (%d) is larger than MaxSinkBlocks (%d)\n", NSinkBlocks, MaxSinkBlocks);
            assertf(NSourceBlocks <= MaxSourceBlocks,
                    "NSourceBlocks (%d) is larger than MaxSourceBlocks (%d)\n", NSourceBlocks, MaxSourceBlocks);
            
            STDLOG(1,"Executing directs for slab %d w = %d k: %d - %d\n",slabID,w,kl,kh);
            SICExecute.Start();
            SlabInteractionCollections[slabID][w*NSplit + n]->GPUExecute(blocking);
            //SlabInteractionCollections[slabID][w*NSplit + n]->CPUExecute();
            SICExecute.Stop();
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
    
    STDLOG(1,"GPU-CPU comparison passed for slab %d\n", slabID);
}

void NearFieldDriver::ExecuteSlabCPU(int slabID){
    ExecuteSlabCPU(slabID,(int *)NULL);
}


void NearFieldDriver::ExecuteSlabCPU(int slabID, int * predicate){
    CPUFallbackTimer.Start();
    if(!LBW->IDPresent(NearAccSlab, slabID))
        LBW->AllocateArena(NearAccSlab,slabID);
    ZeroAcceleration(slabID,NearAccSlab);
    
    #ifdef DIRECTSINGLESPLINE
    FLOAT inv_eps3 = 1./(SofteningLengthInternal*SofteningLengthInternal*SofteningLengthInternal);
    #endif

    uint64 DI_slab = 0;
    uint64 NSink_CPU_slab = 0;
    #pragma omp parallel for schedule(dynamic,1) reduction(+:DI_slab,NSink_CPU_slab)
    for(int y = 0; y < P.cpd; y++){
        int g = omp_get_thread_num();
        //STDLOG(1,"Executing directs on pencil y=%d in slab %d, in OMP thread %d on CPU %d (nprocs: %d)\n", y, slabID, g, sched_getcpu(), omp_get_num_procs());
        for(int z = 0; z < P.cpd; z++){
            if(predicate != NULL && !predicate[y*P.cpd +z]) continue;
            
            // We can use PosCell instead of PosXYZ here because it's for the sink positions
            posstruct * sink_pos = PP->PosCell(slabID,y,z);
            accstruct * sink_acc = PP->NearAccCell(slabID,y,z);
            uint64 np_sink = PP->NumberParticle(slabID,y,z);
            NSink_CPU_slab += np_sink;
            if (np_sink == 0) continue;
            for(int i = slabID - RADIUS; i <= slabID + RADIUS; i++){
                for(int j = y - RADIUS; j <= y + RADIUS; j++){
                    for(int k = z - RADIUS; k <= z + RADIUS; k++){
                        uint64 np_source = PP->NumberParticle(i,j,k);
                        
                        // We assume that PosXYZ is used for all sources
                        // This lets us drift PosSlab much earlier
                        List3<FLOAT> source_pos_xyz = PP->PosXYZCell(i,j,k);
                        posstruct *source_pos = new posstruct[np_source];
                        for(uint64 ii = 0; ii < np_source; ii++){
                            source_pos[ii].x = source_pos_xyz.X[ii];
                            source_pos[ii].y = source_pos_xyz.Y[ii];
                            source_pos[ii].z = source_pos_xyz.Z[ii];
                        }
                        
                        FLOAT3 delta = PP->CellCenter(slabID,y,z)-PP->CellCenter(i,j,k);
                        if(np_source >0) DD[g].AVXExecute(sink_pos,source_pos,np_sink,np_source,
                                delta,eps,sink_acc);
                        delete[] source_pos;
                        
                        DI_slab += np_sink*np_source;
                    }
                }
            }
            #ifdef DIRECTSINGLESPLINE
            for(uint64 i = 0; i < np_sink; i++)
                sink_acc[i] *= inv_eps3;
            #endif
        }
    }
    DirectInteractions_CPU += DI_slab;
    NSink_CPU += NSink_CPU_slab;
    CPUFallbackTimer.Stop();
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
    FinalizeBookkeeping.Start();

    assertf(SlabDone(slab) != 0,
            "Finalize called for slab %d but it is not complete\n",slab);
            
    #ifdef DIRECTSINGLESPLINE
    FLOAT inv_eps3 = 1./(SofteningLengthInternal*SofteningLengthInternal*SofteningLengthInternal);
    #endif

    if(P.ForceOutputDebug)
        CheckInteractionList(slab);
    
    LBW->AllocateArena(NearAccSlab,slab);

    SetInteractionCollection ** Slices = SlabInteractionCollections[slab];
    int NSplit = SlabNSplit[slab];

    int cpd = P.cpd;
    int nfr = P.NearFieldRadius;
    int width = 2*nfr +1;

    // Ideally we would zero-initialize all cell acclerations then add accelerations to that
    // but zeroing every cell is really expensive, so we instead *set* the acceleration
    // instead of add if it's the first time touching the cell
    int *CellAccInit = new int[cpd*cpd]();

    // There's an annoying number of reductions here, would be a lot of overhead to do in parallel
    // Let's do them serially here; hopefully this is quick.
    for(int sliceIdx = 0; sliceIdx < NSplit*WIDTH; sliceIdx++){
        SetInteractionCollection *Slice = Slices[sliceIdx];

        Construction +=Slice->Construction.Elapsed();
            FillSinkLists+=Slice->FillSinkLists.Elapsed();
                CountSinks+=Slice->CountSinks.Elapsed();
                CalcSinkBlocks+=Slice->CalcSinkBlocks.Elapsed();
                AllocAccels += Slice->AllocAccels.Elapsed();
            FillSourceLists+=Slice->FillSourceLists.Elapsed();
                CountSources+=Slice->CountSources.Elapsed();
                CalcSourceBlocks+=Slice->CalcSourceBlocks.Elapsed();
            FillInteractionList+=Slice->FillInteractionList.Elapsed();
        DeviceThreadTimer += Slice->DeviceThreadTimer.Elapsed();
            LaunchDeviceKernels += Slice->LaunchDeviceKernels.Elapsed();
            FillSinks += Slice->FillSinks.Elapsed();
            FillSources += Slice->FillSources.Elapsed();
            WaitForResult += Slice->WaitForResult.Elapsed();
            CopyAccelFromPinned += Slice->CopyAccelFromPinned.Elapsed();

        int g = Slice->AssignedDevice;
        DirectInteractions_GPU[g] += Slice->DirectTotal;
        TotalDirectInteractions_GPU += Slice->DirectTotal;

        GB_to_device[g] += Slice->bytes_to_device/1e9;
        GB_from_device[g] += Slice->bytes_from_device/1e9;
        DeviceSinks[g] += Slice->SinkTotal;
        DeviceSources[g] += Slice->SourceTotal;
    }
    FinalizeBookkeeping.Stop();
    
    CopyPencilToSlab.Start();
    #pragma omp parallel for schedule(static)
    for(int sliceIdx = 0; sliceIdx < NSplit*WIDTH; sliceIdx++){
        int w = sliceIdx / NSplit;
        SetInteractionCollection *Slice = Slices[sliceIdx];
        int Nj = (cpd - w)/width;
        if(Nj * width + w < cpd)
            Nj++;

        for(int sinkIdx = 0; sinkIdx < Slice->NSinkList; sinkIdx++){
            int SinkCount = Slice->SinkSetCount[sinkIdx];
            int j = sinkIdx%Nj;
            int k = sinkIdx/Nj + Slice->K_low;
            int jj = w + j * width;
            int zmid = PP->WrapSlab(jj + P.NearFieldRadius);
            int Start = Slice->SinkSetStart[sinkIdx];

            for(int z=zmid-nfr; z <= zmid+nfr;z++){
                int cid = PP->CellID(k,z);
                int CellNP = PP->NumberParticle(slab,k,z);
                accstruct *a = PP->NearAccCell(slab,k,z);

                if(P.ForceOutputDebug == 1){
                    for(int p =0; p < CellNP; p++){
                        assert(isfinite(Slice->SinkSetAccelerations[Start +p].x));
                        assert(isfinite(Slice->SinkSetAccelerations[Start +p].y));
                        assert(isfinite(Slice->SinkSetAccelerations[Start +p].z));
                    }
                }

                pthread_mutex_lock(&CellLock[cid]);
                // Should rewrite the following to avoid so much repetition
                if(CellAccInit[cid]){
                    #pragma simd assert
                    for(int p = 0; p <CellNP; p++){
                        // It's not good to know about spline here, but it's likely the most efficient way
                        #ifdef DIRECTSINGLESPLINE
                        a[p] += Slice->SinkSetAccelerations[Start +p]*inv_eps3;
                        #else
                        a[p] += Slice->SinkSetAccelerations[Start +p];
                        #endif
                    }
                } else {
                    #pragma simd assert
                    for(int p = 0; p <CellNP; p++){
                        // It's not good to know about spline here, but it's likely the most efficient way
                        #ifdef DIRECTSINGLESPLINE
                        a[p] = Slice->SinkSetAccelerations[Start +p]*inv_eps3;
                        #else
                        a[p] = Slice->SinkSetAccelerations[Start +p];
                        #endif
                    }
                    CellAccInit[cid] = 1;
                }
                pthread_mutex_unlock(&CellLock[cid]);

                Start+=CellNP;
                SinkCount-=CellNP;
            }
            assert(SinkCount == 0);
        }
        delete Slice;
    }
    CopyPencilToSlab.Stop();
    
    delete[] Slices;
    delete[] CellAccInit;
    if(P.ForceOutputDebug) CheckGPUCPU(slab);
    #endif
    FinalizeTimer.Stop();
}

#include <vector>
#include <algorithm>
void NearFieldDriver::CheckInteractionList(int slab){

    vector<uint64> ** il = new std::vector<uint64> *[P.cpd*P.cpd];

    for(int i = 0; i < P.cpd*P.cpd; i++){
        il[i] = new vector<uint64>();
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
    
    vector<uint64> ** il_test = new std::vector<uint64> *[P.cpd*P.cpd];
    for(int i = 0; i < P.cpd*P.cpd; i++){
        il_test[i] = new vector<uint64>();
        il_test[i]->reserve(WIDTH*WIDTH*WIDTH);
    }
    for(int y = 0; y < P.cpd; y++){
        for(int z = 0; z < P.cpd; z++){
            for(int i = slab - RADIUS; i <= slab + RADIUS; i++){
                for(int j = y - RADIUS; j <= y + RADIUS; j++){
                    for(int k = z - RADIUS; k <= z + RADIUS; k++){
                        uint64 sinkCellId = y*P.cpd + z;
                        uint64 sourceCellId = P.cpd*P.cpd * PP->WrapSlab(i) +
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
    
    STDLOG(1,"Checking the interaction list for slab %d passed.\n", slab);
    //Checking the interaction list is time-consuming
    //it should never be done in an actual run
}
    
NearFieldDriver *JJ;
    
#include "PencilOnPencil.cc"
#include "PencilPlan.cc"
