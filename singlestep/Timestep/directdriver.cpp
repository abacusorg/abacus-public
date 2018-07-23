/* directdriver.cpp
 *
 * Abstraction layer that handles everything associated with launching the directs. 
 * All details of direct implementation should be invisible to classes above this one.
 * 
 * FIXME: DWF reccomends refactoring this class to reduce its direct responsibilities. 
 * It has turned into a bit of a God class---NearFieldDriver is currently directly responsible for:
 *  1. Initializing the CPU directs
 *  2. Initializing the GPU
 *  3. Directly executing the procedures for CPU and GPU calculation
 *  4. Calculating and reporting timings for all of the above.
 *  5. Directly cleaning up many parts of the above.
 *
 *  This is basically untestable, and leads to lots of interdependent but widely seperated code.
 *
 * This should be refactored into something like the following hierarchy:
 * class NearFieldDriver:
 *  Top level abstraction layer that sends commands to a specified calculation module. 
 *  Calculation modules must implement 
 *  interface NearFieldCalculator:
 *      Exposes a common set of controls for different types of direct module.
 *      Each implementation should be solely responsible for its own initialization/tear-down
 *      Implemented by:
 *          Class CPUDirectCalculator:
 *              Performs the direct calculation on the CPU. Very simple
 *          Class AVXDirectCalculator:
 *              Performs the direct calculation via AVX on the CPU. 
 *          Class GPUDirectCalculator:
 *              Abstracts away all of the GPU initialization and calculation complexities.
 *              Should probably only be an abstraction layer on top of several classes that
 *              hand GPU memory management, computation, and synchronization
 *      Each implementation that uses PencilOnPencil should take a DirectPencilFactory * in its
 *      constructor that abstracts that away from the rest of the direct computation. The 
 *      DirectPencilFactory should take a slab and its cellinfo, and should return a newly
 *      contructed SetInteractionCollection
 *      Should contain methods like the following:
 *          void ExecuteSlab(..., DirectTimings * timer):
 *              Orders the Direct execution for the specified slab. Uses timer for timings.
 *              Pencil creation is abstracted to the DirectPencilFactory
 *          bool IsNearFieldDoneForSlab(..., DirectTimings *timer):
 *              Queries completion of the slab. Replaces SlabDone
 *          void CleanUpSlab(..., DirectTimings *timer):
 *              Replaces Finalize. Collects timings and frees all memory associated
 *              with calculating the directs for the slab.
 *
 *  This will be much simpler to understand, easier to make mocks for, and much easier to test. 
 *  The common code between directs and microstepping will also be much easier to extract into 
 *  seperate classes. 
 */
#include "StructureOfLists.cc"
#include "SetInteractionCollection.hh"

#include "DeviceFunctions.h"

// Foward declaration; found in kick.cpp
void ZeroAcceleration(int slab,int Slabtype);

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
    
#ifdef CUDADIRECT
    NGPU = GetNGPU();
    GPUMemoryGB = GetDeviceMemory();
    GPUMemoryGB = std::min(5.0e-9*P.np*sizeof(FLOAT3),GPUMemoryGB/(DirectBPD));
    GPUMemoryGB = std::min(GPUMemoryGB,.05/(DirectBPD*NGPU)*P.MAXRAMMB/1024);  // Don't pin more than 5% of host memory
    STDLOG(1, "Using %f GB of GPU memory (per GPU thread)\n", GPUMemoryGB);
    size_t BlockSizeBytes = sizeof(FLOAT) *3 * NFBlockSize;
    MaxSinkBlocks = floor(1e9 * GPUMemoryGB/(6*BlockSizeBytes));
    MaxSourceBlocks = WIDTH * MaxSinkBlocks;
    STDLOG(1,"Initializing GPU with %f x10^6 sink blocks and %f x10^6 source blocks\n",
            MaxSinkBlocks/1e6,MaxSourceBlocks/1e6);
    
    // We want to compute maxkwidth for GPUSetup
    // maxkwidth will occur when we have the fewest splits
    // The fewest splits will occur when we are operating on the smallest slabs
    MinSplits = GetNSplit(WIDTH*Slab->min, Slab->min);
    //TODO: this fudge factor fails sometimes (e.g. Ewald test for CPD 131).  What's the right way to compute this?
    MinSplits = (int) max(1.,floor(.8*MinSplits));  // fudge factor to account for uneven slabs
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
            // We may wish to make these in an order that will alter between GPUs
            int kh = SplitPoint[n];
            // The construction and execution are timed internally in each SIC then reduced in Finalize(slab)
            SlabInteractionCollections[slabID][w*NSplit + n] = 
                new SetInteractionCollection(slabID,w,kl,kh,P.DensityKernelRad2);
            
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
    accstruct * a_gpu = (accstruct *) LBW->ReturnIDPtr(AccSlab,slabID);
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
        acc3struct ai_g = a_gpu[i];
        acc3struct ai_c = a_cpu[i];
        if(ai_g.norm() == 0. && ai_c.norm() == 0.)
            continue;
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
    if(!LBW->IDPresent(AccSlab, slabID))
        LBW->AllocateArena(AccSlab,slabID);
    ZeroAcceleration(slabID,AccSlab);
    
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
			// TODO: At present, the b2 parameter is not passed
			// into the CPU directs, so there is no FOF neighbor
			// computation.
			// TODO: sink_acc may now be a float4, but the CPU routines
			// want a float3.  We'll overload this space and fix it later
                        if(np_source >0) DD[g].AVXExecute(sink_pos,source_pos,np_sink,np_source,
                                delta,eps,(FLOAT3 *)sink_acc);
                        delete[] source_pos;
                        
                        DI_slab += np_sink*np_source;
                    }
                }
            }
	    // All done with this cell.  Fix the float4 to float3 issue
	    acc3struct *acc3 = (acc3struct *)sink_acc;
	    for (uint64 i=np_sink-1; i>=0; i--) 
	    	sink_acc[i] = accstruct(acc3[i]);
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

    if(P.ForceOutputDebug)
        CheckInteractionList(slab);
    
    LBW->AllocateArena(AccSlab,slab);

    SetInteractionCollection ** Slices = SlabInteractionCollections[slab];
    int NSplit = SlabNSplit[slab];

    int cpd = P.cpd;
    int nfr = P.NearFieldRadius;
    int width = 2*nfr +1;

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

    if(P.ForceOutputDebug){
        // Make the accelerations invalid so we can detect improper co-adding
        #pragma omp parallel for schedule(static)
        for(int y = 0; y < cpd; y++){
            for(int z = 0; z < cpd; z++){
                accstruct *acc = PP->NearAccCell(slab, y, z);
                int count = PP->NumberParticle(slab,y,z);
                
                for(int i = 0; i < count; i++)
                    acc[i] = accstruct(std::numeric_limits<float>::infinity());
                    //acc[i] = accstruct(0.);
            }
        }
    }
    FinalizeBookkeeping.Stop();
    
    CopyPencilToSlab.Start();    

    // Compute the Number of j's in each registration
    // This is the same for all y rows (k's).
    int Nj[WIDTH];
    for (int w=0; w<WIDTH; w++) {
        Nj[w] = (cpd - w)/WIDTH;
        if(Nj[w] * WIDTH + w < cpd)
            Nj[w]++;
    }
    // TODO: It's weird that Nj requires computation instead of being
    // stored in the SetInteractionCollection class.  Indeed, shouldn't
    // the class have methods to compute the index number from (y,z)
    // and vice versa?  Too much re-coded math.

    #pragma omp parallel for schedule(static)
    for (int k=0; k<cpd; k++) {
        // We're going to find all the Sets associated with this row.
        // There are WIDTH of them, one for each z registration
        SetInteractionCollection *theseSlices[WIDTH];
        for (int w=0; w<WIDTH; w++)
            theseSlices[w] = NULL;
        for (int sliceIdx=0; sliceIdx< NSplit*WIDTH; sliceIdx++) {
        // Just doing an exhaustive search, so we assume less about
        // the upstream data model.
            if (k >= Slices[sliceIdx]->K_low && k < Slices[sliceIdx]->K_high) {
               // We've found a Slice.
               // Compute which z-registration this is.
               int w = sliceIdx/NSplit; 
               // TODO: Again, why is this not a part of the SetInteractionCollection?
               theseSlices[w] = Slices[sliceIdx];
            }
        }
        for (int w=0; w<width; w++)
            assertf(theseSlices[w] != NULL, "We failed to find all z-registrations");

        // We're going to set the acceleration the first time we encounter
        // a cell.  The cells are encountered in z order, starting from 
        // z=0 and increasing.  That means that when we get to z=cpd,
        // we've already been there once.
        int firsttouch = 0;
        
        // Now we're going to proceed through the pencils
        for (int jj=0; jj<cpd; jj++) {
            int w = jj%width;   // This is which registration we're in
            int j = (jj-w)/width;
            SetInteractionCollection *Slice = theseSlices[w];
            int sinkIdx = (k - Slice->K_low)*Nj[w] + j;

            // And now we can continue with the previous stuff
            int SinkCount = Slice->SinkSetCount[sinkIdx];
            int zmid = PP->WrapSlab(jj + P.NearFieldRadius);
            int Start = Slice->SinkSetStart[sinkIdx];

            for(int z=zmid-nfr; z <= zmid+nfr; z++){
                int CellNP = PP->NumberParticle(slab,k,z);
                accstruct *a = PP->NearAccCell(slab,k,z);

                if(P.ForceOutputDebug == 1){
                    for(int p =0; p < CellNP; p++){
                        assert(isfinite(Slice->SinkSetAccelerations[Start +p].x));
                        assert(isfinite(Slice->SinkSetAccelerations[Start +p].y));
                        assert(isfinite(Slice->SinkSetAccelerations[Start +p].z));
                    }

                    if(z==firsttouch && z<cpd){
                        for(int p = 0; p <CellNP; p++){
                            assert(!isfinite(a[p].x));
                            assert(!isfinite(a[p].y));
                            assert(!isfinite(a[p].z));
                        }
                    } else {
                        for(int p = 0; p <CellNP; p++){
                            assert(isfinite(a[p].x));
                            assert(isfinite(a[p].y));
                            assert(isfinite(a[p].z));
                        }
                    }
                }

                // Set when this is the first touch; accumulate otherwise
                if(z==firsttouch && z<cpd){
                    firsttouch++;
                    #pragma simd assert
                    for(int p = 0; p <CellNP; p++)
                        a[p] = Slice->SinkSetAccelerations[Start +p];
                } else {
                    #pragma simd assert
                    for(int p = 0; p <CellNP; p++)
                        a[p] += Slice->SinkSetAccelerations[Start +p];
                }

                Start+=CellNP;
                SinkCount-=CellNP;
            }
            assert(SinkCount == 0);
        }
    }

    // ********************************************************************* old code //
/*
    // The co-adding in here doesn't really benefit from double precision, as far as Spiral and Ewald indicate
    for(int sliceIdx = 0; sliceIdx < NSplit*WIDTH; sliceIdx++){
        int w = sliceIdx / NSplit;
        SetInteractionCollection *Slice = Slices[sliceIdx];
        int Nj = (cpd - w)/width;
        if(Nj * width + w < cpd)
            Nj++;

        assert(Slice->NSinkList % Nj == 0);
        // Thread over y-rows.  This keeps threads spread out prevents races at the z-wrap
        #pragma omp parallel for schedule(static)
        for(int k = Slice->K_low; k < Slice->NSinkList / Nj + Slice->K_low; k++){  // y-rows
            for(int j = 0; j < Nj; j++){ // z-pencil centers
                int jj = w + j*width;
                int zmid = jj + P.NearFieldRadius;  // this is deliberately not wrapped so we can tell when we cross the wrap

                int sinkIdx = (k - Slice->K_low)*Nj + j;
                int SinkCount = Slice->SinkSetCount[sinkIdx];
                int Start = Slice->SinkSetStart[sinkIdx];

                for(int z=zmid-nfr; z <= zmid+nfr; z++){
                    int CellNP = PP->NumberParticle(slab,k,z);
                    accstruct *a = PP->NearAccCell(slab,k,z);

                    if(P.ForceOutputDebug == 1){
                        for(int p =0; p < CellNP; p++){
                            assert(isfinite(Slice->SinkSetAccelerations[Start +p].x));
                            assert(isfinite(Slice->SinkSetAccelerations[Start +p].y));
                            assert(isfinite(Slice->SinkSetAccelerations[Start +p].z));
                        }

                        if(w != 0 || z >= cpd){
                            for(int p = 0; p <CellNP; p++){
                                assert(isfinite(a[p].x));
                                assert(isfinite(a[p].y));
                                assert(isfinite(a[p].z));
                            }
                        } else {
                            for(int p = 0; p <CellNP; p++){
                                assert(!isfinite(a[p].x));
                                assert(!isfinite(a[p].y));
                                assert(!isfinite(a[p].z));
                            }
                        }
                    }

                    // The first time we see a cell is always when w==0, and we have not yet z-wrapped
                    if(w != 0 || z >= cpd){
                        #pragma simd assert
                        for(int p = 0; p <CellNP; p++)
                            a[p] += Slice->SinkSetAccelerations[Start +p];
                    } else {
                        #pragma simd assert
                        for(int p = 0; p <CellNP; p++)
                            a[p] = Slice->SinkSetAccelerations[Start +p];
                    }

                    Start+=CellNP;
                    SinkCount-=CellNP;
                }
                assert(SinkCount == 0);
            }
        }
        delete Slice;
    }*/

    // Do a final pass to delete all slices
    // TODO: can this be done inline above?
    for(int sliceIdx = 0; sliceIdx < NSplit*WIDTH; sliceIdx++){
        SetInteractionCollection *Slice = Slices[sliceIdx];
        delete Slice;
    }

    CopyPencilToSlab.Stop();
    
    delete[] Slices;
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
            printf("Error: cell %d has %lld =/= %lld interactions\n",i,il[i]->size(),WIDTH*WIDTH*WIDTH);
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
