/* directdriver.cpp
 *
 * Abstraction layer that handles everything associated with launching the directs. 
 * All details of direct implementation should be invisible to classes above this one.
 * 
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

        // Total timings from running the SetInteractionCollections in the host threads that communicate with the GPUs
        // These are gathered in Finalize(slab)
        double DeviceThreadTimer = 0;
        double     LaunchDeviceKernels = 0;
        double     FillSinks = 0;
        double     FillSources = 0;
        double     WaitForResult = 0;
        double     CopyAccelFromPinned = 0;

	int MaxNSplits = 0;
    
        double *GB_to_device, *GB_from_device;
        uint64 *DeviceSinks;
        uint64 *DeviceSources;
        
        STimer CalcSplitDirects;
        STimer SICConstruction;
        STimer SICExecute;
        STimer CPUFallbackTimer;
        STimer FinalizeTimer;

        void AggregateStats();  // Called before shutdown
        // The following totals are filled in by AggregateStats()
        double GPUThroughputTime = 0;
        double total_GB_to = 0, total_GB_from = 0, total_sinks = 0, total_sources = 0, gdi_gpu = 0;
        double mean_splits_per_slab = 0;

    private:
        int *slabcomplete;
        int * SlabNSplit;

        SetInteractionCollection *** SlabInteractionCollections; 
        int WIDTH;
        int RADIUS;
        Direct *DD;
        int NGPU;
	int NBuffers;
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
};



// This constructor calls GPUSetup() to start the GPU threads.
// But first, it has to set up some configurations.
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
    GPUMemoryGB /= DirectBPD;	// Nominal GB per buffer
    NBuffers = NGPU*DirectBPD;
    STDLOG(1, "Running with %d GPUs, each with %d Buffers of max size %f GB\n", NGPU, DirectBPD, GPUMemoryGB);

    // No need to go crazy if the problem is small.  But small 
    // problems can be highly clustered.
    // Even if CPD=WIDTH, we'd be loading all of the positions as sources and 
    // 1/WIDTH as sinks WIDTH-times over.  That's 2*WIDTH FLOAT3 + WIDTH FLOAT4,
    // divided over NBuffers.  
    // Round up by a factor of 1.3 and an extra 4 MB, just to be generous
    GPUMemoryGB = std::min(GPUMemoryGB, P.np*1e-9*sizeof(accstruct)*3*WIDTH/NBuffers*1.3+0.004);

    // GPUMemoryGB = std::min(GPUMemoryGB, 5.0*P.np*1e-9*sizeof(FLOAT3)+0.004);

    // Don't pin more than 10% of the host memory.
    GPUMemoryGB = std::min(GPUMemoryGB,0.10/(NBuffers)*P.MAXRAMMB/1024);  

    STDLOG(1, "Using %f GB of GPU memory (per GPU thread)\n", GPUMemoryGB);
    MinSplits = NBuffers;
    MaxNSplits = MinSplits;

    // Put a floor to insist on using all GPUs
    STDLOG(1,"MinSplits = %d\n", MinSplits);

    GPUSetup(P.cpd, 1.0e9*GPUMemoryGB, NGPU, DirectBPD, P.GPUThreadCoreStart, P.NGPUThreadCores, &MaxSinkBlocks, &MaxSourceBlocks);
    STDLOG(1,"Initializing GPU with %7.3f x10^3 sink blocks and %7.3f x10^3 source blocks\n",
            MaxSinkBlocks/1e3,MaxSourceBlocks/1e3);
    // This returns the number of SinkBlocks and number of SourceBlocks
    // that any future SIC is allowed.  Overheads for other bookkeeping
    // has been included in this estimate, so we just need to obey this.


    SICExecute.Clear();

    NSink_GPU_final = (uint64*) malloc(NGPU*sizeof(uint64));
    for(int g = 0; g < NGPU; g++)
        NSink_GPU_final[g] = 0;
    
    DirectInteractions_GPU = new uint64[NBuffers]();
    GB_to_device = new double[NBuffers]();
    GB_from_device = new double[NBuffers]();
    DeviceSinks = new uint64[NBuffers]();
    DeviceSources = new uint64[NBuffers]();
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

#include "extras_directdriver.cpp"

/** Compute the number of splits for a given number of sources and sinks particles

Now that a SIC includes all of the cells in a Y skewer, I think
that we can make a better estimate.  The skewer will have each
sink particle 5 times.

*/

// ================ Code to queue up work for the GPU ==================

uint64 ComputeSICSize(int cpd, int np, int WIDTH, int NSplit);
int ComputeSinkBlocks(int slab, int j);
int ComputeSourceBlocks(int slab, int j);

void NearFieldDriver::ExecuteSlabGPU(int slabID, int blocking){
    // This will construct the relevant SetInteractionCollections,
    // then submit them to the queue.

    // First, calculate the required subdivisions of the slab to fit in GPU memory
    CalcSplitDirects.Start();

    int *SinkBlocks = new int[P.cpd];
    int *SourceBlocks = new int[P.cpd];
    int totSinkBlocks = 0;
    int totSourceBlocks = 0;

    for (int j=0;j<P.cpd;j++) {
	SinkBlocks[j] = ComputeSinkBlocks(slabID, j);
	totSinkBlocks += SinkBlocks[j];
	SourceBlocks[j] = ComputeSourceBlocks(slabID, j);
	totSourceBlocks += SourceBlocks[j];
	// These contain the number of blocks in each j.
    }
    // We will need at least totBlocks/maxBlocks, rounded up, splits.
    // Given the maximum number of blocks in a SIC, it will be 
    // easy to pack these in.  But the problem is that we want to 
    // use all of the GPU buffers for efficiency, so that sets a 
    // minimum.
    STDLOG(1,"Found %d sink and %d source blocks in slab %d\n", totSinkBlocks, totSourceBlocks, slabID);
    int useMaxSink = totSinkBlocks/NBuffers;
    useMaxSink *= (1+1.1*NBuffers/P.cpd); 
	// Trying to bias up to leave the last one short instead of long
    useMaxSink = std::min(useMaxSink, MaxSinkBlocks);
    int useMaxSource = MaxSourceBlocks;
    STDLOG(2,"Using max sink %d and source %d.\n", useMaxSink, useMaxSource);

    // Now we want to pack these in
    int SplitPoint[1024]; // Just a crazy number
    int NSplit = 0;
    int thisSink = 0, thisSource = 0;

    for (int j=0; j<P.cpd; j++) {
	thisSink += SinkBlocks[j];
	thisSource += SourceBlocks[PP->WrapSlab(j+RADIUS)];
        if (j==0 || thisSink>useMaxSink ||
            thisSource>useMaxSource) {
	    // We've overflowed the previous set; end it 
	    if (NSplit>0) {
		SplitPoint[NSplit-1] = j;
		STDLOG(1,"Split %d ending at y=%d; %d sink and %d source blocks\n", 
		    NSplit-1, j, thisSink-SinkBlocks[j], thisSource-SourceBlocks[PP->WrapSlab(j+RADIUS)]);
	    }
	    // Time to start a new one
	    thisSink = SinkBlocks[j];
	    thisSource = 0;
	    for (int c=-RADIUS; c<=RADIUS; c++) 
		thisSource += SourceBlocks[PP->WrapSlab(j-c)];
	    NSplit++;
	}
    }
    // And always end the last one
    SplitPoint[NSplit-1] = P.cpd;
    STDLOG(1,"Split %d ending at y=%d; %d sink and %d source blocks\n", 
		    NSplit-1, P.cpd, thisSink, thisSource);
    // There is a failure mode here if the last skewer by itself
    // overflows; it will end up assigned without prior question.
    assertf(thisSink<MaxSinkBlocks && thisSource<MaxSourceBlocks,
    	"Sinks or Sources of the last skewer overflow the maxima.");

    uint64 NSink = Slab->size(slabID);
    STDLOG(1,"Using %d direct splits on slab %d, max blocks %d sink and %d source\n", 
    	NSplit, slabID, useMaxSink, useMaxSource);

    delete[] SinkBlocks;
    delete[] SourceBlocks;
	
    SlabNSplit[slabID] = NSplit;
    MaxNSplits = std::max(MaxNSplits, NSplit);
    SlabInteractionCollections[slabID] = new SetInteractionCollection *[NSplit];

    // Now we need to make the NearField_SIC_Slab
    uint64 bsize = ComputeSICSize(P.cpd, NSink, WIDTH, NSplit);
    LBW->AllocateSpecificSize(NearField_SIC_Slab, slabID, bsize);
    char *buffer = LBW->ReturnIDPtr(NearField_SIC_Slab, slabID);

    // And allocate space for the accelerations
    LBW->AllocateArena(AccSlab,slabID);

    CalcSplitDirects.Stop();

    // Initialize the space -- This should no longer be needed
    /*
    FLOAT *p = (FLOAT *)LBW->ReturnIDPtr(AccSlab, slabID);
    for (int j=0; j<LBW->IDSizeBytes(AccSlab,slabID)/sizeof(accstruct)*4; j++) p[j] = 0.0;
    */


    if(P.ForceOutputDebug){
        // Make the accelerations invalid so we can detect improper co-adding
        #pragma omp parallel for schedule(static)
        for(int y = 0; y < P.cpd; y++){
            for(int z = 0; z < P.cpd; z++){
                accstruct *acc = PP->NearAccCell(slabID, y, z);
                int count = PP->NumberParticle(slabID,y,z);
                
                for(int i = 0; i < count; i++)
                    acc[i] = accstruct(std::numeric_limits<float>::infinity());
            }
        }
    }

    int jl =0;
    for(int n = 0; n < NSplit; n++){
        SICConstruction.Start();
	// We may wish to make these in an order that will alter between GPUs
	int jh = SplitPoint[n];
	// The construction and execution are timed internally in each SIC then reduced in Finalize(slab)
	// This is where the SetInteractionCollection is actually constructed
	// STDLOG(1,"Entering SIC Construction with %d bytes\n", bsize);
	SlabInteractionCollections[slabID][n] = 
	    new SetInteractionCollection(slabID,jl,jh,WriteState.DensityKernelRad2, buffer, bsize);
	// STDLOG(1,"Leaving SIC Construction with %d bytes\n", bsize);
	
	// Check that we have enough blocks
	int NSinkBlocks = SlabInteractionCollections[slabID][n]->NSinkBlocks;
	int NSourceBlocks = SlabInteractionCollections[slabID][n]->NSourceBlocks;
	assertf(NSinkBlocks <= MaxSinkBlocks,
		"NSinkBlocks (%d) is larger than MaxSinkBlocks (%d)\n", NSinkBlocks, MaxSinkBlocks);
	assertf(NSourceBlocks <= MaxSourceBlocks,
		"NSourceBlocks (%d) is larger than MaxSourceBlocks (%d)\n", NSourceBlocks, MaxSourceBlocks);

    SICConstruction.Stop();
    SICExecute.Start();
	
	STDLOG(1,"Executing directs for slab %d, y: %d - %d\n",slabID,jl,jh);
	// This SIC is ready; send it to be executed
	SlabInteractionCollections[slabID][n]->GPUExecute(blocking);
	//SlabInteractionCollections[slabID][n]->CPUExecute();
	jl = jh;
    SICExecute.Stop();
    }

    STDLOG(1, "%l bytes remaining after SIC allocation on slab %d (%4.1f%% unused)\n", 
	bsize, slabID, 100.0*bsize/LBW->IDSizeBytes(NearField_SIC_Slab, slabID));
    return;
}



void NearFieldDriver::ExecuteSlab(int slabID, int blocking){
    // This routine is the usual entry point when wanting to queue up
    // a slab for computation.
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


// ================ Code to Handle the return of information ===========

int NearFieldDriver::SlabDone(int slab){
    // Return 1 if all SetInteractionCollections have been completed.
    slab = PP->WrapSlab(slab);

    if (slabcomplete[slab] == 0){

        int NSplit = SlabNSplit[slab];

        int complete = 1;
        for(int i = 0; i < NSplit; i++)
            complete = complete && SlabInteractionCollections[slab][i]->CheckCompletion();

        slabcomplete[slab] = complete;
        return complete;
    }

    return 1;
}


void NearFieldDriver::Finalize(int slab){
#ifndef CUDADIRECT
    return;
#endif

    FinalizeTimer.Start();

    // When all of the SIC are finished, call this routine.
    // It will accumulate the timings and statistics, 
    // and then delete the SIC.
    slab = PP->WrapSlab(slab);

    assertf(SlabDone(slab) != 0,
            "Finalize called for slab %d but it is not complete\n",slab);

    if(P.ForceOutputDebug)
        CheckInteractionList(slab);

    SetInteractionCollection **Slices = SlabInteractionCollections[slab];
    int NSplit = SlabNSplit[slab];

    int cpd = P.cpd;
    int nfr = RADIUS;

    // Collect the statistics and timings
    for(int sliceIdx = 0; sliceIdx < NSplit; sliceIdx++){
        SetInteractionCollection *Slice = Slices[sliceIdx];
        
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

    // Just test that AccSlab is not crazy
    /*
    FLOAT *p = (FLOAT *)LBW->ReturnIDPtr(AccSlab, slab);
    for (int j=0; j<LBW->IDSizeBytes(AccSlab,slab)/sizeof(accstruct)*4; j++) 
    	assertf(isfinite(p[j]) && abs(p[j])<10, "Accelerations appear crazy\n");
    */

    // Do a final pass to delete all slices
    for(int sliceIdx = 0; sliceIdx < NSplit; sliceIdx++){
        SetInteractionCollection *Slice = Slices[sliceIdx];
        delete Slice;
    }
    LBW->DeAllocate(NearField_SIC_Slab, slab);
    delete[] Slices;
    
    if(P.ForceOutputDebug) CheckGPUCPU(slab);
    FinalizeTimer.Stop();
}

/** We want to report some timings and diagnostic information once the sim has completed.
 *  This function computes some useful sums that will be written to the timing file or global
 *  log, so it should be called after timestep() but before the logs are written.
 */
void NearFieldDriver::AggregateStats(){
#ifdef CUDADIRECT
    for(int g = 0; g < NBuffers; g++){
        total_GB_to += GB_to_device[g];
        total_GB_from += GB_from_device[g];
        total_sinks += DeviceSinks[g];
        total_sources += DeviceSources[g];
    }

    // Measures the total amount of time we have at least one GPU thread running
    assertf(!SetInteractionCollection::GPUThroughputTimer.timeron, "GPU throughput timer still on! Atomic thread counting failure?\n");
    GPUThroughputTime = SetInteractionCollection::GPUThroughputTimer.Elapsed();

    gdi_gpu = TotalDirectInteractions_GPU/1e9;

    for(int s = 0; s < P.cpd; s++)
        mean_splits_per_slab += SlabNSplit[s];
    mean_splits_per_slab /= P.cpd;

    STDLOG(0, "Maximum NSplits used in Directs: %d\n", MaxNSplits);
    STDLOG(0, "Mean NSplits used in Directs: %d\n", mean_splits_per_slab);

#endif
}


NearFieldDriver *JJ;
    
#include "PencilOnPencil.cc"
#include "PencilPlan.cc"

