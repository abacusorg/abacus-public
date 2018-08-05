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
        PTimer CopyPencilToSlabSetup;
        PTimer CopyPencilToSlabCopy;
        STimer TearDownPencils;
        STimer ZeroAccel;

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
        int GetNSplit(uint64 NSource, uint64 NSink);
        
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
    STDLOG(1, "Running with %d GPUs, each with %d Buffers of max size %f GB\n", NGPU, DirectBPD, GPUMemoryGB);

    // No need to go crazy if the problem is small.
    // Even if CPD=WIDTH, we'd be loading all of the positions as sources and 
    // 1/WIDTH as sinks.  So if we have 5x the total positions, we can't lose.
    GPUMemoryGB = std::min(GPUMemoryGB, 50.0*P.np*1e-9*sizeof(FLOAT3));

    // Don't pin more than 10% of the host memory.
    GPUMemoryGB = std::min(GPUMemoryGB,0.10/(DirectBPD*NGPU)*P.MAXRAMMB/1024);  

    STDLOG(1, "Using %f GB of GPU memory (per GPU thread)\n", GPUMemoryGB);

    // The fewest splits will occur when we are operating on the smallest slabs
    MinSplits = GetNSplit(WIDTH*Slab->min, Slab->min);
    MinSplits = std::max(NGPU,(int)floor(.8*MinSplits));  
    	// fudge factor to account for uneven slabs
    // Put a floor to insist on using all GPUs
    STDLOG(1,"MinSplits = %d\n", MinSplits);

    NBuffers = NGPU*DirectBPD;
    GPUSetup(P.cpd, 1.0e9*GPUMemoryGB, NGPU, DirectBPD, P.GPUThreadCoreStart, P.NGPUThreadCores, &MaxSinkBlocks, &MaxSourceBlocks);
    STDLOG(1,"Initializing GPU with %f x10^6 sink blocks and %f x10^6 source blocks\n",
            MaxSinkBlocks/1e6,MaxSourceBlocks/1e6);
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

// TODO: Can remove this function?
int NearFieldDriver::GetNSplit(uint64 NSource, uint64 NSink){
    // Pencils are aligned to NFBlocksize boundaries
    // Best case: we will have ~NSource/NFBlocksize blocks
    // Worst case: We have all but one pencil NFBlocksize + 1 long 
    //            with last pencil holding the rest
    int N1 = 2*( P.cpd*(P.cpd + WIDTH) -1 );
    int N2 =std::ceil( (NSource -( NFBlockSize + 1.0) * ( P.cpd*(P.cpd + WIDTH) -1 ))/(1.0*NFBlockSize) + 1);
    if (N2 < 0) N2 = 0; //if we don't have at least GPUBlocksize particles per pencil, this calculation fails
    
    int SourceBlocks = N1 + N2;
    
    N1 = 2*(P.cpd*P.cpd -1);
    N2 =std::ceil( (NSink -( NFBlockSize + 1.0) * ( P.cpd*(P.cpd) -1 ))/(1.0*NFBlockSize) + 1);
    if (N2 < 0) N2 = 0;
    int SinkBlocks = N1 + N2;
    
    int NSplitSource = std::ceil((1.0*SourceBlocks)/MaxSourceBlocks);
    int NSplitSink = std::ceil((1.0*SinkBlocks)/MaxSinkBlocks);
    const int NSplit = std::max(NSplitSource,NSplitSink);
    
    return NSplit;
}


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
    STDLOG(1,"Using max sink %d and source %d.\n", useMaxSink, useMaxSource);

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
    SlabInteractionCollections[slabID] = new SetInteractionCollection *[NSplit];


#ifdef OLDCODE

    //Add up the sources
    uint64 NSource = 0;
    for(int i = slabID-RADIUS; i <= slabID + RADIUS; i++) NSource+=Slab->size(i);

    const int NSplit = std::max(MinSplits,GetNSplit(NSource, NSink));
    assertf(NSplit >= MinSplits, "NSplit (%d) is less than MinSplits (%d)\n", NSplit, MinSplits);
    STDLOG(1,"Splitting slab %d into %d blocks for directs.\n",slabID,NSplit);
    SlabNSplit[slabID] = NSplit;
    
    uint64 NPTarget = NSink/NSplit;
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
#endif

    CalcSplitDirects.Stop();

    // Now we need to make the NearField_SIC_Slab
    uint64 bsize = ComputeSICSize(P.cpd, NSink, WIDTH, NSplit);
    LBW->AllocateSpecificSize(NearField_SIC_Slab, slabID, bsize);
    char *buffer = LBW->ReturnIDPtr(NearField_SIC_Slab, slabID);

    // And allocate space for the accelerations
    LBW->AllocateArena(AccSlab,slabID);

    // Initialize the space
    FLOAT *p = (FLOAT *)LBW->ReturnIDPtr(AccSlab, slabID);
    // TODO: May be able to get rid of this initialization
    for (int j=0; j<LBW->IDSizeBytes(AccSlab,slabID)/sizeof(accstruct)*4; j++) p[j] = 0.0;

    // If we thread over y-splits, would that help with NUMA locality?
    //? for(int k_mod = 0; k_mod < WIDTH; k_mod++){
        int jl =0;
        for(int n = 0; n < NSplit; n++){
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
            
            STDLOG(1,"Executing directs for slab %d, y: %d - %d\n",slabID,jl,jh);
            SICExecute.Start();
	    // This SIC is ready; send it to be executed
            SlabInteractionCollections[slabID][n]->GPUExecute(blocking);
            //SlabInteractionCollections[slabID][n]->CPUExecute();
            SICExecute.Stop();
            jl = jh;
        }
    //}
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



class CoaddPlan {
  // This is a simple class to queue up a list of copies so they can be
  // executed quickly.
  public:
    accstruct *to, *from;
    int np;
    int qinitialize;   // =1 if one is setting the value, =0 if coadding
    CoaddPlan() { }
    ~CoaddPlan() { }
    inline void plan(accstruct *_to, accstruct *_from, 
    			int _np, int _qinitialize) {
	to = _to; from = _from; np = _np; qinitialize = _qinitialize;
    }
    inline int append(CoaddPlan &next) {
	// Consider whether next could be simply appended to this.
	// Return 0 if successful, 1 if not.
	if (qinitialize==next.qinitialize &&
	    from+np == next.from && 
	    to+np == next.to) {
	    np += next.np; return 0;
        } else return 1;
    }
    inline void execute() {
	if(qinitialize) {
	    memcpy(to, from, np*sizeof(accstruct));
	    // #pragma simd assert
	    // for(int p = 0; p <np; p++)
		// to[p] = from[p];
	} else {
	    #pragma simd assert
	    for(int p = 0; p <np; p++)
		to[p] += from[p];
	}
    }
};



void NearFieldDriver::Finalize(int slab){
    // When all of the SIC are finished, call this routine.
    // It will coadd all of the partial accelerations back into AccSlab,
    // accumulate the timings and statistics, 
    // and then delete the SIC.
    slab = PP->WrapSlab(slab);

#ifndef CUDADIRECT
    return;
#endif
    FinalizeTimer.Start();
    FinalizeBookkeeping.Start();

    assertf(SlabDone(slab) != 0,
            "Finalize called for slab %d but it is not complete\n",slab);

    if(P.ForceOutputDebug)
        CheckInteractionList(slab);
    
    // Allocate the space to hold the co-added NF acceleration
    //? LBW->AllocateArena(AccSlab,slab);

    SetInteractionCollection ** Slices = SlabInteractionCollections[slab];
    int NSplit = SlabNSplit[slab];

    int cpd = P.cpd;
    int nfr = RADIUS;

    // Collect the statistics and timings
    for(int sliceIdx = 0; sliceIdx < NSplit; sliceIdx++){
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

#ifdef OLDPARTIAL
#ifdef OLDCODE
    // We're now going to coadd the partial accelerations
    // This first loop is over Y rows, each handled by a different thread
    #pragma omp parallel for schedule(static)
    for (int j=0; j<cpd; j++) {
        // We're going to find all the Sets associated with this row.
        // There are WIDTH of them, one for each z registration
	CopyPencilToSlabSetup.Start();
        SetInteractionCollection *theseSlices[WIDTH];
        for (int w=0; w<WIDTH; w++)
            theseSlices[w] = NULL;
        for (int sliceIdx=0; sliceIdx< NSplit*WIDTH; sliceIdx++) {
        // Just doing an exhaustive search, so we assume less about
        // the upstream data model.
            if (j >= Slices[sliceIdx]->j_low && j < Slices[sliceIdx]->j_high) {
               // We've found a Slice.
               // Compute which z-registration this is.
               theseSlices[Slices[sliceIdx]->k_mod] = Slices[sliceIdx];
            }
        }
        for (int w=0; w<WIDTH; w++)
            assertf(theseSlices[w] != NULL, "We failed to find all z-registrations");
	CoaddPlan copyplan[cpd*(2*nfr+1)];
	int nplan = 0;

        // We're going to set the acceleration the first time we encounter
        // a cell.  The cells are encountered in z order, starting from 
        // z=0 and increasing.  That means that when we get to z=cpd,
        // we've already been there once.
        int firsttouch = 0;
        
        // Now we're going to proceed through the pencils in the Z direction
        for (int kk=0; kk<cpd; kk++) {
            int k_mod = kk%WIDTH;   // This is which registration we're in
            SetInteractionCollection *Slice = theseSlices[k_mod];
            int k = (kk-Slice->k_mod)/Slice->nfwidth;
            int sinkIdx = (j - Slice->j_low)*Slice->Nk + k;

            // And now we can continue with the previous stuff
            int SinkCount = Slice->SinkSetCount[sinkIdx];
            int zmid = PP->WrapSlab(kk + Slice->nfradius);
            int Start = Slice->SinkSetStart[sinkIdx];

            for(int z=zmid-Slice->nfradius; z <= zmid+Slice->nfradius; z++){
                int CellNP = PP->NumberParticle(slab,j,z);
                accstruct *a = PP->NearAccCell(slab,j,z);

                if(P.ForceOutputDebug == 1){
                    for(int p =0; p < CellNP; p++){
                        assert(isfinite(Slice->SinkSetAccelerations[Start +p].x));
                        assert(isfinite(Slice->SinkSetAccelerations[Start +p].y));
                        assert(isfinite(Slice->SinkSetAccelerations[Start +p].z));
                    }
                }

                // Set when this is the first touch; accumulate otherwise
                if(z==firsttouch && z<cpd){
                    firsttouch++;
		    copyplan[nplan++].plan(a, Slice->SinkSetAccelerations+Start, CellNP, 1);
                } else {
		    copyplan[nplan].plan(a, Slice->SinkSetAccelerations+Start, CellNP, 0);
		    // We attempt to concatentate this one to the previous.
		    // nplan=0 always went to the initialization branch, so
		    // there's no risk that nplan-1 is illegal.
		    if (copyplan[nplan-1].append(copyplan[nplan])) nplan++;
                }

                Start+=CellNP;
                SinkCount-=CellNP;
            }
            assert(SinkCount == 0);
        }
	CopyPencilToSlabSetup.Stop();
	CopyPencilToSlabCopy.Start();
	for (int p=0; p<nplan; p++) copyplan[p].execute();
	CopyPencilToSlabCopy.Stop();
    }
#else   // Here's the new code

    // We need to coadd the 5 partial accelerations
    uint64 NSink = Slab->size(slab);
    accstruct *acctot = (accstruct *)LBW->ReturnIDPtr(AccSlab, slab);
    accstruct *accpart[WIDTH];
    accpart[0] = (accstruct *)LBW->ReturnIDPtr(PartialAccSlab, slab);
    for (int s=1; s<WIDTH; s++) accpart[s] = accpart[s-1] + NSink;

// TODO: Restore this pragma
//    #pragma omp parallel for schedule(static)
    for (int j=0; j<cpd; j++) {
	uint64 start = PP->CellInfo(slab,j,0)->startindex;
	uint64 end = NSink;
	if (j<cpd-1) end = PP->CellInfo(slab,j+1,0)->startindex;
	int nbad[5] = { 0,0,0,0,0};
	for (uint64 p=start; p<end; p++) {
	    accstruct tot;
	    for (int s=0; s<WIDTH; s++) 
	    	if( !(fabs(accpart[s][p].x)<1) ) nbad[s]++;
		    // fprintf(stderr, "Acceleration too big %d %lu\n", s, p);
		    // accpart[s][p].x, accpart[s][p].y, accpart[s][p].z, accpart[s][p].w);
	    tot = accpart[0][p];
	    // TODO: Might be worth writing something that is easier to unroll
	    for (int s=1; s<WIDTH; s++) 
		tot += accpart[s][p];
	    acctot[p] = tot;
	}
	fprintf(stderr,"%d %d %d %d %d bad entries in pencil %d of slab %d\n", 
	    nbad[0], nbad[1], nbad[2], nbad[3], nbad[4], j, slab);
    }

#endif  // End new code
#endif  // End OLDPARTIAL; we don't have any coaddition to do!

    CopyPencilToSlab.Stop();

    // Do a final pass to delete all slices
    TearDownPencils.Start();
    for(int sliceIdx = 0; sliceIdx < NSplit; sliceIdx++){
        SetInteractionCollection *Slice = Slices[sliceIdx];
        delete Slice;
    }
    LBW->DeAllocate(NearField_SIC_Slab, slab);
    //?LBW->DeAllocate(PartialAccSlab, slab);
    delete[] Slices;
    TearDownPencils.Stop();
    
    if(P.ForceOutputDebug) CheckGPUCPU(slab);
    FinalizeTimer.Stop();
}
    



NearFieldDriver *JJ;
    
#include "PencilOnPencil.cc"
#include "PencilPlan.cc"





 /*
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
