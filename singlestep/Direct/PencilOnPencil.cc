//code related to direct pencil interaction creation for directdriver

//Collection of Interactions of all 5 cell pencils in a specified region of a slab

// TODO: There is heavy use of signed 32-bit integers in this package.
// For now, DJE has simply added a lot of asserts.

#include "StructureOfLists.cc"
#include "SetInteractionCollection.hh"


// ====================== Helper functions  =========================

inline int SetInteractionCollection::NumPaddedBlocks(int nparticles) {
    // Given the number of particles, compute and return the number of blocks
    // required.  1..N means 1, N+1..2N means 2, etc.  0 means 0.
    return (nparticles+NFBlockSize-1)/NFBlockSize;
}

inline int SetInteractionCollection::PaddedSinkCount(int sinkindex) {
    // The GPU thread needs to figure out how much space has been allocated.
    // No point to allocate a big array when it's so easy to compute.
    return NFBlockSize*NumPaddedBlocks(SinkSetCount[sinkindex]);
}

inline int SetInteractionCollection::PaddedSourceCount(int sourceindex) {
    // The GPU thread needs to figure out how much space has been allocated
    return NFBlockSize*NumPaddedBlocks(SourceSetCount[sourceindex]);
}

inline int SetInteractionCollection::index_to_zcen(int j) {
    // Given a (j) of the internal sink indexing, return
    // the z of the central cell
    int jj = W+j*width+nfradius;   // The central cell
    if (jj<0) jj+=cpd; if (jj>=cpd) jj-=cpd;  // Wrapped
    return jj;
}

// ====================== Constructor =========================

SetInteractionCollection::SetInteractionCollection(int slab, int w, int k_low, int k_high, FLOAT _b2){
    // The constructor of the SIC is where all of the planning happens
    // Given a slab, a range of Y's, and a Z modulus, plus a FOF scale
    // to pass along.
    // Once finished, this contains all of the information needed to
    // gather the particle data and pass it along to the GPU, as well
    // as the space to store the results.
    Construction.Start();

    eps = JJ->eps;
    
    //set known class variables
    CompletionFlag = 0;
    K_low = k_low;
    K_high = k_high;
    cpd = P.cpd;
    W = w;
    SlabId = slab;
    b2 = _b2;
    bytes_to_device = 0, bytes_from_device = 0;
    AssignedDevice = 0;
    
    //useful construction constants
    nfradius = P.NearFieldRadius;
    width = 2*P.NearFieldRadius+1;
    k_width = k_high-k_low;
    Nj = (P.cpd - w)/width;
    if(Nj * width + w < P.cpd) Nj++;


    // Make a bunch of the SinkSet and SourceSet containers
    
    assertf( (uint64)k_width*Nj < INT32_MAX, 
            "The number of sink sets will overflow a 32-bit signed int");
    NSinkList = k_width * Nj;
    assertf(NSinkList <= MaxNSink, "NSinkList (%d) larger than allocated space (MaxNSink = %d)\n", NSinkList, MaxNSink);
    
    assert(posix_memalign((void **) &SinkSetStart, 4096, sizeof(int) * NSinkList) == 0);
    assert(posix_memalign((void **) &SinkSetCount, 4096, sizeof(int) * NSinkList) == 0);
    assert(posix_memalign((void **) &SinkPlan, 4096, sizeof(SinkPencilPlan) * NSinkList) == 0);
    assert(posix_memalign((void **) &SinkSetIdMax, 4096, sizeof(int) * NSinkList) == 0);
    
    int localSinkTotal = 0;

    assertf( (uint64)Nj * (k_width + width) < INT32_MAX,
        "The number of source sets will overflow a 32-bit signed int");
    NSourceSets = Nj * (k_width + width -1);
    
    assertf(NSourceSets <= MaxNSource, "NSourceSets (%d) larger than allocated space (MaxNSource = %d)\n", NSourceSets, MaxNSource);
    assertf(NSourceSets <= P.cpd*(P.cpd+width), "NSourceSets (%d) exceeds SourceSet array allocation (%d)\n", NSourceSets, P.cpd*(P.cpd+width));
    
    assert(posix_memalign((void **) &SourceSetStart, 4096, sizeof(int) * NSourceSets) == 0);  // can this be of size NSourceSets?
    assert(posix_memalign((void **) &SourceSetCount, 4096, sizeof(int) * NSourceSets) == 0);
    assert(posix_memalign((void **) &SourcePlan, 4096, sizeof(SourcePencilPlan) * NSourceSets) == 0);
    
    int localSourceTotal = 0;

    DirectTotal = 0;

    // Next, fill in sink data
    FillSinkLists.Clear(); FillSinkLists.Start();
    // Count the Sinks
    CountSinks.Clear(); CountSinks.Start();

    uint64 skewer_blocks[k_width+width];   // Number of blocks in this k-skewer
            // This is oversized because we'll re-use for Sources

    #pragma omp parallel for schedule(static) reduction(+:localSinkTotal)
    for(int k = 0; k < k_width; k++){
        int this_skewer_blocks = 0;   // Just to have a local variable
        for(int j = 0; j < Nj; j ++) {
            /// int jj = w + j * width;
            /// int zmid = PP->WrapSlab(jj + P.NearFieldRadius);
	    int zmid = index_to_zcen(j);
            int sinkindex = k * Nj + j;
            // This loads all of the position pointers into the Plan,
            // but also returns the total size.
            int pencilsize = SinkPlan[sinkindex].load(slab, k + k_low, zmid);
            SinkSetCount[sinkindex] = pencilsize;
            localSinkTotal += pencilsize;
            this_skewer_blocks += NumPaddedBlocks(pencilsize);
        }
        skewer_blocks[k] = this_skewer_blocks;
    }
    CountSinks.Stop();
    
    CalcSinkBlocks.Clear(); CalcSinkBlocks.Start();

    // Cumulate the number of blocks in each skewer, so we know how 
    // to start enumerating.
    uint64 skewer_blocks_start[k_width+1+width]; 
            // This is oversized because we'll re-use for Sources
    skewer_blocks_start[0] = 0;
    for (int k=0; k<k_width; k++) 
            skewer_blocks_start[k+1] = skewer_blocks_start[k]+skewer_blocks[k];

    NSinkBlocks = skewer_blocks_start[k_width];   // The total for all skewers
    assertf(skewer_blocks_start[k_width]*NFBlockSize < INT32_MAX, 
            "The number of padded sink particles will overflow a 32-bit signed int");

    assertf(NSinkBlocks <= MaxSinkBlocks, "NSinkBlocks (%d) larger than allocated space (MaxSinkBlocks = %d)\n", NSinkBlocks, MaxSinkBlocks);
    assert(posix_memalign((void **) &SinkBlockParentPencil, 4096, sizeof(int) * NSinkBlocks) == 0);

    // Loop over the pencils to set these:
    //   SinkSetStart[set] points to where the padded particles starts
    //   SinkSetIdMax[set] points to where the padded particles end
    //   SinkBlockParentPencil[block] points back to the (k,j) sinkset enumeration
    //
    
    #pragma omp parallel for schedule(static) 
    for(int k = 0; k < k_width; k++){
        // Figure out where to start the enumeration in this skewer
        int block_start = skewer_blocks_start[k];
        int NPaddedSinks = NFBlockSize*block_start; 

        for(int j = 0; j < Nj; j ++) {
            int sinkset = k * Nj + j;
            int sinklength = SinkSetCount[sinkset];    // num particles
            SinkSetStart[sinkset] = NPaddedSinks;
            SinkSetIdMax[sinkset] = SinkSetStart[sinkset] + sinklength;
            int nblocks = NumPaddedBlocks(sinklength);
            NPaddedSinks += nblocks*NFBlockSize;

            int block_end = block_start+nblocks;
            for (int b=block_start; b<block_end; b++)
                SinkBlockParentPencil[b] = sinkset;
            block_start=block_end;
        }
    }

    int NPaddedSinks = NFBlockSize*NSinkBlocks;
    SinkTotal = NPaddedSinks;  // for performance metrics, we always move around the padded amount
            // The total padded number of particles
    assertf(NPaddedSinks <= MaxSinkSize, "NPaddedSinks (%d) larger than allocated space (MaxSinkSize = %d)\n", NPaddedSinks, MaxSinkSize);

//    SinkSetPositions = new List3<FLOAT>(NPaddedSinks);
    CalcSinkBlocks.Stop();
    
    AllocAccels.Start();
    assert(posix_memalign((void **) &SinkSetAccelerations, 4096, sizeof(accstruct) * NPaddedSinks) == 0);
    AllocAccels.Stop();
    
    FillSinkLists.Stop();

    // All done with the Sinks.  Now for the Sources.
    FillSourceLists.Start();
    CountSources.Start();

    // Notice that the Source Pencils go from [k_low-NFR, k_high+NFR).

    // We once again need to precompute the enumeration of the padded
    // blocks, totaled by skewer.  

    #pragma omp parallel for schedule(static) reduction(+:localSourceTotal)
    for(int k = 0; k<k_width+width-1; k++) {
        int this_skewer_blocks = 0;   // Just to have a local variable
        for(int j = 0; j < Nj; j ++) {
            /// int jj = w + j * width;
            /// int zmid = PP->WrapSlab(jj + P.NearFieldRadius);
	    int zmid = index_to_zcen(j);
            int sourceindex = k * Nj + j;
            // This loads the Plan and returns the length
            int sourcelength = SourcePlan[sourceindex].load(slab,
                            k+k_low-nfradius, zmid);
            SourceSetCount[sourceindex] = sourcelength;
            localSourceTotal += sourcelength;
            this_skewer_blocks += NumPaddedBlocks(sourcelength);
        }
        skewer_blocks[k] = this_skewer_blocks;
    }
    CountSources.Stop();

    CalcSourceBlocks.Clear(); CalcSourceBlocks.Start();

    // Cumulate the number of blocks in each skewer, so we know how 
    // to start enumerating.
    skewer_blocks_start[0] = 0;
    for (int k=0; k<k_width+width-1; k++) 
            skewer_blocks_start[k+1] = skewer_blocks_start[k]+skewer_blocks[k];

    NSourceBlocks = skewer_blocks_start[k_width+width-1];   
            // The total for all skewers
    assertf(skewer_blocks_start[k_width+width-1]*NFBlockSize < INT32_MAX, 
            "The number of padded source particles will overflow a 32-bit signed int");

    // SourceSetStart[set] holds the padded particle index
    #pragma omp parallel for schedule(static) 
    for(int k = 0; k<k_width+width-1; k++) {
        int block_start = skewer_blocks_start[k];
        int NPaddedSources = NFBlockSize*block_start; 
        for(int j = 0; j < Nj; j ++) {
            int sourceset = k * Nj + j;
            int sourcelength = SourceSetCount[sourceset];
            SourceSetStart[sourceset] = NPaddedSources;
            int nblocks = NumPaddedBlocks(sourcelength);
            NPaddedSources += nblocks*NFBlockSize;
        }
    }
    
    int NPaddedSources = NFBlockSize*NSourceBlocks;
            // The total number of padded sources
    SourceTotal = NPaddedSources;  // for performance metrics, we always move around the padded amount 
//    SourceSetPositions = new List3<FLOAT>(NPaddedSources);
    assertf(NPaddedSources <= MaxSourceSize, "NPaddedSources (%d) larger than allocated space (MaxSourceSize = %d)\n", NPaddedSources, MaxSourceSize);
    CalcSourceBlocks.Stop();
    FillSourceLists.Stop();


    // Next, we have to pair up the Source and Sinks.  Each sink
    // will be acted on by 5 sources.
    // Fill the interaction lists for the sink sets
    FillInteractionList.Start();

    InteractionCount = width*k_width*Nj;
    assertf(InteractionCount <= MaxNSink*width, "InteractionCount (%d) larger than allocated space (MaxNSink * width = %d)\n", InteractionCount, MaxNSink * width);
    assertf((uint64)width*k_width*Nj < INT32_MAX, 
            "Interaction Count exceeds 32-bit signed int");
    assert(posix_memalign((void **) &SinkSourceInteractionList, 4096, sizeof(int) * InteractionCount) == 0);
    assert(posix_memalign((void **) &SinkSourceYOffset, 4096, sizeof(FLOAT) * InteractionCount) == 0);
    FLOAT cellsize = PP->invcpd;
    
    uint64 localDirectTotal = 0;
    #pragma omp parallel for schedule(static) reduction(+:localDirectTotal)
    for(int k = 0; k < k_width; k++){
        int g = omp_get_thread_num();
        assertf(k*Nj + Nj <= NSinkList, "SinkSetCount array access at %d would exceed allocation %d\n", k*Nj + Nj, NSinkList);
        for(int j=0; j < Nj; j++) {
	    int zmid = index_to_zcen(j);

            int sinkindex = k*Nj + j;
            int l = width * sinkindex;
	    // We have width interactions for each SinkPencil; 
	    // these are packed sequentially

            assertf(l + width <= InteractionCount, "SinkSourceInteractionList array access at %d would exceed allocation %d\n", l + width, InteractionCount);
            assertf((k + width)*(Nj) +  j <= P.cpd*(P.cpd+width), "SourceSetCount array access at %d would exceed allocation %d\n", (k + width)*(Nj), P.cpd*(P.cpd+width));
            for(int y=0;y<width;y++) {
                int sourceindex = (k + y)*(Nj) +  j;
                SinkSourceInteractionList[l+y] = sourceindex;
                #ifndef GLOBAL_POS
                    SinkSourceYOffset[l+y] = (y-nfradius)*cellsize;
                #else
                    // The y coordinate is (k_low+k+y-width/2)
                    int tmpy = k_low+k+y-nfradius;
                    SinkSourceYOffset[l+y] = (tmpy-PP->WrapSlab(tmpy))*cellsize;
                #endif
                localDirectTotal += SinkSetCount[sinkindex] * SourceSetCount[sourceindex];
            }
        }
    }
    DirectTotal = localDirectTotal;

    FillInteractionList.Stop();
    Construction.Stop();
}


SetInteractionCollection::~SetInteractionCollection(){
    // A simple destructor
    free(SinkSetStart);
    free(SinkSetCount);
    free(SinkPlan);
    free(SinkSetIdMax);
    free(SourceSetStart);
    free(SourceSetCount);
    free(SourcePlan);
    free(SinkBlockParentPencil);
    free(SinkSetAccelerations);
    free(SinkSourceInteractionList);
    free(SinkSourceYOffset);
}



void SetInteractionCollection::SetCompleted(){
    // Call this when the Set is detected as done!
    STDLOG(1,"Completed SIC for slab %d w: %d k: %d - %d\n",SlabId,W,K_low,K_high); 
    CompletionFlag = 1;
}


int SetInteractionCollection::CheckCompletion(){
    return CompletionFlag;
}

#include "extras_PencilOnPencil.cc"

#ifndef CUDADIRECT
// If we're not compiling the GPU code, need a stub for this function
void SetInteractionCollection::GPUExecute(int){
    QUIT("Abacus was not compiled with CUDA.  Try running ./configure again?\n");
}
#endif

volatile int SetInteractionCollection::ActiveThreads = 0;
pthread_mutex_t SetInteractionCollection::GPUTimerMutex = PTHREAD_MUTEX_INITIALIZER;
STimer SetInteractionCollection::GPUThroughputTimer;

