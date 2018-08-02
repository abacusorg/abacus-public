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

inline int SetInteractionCollection::index_to_zcen(int k) {
    // Given a (k) of the internal sink indexing, return
    // the z of the central cell
    int kk = k_mod+k*nfwidth+nfradius;   // The central cell
    if (kk<0) kk+=cpd; if (kk>=cpd) kk-=cpd;  // Wrapped
    return kk;
}

// ====================== Constructor =========================

/** Given a slab of size cpd and number of particles np, 
we want to compute the required space for the mallocs for all
SetInteractionCollection in this slab.  This is easier to do for the full
than for a single SIC; no worrying about the Y splits.

NSinkSet = cpd**2
NSourceSet = cpd**2 plus 4*cpd for each split
NPaddedSinks = WIDTH*Nparticles+NFBlocksize*cpd**2
NSinkBlocks = NPaddedSinks/NFBlocksize

The last term is worst case.  Typical case is 0.5 of that;
we use the overage to avoid worrying about the 4*cpd source boundaries,
and just set NSourceSet = NSinkSet.

SinkSetStart int[NSinkSet]
SinkSetCount int[]
SinkPlan     plan[]
SinkSetIdMax int[]

SourceSetStart int[NSourceSet]
SourceSetCount int[]
SourcePlan     plan[]

SinkBlockParentPencil int[NSinkBlocks] 
SinkSetAcceleration   FLOAT4[NPaddedSinks]
SinkSourceInteractionList   int[WIDTH*NSinkSet]
SinkSourceYOffset       FLOAT[]

**/

uint64 MaxSICAllocation(int cpd, int np, int WIDTH) {
    uint64 size = 0;
    size += cpd*cpd*(5*sizeof(int)+2*sizeof(SinkPencilPlan)+WIDTH*sizeof(int)+sizeof(FLOAT));
    size += sizeof(accstruct)*(WIDTH*np + NFBlockSize*cpd*cpd);
    size += sizeof(int)*(WIDTH*np + NFBlockSize*cpd*cpd)/NFBlockSize;
    size += 131072; 	// Just adding in some for alignment and small-problem worst case
    return size;
}

/// Given a pointer to an existing 'buffer' of size 'bsize', return 
/// an aligned pointer 'ptr' that includes space of request bytes.
/// Adjust buffer and bsize on return.
/// Crash if the buffer has run out of space
int posix_memalign_wrap((char *) &buffer, size_t &bsize, (void **)ptr, 
	int alignment, size_t request) {
    int extra = alignment-buffer%alignment;
    assertf(bsize>=request+extra, "Aligned subdivision of buffer has overflowed its space: %l < %l + %d\n", bsize, request, extra);
    buffer += extra; bsize -= extra;
    *ptr = (void *)buffer;   // We're pointing to the correct location
    buffer += request; bsize -= request;
    return 0;
}


/// The constructor of the SIC is where all of the planning happens
/// Given a slab, a range of Y's, and a Z modulus, plus a FOF scale
/// to pass along.
///
/// Once finished, this contains all of the information needed to
/// gather the particle data and pass it along to the GPU, as well
/// as the space to store the results.
///
/// This version takes a pointer to and size of an existing buffer,
/// so that it skips its internal mallocs and instead uses the 
/// buffer.  Important: buffer and bsize are modified and returned
/// with the pointer to and size of the unused space.

SetInteractionCollection::SetInteractionCollection(int slab, int _kmod, int _jlow, int _jhigh, FLOAT _b2, (char *)&buffer, size_t &bsize){
    Construction.Start();

    eps = JJ->eps;
    
    //set known class variables
    CompletionFlag = 0;
    j_low = _jlow;
    j_high = _jhigh;
    cpd = P.cpd;
    k_mod = _kmod;
    SlabId = slab;
    b2 = _b2;
    bytes_to_device = 0, bytes_from_device = 0;
    AssignedDevice = 0;
    
    //useful construction constants
    nfradius = P.NearFieldRadius;
    nfwidth = 2*P.NearFieldRadius+1;
    j_width = j_high-j_low;
    Nk = (P.cpd - k_mod)/nfwidth;
    if(Nk * nfwidth + k_mod < P.cpd) Nk++;


    // Make a bunch of the SinkSet and SourceSet containers
    
    assertf( (uint64)j_width*Nk < INT32_MAX, 
            "The number of sink sets will overflow a 32-bit signed int");
    NSinkSets = j_width * Nk;
    assertf(NSinkSets <= MaxNSink, "NSinkSets (%d) larger than allocated space (MaxNSink = %d)\n", NSinkSets, MaxNSink);
    
    assert(posix_memalign_wrap(buffer, bsize, (void **) &SinkSetStart, 4096, sizeof(int) * NSinkSets) == 0);
    assert(posix_memalign_wrap(buffer, bsize, (void **) &SinkSetCount, 4096, sizeof(int) * NSinkSets) == 0);
    assert(posix_memalign_wrap(buffer, bsize, (void **) &SinkPlan, 4096, sizeof(SinkPencilPlan) * NSinkSets) == 0);
    assert(posix_memalign_wrap(buffer, bsize, (void **) &SinkSetIdMax, 4096, sizeof(int) * NSinkSets) == 0);
    
    int localSinkTotal = 0;

    assertf( (uint64)Nk * (j_width + nfwidth) < INT32_MAX,
        "The number of source sets will overflow a 32-bit signed int");
    NSourceSets = Nk * (j_width + nfwidth -1);
    
    assertf(NSourceSets <= MaxNSource, "NSourceSets (%d) larger than allocated space (MaxNSource = %d)\n", NSourceSets, MaxNSource);
    assertf(NSourceSets <= P.cpd*(P.cpd+nfwidth), "NSourceSets (%d) exceeds SourceSet array allocation (%d)\n", NSourceSets, P.cpd*(P.cpd+nfwidth));
    
    assert(posix_memalign_wrap(buffer, bsize, (void **) &SourceSetStart, 4096, sizeof(int) * NSourceSets) == 0);  
    assert(posix_memalign_wrap(buffer, bsize, (void **) &SourceSetCount, 4096, sizeof(int) * NSourceSets) == 0);
    assert(posix_memalign_wrap(buffer, bsize, (void **) &SourcePlan, 4096, sizeof(SourcePencilPlan) * NSourceSets) == 0);
    
    int localSourceTotal = 0;

    DirectTotal = 0;

    // Next, fill in sink data
    FillSinkLists.Clear(); FillSinkLists.Start();
    // Count the Sinks
    CountSinks.Clear(); CountSinks.Start();

    uint64 skewer_blocks[j_width+nfwidth];   // Number of blocks in this k-skewer
            // This is oversized because we'll re-use for Sources

    #pragma omp parallel for schedule(static) reduction(+:localSinkTotal)
    for(int j = 0; j < j_width; j++){
        int this_skewer_blocks = 0;   // Just to have a local variable
        for(int k = 0; k < Nk; k ++) {
	    int zmid = index_to_zcen(k);
            int sinkindex = j * Nk + k;
            // This loads all of the position pointers into the Plan,
            // but also returns the total size.
            int pencilsize = SinkPlan[sinkindex].load(slab, j + j_low, zmid);
            SinkSetCount[sinkindex] = pencilsize;
            localSinkTotal += pencilsize;
            this_skewer_blocks += NumPaddedBlocks(pencilsize);
        }
        skewer_blocks[j] = this_skewer_blocks;
    }
    CountSinks.Stop();
    
    CalcSinkBlocks.Clear(); CalcSinkBlocks.Start();

    // Cumulate the number of blocks in each skewer, so we know how 
    // to start enumerating.
    uint64 skewer_blocks_start[j_width+1+nfwidth]; 
            // This is oversized because we'll re-use for Sources
    skewer_blocks_start[0] = 0;
    for (int j=0; j<j_width; j++) 
            skewer_blocks_start[j+1] = skewer_blocks_start[j]+skewer_blocks[j];

    NSinkBlocks = skewer_blocks_start[j_width];   // The total for all skewers
    assertf(skewer_blocks_start[j_width]*NFBlockSize < INT32_MAX, 
            "The number of padded sink particles will overflow a 32-bit signed int");

    assertf(NSinkBlocks <= MaxSinkBlocks, "NSinkBlocks (%d) larger than allocated space (MaxSinkBlocks = %d)\n", NSinkBlocks, MaxSinkBlocks);
    assert(posix_memalign_wrap(buffer, bsize, (void **) &SinkBlockParentPencil, 4096, sizeof(int) * NSinkBlocks) == 0);

    // Loop over the pencils to set these:
    //   SinkSetStart[set] points to where the padded particles starts
    //   SinkSetIdMax[set] points to where the padded particles end
    //   SinkBlockParentPencil[block] points back to the (j,k) sinkset enumeration
    //
    
    #pragma omp parallel for schedule(static) 
    for(int j = 0; j < j_width; j++){
        // Figure out where to start the enumeration in this skewer
        int block_start = skewer_blocks_start[j];
        int NPaddedSinks = NFBlockSize*block_start; 

        for(int k = 0; k < Nk; k ++) {
            int sinkset = j * Nk + k;
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
    assert(posix_memalign_wrap(buffer, bsize, (void **) &SinkSetAccelerations, 4096, sizeof(accstruct) * NPaddedSinks) == 0);
    AllocAccels.Stop();
    
    FillSinkLists.Stop();

    // All done with the Sinks.  Now for the Sources.
    FillSourceLists.Start();
    CountSources.Start();

    // Notice that the Source Pencils go from [j_low-NFR, j_high+NFR).

    // We once again need to precompute the enumeration of the padded
    // blocks, totaled by skewer.  

    #pragma omp parallel for schedule(static) reduction(+:localSourceTotal)
    for(int j = 0; j<j_width+nfwidth-1; j++) {
        int this_skewer_blocks = 0;   // Just to have a local variable
        for(int k = 0; k < Nk; k ++) {
	    int zmid = index_to_zcen(k);
            int sourceindex = j * Nk + k;
            // This loads the Plan and returns the length
            int sourcelength = SourcePlan[sourceindex].load(slab,
                            j+j_low-nfradius, zmid);
            SourceSetCount[sourceindex] = sourcelength;
            localSourceTotal += sourcelength;
            this_skewer_blocks += NumPaddedBlocks(sourcelength);
        }
        skewer_blocks[j] = this_skewer_blocks;
    }
    CountSources.Stop();

    CalcSourceBlocks.Clear(); CalcSourceBlocks.Start();

    // Cumulate the number of blocks in each skewer, so we know how 
    // to start enumerating.
    skewer_blocks_start[0] = 0;
    for (int j=0; j<j_width+nfwidth-1; j++) 
            skewer_blocks_start[j+1] = skewer_blocks_start[j]+skewer_blocks[j];

    NSourceBlocks = skewer_blocks_start[j_width+nfwidth-1];   
            // The total for all skewers
    assertf(skewer_blocks_start[j_width+nfwidth-1]*NFBlockSize < INT32_MAX, 
            "The number of padded source particles will overflow a 32-bit signed int");

    // SourceSetStart[set] holds the padded particle index
    #pragma omp parallel for schedule(static) 
    for(int j = 0; j<j_width+nfwidth-1; j++) {
        int block_start = skewer_blocks_start[j];
        int NPaddedSources = NFBlockSize*block_start; 
        for(int k = 0; k < Nk; k ++) {
            int sourceset = j * Nk + k;
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

    InteractionCount = nfwidth*j_width*Nk;
    assertf(InteractionCount <= MaxNSink*nfwidth, "InteractionCount (%d) larger than allocated space (MaxNSink * nfwidth = %d)\n", InteractionCount, MaxNSink * nfwidth);
    assertf((uint64)nfwidth*j_width*Nk < INT32_MAX, 
            "Interaction Count exceeds 32-bit signed int");
    assert(posix_memalign_wrap(buffer, bsize, (void **) &SinkSourceInteractionList, 4096, sizeof(int) * InteractionCount) == 0);
    assert(posix_memalign_wrap(buffer, bsize, (void **) &SinkSourceYOffset, 4096, sizeof(FLOAT) * InteractionCount) == 0);
    FLOAT cellsize = PP->invcpd;
    
    uint64 localDirectTotal = 0;
    #pragma omp parallel for schedule(static) reduction(+:localDirectTotal)
    for(int j = 0; j < j_width; j++){
        int g = omp_get_thread_num();
        assertf(j*Nk + Nk <= NSinkSets, "SinkSetCount array access at %d would exceed allocation %d\n", j*Nk + Nk, NSinkSets);
        for(int k=0; k < Nk; k++) {
	    int zmid = index_to_zcen(k);

            int sinkindex = j*Nk + k;
            int l = nfwidth * sinkindex;
	    // We have nfwidth interactions for each SinkPencil; 
	    // these are packed sequentially

            assertf(l + nfwidth <= InteractionCount, "SinkSourceInteractionList array access at %d would exceed allocation %d\n", l + nfwidth, InteractionCount);
            assertf((j + nfwidth)*(Nk) +  k <= P.cpd*(P.cpd+nfwidth), "SourceSetCount array access at %d would exceed allocation %d\n", (j + nfwidth)*(Nk), P.cpd*(P.cpd+nfwidth));
            for(int y=0;y<nfwidth;y++) {
                int sourceindex = (j + y)*(Nk) +  k;
                SinkSourceInteractionList[l+y] = sourceindex;
                #ifndef GLOBAL_POS
                    SinkSourceYOffset[l+y] = (y-nfradius)*cellsize;
                #else
                    // The y coordinate is (j_low+j+y-nfwidth/2)
                    int tmpy = j_low+j+y-nfradius;
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
    // These are now part of the NearField_SIC_Slab buffer
    //? free(SinkSetStart);
    //? free(SinkSetCount);
    //? free(SinkPlan);
    //? free(SinkSetIdMax);
    //? free(SourceSetStart);
    //? free(SourceSetCount);
    //? free(SourcePlan);
    //? free(SinkBlockParentPencil);
    //? free(SinkSetAccelerations);
    //? free(SinkSourceInteractionList);
    //? free(SinkSourceYOffset);
}



void SetInteractionCollection::SetCompleted(){
    // Call this when the Set is detected as done!
    STDLOG(1,"Completed SIC for slab %d w: %d k: %d - %d\n",SlabId,k_mod,j_low,j_high); 
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

