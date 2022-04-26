/** \file This is the code to load up the SetInteractionCollection
 * instances, which is the CPU's organization of the work involved in
 * computing the near-field pencil-on-pencil directs.
*/

// TODO: There is heavy use of signed 32-bit integers in this package.
// For now, DJE has simply added a lot of asserts.

#include "StructureOfLists.cc"
#include "SetInteractionCollection.hh"


// ====================== Helper functions  =========================

/// Given the number of particles, compute and return the number of blocks
/// required.  1..N means 1, N+1..2N means 2, etc.  0 means 0.
inline int SetInteractionCollection::NumPaddedBlocks(int nparticles) {
    return (nparticles+NFBlockSize-1)/NFBlockSize;
}

/// The GPU thread needs to figure out how much space has been allocated.
/// No point to allocate a big array when it's so easy to compute.
inline int SetInteractionCollection::PaddedSinkCount(int sinkindex) {
    return NFBlockSize*NumPaddedBlocks(SinkSetCount[sinkindex]);
}

/// The GPU thread needs to figure out how much space has been allocated
inline int SetInteractionCollection::PaddedSourceCount(int sourceindex) {
    return NFBlockSize*NumPaddedBlocks(SourceSetCount[sourceindex]);
}

/// Given a (k) of the internal sink indexing, return
/// the z of the central cell
inline int SetInteractionCollection::index_to_zcen(int k, int twoD) {
    // The first sink pencil is the first one with all valid cells in the primary region
    int kk = k + k_low + nfradius;    // The central cell
    int pad = twoD ? nfradius : 0;
    if(kk >= node_z_start + node_z_size + pad) kk -= node_z_size + 2*pad;
    return kk;
}

// ====================== Constructor =========================

/** Given a slab of size cpd and number of particles np, 
we want to compute the required space for the mallocs for all
SetInteractionCollection in this slab.  This is easier to do for the full
slab than for a single SIC; no worrying about the Y splits.

NSinkSet = cpd**2
NSourceSet = cpd**2 plus 4*cpd for each split
NPaddedSinks = WIDTH*Nparticles+NFBlocksize*cpd**2
NSinkBlocks = NPaddedSinks/NFBlocksize

The last term in NPaddedSinks is worst case.  Typical case is 0.5 of that.

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
// REMOVE: SinkSetAcceleration   FLOAT4[NPaddedSinks]
SinkSourceInteractionList   int[WIDTH*NSinkSet]
SinkSourceYOffset       FLOAT[]

This is mildly conservative.  The overage is about 
0.5*NFBlockSize*cpd**2*sizeof(accstruct).  For 64-particle blocks
and 40 particles per cell, this is 32 extra acc's per cell compared
to 200 particles.  The plans are also not tiny: about 400 bytes per
cell, which is roughly like 25 particles.  So this is about 10%
conservative.

TODO: We could consider multiplying that factor by 0.6 and then
adding in a bit more of a floor.  Let's see what's happening in 
practice.

**/

/// Return the number of particles in this (i,j) skewer, summing over k.
uint64 NumParticlesInSkewer(int slab, int j, int ghost) {
    // Source skewers use ghost particles; sink skewers don't
    uint64 start, end;
    if(ghost){
        start = CP->CellInfo(slab, j, node_z_start_ghost)->startindex_with_ghost;
        if (j<P.cpd-1)
            end = CP->CellInfo(slab, j+1, node_z_start_ghost)->startindex_with_ghost;
        else
            end = SS->size_with_ghost(slab);
    }
    else {
        start = CP->CellInfo(slab, j, node_z_start)->startindex;
        if (j<P.cpd-1)
            end = CP->CellInfo(slab, j+1, node_z_start)->startindex;
        else
            end = SS->size(slab);
    }

    return end-start;	// This is the number of particles in this skewer
}

/// For one skewer, compute how many sink blocks will be incurred
/// for this one Skewer in the worst case for all of the SIC.
int ComputeSinkBlocks(int slab, int j, int NearFieldRadius) {
    int ghost = 0;
    int np = NumParticlesInSkewer(slab, j, ghost);
    return (2*NearFieldRadius+1)*(np/NFBlockSize)+P.cpd;
}
/// For one skewer, compute how many source blocks will be incurred
/// for this one Skewer in the worst case for all of the SIC.
int ComputeSourceBlocks(int slab, int j, int NearFieldRadius) {
    int np = 0;
    int ghost = 1;
    for (int c=-NearFieldRadius; c<=NearFieldRadius; c++)
        np += NumParticlesInSkewer(slab+c, j, ghost);
    return (np/NFBlockSize+P.cpd);
}

/// Compute a mildly conservative estimate of the space required
/// for all SIC in a slab.
uint64 ComputeSICSize(int cpd, int np, int WIDTH, int NSplit) {
    uint64 size = 0;
    int64 NSinkSet = cpd*node_z_size;
    int64 NSourceSet = (cpd + (WIDTH-1)*NSplit)*node_z_size_with_ghost;
    int64 NPaddedSinks = (int64)WIDTH*np + NSinkSet*NFBlockSize;
    int64 NSinkBlocks = NPaddedSinks/NFBlockSize;

    size += (3*sizeof(int) + sizeof(SinkPencilPlan)) * NSinkSet;
    size += (2*sizeof(int) + sizeof(SourcePencilPlan)) * NSourceSet;
    size += (sizeof(int) + sizeof(FLOAT)) * WIDTH*NSinkSet;
    size += (sizeof(int)) * NSinkBlocks;

    STDLOG(2,"Nsink = %d, NSource = %d, Nblocks = %d, size = %d\n",
    	NSinkSet, NSourceSet, NSinkBlocks, size);

    size += (uint64) 1024*1024*PAGE_SIZE/4096; 	
    	// Just adding in some for alignment and small-problem worst case

    // TODO: adding 10% extra padding in an attempt to get the 1e12 particle sim to run
    return 1.2*size;
}

/// Given a pointer to an existing 'buffer' of size 'bsize', return 
/// an aligned pointer 'ptr' that includes space of request bytes.
/// Adjust buffer and bsize on return.
/// Crash if the buffer has run out of space
int posix_memalign_wrap(char * &buffer, size_t &bsize, void ** ptr, 
	int alignment, size_t request) {
    int extra = alignment-((uintptr_t)buffer)%alignment;
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

SetInteractionCollection::SetInteractionCollection(int slab, int _jlow, int _jhigh, FLOAT _b2, char * &buffer, size_t &bsize, int NearFieldRadius){
    eps = NFD->eps;
    
    //set known class variables
    CompletionFlag = 0;
    j_low = _jlow;
    j_high = _jhigh;
    cpd = P.cpd;
    SlabId = slab;
    b2 = _b2;
    bytes_to_device = 0, bytes_from_device = 0;
    AssignedDevice = 0;

    //useful construction constants
    nfradius = NearFieldRadius;
    nfwidth = 2*NearFieldRadius+1;
    j_width = j_high-j_low;

    // k_low, k_high refer to the non-ghost sink cells
    // But note we do need sink pencils that are centered on ghost cells!
    k_low = node_z_start;
    k_high = node_z_start + node_z_size;
    Nk = (MPI_size_z > 1) ? (k_high - k_low + nfwidth - 1) : (k_high - k_low);

    // Load the Pointers to the PosXYZ Slabs
    SinkPosSlab = (FLOAT *)SB->GetSlabPtr(PosXYZSlab,slab);
    for (int c=0; c<nfwidth; c++) {
	SourcePosSlab[c] = (FLOAT *)SB->GetSlabPtr(PosXYZSlab,slab+c-nfradius);

        // Nslab[] is used to compute the offset of X,Y,Z in PosXYZSlab, which always includes ghosts
        Nslab[c] = SS->size_with_ghost(slab - nfradius + c);
    }

    // There is a slab that has the WIDTH Partial Acceleration fields.
    // Get a pointer to the appropriate segment of that.
    SinkAccSlab = (accstruct *) SB->GetSlabPtr(AccSlab,slab);

    // Make a bunch of the SinkSet and SourceSet containers
    
    assertf( (uint64)j_width*Nk < INT32_MAX, 
            "The number of sink sets will overflow a 32-bit signed int");
    NSinkSets = j_width * Nk;
    assertf(NSinkSets <= MaxNSink, "NSinkSets (%d) larger than allocated space (MaxNSink = %d)\n", NSinkSets, MaxNSink);
    
    assert(posix_memalign_wrap(buffer, bsize, (void **) &SinkSetStart, PAGE_SIZE, sizeof(int) * NSinkSets) == 0);
    assert(posix_memalign_wrap(buffer, bsize, (void **) &SinkSetCount, PAGE_SIZE, sizeof(int) * NSinkSets) == 0);
    assert(posix_memalign_wrap(buffer, bsize, (void **) &SinkPlan, PAGE_SIZE, sizeof(SinkPencilPlan) * NSinkSets) == 0);
    assert(posix_memalign_wrap(buffer, bsize, (void **) &SinkSetIdMax, PAGE_SIZE, sizeof(int) * NSinkSets) == 0);

    assertf( (uint64)Nk * (j_width + nfwidth) < INT32_MAX,
        "The number of source sets will overflow a 32-bit signed int");
    NSourceSets = Nk * (j_width + nfwidth -1);
    
    assertf(NSourceSets <= MaxNSource, "NSourceSets (%d) larger than allocated space (MaxNSource = %d)\n", NSourceSets, MaxNSource);
    assertf(NSourceSets <= cpd*(cpd+nfwidth), "NSourceSets (%d) exceeds SourceSet array allocation (%d)\n", NSourceSets, cpd*(cpd+nfwidth));
    
    assert(posix_memalign_wrap(buffer, bsize, (void **) &SourceSetStart, PAGE_SIZE, sizeof(int) * NSourceSets) == 0);  
    assert(posix_memalign_wrap(buffer, bsize, (void **) &SourceSetCount, PAGE_SIZE, sizeof(int) * NSourceSets) == 0);
    assert(posix_memalign_wrap(buffer, bsize, (void **) &SourcePlan, PAGE_SIZE, sizeof(SourcePencilPlan) * NSourceSets) == 0);

    DirectTotal = 0;

    // Next, fill in sink data
    // Count the Sinks

    uint64 skewer_blocks[j_width+nfwidth];   // Number of blocks in this k-skewer
            // This is oversized because we'll re-use for Sources

    const int twoD = MPI_size_z > 1;  // Don't do a z-wrap in 2D
    int localSinkTotal = 0;
    #pragma omp parallel for schedule(static) reduction(+:localSinkTotal)
    for(int j = 0; j < j_width; j++){
        int this_skewer_blocks = 0;   // Just to have a local variable
        for(int k = 0; k < Nk; k ++) {
	        int zmid = index_to_zcen(k, twoD);
            int sinkindex = j * Nk + k;
            // This loads all of the position pointers into the Plan,
            // but also returns the total size.
            int pencilsize = SinkPlan[sinkindex].load(slab, j + j_low, zmid, nfradius, twoD);
            SinkSetCount[sinkindex] = pencilsize;
            localSinkTotal += pencilsize;
            this_skewer_blocks += NumPaddedBlocks(pencilsize);
        }
        skewer_blocks[j] = this_skewer_blocks;
    }

    SinkTotal = localSinkTotal;

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
    assert(posix_memalign_wrap(buffer, bsize, (void **) &SinkBlockParentPencil, PAGE_SIZE, sizeof(int) * NSinkBlocks) == 0);

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
    PaddedSinkTotal = NPaddedSinks;  // for performance metrics, we always move around the padded amount
            // The total padded number of particles
    assertf(NPaddedSinks <= MaxSinkSize, "NPaddedSinks (%d) larger than allocated space (MaxSinkSize = %d)\n", NPaddedSinks, MaxSinkSize);
    
    // We only do this in CPU mode
    // assert(posix_memalign_wrap(buffer, bsize, (void **) &SinkSetAccelerations, PAGE_SIZE, sizeof(accstruct) * NPaddedSinks) == 0);
    SinkSetAccelerations = NULL;   // Just to give a value

    // All done with the Sinks.  Now for the Sources.

    // Notice that the Source Pencils go from [j_low-NFR, j_high+NFR).

    // We once again need to precompute the enumeration of the padded
    // blocks, totaled by skewer.  

    int localSourceTotal = 0;
    #pragma omp parallel for schedule(static) reduction(+:localSourceTotal)
    for(int j = 0; j<j_width+nfwidth-1; j++) {
        int this_skewer_blocks = 0;   // Just to have a local variable
        for(int k = 0; k < Nk; k++) {
	    int zmid = index_to_zcen(k, twoD);
            int sourceindex = j * Nk + k;
            // This loads the Plan and returns the length
            int sourcelength = SourcePlan[sourceindex].load(slab,
                            j+j_low-nfradius, zmid, nfradius);
            SourceSetCount[sourceindex] = sourcelength;
            localSourceTotal += sourcelength;
            this_skewer_blocks += NumPaddedBlocks(sourcelength);
        }
        skewer_blocks[j] = this_skewer_blocks;
    }

    SourceTotal = localSourceTotal;

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
    PaddedSourceTotal = NPaddedSources;  // for performance metrics, we always move around the padded amount 
    assertf(NPaddedSources <= MaxSourceSize, "NPaddedSources (%d) larger than allocated space (MaxSourceSize = %d)\n", NPaddedSources, MaxSourceSize);

    // Next, we have to pair up the Source and Sinks.  Each sink
    // will be acted on by 5 sources.
    // Fill the interaction lists for the sink sets

    InteractionCount = nfwidth*j_width*Nk;
    assertf(InteractionCount <= MaxNSink*nfwidth, "InteractionCount (%d) larger than allocated space (MaxNSink * nfwidth = %d)\n", InteractionCount, MaxNSink * nfwidth);
    assertf((uint64)nfwidth*j_width*Nk < INT32_MAX, 
            "Interaction Count exceeds 32-bit signed int");
    assert(posix_memalign_wrap(buffer, bsize, (void **) &SinkSourceInteractionList, PAGE_SIZE, sizeof(int) * InteractionCount) == 0);
    assert(posix_memalign_wrap(buffer, bsize, (void **) &SinkSourceYOffset, PAGE_SIZE, sizeof(FLOAT) * InteractionCount) == 0);
    FLOAT cellsize = CP->invcpd;
    
    assertf(j_width*Nk <= NSinkSets, "SinkSetCount size %d would exceed allocation %d\n", j_width*Nk, NSinkSets);

    uint64 localDirectTotal = 0, localPaddedDirectTotal = 0;
    #pragma omp parallel for schedule(static) reduction(+:localDirectTotal) reduction(+:localPaddedDirectTotal)
    for(int j = 0; j < j_width; j++){
        for(int k=0; k < Nk; k++) {
	    int zmid = index_to_zcen(k, twoD);

            int sinkindex = j*Nk + k;
            int l = nfwidth * sinkindex;
	    // We have nfwidth interactions for each SinkPencil; 
	    // these are packed sequentially

            assertf(l + nfwidth <= InteractionCount, "SinkSourceInteractionList array access at %d would exceed allocation %d\n", l + nfwidth, InteractionCount);
            assertf((j + nfwidth)*Nk +  k <= cpd*(cpd+nfwidth), "SourceSetCount array access at %d would exceed allocation %d\n", (j + nfwidth)*(Nk), cpd*(cpd+nfwidth));
            for(int y=0;y<nfwidth;y++) {
                int sourceindex = (j + y)*Nk +  k;
                SinkSourceInteractionList[l+y] = sourceindex;
                #ifndef GLOBAL_POS
                    SinkSourceYOffset[l+y] = (y-nfradius)*cellsize;
                #else
                    // The y coordinate is (j_low+j+y-nfwidth/2)
                    int tmpy = j_low+j+y-nfradius;
                    SinkSourceYOffset[l+y] = (tmpy-CP->WrapSlab(tmpy))*cellsize;
                #endif
                localDirectTotal += SinkSetCount[sinkindex] * SourceSetCount[sourceindex];
                localPaddedDirectTotal += PaddedSinkCount(sinkindex) * SourceSetCount[sourceindex];
            }
        }
    }
    DirectTotal = localDirectTotal;
    PaddedDirectTotal = localPaddedDirectTotal;
}


/// A simple destructor
SetInteractionCollection::~SetInteractionCollection(){
}


/// Call this when the Set is detected as done!
void SetInteractionCollection::SetCompleted(){
    STDLOG(2,"Completed SIC for slab %d, y: %d - %d\n",SlabId,j_low,j_high); 
    CompletionFlag = 1;
}


/// Tell us if the Set is complete
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

