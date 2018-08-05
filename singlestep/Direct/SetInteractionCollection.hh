// Holds information about an group of interactions (i.e. pencil on pencil) to compute directs for 
// These sets can be any arbitrary collections of particles (e.g. groups), but frequently are 5 cell pencils.

/* We arrange the global calculation into sets of cells in the X
direction acting on sets of cells in the Z direction.  If NR=2, then
there are 5 cells in each Pencil.  A given sink cell will appear in 5
sink pencils, and we compute these 5 partial accelerations in 5 different
SetInteractionCollections, coadding at the end.  The source pencil and
sink pencil overlap at the central x & z, but by forming these pencils
for a range of Y, we can have 5 source pencils act on each sink pencil,
thereby limiting the data transfer to the GPU.  

Using pencils is more GPU efficient than having a single sink cell
because we want to have the sinks fill the GPU blocks.  Further, it
is important that a given source be used on many sinks when it is
loaded: if the GPU ram can load 240 GB/sec of source positions, 
this is "only" 20 Gsources/sec.  But the GPU can compute ~200 Gdirect/sec,
so if we can't provide at least 10 sinks, we will go lacking.
And the load across the PCIe bus is >20x slower than that, so we
want to have >200 sink particles in the 5 sink pencils that a source
pencil will act on.  

The particles in the pencils have been shifted to the cell-centered
coordinate system of the central cell.  In this way, a sink and source
pencils' coordinates differ only in a Y offset, which is added on
the fly.
*/

#ifndef __SIC_HH
#define __SIC_HH

#include <vector>

#ifdef CUDADIRECT
namespace cuda {
#include "driver_types.h"
}
#endif

#include "StructureOfLists.cc"

/// Cells are assumed to be x[0..N), y[0..N), z[0..N), contiguous.
/// NOTE: This routine assumes that start can be held as a 32-bit
/// unsigned integer.  Other places in the code allow for a uint64.
/// TODO: Is that loss of generality ok?
class CellPencilPlan {
public:
    unsigned int start;   //< The starting position of the posXYZ for this cell
    int N;	//< The number of particles in this cell

    //? List3<FLOAT> pos;
    
    #ifdef GLOBAL_POS
    FLOAT offset;
    // The offset to be applied to x or z, relative to the center cell
    // With cell-centered coordinates, we can compute this on the fly,
    // so we'll save the space
    #endif
};

class SinkPencilPlan {
  public:
    CellPencilPlan cell[2*NFRADIUS+1];
    // The cells are not assumed to be contiguous (e.g., periodic wraps)

    void copy_into_pinned_memory(List3<FLOAT> &pinpos, int start, int total,
    	void *SinkPosSlab);
    void copy_from_pinned_memory(void *pinacc, int start, int total, 
    	void *SinkAccSlab);
    int load(int x, int y, int z);
};


class SourcePencilPlan {
    // The same as above, but with motion in the x direction
  public:
    CellPencilPlan cell[2*NFRADIUS+1];
    // The cells are not assumed to be contiguous (e.g., periodic wraps)

    void copy_into_pinned_memory(List3<FLOAT> &pinpos, int start, int total,
    	void **SourcePosSlab);
    int load(int x, int y, int z);
};

/* The Sink particle positions and accelerations, as well as the Source
particle positions, are in large single vectors.  Blocks of NFBlockSize
particles are the unit of computation; individual pencils (not cells)
are always padded out to a multiple of NFBlockSize.

Each Sink Pencil is given an individual index.  There are simple
vectors to hold the particle count (SinkSetCount) and the starting
particle index (SinkSetStart), as well as a plan (SinkPlan) that 
tells how the particles should be loaded from external memory.

Then there are parallel items for the Sources.

There is also a vector SinkSetIdMax that gives the end of each Sink
Pencil.

Next, we have to describe all of the pairwise interactions between
Sinks and Sources.  These are pairs of Pencils, not Blocks. 
There is a large vector of these pairs, organized so that each 
Sink Pencil has its 2*NFR+1 Source pairs contiguous.
The SinkSourceInteractionList holds the Source index.
The SinkSourceYOffset holds the delta(Y) position shift to be applied
in each Pencil on Pencil direct interaction.

Meanwhile the computation itself is organized by thread blocks of
Sink Blocks.  We therefore store a vector SinkBlockParentPencil
that holds for each block what its Sink Pencil id is.  From
that, one can look up the Source Pencil in SinkSourceInteractionList.

The calculation then loops over the Source Pencils, using each 
Source Block.

*/

class SetInteractionCollection{
    public:
        // Synchronous Timers
        STimer  Construction;
        STimer      FillSinkLists;
        STimer          CountSinks;
        STimer          CalcSinkBlocks;
        STimer          AllocAccels;
        STimer      FillSourceLists;
        STimer          CountSources;
        STimer          CalcSourceBlocks;
        STimer      FillInteractionList;
    
        // Asynchronous, CPU side
        STimer DeviceThreadTimer;
        STimer     LaunchDeviceKernels;
        STimer     FillSinks;
        STimer     FillSources;
        STimer     WaitForResult;
        STimer     CopyAccelFromPinned;
    
        // Data transfer metrics
        uint64 bytes_to_device, bytes_from_device;
    
        // These are shared by all threads
        static volatile int ActiveThreads;
        static pthread_mutex_t GPUTimerMutex;
        static STimer GPUThroughputTimer;

	void  *SinkPosSlab;    //< Ptr to the start of the PosXYZ slab
	void  **SourcePosSlab;    //< Ptr to the start of the PosXYZ slabs for Sources
	void *SinkAccSlab;  //< Ptr to the start of the slab of accelerations for this sink slab
	//? void *SinkPartialAccSlab;  //< Ptr to the start of the slab of accelerations for this Z modulus.

        int *           SinkSetStart; //The index in the Sink Pos/Acc lists where this set begins
        int *           SinkSetCount; //The number of particles in the SinkSet
        int *           SinkSetIdMax; //The sum of the above, i.e., the end of the Sink Pencil.  
		// TODO: Might prefer to code this out, but it is vastly less than the pos/acc data
        SinkPencilPlan *           SinkPlan; // The plan for this pencil
        accstruct *        SinkSetAccelerations; //Where the computed accelerations for the collection will be stored

        int *           SourceSetStart;  // The index where the Source Pencils start
        int *           SourceSetCount;  // The number of particles in the Source pencils
        SourcePencilPlan *           SourcePlan;   // The plans to load the Source pencils

	FLOAT b2;	// The square radius for FOF density computation

        volatile int CompletionFlag;

        int         SlabId;	// The X slab number
        //? int         k_mod;	// The modulus of the Z cells, range [0,nfwidth)
	int	    Nk;		// The number of Z cells
        int         j_low;	// The minimum sink Y
        int         j_high;	// The maximum sink Y+1 [j_low,j_high)
	int	    j_width;	// The number of Y cells

        int         cpd;
	int 	    nfradius;	// The NearField Radius
	int 	    nfwidth;	// The NearField Diameter (2*R+1)

        int         InteractionCount;//How many source cell on sink cell interactions are in this set

        int         NSinkBlocks; //The number of gpu blocks the sink sets translate into
        int *       SinkBlockParentPencil; //What sink pencil each sink block come from

        int         NSourceBlocks;   // The number of Source blocks


        int         NSinkSets;	// The number of Sink Pencils
        int         NSourceSets;	// The number of Source Pencils

        int *       SinkSourceInteractionList;
		// Each SinkPencil will interact with nfwidth SourcePencil.
		// This is where we enumerate those interactions.
        FLOAT *       SinkSourceYOffset;
		// Each interaction requires a different Delta(y) to
		// shift the particle positions from the Pencil center.

        uint64 DirectTotal; //Number of directs for this interection collection
		// This is the unpadded, science-useful count of directs
        uint64 SinkTotal;	// The number of total sinks (padded out to blocks)
        uint64 SourceTotal;	// The number of total sources (padded)
    
        // Different softenings use different eps
        FLOAT eps;

        int AssignedDevice;
        int Blocking;

        //Methods

	// Constructor
        SetInteractionCollection(int slab, int _kmod, int _jlow, int _jhigh, FLOAT _b2, char * &buffer, size_t &bsize);
        ~SetInteractionCollection();

	int NumPaddedBlocks(int nparticles);
		// Returns the number of padded blocks for a given true number of particles
	int PaddedSinkCount(int sinkindex);
		// Returns the padded length of the array for this pencil
	int PaddedSourceCount(int sourceindex);
		// Returns the padded length of the array for this pencil

	int index_to_zcen(int j);
		// Returns the zcentral cell given the internal j indexing


        //execute this collection on the GPU
        void GPUExecute(int blocking);
        void CPUExecute();
        void Unpin();

        //check if this collection has been completed;
        int CheckCompletion();

        //Mark this collection as completed and clean up everything but results
        void SetCompleted();

        //Prints debugging information for this interaction set to stdout
        void PrintInteractions();

        //Fill a cell on cell interaction list with the interactions in this SIC (debugging only)
        void  AddInteractionList( std::vector<uint64> ** il);

};

#endif
