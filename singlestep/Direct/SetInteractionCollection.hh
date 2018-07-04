//holds information about an group of interactions (i.e. pencil on pencil) to compute directs for 
//These sets can be any arbitrary collections of particles (e.g. groups), but frequently are 5 cell pencils

#ifndef __SIC_HH
#define __SIC_HH

#include <vector>

#ifdef CUDADIRECT
namespace cuda {
#include "driver_types.h"
}
#endif

#include "StructureOfLists.cc"

class CellPencilPlan {
public:
    List3<FLOAT> pos;
    // Cells are assumed to be x[0..N), y[0..N), z[0..N), contiguous,
    // with x[0] being the given position.
    // pos.N holds the count
    
    FLOAT offset;
    // The offset to be applied to x or z, relative to the center cell
};

class SinkPencilPlan {
  public:
    CellPencilPlan cell[2*NFRADIUS+1];
    // The cells are not assumed to be contiguous (e.g., periodic wraps)

    void copy_into_pinned_memory(List3<FLOAT> &pinpos, int start, int total);
    int load(int x, int y, int z);
};


class SourcePencilPlan {
    // The same as above, but with motion in the x direction
  public:
    CellPencilPlan cell[2*NFRADIUS+1];
    // The cells are not assumed to be contiguous (e.g., periodic wraps)

    void copy_into_pinned_memory(List3<FLOAT> &pinpos, int start, int total);
    int load(int x, int y, int z);
};

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

        int *           SinkSetStart; //The index in the Sink Pos/Acc lists where this set begins
        int *           SinkSetCount; //The number of particles in the SinkSet
        int *           SinkSetIdMax; //The sum of the above.  We may even be able to get rid of them and just send this to the GPU.
        SinkPencilPlan *           SinkPlan; // The plan for this pencil
        FLOAT3 *        SinkSetAccelerations; //Where the computed accelerations for the collection will be stored

        volatile int CompletionFlag;

        int         SlabId;
        int         W;
        int         K_low;
        int         K_high;
        int         cpd;
        int         InteractionCount;//How many source cell on sink cell interactions are in this set

        // List3<FLOAT> *  SinkSetPositions; //Position data for particles in all sink sets in the collection
        int             NSinkBlocks; //The number of gpu blocks the sink sets translate into
        int *           SinkBlockParentPencil; //What sink pencil each sink block come from


        int *           SourceSetStart;
        int *           SourceSetCount;
        SourcePencilPlan *           SourcePlan;
        // List3<FLOAT>*   SourceSetPositions;
        int             NSourceBlocks;


        int         NSinkList;
        int         NSourceSets;

        int *       SinkSourceInteractionList;
		// Each SinkPencil will interact with width SourcePencil.
		// This is where we enumerate those interactions.
        FLOAT *       SinkSourceYOffset;
		// Each interaction requires a different Delta(y) to
		// shift the particle positions from the Pencil center.

        uint64 DirectTotal; //Number of directs for this interection collection
        uint64 SinkTotal;
        uint64 SourceTotal;
    
        // Different softenings use different eps
        FLOAT eps;

        int AssignedDevice;
        int Blocking;

        //Methods
        //Count the number of particles in the specified pencil
        // int SourcePencilCount(int slab, int ymid, int zz);
        // int SinkPencilCount(int slab, int ymid, int zz);

        //Copy the particle data for a pencil into the given memory
        // void CreateSinkPencil(int sinkx, int sinky, int sinkz, uint64 OutputOffset);
        // void CreateSourcePencil(int sx, int sy, int nz, uint64 OutputOffset);

	int NumPaddedBlocks(int nparticles);
		// Returns the number of padded blocks for a given true number of particles
	int PaddedSinkCount(int sinkindex);
		// Returns the padded length of the array for this pencil
	int PaddedSourceCount(int sourceindex);
		// Returns the padded length of the array for this pencil

        SetInteractionCollection(int slab, int w, int k_low, int k_high);
        ~SetInteractionCollection();

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
