//holds information about an group of interactions (i.e. pencil on pencil) to compute directs for 
//These sets can be any arbitrary collections of particles (e.g. groups), but frequently are 5 cell pencils

#ifndef __SIC_HH
#define __SIC_HH

#include <vector>

#ifdef CUDADIRECT
#include "driver_types.h"
#endif

class SetInteractionCollection{
    public:
        // Synchronous Timers
        STimer  Construction;
        STimer      FillSinkLists;
        STimer          CountSinks;
        STimer          CalcSinkBlocks;
        STimer          FillSinks;
        STimer      FillSourceLists;
        STimer          CountSources;
        STimer          CalcSourceBlocks;
        STimer          FillSources;
        STimer      FillInteractionList;
    
        // Asynchronous, CPU side
        STimer LaunchDeviceKernels;
    
        // Asynchronous, GPU side
#ifdef CUDADIRECT
        cudaEvent_t CopyStart,CopyStopExecStart,ExecStopResStart,ResStop;
#endif
        float CopyTime;
        float ExecutionTime;
        float CopybackTime;
        float TotalTime;
    
        // Data transfer metrics
        uint64 bytes_to_device, bytes_from_device;
    
        // These are shared by all threads
        static int ActiveThreads;
        static pthread_mutex_t GPUTimerMutex;
        static STimer GPUThroughputTimer;

        int *           SinkSetStart; //The index in the Sink Pos/Acc lists where this set begins
        int *           SinkSetCount; //The number of particles in the SinkSet
        FLOAT3 *        SinkSetAccelerations; //Where the computed accelerations for the collection will be stored

        volatile int CompletionFlag;

        int         SlabId;
        int         W;
        int         K_low;
        int         K_high;
        int         InteractionCount;//How many source cell on sink cell interactions are in this set

        List3<FLOAT> *  SinkSetPositions; //Position data for particles in all sink sets in the collection
        int             NSinkBlocks; //The number of gpu blocks the sink sets translate into
        int *           SinkBlockParentPencil; //What sink pencil each sink block come from


        int *           SourceSetStart;
        int *           SourceSetCount;
        List3<FLOAT>*   SourceSetPositions;
        int             NSourceBlocks;


        int         NSinkList;
        int         NSourceSets;

        int *       SinkSourceInteractionList;

        uint64 DirectTotal; //Number of directs for this interection collection
        uint64 SinkTotal;
        FLOAT eps2;

        int AssignedDevice;
        int Blocking;

        //Methods
        //Count the number of particles in the specified pencil
        int SourcePencilCount(int slab, int ymid, int zz);
        int SinkPencilCount(int slab, int ymid, int zz);

        //Copy the particle data for a pencil into the given memory
        void CreateSinkPencil(int sinkx, int sinky, int sinkz, uint64 OutputOffset);
        void CreateSourcePencil(int sx, int sy, int nz, uint64 OutputOffset);

        SetInteractionCollection(int slab, int w, int k_low, int k_high);
        ~SetInteractionCollection();

        //execute this collection on the GPU
        void GPUExecute(int blocking);
        void CPUExecute();

        //check if this collection has been completed;
        int CheckCompletion();

        //Mark this collection as completed and clean up everything but results
        void SetCompleted();

        //Prints debugging information for this interaction set to stdout
        void PrintInteractions();

        //Fill a cell on cell interaction list with the interactions in this SIC (debugging only)
        void  AddInteractionList( std::vector<int> ** il);

};

#endif