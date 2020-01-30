/* Parameters.cpp

The Parameters class contains the time-independent global variables
for the simulation.  These get read from an ASCII parameter file
via ParseHeader.

NB: When adding parameters, you should add:

1) a variable to the class definition,
2) an installscalar/vector line to the constructor, 
3) (optional) a validation check in ValidateParameters.

*/

#ifndef __PARAMETERS_CPP
#define __PARAMETERS_CPP

#define MAX_LINE_LENGTH 1024
#define STRUNDEF "NONE"

#include "ParseHeader.hh"
#include "file.cpp"

#include <sys/sysinfo.h>

class Parameters: public ParseHeader {
public:
    
    char SimName[1024]; //What to call this run
    long long int np;
    int cpd;
    int order;

    int NearFieldRadius;    // Radius of cells in the near-field
    double SofteningLength; // Softening length in the same units as BoxSize

    int  DerivativeExpansionRadius;
    int  MAXRAMMB;
    float  ConvolutionCacheSizeMB; // Set to manually override the detected cache size; this is for L3
    float  ConvolutionL1CacheSizeMB; // Set to manually override the detected cache size
    int RamDisk;        // ==0 for a normal disk, ==1 for a ramdisk (which don't have DIO support)  // TODO: automatically detect this, or at least provide per-directory options
    int ForceBlockingIO;   // ==1 if you want to force all IO to be blocking.
    char StateIOMode[64];  //  "normal", "slosh", "overwrite", "stripe"
    char Conv_IOMode[64];  //  "normal", "slosh", "overwrite", "stripe"
    
    int OMP_NUM_THREADS;  // Number of OpenMP threads.  0 does not modify the system value (usually OMP_NUM_THREADS, or all threads).
                        // Negative values use that many fewer than the max.

    int  DirectNewtonRaphson;  // 0 or 1 

    char DerivativesDirectory[1024];

    char InitialConditionsDirectory[1024];   // The initial condition file name
    char ICFormat[1024];                // The format of the IC files
    int FlipZelDisp;                    // If non-zero and using ICFormat = Zeldovich, flip the Zeldovich displacements
    double ICPositionRange;                // The box size of the IC positions, 
            // in file units.  If ==0, then will default to BoxSize;
    double ICVelocity2Displacement;        // The conversion factor from file velocities
            // to redshift-space comoving displacements (at the IC redshift!).
        // =-1 to instead supply file velocities in km/s.
        // =0 to force the initial velocities to zero!
    double NumSlabsInsertList;                
         // The amount of space to allocate for the insert list, in units 
         // of np/cpd particles.  Set =0 to allocate the full np particles.
         // Default is 2.
    double NumSlabsInsertListIC;                
         // The amount of space to allocate for the insert list, in units 
         // of np/cpd particles.  Set =0 to allocate the full np particles.
         // This parameter applies only to the IC step.
         // Recommend either 4 or 0.

    char ReadStateDirectory[1024];  // Where the input State lives
    char WriteStateDirectory[1024]; // Where the output State lives

    // The node-local directories in the parallel verison
    // These will be equal to the global directories if not given
    char LocalWorkingDirectory[1024];        // node-local working directory
    char LocalReadStateDirectory[1024];
    char LocalWriteStateDirectory[1024];

    char MultipoleDirectory[1024];
    char TaylorDirectory[1024];
    char MultipoleDirectory2[1024];  // for splitting even/odd multipoles
    char TaylorDirectory2[1024];
    char WorkingDirectory[1024];        // If Read/Write/Past not specified, where to put the states
    char LogDirectory[1024];
    char OutputDirectory[1024];     // Where the outputs go
    char GroupDirectory[1024]; //Where the Group files go
    char BackupDirectory[1024]; // The directory from which to restore backups (the Python code also writes backups here)
    
    int OutputEveryStep; //Force timeslices to be output every step if 1
    char OutputFormat[1024];                // The format of the Output files
    int  OmitOutputHeader;                // =1 if you want to skip the ascii header

    double FinalRedshift;        // When to stop.  This will override TimeSliceRedshifts.
	#define MAX_TIMESLICE_REDSHIFTS 1024
    int    nTimeSlice;
    int    nTimeSliceSubsample;
    double TimeSliceRedshifts[MAX_TIMESLICE_REDSHIFTS];
	double TimeSliceRedshifts_Subsample[MAX_TIMESLICE_REDSHIFTS];

    #define NUM_SUBSAMPLES 2
	double ParticleSubsampleA; //a consistently sampled small fraction of particles to output during some output steps. 
	double ParticleSubsampleB; //an additional, unique fraction of particles to output during some output steps. 
	
    double H0;          // The Hubble constant in km/s/Mpc
    double Omega_M;
    double Omega_Smooth;
    double Omega_DE;
    double Omega_K;
    double w0;          // w(z) = w_0 + (1-a)*w_a
    double wa;

    double BoxSize;
    int hMpc;           // =1 if we're using Mpc/h units.  =0 if Mpc units
    double InitialRedshift;
    int LagrangianPTOrder;  // =1 for Zel'dovich, =2 for 2LPT, =3 for 3LPT

    int GroupRadius;        // Maximum size of a group, in units of cell sizes
    double TimeStepAccel;         // Time-step parameter based on accelerations
    double TimeStepDlna;        // Maximum time step in d(ln a)

    // Could have microstepping instructions
    // Could have group finding or coevolution set instructions

    int  NLightCones;
    int  OutputFullLightCones;

    double LightConeOrigins[24];  // Same units as BoxSize

    char LightConeDirectory[1024];

    int PowerSpectrumStepInterval;
    int PowerSpectrumN1d; //1D number of bins to use in the powerspectrum

    int LogVerbosity;   // If 0, production-level log; higher numbers are more verbose
    int StoreForces; // If 1, store the accelerations
    int ForceOutputDebug; // If 1, output near and far forces seperately. 

    int MakeGlass; // Reverse the sign of the acceleration

    int ForceCPU; //IF 1, force directs to be executed exclusively on cpu even if cuda is available
    int GPUMinCellSinks;// If AVX directs are compiled, cells with less than this many particles go to cpu
    int ProfilingMode;//If 1, enable profiling mode, i.e. delete the write-state after creating it to run repeatedly on same dat
    
    #define MAX_IO_THREADS 16
    int IOCores[MAX_IO_THREADS];  // The cores that the IO threads will be bound to.  -1 means don't bind
    #define MAX_IODIRS 100
    char **IODirs; //[MAX_IODIRS][1024];
    int nIODirs;
    int IODirThreads[MAX_IODIRS];
    
    int Conv_OMP_NUM_THREADS;
    int Conv_IOCores[MAX_IO_THREADS];
    int Conv_zwidth;
    
    // TODO: this scheme doesn't account for more complicated NUMA architectures
    int GPUThreadCoreStart[MAX_GPUS];  // The core on which to start placing GPU device threads.
    int NGPUThreadCores;  // The number of free cores on which to place GPU device threads.
    int GPUQueueAssignments[MAX_GPUS];  // The work queue assignments
    int DirectBPD;

    double DensityKernelRad;  // The kernel Radius for the density computation, specified in units of the interparticle spacing.  0 will default to FoFLinkingLength[0]
    double L0DensityThreshold;  // The kernel density required for a particle to be eligible to be in a L0 group; specified in units of the cosmic mean density. This is ignored (uses 0) if DensityKernelRad==0.  Value = 0 triggers code to make a particle eligible if it has any non-self neighbor within DensityKernelRad

    int AllowGroupFinding;
    double FoFLinkingLength[3]; //Linking lengths for level 0,1,2 groupfinding in fractional interparticle spacing 
    double SODensity[2];  // Overdensities for SO groupfinding level 1 and 2
    int MinL1HaloNP; // minimum L1 halo size to output
	float L1Output_dlna;  // minimum delta ln(a) between L1 halo outputs
    float SO_RocheCoeff; 
    float SO_NPForMinDensity; 
    int SO_EvolvingThreshold; //allow evolving (redshift-dependent) density threshold 

    #define MAX_L1OUTPUT_REDSHIFTS 1024
    int nTimeSliceL1; 
    float L1OutputRedshifts[MAX_L1OUTPUT_REDSHIFTS];
    int OutputAllHaloParticles;  // ==0 normally, to output only taggable L1 particles.  If non-zero, output all particles

    double MicrostepTimeStep; // Timestep parameter that controls microstep refinement

    long long int MaxPID;  // Maximum PID to expect.  A PID equal or larger than this indicates corruption of some sort.  0 means NP; -1 means don't check.

    int ProperSoftening;  // Keep the softening length fixed in proper coordinates.  SofteningLength is specified at z=0.
    double SofteningMax;  // The maximum comoving softening to allow when using ProperSoftening

    int NoChecksum;  // Explicitly disable output checksumming

    int UsePinnedGPUMemory;  // Whether to pin the CPU-side GPU staging buffers

    double MPICallRateLimit_ms;  // Enforce a delay between MPI_Test calls

    int UseMunmapThread;  // dedicated munmap() thread in ArenaAllocator

    int MunmapThreadCore;  // Core binding for disposal thread

    // Return the L{tier} size in MB
    float getCacheSize(int tier){
        int cache_size = 0;
        FILE *fp = NULL;
        char fn[1024];

        // find the last-level cache
        for(int i = 0; i<=tier; i++){
            sprintf(fn, "/sys/devices/system/cpu/cpu0/cache/index%d/size", i);
            fp = fopen(fn, "r");
            if(fp == NULL)
                break;
            int nscan = fscanf(fp, "%dK", &cache_size);  // cache size in KB
            fclose(fp);
            assertf(nscan == 1, "Unexpected cache size file format (\"%s\")\n", fn);
        }
        return (float)cache_size/1024.0;    // to MB
    }
    
    // in MB
    int getRAMSize(){
        struct sysinfo info;
        int status = sysinfo(&info);
        assert(status == 0);
        return (int) ((double) info.totalram*info.mem_unit / 1024./1024.);
    }

    Parameters() {
        installscalar("NP",np, MUST_DEFINE);
        installscalar("CPD",cpd,MUST_DEFINE);
        installscalar("Order",order,MUST_DEFINE);

        installscalar("NearFieldRadius",NearFieldRadius,MUST_DEFINE);    // Radius of cells in the near-field
        installvector("SofteningLength", &SofteningLength, 2, 0, MUST_DEFINE); // Softening length in the same units as BoxSize
        installscalar("DerivativeExpansionRadius", DerivativeExpansionRadius,MUST_DEFINE);
        MAXRAMMB = getRAMSize();
        installscalar("MAXRAMMB", MAXRAMMB, DONT_CARE);
        ConvolutionCacheSizeMB = getCacheSize(4);
        installscalar("ConvolutionCacheSizeMB", ConvolutionCacheSizeMB, DONT_CARE);
        ConvolutionL1CacheSizeMB = getCacheSize(1);
        installscalar("ConvolutionL1CacheSizeMB", ConvolutionL1CacheSizeMB, DONT_CARE);
        RamDisk = 0;
        installscalar("RamDisk",RamDisk,DONT_CARE);
        ForceBlockingIO = 0;
        installscalar("ForceBlockingIO",ForceBlockingIO,DONT_CARE);

        sprintf(StateIOMode, "normal");
        installscalar("StateIOMode", StateIOMode, DONT_CARE);
        sprintf(Conv_IOMode, "normal");
        installscalar("Conv_IOMode", Conv_IOMode, DONT_CARE);

        OMP_NUM_THREADS = 0;
        installscalar("OMP_NUM_THREADS",OMP_NUM_THREADS,DONT_CARE);

        DirectNewtonRaphson = 1;
        installscalar("DirectNewtonRaphson",DirectNewtonRaphson,DONT_CARE);  // 0 or 1

        installscalar("DerivativesDirectory",DerivativesDirectory,MUST_DEFINE);

        installscalar("InitialConditionsDirectory",InitialConditionsDirectory,DONT_CARE);   // The initial condition file name
        installscalar("ICFormat",ICFormat,DONT_CARE);   // The initial condition file format
        ICPositionRange = -1.;
        installscalar("ICPositionRange",ICPositionRange,DONT_CARE);   // The initial condition file position convention
        ICVelocity2Displacement = -1;
        installscalar("ICVelocity2Displacement",ICVelocity2Displacement,DONT_CARE);   // The initial condition file velocity convention
        FlipZelDisp = 0;
        installscalar("FlipZelDisp",FlipZelDisp,DONT_CARE);   // Flip Zeldovich ICs

        NumSlabsInsertList = 2.0;
        installscalar("NumSlabsInsertList",NumSlabsInsertList,DONT_CARE);   
        NumSlabsInsertListIC = 4.0;
        installscalar("NumSlabsInsertListIC",NumSlabsInsertListIC,DONT_CARE);   

        sprintf(ReadStateDirectory,STRUNDEF);
        sprintf(WriteStateDirectory,STRUNDEF);
        sprintf(WorkingDirectory,STRUNDEF);
        sprintf(MultipoleDirectory2,STRUNDEF);
        sprintf(TaylorDirectory2,STRUNDEF);
        sprintf(BackupDirectory,STRUNDEF);

        sprintf(LocalWorkingDirectory,STRUNDEF);
        sprintf(LocalReadStateDirectory,STRUNDEF);
        sprintf(LocalWriteStateDirectory,STRUNDEF);

        installscalar("ReadStateDirectory",ReadStateDirectory,DONT_CARE);  // Where the input State lives
        installscalar("WriteStateDirectory",WriteStateDirectory,DONT_CARE); // Where the output State lives
        installscalar("WorkingDirectory",WorkingDirectory,DONT_CARE);
        installscalar("MultipoleDirectory",MultipoleDirectory,MUST_DEFINE);
        installscalar("TaylorDirectory",TaylorDirectory,MUST_DEFINE);
        installscalar("MultipoleDirectory2",MultipoleDirectory2,DONT_CARE);
        installscalar("TaylorDirectory2",TaylorDirectory2,DONT_CARE);
    	installscalar("LogDirectory",LogDirectory,MUST_DEFINE);
    	installscalar("OutputDirectory",OutputDirectory,MUST_DEFINE);     // Where the outputs go
        installscalar("GroupDirectory",GroupDirectory,MUST_DEFINE);
        installscalar("BackupDirectory",BackupDirectory,DONT_CARE);

        installscalar("LocalWorkingDirectory",LocalWorkingDirectory,DONT_CARE);
        installscalar("LocalReadStateDirectory",LocalReadStateDirectory,DONT_CARE);
        installscalar("LocalWriteStateDirectory",LocalWriteStateDirectory,DONT_CARE);

    	OutputEveryStep = 0;
    	installscalar("OutputEveryStep",OutputEveryStep,DONT_CARE);

        installscalar("LightConeDirectory",LightConeDirectory,MUST_DEFINE); //Where the lightcones go. Generally will be the same as the Output directory
        installscalar("NLightCones",NLightCones,DONT_CARE); //if not set, we assume 0
        OutputFullLightCones = 0;
        installscalar("OutputFullLightCones",OutputFullLightCones,DONT_CARE); //if not set, we assume 0
        installvector("LightConeOrigins",LightConeOrigins,24,1,DONT_CARE);

        FinalRedshift = -2.0;        // If <-1, then we will cascade back to the minimum of the TimeSliceRedshifts list
        installscalar("FinalRedshift",FinalRedshift,DONT_CARE);
		
        for (int i = 0; i < MAX_TIMESLICE_REDSHIFTS; i++)
            TimeSliceRedshifts[i] = -2;		
        installvector("TimeSliceRedshifts",TimeSliceRedshifts,MAX_TIMESLICE_REDSHIFTS,1,DONT_CARE);
		
        for (int i = 0; i < MAX_TIMESLICE_REDSHIFTS; i++)
            TimeSliceRedshifts_Subsample[i] = -2;		
        installvector("TimeSliceRedshifts_Subsample",TimeSliceRedshifts_Subsample,MAX_TIMESLICE_REDSHIFTS,1,DONT_CARE);

        for (int i = 0; i < MAX_L1OUTPUT_REDSHIFTS; i++)
            L1OutputRedshifts[i] = -2;
        installvector("L1OutputRedshifts", L1OutputRedshifts, MAX_L1OUTPUT_REDSHIFTS, 1, DONT_CARE);

        ParticleSubsampleA = 0.;
        ParticleSubsampleB = 0.;
        installscalar("ParticleSubsampleA", ParticleSubsampleA, DONT_CARE);
        installscalar("ParticleSubsampleB", ParticleSubsampleB, DONT_CARE); 
		
        strcpy(OutputFormat,"RVdouble");
        // strcpy(OutputFormat,"Packed");
        installscalar("OutputFormat",OutputFormat,DONT_CARE);
        OmitOutputHeader = 0;
        installscalar("OmitOutputHeader",OmitOutputHeader,DONT_CARE);

        installscalar("H0", H0, MUST_DEFINE);
        installscalar("Omega_M", Omega_M, MUST_DEFINE);
        Omega_Smooth = 0.0;
        installscalar("Omega_Smooth", Omega_Smooth, DONT_CARE);
        installscalar("Omega_DE", Omega_DE, MUST_DEFINE);
        installscalar("Omega_K", Omega_K, MUST_DEFINE);
        installscalar("w0", w0, MUST_DEFINE);
        installscalar("wa", wa, MUST_DEFINE);

        installscalar("BoxSize",BoxSize,MUST_DEFINE);
        installscalar("hMpc",hMpc,MUST_DEFINE);           // =1 if we're using Mpc/h units.  =0 if Mpc units
        installscalar("InitialRedshift",InitialRedshift,MUST_DEFINE);
        installscalar("LagrangianPTOrder",LagrangianPTOrder,MUST_DEFINE);  // =1 for Zel'dovich, =2 for 2LPT, =3 for 3LPT

        installscalar("GroupRadius",GroupRadius,MUST_DEFINE);        // Maximum size of a group, in units of cell sizes
        installscalar("TimeStepAccel",TimeStepAccel,MUST_DEFINE);         // Time-step parameter based on accelerations
        installscalar("TimeStepDlna",TimeStepDlna,MUST_DEFINE);        // Maximum time step in d(ln a)

        LogVerbosity = 1;
        installscalar("LogVerbosity",LogVerbosity, DONT_CARE);
        StoreForces = 0;
        installscalar("StoreForces",StoreForces, DONT_CARE);
        
        ForceOutputDebug = 0;
        installscalar("ForceOutputDebug",ForceOutputDebug,DONT_CARE);

        MakeGlass = 0;
        installscalar("MakeGlass",MakeGlass,DONT_CARE);

        ForceCPU = 0;
        installscalar("ForceCPU",ForceCPU,DONT_CARE);
        GPUMinCellSinks = 0;
        installscalar("GPUMinCellSinks",GPUMinCellSinks,DONT_CARE);
        ProfilingMode=0;
        installscalar("ProfilingMode",ProfilingMode,DONT_CARE);

        installscalar("SimName",SimName,MUST_DEFINE);
        installscalar("PowerSpectrumStepInterval",PowerSpectrumStepInterval,DONT_CARE);
        installscalar("PowerSpectrumN1d",PowerSpectrumN1d,DONT_CARE);
        PowerSpectrumStepInterval = -1; //Do not calculate OTF powerspectra
        PowerSpectrumN1d = 1;
        hs = NULL;

        // default means don't bind to core
        for (int i = 0; i < MAX_IO_THREADS; i++)
            IOCores[i] = -1;
        installvector("IOCores", IOCores, MAX_IO_THREADS, 1, DONT_CARE);

        for (int i = 0; i < MAX_IO_THREADS; i++)
            Conv_IOCores[i] = -1;
        installvector("Conv_IOCores", Conv_IOCores, MAX_IO_THREADS, 1, DONT_CARE);

        Conv_OMP_NUM_THREADS = 0;
        installscalar("Conv_OMP_NUM_THREADS", Conv_OMP_NUM_THREADS, DONT_CARE);

        Conv_zwidth = -1;
        installscalar("Conv_zwidth", Conv_zwidth, DONT_CARE);

        // Using staticly allocated memory didn't seem to work with installvector
        IODirs = (char**) malloc(MAX_IODIRS*sizeof(char*));
        char *block = (char *) malloc(MAX_IODIRS*1024*sizeof(char));
        for(int i = 0; i < MAX_IODIRS; i++){
            IODirs[i] = block + 1024*i;
            strcpy(IODirs[i], STRUNDEF);
            IODirThreads[i] = -1;
        }
        installvector("IODirs", IODirs, MAX_IODIRS, 1024, DONT_CARE);
        nIODirs = 0;
        installscalar("nIODirs", nIODirs, DONT_CARE);
        installvector("IODirThreads", IODirThreads, MAX_IODIRS, 1, DONT_CARE);

        // If GPUThreadCoreStart is undefined, GPU threads will not be bound to cores
        for(int i = 0; i < MAX_GPUS; i++){
            GPUThreadCoreStart[i] = -1;
            GPUQueueAssignments[i] = i;  // one per GPU
        }
        installvector("GPUThreadCoreStart", GPUThreadCoreStart, MAX_GPUS, 1, DONT_CARE);
        installvector("GPUQueueAssignments", GPUQueueAssignments, MAX_GPUS, 1, DONT_CARE);
        NGPUThreadCores = -1;
        installscalar("NGPUThreadCores", NGPUThreadCores, DONT_CARE);

        DirectBPD = 3;
        installscalar("DirectBPD", DirectBPD, DONT_CARE);

	    DensityKernelRad = 0.0; // if no value is given, this will default to L0 linking length. 
        installscalar("DensityKernelRad",DensityKernelRad, DONT_CARE);
	    L0DensityThreshold = 75.0; 
        installscalar("L0DensityThreshold",L0DensityThreshold, DONT_CARE); 

        AllowGroupFinding = 1;
        installscalar("AllowGroupFinding",AllowGroupFinding, DONT_CARE);
        FoFLinkingLength[0] = .25;
        FoFLinkingLength[1] = .186;
        FoFLinkingLength[2] = .138;
        installvector("FoFLinkingLength",FoFLinkingLength,3,1,DONT_CARE);
        SODensity[0] = 200.0;
        SODensity[1] = 800.0;
        installvector("SODensity",SODensity,2,1,DONT_CARE);
        MinL1HaloNP = 40;
        installscalar("MinL1HaloNP", MinL1HaloNP, DONT_CARE);
		L1Output_dlna = -1;
		installscalar("L1Output_dlna", L1Output_dlna, DONT_CARE);

        SO_RocheCoeff = 2.0; 
        installscalar("SO_RocheCoeff", SO_RocheCoeff, DONT_CARE);
        SO_NPForMinDensity = 35.0; 
        installscalar("SO_NPForMinDensity", SO_NPForMinDensity, DONT_CARE);
        SO_EvolvingThreshold = 1; 
        installscalar("SO_EvolvingThreshold", SO_EvolvingThreshold, DONT_CARE); 

        OutputAllHaloParticles = 0;
        installscalar("OutputAllHaloParticles", OutputAllHaloParticles, DONT_CARE);

        MicrostepTimeStep = 0.;  // no microstepping
        installscalar("MicrostepTimeStep", MicrostepTimeStep, DONT_CARE);

        MaxPID = -1;
        installscalar("MaxPID", MaxPID, DONT_CARE);

        ProperSoftening = 0;
        installscalar("ProperSoftening", ProperSoftening, DONT_CARE);

        SofteningMax = DBL_MAX;
        installscalar("SofteningMax", SofteningMax, DONT_CARE);

        NoChecksum = 0;
        installscalar("NoChecksum", NoChecksum, DONT_CARE);

        UsePinnedGPUMemory = -1;  // auto
        installscalar("UsePinnedGPUMemory", UsePinnedGPUMemory, DONT_CARE);

        MPICallRateLimit_ms = 1;
        installscalar("MPICallRateLimit_ms", MPICallRateLimit_ms, DONT_CARE);

        UseMunmapThread = 1;
        installscalar("UseMunmapThread", UseMunmapThread, DONT_CARE);

        MunmapThreadCore = -1;
        installscalar("MunmapThreadCore", MunmapThreadCore, DONT_CARE);
    }

    // We're going to keep the HeaderStream, so that we can output it later.
    HeaderStream *hs;
    ~Parameters(void) {
        delete hs;
        free(IODirs[0]);
        free(IODirs);
    }
    char *header() { 
        assert(hs!=NULL); assert(hs->buffer!=NULL);
        return hs->buffer;        // This is just a standard C-style string.
    }


    void ReadParameters(char *paramaterfile, int icflag);
    void ValidateParameters(void);
    void register_vars();

    double ppd() {
        // return the cube root of np, but be careful to avoid round-off 
        // of perfect cubes.
        double _ppd = pow((double)np, 1.0/3.0);
        if ( fabs(_ppd-floor(_ppd+0.1))<1e-10 ) _ppd = floor(_ppd+0.1);
        return _ppd;
    }
    int is_np_perfect_cube() {
        // Return 1 if np is a perfect cube.
        long long int n = floor(ppd());
        if (n*n*n==np) return 1; else return 0;
    }

    double FinishingRedshift() {
        // Return the redshift where the code should halt.
        // FinalRedshift is the controlling item, but if that's
        // not set (value<=-1), then we will run to the minimum of the TimeSliceRedshifts list.
        // If no TimeSlices are requested, then z=0.
        if (FinalRedshift>-1) return FinalRedshift;
        double minz = 1e10;
        for (int i=0; i<MAX_TIMESLICE_REDSHIFTS; i++)
            if (TimeSliceRedshifts[i] <= -1) break; 
            else if (TimeSliceRedshifts[i]<minz) minz = TimeSliceRedshifts[i];
        return minz;
    
        return 0.0;
    }

private:
    void ProcessStateDirectories();
    void CountTimeSlices();
};

// Convert a whole string to lower case, in place.
void strlower(char* str){
    for ( ; *str; ++str)
        *str = tolower(*str);
}

void Parameters::CountTimeSlices(){
    nTimeSlice = -1; 
    for (int i = 0; i < MAX_TIMESLICE_REDSHIFTS; i++){
        if (TimeSliceRedshifts[i] <= -1) {
            nTimeSlice = i; 
            break; 
        }
    }   

    nTimeSliceSubsample = -1; 
    for (int i = 0; i < MAX_TIMESLICE_REDSHIFTS; i++){
        if (TimeSliceRedshifts_Subsample[i] <= -1) {
            nTimeSliceSubsample = i; 
            break; 
        }
    } 

    nTimeSliceL1 = -1; 
    for (int i = 0; i < MAX_L1OUTPUT_REDSHIFTS; i++){
        if (L1OutputRedshifts[i] <= -1) {
            nTimeSliceL1 = i; 
            break; 
        }
    } 

    if (not (nTimeSlice > 0 || nTimeSliceSubsample > 0 || nTimeSliceL1 > 0 || OutputEveryStep))
        printf("Warning! No output requested. Are you sure you want this?\n");  //must request at least one kind of output. 
}
void Parameters::ProcessStateDirectories(){
    strlower(StateIOMode);
    strlower(Conv_IOMode);

    // Set read dir and write dir from working dir if they were not given
    if (strcmp(WorkingDirectory,STRUNDEF) != 0){
        if(strcmp(ReadStateDirectory,STRUNDEF) == 0){
            int ret = snprintf(ReadStateDirectory, 1024, "%s/read",WorkingDirectory);
            assert(ret >= 0 && ret < 1024);
        }
        if(strcmp(WriteStateDirectory,STRUNDEF) == 0){
            int ret = snprintf(WriteStateDirectory, 1024, "%s/write",WorkingDirectory);
            assert(ret >= 0 && ret < 1024);
            if(strcmp(StateIOMode, "overwrite") == 0) {  // later, we will set WriteState.OverwriteState
                strcpy(WriteStateDirectory, ReadStateDirectory);
            }
        }
    }

    // Set read dir and write dir from local working dir if they were not given
    if (strcmp(LocalWorkingDirectory,STRUNDEF) != 0){
        if(strcmp(LocalReadStateDirectory,STRUNDEF) == 0){
            int ret = snprintf(LocalReadStateDirectory, 1024, "%s/read",LocalWorkingDirectory);
            assert(ret >= 0 && ret < 1024);
        }
        if(strcmp(LocalWriteStateDirectory,STRUNDEF) == 0){
            int ret = snprintf(LocalWriteStateDirectory, 1024, "%s/write",LocalWorkingDirectory);
            assert(ret >= 0 && ret < 1024);
            if(strcmp(StateIOMode, "overwrite") == 0) {  // later, we will set WriteState.OverwriteState
                strcpy(LocalWriteStateDirectory, LocalReadStateDirectory);
            }
        }
    } else {
        // LocalWorkingDirectory not given; copy the global values (unless the local values were explicitly given)
        if(strcmp(LocalReadStateDirectory,STRUNDEF) == 0){
            strcpy(LocalReadStateDirectory, ReadStateDirectory);
        }
        if(strcmp(LocalWriteStateDirectory,STRUNDEF) == 0){
            strcpy(LocalWriteStateDirectory, WriteStateDirectory);
        }
    }
}


void Parameters::ReadParameters(char *parameterfile, int icflag) {
    hs = new HeaderStream(parameterfile);
    ReadHeader(*hs);
    hs->Close();
    CountTimeSlices();
    ProcessStateDirectories();
    if(!icflag) ValidateParameters();
}

void Parameters::ValidateParameters(void) {
    // Warning: Can't use STDLOG(), QUIT(), or assertf() calls in here:
    // We haven't opened the stdlog file yet!

    if(MaxPID == 0){
        MaxPID = np;
    }

    if(np<0) {
        fprintf(stderr,
            "[ERROR] np = %lld must be greater than zero!\n", np);
        assert(1==0);
    }


    if(cpd<0) {
        fprintf(stderr,
            "[ERROR] cpd = %d must be greater than zero!\n", cpd);
        assert(1==0);
    }

    if( !( (DirectNewtonRaphson == 1) ||  (DirectNewtonRaphson == 0) ) ) {
        fprintf(stderr,"DirectNewtonRapson must be 0 or 1\n");
        assert(1==0);
    }

    if(cpd%2==0) {
        fprintf(stderr,
            "[ERROR] cpd = %d  must be odd!\n", cpd);
        assert(1==0);
    }

    if(NearFieldRadius<=0) {
        fprintf(stderr, "[ERROR] NearFieldRadius = %d  must be greater than 0\n",
            NearFieldRadius );
        assert(1==0);
    }

    if(NearFieldRadius>(cpd-1)/2) {
        fprintf(stderr,
            "[ERROR] NearFieldRadius = %d must be less than (cpd-1)/2 = %d\n",
                    NearFieldRadius, (cpd-1)/2  );
        assert(1==0);
    }

    if(GroupRadius<=0) {
        fprintf(stderr, "[ERROR] GroupRadius = %d  must be greater than 0\n",
            GroupRadius );
        assert(1==0);
    }

    if(order>16 || order < 2 ) {
        fprintf(stderr,
            "[ERROR] order = %d must be less than or equal to 16 and greater than 1\n",
            order);
        assert(1==0);
    }

    if(MAXRAMMB<0) {
        fprintf(stderr,
            "[ERROR] MAXRAMMB = %d must be greater than 0 \n",
                MAXRAMMB);
        assert(1==0);
    }

    if(ConvolutionCacheSizeMB<=0) {
        fprintf(stderr,
            "[ERROR] ConvolutionCacheSizeMB = %f must be greater than 0 \n",
                ConvolutionCacheSizeMB);
        assert(1==0);
    }

    if(ConvolutionL1CacheSizeMB<=0) {
        fprintf(stderr,
            "[ERROR] ConvolutionL1CacheSizeMB = %f must be greater than 0 \n",
                ConvolutionL1CacheSizeMB);
        assert(1==0);
    }


    if( !((DerivativeExpansionRadius >= 1) && (DerivativeExpansionRadius <= 8)) &&
        (DerivativeExpansionRadius!=16) &&
        (DerivativeExpansionRadius!=32) ) {

        fprintf(stderr,
            "[ERROR] DerivativeExpansionRadius = %d has to be 1-8, 16, or 32\n",
                DerivativeExpansionRadius);
        assert(1==0);
    }

    if( (SofteningLength < 0) || (SofteningLength>BoxSize) ) {
        fprintf(stderr,
            "[ERROR] SofteningLength = %e has to be in [0,BoxSize)\n",
                SofteningLength);
        assert(1==0);
    }

    if (NumSlabsInsertList<0.0 || NumSlabsInsertList>cpd) {
        fprintf(stderr,
            "[ERROR] NumslabsInsertList = %e must be in range [0..CPD]\n",
                NumSlabsInsertList);
        assert(1==0);
    }

    if (NumSlabsInsertListIC<0.0 || NumSlabsInsertListIC>cpd) {
        fprintf(stderr,
            "[ERROR] NumslabsInsertListIC = %e must be in range [0..CPD]\n",
                NumSlabsInsertListIC);
        assert(1==0);
    }

    if (nTimeSlice<0) {
        fprintf(stderr,"nTimeSlice must be >=0\n");
        assert(1==0);
    }
    
    if(Omega_Smmoth<0.0 || Omega_Smooth>Omega_M){
        fprintf(stderr,"Must have 0<=Omega_Smooth<Omega_M, but told Omega_Smooth = %g\n", Omega_Smooth);
        assert(1==0);
    }
    
    if(abs(Omega_M + Omega_DE + Omega_K - 1.) > 1e-6){
        fprintf(stderr,"Omega_M + Omega_DE + Omega_K must equal 1, but is %g\n", Omega_M + Omega_DE + Omega_K);
        assert(1==0);
    }

    // Illegal ICFormat's will crash in loadIC.cpp.
    // But don't let the IC step run if it's going to crash when we try to do LPT

    // If invoking LPT, must have displacement-oriented format
    if(LagrangianPTOrder > 1){
        if(!(strcasecmp(ICFormat, "Zeldovich") == 0 ||
            strcasecmp(ICFormat, "RVZel") == 0 ||
            strcasecmp(ICFormat, "RVdoubleZel") == 0)) {
            fprintf(stderr, "ICFormat %s is not displacement-oriented and LagrangianPTOrder = %d\n",
                    ICFormat, LagrangianPTOrder);
        }
    }

    /*
    ExpandPathName(DerivativesDirectory);
    ExpandPathName(ReadStateDirectory);
    ExpandPathName(WriteStateDirectory);
    ExpandPathName(OutputDirectory);
    ExpandPathName(MultipoleDirectory);
    ExpandPathName(LogDirectory);
    ExpandPathName(InitialConditionsDirectory);

    CheckDirectoryExists(DerivativesDirectory);
    CheckDirectoryExists(ReadStateDirectory);
    CheckDirectoryExists(WriteStateDirectory);
    CheckDirectoryExists(OutputDirectory);
    CheckDirectoryExists(LogDirectory);

    CheckDirectoryExists(InitialConditionsDirectory);
    */

    if (ForceOutputDebug && StoreForces) {
            fprintf(stderr,"ForcesOutputDebug and StoreForces both set. This is not supported.\n");
        assert(1==0);
    }

    if (LogVerbosity<0) {
        fprintf(stderr,"LogVerbosity must be >=0\n");
        assert(1==0);
    }

    assert(nIODirs < MAX_IODIRS);
    for (int i = 0; i < nIODirs; i++)
        assert(IODirThreads[i] >= 1);

    assertf(
        strcmp(StateIOMode, "normal") == 0 ||
        strcmp(StateIOMode, "overwrite") == 0 ||
        strcmp(StateIOMode, "slosh") == 0 ||
        strcmp(StateIOMode, "stripe") == 0,
        "StateIOMode = \"%s\" must be one of normal, overwrite, slosh, stripe.",
        StateIOMode
        );

    assertf(
        strcmp(Conv_IOMode, "normal") == 0 ||
        strcmp(Conv_IOMode, "overwrite") == 0 ||
        strcmp(Conv_IOMode, "slosh") == 0 ||
        strcmp(Conv_IOMode, "stripe") == 0,
        "Conv_IOMode = \"%s\" must be one of normal, overwrite, slosh, stripe.",
        Conv_IOMode
        );
}
Parameters P;

#endif
