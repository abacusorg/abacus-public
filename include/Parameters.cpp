// Copyright 2012-2025 The Abacus Developers
// SPDX-License-Identifier: GPL-3.0-or-later

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
    
    std::string SimName; //What to call this run
    uint64 np;
    int cpd;
    int order;

    int NearFieldRadius;    // Radius of cells in the near-field
    double SofteningLength; // Softening length in the same units as BoxSize

    int  DerivativeExpansionRadius;
    int  MAXRAMMB;
    double MemPerGPUBufferGB;  // Max allocation for each GPU buffer, of which there will be NGPU * DirectBPD
    double  ConvolutionCacheSizeMB; // Set to manually override the detected cache size; this is for L3
    double  ConvolutionL1CacheSizeMB; // Set to manually override the detected cache size
    int AllowDirectIO;        // ==1 for a normal disk, ==0 for a ramdisk or sometimes network file system
    int ForceBlockingIO;   // ==1 if you want to force all IO to be blocking.
    std::string StateIOMode;  //  "normal", "slosh", "overwrite", "stripe"
    std::string Conv_IOMode;  //  "normal", "slosh", "overwrite", "stripe"
    
    int OMP_NUM_THREADS;  // Number of OpenMP threads.  0 does not modify the system value (usually OMP_NUM_THREADS, or all threads).
                        // Negative values use that many fewer than the max.

    int  DirectNewtonRaphson;  // 0 or 1 

    fs::path DerivativesDirectory;

    fs::path InitialConditionsDirectory;   // The initial condition file name
    std::string ICFormat;                // The format of the IC files
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

    fs::path ReadStateDirectory;  // Where the input State lives
    fs::path WriteStateDirectory; // Where the output State lives

    // The node-local directories in the parallel verison
    // These will be equal to the global directories if not given
    fs::path LocalWorkingDirectory;        // node-local working directory
    fs::path LocalReadStateDirectory;
    fs::path LocalWriteStateDirectory;

    fs::path MultipoleDirectory;
    fs::path TaylorDirectory;
    fs::path MultipoleDirectory2;  // for splitting even/odd multipoles
    fs::path TaylorDirectory2;
    fs::path WorkingDirectory;        // If Read/Write/Past not specified, where to put the states
    fs::path LogDirectory;
    fs::path OutputDirectory;     // Where the outputs go
    fs::path GroupDirectory; //Where the Group files go
    fs::path BackupDirectory; // The directory from which to restore backups (the Python code also writes backups here)
    
    int OutputEveryStep; //Force timeslices to be output every step if 1
    std::string OutputFormat;                // The format of the Output files
    int  OmitOutputHeader;                // =1 if you want to skip the ascii header

    double FinalRedshift;        // When to stop.  This will override TimeSliceRedshifts.
    std::vector<double> TimeSliceRedshifts;
    std::vector<double> TimeSliceRedshifts_Subsample;

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

    // Could have group finding or coevolution set instructions

    int  OutputFullLightCones;

    std::vector<double> LCorigins;  // Same units as BoxSize

    fs::path LCDirectory;
    int LCCheckAcrossWrap;  // If 1, check for light cone particles that cross the periodic wrap
    int LCBoxRepeats;  // The number of times to repeat the box in +/- x, y, z directions (0 = normal light cone)

    int PowerSpectrumStepInterval;
    int PowerSpectrumN1d; //1D number of bins to use in the powerspectrum

    int LogVerbosity;   // If 0, production-level log; higher numbers are more verbose
    int StoreForces; // If 1, store the accelerations. 2 = separate near and far forces. 3 = only on output steps
    int ForceOutputDebug; // If 1, output near and far forces seperately. 

    int MakeGlass; // Reverse the sign of the acceleration

    int ForceCPU; //IF 1, force directs to be executed exclusively on cpu even if cuda is available
    int GPUMinCellSinks;// If AVX directs are compiled, cells with less than this many particles go to cpu
    int ProfilingMode;//If 1, enable profiling mode, i.e. delete the write-state after creating it to run repeatedly on same dat
    
    std::vector<int> IOCores;  // The cores that the IO threads will be bound to.  -1 means don't bind
    std::vector<fs::path> IODirs;
    std::vector<int> IODirThreads;
    
    int Conv_OMP_NUM_THREADS;
    std::vector<int> Conv_IOCores;
    int Conv_zwidth;
    
    // TODO: this scheme doesn't account for more complicated NUMA architectures
    std::vector<int> GPUThreadCoreStart;  // The core on which to start placing GPU device threads.
    int NGPUThreadCores;  // The number of free cores on which to place GPU device threads.
    std::vector<int> GPUQueueAssignments;  // The work queue assignments
    int DirectBPD;

    double DensityKernelRad;  // The kernel Radius for the density computation, specified in units of the interparticle spacing.  0 will default to FoFLinkingLength[0]
    double L0DensityThreshold;  // The kernel density required for a particle to be eligible to be in a L0 group; specified in units of the cosmic mean density. This is ignored (uses 0) if DensityKernelRad==0.  Value = 0 triggers code to make a particle eligible if it has any non-self neighbor within DensityKernelRad

    int AllowGroupFinding;
    
    std::vector<double> FoFLinkingLength = std::vector<double>(3); //Linking lengths for level 0,1,2 groupfinding in fractional interparticle spacing 
    std::vector<double> SODensity = std::vector<double>(2);  // Overdensities for SO groupfinding level 1 and 2
    
    int MinL1HaloNP; // minimum L1 halo size to output
	double L1Output_dlna;  // minimum delta ln(a) between L1 halo outputs
    double SO_RocheCoeff; 
    double SO_alpha_eligible;   // New centers must be outside alpha*R_Delta
    double SO_NPForMinDensity; 
    int SO_EvolvingThreshold; //allow evolving (redshift-dependent) density threshold 

    std::vector<double> L1OutputRedshifts;
    int OutputAllHaloParticles;  // ==0 normally, to output only taggable L1 particles.  If non-zero, output all particles

    long long int MaxPID;  // Maximum PID to expect.  A PID equal or larger than this indicates corruption of some sort.  0 means NP; -1 means don't check.

    int ProperSoftening;  // Keep the softening length fixed in proper coordinates.  SofteningLength is specified at z=0.
    double SofteningMax;  // The maximum comoving softening to allow when using ProperSoftening

    int NoChecksum;  // Explicitly disable output checksumming

    int UsePinnedGPUMemory;  // Whether to pin the CPU-side GPU staging buffers

    double MPICallRateLimit_ms;  // Enforce a delay between MPI_Test calls

    int UseMunmapThread;  // dedicated munmap() thread in ArenaAllocator

    int MunmapThreadCore;  // Core binding for disposal thread

    int NumZRanks;  // Number of ranks in the Z-dimension in the 2D domain decomposition

    int InsertListGapElems;  // Size of the MAL gaps in the Insert List, in number of ilstruct elems

    int ForceAuxDensity;  // Use densities from the aux in group finding, even if we have acc densities

    uint64_t LCHealpixNside;  // Healpix lightcone resolution

    int LCHealpixOutputSparseMap;  // Output healpix "structs" that can later be binned into a map

    int LCOutputRVPID;  // Whether to output the RV & PID on the lightcone

    int OutputRVAtL1RedshiftsA;  // Whether to output subsample A RVs at L1 output redshifts

    int ReleaseFreeMemory;  // Whether to intermittently release free memory to the system

    // Return the L{tier} size in MiB
    double getCacheSize(int tier){
        int cache_size = 0;
        FILE *fp = NULL;
        fs::path fname;

        // find the last-level cache
        for(int i = 0; i<=tier; i++){
            fname = fmt::format("/sys/devices/system/cpu/cpu0/cache/index{:d}/size", i);
            fp = fopen(fname.c_str(), "r");
            if(fp == NULL)
                break;
            int nscan = fscanf(fp, "%dK", &cache_size);  // cache size in KiB
            fclose(fp);
            assertf(nscan == 1, "Unexpected cache size file format (\"{}\")\n", fname);
        }
        return cache_size/1024.0;    // to MiB
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
        installscalar("SofteningLength", SofteningLength, MUST_DEFINE); // Softening length in the same units as BoxSize
        installscalar("DerivativeExpansionRadius", DerivativeExpansionRadius,MUST_DEFINE);
        MAXRAMMB = getRAMSize();
        installscalar("MAXRAMMB", MAXRAMMB, DONT_CARE);
        MemPerGPUBufferGB = 0.0;  // auto
        installscalar("MemPerGPUBufferGB", MemPerGPUBufferGB, DONT_CARE);
        ConvolutionCacheSizeMB = getCacheSize(4);
        installscalar("ConvolutionCacheSizeMB", ConvolutionCacheSizeMB, DONT_CARE);
        ConvolutionL1CacheSizeMB = getCacheSize(1);
        installscalar("ConvolutionL1CacheSizeMB", ConvolutionL1CacheSizeMB, DONT_CARE);
        AllowDirectIO = 1;
        installscalar("AllowDirectIO",AllowDirectIO,DONT_CARE);
        ForceBlockingIO = 0;
        installscalar("ForceBlockingIO",ForceBlockingIO,DONT_CARE);

        StateIOMode = "normal";
        installscalar("StateIOMode", StateIOMode, DONT_CARE);
        Conv_IOMode = "normal";
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

        NumSlabsInsertList = 0.;  // auto
        installscalar("NumSlabsInsertList",NumSlabsInsertList,DONT_CARE);   
        NumSlabsInsertListIC = 0.;
        installscalar("NumSlabsInsertListIC",NumSlabsInsertListIC,DONT_CARE);   

        ReadStateDirectory = STRUNDEF;
        WriteStateDirectory = STRUNDEF;
        WorkingDirectory = STRUNDEF;
        MultipoleDirectory = STRUNDEF;
        TaylorDirectory = STRUNDEF;
        MultipoleDirectory2 = STRUNDEF;
        TaylorDirectory2 = STRUNDEF;
        BackupDirectory = STRUNDEF;
        LocalWorkingDirectory = STRUNDEF;
        LocalReadStateDirectory = STRUNDEF;
        LocalWriteStateDirectory = STRUNDEF;

        installscalar("ReadStateDirectory",ReadStateDirectory,DONT_CARE);  // Where the input State lives
        installscalar("WriteStateDirectory",WriteStateDirectory,DONT_CARE); // Where the output State lives
        installscalar("WorkingDirectory",WorkingDirectory,DONT_CARE);
        installscalar("MultipoleDirectory",MultipoleDirectory,DONT_CARE);
        installscalar("TaylorDirectory",TaylorDirectory,DONT_CARE);
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

        installscalar("LCDirectory",LCDirectory,MUST_DEFINE); //Where the lightcones go. Generally will be the same as the Output directory
        OutputFullLightCones = 0;
        installscalar("OutputFullLightCones",OutputFullLightCones,DONT_CARE); //if not set, we assume 0

        installvector("LCOrigins",LCorigins,DONT_CARE);
        LCCheckAcrossWrap = 0;
        installscalar("LCCheckAcrossWrap",LCCheckAcrossWrap,DONT_CARE);
        LCBoxRepeats = 0;
        installscalar("LCBoxRepeats",LCBoxRepeats,DONT_CARE);

        FinalRedshift = -2.0;        // If <-1, then we will cascade back to the minimum of the TimeSliceRedshifts list
        installscalar("FinalRedshift",FinalRedshift,DONT_CARE);
		
        installvector("TimeSliceRedshifts", TimeSliceRedshifts, DONT_CARE);
        installvector("TimeSliceRedshifts_Subsample", TimeSliceRedshifts_Subsample, DONT_CARE);
        installvector("L1OutputRedshifts", L1OutputRedshifts, DONT_CARE);

        ParticleSubsampleA = 0.;
        ParticleSubsampleB = 0.;
        installscalar("ParticleSubsampleA", ParticleSubsampleA, DONT_CARE);
        installscalar("ParticleSubsampleB", ParticleSubsampleB, DONT_CARE); 
		
        OutputFormat = "RVdouble";
        // OutputFormat = "Packed";
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

        GroupRadius = 0;
        installscalar("GroupRadius",GroupRadius,DONT_CARE);        // Maximum size of a group, in units of cell sizes
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

        installvector("IOCores", IOCores, DONT_CARE);
        installvector("Conv_IOCores", Conv_IOCores, DONT_CARE);

        Conv_OMP_NUM_THREADS = 0;
        installscalar("Conv_OMP_NUM_THREADS", Conv_OMP_NUM_THREADS, DONT_CARE);

        Conv_zwidth = -1;
        installscalar("Conv_zwidth", Conv_zwidth, DONT_CARE);

        installvector("IODirs", IODirs, DONT_CARE);
        installvector("IODirThreads", IODirThreads, DONT_CARE);

        // If GPUThreadCoreStart is undefined, GPU threads will not be bound to cores
        installvector("GPUThreadCoreStart", GPUThreadCoreStart, DONT_CARE);
        installvector("GPUQueueAssignments", GPUQueueAssignments, DONT_CARE);
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
        installvector("FoFLinkingLength",FoFLinkingLength,DONT_CARE);
        SODensity[0] = 200.0;
        SODensity[1] = 800.0;
        installvector("SODensity",SODensity,DONT_CARE);
        MinL1HaloNP = 40;
        installscalar("MinL1HaloNP", MinL1HaloNP, DONT_CARE);
		L1Output_dlna = -1;
		installscalar("L1Output_dlna", L1Output_dlna, DONT_CARE);

        SO_RocheCoeff = 2.0; 
        installscalar("SO_RocheCoeff", SO_RocheCoeff, DONT_CARE);
        SO_alpha_eligible = 0.8; 
        installscalar("SO_alpha_eligible", SO_alpha_eligible, DONT_CARE);
        SO_NPForMinDensity = 35.0; 
        installscalar("SO_NPForMinDensity", SO_NPForMinDensity, DONT_CARE);
        SO_EvolvingThreshold = 1; 
        installscalar("SO_EvolvingThreshold", SO_EvolvingThreshold, DONT_CARE); 

        OutputAllHaloParticles = 0;
        installscalar("OutputAllHaloParticles", OutputAllHaloParticles, DONT_CARE);

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

        NumZRanks = 1;
        installscalar("NumZRanks", NumZRanks, DONT_CARE);

        InsertListGapElems = PAGE_SIZE/40;
        installscalar("InsertListGapElems", InsertListGapElems, DONT_CARE);

        ForceAuxDensity = 0;
        installscalar("ForceAuxDensity", ForceAuxDensity, DONT_CARE);

        LCHealpixNside = 16384;
        installscalar("LCHealpixNside", LCHealpixNside, DONT_CARE);

        LCOutputRVPID = 1;
        installscalar("LCOutputRVPID", LCOutputRVPID, DONT_CARE);

        LCHealpixOutputSparseMap = 0;
        installscalar("LCHealpixOutputSparseMap", LCHealpixOutputSparseMap, DONT_CARE);

        OutputRVAtL1RedshiftsA = 0;
        installscalar("OutputRVAtL1RedshiftsA", OutputRVAtL1RedshiftsA, DONT_CARE);

        ReleaseFreeMemory = 1;
        installscalar("ReleaseFreeMemory", ReleaseFreeMemory, DONT_CARE);
    }

    // We're going to keep the HeaderStream, so that we can output it later.
    HeaderStream *hs;
    ~Parameters(void) {
        delete hs;
    }
    std::string header() { 
        assert(hs!=NULL); assert(hs->buffer!=NULL);
        return hs->buffer;        // This is just a standard C-style string.
    }


    void ReadParameters(const fs::path &paramaterfile, int icflag);
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
        uint64 n = floor(ppd());
        if (n*n*n==np) return 1; else return 0;
    }

    double FinishingRedshift(bool *finalz_is_L1 = NULL) {
        // Return the redshift where the code should halt.
        // FinalRedshift is the controlling item, but if that's
        // not set (value<=-1), then we will the redshifts of
        // the requested outputs.
        // If no TimeSlices are requested, then z=0.
        if (FinalRedshift > -1) return FinalRedshift;

        double minz = 1e100;
        bool have_minz = false;

        for (const auto& redshift : TimeSliceRedshifts) {
            minz = std::min(minz, redshift);
            have_minz = true;
        }

        for (const auto& redshift : TimeSliceRedshifts_Subsample) {
            minz = std::min(minz, redshift);
            have_minz = true;
        }

        for (const auto& redshift : L1OutputRedshifts) {
            minz = std::min(minz, redshift);
            have_minz = true;
        }

        if (have_minz) {
            if(finalz_is_L1 != NULL) *finalz_is_L1 = contains(L1OutputRedshifts.begin(), L1OutputRedshifts.end(), minz);
            return minz;
        }
        return 0.0;
    }

    void ProcessStateDirectories();

private:
    void SortTimeSlices();
};

void strlower(std::string &str){
    for (size_t i = 0; i < str.size(); i++)
        str[i] = tolower(str[i]);
}

void Parameters::SortTimeSlices(){
    std::sort(TimeSliceRedshifts.begin(), TimeSliceRedshifts.end(), std::greater<double>());
    std::sort(TimeSliceRedshifts_Subsample.begin(), TimeSliceRedshifts_Subsample.end(), std::greater<double>());
    std::sort(L1OutputRedshifts.begin(), L1OutputRedshifts.end(), std::greater<double>());

    if(! (TimeSliceRedshifts.size() > 0 || TimeSliceRedshifts_Subsample.size() > 0 || L1OutputRedshifts.size() > 0 || OutputEveryStep || StoreForces || LCorigins.size()) )
        fmt::print("Warning! No output requested. Are you sure you want this?\n");
}

void Parameters::ProcessStateDirectories(){
    strlower(StateIOMode);
    strlower(Conv_IOMode);

    // Set read dir and write dir from working dir if they were not given
    if (WorkingDirectory != STRUNDEF){
        if(ReadStateDirectory == STRUNDEF){
            ReadStateDirectory = WorkingDirectory / "read";
        }
        if(WriteStateDirectory == STRUNDEF){
            WriteStateDirectory = WorkingDirectory / "write";
        }
    }

    // Set read dir and write dir from local working dir if they were not given
    if (LocalWorkingDirectory != STRUNDEF){
        // append the rank to the LocalWorkingDirectory
        // LocalWorkingDirectory should not have a trailing slash
        LocalWorkingDirectory += NodeString;
    } else {
        LocalWorkingDirectory = WorkingDirectory.string() + NodeString;
    }
    

    if(LocalReadStateDirectory == STRUNDEF){
        LocalReadStateDirectory = LocalWorkingDirectory / "read";
    }
    if(LocalWriteStateDirectory == STRUNDEF){
        LocalWriteStateDirectory = LocalWorkingDirectory / "write";
    }

    if (MultipoleDirectory == STRUNDEF){
        // not given
        MultipoleDirectory = LocalWorkingDirectory / "multipole";
    }

    if (TaylorDirectory == STRUNDEF){
        // not given
        TaylorDirectory = LocalWorkingDirectory / "taylor";
    }
}


void Parameters::ReadParameters(const fs::path &parameterfile, int icflag) {
    hs = new HeaderStream(parameterfile);
    ReadHeader(*hs);
    hs->Close();
    SortTimeSlices();
    // ProcessStateDirectories();  // this is now done in singlestep.cpp (after MPI)
    if(!icflag) ValidateParameters();
}

void Parameters::ValidateParameters(void) {
    // Warning: Can't use STDLOG(), QUIT(), or assertf() calls in here:
    // We haven't opened the stdlog file yet!

    if(MaxPID == 0){
        MaxPID = np;
    }

    if(cpd<0) {
        fmt::print(stderr,
            "[ERROR] cpd = {:d} must be greater than zero!\n", cpd);
        assert(1==0);
    }

    if( !( (DirectNewtonRaphson == 1) ||  (DirectNewtonRaphson == 0) ) ) {
        fmt::print(stderr,"DirectNewtonRapson must be 0 or 1\n");
        assert(1==0);
    }

    if(cpd%2==0) {
        fmt::print(stderr,
            "[ERROR] cpd = {:d}  must be odd!\n", cpd);
        assert(1==0);
    }

    if(NearFieldRadius<=0) {
        fmt::print(stderr, "[ERROR] NearFieldRadius = {:d}  must be greater than 0\n",
            NearFieldRadius );
        assert(1==0);
    }

    if(NearFieldRadius>(cpd-1)/2) {
        fmt::print(stderr,
            "[ERROR] NearFieldRadius = {:d} must be less than (cpd-1)/2 = {:d}\n",
                    NearFieldRadius, (cpd-1)/2  );
        assert(1==0);
    }

    if(AllowGroupFinding && GroupRadius<0) {
        fmt::print(stderr, "[ERROR] GroupRadius = {:d}  must be >= 0\n",
            GroupRadius );
        assert(1==0);
    }

    if(! (order <= 16 && order >= 0) ) {
        fmt::print(stderr,
            "[ERROR] order = {:d} must be less than or equal to 16 and greater than 1\n",
            order);
        assert(1==0);
    }

    if(MAXRAMMB<0) {
        fmt::print(stderr,
            "[ERROR] MAXRAMMB = {:d} must be greater than 0 \n",
                MAXRAMMB);
        assert(1==0);
    }

    if(MemPerGPUBufferGB<0) {
        fmt::print(stderr,
            "[ERROR] MemPerGPUBufferGB = {:f} must be greater than 0 \n",
                MemPerGPUBufferGB);
        assert(1==0);
    }

    if(ConvolutionCacheSizeMB<=0) {
        fmt::print(stderr,
            "[ERROR] ConvolutionCacheSizeMB = {:f} must be greater than 0 \n",
                ConvolutionCacheSizeMB);
        assert(1==0);
    }

    if(ConvolutionL1CacheSizeMB<=0) {
        fmt::print(stderr,
            "[ERROR] ConvolutionL1CacheSizeMB = {:f} must be greater than 0 \n",
                ConvolutionL1CacheSizeMB);
        assert(1==0);
    }


    if( !((DerivativeExpansionRadius >= 1) && (DerivativeExpansionRadius <= 8)) &&
        (DerivativeExpansionRadius!=16) &&
        (DerivativeExpansionRadius!=32) ) {

        fmt::print(stderr,
            "[ERROR] DerivativeExpansionRadius = {:d} has to be 1-8, 16, or 32\n",
                DerivativeExpansionRadius);
        assert(1==0);
    }

    if( (SofteningLength < 0) || (SofteningLength>BoxSize) ) {
        fmt::print(stderr,
            "[ERROR] SofteningLength = {:e} has to be in [0,BoxSize)\n",
                SofteningLength);
        assert(1==0);
    }

    if (NumSlabsInsertList<0.0) {
        fmt::print(stderr,
            "[ERROR] NumslabsInsertList = {:e} must be >= 0\n",
                NumSlabsInsertList);
        assert(1==0);
    }

    if (NumSlabsInsertListIC<0.0) {
        fmt::print(stderr,
            "[ERROR] NumslabsInsertListIC = {:e} must be >= 0\n",
                NumSlabsInsertListIC);
        assert(1==0);
    }
    
    if(Omega_Smooth<0.0 || Omega_Smooth>Omega_M){
        fmt::print(stderr,"Must have 0<=Omega_Smooth<Omega_M, but told Omega_Smooth = {:g}\n", Omega_Smooth);
        assert(1==0);
    }
    
    if(abs(Omega_M + Omega_DE + Omega_K - 1.) > 1e-6){
        fmt::print(stderr,"Omega_M + Omega_DE + Omega_K must equal 1, but is {:g}\n", Omega_M + Omega_DE + Omega_K);
        assert(1==0);
    }

    // Illegal ICFormat's will crash in loadIC.cpp.
    // But don't let the IC step run if it's going to crash when we try to do LPT

    // If invoking LPT, must have displacement-oriented format
    if(LagrangianPTOrder > 1){
        if(!((ICFormat == "Zeldovich") ||
            (ICFormat == "RVZel") ||
            (ICFormat == "RVdoubleZel"))) {
            fmt::print(stderr, "Warning! ICFormat {:s} is not displacement-oriented and LagrangianPTOrder = {:d}. Forcing LagrangianPTOrder = 0.\n",
                    ICFormat, LagrangianPTOrder);
            LagrangianPTOrder = 0;
        }
    }

    /*
    DerivativesDirectory = fs::absolute(DerivativesDirectory);
    ReadStateDirectory = fs::absolute(ReadStateDirectory);
    WriteStateDirectory = fs::absolute(WriteStateDirectory);
    OutputDirectory = fs::absolute(OutputDirectory);
    MultipoleDirectory = fs::absolute(MultipoleDirectory);
    LogDirectory = fs::absolute(LogDirectory);
    InitialConditionsDirectory = fs::absolute(InitialConditionsDirectory);

    assert(fs::is_directory(DerivativesDirectory));
    assert(fs::is_directory(ReadStateDirectory));
    assert(fs::is_directory(WriteStateDirectory));
    assert(fs::is_directory(OutputDirectory));
    assert(fs::is_directory(LogDirectory));
    assert(fs::is_directory(InitialConditionsDirectory));
    */

    if (ForceOutputDebug) {
            StoreForces = 2;  // Output near and far separately
    }

    assertf(
        (StateIOMode == "normal") ||
        (StateIOMode == "overwrite") ||
        (StateIOMode == "slosh") ||
        (StateIOMode == "stripe"),
        "StateIOMode = \"{:s}\" must be one of normal, overwrite, slosh, stripe.",
        StateIOMode
        );

    assertf(
        (Conv_IOMode == "normal") ||
        (Conv_IOMode == "overwrite") ||
        (Conv_IOMode == "slosh") ||
        (Conv_IOMode == "stripe"),
        "Conv_IOMode = \"{:s}\" must be one of normal, overwrite, slosh, stripe.",
        Conv_IOMode
        );

#ifdef PARALLEL
    if(NumZRanks < 1){
        fmt::print(stderr,"NumZRanks={:d} must be >= 1!\n", NumZRanks);
        assert(1==0);
    }
#else
    if(NumZRanks > 1){
        fmt::print(stderr,"Warning: NumZRanks={:d} will have no effect because code is not compiled for parallel. Forcing NumZRanks=1.\n", NumZRanks);
        NumZRanks = 1;
    }
#endif
}

Parameters P;

#endif
