// Copyright 2012-2025 The Abacus Developers
// SPDX-License-Identifier: GPL-3.0-or-later

//#include "read_dio.h"
//#include "write_dio.h"

// used by both serial and parallel io
ReadDirect *RD_RDD;
ReadDirect *RD_RDM;
WriteDirect *WD_WDT;

typedef struct { 
    double ReadDerivatives = 0.;
    double ReadMultipoles = 0.;
	double TransposeBuffering = 0.;
	double TransposeAlltoAllv = 0.;
    double WriteTaylor = 0.;
    double ForwardZFFTMultipoles = 0.;
    double InverseZFFTTaylor = 0.;
    double ConvolutionArithmetic = 0.;
    double ArraySwizzle = 0.;
	
		
#ifdef PARALLEL
	double Constructor = 0.0;
	double AllocMT = 0.0;
	double AllocDerivs = 0.0;
	double QueueTaylors = 0.0; 
	double FFTPlanning = 0.0; 
	double Destructor = 0.0; 
	double ThreadCleanUp = 0.0; 
	

	double MmapMT = 0.0; 
	double MunmapMT = 0.0; 
	double MmapDerivs = 0.0; 
	double MunmapDerivs = 0.0;
	
	
#endif	
	
#ifdef CONVIOTHREADED
    double WaitForIO = 0;
#endif
    
    double ConvolveWallClock = 0.0;
    double Discrepency = 0.;

    uint64_t ReadDerivativesBytes = 0, ReadMultipolesBytes = 0, TransposeBufferingBytes = 0, TransposeAlltoAllvBytes = 0, WriteTaylorBytes = 0;
    uint64_t ops = 0;
    uint64_t totalMemoryAllocated = 0;

    float runtime_ConvolutionCacheSizeMB = 0.;
    float runtime_ConvolutionL1CacheSizeMB = 0.;
    int ComputeCores = 0;
} ConvolutionStatistics;

class ConvolutionParameters{
public:
    int runtime_cpd = 0;
    int runtime_order = 0;
    int runtime_NearFieldRadius = 0;
    int runtime_DerivativeExpansionRadius = 0;

    int runtime_AllowDIO = 0;
    int runtime_DIOBufferSizeKB = 0;
    float runtime_ConvolutionCacheSizeMB = 0;
    float runtime_ConvolutionL1CacheSizeMB = 0;
    float runtime_MaxConvolutionRAMMB = 0;

    // These directories should be accessed through the functions below
    fs::path runtime_MultipoleDirectory;
    fs::path runtime_TaylorDirectory;
    fs::path runtime_DerivativesDirectory;
    fs::path runtime_MultipoleDirectory2;
    fs::path runtime_TaylorDirectory2;
    fs::path runtime_MultipolePrefix;
    fs::path runtime_TaylorPrefix;
    
    uint64_t blocksize = 0, zwidth = 0, rml = 0, CompressedMultipoleLengthXY = 0;

	uint64_t z_slabs_per_node = 0;

    int *first_slabs_all = NULL;
    int *total_slabs_all = NULL;
	
	std::vector<int> io_cores;
    int niothreads = 0;

    int ProfilingMode = 0;

    int StripeConvState = 0;  // in analogy with WriteState.StripeConvState
    int OverwriteConvState = 0;  // in analogy with WriteState.OverwriteConvState

    ConvolutionParameters(int mpi_nranks){
        first_slabs_all = new int[mpi_nranks];
        total_slabs_all = new int[mpi_nranks];
    }

    ~ConvolutionParameters(){
        delete[] first_slabs_all;
        delete[] total_slabs_all;
    }

    fs::path MultipoleFN(int slab){
        // We elsewhere generically support N threads, but here is where we assume 2
        if(StripeConvState && slab % 2 == 1)
            return runtime_MultipoleDirectory2 / fmt::format("{}_{:04d}", runtime_MultipolePrefix, slab);
        else
            return runtime_MultipoleDirectory / fmt::format("{}_{:04d}", runtime_MultipolePrefix, slab);
    }

    fs::path TaylorFN(int slab){
        // We elsewhere generically support N threads, but here is where we assume 2
        if(StripeConvState && slab % 2 == 1)
            return runtime_TaylorDirectory2 / fmt::format("{}_{:04d}", runtime_TaylorPrefix, slab);
        else
            return runtime_TaylorDirectory / fmt::format("{}_{:04d}", runtime_TaylorPrefix, slab);
    }

    // For zwidth purposes, we'll need to know if the MT are on ramdisk
    int is_ramdisk(){
        // For simplicity, multipoles and taylors must be either both or neither on ramdisk
        int mramdisk = is_path_on_ramdisk(runtime_MultipoleDirectory);
        assert(mramdisk == is_path_on_ramdisk(runtime_TaylorDirectory));

        return mramdisk;
    }
};
