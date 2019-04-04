//#include "read_dio.h"
//#include "write_dio.h"

// used by both serial and parallel io
ReadDirect *RD_RDD;
ReadDirect *RD_RDM;
WriteDirect *WD_WDT;

typedef struct { 
    double ReadDerivatives;
    double ReadMultipoles;
	double TransposeBuffering;
	double TransposeAlltoAllv;
    double WriteTaylor;
    double ForwardZFFTMultipoles;
    double InverseZFFTTaylor;
    double ConvolutionArithmetic;
    double ArraySwizzle;
	
		
	
	
#ifdef CONVIOTHREADED
    double WaitForIO;
#endif
    
    double ConvolveWallClock;
    double Discrepency;

    uint64_t ReadDerivativesBytes, ReadMultipolesBytes, TransposeBufferingBytes, TransposeAlltoAllvBytes, WriteTaylorBytes;
    uint64_t ops;
    uint64_t totalMemoryAllocated;

    int runtime_ConvolutionCacheSizeMB;
    int ComputeCores;
} ConvolutionStatistics;

class ConvolutionParameters{
public:
    int runtime_cpd;
    int runtime_order;
    int runtime_NearFieldRadius;
    int runtime_DerivativeExpansionRadius;

    int runtime_IsRamDisk;
    int runtime_DIOBufferSizeKB;
    int runtime_ConvolutionCacheSizeMB;
    int runtime_MaxConvolutionRAMMB;

    // These directories should be accessed through the functions below
    char runtime_MultipoleDirectory[1024];
    char runtime_TaylorDirectory[1024];
    char runtime_DerivativesDirectory[1024];
    char runtime_MultipoleDirectory2[1024];
    char runtime_TaylorDirectory2[1024];
    char runtime_MultipolePrefix[1024];
    char runtime_TaylorPrefix[1024];
    
    uint64_t blocksize, zwidth, rml, CompressedMultipoleLengthXY;

	uint64_t z_slabs_per_node;

    int *first_slabs_all;
    int *total_slabs_all;
	
	int io_cores[MAX_IO_THREADS];
    int niothreads;

    int ProfilingMode;

    int StripeConvState;  // in analogy with WriteState.StripeConvState
    int OverwriteConvState;  // in analogy with WriteState.OverwriteConvState

    ConvolutionParameters(int mpi_nranks){
        first_slabs_all = new int[mpi_nranks];
        total_slabs_all = new int[mpi_nranks];
    }

    ~ConvolutionParameters(){
        delete[] first_slabs_all;
        delete[] total_slabs_all;
    }

    void MultipoleFN(int slab, char * const fn){
        // We elsewhere generically support N threads, but here is where we assume 2
        if(StripeConvState && slab % 2 == 1)
            sprintf(fn, "%s/%s_%04d", runtime_MultipoleDirectory2, runtime_MultipolePrefix, slab);
        else
            sprintf(fn, "%s/%s_%04d", runtime_MultipoleDirectory, runtime_MultipolePrefix, slab);
    }

    void TaylorFN(int slab, char * const fn){
        // We elsewhere generically support N threads, but here is where we assume 2
        if(StripeConvState && slab % 2 == 1)
            sprintf(fn, "%s/%s_%04d", runtime_TaylorDirectory2, runtime_TaylorPrefix, slab);
        else
            sprintf(fn, "%s/%s_%04d", runtime_TaylorDirectory, runtime_TaylorPrefix, slab);
    }

    // For zwidth purposes, we'll need to know if the MT are on ramdisk
    int is_ramdisk(){
        // For simplicity, multipoles and taylors must be either both or neither on ramdisk
        int mramdisk = is_path_on_ramdisk(runtime_MultipoleDirectory);
        assert(mramdisk == is_path_on_ramdisk(runtime_TaylorDirectory));

        return mramdisk;
    }
};

#include "block_io_utils.cpp"

#ifdef CONVIOTHREADED
#include "ConvolutionIOThread.cpp"
#endif

class OutofCoreConvolution { 
public: 

    ConvolutionStatistics CS; 
    ConvolutionParameters CP;
    
    OutofCoreConvolution(ConvolutionParameters &_CP);
    ~OutofCoreConvolution(void) { }     
    void Convolve();

    uint64_t blocksize, zwidth;
	
#ifdef PARALLEL
	uint64_t z_slabs_per_node; 
#endif


private:

    STimer ForwardZFFTMultipoles;
    STimer InverseZFFTTaylor;

    STimer ConvolutionArithmetic;
    STimer ArraySwizzle;

    STimer ConvolveWallClock;

#ifdef CONVIOTHREADED
    ConvIOThread **iothreads;
    STimer WaitForIO;
#endif

    uint64_t cpd,order,rml,CompressedMultipoleLengthXY;

    void BlockConvolve(void);
    void WriteDiskTaylor(int z);
    void ReadDiskMultipolesAndDerivs(int z);
    void SwizzleMultipoles(int z);
    void SwizzleTaylors(int z);

    void RenameMultipolesToTaylors();

    Complex *PlaneBuffer;
    DFLOAT **CompressedDerivatives;
    MTCOMPLEX **DiskBuffer;
    Block *CurrentBlock;
    
    double invcpd3;
};
