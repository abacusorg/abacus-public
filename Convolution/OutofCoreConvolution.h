//#include "read_dio.h"
//#include "write_dio.h"

// used by both serial and parallel io
ReadDirect *RD_RDD;
ReadDirect *RD_RDM;
WriteDirect *WD_WDT;

typedef struct { 
    double ReadDerivatives;
    double ReadMultipoles;
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

    uint64_t ReadDerivativesBytes, ReadMultipolesBytes, WriteTaylorBytes;
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
    int runtime_DiskBufferSizeKB;
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
    
    int delete_multipoles_after_read;
    
    uint64_t blocksize, zwidth, rml, CompressedMultipoleLengthXY;
    int io_cores[MAX_IO_THREADS];
    int niothreads;

    int ProfilingMode;

    int StripeConvState;  // in analogy with WriteState.StripeConvState

    void MultipoleDirectory(int slab, char * const fn){
        // We elsewhere generically support N threads, but here is where we assume 2
        if(StripeConvState && slab % 2 == 1)
            sprintf(fn, "%s/%s_%04d", runtime_MultipoleDirectory2, runtime_MultipolePrefix, slab);
        else
            sprintf(fn, "%s/%s_%04d", runtime_MultipoleDirectory, runtime_MultipolePrefix, slab);
    }

    void TaylorDirectory(int slab, char * const fn){
        // We elsewhere generically support N threads, but here is where we assume 2
        if(StripeConvState && slab % 2 == 1)
            sprintf(fn, "%s/%s_%04d", runtime_TaylorDirectory2, runtime_TaylorPrefix, slab);
        else
            sprintf(fn, "%s/%s_%04d", runtime_TaylorDirectory, runtime_TaylorPrefix, slab);
    }
};

#include "block_io_utils.cpp"

#ifdef CONVIOTHREADED
#include "ConvolutionIOThread.cpp"
#endif

class OutofCoreConvolution { 
public: 
    
    OutofCoreConvolution(void) {
        memset(&CS, 0, sizeof(ConvolutionStatistics));
    } 
    ~OutofCoreConvolution(void) { } 

    ConvolutionParameters CP;
    void Convolve( ConvolutionParameters CP );

    uint64_t blocksize, zwidth;

    ConvolutionStatistics CS; 

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

    Complex *DiskBuffer;
    DFLOAT **CompressedDerivatives;
    MTCOMPLEX **PlaneBuffer;
    Block *CurrentBlock;
    
    double invcpd3;
};
