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

typedef struct {
    int runtime_cpd;
    int runtime_order;
    int runtime_NearFieldRadius;
    int runtime_DerivativeExpansionRadius;

    int runtime_IsRamDisk;
    int runtime_DiskBufferSizeKB;
    int runtime_ConvolutionCacheSizeMB;
    int runtime_MaxConvolutionRAMMB;

    char runtime_DerivativesDirectory[1024];
    char runtime_MultipoleDirectory[1024];
    char runtime_TaylorDirectory[1024];

    char runtime_MultipolePrefix[1024];
    char runtime_TaylorPrefix[1024];
    
    int delete_multipoles_after_read;
    
    uint64_t blocksize, zwidth, rml, CompressedMultipoleLengthXY;
    int io_core;
} ConvolutionParameters;

#include "block_io_utils.cpp"

#ifdef CONVIOTHREADED
#include "ConvolutionIOThread.cpp"
#endif

class OutofCoreConvolution { 
public: 
    
    OutofCoreConvolution(void) { } 
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
    // Can use more of these for simultaneous IO on multiple devices
    ConvIOThread *iothread;
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
