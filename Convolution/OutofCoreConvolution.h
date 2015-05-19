//#include "read_dio.h"
//#include "write_dio.h"

typedef struct { 
    double ReadDerivatives;
    double ReadMultipoles;
    double WriteTaylor;
    double ForwardZFFTMultipoles;
    double InverseZFFTTaylor;
    double ConvolutionArithmetic;
    double ArraySwizzle;
    
    double ConvolveWallClock;
    double Discrepency;

    unsigned long long int ReadDerivativesBytes;
    unsigned long long int ReadMultipolesBytes;
    unsigned long long int WriteTaylorBytes;
    unsigned long long int ops;
    unsigned long long int totalMemoryAllocated;

    int runtime_ConvolutionCacheSizeMB;
    int blocksize, zwidth;
    int runtime_order;
    int runtime_cpd;

} ConvolutionStatistics;

typedef struct {
    int runtime_cpd;
    int runtime_order;
    int runtime_NearFieldRadius;
    int runtime_DerivativeExpansionRadius;

    int runtime_IsRamDisk;
    int runtime_DiskBufferSizeKB;
    int runtime_ConvolutionCacheSizeMB;
    long long unsigned int runtime_MaxConvolutionRAMMB;

    char runtime_DerivativesDirectory[1024];
    char runtime_MultipoleDirectory[1024];
    char runtime_TaylorDirectory[1024];

    char runtime_MultipolePrefix[1024];
    char runtime_TaylorPrefix[1024];
} ConvolutionParameters;

class OutofCoreConvolution { 
public: 
    
    OutofCoreConvolution(void) { } 
    ~OutofCoreConvolution(void) { } 

    ConvolutionParameters CP;
    void Convolve( ConvolutionParameters CP );

    long long unsigned int blocksize, zwidth;

    ConvolutionStatistics CS; 

private:

    STimer ReadDerivatives;
    STimer ReadMultipoles;
    STimer WriteTaylor;

    STimer ForwardZFFTMultipoles;
    STimer InverseZFFTTaylor;

    STimer ConvolutionArithmetic;
    STimer ArraySwizzle;

    STimer ConvolveWallClock;


    long long unsigned int cpd,order,rml,CompressedMultipoleLengthXY;

    void BlockConvolve(void);
    void WriteDiskTaylor(int z);
    void ReadDiskDerivatives(int z);
    void ReadDiskMultipoles(int z);

    int mapM[8192];
    int remap[8192];

    ReadDirect *RD_RDD;
    ReadDirect *RD_RDM;
    WriteDirect *WD_WDT;

    Complex *DiskBuffer;
    double *CompressedDerivatives;
    Complex *TemporarySpace;
};
