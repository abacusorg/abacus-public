#include "ConvolutionParametersStatistics.cpp" //NAM TODO haven't tested if this compiles. 
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

    int64_t blocksize, zwidth;
	
#ifdef PARALLEL
	int64_t z_slabs_per_node; 
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

    int64_t cpd,order,rml,CompressedMultipoleLengthXY;

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
