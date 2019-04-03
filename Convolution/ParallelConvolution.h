#define CPD2PAD 1024    // Byte padding



/// The ParallelConvolution class is responsible for performing the Fourier space convolution of the multipole moments with the derivatives tensor to obtain the Taylor coefficients required to compute the far field, in the use case where Abacus is running on multipole nodes. 
/// ParallelConvolution:
/// 1) handles the MPI work that re-distributes existing multipole moment data to all nodes from each node (which has multipoles generated during singlestep for a limited set of x and for all z, but requires multipoles for all x and for a limitet set of z to perform the convolution). 
/// 2) performs array manipulation, the FFT, the convolution (via InCoreConvolution.cpp),  the inverse FFT, and more array manipulation. 
/// 3) does the MPI work to send the Taylors back to the nodes they live on (each node needs, to run singlestep, Taylors for a set of x and for all z). 

class ParallelConvolution { 
public: 

    //ConvolutionStatistics CS; 
    
	// ParallelConvolution();
    ParallelConvolution(int _cpd, int _order, char *MTfile, ConvolutionParameters &_CP);
    ~ParallelConvolution(void);
	
    /// The zstart for the node of the given rank.
    int Zstart(int rank);  
	
    void Convolve();
	
	void RecvMultipoleSlab(int first_slab_finished);
	void SendMultipoleSlab(int slab);
	void WaitForMultipoleTransferComplete();
	
	void RecvTaylorSlab(int slab);
	void SendTaylorSlab(int slab); //TODO do we need this?
	int  CheckRecvTaylorComplete(int slab);
	
    uint64_t blocksize, zwidth, z_slabs_per_node; 

private:
	ConvolutionParameters ConvolveParams; //NAM TODO do we really need this? 
	
    uint64_t cpd;
    uint64_t cpd2p1;   // (CPD+1)/2 is the length of the z array
    uint64_t cpd2pad;  // We might want to pad CPD**2 to obtain better FFTW behavior
    int order;    // multipole order
    uint64_t rml;      // (order+1)**2
	uint64_t this_node_size; 

    int zstart;   // The first z for this node
    int znode;    // The number of z's on this node
	
	int * node_start;
	int * node_size; 

    MTCOMPLEX * MTdisk;   // This is the pointer to the [x][znode][m][y] buffer on disk //NAM TODO * or ** ?
	DFLOAT ** Ddisk; //This is the pointer to the derivatives buffer. 
    Complex * MTzmxy;   // This is the pointer to the [znode][m][x][y] work buffer 

    #define M_TAG (10000)     // An offset for the multipole tags
    MPI_Request *Mrecv_requests;    // We'll set up a [CPD] array to receive one packet per x
    MPI_Request **Msend_requests;    // We'll set up a [CPD][MPI_size] array even though each node will only use the x's in its NodeSlabs range
    int Msend_active;

    #define T_TAG (20000)     // An offset for the taylor tags
    MPI_Request *Tsend_requests;    // We'll set up a [CPD] array to send one packet per x
    MPI_Request **Trecv_requests;    // We'll set up a [CPD][MPI_size] array even though each node will only use the x's in its NodeSlabs range
	int Trecv_active; 
	
	void AllocMT(char * MTfile);
	void AllocDerivs(); 
	
	void LoadDerivatives(int zstart, int zwidth = 1);
	void Swizzle_to_zmxy();
	void Swizzle_to_xzmy();
	void FFT();
	void iFFT();
	
	int  CheckSendMultipoleComplete();
	void WaitRecvMultipoleComplete();
	
	int GetTaylorRecipient(int slab);
	
    // STimer ForwardZFFTMultipoles;
//     STimer InverseZFFTTaylor;
//
//     STimer ConvolutionArithmetic;
//     STimer ArraySwizzle;
//
//     STimer ConvolveWallClock;
//
//     uint64_t cpd,order,rml,CompressedMultipoleLengthXY;
//
//     void BlockConvolve(void);
//     void WriteDiskTaylor(int z);
//     void ReadDiskMultipolesAndDerivs(int z);
//     void SwizzleMultipoles(int z);
//     void SwizzleTaylors(int z);
//
//     void RenameMultipolesToTaylors();

    // Complex *PlaneBuffer;
 //    DFLOAT **CompressedDerivatives;
 //    MTCOMPLEX **DiskBuffer;
 //    Block *CurrentBlock;
    
    // double invcpd3;
};

ParallelConvolution *ParallelConvolveDriver;

