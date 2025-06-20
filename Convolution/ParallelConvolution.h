/*
 * Copyright 2012-2025 The Abacus Developers
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#define CPD2PADDING 1024    // Byte padding

/// The ParallelConvolution class is responsible for performing the Fourier space convolution of the multipole moments with the derivatives tensor to obtain the Taylor coefficients required to compute the far field, in the use case where Abacus is running on multipole nodes. 
/// ParallelConvolution:
/// 1) handles the MPI work that re-distributes existing multipole moment data to all nodes from each node (which has multipoles generated during singlestep for a limited set of x and for all z, but requires multipoles for all x and for a limitet set of z to perform the convolution). 
/// 2) performs array manipulation, the FFT, the convolution (via InCoreConvolution.cpp),  the inverse FFT, and more array manipulation. 
/// 3) does the MPI work to send the Taylors back to the nodes they live on (each node needs, to run singlestep, Taylors for a set of x and for all z). 

#include "mpi_limiter.h"

#ifdef DOUBLEPRECISION
#define MPI_MTCOMPLEX MPI_C_DOUBLE_COMPLEX
#else
#define MPI_MTCOMPLEX MPI_C_COMPLEX
#endif

class ParallelConvolution { 
public:     
	// ParallelConvolution();
    ParallelConvolution(int _cpd, int _order, const fs::path &MultipoleDirectory, int _create_MT_file = 0);
    ~ParallelConvolution(void);
		
	void AllocMT_two();	
	void AllocMT();
	void AllocDerivs(); 
	
    /// The zstart for the node of the given rank.
    int Zstart(int rank);  
	
    void Convolve();
	
	void RecvMultipoleSlab(int first_slab_finished);
	void SendMultipoleSlab(int slab);
	void WaitForMultipoleTransferComplete(int offset);
	
	int  CheckForMultipoleTransferComplete(int slab);
	
	
	void RecvTaylorSlab(int slab);
	void SendTaylors(int offset); //TODO do we need this?
	int  CheckTaylorSlabReady(int slab);
	
	//timers and stats to log:
    ConvolutionStatistics CS; 
    int64_t blocksize, znode; // The number of z's on this node
	void dumpstats();

    STimer *MsendTimer;   // A timer for the sending of each Multipole Slab

private:
	
	fs::path mt_file; 
	int create_MT_file = 0; 
	
	//int OverwriteConvState;
	//int StripeConvState;
	
    int64_t cpd;
	int64_t node_ky_size;  // length of ky, which is split over nodes in 2D
	int64_t kz_size;  // length of kz, different in 1D and 2D
    int64_t cpdky_pad;  // We might want to pad CPD*node_ky_size to obtain better FFTW behavior
    int order;    // multipole order
    int64_t rml;      // (order+1)**2
	int64_t node_slab_elem;  // number of elements in each slab this node handles
	int64_t CompressedMultipoleLengthXY; 
	
	Complex invcpd3; 

	size_t mt_offset;
	int ramdisk_derivs = -1; 
	
    int zstart;   // The first z for this node
	
	int * node_start;
	int * node_size; 
	
	int * first_slabs_all;
	int * total_slabs_all; 
	
	int first_slab_on_node;
	int total_slabs_on_node;

    MTCOMPLEX * MTdisk;   // This is the pointer to the [x][znode][m][y] buffer on disk //NAM TODO * or ** ?
	DFLOAT ** Ddisk; //This is the pointer to the derivatives buffer. 
    Complex * MTzmxy;   // This is the pointer to the [znode][m][x][y] work buffer 

	size_t Ddisk_bytes;  // bytes per derivatives slab
	size_t MTdisk_bytes;  // MTdisk bytes
	off_t dfile_offset;  // where to start reading the derivatives file (2D)

    #define M_TAG (100000)     // An offset for the multipole tags
    MPI_Request *Mrecv_requests;    // We'll set up a [CPD] array to receive one packet per x
    MPI_Request **Msend_requests;    // We'll set up a [CPD][MPI_size] array even though each node will only use the x's in its NodeSlabs range
	int *Mrecv_flags;
	int **Msend_flags;
    int Msend_active;

    #define T_TAG (105000)     // An offset for the taylor tags; we step by 10000 per rank, so this is good to CPD=5000
    MPI_Request *Tsend_requests;    // We'll set up a [CPD] array to send one packet per x
    MPI_Request **Trecv_requests;    // We'll set up a [CPD][MPI_size] array even though each node will only use the x's in its NodeSlabs range
	int *Tsend_flags;
	int **Trecv_flags;
	int Trecv_active; 

    AbacusMPILimiter mpi_limiter;
    int *TaylorSlabAllMPIDone;
    int *MultipoleSlabAllMPIDone;

	ReadDirect *RD;
	
	
	fs::path MultipoleFN(int slab);
	fs::path TaylorFN(int slab);
	
	void LoadDerivatives(int z);
	void Swizzle_to_zmxy();
	void Swizzle_to_xzmy();
	
	fftw_plan PlanFFT(int sign);
	void FFT(fftw_plan plan);
		
	int  CheckSendMultipoleComplete(int slab);
	void WaitRecvMultipoleComplete(int slab);
	
	
	int CheckRecvMultipoleComplete(int slab);
	
	
	int GetTaylorRecipient(int slab, int offset);
	int CheckTaylorRecvReady(int slab);
	int CheckTaylorSendComplete(int slab);
		
};

ParallelConvolution *ParallelConvolveDriver;

