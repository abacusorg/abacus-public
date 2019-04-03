/* This class will handle the Multipoles and Taylors in a parallel
in-memory setting.  The timestep code should compute MultipoleSlab
and then hand it to this code.  Here, we trigger MPI_Isends() and
MPI_Irecvs() to disperse the information out to the nodes that will handle
the individual z's.  This gets stored in a large buffer on a RAMdisk,
in [x][znode][m][y] order.  We provide a method to check that a 
given Multipole's transfer has been complete, so that MultipoleSlab
can be freed.

Then during the convolve stage, the information is transposed from
this buffer into a second buffer in [znode][m][y][x] order.  The
FFTs are called, then the InCoreConvolution, then the iFFTs.  The
information is then transposed back into the RAMdisk buffer.  This
routine also must load the appropriate derivatives from RAMdisk.

During timestep, we provide a method by which the TaylorSlab can
be requested, which allocates space and triggers the MPI call
to gather the information.  And a method to check that the MPI
work has been completed.
*/

#include "ParallelConvolution.h"

/// The zstart for the node of the given rank.
int ParallelConvolution::Zstart(int rank) {
    return (cpd2p1*rank/MPI_size); //TODO NAM sanity check: check that this returns integer division answers. 
	//TODO NAM does below code handle the case where some nodes have nothing to do or not a full range of z? 
}

// ParallelConvolution::ParallelConvolution(){}
/// The constructor.  This should open the MTfile as read/write or allocate
/// space if it doesn't exist.
ParallelConvolution::ParallelConvolution(int _cpd, int _order, char *MTfile, ConvolutionParameters &_ConvolveParams) : ConvolveParams(_ConvolveParams) { //TODO NAM does this always assume convolution overwrite mode? 
	cpd = _cpd;
	cpd2p1 = (cpd+1)/2; //NAM TODO why is range of z only half the cpd length? Definitely introduced some inconsistencies further down. Fix. 
	
	int pad = CPD2PAD/sizeof(Complex);
	cpd2pad = floor((cpd*cpd+pad)/pad)*pad;
	
	order = _order;
	rml = (order+1)*(order+1);
	zstart = Zstart(MPI_rank);
	znode = Zstart(MPI_rank+1)-zstart;
	this_node_size = znode*rml*cpd;
	
	assertf(sizeof(fftw_complex)==sizeof(Complex), "Mismatch of complex type sizing between FFTW and Complex\n");

	// Here we compute the start and number of z's in each node
	node_start = new int[MPI_size];
	node_size  = new int[MPI_size];
	for (int j=0; j<MPI_size; j++) {
	    node_start[j] = Zstart(j)*rml*cpd;
	    node_size[j]  = (Zstart(j+1)-Zstart(j))*rml*cpd;
	    assertf(node_size[j]*sizeof(MTCOMPLEX)<(((int64)1)<<31),
	    	"Amount of M/T data in a z range of one slab exceeds MPI 2 GB limit");
	}
	
	// create arrays to monitor statuses of Multipole MPI sends and receives. 
	Tsend_requests = new MPI_Request[cpd];
	for (int x = 0; x < cpd; x++) Tsend_requests[x] = MPI_REQUEST_NULL; 
	
	Msend_requests = new MPI_Request * [cpd];
	for (int x = 0; x < cpd; x++) Msend_requests[x] = NULL;  // This indicates no pending req
	Msend_active = 0;
	for (int x = 0; x < cpd; x++) Mrecv_requests[x] = MPI_REQUEST_NULL;  // This indicates no pending req

	Trecv_requests = new MPI_Request * [cpd];
	for (int x = 0; x < cpd; x++) Trecv_requests[x] = NULL;  // This indicates no pending req
	Trecv_active = 0;
	for (int x = 0; x < cpd; x++) Tsend_requests[x] = MPI_REQUEST_NULL;  // This indicates no pending req
	


	// However, a challenge is that on the receiving end, the x locations
	// are sufficiently discontinuous that we'll overflow the 32-bit indexing.
	// We might prefer to have the x's not monotonic, but rather node-cycling.
	// That could get fixed.
	
	//TODO NAM this assumes we're doing a convolution IO overwrite, since constructor only took one *MTfile! double check. 
	// mmap the MTfile to *MTdisk
    // For simplicity, multipoles and taylors must both be on ramdisk
	
	//Allocate MT and derivatives buffers and map to shared memory. 
	AllocMT(MTfile);
	AllocDerivs();	
}

/// The destructor.  This should close the MTfile. NAM: MTfile can be closed right after mapping, I think. 
/// Also can do any checks to make sure that all work was completed.
ParallelConvolution::~ParallelConvolution() {
	
	//NAM how to check work was completed, ideas:
	// assert all values of the four arrays below have been set to MPI_REQUEST_NULL? 
	// add class variable that keeps track of the number of FFTS, their success? check this? 
	
	assertf(Trecv_active == 0, "Error: active taylor receive requests ! = 0 \n");
	assertf(Msend_active == 0, "Error: active multipole send requests ! = 0 \n");
	
	delete[] Msend_requests;
	delete[] Trecv_requests; 
	delete[] Mrecv_requests; 
	delete[] Tsend_requests; 
	
	delete[] node_start; 
	delete[] node_size; 
	
	
    // TODO: do we need to munmap?  We do want the file to persist. NAM: i think so. munmap destroys the mapped virtual memory but not the file on disk. 
	size_t size = sizeof(MTCOMPLEX)*cpd*znode*cpd*rml; //size of the [x][znode][m][y] buffer on disk
    int res = munmap((void *) MTdisk, size);
    assertf(res == 0, "Failed to munmap MTdisk\n");
	
	size = sizeof(DFLOAT)*rml*ConvolveParams.CompressedMultipoleLengthXY; //size of the Derivatives buffer on disk. Double check this! Also, add CompressedMultipoleLengthXY to 
    res = munmap((void *) Ddisk, size);
    assertf(res == 0, "Failed to munmap Ddisk\n");

	
	for (int z = 0; z < znode; z++){
		delete[] Ddisk[z];
	}
	
	delete[] MTdisk; 
	delete[] Ddisk; 
			
	
   
}

/* =========================   ALLOC BUFFERS   ================== */ 


void ParallelConvolution::AllocMT(char * MTfile){
	//MTdisk has shape [x][znode][m][y]

    size_t size = sizeof(MTCOMPLEX)*znode*cpd*rml*cpd; //size of each x slice of MTdisk. 
    MTdisk = (MTCOMPLEX *) malloc(size);
	assert(MTdisk != NULL);
    
	//check that multipoles and taylors are on ramdisk. 
    assert(1 == is_path_on_ramdisk(ConvolveParams.runtime_MultipoleDirectory) == is_path_on_ramdisk(ConvolveParams.runtime_TaylorDirectory));
	
    size_t file_offset = cpd*zstart*cpd*rml*sizeof(MTCOMPLEX);
	
    int shm_fd_flags = O_RDWR;

    // map the shared memory fd to an address
    // If we're doing a ramdisk overwrite, this maps the multipoles directly into memory
	int fd = open(MTfile, shm_fd_flags, S_IRUSR | S_IWUSR);
    assertf(fd != -1, "Failed to open shared memory file at \"%s\"\n", MTfile);

	MTdisk = (MTCOMPLEX *) mmap(NULL, size, PROT_READ | PROT_WRITE, MAP_SHARED, fd, file_offset);

    int res = close(fd);
    assertf((void *) MTdisk != MAP_FAILED, "%d mmap shared memory from fd = %d of size = %d at offset = %d failed\n", MPI_rank, fd, size, file_offset);
    assertf(MTdisk != NULL, "mmap shared memory from fd = %d of size = %d at offset = %d failed\n", fd, size, file_offset);
    assertf(res == 0, "Failed to close fd %d\n", fd);
	
	return;
}

void ParallelConvolution::AllocDerivs(){
    Ddisk = new DFLOAT*[znode];
    for(int z = 0; z < znode; z++)
        Ddisk[z] = NULL; //NAM TODO  why not another malloc? taken from block_io utils. Double check. 
}

/* =========================   DERIVATIVES   ================== */ 

void ParallelConvolution::LoadDerivatives(int zstart , int zwidth) {
	
    //ReadDerivatives.Start(); //TODO NAM this currently exists nowhere else, fix. 
	
	size_t size = sizeof(DFLOAT)*rml*ConvolveParams.CompressedMultipoleLengthXY; //TODO NAM DFLOAT defined? What about ConvolveParams? 

    const char *fnfmt;
    if(sizeof(DFLOAT) == sizeof(float))
        fnfmt = "%s/fourierspace_float32_%d_%d_%d_%d_%d";
    else
        fnfmt = "%s/fourierspace_%d_%d_%d_%d_%d";
	
    // note the derivatives are stored in z-slabs, not x-slabs
    for(int z = zstart; z < zstart + zwidth; z++) {

        char fn[1024];
        sprintf(fn, fnfmt,
                ConvolveParams.runtime_DerivativesDirectory,
                (int) cpd, ConvolveParams.runtime_order, ConvolveParams.runtime_NearFieldRadius,
                ConvolveParams.runtime_DerivativeExpansionRadius, z);

        if(Ddisk[z-zstart] != NULL){ //NAM TODO why munmap here? confused... if it was previously mapped and we forgot to munmap? 
            int res = munmap(Ddisk[z-zstart], size);
            assertf(res == 0, "Failed to munmap derivs\n"); 
        }

        int fd = open(fn, O_RDWR, S_IRUSR | S_IWUSR);
        assertf(fd != -1, "Failed to open shared memory file at \"%s\"\n", fn);

        // map the shared memory fd to an address
        Ddisk[z-zstart] = (DFLOAT *) mmap(NULL, size, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
        int res = close(fd);
        assertf((void *) Ddisk[z-zstart] != MAP_FAILED, "mmap shared memory from fd = %d of size = %d failed\n", fd, size);
        assertf(Ddisk[z-zstart] != NULL, "mmap shared memory from fd = %d of size = %d failed\n", fd, size);
        assertf(res == 0, "Failed to close fd %d\n", fd);

    }

	//ReadDerivatives.Stop();
}
 

/* =========================  MPI Multipoles ================== */


/// Launch the CPD MPI_Irecv jobs.
/// Call this when the first slab is being finished.
/// We need to know the slabs that will be computed on each node.
void ParallelConvolution::RecvMultipoleSlab(int first_slab_finished) {
	int offset = first_slab_finished - first_slab_on_node;
	for (int r = 0; r < MPI_size; r++) {
	    for (int xraw = first_slabs_all[r] + offset; xraw < node_size[r]; xraw++) {
		int x = CP->WrapSlab(xraw);
		int tag = x + M_TAG;
		// We're receiving the full range of z's in each transfer.
		// Key thing is to route this to the correct x location
		MPI_Irecv(&MTdisk[x*this_node_size], this_node_size,
		    MPI_COMPLEX, r, tag, MPI_COMM_WORLD, &Mrecv_requests[x]);
	    }
	}
	STDLOG(1,"MPI_Irecv set for incoming Multipoles\n");
	return;
}


/// Take 1 MultipoleSlab and trigger the MPI calls to disperse
/// it to the other nodes.
void ParallelConvolution::SendMultipoleSlab(int slab) {
	slab = CP->WrapSlab(slab);
    MTCOMPLEX * mt = (MTCOMPLEX *) SB->GetSlabPtr(MultipoleSlab, slab);
	Msend_requests[slab] = new MPI_Request[MPI_size];
	Msend_active++;
	for (int r = 0; r < MPI_size; r++) {
    	int tag = slab + M_TAG;
    	MPI_Isend(mt + node_start[r], node_size[r], MPI_COMPLEX, r, tag, MPI_COMM_WORLD, &Msend_requests[slab][r]);
	}
	STDLOG(1,"Multipole slab %d has been queued for MPI_Isend\n", slab);
}


/// Return only once all of the multipole slabs have been received.
/// We'll only check this at the end of the timestep loop, so 
/// this will just spin until done.
void ParallelConvolution::WaitRecvMultipoleComplete() {
	while (1) //spin until done. 
	{
		int err = MPI_Waitall(cpd, Mrecv_requests, MPI_STATUSES_IGNORE);  // TODO: NAM okay to have MPI_STATUSES_IGNORE, or need array_of_statuses? 
		if (err == 0) return; 
	}
}

/// Call periodically to check on the status of the MPI calls
/// and delete the MultipoleSlab when a send is complete.
/// Probably best if this is not expecting the completion to be 
/// in order, but we could discuss this.
int ParallelConvolution::CheckSendMultipoleComplete() {
    for (int x = 0; x < cpd; x++) {
		
	    if (Msend_requests[x]==NULL) continue;  // Nothing to do for this slab
		
	    // This slab has pending MPI calls
	    int err, sent = 0, done = 1;
		
	    for (int r = 0; r < MPI_size; r++) {
			if (Msend_requests[x][r] == MPI_REQUEST_NULL) continue;  // Already done
			
			err = MPI_Test(Msend_requests[x] + r, &sent, MPI_STATUS_IGNORE);
			if (not sent) {
			    // Found one that is not done
			    done=0; break;
			}
	    }
	    
		if (done) {
			STDLOG(1,"MPI_Isend complete for MultipoleSlab %d\n", x);
	        delete[] Msend_requests[x];
			Msend_requests[x] = NULL;   // This will allow us to check faster
			Msend_active--;
			SB->DeAllocate(MultipoleSlab, x);   // Can delete it ---- NAM TODO if this gets deleted, never written, then where did it come from in the first place? 
			// TODO: May need to force this to skip writing out?
	    }
	}
}

/// Call after the timestep loop to spin until all of the
/// MPI calls have been finished.
void ParallelConvolution::WaitForMultipoleTransferComplete() {
	STDLOG(1,"Begin waiting to determine that all incoming Multipoles received\n");
	WaitRecvMultipoleComplete(); //this spins until recvs are done. 
	STDLOG(1,"End waiting to determine that all incoming Multipoles received\n");
	STDLOG(1,"Begin waiting to determine that all outgoing Multipoles sent\n");
	while (Msend_active>0) 
		CheckSendMultipoleComplete(); //this spins until sends are done. 
	STDLOG(1,"End waiting to determine that all outgoing Multipoles sent\n");
}

/* ======================= MPI Taylors ========================= */


// TODO: Need to figure out how to pull the Taylors over.
// Need to be careful about the ordering of the Isend's.
// If we queue up all of the Isend's, will we flood the buffers and
// block the Manifest?  But we have no guarantee that the nodes will
// stay in timing agreement.  Do we need to provision to hold all of the 
// TaylorSlabs at once, as opposed to providing on request?
// Would having a separate thread to field requests be safer?  but multi-threading MPI
// Would using a different COMM for the Manifest help?
// If we put tags on every message, might it help to avoid flooding a buffer?

///for a given x slab = slab, figure out which node rank is responsible for it and should receive its Taylors. 
int ParallelConvolution::GetTaylorRecipient(int slab){
	for (int r = 0; r < MPI_size; r++){
		if (ConvolveParams.first_slabs_all[r] <= slab < ConvolveParams.first_slabs_all[r] + ConvolveParams.total_slabs_all[r])
			return r; 
	}
}

void ParallelConvolution::SendTaylorSlab(int slab) {
	//Each node has Taylors for a limited range of z but for every x slab. 
	//figure out who the receipient should be based on x slab and send to them. 
	// Take from MTdisk. Set Tsend_requests[x] as request. 
	
	slab = CP->WrapSlab(slab);
	int r = GetTaylorRecipient(slab); //x-slab slab is in node r's domain. Send to node r. 
	
	int tag = slab + T_TAG; 
	MPI_Isend(MTdisk + slab*this_node_size, this_node_size, MPI_COMPLEX, r, tag, MPI_COMM_WORLD, &Tsend_requests[slab]);
	
	STDLOG(1,"Taylor slab %d has been queued for MPI_Isend to rank %d\n", slab, r);
}
	
/// Allocate space for this TaylorSlab, and trigger the MPI calls
/// to gather the information from other nodes.  To be used in FetchSlab.
void ParallelConvolution::RecvTaylorSlab(int slab) {	
	//We are receiving a specified x slab. For each node, receive its z packet for this x. For each node, use Trecv_requests[x][rank] as request. 
	slab = CP->WrapSlab(slab);
    MTCOMPLEX *mt = (MTCOMPLEX *) SB->GetSlabPtr(TaylorSlab, slab);
	
	Trecv_requests[slab] = new MPI_Request[MPI_size];
	Trecv_active++; 
	
	for (int r = 0; r < MPI_size; r++) {
		int tag = slab + T_TAG; 
		MPI_Irecv(mt + node_start[r], node_size[r], MPI_COMPLEX, r, tag, MPI_COMM_WORLD, &Trecv_requests[slab][r]);
	}
	STDLOG(1,"MPI_Irecv set for incoming Taylors\n");
	return;
}

/// Check whether the MPI calls are complete for this TaylorSlab,
/// returning 1 if yes.  To be used in the TaylorForcePrecondition.
/// Also, have to be careful about the SendManifest.
int ParallelConvolution::CheckRecvTaylorComplete(int slab) {
	slab = CP->WrapSlab(slab);
	
    if (Trecv_requests[slab]==NULL) return 1;  // Nothing to do for this slab
	
    // This slab has pending MPI calls
    int err, sent=0, done=1;
	
    for (int r=0; r<MPI_size; r++) {
		if (Trecv_requests[slab][r]==MPI_REQUEST_NULL) continue;  // Already done
		
		err = MPI_Test(Trecv_requests[slab]+r,&sent,MPI_STATUS_IGNORE);
		if (not sent) {
		    // Found one that is not done
		    done=0; break;
		}
    }
    
	if (done) {
		STDLOG(1,"MPI_Irecv complete for TaylorSlab %d\n", slab);
        delete[] Trecv_requests[slab];
		Trecv_requests[slab] = NULL;   // This will allow us to check faster
		Trecv_active--;
    }
	
	return done; 
}



/* =========================  Convolve functions ================== */

/// Create the MTzmxy buffer and copy into it from MTdisk.
// Note that this is also a type-casting from float to double.
// I wonder if this is why this step has been slow?
void ParallelConvolution::Swizzle_to_zmxy() {
    size_t bufsize = sizeof(Complex)*znode*rml*cpd2pad; //this implies that cpd2pad is the size of x * y for a given z. 
	STDLOG(1,"Allocating %d bytes for MTzmxy buffer\n", bufsize);
	int ret = posix_memalign((void **)MTzmxy, 4096, bufsize);
	assertf(ret==0, "Allocation of MTzmxy failed\n");
	
	// TODO: Do the copies, including the -1 x offset
	// Converting from [x][znode][m][y] to [znode][m][x][y]
	// Loop over the combination of [znode][m]
	#pragma omp parallel for schedule(static)
	for (int zm = 0; zm < cpd2p1*rml; zm++) {
	    int zoffset = zm/rml;   // This is not actually z, but z-zstart
	    int m = zm-zoffset*rml; 
		
	    for (int x=0; x<cpd; x++) {
			int xuse = x+1; if (xuse>=cpd) xuse=0;   
			
			for (int y=0; y<cpd; y++) {
			    MTzmxy[(int64)zm*cpd2pad + xuse*cpd + y] =
			        MTdisk[(((int64)x*znode + zoffset)*rml + m)*cpd + y];
			}
	    }
	}
	return;
}

/// Copy from MTzmxy buffer back to MTdisk, then free buffer
void ParallelConvolution::Swizzle_to_xzmy() {
	// TODO: Do the copies, including the -1 x offset
	// TODO: Might prefer a different pragma structure than the other Swizzle
	// e.g., to combine the [x][znode] loop, so that each thread writes contiguous [m][y].

	#pragma omp parallel for schedule(static)
	for (int xz = 0; xz < cpd * cpd2p1; xz ++) {
		int xoffset = xz / cpd2p1; 
		int z = xz - xoffset * cpd2p1; 
		
		int x = xoffset + 1; //NAM. i think? 
		
		for (int m = 0; m < rml; m ++){
			for (int y = 0; y < cpd; y++){
				MTdisk[ (int64)((xz  * rml + m)) * cpd + y ] = 
					MTzmxy[ z * cpd2pad * rml + m * cpd2pad + x * cpd + y ]; 
			}
		}
	}
	
	
    free(MTzmxy);
	return;
}

/// Create the plan and do the FFT
void ParallelConvolution::FFT() {
	// TODO: Import wisdom
	// We're going to plan to process each [x][y] slab separately,
	// doing the cpd FFTs of cpd length.
	// stride = cpd, dist = 1, nembed = NULL
	
	int err = fftw_init_threads(); 
	assertf(err == 0, "Failed to fftw_init_threads\n");
	
	fftw_plan_with_nthreads(omp_get_max_threads()); //NAM TODO how many threads in here? 
	
	fftw_plan plan;
	int * n = new int[cpd]; 
	for(int i = 0; i < cpd; i ++) n[i] = cpd; 
	
	plan = fftw_plan_many_dft(1, n, cpd, 
		(fftw_complex *) MTzmxy, NULL, cpd, 1,
		(fftw_complex *) MTzmxy, NULL, cpd, 1,
		FFTW_FORWARD, FFTW_PATIENT);
	// TODO: Is it ok to have only one plan but use many threads?

	// Now we loop over the combination of [znode][m]
	#pragma omp parallel for schedule(static)
	for (int zm = 0; zm < cpd2p1*rml; zm++) {
	    fftw_execute_dft(plan, (fftw_complex *)MTzmxy+zm*cpd2pad,
			       (fftw_complex *)MTzmxy+zm*cpd2pad);
	}
	
	fftw_destroy_plan(plan);
	delete[] n; 
	
	//fftw_cleanup_threads() NAM TODO do we need this?  or maybe in destructor? 
	return;
}

void ParallelConvolution::iFFT() {
	// Do we need to re-import WISDOM?
	// We're going to plan to process each [x][y] slab separately,
	// doing the cpd FFTs of cpd length.
	// stride = cpd, dist = 1, nembed = NULL
	
	int err = fftw_init_threads(); 
	assertf(err == 0, "Failed to fftw_init_threads\n");
	
	fftw_plan_with_nthreads(omp_get_max_threads()); //NAM TODO how many threads in here? 
	
	fftw_plan plan;
	int * n = new int[cpd]; 
	for(int i = 0; i < cpd; i ++) n[i] = cpd; 
	
	plan = fftw_plan_many_dft(1, n, cpd, 
		(fftw_complex *) MTzmxy, NULL, cpd, 1,
		(fftw_complex *) MTzmxy, NULL, cpd, 1,
		FFTW_BACKWARD, FFTW_PATIENT);

	// Now we loop over the combination of [znode][m]
	#pragma omp parallel for schedule(static)
	for (int zm = 0; zm < cpd2p1*rml; zm++) {
	    fftw_execute_dft(plan, (fftw_complex *)MTzmxy+zm*cpd2pad,
			       (fftw_complex *)MTzmxy+zm*cpd2pad);
	}
	fftw_destroy_plan(plan);
	delete[] n; 
	
	//fftw_cleanup_threads() NAM TODO do we need this?  or maybe in destructor? 
	return;
}

/// This executes the convolve on the MTfile info.
void ParallelConvolution::Convolve() {
	
    InCoreConvolution *ICC = new InCoreConvolution(order,cpd,blocksize);
	
	
	// We're beginning in [x][znode][m][y] order on the RAMdisk
	// Copy to the desired order in a new buffer
	Swizzle_to_zmxy();   // Probably apply the -1 x offset here
	FFT();

	// Two options: 
	// a) Swizzle to zmyx, then do a bunch of already-contiguous FFTs, and alter ICC
	// b) Swizzle to zmxy, then do a bunch of FFTs with a stride; use ICC as is
	// If option a, then the swizzle should be a smarter transpose.
	// Option b can be boring code.

	for (int z = zstart; z < zstart + znode; z++) {
	    LoadDerivatives(z);
		
		//void InCoreConvolution::InCoreConvolve(Complex *FFTM, DFLOAT *CompressedD) {
		// InCoreConvolution(int order, int cpd, int blocksize) 
		//NAM TODO careful! does ICC know about doing one z at a time? 

	    ICC->InCoreConvolve(MTzmxy, Ddisk[z]); // TODO NAM add z argument. , z);
	}
	iFFT();
	Swizzle_to_xzmy();   // Restore to the RAMdisk
	
	delete ICC;
}


