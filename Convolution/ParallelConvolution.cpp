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
//#define DO_NOTHING

#include "ParallelConvolution.h"

/// The zstart for the node of the given rank.
int ParallelConvolution::Zstart(int rank) {
    return (cpd2p1*rank/MPI_size); //TODO NAM sanity check: check that this returns integer division answers. 
	//TODO NAM does below code handle the case where some nodes have nothing to do or not a full range of z? 
}

// ParallelConvolution::ParallelConvolution(){}
/// The constructor.  This should open the MTfile as read/write or allocate
/// space if it doesn't exist.
ParallelConvolution::ParallelConvolution(int _cpd, int _order, const char *MTfile, int _create_MT_file) { //TODO NAM does this always assume convolution overwrite mode? 
	
	STDLOG(1, "Calling constructor for ParallelConvolution\n");
	
    StripeConvState = strcmp(P.Conv_IOMode, "stripe") == 0;
    OverwriteConvState = strcmp(P.Conv_IOMode, "overwrite") == 0;	
	
	mt_file = MTfile;
	create_MT_file = _create_MT_file; 
	cpd = _cpd;
	cpd2p1 = (cpd+1)/2; //NAM TODO why is range of z only half the cpd length? Definitely introduced some inconsistencies further down. Fix. 
	
	int pad = CPD2PAD/sizeof(Complex);
	cpd2pad = floor((cpd*cpd+pad)/pad)*pad;
	
	order = _order;
	rml = (order+1)*(order+1);
	zstart = Zstart(MPI_rank);
	znode = Zstart(MPI_rank+1)-zstart;
	this_node_size = znode*rml*cpd;
	
	STDLOG(1, "Doing zstart = %d and znode = %d\n", zstart, znode);
	
	
	CompressedMultipoleLengthXY = ((1+P.cpd)*(3+P.cpd))/8;
		 
    int cml = ((order+1)*(order+2)*(order+3))/6;
    int nprocs = omp_get_max_threads();
	size_t cacherambytes = P.ConvolutionCacheSizeMB*(1024LL*1024LL);

    blocksize = 0;
    for(blocksize=cpd*cpd;blocksize>=2;blocksize--) 
        if((cpd*cpd)%blocksize==0)
            if(nprocs*2.5*cml*blocksize*sizeof(Complex) < cacherambytes) break;
                // 2.5 = 2 Complex (mcache,tcache) 1 double dcache

    
	STDLOG(1, "Conv params: %d %d %d %d %d\n", cml, nprocs, cacherambytes, blocksize, CompressedMultipoleLengthXY);
	
	first_slabs_all = new int[MPI_size];
    total_slabs_all = new int[MPI_size]; 
	
	ReadNodeSlabs(1, first_slabs_all, total_slabs_all);
	
	first_slab_on_node  = first_slabs_all[MPI_rank];
	total_slabs_on_node = total_slabs_all[MPI_rank];
	
	
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
	
	Msend_requests = new MPI_Request * [cpd];	
	Mrecv_requests = new MPI_Request[cpd];
	Trecv_requests = new MPI_Request * [cpd];
	Tsend_requests = new MPI_Request[cpd];
	
	Msend_active = 0;
	Trecv_active = 0;


	for (int x = 0; x < cpd; x++) {
		Msend_requests[x] = NULL;  // This indicates no pending req
		Mrecv_requests[x] = MPI_REQUEST_NULL;  // This indicates no pending req
		Trecv_requests[x] = NULL;  // This indicates no pending req
		Tsend_requests[x] = MPI_REQUEST_NULL;  // This indicates no pending req
	} 
	
	AllocMT();
	AllocDerivs();
}

/// The destructor.  This should close the MTfile. NAM: MTfile can be closed right after mapping, I think. 
/// Also can do any checks to make sure that all work was completed.
ParallelConvolution::~ParallelConvolution() {
	
	STDLOG(1, "Entering desctructor\n")
	//NAM how to check work was completed, ideas:
	// assert all values of the four arrays below have been set to MPI_REQUEST_NULL? 
	// add class variable that keeps track of the number of FFTS, their success? check this? 
	
	assertf(Trecv_active == 0, "Error: active taylor receive requests ! = 0, = %d\n", Trecv_active);
	assertf(Msend_active == 0, "Error: active multipole send requests ! = 0, = %d\n", Msend_active);
	
	delete[] Msend_requests;
	delete[] Trecv_requests; 
	delete[] Mrecv_requests; 
	delete[] Tsend_requests; 
	
	delete[] node_start; 
	delete[] node_size; 
	
	
	delete[] first_slabs_all;
	delete[] total_slabs_all; 
	
    size_t size = sizeof(MTCOMPLEX)*this_node_size;
    int res = munmap(MTdisk, size);
    assertf(res == 0, "munmap failed\n");
	
    for (int z = 0; z < znode; z++){
		delete[] Ddisk[z];
	}
	
	delete[] Ddisk; 
	
	fftw_cleanup_threads();
	
	STDLOG(1, "Exiting destructor\n")
			
	
   
}

/* TODO move these back to CP.cpp */

void ParallelConvolution::MultipoleFN(int slab, char * const fn){
    // We elsewhere generically support N threads, but here is where we assume 2
    if(StripeConvState && slab % 2 == 1)
        sprintf(fn, "%s/%s_%04d", P.MultipoleDirectory2, "Multipoles", slab);
    else
        sprintf(fn, "%s/%s_%04d", P.MultipoleDirectory, "Multipoles", slab);
}

void ParallelConvolution::TaylorFN(int slab, char * const fn){
    // We elsewhere generically support N threads, but here is where we assume 2
    if(StripeConvState && slab % 2 == 1)
        sprintf(fn, "%s/%s_%04d", P.TaylorDirectory2, "Taylor", slab);
    else
        sprintf(fn, "%s/%s_%04d", P.TaylorDirectory, "Taylor", slab);
}

/* =========================   ALLOC BUFFERS   ================== */ 


void ParallelConvolution::AllocMT(){ //could add create flag to this to handle file creation and pass in 1 from timestep_IC and default to 0 otherwise. 
	//MTdisk has shape [x][znode][m][y]
	
	//TODO this appears almost exactly as in arena allocator. Consider factoring out to ramdisk util. 
	
	assert(mt_file != NULL);
    assert(strnlen(mt_file,1) > 0);
    STDLOG(1,"Mapping MTdisk from shared memory\n");

	if(create_MT_file)
		unlink(mt_file);
	
    int shm_fd_flags = create_MT_file ? (O_RDWR | O_CREAT | O_EXCL) : O_RDWR;
    int fd = open(mt_file, shm_fd_flags, S_IRUSR | S_IWUSR);
    assertf(fd != -1, "Failed to open shared memory file at \"%s\"\n", mt_file);
    
	size_t size = sizeof(MTCOMPLEX)*this_node_size*cpd; //size of MTdisk.
	
    if(create_MT_file){
	    // set the size
	    int res = ftruncate(fd, size);
	    assertf(res == 0, "ftruncate on shared memory mt_file = %s to size = %d failed\n", mt_file, size);
    } else {
        // check that the size is as expected
        struct stat shmstat;
        int res = fstat(fd, &shmstat);
        assertf(res == 0, "fstat on shared memory ramdisk_fn = %s\n", mt_file);
        assertf(shmstat.st_size == size, "Found shared memory size %d; expected %d (ramdisk_fn = %s)\n", shmstat.st_size, size, mt_file);
    }
    // map the shared memory fd to an address
    int mmap_flags = create_MT_file ? (PROT_READ | PROT_WRITE) : (PROT_READ | PROT_WRITE);  // the same
	MTdisk = (MTCOMPLEX *) mmap(NULL, size, mmap_flags, MAP_SHARED, fd, 0);
    int res = close(fd);
    assertf((void *) MTdisk != MAP_FAILED, "mmap shared memory from fd = %d of size = %d failed\n", fd, size);
    assertf(MTdisk != NULL, "mmap shared memory from fd = %d of size = %d failed\n", fd, size);
    assertf(res == 0, "Failed to close fd %d\n", fd);
	
	STDLOG(1, "Successfully mapped MTdisk of size %d bytes.\n", size)
		
    size_t bufsize = sizeof(Complex)*znode*rml*cpd2pad; 	

	// int ret = posix_memalign((void **)&MTzmxy, 4096, bufsize);
// 	assertf(ret==0, "Failed to memalign MTzmxy\n");
// 	STDLOG(1,"Allocated %d bytes for MTzmxy buffer\n", bufsize);


	MTzmxy = (Complex *) fftw_malloc(bufsize);
	assertf(MTzmxy!=NULL, "Failed fftw_malloc for MTzmxy\n");
	STDLOG(1, "Successfully malloced MTzmxy\n");
	return;
}

void ParallelConvolution::AllocDerivs(){
	
	STDLOG(1,"AllocDerivs a, %d\n", MPI_rank);
	
    Ddisk = new DFLOAT*[znode];
	
	STDLOG(1,"AllocDerivs b, %d\n", MPI_rank);
	
	
    for(int z = 0; z < znode; z++)
        Ddisk[z] = NULL; //NAM TODO  why not another malloc? taken from block_io utils. Double check. 
	
	STDLOG(1,"AllocDerivs c, %d\n", MPI_rank);
	
}

/* =========================   DERIVATIVES   ================== */ 

void ParallelConvolution::LoadDerivatives(int z) {
	int z_file = z + zstart;
	STDLOG(1,"Entering LoadDerivatives z %d from file z_file = %d, zstart = %d\n", z, z_file, zstart);
	
    //ReadDerivatives.Start(); 
	
	size_t size = sizeof(DFLOAT)*rml*CompressedMultipoleLengthXY;

    const char *fnfmt;
    if(sizeof(DFLOAT) == sizeof(float))
        fnfmt = "%s/fourierspace_float32_%d_%d_%d_%d_%d";
    else
        fnfmt = "%s/fourierspace_%d_%d_%d_%d_%d";
	
    // note the derivatives are stored in z-slabs, not x-slabs
    char fn[1024];
    sprintf(fn, fnfmt,
            P.DerivativesDirectory,
            (int) cpd, order, P.NearFieldRadius,
            P.DerivativeExpansionRadius, z_file);

    if(Ddisk[z] != NULL){ //NAM TODO why munmap here? 
        int res = munmap(Ddisk[z], size);
        assertf(res == 0, "Failed to munmap derivs\n"); 
    }

    int fd = open(fn, O_RDWR, S_IRUSR | S_IWUSR);
    assertf(fd != -1, "Failed to open shared memory file at \"%s\"\n", fn);

    // map the shared memory fd to an address
    Ddisk[z] = (DFLOAT *) mmap(NULL, size, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
    int res = close(fd);
    assertf((void *) Ddisk[z] != MAP_FAILED, "mmap shared memory from fd = %d of size = %d failed\n", fd, size);
    assertf(Ddisk[z] != NULL, "mmap shared memory from fd = %d of size = %d failed\n", fd, size);
    assertf(res == 0, "Failed to close fd %d\n", fd);


	STDLOG(1,"Exiting LoadDerivatives z %d at z_file %d\n", z, z_file);
	
	//ReadDerivatives.Stop();
}
 

/* =========================  MPI Multipoles ================== */


/// Launch the CPD MPI_Irecv jobs.
/// Call this when the first slab is being finished.
/// We need to know the slabs that will be computed on each node.
void ParallelConvolution::RecvMultipoleSlab(int first_slab_finished) {
	
	STDLOG(1,"Setting MPI_Irecv for incoming Multipoles\n");
	
	STDLOG(1,"finished: %d on node with first slab %d \n", first_slab_finished, first_slab_on_node);
	
	
	int offset = first_slab_finished - first_slab_on_node;
	for (int r = 0; r < MPI_size; r++) {
		
	    for (int xraw = first_slabs_all[r] + offset; xraw < first_slabs_all[r] + offset + total_slabs_all[r]; xraw++) {
			int x = CP->WrapSlab(xraw);
			int tag = x + M_TAG;
			
			STDLOG(1, "recv tag = %d, x= %d\n", tag, x); 
			// We're receiving the full range of z's in each transfer.
			// Key thing is to route this to the correct x location
			
			STDLOG(1, "%p\n", MTdisk + x * this_node_size); 
			
			MPI_Irecv(MTdisk + x * this_node_size, this_node_size,
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
	
	STDLOG(1,"slab %d wrapped slab\n", slab);
	
    MTCOMPLEX * mt = (MTCOMPLEX *) SB->GetSlabPtr(MultipoleSlab, slab);
	
	STDLOG(1,"slab %d got slab ptr %p \n", slab, mt);
	
	Msend_requests[slab] = new MPI_Request[MPI_size];
	
	STDLOG(1,"slab %d init Msend_requests\n", slab);
	
	
	Msend_active++;
	for (int r = 0; r < MPI_size; r++) {
    	int tag = slab + M_TAG;
		
		STDLOG(1,"slab %d about to send %d, node_start[r] = %d, send tag  = %d\n", slab, r, node_start[r], tag);
				
		STDLOG(1, "%p\n", mt + node_start[r]); 
				
    	MPI_Isend(mt + node_start[r], node_size[r], MPI_COMPLEX, r, tag, MPI_COMM_WORLD, &Msend_requests[slab][r]);
		
		STDLOG(1,"slab %d sent %d\n", slab, r);
		
	}
	STDLOG(1,"Multipole slab %d has been queued for MPI_Isend, Msend_active = %d\n", slab, Msend_active);
}


// /// Return only once all of the multipole slabs have been received.
// /// We'll only check this at the end of the timestep loop, so
// /// this will just spin until done.
// void ParallelConvolution::WaitRecvMultipoleComplete(int slab) {
//
// 	STDLOG(1,"Entering WaitRecvMultipoleComplete for slab = %d\n", slab);
//
//
// 	while (1) //spin until done.
// 	{
// 		STDLOG(1,"WaitRecvMultipoleComplete spinning for slab = %d\n", slab);
// 		int err = MPI_Wait(&Mrecv_requests[slab], MPI_STATUS_IGNORE);
// 		//int err = MPI_Waitall(cpd, Mrecv_requests[slab], MPI_STATUSES_IGNORE);  // TODO: NAM okay to have MPI_STATUSES_IGNORE, or need array_of_statuses?
// 		STDLOG(1,"WaitRecvMultipoleComplete spinning, err = %d for slab = %d\n", err, slab);
//
// 		if (err == 0) return;
// 	}
// }

int ParallelConvolution::CheckRecvMultipoleComplete(int slab) {
	
	STDLOG(1,"Entering WaitRecvMultipoleComplete for slab = %d\n", slab);
	int done = 0; 
	
	STDLOG(1,"WaitRecvMultipoleComplete spinning for slab = %d\n", slab);
	int err = MPI_Wait(&Mrecv_requests[slab], MPI_STATUS_IGNORE);
	//int err = MPI_Waitall(cpd, Mrecv_requests[slab], MPI_STATUSES_IGNORE);  // TODO: NAM okay to have MPI_STATUSES_IGNORE, or need array_of_statuses? 
	STDLOG(1,"WaitRecvMultipoleComplete spinning, err = %d for slab = %d\n", err, slab);
	
	if (err == 0) done = 1;
	
	return done;  
}

/// Call periodically to check on the status of the MPI calls
/// and delete the MultipoleSlab when a send is complete.
/// Probably best if this is not expecting the completion to be 
/// in order, but we could discuss this.
int ParallelConvolution::CheckSendMultipoleComplete(int slab) {
	STDLOG(1,"Entering CheckSendMultipoleComplete for slab %d\n", slab);
	
    int x = slab; 	
    if (Msend_requests[x] == NULL) return 1;  // Nothing to do for this slab
	
	STDLOG(1,"2 CheckSendMultipoleComplete for slab %d\n", slab);
	
	
    // This slab has pending MPI calls
    int err, sent = 0, done = 1;
	
    for (int r = 0; r < MPI_size; r++) {		
		if (Msend_requests[x][r] == MPI_REQUEST_NULL) continue;  // Already done
		
		STDLOG(1,"4 %d CheckSendMultipoleComplete for slab %d\n", r, slab);
		
		
		err = MPI_Test(Msend_requests[x] + r, &sent, MPI_STATUS_IGNORE);
		STDLOG(1,"CheckSendMultipoleComplete, r= %d, err= %d, sent = %d\n", r, err, sent);
		
		if (not sent) {
		    // Found one that is not done
		    done=0; break;
		}
    }
    
	if (done) {
        delete[] Msend_requests[x];
		Msend_requests[x] = NULL;   // This will allow us to check faster
		Msend_active--;
		
		STDLOG(1,"MPI_Isend complete for MultipoleSlab %d, Msend_active = %d\n", x, Msend_active);
		
		// SB->DeAllocate(MultipoleSlab, x);
		// TODO: May need to force this to skip writing out?
    }
	
	return done; 
}

// /// Call after the timestep loop to spin until all of the
// /// MPI calls have been finished.
// void ParallelConvolution::WaitForMultipoleTransferComplete(int offset) {
// 	for (int _slab = first_slab_on_node + offset; _slab < first_slab_on_node + offset + total_slabs_on_node; _slab ++){
// 		int slab = _slab % cpd;
// 		STDLOG(1,"Begin waiting to determine that all incoming Multipoles received for slab %d, first %d, total %d\n", slab, first_slab_on_node, total_slabs_on_node);
// 		WaitRecvMultipoleComplete(slab); //this spins until recvs are done.
// 		STDLOG(1,"End waiting to determine that all incoming Multipoles received for slab %d\n", slab);
// 		STDLOG(1,"Begin waiting to determine that all outgoing Multipoles sent for slab %d\n", slab);
// 		int done = 0;
// 		while (not done)
// 			done = CheckSendMultipoleComplete(slab); //this spins until sends are done.
// 		STDLOG(1,"End waiting to determine that all outgoing Multipoles sent for slab %d\n", slab);
// 	}
// }

int ParallelConvolution::CheckForMultipoleTransferComplete(int _slab) {
	int slab = _slab % cpd;
	int done_recv = 0; 
	STDLOG(1,"Begin waiting to determine that all incoming Multipoles received for slab %d, first %d, total %d\n", slab, first_slab_on_node, total_slabs_on_node);
	done_recv = CheckRecvMultipoleComplete(slab); //this spins until recvs are done. 
	STDLOG(1,"End waiting to determine that all incoming Multipoles received for slab %d\n", slab);
	STDLOG(1,"Begin waiting to determine that all outgoing Multipoles sent for slab %d\n", slab);
	int done_send = 0; 
	done_send = CheckSendMultipoleComplete(slab); //this spins until sends are done. 
	STDLOG(1,"End waiting to determine that all outgoing Multipoles sent for slab %d\n", slab);
	
	return done_send and done_recv; 
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
	STDLOG(1,"Entering GetTaylorRecipient\n");
	
	for (int r = 0; r < MPI_size; r++){
		if (first_slabs_all[r] <= slab and slab < first_slabs_all[r] + total_slabs_all[r]){
			STDLOG(1, "For slab %d, identified node %d as Taylor target\n", slab, r);
			return r;
		} 
	}
}

void ParallelConvolution::SendTaylorSlab(int slab) {
	//Each node has Taylors for a limited range of z but for every x slab. 
	//figure out who the receipient should be based on x slab and send to them. 
	// Take from MTdisk. Set Tsend_requests[x] as request. 
	STDLOG(1,"Entering SendTaylorSlab\n");
	
	slab = CP->WrapSlab(slab);
	int r = GetTaylorRecipient(slab); //x-slab slab is in node r's domain. Send to node r. 
	
	int tag = slab + T_TAG; 
	STDLOG(1,"About to Isend SendTaylorSlab with tag %d to node %d for slab %d\n", tag, r, slab);
	
	MPI_Isend(MTdisk + slab * this_node_size, this_node_size, MPI_COMPLEX, r, tag, MPI_COMM_WORLD, &Tsend_requests[slab]);
	
	STDLOG(1,"Taylor slab %d has been queued for MPI_Isend to rank %d\n", slab, r);
}
	
/// Allocate space for this TaylorSlab, and trigger the MPI calls
/// to gather the information from other nodes.  To be used in FetchSlab.
void ParallelConvolution::RecvTaylorSlab(int slab) {	
	//We are receiving a specified x slab. For each node, receive its z packet for this x. For each node, use Trecv_requests[x][rank] as request. 
	
	STDLOG(1,"Entering RecvTaylorSlab\n");
	
	slab = CP->WrapSlab(slab);
    MTCOMPLEX *mt = (MTCOMPLEX *) SB->GetSlabPtr(TaylorSlab, slab);
	
	Trecv_requests[slab] = new MPI_Request[MPI_size];
	Trecv_active++; 
	
	for (int r = 0; r < MPI_size; r++) {
		int tag = slab + T_TAG; 
		
		STDLOG(1,"About to Irecv RecvTaylorSlab rank = %d, tag = %d\n",r, tag);
		
		MPI_Irecv(mt + node_start[r], node_size[r], MPI_COMPLEX, r, tag, MPI_COMM_WORLD, &Trecv_requests[slab][r]);
	}
	STDLOG(1,"MPI_Irecv set for incoming Taylors\n");
	return;
}

/// Check whether the MPI calls are complete for this TaylorSlab,
/// returning 1 if yes.  To be used in the TaylorForcePrecondition.
/// Also, have to be careful about the SendManifest.
int ParallelConvolution::CheckRecvTaylorComplete(int slab) {
	
	STDLOG(1,"Entering CheckRecvTaylorComplete for slab %d\n", slab);
	
	slab = CP->WrapSlab(slab);
	
    if (Trecv_requests[slab]==NULL) return 1;  // Nothing to do for this slab
	
    // This slab has pending MPI calls
    int err, sent=0, done=1;
	
    for (int r=0; r<MPI_size; r++) {
		if (Trecv_requests[slab][r]==MPI_REQUEST_NULL) continue;  // Already done
		
		err = MPI_Test(Trecv_requests[slab]+r, &sent,MPI_STATUS_IGNORE);
		STDLOG(1,"err = %d, sent = %d CheckRecvTaylorComplete\n", err, sent);
		
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
	STDLOG(1,"Entering Swizzle_to_zmxy\n");
	
  
	// TODO: Do the copies, including the -1 x offset
	// The convolve code expects the data from file x (MTdisk[x]) to be in MTzmxy[x+1]
	// Converting from [x][znode][m][y] to [znode][m][x][y]
	// Loop over the combination of [znode][m]
	#pragma omp parallel for schedule(static) 
	for (int zm = 0; zm < znode*rml; zm++) {
	    int zoffset = zm/rml;   // This is not actually z, but z-zstart
	    int m = zm-zoffset*rml; 
		for (int x=0; x<cpd; x++) {
			int xuse = x+1; if (xuse>=0) xuse=0;   
			//int xuse = x; 
			
			for (int y=0; y<cpd; y++) {
			    MTzmxy[(int64)zm*cpd2pad + xuse*cpd + y] = //MTdisk[x][zoffset*rml*cpd + m*cpd + y];
				       MTdisk[(int64)x*znode*rml*cpd + (int64)zoffset*rml*cpd + m*cpd + y];
				   }
	    }
	}

	STDLOG(1,"Exiting Swizzle_to_zmxy\n");
	
	return;
}

/// Copy from MTzmxy buffer back to MTdisk, then free buffer
void ParallelConvolution::Swizzle_to_xzmy() {
	// TODO: Do the copies, including the -1 x offset
	// TODO: Might prefer a different pragma structure than the other Swizzle
	// e.g., to combine the [x][znode] loop, so that each thread writes contiguous [m][y].
	STDLOG(1,"Entering Swizzle_to_xzmy\n");

	#pragma omp parallel for schedule(static)
	for (int xz = 0; xz < cpd * znode; xz ++) {
		int xoffset = xz / znode; 
		int z = xz - xoffset * znode; 
		
		//int x = xoffset;
		int x = xoffset+1; if (x>=0) x = 0;
		
		for (int m = 0; m < rml; m ++){
			for (int y = 0; y < cpd; y++){	 
				
#ifdef DO_NOTHING
				STDLOG(1, "%f %f %f %f\n", real(MTdisk[ (int64)((xz  * rml + m)) * cpd + y ]), real(MTzmxy[ z * cpd2pad * rml + m * cpd2pad + x * cpd + y ]/((Complex)cpd)), imag(MTdisk[ (int64)((xz  * rml + m)) * cpd + y ]), imag(MTzmxy[ z * cpd2pad * rml + m * cpd2pad + x * cpd + y ]/((Complex)cpd)));
										
#else
				MTdisk[ (int64)((xz  * rml + m)) * cpd + y ] = 
					MTzmxy[ z * cpd2pad * rml + m * cpd2pad + x * cpd + y ]/((Complex)cpd); //normalization constant for inverse fftw. 
 #endif	  
				
			}
		}
	}
#ifdef DO_NOTHING
	STDLOG(1, "Do nothing swizzle passed.\n");
#endif
	
	STDLOG(1,"Exiting Swizzle_to_xzmy\n");

	
    fftw_free(MTzmxy);
	
	STDLOG(1,"Exiting (freed) Swizzle_to_xzmy\n");
	
	return;
}



fftw_plan ParallelConvolution::PlanFFT(int sign){
	int err = fftw_init_threads(); 
	assertf(err != 0, "Failed to fftw_init_threads\n");
	fftw_plan_with_nthreads(omp_get_max_threads()); //NAM TODO how many threads in here? 
	
	fftw_plan plan;
	
	int  n[] = {cpd};	
	plan = fftw_plan_many_dft(1, n, cpd, //we will do znode*rml of these FFTS
		(fftw_complex *) MTzmxy, NULL, 1, cpd, //each new [x][y] chunk is located at strides of cpd. 
		(fftw_complex *) MTzmxy, NULL, 1, cpd,
		sign, FFTW_PATIENT);
		
	STDLOG(1,"Planned FFT with sign %d\n", sign);
	
	return plan;
		
}


/// Create the plan and do the FFT
void ParallelConvolution::FFT(fftw_plan plan) {
	// TODO: Import wisdom
	// We're going to plan to process each [x][y] slab separately,
	// doing the cpd FFTs of cpd length.
	// stride = cpd, dist = 1, nembed = NULL
				
	// Now we loop over the combination of [znode][m]
	#pragma omp parallel for schedule(static) //pragma appears okay 
	for (int zm = 0; zm < znode*rml; zm++) {
		    fftw_execute_dft(plan, (fftw_complex *)(MTzmxy+zm*cpd2pad),
				       (fftw_complex *)(MTzmxy+zm*cpd2pad));
	}
		
	fftw_destroy_plan(plan);

	STDLOG(1,"Exiting FFT\n");
	
	//fftw_cleanup_threads() NAM TODO do we need this?  or maybe in destructor? 
	return;
}



/// This executes the convolve on the MTfile info.
void ParallelConvolution::Convolve() {
	STDLOG(1,"Entering convolve\n");
	
	STDLOG(1,"order %d cpd %d\n", order, cpd);
	
	STDLOG(1,"blocksize %d\n", blocksize);
	
    InCoreConvolution *ICC = new InCoreConvolution(order,cpd,blocksize);
	
	STDLOG(1,"Made ICC convolve\n");
	
	// We're beginning in [x][znode][m][y] order on the RAMdisk
	// Copy to the desired order in a new buffer
	fftw_plan forward_plan  = PlanFFT(FFTW_FORWARD);
	fftw_plan backward_plan = PlanFFT(FFTW_BACKWARD);
		
	Swizzle_to_zmxy();   // Probably apply the -1 x offset here
	
	STDLOG(1, "Beginning forward FFT\n");
	FFT(forward_plan);
	
#ifndef DO_NOTHING
	
	
	

	// Two options: 
	// a) Swizzle to zmyx, then do a bunch of already-contiguous FFTs, and alter ICC
	// b) Swizzle to zmxy, then do a bunch of FFTs with a stride; use ICC as is
	// If option a, then the swizzle should be a smarter transpose.
	// Option b can be boring code.
	for (int z = 0; z < znode; z++) {
	    LoadDerivatives(z);
		
		//void InCoreConvolution::InCoreConvolve(Complex *FFTM, DFLOAT *CompressedD) {
		// InCoreConvolution(int order, int cpd, int blocksize) 
		//NAM TODO careful! does ICC know about doing one z at a time? 

	    ICC->InCoreConvolve(MTzmxy+z*rml*cpd2pad, Ddisk[z]); // TODO NAM add z argument. , z);
	}
	
	
	
#endif 
	STDLOG(1, "Beginning inverse FFT\n");
	FFT(backward_plan);
	Swizzle_to_xzmy();   // Restore to the RAMdisk
	
#ifdef DO_NOTHING
	exit(1);
#endif
	
	delete ICC;
	
    // if(OverwriteConvState){
//         RenameMultipolesToTaylors();
//     }
	
	STDLOG(1,"Leaving convolve\n");
	
}

// void ParallelConvolution::RenameMultipolesToTaylors(){
//     assert(OverwriteConvState);
//
// 	STDLOG(1, "Renaming multipole files to taylors\n");
//
//     char mfn[1024], tfn[1024];
//     for(int x = first_slab_on_node; x < first_slab_on_node + total_slabs_on_node; x++){
//         MultipoleFN(x % cpd, mfn);
//         TaylorFN(x % cpd, tfn);
//
//         int res = rename(mfn, tfn);
//         assertf(res == 0, "Failed to rename multipoles file \"%s\" to taylors file \"%s\".", mfn, tfn);
//     }
// }
//
//
