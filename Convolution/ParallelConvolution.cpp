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
    return (cpd2p1*rank/MPI_size); 
}

// ParallelConvolution::ParallelConvolution(){}
/// The constructor.  This should open the MTfile as read/write or allocate
/// space if it doesn't exist.
ParallelConvolution::ParallelConvolution(int _cpd, int _order, const char *MTfile, int _create_MT_file) { //TODO NAM does this always assume convolution overwrite mode? 
	
	STDLOG(1, "Calling constructor for ParallelConvolution\n");

    //StripeConvState = strcmp(P.Conv_IOMode, "stripe") == 0;
    //OverwriteConvState = strcmp(P.Conv_IOMode, "overwrite") == 0;	
	
	mt_file = MTfile;
	create_MT_file = _create_MT_file; 
	cpd = _cpd;
	cpd2p1 = (cpd+1)/2; 
	
	int pad = CPD2PAD/sizeof(Complex);
	cpd2pad = floor((cpd*cpd+pad)/pad)*pad;
	
	order = _order;
	rml = (order+1)*(order+1);
	zstart = Zstart(MPI_rank);
	znode = Zstart(MPI_rank+1)-zstart;
	this_node_size = znode*rml*cpd;
	
	STDLOG(1, "Doing zstart = %d and znode = %d\n", zstart, znode);
	
	CompressedMultipoleLengthXY = ((1+P.cpd)*(3+P.cpd))/8;
	invcpd3 = (Complex) (pow(cpd * cpd * cpd, -1.0)); //NAM TODO get rid of Complex typecasting. 
		 
    int cml = ((order+1)*(order+2)*(order+3))/6;
    int nprocs = omp_get_max_threads();
	size_t cacherambytes = P.ConvolutionCacheSizeMB*(1024LL*1024LL);

    blocksize = 0;
    for(blocksize=cpd*cpd;blocksize>=2;blocksize--) 
        if((cpd*cpd)%blocksize==0)
            if(nprocs*2.5*cml*blocksize*sizeof(Complex) < cacherambytes) break;
                // 2.5 = 2 Complex (mcache,tcache) 1 double dcache
	
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
	
	int DIOBufferSizeKB = 1LL<<11;
    size_t sdb = DIOBufferSizeKB;
    sdb *= 1024LLU;

    int direct = P.RamDisk; 
	assert( (direct == 1) || (direct==0) );
	STDLOG(1, "Making RD object with %d %f\n", direct, sdb);
	ReadDirect * RD_RDD = new ReadDirect(direct,sdb);
	
	AllocMT();
	AllocDerivs();
}

/// The destructor.  This should close the MTfile. NAM: MTfile can be closed right after mapping, I think. 
/// Also can do any checks to make sure that all work was completed.
ParallelConvolution::~ParallelConvolution() {
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
	    if(ramdisk_derivs){
		    if(Ddisk[z] != NULL){
		        int res = munmap(Ddisk[z], size);
		        assertf(res == 0, "Failed to munmap derivs %d\n", z); 
		    }
		} else {
			free(Ddisk[z]); 
		}
	}
	
	delete[] Ddisk; 
	delete RD_RDD;
	fftw_cleanup_threads(); //NAM check usage. 	
}



/* =========================   ALLOC BUFFERS   ================== */ 


void ParallelConvolution::AllocMT(){ 
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
	
	STDLOG(1, "Successfully mapped MTdisk of size %d bytes at address %p.\n", size, MTdisk);
		
    size_t bufsize = sizeof(Complex)*znode*rml*cpd2pad; 	
	
	MTzmxy = (Complex *) fftw_malloc(bufsize);
	assertf(MTzmxy!=NULL, "Failed fftw_malloc for MTzmxy\n");
	STDLOG(1, "Successfully malloced MTzmxy at address %p\n", MTzmxy);
	return;
}

void ParallelConvolution::AllocDerivs(){
    if(is_path_on_ramdisk(P.DerivativesDirectory))
        ramdisk_derivs = 1;
	
	Ddisk = new DFLOAT*[znode];
	for(int z = 0; z < znode; z++)
        Ddisk[z] = NULL;
	
    if(ramdisk_derivs)
        return;
		
    size_t s = sizeof(DFLOAT)*(rml*CompressedMultipoleLengthXY);
    for (int z = 0; z < znode; z++){
        int memalign_ret = posix_memalign((void **) (Ddisk + z), 4096, s);
        assert(memalign_ret == 0);
        //alloc_bytes += s;
    }
}

/* =========================   DERIVATIVES   ================== */ 

void ParallelConvolution::LoadDerivatives(int z) {
	int z_file = z + zstart;	
    //ReadDerivatives.Start(); 
	
	size_t size = sizeof(DFLOAT)*rml*CompressedMultipoleLengthXY;

    const char *fnfmt;
    if(sizeof(DFLOAT) == sizeof(float))
        fnfmt = "%s/fourierspace_float32_%d_%d_%d_%d_%d";
    else
        fnfmt = "%s/fourierspace_%d_%d_%d_%d_%d";
	
    // note the derivatives are stored in z-slabs, not x-slabs
    char fn[1024]; //abacus.py has code to copy Derivs to DerivativesDirectory, and this looks for them there. Change DerivsDirec = RAMDISK?
    sprintf(fn, fnfmt,
            P.DerivativesDirectory, 
            (int) cpd, order, P.NearFieldRadius,
            P.DerivativeExpansionRadius, z_file);

    if(!ramdisk_derivs){						
		assert(Ddisk[z] != NULL); 
        RD_RDD->BlockingRead( fn, (char *) Ddisk[z], size, 0, 1); //NAM TODO the last arg (1) should be ReadDirect's ramdisk flag set during constructor, but that seg faults right now. Fix! 
		
		        // ReadDerivativeBytes += size;
		
    } else {
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
    }

	
	STDLOG(1, "Mapped derivatives for z = %d from file %c\n", z_file, fn);
	//ReadDerivatives.Stop();
}
 

/* =========================  MPI Multipoles ================== */


/// Launch the CPD MPI_Irecv jobs.
/// Call this when the first slab is being finished.
/// We need to know the slabs that will be computed on each node.
void ParallelConvolution::RecvMultipoleSlab(int first_slab_finished) {
	int offset = first_slab_finished - first_slab_on_node;
	for (int r = 0; r < MPI_size; r++) {
		
	    for (int xraw = first_slabs_all[r] + offset; xraw < first_slabs_all[r] + offset + total_slabs_all[r]; xraw++) {
			int x = CP->WrapSlab(xraw);
			int tag = (MPI_rank+1) * 10000 + x + M_TAG; 
			
			// We're receiving the full range of z's in each transfer.
			// Key thing is to route this to the correct x location
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
	
	MTCOMPLEX * mt = (MTCOMPLEX *) SB->GetSlabPtr(MultipoleSlab, slab);
	
#ifdef DO_NOTHING	
	//fill up multipole slab with predictable values and pass them through MPI piping without convolving.
	int zcpd = (cpd+1)/2;
	for (int z = 0; z < zcpd; z++){
		for (int m = 0; m < rml; m ++){
			for (int y = 0; y < cpd; y ++){ 
				int i = z*cpd*rml + m*cpd + y; 
				MTCOMPLEX value = (slab+1)*1000000 + y * 1000 + z + 0.01 * m ; 
				mt[i] = value; 				
			}
		}
	}
#endif	
	
	Msend_requests[slab] = new MPI_Request[MPI_size];
	Msend_active++;
	for (int r = 0; r < MPI_size; r++) {
		int tag = (r+1) * 10000 + slab + M_TAG; 
		MPI_Isend(mt + node_start[r], node_size[r], MPI_COMPLEX, r, tag, MPI_COMM_WORLD, &Msend_requests[slab][r]);		
	}
	STDLOG(1,"Multipole slab %d has been queued for MPI_Isend, Msend_active = %d\n", slab, Msend_active);
}


int ParallelConvolution::CheckRecvMultipoleComplete(int slab) {
	STDLOG(1,"CheckRecvMultipoleComplete testing for slab = %d\n", slab);
	int received = 0; 
	int err = MPI_Test(&Mrecv_requests[slab], &received, MPI_STATUS_IGNORE);
	STDLOG(1,"CheckRecvMultipoleComplete tested, err = %d  received = %d for slab = %d\n", err, received, slab);
		
	return received;  
}

/// Call periodically to check on the status of the MPI calls
/// and delete the MultipoleSlab when a send is complete.
/// Probably best if this is not expecting the completion to be 
/// in order, but we could discuss this.
int ParallelConvolution::CheckSendMultipoleComplete(int slab) {	
	STDLOG(1,"Entering CheckSendMultipoleComplete for slab %d\n", slab);
	
    int x = slab; 	
    if (Msend_requests[x] == NULL) return 1;  // Nothing to do for this slab

    // This slab has pending MPI calls
    int err = 0, sent = 0, done = 1;
	
    for (int r = 0; r < MPI_size; r++) {		
		if (Msend_requests[x][r]==MPI_REQUEST_NULL) continue;  // Already done
		err = MPI_Test(Msend_requests[x] + r, &sent, MPI_STATUS_IGNORE);
	
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
    }
	
	return done; 
}

int ParallelConvolution::CheckForMultipoleTransferComplete(int _slab) {
	int slab = _slab % cpd;
	int done_recv = 0; 
	done_recv = CheckRecvMultipoleComplete(slab); //this spins until recvs are done. 
	int done_send = 0; 
	done_send = CheckSendMultipoleComplete(slab); //this spins until sends are done. 
	STDLOG(1,"Multipoles for slab %d: send_done %d, recv_done %d\n", slab, done_send, done_recv);
	
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
int ParallelConvolution::GetTaylorRecipient(int slab, int offset){
	STDLOG(1,"Entering GetTaylorRecipient\n");
	
	for (int r = 0; r < MPI_size; r++) {
		for (int _x = first_slabs_all[r] + offset; _x < first_slabs_all[r] + offset + total_slabs_all[r]; _x++) {
			int x = _x % cpd;
			if (x == slab){
				STDLOG(1, "For slab %d, identified node %d as Taylor target\n", slab, r);
				return r;
			}
		} 
	}
}

void ParallelConvolution::SendTaylorSlab(int slab, int offset) {
	//Each node has Taylors for a limited range of z but for every x slab. 
	//figure out who the receipient should be based on x slab and send to them. 
	// Take from MTdisk. Set Tsend_requests[x] as request. 	
	slab = CP->WrapSlab(slab);
	int r = GetTaylorRecipient(slab, offset); //x-slab slab is in node r's domain. Send to node r. 
	
	int tag = (MPI_rank+1) * 10000 + slab + T_TAG; 
	
#ifdef DO_NOTHING
	for (int z = 0; z < znode; z++){
		for(int m = 0; m < rml; m ++){ ;
			for(int y = 0; y < cpd; y++){
				
				int i = z*cpd*rml + m*cpd + y;
				MTCOMPLEX expected = (slab+1)*1000000 + y * 1000 + (z+zstart) + 0.01 * m ; 
				MTCOMPLEX observed = MTdisk[slab*this_node_size + i]; 
				if (abs(real(expected - observed)) > 0.00001 or abs(imag(expected - observed)) > 0.00001){
					STDLOG(1, "ERROR SENDING MTDISK SLAB %d, i %d, z %d, m %d, y %d, obs = %1.7f %1.7f, exp = %1.7f %1.7f\n", slab, i, z, m, y, real(MTdisk[slab*this_node_size + i]), imag(MTdisk[slab*this_node_size + i]), real(expected), imag(expected));
					assertf(0 == 99, "Failed piping test during SendTaylorSlab!\n");		
				}
			}
		}
	}

#endif
	
	
	MPI_Isend(MTdisk + slab * this_node_size, this_node_size, MPI_COMPLEX, r, tag, MPI_COMM_WORLD, &Tsend_requests[slab]);
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
		int tag = (r+1) * 10000 + slab + T_TAG; 
		MPI_Irecv(mt + node_start[r], node_size[r], MPI_COMPLEX, r, tag, MPI_COMM_WORLD, &Trecv_requests[slab][r]);
		
		
	}
	STDLOG(1,"MPI_Irecv set for incoming Taylors\n");
	return;
}

int ParallelConvolution::CheckTaylorRecvReady(int slab){
    if (Trecv_requests[slab]==NULL) return 1;  // Nothing to do for this slab
	
    int err, received=0, done=1;
    for (int r=0; r<MPI_size; r++) {
		if (Trecv_requests[slab][r]==MPI_REQUEST_NULL) continue;  // Already done
		
		received = 0 ;
		err = MPI_Test(&Trecv_requests[slab][r], &received, MPI_STATUS_IGNORE);
		
		if (not received) {
		    // Found one that is not done
		    done=0; break;
		}		
    }
	
	if (done) {
		delete[] Trecv_requests[slab];
		Trecv_requests[slab] = NULL; 
		Trecv_active--; 
	}
	
	return done;
}

int ParallelConvolution::CheckTaylorSendComplete(int slab){
	int sent = 0; 
	int err = MPI_Test(&Tsend_requests[slab], &sent, MPI_STATUS_IGNORE);
	return sent; 
}

/// Check whether the MPI calls are complete for this TaylorSlab,
/// returning 1 if yes.  To be used in the TaylorForcePrecondition.
/// Also, have to be careful about the SendManifest.
int ParallelConvolution::CheckTaylorSlabReady(int slab) {
	slab = CP->WrapSlab(slab);	
	int send_done = 0; int recv_done = 0; 
	send_done = CheckTaylorSendComplete(slab); 
	
	if (send_done) {
		recv_done = CheckTaylorRecvReady(slab);
	}
	
	STDLOG(1, "Taylor slab %d send_done %d, recv_done %d\n", slab, send_done, recv_done); 
	
#ifdef DO_NOTHING
	if (send_done and recv_done){
	    MTCOMPLEX *mt = (MTCOMPLEX *) SB->GetSlabPtr(TaylorSlab, slab);
		
		int zcpd = (P.cpd+1)/2;
		int oops = 0; 
		MTCOMPLEX expected; 
		MTCOMPLEX observed; 
	
		for (int z = 0; z < zcpd; z++){
			for (int m = 0; m < rml; m ++){
				for (int y = 0; y < cpd; y ++){ //lcpd*(lcpd+1)/2*rml*sizeof(MTCOMPLEX)
				
					int i = z*cpd*rml + m*cpd + y; 
			
					MTCOMPLEX value = (slab+1)*1000000 + y * 1000 + z + 0.01 * m ; 
			
					if (oops == 0  and (abs(real(value - mt[i])) > 0.00001 or abs(imag(value - mt[i])) > 0.00001) ){
						if (oops) STDLOG(1, "Slab %d, i %d y %d z %d m %d received %1.7f %1.7f, expected %1.7f %1.7f\n", 
								slab, i, y, z, m, real(observed), imag(observed), real(expected), imag(expected));
						assertf(0 == 99, "Failed piping test during CheckTaylorSlabReady!\n");		
					}
					mt[i] = 0.0; 		
				}
			}
		}
	}
	

#endif
	
	
	return (send_done and recv_done); 

}



/* =========================  Convolve functions ================== */

/// Create the MTzmxy buffer and copy into it from MTdisk.
// Note that this is also a type-casting from float to double.
// I wonder if this is why this step has been slow?
void ParallelConvolution::Swizzle_to_zmxy() {
	STDLOG(1,"Entering Swizzle_to_zmxy\n");
	  
	// The convolve code expects the data from file x (MTdisk[x]) to be in MTzmxy[x+1]
	// Converting from [x][znode][m][y] to [znode][m][x][y]
	// Loop over the combination of [znode][m]
	#pragma omp parallel for schedule(static) 
	for (int zm = 0; zm < znode*rml; zm++) {
	    int zoffset = zm/rml;   // This is not actually z, but z-zstart
	    int m = zm-zoffset*rml; 
		for (int x=0; x<cpd; x++) {
			int xuse = x + 1;
			if (xuse>=cpd) xuse=0;
			
			for (int y=0; y<cpd; y++) {
			    MTzmxy[(int64)zm*cpd2pad + xuse*cpd + y] =  
				       MTdisk[(int64)x*znode*rml*cpd + (int64)zoffset*rml*cpd + m*cpd + y] ;
				   }
	    }
	}

	STDLOG(1,"Exiting Swizzle_to_zmxy\n");
	
	return;
}

/// Copy from MTzmxy buffer back to MTdisk, then free buffer
void ParallelConvolution::Swizzle_to_xzmy() {
	STDLOG(1,"Swizzling to xzmy\n");

	#pragma omp parallel for schedule(static)
	for (int xz = 0; xz < cpd * znode; xz ++) {
		int xoffset = xz / znode; 
		int z = xz - xoffset * znode; 
		
		int x = xoffset + 1;
		if (x>=cpd) x = 0;
		
		for (int m = 0; m < rml; m ++){
			for (int y = 0; y < cpd; y++){	 
				
#ifdef DO_NOTHING
				MTdisk[ (int64)((xz  * rml + m)) * cpd + y ] = 
					MTzmxy[ z * cpd2pad * rml + m * cpd2pad + x * cpd + y ] / ((Complex) cpd);	
#else
				MTdisk[ (int64)((xz  * rml + m)) * cpd + y ] = 
					MTzmxy[ z * cpd2pad * rml + m * cpd2pad + x * cpd + y ] * invcpd3; //one factor of 1/CPD comes from the FFT below. Two more are needed for the FFTs in slab taylor -- we choose to throw them in here. 
#endif 
			}
		}
	}

    fftw_free(MTzmxy);
	return;
}



fftw_plan ParallelConvolution::PlanFFT(int sign ){
	fftw_plan plan;
	
	// TODO: Import wisdom
	// We're going to plan to process each [x][y] slab separately,
	// doing the cpd FFTs of cpd length.
	// stride = cpd, dist = 1, nembed = NULL
	
	
	int  n[] = {cpd};	
	plan = fftw_plan_many_dft(1, n, cpd, //we will do znode*rml of sets of cpd FFTs.
		(fftw_complex *) MTzmxy, NULL, cpd, 1, //each new [x][y] chunk is located at strides of cpd.
		(fftw_complex *) MTzmxy, NULL, cpd, 1,
		sign, FFTW_PATIENT);
		

		
	return plan;
		
}


/// Create the plan and do the FFT
void ParallelConvolution::FFT(fftw_plan plan) {			
	// Now we loop over the combination of [znode][m]
	#pragma omp parallel for schedule(static) //pragma appears okay
	for (int zm = 0; zm < znode*rml; zm++) {
		fftw_execute_dft(plan, (fftw_complex *)(MTzmxy+zm*cpd2pad),
				       (fftw_complex *)(MTzmxy+zm*cpd2pad));
	}
	
		
	fftw_destroy_plan(plan);
	return;
}



/// This executes the convolve on the MTfile info.
void ParallelConvolution::Convolve() {
	STDLOG(1,"Beginning Convolution.\n");
	
    InCoreConvolution *ICC = new InCoreConvolution(order,cpd,blocksize, cpd2pad);

	// We're beginning in [x][znode][m][y] order on the RAMdisk
	// Copy to the desired order in a new buffer
	int err = fftw_init_threads(); 
	assertf(err != 0, "Failed to fftw_init_threads\n");
	fftw_plan_with_nthreads(omp_get_max_threads()); //NAM TODO how many threads in here? 
	
	fftw_plan forward_plan  = PlanFFT(FFTW_FORWARD);
	fftw_plan backward_plan = PlanFFT(FFTW_BACKWARD);
		
	Swizzle_to_zmxy();   // Probably apply the -1 x offset here
	
	STDLOG(1, "Beginning forward FFT\n");
	FFT(forward_plan);
	
#ifndef DO_NOTHING

	//Swizzle to zmxy, then do a bunch of FFTs with a stride; use ICC as is
	for (int z = 0; z < znode; z++) {
	    LoadDerivatives(z);
		ICC->InCoreConvolve(MTzmxy+z*rml*cpd2pad, Ddisk[z]); 
		
	}
	
#endif 
	STDLOG(1, "Beginning inverse FFT\n");
	FFT(backward_plan);
	Swizzle_to_xzmy();   // Restore to the RAMdisk
	
	delete ICC;
	
	STDLOG(1,"Convolution complete.\n");
	
}
