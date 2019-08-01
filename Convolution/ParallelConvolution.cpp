/* This class will handle the Multipoles and Taylors in a parallel
in-memory setting.  The timestep code should compute MultipoleSlab
and then hand it to this code.  Here, we trigger MPI_Issends() and
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

STimer  Constructor;
STimer  AllocateMT;
STimer  AllocateDerivs; 
STimer  QueueTaylors; 
STimer  FFTPlanning;
STimer 	ForwardZFFTMultipoles;
STimer  InverseZFFTTaylor; 
STimer  ArraySwizzle; 
STimer  ReadDerivatives;
STimer  ConvolutionArithmetic; 
STimer  Destructor; 
STimer  ThreadCleanUp; 

STimer  MmapMT; 
STimer  MunmapMT; 
STimer  MmapDerivs; 
STimer  MunmapDerivs; 

/// The zstart for the node of the given rank.
int ParallelConvolution::Zstart(int rank) {
    return (cpd2p1*rank/MPI_size); 
}

// ParallelConvolution::ParallelConvolution(){}
/// The constructor.  This should open the MTfile as read/write or allocate
/// space if it doesn't exist.
ParallelConvolution::ParallelConvolution(int _cpd, int _order, char MultipoleDirectory[1024], int _create_MT_file) { //TODO NAM does this always assume convolution overwrite mode? 
	
	STDLOG(2, "Calling constructor for ParallelConvolution\n");
	
	Constructor.Clear(); Constructor.Start(); 
	
    std::stringstream ss;
    std::string s;
	ss << MultipoleDirectory << "/Multipoles";
	s = ss.str();
	mt_file = s.c_str(); 
    STDLOG(2,"Working with MTfile %s\n", mt_file);

	//in previous timestep's finish action, multipoles were distributed to nodes such that each node now stores multipole moments for all x and its given range of z that it will convolve. Fetch these from shared memory.
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
	
	STDLOG(2, "Doing zstart = %d and znode = %d\n", zstart, znode);
	
	CompressedMultipoleLengthXY = ((1+P.cpd)*(3+P.cpd))/8;
	invcpd3 = (Complex) (pow(cpd * cpd * cpd, -1.0)); //NAM TODO get rid of Complex typecasting. 
		 
    int cml = ((order+1)*(order+2)*(order+3))/6;
    int nprocs = omp_get_max_threads();
	size_t cacherambytes = P.ConvolutionCacheSizeMB*(1024LL*1024LL);

/* Old code: better to count up than count down!
    blocksize = 0;
    for(blocksize=cpd*cpd;blocksize>=2;blocksize--) 
        if((cpd*cpd)%blocksize==0)
            if(nprocs*2.5*cml*blocksize*sizeof(Complex) < cacherambytes) break;
                // 2.5 = 2 Complex (mcache,tcache) 1 double dcache
*/
	
    int blocksize = 1;
    for (int b=2; b<P.cpd*P.cpd; b++) {
        if (nprocs*2.5*cml*b*sizeof(Complex)>=cacherambytes) break;
            // Too much memory, so stop looking
        if ((P.cpd*P.cpd)%b == 0) blocksize = b;  // Could use this value
    }
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
	STDLOG(2, "Making RD object with %d %f\n", direct, sdb);
	ReadDirect * RD_RDD = new ReadDirect(direct,sdb);
	
	Constructor.Stop(); 
	CS.Constructor = Constructor.Elapsed();
	
	AllocMT();
	AllocDerivs();
	
	CS.ops = 0;
	CS.ComputeCores = omp_get_max_threads();
}

/// The destructor.  This should close the MTfile. NAM: MTfile can be closed right after mapping, I think. 
/// Also can do any checks to make sure that all work was completed.
ParallelConvolution::~ParallelConvolution() {
	//NAM how to check work was completed, ideas:
	// assert all values of the four arrays below have been set to MPI_REQUEST_NULL? 
	// add class variable that keeps track of the number of FFTS, their success? check this? 
	
	Destructor.Clear(); Destructor.Start(); 
		
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
	
	MunmapMT.Clear(); MunmapMT.Start(); 
    int res = munmap(MTdisk, size);
    assertf(res == 0, "munmap failed\n");
	MunmapMT.Stop(); 
	CS.MunmapMT += MunmapMT.Elapsed(); 

	assertf(ramdisk_derivs == 0 || ramdisk_derivs == 1, "Ramdisk_derivs detected bad value: %d.\n", ramdisk_derivs);
	
    for (int z = 0; z < znode; z++){
	    if(ramdisk_derivs){
		    if(Ddisk[z] != NULL){
				MunmapDerivs.Clear(); MunmapDerivs.Start(); 
		        int res = munmap(Ddisk[z], size);
		        assertf(res == 0, "Failed to munmap derivs %d\n", z); 
				MunmapDerivs.Stop(); 
				CS.MunmapDerivs += MunmapDerivs.Elapsed(); 
		    }
		} else {
			free(Ddisk[z]); 
		}
	}
	
	delete[] Ddisk; 
	delete RD_RDD;
	
	Destructor.Stop(); CS.Destructor = Destructor.Elapsed();
	
	ThreadCleanUp.Clear(); ThreadCleanUp.Start();
	//fftw_cleanup_threads(); //NAM check usage. 	
	ThreadCleanUp.Stop(); CS.ThreadCleanUp = ThreadCleanUp.Elapsed();
	
	CS.ConvolveWallClock += Destructor.Elapsed() + ThreadCleanUp.Elapsed(); 
	
	
	STimer dumpstats_timer; 
	
	dumpstats_timer.Start(); 
	
	if (not create_MT_file){ //if this is not the 0th step, dump stats. 
	    char timingfn[1050];
	    sprintf(timingfn,"%s/last%s.convtime",P.LogDirectory,NodeString);
	    dumpstats(timingfn); 
	}
	
	dumpstats_timer.Stop(); 
	
	STDLOG(1, "Dumpstats took %f seconds\n", dumpstats_timer.Elapsed());
}




/* =========================   ALLOC BUFFERS   ================== */ 


void ParallelConvolution::AllocMT(){ 
	//MTdisk has shape [x][znode][m][y]
	
	//TODO this appears almost exactly as in arena allocator. Consider factoring out to ramdisk util. 
	AllocateMT.Clear(); AllocateMT.Start(); 
	
	
	assert(mt_file != NULL);
    assert(strnlen(mt_file,1) > 0);
    STDLOG(2,"Mapping MTdisk from shared memory\n");

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
	
	MmapMT.Clear(); MmapMT.Start(); 
	
	MTdisk = (MTCOMPLEX *) mmap(NULL, size, mmap_flags, MAP_SHARED, fd, 0);
    int res = close(fd);
    assertf((void *) MTdisk != MAP_FAILED, "mmap shared memory from fd = %d of size = %d failed\n", fd, size);
    assertf(MTdisk != NULL, "mmap shared memory from fd = %d of size = %d failed\n", fd, size);
    assertf(res == 0, "Failed to close fd %d\n", fd);
	
	STDLOG(2, "Successfully mapped MTdisk of size %d bytes at address %p.\n", size, MTdisk);
	
	MmapMT.Stop(); 
		
    size_t bufsize = sizeof(Complex)*znode*rml*cpd2pad; 	
	
	MTzmxy = (Complex *) fftw_malloc(bufsize);
	assertf(MTzmxy!=NULL, "Failed fftw_malloc for MTzmxy\n");
	STDLOG(2, "Successfully malloced MTzmxy at address %p\n", MTzmxy);
	
	AllocateMT.Stop(); 
	CS.AllocMT = AllocateMT.Elapsed();
	CS.MmapMT  = MmapMT.Elapsed();
	
	return;
}

void ParallelConvolution::AllocDerivs(){
	AllocateDerivs.Clear(); AllocateDerivs.Start();
	
    if(is_path_on_ramdisk(P.DerivativesDirectory)) ramdisk_derivs = 1;
	else ramdisk_derivs = 0; 
	assertf(ramdisk_derivs == 0 || ramdisk_derivs == 1, "Ramdisk_derivs detected bad value: %d.\n", ramdisk_derivs);
	
	Ddisk = new DFLOAT*[znode];
	for(int z = 0; z < znode; z++)
        Ddisk[z] = NULL;
	
    if(ramdisk_derivs){
		AllocateDerivs.Stop();
		CS.AllocDerivs = AllocateDerivs.Elapsed(); 
        return;
	}
	
    size_t s = sizeof(DFLOAT)*(rml*CompressedMultipoleLengthXY);
    for (int z = 0; z < znode; z++){
        int memalign_ret = posix_memalign((void **) (Ddisk + z), 4096, s);
        assert(memalign_ret == 0);
        //alloc_bytes += s;
    }
	
	CS.ReadDerivativesBytes = 0.0; 
	AllocateDerivs.Stop();
	CS.AllocDerivs = AllocateDerivs.Elapsed(); 
}

/* =========================   DERIVATIVES   ================== */ 

void ParallelConvolution::LoadDerivatives(int z) {
	int z_file = z + zstart;	
		
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
	
	assertf(ramdisk_derivs == 0 || ramdisk_derivs == 1, "Ramdisk_derivs detected bad value: %d.\n", ramdisk_derivs);
	
    if(!ramdisk_derivs){						
		assert(Ddisk[z] != NULL); 
        RD_RDD->BlockingRead( fn, (char *) Ddisk[z], size, 0, 1); //NAM TODO the last arg (1) should be ReadDirect's ramdisk flag set during constructor, but that seg faults right now. Fix! 
        CS.ReadDerivativesBytes += size;
		STDLOG(2, "BlockingRead derivatives for z = %d from file %c\n", z_file, fn);
		
		
    } else {
	    if(Ddisk[z] != NULL){ 
	        int res = munmap(Ddisk[z], size);
	        assertf(res == 0, "Failed to munmap derivs\n"); 
	    }

	    int fd = open(fn, O_RDWR, S_IRUSR | S_IWUSR);
	    assertf(fd != -1, "Failed to open shared memory file at \"%s\"\n", fn);

		MmapDerivs.Clear(); MmapDerivs.Start(); 
		
	    // map the shared memory fd to an address
	    Ddisk[z] = (DFLOAT *) mmap(NULL, size, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
	    int res = close(fd);
	    assertf((void *) Ddisk[z] != MAP_FAILED, "mmap shared memory from fd = %d of size = %d failed\n", fd, size);
	    assertf(Ddisk[z] != NULL, "mmap shared memory from fd = %d of size = %d failed\n", fd, size);
	    assertf(res == 0, "Failed to close fd %d\n", fd);
		STDLOG(2, "Mapped derivatives for z = %d from file %c\n", z_file, fn);
		
		MmapDerivs.Stop(); 
		CS.MmapDerivs += MmapDerivs.Elapsed(); 
    }
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
	STDLOG(2,"MPI_Irecv set for incoming Multipoles\n");
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
		int tag = (r+1) * 10000 + slab + M_TAG; 
		MPI_Issend(mt + node_start[r], node_size[r], MPI_COMPLEX, r, tag, MPI_COMM_WORLD, &Msend_requests[slab][r]);		
	}
	STDLOG(2,"Multipole slab %d has been queued for MPI_Issend, Msend_active = %d\n", slab, Msend_active);
}


int ParallelConvolution::CheckRecvMultipoleComplete(int slab) {
	int received = 0; 
	int err = MPI_Test(&Mrecv_requests[slab], &received, MPI_STATUS_IGNORE);	
	if (not received) STDLOG(4, "Multipole slab %d not received yet...\n", slab);	
	return received;  
}

/// Call periodically to check on the status of the MPI calls
/// and delete the MultipoleSlab when a send is complete.
/// Probably best if this is not expecting the completion to be 
/// in order, but we could discuss this.
int ParallelConvolution::CheckSendMultipoleComplete(int slab) {	
	
    int x = slab; 	
    if (Msend_requests[x] == NULL) return 1;  // Nothing to do for this slab

    // This slab has pending MPI calls
    int err = 0, sent = 0, done = 1;
	
    for (int r = 0; r < MPI_size; r++) {		
		if (Msend_requests[x][r]==MPI_REQUEST_NULL) continue;  // Already done
		err = MPI_Test(Msend_requests[x] + r, &sent, MPI_STATUS_IGNORE);
	
		if (not sent) {
		    // Found one that is not done
			STDLOG(4, "Multipole slab %d not sent yet...\n", slab);
			
		    done=0; break;
		}
		
    }
    
	if (done) {
        delete[] Msend_requests[x];
 		Msend_requests[x] = NULL;   // This will allow us to check faster

		Msend_active--;
		STDLOG(2,"MPI_Issend complete for MultipoleSlab %d, Msend_active = %d\n", x, Msend_active);
    }
	
	return done; 
}

int ParallelConvolution::CheckForMultipoleTransferComplete(int _slab) {
	int slab = _slab % cpd;
	int done_recv = 0; 
	done_recv = CheckRecvMultipoleComplete(slab); //this spins until recvs are done. 
	int done_send = 0; 
	done_send = CheckSendMultipoleComplete(slab); //this spins until sends are done. 
	if (done_send and done_recv) STDLOG(1,"Multipole MPI work complete for slab %d.\n", slab);
	
	return done_send and done_recv; 
}

/* ======================= MPI Taylors ========================= */


// TODO: Need to figure out how to pull the Taylors over.
// Need to be careful about the ordering of the Issend's.
// If we queue up all of the Issend's, will we flood the buffers and
// block the Manifest?  But we have no guarantee that the nodes will
// stay in timing agreement.  Do we need to provision to hold all of the 
// TaylorSlabs at once, as opposed to providing on request?
// Would having a separate thread to field requests be safer?  but multi-threading MPI
// Would using a different COMM for the Manifest help?
// If we put tags on every message, might it help to avoid flooding a buffer?

///for a given x slab = slab, figure out which node rank is responsible for it and should receive its Taylors. 
int ParallelConvolution::GetTaylorRecipient(int slab, int offset){	
	for (int r = 0; r < MPI_size; r++) {
		for (int _x = first_slabs_all[r] + offset; _x < first_slabs_all[r] + offset + total_slabs_all[r]; _x++) {
			int x = _x % cpd;
			if (x == slab){
				STDLOG(2, "For slab %d, identified node %d as Taylor target\n", slab, r);
				return r;
			}
		} 
	}
    QUIT("Did not find target node for slab %d!", slab);
    return -1;
}

void ParallelConvolution::SendTaylors(int offset) {
	QueueTaylors.Clear(); QueueTaylors.Start();
	for (int slab = 0; slab < cpd; slab ++){ 
		//Each node has Taylors for a limited range of z but for every x slab. 
		//figure out who the receipient should be based on x slab and send to them. 
		// Take from MTdisk. Set Tsend_requests[x] as request. 	
		slab = CP->WrapSlab(slab);
		STDLOG(2, "About to SendTaylor Slab %d with offset %d\n", slab, FORCE_RADIUS); 
		
		int r = GetTaylorRecipient(slab, offset); //x-slab slab is in node r's domain. Send to node r. 
	
		int tag = (MPI_rank+1) * 10000 + slab + T_TAG; 
		
		MPI_Issend(MTdisk + slab * this_node_size, this_node_size, MPI_COMPLEX, r, tag, MPI_COMM_WORLD, &Tsend_requests[slab]);
		STDLOG(2,"Taylor slab %d has been queued for MPI_Issend to rank %d\n", slab, r);
	}
	
	QueueTaylors.Stop(); CS.SendTaylors = QueueTaylors.Elapsed(); 
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
	STDLOG(2,"MPI_Irecv set for incoming Taylors\n");
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
			STDLOG(4, "Taylor slab %d not received yet...\n", slab);
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
	
	if (not sent) STDLOG(4, "Taylor slab %d not sent yet...\n", slab);
	
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
	
	if(send_done and recv_done) STDLOG(1, "Taylor MPI work complete for slab %d.\n", slab); 

	return (send_done and recv_done); 

}



/* =========================  Convolve functions ================== */

/// Create the MTzmxy buffer and copy into it from MTdisk.
// Note that this is also a type-casting from float to double.
// I wonder if this is why this step has been slow?
void ParallelConvolution::Swizzle_to_zmxy() {
	STDLOG(2,"Entering Swizzle_to_zmxy\n");
	  
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

	STDLOG(2,"Exiting Swizzle_to_zmxy\n");
	
	return;
}

/// Copy from MTzmxy buffer back to MTdisk, then free buffer
void ParallelConvolution::Swizzle_to_xzmy() {
	STDLOG(2,"Swizzling to xzmy\n");

	#pragma omp parallel for schedule(static)
	for (int xz = 0; xz < cpd * znode; xz ++) {
		int xoffset = xz / znode; 
		int z = xz - xoffset * znode; 
		
		int x = xoffset + 1;
		if (x>=cpd) x = 0;
		
		for (int m = 0; m < rml; m ++){
			for (int y = 0; y < cpd; y++){	 
				MTdisk[ (int64)((xz  * rml + m)) * cpd + y ] = 
					MTzmxy[ z * cpd2pad * rml + m * cpd2pad + x * cpd + y ] * invcpd3; //one factor of 1/CPD comes from the FFT below. Two more are needed for the FFTs in slab taylor -- we choose to throw them in here. 
			}
		}
	}

    fftw_free(MTzmxy);
	return;
}



fftw_plan ParallelConvolution::PlanFFT(int sign ){
	char wisdom_file[sizeof "parallel_convolve_fft_1.wisdom"];
	sprintf(wisdom_file, "parallel_convolve_fft_%d.wisdom", sign);
	int wisdom_exists = fftw_import_wisdom_from_filename(wisdom_file);
	
	STDLOG(1, "Wisdom import returned %d. 1 indicates successful import, 0 otherwise.\n", wisdom_exists);
	
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
		
	if (wisdom_exists == 0) {
		STDLOG(1, "Exporting wisdom to file.\n");
		fftw_export_wisdom_to_filename(wisdom_file);
	}	
		
		
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
	FFTPlanning.Clear(); FFTPlanning.Start();
		int err = fftw_init_threads(); 
		assertf(err != 0, "Failed to fftw_init_threads\n");
		fftw_plan_with_nthreads(omp_get_max_threads()); //NAM TODO how many threads in here? 
	
		fftw_plan forward_plan  = PlanFFT(FFTW_FORWARD);
		fftw_plan backward_plan = PlanFFT(FFTW_BACKWARD);
	FFTPlanning.Stop();
	CS.FFTPlanning = FFTPlanning.Elapsed();	
	
		
	ArraySwizzle.Clear(); ArraySwizzle.Start(); 	
		Swizzle_to_zmxy();   // Probably apply the -1 x offset here
	ArraySwizzle.Stop(); 
	
	STDLOG(2, "Beginning forward FFT\n");
	ForwardZFFTMultipoles.Clear(); 
	ForwardZFFTMultipoles.Start(); 
		FFT(forward_plan);
	ForwardZFFTMultipoles.Stop(); 
	CS.ForwardZFFTMultipoles = ForwardZFFTMultipoles.Elapsed();
	
	ConvolutionArithmetic.Clear(); 
	ReadDerivatives.Clear(); 
	//Swizzle to zmxy, then do a bunch of FFTs with a stride; use ICC as is
	for (int z = 0; z < znode; z++) {
		ReadDerivatives.Start();
	    	LoadDerivatives(z);
		ReadDerivatives.Stop(); 
		
		ConvolutionArithmetic.Start(); 
			ICC->InCoreConvolve(MTzmxy+z*rml*cpd2pad, Ddisk[z]); 
		ConvolutionArithmetic.Stop();
	}
	
	CS.ConvolutionArithmetic = ConvolutionArithmetic.Elapsed();
	CS.ReadDerivatives       =       ReadDerivatives.Elapsed(); 

	STDLOG(2, "Beginning inverse FFT\n");
	InverseZFFTTaylor.Clear(); 
	InverseZFFTTaylor.Start(); 
		FFT(backward_plan);
	InverseZFFTTaylor.Stop(); 
	CS.InverseZFFTTaylor = InverseZFFTTaylor.Elapsed();

	ArraySwizzle.Start();
		Swizzle_to_xzmy();   // Restore to the RAMdisk
	ArraySwizzle.Stop(); 
	CS.ArraySwizzle = ArraySwizzle.Elapsed(); 
	
	
	CS.ops = ICC->ConvolutionArithmeticCount();
	delete ICC;
	
	STDLOG(1,"Convolution complete.\n");
	
}



/* ======================== REPORTING ======================== */ 

void ParallelConvolution::dumpstats(char *fn) {

    FILE *fp;
    fp = fopen(fn,"w");
    assert(fp!=NULL);


    double accountedtime  = CS.ConvolutionArithmetic;
           accountedtime += CS.ForwardZFFTMultipoles + CS.InverseZFFTTaylor;
           accountedtime += CS.ReadDerivatives;
		   accountedtime += CS.Constructor + CS.AllocMT + CS.AllocDerivs + CS.SendTaylors + CS.FFTPlanning + CS.Destructor + CS.ThreadCleanUp; 
           accountedtime += CS.ArraySwizzle;
    double discrepency = CS.ConvolveWallClock - accountedtime;
	

    int computecores = CS.ComputeCores;
    fprintf(fp,"Convolution parameters:  RamAllocated = %dMB CacheSizeMB = %dMB nreal_cores=%d blocksize=%d zwidth=%d cpd=%d order=%d",
	         (int) (CS.totalMemoryAllocated/(1<<20)), P.ConvolutionCacheSizeMB, computecores, (int) blocksize, (int) zwidth, cpd, order);
			 

     fprintf(fp,"\n\n");
	 
     fprintf(fp,"\t ConvolutionWallClock:  %2.2e seconds \n", CS.ConvolveWallClock );
     fprintf(fp,"\t \t %50s : %1.2e seconds\n", "Constructor", CS.Constructor );
     fprintf(fp,"\t \t %50s : %1.2e seconds\n", "Allocating Multipole/Taylor", CS.AllocMT);
	 
     fprintf(fp,"\t \t %50s : %1.2e seconds\n", "Mmaping Multipole/Taylor (part of alloc)", CS.MmapMT);
     fprintf(fp,"\t \t %50s : %1.2e seconds\n", "Munmaping Multipole/Taylor (part of alloc)", CS.MunmapMT);
	 
     fprintf(fp,"\t \t %50s : %1.2e seconds\n", "Allocating Derivatives", CS.AllocDerivs);	 
	 
     fprintf(fp,"\t \t %50s : %1.2e seconds\n", "Mmaping Derivatives (part of alloc)", CS.MmapDerivs);
     fprintf(fp,"\t \t %50s : %1.2e seconds\n", "Munmaping Derivatives (part of alloc)", CS.MunmapDerivs);
	 
     fprintf(fp,"\t \t %50s : %1.2e seconds\n", "Array Swizzling", CS.ArraySwizzle );

     double e = CS.ReadDerivativesBytes/CS.ReadDerivatives/(1.0e+6);
     fprintf(fp,"\t \t %50s : %1.2e seconds --> rate was %4.0f MB/s\n", "ReadDiskDerivatives", CS.ReadDerivatives, e );


     double Gops = ((double) CS.ops)/(1.0e+9);
     fprintf(fp,"\t \t %50s : %1.2e seconds for %5.3f billion double precision operations\n", "Convolution Arithmetic", CS.ConvolutionArithmetic, Gops );
	 
     fprintf(fp,"\t \t %50s : %1.2e seconds\n", "FFT Planning", CS.FFTPlanning );
     fprintf(fp,"\t \t %50s : %1.2e seconds\n", "Forward FFT Z Multipoles", CS.ForwardZFFTMultipoles );
     fprintf(fp,"\t \t %50s : %1.2e seconds\n", "Inverse FFT Z Taylor",         CS.InverseZFFTTaylor );
	 
     fprintf(fp,"\t \t %50s : %1.2e seconds\n", "Queuing MPI Issends for Taylors", CS.SendTaylors );
	 
     fprintf(fp,"\t \t %50s : %1.2e seconds\n", "Destructor", CS.Destructor );
     fprintf(fp,"\t \t %50s : %1.2e seconds\n", "fftw thread clean up (part of destructor)", CS.ThreadCleanUp );
	 
	 
	 	 

     fprintf(fp,"\t %50s : %1.2e seconds which is %d%% \n", "Unaccounted remaining wallclock time", discrepency, (int) (discrepency/CS.ConvolveWallClock*100) );

     double cae       = CS.ConvolutionArithmetic;
     double farithp   = cae/CS.ConvolveWallClock*100;
	 double setupp    = (CS.Constructor + CS.Destructor + CS.AllocMT + CS.AllocDerivs)/CS.ConvolveWallClock*100;
	 double planp     = (CS.FFTPlanning)/CS.ConvolveWallClock*100;
	 double sendp     = (CS.SendTaylors)/CS.ConvolveWallClock*100;
	 double ffftp     = (CS.ForwardZFFTMultipoles + CS.InverseZFFTTaylor)/CS.ConvolveWallClock*100;
     double fiop      = (CS.ReadDerivatives)/CS.ConvolveWallClock*100;
	 double fclean    = (CS.ThreadCleanUp)/CS.ConvolveWallClock*100;

     double swzp      = CS.ArraySwizzle/CS.ConvolveWallClock*100;

     fprintf(fp,"\n \t Summary: Fourier Transforms = %2.0f%%     Convolution Arithmetic = %2.0f%%     Array Swizzle = %2.0f%%    Disk IO = %2.0f%%     ", ffftp, farithp, swzp, fiop );

     fprintf(fp,"\n \t Set up (con/destructor + allocs/mmaps) = %2.0f%%     FFT planning = %2.0f%%     FFT threan clean up = %2.0f%%   Queuing Taylor Issends = %2.0f%%    \n", setupp, planp, fclean, sendp);
     fprintf(fp,"\t          Arithmetic rate = %2.0f DGOPS --> rate per core = %1.1f DGOPS\n", Gops/cae, Gops/cae/computecores );
     fprintf(fp,"\t          [DGOPS == Double Precision Billion operations per second]\n");
     fprintf(fp,"\n");
     
    fclose(fp);
}
