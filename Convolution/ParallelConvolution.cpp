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

2D:

With regards to the 2D code, ParallelConvolution is mostly agnostic.
Each node has already done the Y-Z FFT, with its multipole (sub)slabs
ordered as [kz][m][local_ky].  ky is split over nodes, but ParallelConvolution
doesn't care; it just needs to know this node is doing fewer transforms.

One bookkeepping difference is that the kz dimension in the 1D code
has length (CPD+1)/2, while in the 2D code it has length CPD. This is
because the z-FFT is done second in the 2D code.  We'll use the variables
`kz_size` and `node_ky_size` to track the relevant lengths.

*/

//#define DO_NOTHING 1 

#include "ParallelConvolution.h"

#include "mpi_limiter.h"

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

#define CONVTIMEBUFSIZE 65536
char *convtimebuffer;    // This has to be allocated and freed outside of this class

/// The zstart for the node of the given rank.
int ParallelConvolution::Zstart(int rank) {
    return kz_size*rank/MPI_size_x;
}

// ParallelConvolution::ParallelConvolution(){}
/// The constructor.  This should open the MTfile as read/write or allocate
/// space if it doesn't exist.
ParallelConvolution::ParallelConvolution(int _cpd, int _order, char MultipoleDirectory[1024], int _create_MT_file)
    : mpi_limiter(P.MPICallRateLimit_ms)

    { //TODO NAM does this always assume convolution overwrite mode? 
	
	STDLOG(2, "Calling constructor for ParallelConvolution\n");
	
	Constructor.Clear(); Constructor.Start(); 
	
	mt_file = MultipoleDirectory / "Multipoles"; 
    STDLOG(2,"Working with MTfile {}\n", mt_file);

	//in previous timestep's finish action, multipoles were distributed to nodes such that each node now stores multipole moments for all x and its given range of z that it will convolve. Fetch these from shared memory.
	create_MT_file = _create_MT_file; 
	cpd = _cpd;
	
	if(MPI_size_z > 1){
		node_ky_size = node_cpdp1half;
		kz_size = cpd;
	} else {
		node_ky_size = cpd;
		kz_size = (cpd+1)/2;
	}
	
	int pad = CPD2PADDING/sizeof(Complex);
	cpdky_pad = ((cpd*node_ky_size+pad)/pad)*pad;
	
	order = _order;
	rml = (order+1)*(order+1);
	zstart = Zstart(MPI_rank_x);
	znode = Zstart(MPI_rank_x+1)-zstart;

	node_slab_elem = znode*rml*node_ky_size;  // number of elements in each slab this node handles
    convtimebuffer = NULL;
	
	STDLOG(2, "Doing zstart = {:d} and znode = {:d}\n", zstart, znode);
	
	CompressedMultipoleLengthXY = ((1+cpd)*(3+cpd))/8;
	invcpd3 = (Complex) (pow(cpd * cpd * cpd, -1.0)); //NAM TODO get rid of Complex typecasting. 
		 
    int cml = ((order+1)*(order+2)*(order+3))/6;
    int nprocs = omp_get_max_threads();
	size_t cacherambytes   = P.ConvolutionCacheSizeMB*(1024LL*1024LL);
    size_t L1cacherambytes = P.ConvolutionL1CacheSizeMB*(1024LL*1024L);
	
    blocksize = 1;
    for (int b=2; b<cpd*node_ky_size; b++) {
        if (2.5*b*sizeof(Complex)>=L1cacherambytes) break;
            // Too much L1 memory: can't hold one example of D,M,T
        if (nprocs*2.5*cml*b*sizeof(Complex)>=cacherambytes) break;
            // Too much L3 memory, can't hold all D,M,T, so stop looking
        if ((cpd*node_ky_size)%b == 0) blocksize = b;  // Could use this value
    }
        // 1.5 = 1 Complex (mcache) 1 double dcache
        // 2.5 = 2 Complex (mcache,tcache) 1 double dcache
	
	first_slabs_all = new int[MPI_size_x];
    total_slabs_all = new int[MPI_size_x]; 
	
	ReadNodeSlabs(1, first_slabs_all, total_slabs_all);
	
	first_slab_on_node  = first_slabs_all[MPI_rank_x];
	total_slabs_on_node = total_slabs_all[MPI_rank_x];
	
	
	assertf(sizeof(fftw_complex)==sizeof(Complex), "Mismatch of complex type sizing between FFTW and Complex\n");
	
	// Here we compute the start and number of z's in each node
	node_start = new int[MPI_size_x];
	node_size  = new int[MPI_size_x];

	for (int j=0; j<MPI_size_x; j++) {
	    node_start[j] = Zstart(j)*rml*node_ky_size;
	    node_size[j]  = (Zstart(j+1)-Zstart(j))*rml*node_ky_size;
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

    TaylorSlabAllMPIDone = new int[cpd];
    MultipoleSlabAllMPIDone = new int[cpd];

	for (int x = 0; x < cpd; x++) {
		Msend_requests[x] = NULL;  // This indicates no pending req
		Mrecv_requests[x] = MPI_REQUEST_NULL;  // This indicates no pending req
		Trecv_requests[x] = NULL;  // This indicates no pending req
		Tsend_requests[x] = MPI_REQUEST_NULL;  // This indicates no pending req	

        TaylorSlabAllMPIDone[x] = 0;
        MultipoleSlabAllMPIDone[x] = 0;
	} 
    
    MsendTimer = new STimer[cpd];


	int DIOBufferSizeKB = 1LL<<11;
    size_t sdb = DIOBufferSizeKB;
    sdb *= 1024LLU;

    int no_dio = !P.AllowDirectIO;
	STDLOG(2, "Making RD object with no_dio={:d} bufsize={:d}\n", no_dio, sdb);
	RD = new ReadDirect(no_dio, sdb);
	
	Constructor.Stop(); 
	CS.Constructor = Constructor.Elapsed();
	
	CS.totalMemoryAllocated = 0;
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
		
	assertf(Trecv_active == 0, "Error: active taylor receive requests ! = 0, = {:d}\n", Trecv_active);
	assertf(Msend_active == 0, "Error: active multipole send requests ! = 0, = {:d}\n", Msend_active);
	
	delete[] Msend_requests;
	delete[] Trecv_requests; 
	delete[] Mrecv_requests; 
	delete[] Tsend_requests; 

    delete[] TaylorSlabAllMPIDone;
    delete[] MultipoleSlabAllMPIDone;
    delete[] MsendTimer;

	delete[] node_start; 
	delete[] node_size; 

	delete[] first_slabs_all;
	delete[] total_slabs_all;
	
	MunmapMT.Clear(); MunmapMT.Start(); 
    int res = munmap(MTdisk, MTdisk_bytes);
    assertf(res == 0, "munmap failed\n");
	MunmapMT.Stop(); 
	CS.MunmapMT += MunmapMT.Elapsed(); 

    for (int z = 0; z < znode; z++){
	    if(ramdisk_derivs){
		    if(Ddisk[z] != NULL){
				MunmapDerivs.Clear(); MunmapDerivs.Start(); 
		        int res = munmap(Ddisk[z], Ddisk_bytes);
		        assertf(res == 0, "Failed to munmap derivs {:d}\n", z); 
				MunmapDerivs.Stop(); 
				CS.MunmapDerivs += MunmapDerivs.Elapsed(); 
		    }
		} else {
			free(Ddisk[z]); 
		}
	}
	
	delete[] Ddisk; 
	delete RD;
	
	Destructor.Stop(); CS.Destructor = Destructor.Elapsed();
	
	ThreadCleanUp.Clear(); ThreadCleanUp.Start();
	ThreadCleanUp.Stop(); CS.ThreadCleanUp = ThreadCleanUp.Elapsed();
	
	CS.ConvolveWallClock += Destructor.Elapsed() + ThreadCleanUp.Elapsed(); 
	
	
	STimer dumpstats_timer; 
	
	dumpstats_timer.Start(); 
	
	if (!create_MT_file){ //if this is not the 0th step, dump stats. 
#ifndef DO_NOTHING
	    dumpstats(); 
#endif
	}
	
	dumpstats_timer.Stop(); 
	
	STDLOG(1, "Dumpstats took {:f} seconds\n", dumpstats_timer.Elapsed());
}




/* ==================   ALLOC BUFFERS   ================== */ 


void ParallelConvolution::AllocMT(){ 
	//MTdisk has shape [x][znode][m][y]
	
	//TODO this appears almost exactly as in arena allocator. Consider factoring out to ramdisk util. 
	AllocateMT.Clear(); AllocateMT.Start(); 
	
    assert(!mt_file.empty());
    STDLOG(2,"Mapping MTdisk from shared memory\n");

	if(create_MT_file)
		fs::remove(mt_file);
	
    int shm_fd_flags = create_MT_file ? (O_RDWR | O_CREAT | O_EXCL) : O_RDWR;
    int fd = open(mt_file.c_str(), shm_fd_flags, S_IRUSR | S_IWUSR);
    assertf(fd != -1, "Failed to open shared memory file at \"{}\"\n", mt_file);
    
	MTdisk_bytes = sizeof(MTCOMPLEX)*node_slab_elem*cpd;
	
    if(create_MT_file){
	    // set the size
	    int res = ftruncate(fd, MTdisk_bytes);
	    assertf(res == 0, "ftruncate on shared memory mt_file = {} to size = {:d} failed\n", mt_file, MTdisk_bytes);
    } else {
        // check that the size is as expected
        assertf(res == 0, "fstat on shared memory ramdisk_fn = {}\n", mt_file);
        assertf(fs::file_size(mt_file) == MTdisk_bytes, "Found shared memory size {:d}; expected {:d} (ramdisk_fn = {})\n", shmstat.st_size, MTdisk_bytes, mt_file);
    }
    // map the shared memory fd to an address
    int mmap_flags = create_MT_file ? (PROT_READ | PROT_WRITE) : (PROT_READ | PROT_WRITE);  // the same
	
	MmapMT.Clear(); MmapMT.Start(); 
	
	MTdisk = (MTCOMPLEX *) mmap(NULL, MTdisk_bytes, mmap_flags, MAP_SHARED, fd, 0);
    int res = close(fd);
    assertf((void *) MTdisk != MAP_FAILED, "mmap shared memory from fd = {:d} of size = {:d} failed\n", fd, MTdisk_bytes);
    assertf(MTdisk != NULL, "mmap shared memory from fd = {:d} of size = {:d} failed\n", fd, MTdisk_bytes);
    assertf(res == 0, "Failed to close fd {:d}\n", fd);
	
	STDLOG(2, "Successfully mapped MTdisk of size {:d} bytes at address {:p}.\n", MTdisk_bytes, MTdisk);
	
	MmapMT.Stop(); 
		
    size_t bufsize = sizeof(Complex)*znode*rml*cpdky_pad;
	
	MTzmxy = (Complex *) fftw_malloc(bufsize);
	assertf(MTzmxy!=NULL, "Failed fftw_malloc for MTzmxy of size {:d}\n", bufsize);
	STDLOG(2, "Successfully malloced MTzmxy at address {:p}\n", MTzmxy);
	
	AllocateMT.Stop(); 
	CS.AllocMT = AllocateMT.Elapsed();
	CS.MmapMT  = MmapMT.Elapsed();

	CS.totalMemoryAllocated += MTdisk_bytes + bufsize;
	
	return;
}

void ParallelConvolution::AllocDerivs(){
	AllocateDerivs.Clear(); AllocateDerivs.Start();

	if(MPI_size_z > 1){
		Ddisk_bytes = sizeof(DFLOAT)*rml*(cpd+1)/2*node_ky_size;
	} else {
		Ddisk_bytes = sizeof(DFLOAT)*rml*CompressedMultipoleLengthXY;
	}

	// 2D: deriv files are [ky][m][kx], read a subset of ky (kx has length cpdp1half)
	dfile_offset = MPI_size_z > 1 ? sizeof(DFLOAT)*rml*(cpd+1)/2*all_node_ky_start[MPI_rank_z] : 0;
	
    ramdisk_derivs = is_path_on_ramdisk(P.DerivativesDirectory) &&
						MPI_size_z == 1;  // TODO: need page padding for ramdisk mapping
	
	Ddisk = new DFLOAT*[znode];
	for(int z = 0; z < znode; z++)
        Ddisk[z] = NULL;
	
    if(ramdisk_derivs){
		AllocateDerivs.Stop();
		CS.AllocDerivs = AllocateDerivs.Elapsed(); 
        return;
	}
	
    for (int z = 0; z < znode; z++){
        assert(posix_memalign((void **) (Ddisk + z), PAGE_SIZE, Ddisk_bytes) == 0);
    }
	
	CS.ReadDerivativesBytes = 0.0; 
	AllocateDerivs.Stop();
	CS.AllocDerivs = AllocateDerivs.Elapsed(); 
	CS.totalMemoryAllocated += Ddisk_bytes;
}

/* ==================   DERIVATIVES   ================== */ 

void ParallelConvolution::LoadDerivatives(int z) {
	int z_file = z + zstart;	

    const char *f32str = sizeof(DFLOAT) == 4 ? "_float32" : "";
	const char *twoDstr = MPI_size_z > 1 ? "_2D" : "";
	const char *fnfmt = "%s/fourierspace%s%s_%d_%d_%d_%d_%d";
	
    // note the derivatives are stored in z-slabs, not x-slabs
    char fn[1024];
    sprintf(fn, fnfmt,
            P.DerivativesDirectory,
			f32str, twoDstr,
            (int) cpd, order, P.NearFieldRadius,
            P.DerivativeExpansionRadius, z_file);
	
	assertf(ramdisk_derivs == 0 || ramdisk_derivs == 1, "Ramdisk_derivs detected bad value: {:d}.\n", ramdisk_derivs);
	
    if(!ramdisk_derivs){
		assert(Ddisk[z] != NULL); 
        RD->BlockingRead( fn, (char *) Ddisk[z], Ddisk_bytes, dfile_offset, ramdisk_derivs);
        CS.ReadDerivativesBytes += Ddisk_bytes;
		STDLOG(2, "BlockingRead derivatives for z = {:d} from file {}\n", z_file, fn);
    } else {
	    if(Ddisk[z] != NULL){
	    	MunmapDerivs.Clear(); MunmapDerivs.Start();
	        int res = munmap(Ddisk[z], Ddisk_bytes);
	        MunmapDerivs.Stop();
	        CS.MunmapDerivs += MunmapDerivs.Elapsed();
	        assertf(res == 0, "Failed to munmap derivs\n");
	    }

	    int fd = open(fn.c_str(), O_RDWR, S_IRUSR | S_IWUSR);
	    assertf(fd != -1, "Failed to open shared memory file at \"{}\"\n", fn);

		MmapDerivs.Clear(); MmapDerivs.Start(); 
		
	    // map the shared memory fd to an address
		// TODO 2D: deriv files may need padding such that dfile_offset is page-aligned
		// Alternatively, map the whole thing (need to keep original pointers for munmap).
		// But that's slower.
	    Ddisk[z] = (DFLOAT *) mmap(NULL, Ddisk_bytes, PROT_READ | PROT_WRITE, MAP_SHARED, fd, dfile_offset);
	    int res = close(fd);
	    assertf((void *) Ddisk[z] != MAP_FAILED, "mmap shared memory from fd = {:d} of size = {:d} failed\n", fd, Ddisk_bytes);
	    assertf(Ddisk[z] != NULL, "mmap shared memory from fd = {:d} of size = {:d} failed\n", fd, Ddisk_bytes);
	    assertf(res == 0, "Failed to close fd {:d}\n", fd);
		STDLOG(2, "Mapped derivatives for z = {:d} from file {}\n", z_file, fn);
		
		MmapDerivs.Stop(); 
		CS.MmapDerivs += MmapDerivs.Elapsed(); 
    }
}
 

/* ==================  MPI Multipoles ================== */

/// Launch the CPD MPI_Irecv jobs.
/// Call this when the first slab is being finished.
/// We need to know the slabs that will be computed on each node.
void ParallelConvolution::RecvMultipoleSlab(int first_slab_finished) {
	int offset = first_slab_finished - first_slab_on_node;
	for (int r = 0; r < MPI_size_x; r++) {
		
	    for (int xraw = first_slabs_all[r] + offset; xraw < first_slabs_all[r] + offset + total_slabs_all[r]; xraw++) {
			int x = CP->WrapSlab(xraw);

		    STDLOG(4,"MPI_Irecv set for incoming Multipoles on slab {:d}\n", x);

			int tag = (MPI_rank_x+1) * 10000 + x + M_TAG; 
			
			// We're receiving the full range of z's in each transfer.
			// Key thing is to route this to the correct x location
			MPI_Irecv(MTdisk + x * node_slab_elem, node_slab_elem,
			    MPI_MTCOMPLEX, r, tag, comm_multipoles, &Mrecv_requests[x]);
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

	Msend_requests[slab] = new MPI_Request[MPI_size_x];
	Msend_active++;
	for (int r = 0; r < MPI_size_x; r++) {
		int tag = (r+1) * 10000 + slab + M_TAG; 
		MPI_Issend(mt + node_start[r], node_size[r], MPI_MTCOMPLEX, r, tag, comm_multipoles, &Msend_requests[slab][r]);		
	}
    MsendTimer[slab].Start();
	STDLOG(2,"Multipole slab {:d} has been queued for MPI_Issend, Msend_active = {:d}\n", slab, Msend_active);
}


int ParallelConvolution::CheckRecvMultipoleComplete(int slab) {
    if(Mrecv_requests[slab] == MPI_REQUEST_NULL)
        return 1;

    int received = 0;
    int err = MPI_Test(&Mrecv_requests[slab], &received, MPI_STATUS_IGNORE);

	if (!received) STDLOG(4, "Multipole slab {:d} not received yet...\n", slab);
    else assert(Mrecv_requests[slab] == MPI_REQUEST_NULL);
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
	
    for (int r = 0; r < MPI_size_x; r++) {
		if (Msend_requests[x][r]==MPI_REQUEST_NULL) continue;  // Already done

	    err = MPI_Test(Msend_requests[x] + r, &sent, MPI_STATUS_IGNORE);
	
		if (!sent) {
		    done=0;
			break;
		} else {
            assert(Msend_requests[x][r] == MPI_REQUEST_NULL);
        }
    }
    
	if (done) {
        delete[] Msend_requests[x];
 		Msend_requests[x] = NULL;   // This will allow us to check faster

		Msend_active--;
        MsendTimer[slab].Stop();
		STDLOG(1,"MPI_Issend complete for MultipoleSlab {:d} after {:6.3f} sec, Msend_active = {:d}\n", x, MsendTimer[slab].Elapsed(), Msend_active);
    }
	
	return done; 
}

int ParallelConvolution::CheckForMultipoleTransferComplete(int _slab) {
	int slab = _slab % cpd;

    if(MultipoleSlabAllMPIDone[slab])
        return 1;

    // Has it been long enough since the last MPI query?
    if(!mpi_limiter.Try())
        return 0;

	int done = CheckRecvMultipoleComplete(slab) &&
			   CheckSendMultipoleComplete(slab); 
	if (done){
        STDLOG(1,"Multipole MPI work complete for slab {:d}.\n", slab);
        MultipoleSlabAllMPIDone[slab] = 1;
    }
	
	return done;
}

/* ================== MPI Taylors ================== */

///for a given x slab = slab, figure out which node rank is responsible for it and should receive its Taylors. 
int ParallelConvolution::GetTaylorRecipient(int slab, int offset){	
	for (int r = 0; r < MPI_size_x; r++) {
		for (int _x = first_slabs_all[r] + offset; _x < first_slabs_all[r] + offset + total_slabs_all[r]; _x++) {
			int x = _x % cpd;
			if (x == slab){
				STDLOG(4, "For slab {:d}, identified node {:d} as Taylor target\n", slab, r);
				return r;
			}
		} 
	}
    QUIT("Did not find target node for slab {:d}!", slab);
    return -1;
}

void ParallelConvolution::SendTaylors(int offset) {
	QueueTaylors.Start();

	//for (int slab = 0; slab < cpd; slab ++){ 
	for (int _slab = -cpd/2; _slab < cpd/2+1; _slab ++){ 
		//Each node has Taylors for a limited range of z but for every x slab. 
		//figure out who the receipient should be based on x slab and send to them. 
		// Take from MTdisk. Set Tsend_requests[x] as request. 	
		int slab = CP->WrapSlab(_slab);
		STDLOG(4, "About to SendTaylor Slab {:d} with offset {:d}\n", slab, offset); 
		
		int r = GetTaylorRecipient(slab, offset); //x-slab slab is in node r's domain. Send to node r. 
	
		int tag = (MPI_rank_x+1) * 10000 + slab + T_TAG;
		
		MPI_Issend(MTdisk + slab * node_slab_elem, node_slab_elem, MPI_MTCOMPLEX, r, tag, comm_taylors, &Tsend_requests[slab]);
		STDLOG(3,"Taylor slab {:d} has been queued for MPI_Issend to rank {:d}, offset {:d}\n", slab, r, offset);
	}
	STDLOG(2, "MPI_Issend set for outgoing Taylors.\n");
	QueueTaylors.Stop(); CS.QueueTaylors = QueueTaylors.Elapsed(); 
}
	
/// Allocate space for this TaylorSlab, and trigger the MPI calls
/// to gather the information from other nodes.  To be used in FetchSlab.
void ParallelConvolution::RecvTaylorSlab(int slab) {
	QueueTaylors.Start();
	// We are receiving a specified x slab. For each node, receive its z packet for
	// this x. For each node, use Trecv_requests[x][rank] as request. 
	slab = CP->WrapSlab(slab);

    MTCOMPLEX *mt = (MTCOMPLEX *) SB->GetSlabPtr(TaylorSlab, slab);
	
	Trecv_requests[slab] = new MPI_Request[MPI_size_x];
	Trecv_active++; 
	for (int r = 0; r < MPI_size_x; r++) {
		int tag = (r+1) * 10000 + slab + T_TAG;
		MPI_Irecv(mt + node_start[r], node_size[r], MPI_MTCOMPLEX, r, tag, comm_taylors, &Trecv_requests[slab][r]);
	}
	STDLOG(2,"MPI_Irecv setup for incoming Taylors on slab {:d}.\n", slab);
	
	QueueTaylors.Stop(); CS.QueueTaylors = QueueTaylors.Elapsed(); 
}

int ParallelConvolution::CheckTaylorRecvReady(int slab){
    if (Trecv_requests[slab]==NULL) return 1;  // Nothing to do for this slab
	
    int done=0;
    MPI_Testall(MPI_size_x, Trecv_requests[slab], &done, MPI_STATUSES_IGNORE);
	if (done) {
		delete[] Trecv_requests[slab];
		Trecv_requests[slab] = NULL; 
		Trecv_active--;
		STDLOG(2, "Taylor slab {:d} receives completed\n", slab);
	}
	
	return done;
}

int ParallelConvolution::CheckTaylorSendComplete(int slab){
    if(Tsend_requests[slab] == MPI_REQUEST_NULL)
        return 1;

	int sent = 0;
    int err = MPI_Test(&Tsend_requests[slab], &sent, MPI_STATUS_IGNORE);
	
    if (sent) {
    	assert(Tsend_requests[slab] == MPI_REQUEST_NULL);
    	STDLOG(1, "Taylor slab {:d} send completed\n", slab);
    }
	
	return sent; 
}

/// Check whether the MPI calls are complete for this TaylorSlab,
/// returning 1 if yes.  To be used in the TaylorForcePrecondition.
/// Also, have to be careful about the SendManifest.
int ParallelConvolution::CheckTaylorSlabReady(int slab) {
	slab = CP->WrapSlab(slab);

    if(TaylorSlabAllMPIDone[slab])
        return 1;

    // Has it been long enough since the last MPI query?
    if(!mpi_limiter.Try())
        return 0;

	int send_done = 0; int recv_done = 0; 
	send_done = CheckTaylorSendComplete(slab); 
	
	if (send_done) {
		recv_done = CheckTaylorRecvReady(slab);
	}
	
	if(send_done and recv_done){
        STDLOG(1, "Taylor MPI work complete for slab {:d}.\n", slab);
        TaylorSlabAllMPIDone[slab] = 1;
    }

	return (send_done and recv_done); 

}



/* ==================  Convolve functions ================== */

/// Create the MTzmxy buffer and copy into it from MTdisk.
// Note that this is also a type-casting from float to double.
void ParallelConvolution::Swizzle_to_zmxy() {
	STDLOG(2,"Entering Swizzle_to_zmxy\n");

	// in: [cpd][znode][rml][node_ky_size]
	// out: [znode][rml][cpd][node_ky_size]
	#pragma omp parallel for schedule(static) collapse(2)
	for(int64_t zoff = 0; zoff < znode; zoff++){
		for(int64_t m = 0; m < rml; m++){
			for(int64_t x = 0; x < cpd; x++){
				// The convolve code expects the data from file x (MTdisk[x]) to be in MTzmxy[x+1]
				int xout = x + 1;
				if (xout>=cpd) xout=0;
				for(int64_t ky = 0; ky < node_ky_size; ky++){
					MTzmxy[ zoff*rml*cpdky_pad + m*cpdky_pad + xout*node_ky_size + ky] =
						MTdisk[ x*znode*rml*node_ky_size + zoff*rml*node_ky_size + m*node_ky_size + ky ];
				}
			}
		}
	}
	
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
			for (int y = 0; y < node_ky_size; y++){	 
				// One factor of 1/CPD comes from the FFT below.
				// Two more are needed for the FFTs in slab taylor -- we choose to throw them in here. 
				MTdisk[ (int64)((xz  * rml + m)) * node_ky_size + y ] = 
					MTzmxy[ z * cpdky_pad * rml + m * cpdky_pad + x * node_ky_size + y ] * invcpd3;
			}
		}
	}

    fftw_free(MTzmxy);
	return;
}



fftw_plan ParallelConvolution::PlanFFT(int sign){
	fftw_plan plan = NULL;
	
	// We're going to plan to process each [x][y] slab separately,
	// doing the cpd FFTs of cpd length.
	// stride = cpd, dist = 1, nembed = NULL
	// We imported any wisdom in Prologue

    STDLOG(2, "MTzmxy byte alignment is {:d}\n", 1 << __builtin_ctzll((unsigned long long) MTzmxy));
	
	int  n[] = {(int) cpd};
    if(wisdom_exists){
        // can we enforce better alignment between rows...?
    	plan = fftw_plan_many_dft(1, n, node_ky_size, //we will do znode*rml of sets of node_ky_size FFTs.
    		(fftw_complex *) MTzmxy, NULL, node_ky_size, 1, //each new [x][y] chunk is located at strides of node_ky_size.
    		(fftw_complex *) MTzmxy, NULL, node_ky_size, 1,
    		sign, FFTW_PATIENT | FFTW_WISDOM_ONLY);
        if(plan == NULL){
            if(ReadState.FullStepNumber > 1)  // Haven't done any convolutions yet, don't expect wisdom!
                WARNING("Wisdom was imported but wisdom planning failed. Generating new plans...\n");
            else
                STDLOG(1,"Wisdom was imported but wisdom planning failed (probably okay because this is first convolution). Generating new plans...\n");
        }
        else{
            STDLOG(1,"Wisdom planning succeeded.\n");
        }
    }

    if(plan == NULL){
        plan = fftw_plan_many_dft(1, n, node_ky_size, //we will do znode*rml of sets of node_ky_size FFTs.
            (fftw_complex *) MTzmxy, NULL, node_ky_size, 1, //each new [x][y] chunk is located at strides of node_ky_size.
            (fftw_complex *) MTzmxy, NULL, node_ky_size, 1,
            sign, FFTW_PATIENT);
    }
    assertf(plan != NULL, "Failed to generate ParallelConvolve FFTW plan for sign {:d}\n", sign);
	
	return plan;
}


/// Create the plan and do the FFT
void ParallelConvolution::FFT(fftw_plan plan) {			
	// Now we loop over the combination of [znode][m]
	#pragma omp parallel for schedule(static) //pragma appears okay
	for (int zm = 0; zm < znode*rml; zm++) {
		fftw_execute_dft(plan, (fftw_complex *)(MTzmxy+zm*cpdky_pad),
				       (fftw_complex *)(MTzmxy+zm*cpdky_pad));
	}
	
	fftw_destroy_plan(plan);
	return;
}



/// This executes the convolve on the MTfile info.
void ParallelConvolution::Convolve() {
	STDLOG(1,"Beginning Convolution.\n");
	
    InCoreConvolution *ICC = new InCoreConvolution(order, cpd, node_ky_size, blocksize, MPI_size_z > 1, cpdky_pad);

	// We're beginning in [x][znode][m][y] order on the RAMdisk
	// Copy to the desired order in a new buffer
	FFTPlanning.Clear(); FFTPlanning.Start();
		fftw_plan forward_plan  = PlanFFT(FFTW_FORWARD);
		fftw_plan backward_plan = PlanFFT(FFTW_BACKWARD);
	FFTPlanning.Stop();
	CS.FFTPlanning = FFTPlanning.Elapsed();	
	
		
	ArraySwizzle.Clear(); ArraySwizzle.Start(); 	
		Swizzle_to_zmxy();   // The -1 x offset is applied here
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
			ICC->InCoreConvolve(MTzmxy+z*rml*cpdky_pad, Ddisk[z]);
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
	
#ifdef DO_NOTHING
	STDLOG(4,"DO_NOTHING test passed.\n");
#endif
	
	
	CS.ops = ICC->ConvolutionArithmeticCount();
	delete ICC;
	
	STDLOG(1,"Convolution complete.\n");
	
}



/* ================== REPORTING ================== */ 

void ParallelConvolution::dumpstats() {
    if (convtimebuffer==NULL) return;

    FILE *fp = fmemopen(convtimebuffer, CONVTIMEBUFSIZE, "w");

    //FILE *fp;
    //fp = fopen(fn,"w");
    //assert(fp!=NULL);


    double accountedtime  = CS.ConvolutionArithmetic;
           accountedtime += CS.ForwardZFFTMultipoles + CS.InverseZFFTTaylor;
           accountedtime += CS.ReadDerivatives;
		   accountedtime += CS.Constructor + CS.AllocMT + CS.AllocDerivs + CS.QueueTaylors + CS.FFTPlanning + CS.Destructor + CS.ThreadCleanUp; 
           accountedtime += CS.ArraySwizzle;
    double discrepancy = CS.ConvolveWallClock - accountedtime;
	
    int computecores = CS.ComputeCores;
    fprintf(fp,"Convolution parameters:  RamAllocated = %4.1fMiB CacheSize = %4.1fMiB nreal_cores=%d blocksize=%d zwidth=%d cpd=%d order=%d",
	         (double) CS.totalMemoryAllocated/(1<<20), P.ConvolutionCacheSizeMB, computecores, (int) blocksize, (int) znode, (int) cpd, (int) order);
			 

     fprintf(fp,"\n\n");
	 
     fprintf(fp,"ConvolutionWallClock:  %2.2e seconds \n", CS.ConvolveWallClock );
     fprintf(fp,"\t %50s : %1.2e seconds\n", "Constructor", CS.Constructor );
     fprintf(fp,"\t %50s : %1.2e seconds\n", "Allocating Multipole/Taylor", CS.AllocMT);
	 
     fprintf(fp,"\t %50s : %1.2e seconds\n", "Mmaping Multipole/Taylor (part of alloc)", CS.MmapMT);
     fprintf(fp,"\t %50s : %1.2e seconds\n", "Munmaping Multipole/Taylor (part of alloc)", CS.MunmapMT);
	 
     fprintf(fp,"\t %50s : %1.2e seconds\n", "Allocating Derivatives", CS.AllocDerivs);	 
	 
     fprintf(fp,"\t %50s : %1.2e seconds\n", "Mmaping Derivatives (part of alloc)", CS.MmapDerivs);
     fprintf(fp,"\t %50s : %1.2e seconds\n", "Munmaping Derivatives (part of alloc)", CS.MunmapDerivs);
	 
     fprintf(fp,"\t %50s : %1.2e seconds\n", "Array Swizzling", CS.ArraySwizzle );

     double e = CS.ReadDerivatives ? CS.ReadDerivativesBytes/CS.ReadDerivatives/(1.0e+6) : 0.;
     fprintf(fp,"\t %50s : %1.2e seconds --> rate was %4.0f MB/s\n", "ReadDiskDerivatives", CS.ReadDerivatives, e );


     double Gops = ((double) CS.ops)/(1.0e+9);
     fprintf(fp,"\t %50s : %1.2e seconds for %5.3f billion double precision operations per node\n", "Convolution Arithmetic", CS.ConvolutionArithmetic, Gops/MPI_size );
	 
     fprintf(fp,"\t %50s : %1.2e seconds\n", "FFT Planning", CS.FFTPlanning );
     fprintf(fp,"\t %50s : %1.2e seconds\n", "Forward FFT Z Multipoles", CS.ForwardZFFTMultipoles );
     fprintf(fp,"\t %50s : %1.2e seconds\n", "Inverse FFT Z Taylor",         CS.InverseZFFTTaylor );
	 
     fprintf(fp,"\t %50s : %1.2e seconds\n", "Queuing MPI send/recv for Taylors", CS.QueueTaylors );
	 
     fprintf(fp,"\t %50s : %1.2e seconds\n", "Destructor", CS.Destructor );
     fprintf(fp,"\t %50s : %1.2e seconds\n", "fftw thread clean up (part of destructor)", CS.ThreadCleanUp );
	 
	 

     fprintf(fp,"%50s : %1.2e seconds which is %d%% \n", "Unaccounted remaining wallclock time", discrepancy, (int) (discrepancy/CS.ConvolveWallClock*100) );

     double cae       = CS.ConvolutionArithmetic;
     double farithp   = cae/CS.ConvolveWallClock*100;
	 double setupp    = (CS.Constructor + CS.Destructor + CS.AllocMT + CS.AllocDerivs)/CS.ConvolveWallClock*100;
	 double planp     = (CS.FFTPlanning)/CS.ConvolveWallClock*100;
	 double sendp     = (CS.QueueTaylors)/CS.ConvolveWallClock*100;
	 double ffftp     = (CS.ForwardZFFTMultipoles + CS.InverseZFFTTaylor)/CS.ConvolveWallClock*100;
     double fiop      = (CS.ReadDerivatives)/CS.ConvolveWallClock*100;
	 double fclean    = (CS.ThreadCleanUp)/CS.ConvolveWallClock*100;

     double swzp      = CS.ArraySwizzle/CS.ConvolveWallClock*100;

     fprintf(fp,"\nSummary: Fourier Transforms = %2.0f%%     Convolution Arithmetic = %2.0f%%     Array Swizzle = %2.0f%%    Disk IO = %2.0f%%     ", ffftp, farithp, swzp, fiop );

     fprintf(fp,"\nSet up (con/destructor + allocs/mmaps) = %2.0f%%     FFT planning = %2.0f%%     FFT thread clean up = %2.0f%%   Queuing Taylor MPI send/recv = %2.0f%%    \n", setupp, planp, fclean, sendp);
     fprintf(fp,"          Arithmetic rate = %2.0f DGOPS --> rate per core = %1.1f DGOPS\n", cae ? Gops/cae : 0, cae ? Gops/cae/computecores/MPI_size : 0);
     fprintf(fp,"          [DGOPS == Double Precision Billion operations per second]\n");
     fprintf(fp,"\n");
     
    fclose(fp);
}
