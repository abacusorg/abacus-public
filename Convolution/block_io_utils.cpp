class Block {
    /* All of our IO is done on z-blocks.  The block class is a wrapper
       to manage memory and IO on these blocks.
    */
    
private:
    // this is what we alloc and free.
    // "mtblock" will be offset from this by up to 4096 to ensure DIO
    MTCOMPLEX **raw_mtblock = NULL;
    ConvolutionParameters CP;
    uint64_t cpd = 0, rml = 0, alloc_zwidth = 0;

    int ramdisk_MT = 0;
    int ramdisk_derivs = 0;

    int derivs_file_size = 0;
    
public:
    MTCOMPLEX **mtblock = NULL; // multipoles/taylors
    DFLOAT **dblock = NULL;  // derivatives
    
    size_t alloc_bytes = 0;
    size_t ReadMultipoleBytes = 0, WriteTaylorBytes = 0, ReadDerivativeBytes = 0;
    PTimer ReadDerivatives, ReadMultipoles, TransposeMultipoles, WriteTaylor;
    
    Block(ConvolutionParameters &_CP) : ReadDerivatives(_CP.niothreads),
                                        ReadMultipoles(_CP.niothreads),
                                        WriteTaylor(_CP.niothreads) {
        CP = _CP;
        cpd = CP.runtime_cpd; rml = CP.rml; alloc_zwidth = CP.zwidth;

        derivs_file_size = sizeof(DFLOAT)*rml*CP.CompressedMultipoleLengthXY;

        // For simplicity, multipoles and taylors must be either both or neither on ramdisk
        assert(is_path_on_ramdisk(CP.runtime_MultipoleDirectory) == is_path_on_ramdisk(CP.runtime_TaylorDirectory));

        if(is_path_on_ramdisk(CP.runtime_MultipoleDirectory))
            ramdisk_MT = 1;

        if(is_path_on_ramdisk(CP.runtime_DerivativesDirectory))
            ramdisk_derivs = 1;

        alloc();
        alloc_derivs();
    }
    
    ~Block(){
        if(!ramdisk_MT){
            for(int x = 0; x < cpd; x++)
                if(raw_mtblock[x] != NULL)
                    free(raw_mtblock[x]);
            delete[] raw_mtblock;
        }
        delete[] mtblock;

        for (int z = 0; z < alloc_zwidth; z++){
            if(ramdisk_derivs){
                if(dblock[z] != NULL){
                    int res = munmap(dblock[z], derivs_file_size);
                    assertf(res == 0, "Failed to munmap derivs\n");
                }
            } else{
                free(dblock[z]);
            }
        }
        delete[] dblock;
    }

    void read(int zstart, int zwidth, int thread_num){
        ReadMultipoles.Start(thread_num);

        size_t size = sizeof(MTCOMPLEX)*zwidth*cpd*rml;
        size_t file_offset = zstart*cpd*rml*sizeof(MTCOMPLEX);
        // The IO module can only do DIO if the file and memory buffer have the same alignment
        // We give the offset in units of MTCOMPLEXes
        int buffer_start_offset = (int)((file_offset%4096)/sizeof(MTCOMPLEX));
				
       // for(int x = first_slab_on_node; x < first_slab_on_node + total_slabs_on_node; x++) {
	
        for(int x = 0; x < P.cpd; x++) {
		
            // Different threads are responsible for different files (but they all read into one block)
            if (x % CP.niothreads != thread_num)
                continue;

            char fn[1024];
            int mapM_slab = (x+cpd-1)%cpd;
            CP.MultipoleFN(mapM_slab, fn);
			
            if(ramdisk_MT){
                char tfn[1024];    // used by non-overwriting ramdisk
                int shm_fd_flags;
                char *ramdisk_fn;

                if(CP.OverwriteConvState){
                    shm_fd_flags = O_RDWR;
                    ramdisk_fn = fn;
                } else {
                    shm_fd_flags = O_RDWR | O_CREAT;
                    CP.TaylorFN(mapM_slab, tfn);
                    ramdisk_fn = tfn;
                }

                int fd = open(ramdisk_fn, shm_fd_flags, S_IRUSR | S_IWUSR);
                assertf(fd != -1, "Failed to open shared memory file at \"%s\"\n", ramdisk_fn);

                if(!CP.OverwriteConvState){
                    // expand the Taylors file
                    // TODO: page alignment?
                    int res = ftruncate(fd, file_offset + size);
                    assertf(res == 0, "ftruncate on shared memory ramdisk_fn = %s to size = %d failed\n", ramdisk_fn, file_offset + size);
                }

                // map the shared memory fd to an address
                // If we're doing a ramdisk overwrite, this maps the multipoles directly into memory
                // If it's not an overwrite, this maps our new Taylors write arena into memory (so we still have a read/memcpy to do)
                // TODO: file_offset must be page-aligned. So we could back up to a page boundary and then return an offset pointer, not dissimilar to the DIO library
                // Or is it that bad to require zwidth = full?
                mtblock[x] = (MTCOMPLEX *) mmap(NULL, size, PROT_READ | PROT_WRITE, MAP_SHARED, fd, file_offset);
				
                int res = close(fd);
                assertf((void *) mtblock[x] != MAP_FAILED, "mmap shared memory from fd = %d of size = %d at offset = %d failed\n", fd, size, file_offset);
                assertf(mtblock[x] != NULL, "mmap shared memory from fd = %d of size = %d at offset = %d failed\n", fd, size, file_offset);
                assertf(res == 0, "Failed to close fd %d\n", fd);
            } else {
                mtblock[x] = raw_mtblock[x] + buffer_start_offset;
            }

            if(!(CP.OverwriteConvState && ramdisk_MT)){
                // Either we're not on the ramdisk, in which case we're doing an ordinary DIO read
                // or we are on the ramdisk but we aren't overwriting, so we need to read into the newly-allocated Taylors block
                RD_RDM->BlockingRead( fn, (char *) mtblock[x], size, file_offset, ramdisk_MT);
                ReadMultipoleBytes += size;
            }
        }
        
        ReadMultipoles.Stop(thread_num);
    }
	
#ifdef PARALLEL
	void transpose(int zstart, int zwidth, int thread_num){
		
		int rml_times_cpd = rml * cpd; 

		MTCOMPLEX * sendbuf;
		MTCOMPLEX * recvbuf;
		int sendcounts[MPI_size];
		int    sdispls[MPI_size];
		int recvcounts[MPI_size];
		int    rdispls[MPI_size];

		sendbuf = (MTCOMPLEX *) malloc(zwidth * total_slabs_on_node * rml_times_cpd * sizeof(MTCOMPLEX));
		recvbuf = (MTCOMPLEX *) malloc(rml_times_cpd * cpd * sizeof(MTCOMPLEX));


		for (int i = 0; i < MPI_size; i++)
		{
			sendcounts[i] = total_slabs_on_node * rml_times_cpd; // send total_slabs_on_node * rml * cpd complex numbers)
			sdispls[i]    = i * total_slabs_on_node * rml_times_cpd; // contiguous chunk starts after sdispls[i] complex numbers.
			recvcounts[i] = total_slabs_all[i] * rml_times_cpd; //*zwidth
			rdispls[i]    = first_slabs_all[i] * rml_times_cpd; //doesn't work for zwidth/node > 1
		}
		
		
		
		for(int z=zstart; z<zstart + zwidth; z++){
			for(int x=0; x<total_slabs_on_node;x++){
		  		for(int m=0;m<rml;m++){
					for(int y=0;y<cpd;y++){
						int i = z*total_slabs_on_node*rml_times_cpd + x*rml_times_cpd + m*cpd + y;
						sendbuf[i] = mtblock[x + first_slab_on_node][z*rml_times_cpd + m*cpd + y ];
					}
				}
			}
		}
		
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Alltoallv(sendbuf, sendcounts, sdispls, MPI_COMPLEX, recvbuf, recvcounts, rdispls, MPI_COMPLEX, MPI_COMM_WORLD);
		MPI_Barrier(MPI_COMM_WORLD);
		
		int z = zstart + MPI_rank; 

		for(int x=0; x<P.cpd;x++){
		  		for(int m=0;m<rml;m++){
					for(int y=0;y<cpd;y++){
						
						mtblock[x][z*rml_times_cpd + m*cpd + y] = recvbuf[x*rml_times_cpd + m*cpd + y]; 
						
						// int i = x*rml_times_cpd + m*cpd + y;
 						//printf("%d %d %d %d %f \n", MPI_rank, x, m, y, recvbuf[i]);
					}
				}
			}
			
		free(sendbuf);
		free(recvbuf);
	
	}
#endif

    void write(int zstart, int zwidth, int thread_num){
        // zstart is only used to determine if this is the last iteration
        // and thus memory is eligible to be freed
        WriteTaylor.Start(thread_num);
        
        size_t size = sizeof(MTCOMPLEX)*zwidth*cpd*rml;
        
        for(int x = 0; x < cpd; x++) {
            if (x % CP.niothreads != thread_num)
                continue;

            char fn[1024];
            int remap_slab = (x+cpd-1)%cpd;
            CP.TaylorFN(remap_slab, fn);
            
            if(ramdisk_MT){
                // If mtblock is on the ramdisk, then the Taylors swizzle was effectively our "write"
                // That's because in read() we mapped the multipoles file directly if overwriting,
                // or we created and mapped a new (section of a) Taylors file if not overwriting.
                // All that remains is to munmap.
                int res = munmap((void *) mtblock[x], size);
                assertf(res == 0, "Failed to munmap\n");
            } else {
                if(CP.OverwriteConvState){
                    // TODO: add DIO feature to overwrite files
                    assertf(!CP.OverwriteConvState, "DIO overwriting not yet implemented\n");
                } else{
                    // the mtblock alignment from the read should be valid for the write
                    WD_WDT->BlockingAppend(fn, (char *) mtblock[x], size, 0);
                    WriteTaylorBytes += size;
                }

                if(zstart + zwidth >= (cpd+1)/2){
                    free(raw_mtblock[x]);
                    raw_mtblock[x] = NULL;
                }
            }
        }
        
        WriteTaylor.Stop(thread_num);
    }

    void read_derivs(int zstart, int zwidth, int thread_num){
        ReadDerivatives.Start(thread_num);
     
        size_t size = derivs_file_size;
        
        const char *fnfmt;
        if(sizeof(DFLOAT) == sizeof(float))
            fnfmt = "%s/fourierspace_float32_%d_%d_%d_%d_%d";
        else
            fnfmt = "%s/fourierspace_%d_%d_%d_%d_%d";
        // note the derivatives are stored in z-slabs, not x-slabs
        for(int z=zstart; z < zstart+zwidth; z++) {
            // TODO: might be more efficient to chunk in blocks instead of stripes
            if (z % CP.niothreads != thread_num)
                continue;
            char fn[1024];
            sprintf(fn, fnfmt,
                    CP.runtime_DerivativesDirectory, 
                    (int) cpd, CP.runtime_order, CP.runtime_NearFieldRadius, 
                    CP.runtime_DerivativeExpansionRadius, z);

            if(!ramdisk_derivs){
                RD_RDD->BlockingRead( fn, (char *) dblock[z-zstart], size, 0);
                ReadDerivativeBytes += size;
            } else {
                if(dblock[z-zstart] != NULL){
                    int res = munmap(dblock[z-zstart], size);
                    assertf(res == 0, "Failed to munmap derivs\n");
                }
                int fd = open(fn, O_RDWR, S_IRUSR | S_IWUSR);
                assertf(fd != -1, "Failed to open shared memory file at \"%s\"\n", fn);
                
                // map the shared memory fd to an address
                dblock[z-zstart] = (DFLOAT *) mmap(NULL, size, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
                int res = close(fd);
                assertf((void *) dblock[z-zstart] != MAP_FAILED, "mmap shared memory from fd = %d of size = %d failed\n", fd, size);
                assertf(dblock[z-zstart] != NULL, "mmap shared memory from fd = %d of size = %d failed\n", fd, size);
                assertf(res == 0, "Failed to close fd %d\n", fd);
            }
        }
        
        ReadDerivatives.Stop(thread_num);
    }

private:
    // Allocate CPD separate chunks of space for multipoles.  This ensures they all have the same alignment offset
    void alloc(){
        mtblock = new MTCOMPLEX*[cpd];
        if(ramdisk_MT)
            return;

        raw_mtblock = new MTCOMPLEX*[cpd];
        size_t s = sizeof(MTCOMPLEX) * alloc_zwidth * rml * cpd + 4096;  // wiggle room to adjust start to align with file
        for (int x = 0; x < cpd; x++){
            int memalign_ret = posix_memalign((void **) (raw_mtblock + x), 4096, s);
            assert(memalign_ret == 0);
            alloc_bytes += s;
        }
    }

    void alloc_derivs(){
        dblock = new DFLOAT*[alloc_zwidth];
        for(int z = 0; z < alloc_zwidth; z++)
            dblock[z] = NULL;

        if(ramdisk_derivs)
            return;

        size_t s = sizeof(DFLOAT)*(rml*CP.CompressedMultipoleLengthXY);
        for (int z = 0; z < alloc_zwidth; z++){
            int memalign_ret = posix_memalign((void **) (dblock + z), 4096, s);
            assert(memalign_ret == 0);
            alloc_bytes += s;
        }
    }
};
