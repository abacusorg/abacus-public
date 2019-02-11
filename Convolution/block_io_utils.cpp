#ifdef DOUBLEPRECISION
#define MPI_MTCOMPLEX MPI_DOUBLE_COMPLEX
#else
#define MPI_MTCOMPLEX MPI_COMPLEX
#endif


//#define DEBUG

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
	
#ifdef PARALLEL
	int first_slab_on_disk; 
#endif

    int ramdisk_MT = 0;
    int ramdisk_derivs = 0;

    int derivs_file_size = 0;
    
public:
    MTCOMPLEX **mtblock = NULL; // multipoles/taylors
    DFLOAT **dblock = NULL;  // derivatives
    
    size_t alloc_bytes = 0;
    size_t ReadMultipoleBytes = 0, WriteTaylorBytes = 0, ReadDerivativeBytes = 0;
    PTimer ReadDerivatives, ReadMultipoles, TransposeMultipoles, WriteTaylor;
    
    Block(ConvolutionParameters &_CP) : CP(_CP),
                                        ReadDerivatives(_CP.niothreads),
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
#ifdef PARALLEL
		int last_slab = first_slab_on_node + total_slabs_on_node + 1;
		for (int _x = last_slab; _x < last_slab + cpd - total_slabs_on_node; _x++){
			int x = _x % cpd; 
			free(mtblock[x]);
		}
#endif
		
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
				

		for (int _x = first_slab_on_node; _x < first_slab_on_node + total_slabs_on_node; _x++) //each node must mmap full size of mtblock for all x. if _x in this loop is one of the nodes files, it will load it, otherwise it will initialize mmap to a bunch of zeros. 
		{			
			
	            int x = _x % cpd;

				STDLOG(1,"Reading multipoles for x-slab %d\n", x);

	            // Different threads are responsible for different files (but they all read into one block)
	            if (x % CP.niothreads != thread_num)
	                continue;

	            char fn[1024];

	            char tfn[1024];    // used by non-overwriting ramdisk
	            CP.MultipoleFN(x, fn);
	            CP.TaylorFN(x, tfn);


	            // The convolve code expects the data from file x to be in mtblock[x+1]
	            x = (x+1)%cpd;

	            if(ramdisk_MT){
	                int shm_fd_flags;
	                char *ramdisk_fn;

	                if(CP.OverwriteConvState){
	                    shm_fd_flags = O_RDWR;
	                    ramdisk_fn = fn;
	                } else {
	                    shm_fd_flags = O_RDWR | O_CREAT;
	                    ramdisk_fn = tfn;
	                }

	                // map the shared memory fd to an address
	                // If we're doing a ramdisk overwrite, this maps the multipoles directly into memory
	                // If it's not an overwrite, this maps our new Taylors write arena into memory (so we still have a read/memcpy to do)
	                // TODO: file_offset must be page-aligned. So we could back up to a page boundary and then return an offset pointer, not dissimilar to the DIO library
	                // Or is it that bad to require zwidth = full?
					
					
	                int fd = open(ramdisk_fn, shm_fd_flags, S_IRUSR | S_IWUSR);
	                assertf(fd != -1, "Failed to open shared memory file at \"%s\"\n", ramdisk_fn);
					

	                if(!CP.OverwriteConvState){
	                    // expand the Taylors file
	                    // TODO: page alignment?
	                    int res = ftruncate(fd, file_offset + size);
	                    assertf(res == 0, "ftruncate on shared memory ramdisk_fn = %s to size = %d failed\n", ramdisk_fn, file_offset + size);
	                }
                	mtblock[x] = (MTCOMPLEX *) mmap(NULL, size, PROT_READ | PROT_WRITE, MAP_SHARED, fd, file_offset);
					
					
	                int res = close(fd);
	                assertf((void *) mtblock[x] != MAP_FAILED, "%d mmap shared memory from fd = %d of size = %d at offset = %d failed\n", MPI_rank, fd, size, file_offset);
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
				
				
				STDLOG(1,"Finished reading multipoles for x-slab %d\n", x);
			
			
        }
        
		
		
        ReadMultipoles.Stop(thread_num);
    }
	
#ifdef PARALLEL
	void transpose_z_to_x(int zstart, int zwidth, int thread_num, int z_slabs_per_node, MTCOMPLEX * sendbuf, MTCOMPLEX * recvbuf, int * first_slabs_all, int * total_slabs_all, uint64_t sendbufsize, uint64_t recvbufsize){
		uint64_t rml_times_cpd = rml * cpd; 
		
		 int size_skewer;
		 MPI_Type_size(MPI_skewer, &size_skewer);
		 assert(size_skewer == sizeof(MTCOMPLEX)*rml_times_cpd);
		 
		STDLOG(1,"Beginning first transpose for zstart with int64 %d\n", zstart);
		 
	 

		int sendcounts[MPI_size];
		int    sdispls[MPI_size];
		int recvcounts[MPI_size];
		int    rdispls[MPI_size];
		
		STDLOG(1,"Allocated send/recvcounts/displs. %d\n", zstart);
		
		
		for (int i = 0; i < MPI_size; i++)
		{
			
#ifdef DEBUG	
			sendcounts[i] = total_slabs_on_node * rml_times_cpd; 
			sdispls[i]    = i * total_slabs_on_node * rml_times_cpd; 
			recvcounts[i] = total_slabs_all[i] * rml_times_cpd; 
			rdispls[i]    = (first_slabs_all[i] - first_slabs_all[0]) * rml_times_cpd;
			if (rdispls[i] < 0) rdispls[i]  += cpd;
			
#else		
            // All of these are in units of rml_times_cpd.
            // We only are sending one z per pairwise interaction, so we just need to know the range of x's.
            // send: contains the indexing of the outbound information for each rank
			sendcounts[i] = total_slabs_on_node;
			sdispls[i]    = i * total_slabs_on_node;
            // recv: contains the indexing of the inbound information for this rank
			recvcounts[i] = total_slabs_all[i];
			rdispls[i]    = (first_slabs_all[i] - first_slabs_all[0]);
            // DJE TODO: Is this assuming that first_slabs_all[] is monotonically increasing?
#endif					
		}
		
		STDLOG(1,"Populated send/recvcounts/displs. %d\n", zstart);
		
        for (int zbig=0; zbig<z_slabs_per_node; zbig++) {
          // We're going to have to do z_slabs_per_node separate MPI_Alltoallv() calls.
          // Each one will send one z (and all x's) to each node.
          // The z's being sent in each call are spaced out, so that after all of the
          // MPI calls, each node will have a contiguous range of z's.

		//for(int z=0; z< z_slabs_per_node * MPI_size; z++)
          for (int zr=0; zr<MPI_size; zr++) {
            int z = zbig+zr*z_slabs_per_node;    // The z that is intended for rank zr
			for(int x=0; x<total_slabs_on_node;x++){	
                if (z < zwidth) memcpy(
						&(sendbuf[rml_times_cpd * zr*total_slabs_on_node + rml_times_cpd * x + 0*cpd + 0]),
                        &(mtblock[(x + first_slab_on_node + 1) % cpd][rml_times_cpd*z + 0*cpd + 0 ]),
                        sizeof(MTCOMPLEX)*rml_times_cpd);
                else memset(
						&(sendbuf[rml_times_cpd * zr*total_slabs_on_node + rml_times_cpd * x + 0*cpd + 0]),
                        0, sizeof(MTCOMPLEX)*rml_times_cpd);

                /* // Equivalent to this:
		  		for(int m=0;m<rml;m++){
					for(int y=0;y<cpd;y++){
						uint64_t i = rml_times_cpd * zr*total_slabs_on_node + rml_times_cpd * x + m*cpd + y;
						
						assertf(sizeof(MTCOMPLEX)*i<sendbufsize, "%ld, %d %d %d %d, %ld", i, z, x, m, y, sendbufsize);
						
						
						if (z < zwidth) sendbuf[i] = mtblock[(x + first_slab_on_node + 1) % cpd][rml_times_cpd*z + m*cpd + y ];
						else sendbuf[i] = 0.0;
					}
				}
                */
			}
		  } // End zr loop
		
		STDLOG(1,"Populated sendbuf. %d\n", zbig);
	
        // TODO: are these barriers needed?
		MPI_Barrier(MPI_COMM_WORLD);
		
#ifdef DEBUG
		MPI_Alltoallv(sendbuf, sendcounts, sdispls, MPI_MTCOMPLEX, recvbuf, recvcounts, rdispls, MPI_MTCOMPLEX, MPI_COMM_WORLD);
#else
		MPI_Alltoallv(sendbuf, sendcounts, sdispls, MPI_skewer, recvbuf, recvcounts, rdispls, MPI_skewer, MPI_COMM_WORLD);
#endif
		MPI_Barrier(MPI_COMM_WORLD);
		
		
		STDLOG(1,"Did Alltoallv. %d\n", zbig);
		
	
		uint64_t r = 0; 
        int z = zbig+z_slabs_per_node * MPI_rank;  // This is the z for this node
        if (z<zwidth) {     // Skip if it's out of bounds; we received filler data
            for (int i = 0; i < MPI_size; i ++){			
				for(int x=0; x<total_slabs_all[i]; x++){
                    memcpy( &(mtblock[(x + first_slabs_all[i] + 1) % cpd][rml_times_cpd*z + 0*cpd + 0]),
                            &(recvbuf[r]),
                            sizeof(MTCOMPLEX)*rml_times_cpd);
                    r+=rml_times_cpd;
                    /* // Equivalent to this:
					for(int m=0;m<rml;m++){
						for(int y=0;y<cpd;y++){
							mtblock[(x + first_slabs_all[i] + 1) % cpd][rml_times_cpd*z + m*cpd + y] = recvbuf[r++];
						 }
	 				}
                    */
				}
			}
		}
		
        }
	
		STDLOG(1,"Finishing first transpose for zstart %d\n", zstart);
	}
	
	
	
	void transpose_x_to_z(int zstart, int zwidth, int thread_num, int z_slabs_per_node,  MTCOMPLEX * sendbuf, MTCOMPLEX * recvbuf, int * first_slabs_all, int * total_slabs_all){
		
		
		STDLOG(1,"Beginning second transpose for zstart %d\n", zstart);
		
		uint64_t rml_times_cpd = rml * cpd; 
		
		
		 int size_skewer;
		 MPI_Type_size(MPI_skewer, &size_skewer);
		 assert(size_skewer == sizeof(MTCOMPLEX)*rml_times_cpd);
 

		int sendcounts[MPI_size];
		int    sdispls[MPI_size];
		int recvcounts[MPI_size];
		int    rdispls[MPI_size];

		for (int i = 0; i < MPI_size; i++)
		{
#ifdef DEBUG
			sendcounts[i] = total_slabs_all[i] * rml_times_cpd; 
			sdispls[i]    = (first_slabs_all[i] - first_slabs_all[0]) * rml_times_cpd; 			
			recvcounts[i] = total_slabs_on_node * rml_times_cpd; 
			rdispls[i]    = i * total_slabs_on_node * rml_times_cpd; 
#else
			
			sendcounts[i] = total_slabs_all[i];
			sdispls[i]    = (first_slabs_all[i] - first_slabs_all[0]);
			if (sdispls[i] < 0) sdispls[i]  += cpd;
			recvcounts[i] = total_slabs_on_node;
			rdispls[i]    = i * total_slabs_on_node;
#endif

		}

        for (int zbig=0; zbig<z_slabs_per_node; zbig++) {
          // We're going to have to do z_slabs_per_node separate MPI_Alltoallv() calls.
          // Each one will send one z (and all x's) to each node.
          // The z's being sent in each call are spaced out, so that after all of the
          // MPI calls, each node will have a contiguous range of z's.
		
            uint64_t r = 0; 
            int z = zbig+z_slabs_per_node * MPI_rank;  // This is the z for this node
            for (int i = 0; i < MPI_size; i++){			
				for(int x=0; x<total_slabs_all[i]; x++,r+=rml_times_cpd){
                    if (z < zwidth) 
                        memcpy( &(sendbuf[r]),
                            &(mtblock[(x + first_slabs_all[i] + 1) % cpd][rml_times_cpd*z + 0*cpd + 0]),
                            sizeof(MTCOMPLEX)*rml_times_cpd);
                     else memset( &(sendbuf[r]), 0, sizeof(MTCOMPLEX)*rml_times_cpd);

                    /* // Equivalent to this:
					for(int m=0;m<rml;m++){
						for(int y=0;y<cpd;y++){
                            // TODO: Would be better to get this conditional out of the inner loop
							if (z < zwidth) sendbuf[r] = mtblock[(x + first_slabs_all[i] + 1) % cpd][rml_times_cpd*z + m*cpd + y];
							else sendbuf[r] = 0.0; 
							r++;
						}
					}
                    */
				}
			}
		
		
		MPI_Barrier(MPI_COMM_WORLD);
#ifdef DEBUG
		MPI_Alltoallv(sendbuf, sendcounts, sdispls, MPI_MTCOMPLEX, recvbuf, recvcounts, rdispls, MPI_MTCOMPLEX, MPI_COMM_WORLD);
#else
		MPI_Alltoallv(sendbuf, sendcounts, sdispls, MPI_skewer, recvbuf, recvcounts, rdispls, MPI_skewer, MPI_COMM_WORLD);
#endif
		
		MPI_Barrier(MPI_COMM_WORLD);
		
		// for(int z=0; z<z_slabs_per_node * MPI_size; z++)
          for (int zr=0; zr<MPI_size; zr++) {
            int z = zbig+zr*z_slabs_per_node;    // The z that is intended for rank zr
            if (z>=zwidth) continue;    // Skip writing any information; we received filler
            assertf(z>=0 && z<(cpd+1)/2, "z is out of bounds\n", z);
			for(int x=0; x<total_slabs_on_node;x++){
                memcpy( &( mtblock[(x + first_slab_on_node + 1) % cpd][rml_times_cpd*z + 0*cpd + 0 ]),
                        &(recvbuf[rml_times_cpd*zr*total_slabs_on_node + rml_times_cpd*x + 0*cpd + 0]),
                        sizeof(MTCOMPLEX)*rml_times_cpd);
                /* // Equivalent to this:
		  		for(int m=0;m<rml;m++){
					for(int y=0;y<cpd;y++){
						uint64_t i = rml_times_cpd*z*total_slabs_on_node + rml_times_cpd*x + m*cpd + y;
						mtblock[(x + first_slab_on_node + 1) % cpd][rml_times_cpd*z + m*cpd + y ] = recvbuf[i++];
					}
				}
                */
			}
		  }			
		
        }  // zbig
		STDLOG(1,"Finishing second transpose for zstart %d\n", zstart);
		
	}
	
#endif

    void write(int zstart, int zwidth, int thread_num){
        // zstart is only used to determine if this is the last iteration
        // and thus memory is eligible to be freed
		
		
        WriteTaylor.Start(thread_num);
        
        size_t size = sizeof(MTCOMPLEX)*zwidth*cpd*rml;

        for(int _x = first_slab_on_node; _x < first_slab_on_node + total_slabs_on_node; _x++) {
            int x = _x % cpd;
			
			STDLOG(1,"Beginning write for x %d\n", x);
			
			
            if (x % CP.niothreads != thread_num)
                continue;

            char fn[1024];
            CP.TaylorFN(x, fn);

			STDLOG(1, "Writing out to Taylor x-slab %d \n", x);	
			
            // The convolve code expects the data from file x to be in mtblock[x+1]
            x = (x+1)%cpd;			
            
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
        
			STDLOG(1,"Finishing write for x %d\n", x);
			
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
        mtblock = new MTCOMPLEX*[cpd]();
		
#ifdef PARALLEL
        size_t size = sizeof(MTCOMPLEX)*CP.z_slabs_per_node*cpd*rml;

    	//for x slabs not in this nodes domain, allocate mtblock room of z_slabs_per_node to receive transposes. 
		int last_slab = first_slab_on_node + total_slabs_on_node + 1;
		for (int _x = last_slab; _x < last_slab + cpd - total_slabs_on_node; _x++){
			int x = _x % cpd; 
			mtblock[x] = (MTCOMPLEX*) malloc(size);
			assert(mtblock[x]!=NULL);
		}
#endif
		
        if(ramdisk_MT){
			return;
        }

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
