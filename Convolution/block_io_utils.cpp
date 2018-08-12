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
    
public:
    MTCOMPLEX **mtblock = NULL; // multipoles/taylors
    DFLOAT **dblock = NULL;  // derivatives
    
    size_t alloc_bytes = 0;
    size_t ReadMultipoleBytes = 0, WriteTaylorBytes = 0, ReadDerivativeBytes = 0;
    PTimer ReadDerivatives, ReadMultipoles, WriteTaylor;
    
    Block(ConvolutionParameters &_CP) : ReadDerivatives(_CP.niothreads),
                                        ReadMultipoles(_CP.niothreads),
                                        WriteTaylor(_CP.niothreads) {
        CP = _CP;
        cpd = CP.runtime_cpd; rml = CP.rml; alloc_zwidth = CP.zwidth;
        alloc();
        alloc_derivs();
    }
    
    ~Block(){
        for(int x = 0; x < cpd; x++)
            free(raw_mtblock[x]);
        for (int z = 0; z < alloc_zwidth; z++)
            free(dblock[z]);
        delete[] raw_mtblock;
        delete[] dblock;
    }

    void read(int zstart, int zwidth, int thread_num){
        ReadMultipoles.Start(thread_num);
        
        size_t size = sizeof(MTCOMPLEX)*zwidth*cpd*rml;
        size_t file_offset = zstart*cpd*rml*sizeof(MTCOMPLEX);
        // The IO module can only do DIO if the file and memory buffer have the same alignemnt
        // We give the offset in units of MTCOMPLEXes
        int buffer_start_offset = (int)((file_offset%4096)/sizeof(MTCOMPLEX));
        
        for(int x=0;x<cpd;x++) {
            // Different threads are responsible for different files (but they all read into one block)
            if (x % CP.niothreads != thread_num)
                continue;

            mtblock[x] = raw_mtblock[x] + buffer_start_offset;
            
            char fn[1024];
            int mapM_slab = (x+cpd-1)%cpd;
            CP.MultipoleDirectory(mapM_slab, fn);
            
            RD_RDM->BlockingRead( fn, (char *) mtblock[x], size, file_offset);
            
            // If this is the last loop iteration, we're safe to delete the multipoles
            if(CP.delete_multipoles_after_read
                && zstart + zwidth >= (cpd+1)/2){
                if (remove(fn) != 0)
                    perror("Error deleting file");
            }
            
            ReadMultipoleBytes += size;
        }
        
        ReadMultipoles.Stop(thread_num);
    }

    void write(int zwidth, int thread_num){
        WriteTaylor.Start(thread_num);
        
        size_t size = sizeof(MTCOMPLEX)*zwidth*cpd*rml;
        
        for(int x=0;x<cpd;x++) {
            if (x % CP.niothreads != thread_num)
                continue;

            char fn[1024];
            int remap_slab = (x+cpd-1)%cpd;
            CP.TaylorDirectory(remap_slab, fn);
            
            // the mtblock alignment from the read should be valid for the write
            WD_WDT->BlockingAppend( fn, (char *) mtblock[x], size);
            WriteTaylorBytes += size;
        }
        
        WriteTaylor.Stop(thread_num);
    }

    void read_derivs(int zstart, int zwidth, int thread_num){
        ReadDerivatives.Start(thread_num);
     
        size_t size = sizeof(DFLOAT)*rml*CP.CompressedMultipoleLengthXY;
        
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

            RD_RDD->BlockingRead( fn, (char *) dblock[z-zstart], size, 0);
            
            ReadDerivativeBytes += size;
        }
        
        ReadDerivatives.Stop(thread_num);
    }

private:
    // Allocate CPD separate chunks of space for multipoles.  This ensures they all have the same alignment offset
    void alloc(){
        raw_mtblock = new MTCOMPLEX*[cpd];
        mtblock = new MTCOMPLEX*[cpd];
        size_t s = sizeof(MTCOMPLEX) * alloc_zwidth * rml * cpd + 4096;  // wiggle room to adjust start to align with file
        for (int x = 0; x < cpd; x++){
            int memalign_ret = posix_memalign((void **) (raw_mtblock + x), 4096, s);
            assert(memalign_ret == 0);
            alloc_bytes += s;
        }
    }

    void alloc_derivs(){
        dblock = new DFLOAT*[alloc_zwidth];
        size_t s = sizeof(DFLOAT)*(rml*CP.CompressedMultipoleLengthXY);
        for (int z = 0; z < alloc_zwidth; z++){
            int memalign_ret = posix_memalign((void **) (dblock + z), 4096, s);
            assert(memalign_ret == 0);
            alloc_bytes += s;
        }
    }
};
