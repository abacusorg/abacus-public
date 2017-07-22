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
    STimer ReadDerivatives, ReadMultipoles, WriteTaylor;
    
    Block(ConvolutionParameters &_CP){
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

    void read(int zstart, int zwidth){
        ReadMultipoles.Start();
        
        size_t size = sizeof(MTCOMPLEX)*zwidth*cpd*rml;
        size_t file_offset = zstart*cpd*rml*sizeof(MTCOMPLEX);
        // The IO module can only do DIO if the file and memory buffer have the same alignemnt
        // We give the offset in units of MTCOMPLEXes
        int buffer_start_offset = (int)((file_offset%4096)/sizeof(MTCOMPLEX));
        
        for(int x=0;x<cpd;x++) {
            mtblock[x] = raw_mtblock[x] + buffer_start_offset;
            
            char fn[1024];
            int mapM_slab = (x+cpd-1)%cpd;
            sprintf(fn, "%s/%s_%04d", CP.runtime_MultipoleDirectory, CP.runtime_MultipolePrefix, mapM_slab);
            
            RD_RDM->BlockingRead( fn, (char *) mtblock[x], size, file_offset);
            
            // If this is the last loop iteration, we're safe to delete the multipoles
            if(CP.delete_multipoles_after_read
                && zstart + zwidth >= (cpd+1)/2){
                if (remove(fn) != 0)
                    perror("Error deleting file");
            }
            
            ReadMultipoleBytes += size;
        }
        
        ReadMultipoles.Stop();
    }

    void write(int zwidth){
        WriteTaylor.Start();
        
        size_t size = sizeof(MTCOMPLEX)*zwidth*cpd*rml;
        
        for(int x=0;x<cpd;x++) {
            char fn[1024];
            int remap_slab = (x+cpd-1)%cpd;
            sprintf(fn,"%s/%s_%04d",  CP.runtime_TaylorDirectory, CP.runtime_TaylorPrefix, remap_slab);
            
            // the mtblock alignment from the read should be valid for the write
            WD_WDT->BlockingAppend( fn, (char *) mtblock[x], size);
            WriteTaylorBytes += size;
        }
        
        WriteTaylor.Stop();
    }

    void read_derivs(int zstart, int zwidth){
        ReadDerivatives.Start();
     
        size_t size = sizeof(DFLOAT)*rml*CP.CompressedMultipoleLengthXY;
        
        char *fnfmt;
        if(sizeof(DFLOAT) == sizeof(float))
            fnfmt = "%s/fourierspace_float32_%d_%d_%d_%d_%d";
        else
            fnfmt = "%s/fourierspace_%d_%d_%d_%d_%d";
        // note the derivatives are stored in z-slabs, not x-slabs
        for(int z=zstart; z < zstart+zwidth; z++) {
            char fn[1024];
            sprintf(fn, fnfmt,
                    CP.runtime_DerivativesDirectory, 
                    (int) cpd, CP.runtime_order, CP.runtime_NearFieldRadius, 
                    CP.runtime_DerivativeExpansionRadius, z);

            RD_RDD->BlockingRead( fn, (char *) dblock[z-zstart], size, 0);
            
            ReadDerivativeBytes += size;
        }
        
        ReadDerivatives.Stop();
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
