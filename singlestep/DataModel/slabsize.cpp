/* slabsize.cpp
*
* \brief Provides a class to load and store the number of particles per slab from disk.
*
* The format of the file on disk is:
* ```
* cpd,num_zsplit
* primary_(x=0,z=0),with_ghost_(x=0,z=0)
* primary_(x=0,z=1),with_ghost_(x=0,z=1)
* ...
* primary_(x=cpd-1,z=num_zsplit-1),with_ghost_(x=cpd-1,z=num_zsplit-1)
* ```
*/

#include <stdio.h> 

// The global instantiation of this object is called `SS`.
class SlabSize {
    uint64 *sizes;  // the slab sizess in the read state
    uint64 *newsizes;  // the slab sizess in the write state

    uint64 *sizes_with_ghost;  // in the read state
    uint64 *newsizes_with_ghost;  // in the write state

    int cpd;
    int num_zsplit;
    int zsplit;
    
    public:
    
    SlabSize(int _cpd, int _num_zsplit, int _zsplit) {
        // Each node will hold the full [cpd*num_zsplit] SlabSize arrays.
        // But all the setters and getters will just touch the [slab,this_zrank] row.

        cpd = _cpd;
        num_zsplit = _num_zsplit;
        zsplit = _zsplit;

        sizes = new uint64[cpd*num_zsplit];
        newsizes = new uint64[cpd*num_zsplit];
        sizes_with_ghost = new uint64[cpd*num_zsplit];
        newsizes_with_ghost = new uint64[cpd*num_zsplit];
        
        for (int j=0;j<cpd*num_zsplit;j++){
            sizes[j] = 0;
            newsizes[j] = 0;
            sizes_with_ghost[j] = 0;
            newsizes_with_ghost[j] = 0;
        }
    }

    SlabSize(Parameters &P) 
        : SlabSize(P.cpd, P.NumZRanks, MPI_rank_z) {
        load_from_params(P);
    }

    ~SlabSize(void) { 
        delete[] sizes; 
        delete[] newsizes;
        delete[] sizes_with_ghost;
        delete[] newsizes_with_ghost;
    }

    // Provide access with a wrapped index.
    void set(int slab, uint64 _size, uint64 _size_with_ghost) {
        int ws = Grid->WrapSlab(slab);
        newsizes[ws*num_zsplit + zsplit] = _size;
        newsizes_with_ghost[ws*num_zsplit + zsplit] = _size_with_ghost;
    }
    void setold(int slab, uint64 _size, uint64 _size_with_ghost) {
        int ws = Grid->WrapSlab(slab);
        sizes[ws*num_zsplit + zsplit] = _size;
        sizes_with_ghost[ws*num_zsplit + zsplit] = _size_with_ghost;
    }
    
    uint64 size(int slab) { return sizes[Grid->WrapSlab(slab)*num_zsplit + zsplit]; }
    uint64 size_with_ghost(int slab) { return sizes_with_ghost[Grid->WrapSlab(slab)*num_zsplit + zsplit]; }
    
    // For parallel codes, we want to gather the newsizes information
    void parallel_gather() {
        #ifdef PARALLEL
            // Send all non-zero values of newsizes to node 0
            // Since newsizes is set only when a slab is finished,
            // we can economize our code and just add the vectors.
            MPI_REDUCE_TO_ZERO(newsizes, cpd*num_zsplit, MPI_UINT64_T, MPI_SUM);
            MPI_REDUCE_TO_ZERO(newsizes_with_ghost, cpd*num_zsplit, MPI_UINT64_T, MPI_SUM);
        #endif
        return;
    }

    // For echoing the ghost sizes to the WriteState. Should be called after parallel_gather().
    uint64 total_new_size_with_ghost(){
        uint64 sum =0;
        for(uint64 i = 0; i < cpd*num_zsplit; i++) sum += newsizes_with_ghost[i];
        return sum;
    }

    void load_from_params(Parameters &P){
        fs::path filename = P.ReadStateDirectory / "slabsize";
        STDLOG(2,"Reading SlabSize file from {}\n", filename);
        read(filename);
    }

    void store_from_params(Parameters &P){
        parallel_gather();
        if (MPI_rank==0) {
            fs::path filename = P.WriteStateDirectory / "slabsize";
            STDLOG(2,"Writing SlabSize file to {}\n", filename);
            write(filename);
        }
    }

    // Read and write from files.  Return 0 if ok.
    int read(char *fname) {
        FILE *fp; 
        fp = fopen(fname.c_str(), "r");
        assertf(fp!=NULL, "Couldn't open SlabSize read file {}\n", fname);
        int nread;
        uint64 cpdval,zval;
        nread = fscanf(fp, "%ld,%ld", &cpdval, &zval);
        assertf(nread==2,
                "Couldn't read first entry from SlabSize file {}\n", fname);
        assertf(cpdval==cpd, 
                "SlabSize file opens with incorrect CPD: {:d}\n", cpdval);
        assertf(zval==num_zsplit, 
                "SlabSize file opens with incorrect num_zsplit: {:d}\n", zval);
        for (int j = 0; j<cpd*num_zsplit; j++) {
            nread = fscanf(fp, "%ld,%ld", sizes+j, sizes_with_ghost+j);
            assertf(nread==2, "SlabSize fie ended prematurely\n");
        }
        fclose(fp);
        return 0;
    }

    int write(char *fname) {
        FILE *fp;
        fp = fopen(fname, "w");
        assertf(fp!=NULL, "Couldn't open SlabSize write file %s\n", fname);
        fprintf(fp, "%ld,%ld\n", cpd, num_zsplit);
        for (int j = 0; j < cpd*num_zsplit; j++) {
            fprintf(fp, "%ld,%ld\n", newsizes[j], newsizes_with_ghost[j]);
        }
        fclose(fp);
        return 0;
    }
};
