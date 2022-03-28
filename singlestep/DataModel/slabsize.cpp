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
    uint64 *_size;  // the slab sizes in the read state
    uint64 *_newsize;  // the slab sizes in the write state

    uint64 *_size_with_ghost;  // in the read state
    uint64 *_newsize_with_ghost;  // in the write state

    int _cpd;
    int _num_zsplit;
    int _zsplit;
    
    public:
    
    SlabSize(int cpd, int num_zsplit, int zsplit) {
        // Each node will hold the full [cpd*num_zsplit] SlabSize arrays.
        // But all the setters and getters will just touch the [slab,this_zrank] row.

        _cpd = cpd;
        _num_zsplit = num_zsplit;
        _zsplit = zsplit;

        _size = new uint64[cpd];
        _newsize = new uint64[cpd*num_zsplit];
        _size_with_ghost = new uint64[cpd*num_zsplit];
        _newsize_with_ghost = new uint64[cpd*num_zsplit];
        
        for (int j=0;j<cpd*num_zsplit;j++){
            _size[j] = 0;
            _newsize[j] = 0;
            _size_with_ghost[j] = 0;
            _newsize_with_ghost[j] = 0;
        }
    }

    SlabSize(Parameters &P) 
        : SlabSize(P.cpd, P.NumZRanks, MPI_rank_z) {
        load_from_params(P);
    }

    ~SlabSize(void) { 
        delete[] _size; 
        delete[] _newsize;
        delete[] _size_with_ghost;
        delete[] _newsize_with_ghost;
    }

    // Provide access with a wrapped index.
    void set(int slab, uint64 size, uint64 size_with_ghost) {
        _newsize[Grid->WrapSlab(slab) + _zsplit] = size;
        _newsize_with_ghost[Grid->WrapSlab(slab) + _zsplit] = size_with_ghost;
    }
    void setold(int slab, uint64 size, uint64 size_with_ghost) {
        _size[Grid->WrapSlab(slab) + _zsplit] = size;
        _size_with_ghost[Grid->WrapSlab(slab) + _zsplit] = size_with_ghost;
    }
    uint64 size(int slab) { return _size[Grid->WrapSlab(slab) + _zsplit]; }
    uint64 size_with_ghost(int slab) { return _size_with_ghost[Grid->WrapSlab(slab) + _zsplit]; }

    // For parallel codes, we want to gather the newsize information
    void parallel_gather() {
        #ifdef PARALLEL
            // Send all non-zero values of _newsize to node 0
            // Since _newsize is set only when a slab is finished,
            // we can economize our code and just add the vectors.
            MPI_REDUCE_TO_ZERO(_newsize, _cpd*_num_zsplit, MPI_UINT64_T, MPI_SUM);
            MPI_REDUCE_TO_ZERO(_newsize_with_ghost, _cpd*_num_zsplit, MPI_UINT64_T, MPI_SUM);
        #endif
        return;
    }

    void load_from_params(Parameters &P){
        char filename[1024];
        int ret = snprintf(filename, 1024, "%s/slabsize", P.ReadStateDirectory);
        assert(ret >= 0 && ret < 1024);
        STDLOG(2,"Reading SlabSize file from %s\n", filename);
        read(filename);
    }

    void store_from_params(Parameters &P){
        parallel_gather();
        if (MPI_rank==0) {
            char filename[1024];
            int ret = snprintf(filename, 1024, "%s/slabsize",P.WriteStateDirectory);
            assert(ret >= 0 && ret < 1024);
            STDLOG(2,"Writing SlabSize file to %s\n", filename);
            write(filename);
        }
    }

    // Read and write from files.  Return 0 if ok.
    int read(char *fname) {
        FILE *fp; 
        fp = fopen(fname, "r");
        assertf(fp!=NULL, "Couldn't open SlabSize read file %s\n", fname);
        int nread;
        uint64 cpdval,zval;
        nread = fscanf(fp, "%ld,%ld", &cpdval, &zval);
        assertf(nread==2,
                "Couldn't read first entry from SlabSize file %s\n", fname);
        assertf(cpdval==_cpd, 
                "SlabSize file opens with incorrect CPD: %d\n", cpdval);
        assertf(zval==_num_zsplit, 
                "SlabSize file opens with incorrect num_zsplit: %d\n", zval);
        for (int j = 0; j<_cpd*_num_zsplit; j++) {
            nread = fscanf(fp, "%ld,%ld", _size+j, _size_with_ghost+j);
            assertf(nread==2, "SlabSize file ended prematurely\n");
        }
        fclose(fp);
        return 0;
    }

    int write(char *fname) {
        FILE *fp;
        fp = fopen(fname, "w");
        assertf(fp!=NULL, "Couldn't open SlabSize write file %s\n", fname);
        fprintf(fp, "%ld,%ld\n", _cpd, _num_zsplit);
        for (int j = 0; j < _cpd*_num_zsplit; j++) {
            fprintf(fp, "%ld,%ld\n", _newsize[j], _newsize_with_ghost[j]);
        }
        fclose(fp);
        return 0;
    }
};
