/* slabsize.cpp
*
* \brief Provides a class to load and store the number of particles per slab from disk.
*/

#include <stdio.h> 

// The global instantiation of this object is called `SS`.
class SlabSize {
    uint64 *_size;  // the slab sizes in the read state
    uint64 *_newsize;  // the slab sizes in the write state
    int _cpd;
    
    public:
    uint64 max;
    uint64 min;
    
    SlabSize(int cpd) {
        _cpd = cpd;
        _size = new uint64[cpd];
        _newsize = new uint64[cpd];
        STDLOG(2,"SlabSize vectors %p %p\n", (void *)_size, (void *)_newsize);
        for (int j=0;j<cpd;j++){
            _size[j] = 0;
            _newsize[j] = 0;
        }
        max = 0;
        min = UINT64_MAX;
    }

    SlabSize(Parameters &P) 
        : SlabSize(P.cpd) {
        load_from_params(P);
    }

    ~SlabSize(void) { 
        delete[] _size; 
        delete[] _newsize; 
    }

    // Provide access with a wrapped index.
    void set(int slab, uint64 size) { _newsize[Grid->WrapSlab(slab)] = size; }
    void setold(int slab, uint64 size) { _size[Grid->WrapSlab(slab)] = size; }
    uint64 size(int slab) { return _size[Grid->WrapSlab(slab)]; }

    // For parallel codes, we want to gather the newsize information
    void parallel_gather() {
        #ifdef PARALLEL
            // Send all non-zero values of _newsize to node 0
            // Since _newsize is set only when a slab is finished,
            // we can economize our code and just add the vectors.
            MPI_REDUCE_TO_ZERO(_newsize, _cpd, MPI_UINT64_T, MPI_SUM);
        #endif
        return;
    }

    void load_from_params(Parameters &P){
        char filename[1024];
        int ret = snprintf(filename, 1024, "%s/slabsize", P.ReadStateDirectory);
        assert(ret >= 0 && ret < 1024);
        read(filename);
        STDLOG(2,"Reading SlabSize file from %s\n", filename);
    }

    void store_from_params(Parameters &P){
        parallel_gather();
        if (MPI_rank==0) {
            char filename[1024];
            int ret = snprintf(filename, 1024, "%s/slabsize",P.WriteStateDirectory);
            assert(ret >= 0 && ret < 1024);
            write(filename);
            STDLOG(2,"Writing SlabSize file to %s\n", filename);
        }
    }

    // Read and write from files.  Return 0 if ok.
    int read(char *fname) {
        FILE *fp; 
        fp = fopen(fname, "r");
        assertf(fp!=NULL, "Couldn't open SlabSize read file %s\n", fname);
        int nread;
        uint64 value;
        nread = fscanf(fp, "%ld", &value); 
        assertf(nread==1,
                "Couldn't read first entry from SlabSize file %s\n", fname);
        assertf(value==Grid->cpd, 
                "SlabSize file opens with incorrect CPD: %d\n", value);
        for (int j = 0; j<Grid->cpd; j++) {
            int ret = fscanf(fp, "%ld", _size+j);
            if(_size[j] > max) max = _size[j];
            if(_size[j] < min) min = _size[j];
            assertf(nread==1, "SlabSize file ended prematurely\n");
        }
        fclose(fp);
        return 0;
    }

    int write(char *fname) {
        FILE *fp;
        fp = fopen(fname, "w");
        assertf(fp!=NULL, "Couldn't open SlabSize write file %s\n", fname);
        fprintf(fp, "%d\n", Grid->cpd);
        for (int j = 0; j<Grid->cpd; j++) {
            fprintf(fp, "%ld\n", _newsize[j]);
        }
        fclose(fp);
        return 0;
    }
};
