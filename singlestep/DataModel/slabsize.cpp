/* slabsize.cpp
*
* \brief Provides a class to load and store the number of particles per slab from disk.
* 
* FIXME: DWF believes it is unclear why this is a seperate class that inherits from grid, rather
* than a property of the Grid itself. DWF reccomends refactoring so it is a method of Grid, with 
* the read/write components broken off into SlabSizeFile
*/

#include <stdio.h> 

class SlabSize: public grid {
    uint64 *_size;
    public:
    uint64 max;
    uint64 min;
    SlabSize(int cpd): grid(cpd) {
        _size = new uint64[cpd];
        max = 0;
        min = UINT64_MAX;
    }
    ~SlabSize(void) { delete[] _size; }

    // Provide access with a wrapped index.
    void set(int slab, uint64 size) { _size[WrapSlab(slab)] = size; }
    uint64 size(int slab) { return _size[WrapSlab(slab)]; }

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
        assertf(value==cpd, 
                "SlabSize file opens with incorrect CPD: %d\n", value);
        for (int j = 0; j<cpd; j++) {
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
        fprintf(fp, "%d\n", cpd);
        for (int j = 0; j<cpd; j++) {
            fprintf(fp, "%ld\n", _size[j]);
        }
        fclose(fp);
        return 0;
    }
};
