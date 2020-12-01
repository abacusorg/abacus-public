#include "cell_header.h"
#include "packN_storage.cpp"
#include "threevector.hh"
#include "particle_subsample.cpp"

#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <fcntl.h>
#include <string.h>
#define __STDC_FORMAT_MACROS
#include <inttypes.h>
#include <stdint.h>
#include <errno.h>
#include "iolib.cpp"


// Reads from byte `read_start` to the end of the file.
// Allocates a buffer and loads the data into it.  Returns the buffer.
// Also sets 'datasize' to the size of the buffer.
// 'ramdisk' is a flag that indicates whether 'fp' is on a ramdisk or not.
void* loadData(char *fn, size_t read_start, int ramdisk, size_t* datasize){
    struct stat filestatus;
    stat(fn,&filestatus);
    size_t filesize = filestatus.st_size;
    *datasize = filesize - read_start;
    
    void *data = NULL;
    int ret = posix_memalign(&data, 4096, *datasize);
    assert(ret == 0 && data != NULL && "Malloc failed");
    
    ReadDirect rd(ramdisk, 1024*1024);

    rd.BlockingRead(fn, (char *) data, *datasize, read_start);
    
    return data;
}

// Read from `fn`, starting at `offset`
// The particles are written into `out`
// number of read particles is returned
template <int N, typename T>
uint64_t read_packN(char *fn, size_t offset, int zspace, double subsample_frac, ThreeVector<T> *posout, ThreeVector<T> *velout, uint64 *pidout){
    if(subsample_frac <= 0)
        return 0;

    size_t datasize;
    // Is reading the whole file then parsing it more efficient than doing fread from disk?
    // This way certainly uses more memory.
    packN<N> *data = (packN<N> *) loadData(fn, offset, 1, &datasize);
    if(datasize%N != 0){ //ensure the file is sensible
        fprintf(stderr, "Datasize %zd of file \"%s\" not divisible by %d.  Is this a pack%d file?\n", datasize, fn, N, N);
        free(data);
        exit(1);
    }
    uint64_t max_NP = datasize/N;  // upper limit on the number of particles
    
    if (datasize == 0){
        printf("Empty file encountered: %s\n", fn);
        free(data);
        return 0;
    }

    int do_subsample = subsample_frac > 0 && subsample_frac < 1;  // do we need to bother hashing?
        
    cell_header current_cell;
    current_cell.vscale = 0;   // Make this illegal

    int return_pos = posout != NULL;
    int return_vel = velout != NULL;
    int return_pid = pidout != NULL;

    uint64_t i = 0;
    for(uint64 j = 0; j < max_NP; j++){
        packN<N> p = data[j];

        if (p.iscell()) {
            current_cell = p.unpack_cell();
            if(!current_cell.islegal()){
                fprintf(stderr, "Illegal pack14 cell encountered.\n");
                free(data);
                exit(1);
            }
        } else {
            ThreeVector<T> pos;
            ThreeVector<T> vel;
            uint64_t id;
            p.unpack(pos, vel, id, current_cell);

            if(do_subsample && is_subsample_particle(id, subsample_frac, 0.) != 0)
                continue;

            if(return_pos){
                posout[i] = pos;
                if(zspace)
                    posout[i][0] += vel[0];
            }

            if(return_vel){
                velout[i] = vel;
            }

            if(return_pid){
                pidout[i] = id;
            }
            i++;
        }
    }
    
    free(data);
    
    return i;
}

// Traditional C++ "using" aliases didn't seem to work in the extern "C" block
#define ALIAS(NAME,N,T) \
uint64_t read_##NAME(char *fn, size_t offset, int zspace, double subsample_frac, ThreeVector<T> *posout, ThreeVector<T> *velout, uint64 *pidout){ \
    return read_packN<N,T>(fn, offset, zspace, subsample_frac, posout, velout, pidout); \
}

extern "C" {
    ALIAS(pack14, 14,double)
    ALIAS(pack14f,14,float)
    ALIAS(pack9, 9,double)
    ALIAS(pack9f,9,float)
}
