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

extern "C" {
#define RA_DOUBLEPRECISION
#include "readabacus_functions.cpp"

#undef RA_DOUBLEPRECISION
#include "readabacus_functions.cpp"

}