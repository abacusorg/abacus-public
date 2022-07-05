#include "write_dio.h"

// Some compilers didn't like when this was a macro
inline int is_aligned(const void* pointer, unsigned int byte_count){
    return (((uintptr_t)(const void *)(pointer)) % (byte_count) == 0);
}

int WriteDirect::wropenflags(void) {
    if(!allow_directio)
        return O_CREAT|O_APPEND|O_WRONLY;
    else
        return O_CREAT|O_WRONLY|O_DIRECT|O_ASYNC|O_LARGEFILE;
}

int WriteDirect::wropenfd(char *fn, int writeflags) {
    errno = 0;
    int fd = open(fn, writeflags, 0666);
    if(fd == -1)
        fprintf(stderr, "Failed to open \"%s\" with writeflags %d\n", fn, writeflags);
    FAILERRNO;
    assert(fd!=-1);
    return fd;
}

void WriteDirect::BlockingWrite_Append_Aligned(char *fn, char *x, size_t length) {
    size_t lengthblocks = length/4096;
    assert( lengthblocks * 4096 == length );

    size_t filesize = fsize(fn);
    size_t filesizeblocks = filesize/4096;
    assert( filesizeblocks * 4096 == filesize );

    int fd = wropenfd(fn, wropenflags() );

    errno = 0;
    int rval = lseek(fd,filesize,SEEK_SET); // seek to the end of the file
    FAILERRNO;
    assert(rval!=-1);

    size_t byteswritten = 0;
    if(is_aligned(x,4096)){
    	while(byteswritten<length) {
    		byteswritten += write(fd, &(x[byteswritten]), length-byteswritten);
    	}
    }
    else{
        fprintf(stderr, "write_dio: File aligned, but memory not!\n"); assert(false);

        /*
        // One should strive to always pass aligned buffers to Direct IO
        // but this provides a way around it by memcpy-ing into an aligned buffer
        while(byteswritten<length) {
            size_t bytes2write = min(alignedbytes, length-byteswritten);
            memcpy(alignedbuffer, &(x[byteswritten]), bytes2write);
            errno = 0;
            size_t bw = write(fd, alignedbuffer, bytes2write);
            FAILERRNO;
            assert(bw == bytes2write);
            byteswritten += bytes2write;
        }
        */
    }
    assert(byteswritten==length);

    //fdatasync(fd);
    //posix_fadvise(fd, filesize, length, POSIX_FADV_DONTNEED);

    close(fd);
}

void WriteDirect::BlockingAppendfwrite( char *fn, char *x, size_t length ) {
    //size_t filesize = fsize(fn);
    FILE *fp = fopen(fn,"ab"); // appending
    assert(fp!=NULL);
    errno = 0;
    size_t byteswritten = fwrite( x,  sizeof(char), length, fp );
    FAILERRNO;
    assert(byteswritten == length);

    //int fd = fileno(fp);
    //fdatasync(fd);
    //posix_fadvise(fd, filesize, length, POSIX_FADV_DONTNEED);

    fclose(fp);
}

void WriteDirect::BlockingAppendPointer( FILE *f, char *x, size_t length) {
    // Check that file is not null
    assert(f != nullptr);

    // Write to file keeping track of bytes
    size_t byteswritten = fwrite(x, sizeof(char), length, f);

    // Ensure that the correct number of bytes were written
    FAILERRNO;
    assert(byteswritten == length);
}

void WriteDirect::BlockingAppendDirect(char *fn, char *x, size_t length) {
    size_t bytesleft = length;
    size_t offset = 0;

    // How many bytes do we need to write to bring the file to a 4096 boundary?
    size_t fs = fsize(fn);
    size_t nonaligneddiskbytes = (4096-(fs%4096))%4096;

    if(nonaligneddiskbytes) {
        size_t bytes2write = min( nonaligneddiskbytes, length);
        BlockingAppendfwrite(fn, x, bytes2write);
        offset += bytes2write;
        bytesleft -= bytes2write;
    }

    if(bytesleft) {
        // if there's more to write, the disk file should already have been aligned
        assert( fsize(fn)%4096 == 0 );

        size_t remainingblocks = bytesleft/4096;
        size_t alignedbytes = 4096*remainingblocks;

        BlockingWrite_Append_Aligned(fn, &(x[offset]), alignedbytes);
        bytesleft -= alignedbytes;
        offset += alignedbytes;

        if(bytesleft) BlockingAppendfwrite(fn, &(x[offset]) , bytesleft );
    }
}

void WriteDirect::BlockingAppend(char *fn, char *x, size_t length) {
    WriteDirect::BlockingAppend(fn, x, length, !allow_directio);
}

void WriteDirect::BlockingAppend(char *fn, char *x, size_t length, int no_dio) {
    if(no_dio || !allow_directio)
        BlockingAppendfwrite(fn,x,length);
    else
        BlockingAppendDirect(fn,x,length);
}

void WriteDirect::BlockingAppend( FILE *f, char *x, size_t length) {
    WriteDirect::BlockingAppendPointer(f, x, length);
}
