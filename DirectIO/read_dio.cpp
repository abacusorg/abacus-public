#include "read_dio.h"

#define is_aligned(POINTER, BYTE_COUNT) \
    (((uintptr_t)(const void *)(POINTER)) % (BYTE_COUNT) == 0)
int ReadDirect::rdopenflags(void) {
    if(!allow_directio) 
        return O_RDONLY;
    else
        return O_RDONLY|O_DIRECT|O_LARGEFILE;
}

int ReadDirect::rdopenfd(const fs::path &fn, int readflags) {
    errno = 0;
    int fd = open(fn.c_str(), readflags, 0666);
    FAILERRNO;
    assert(fd!=-1);
    return fd;
}

void ReadDirect::BlockingDirectReadAligned(const fs::path &fn, char *x, size_t length, off_t fileoffsetbytes) {
    size_t lengthblocks = length/4096;          assert( lengthblocks * 4096 == length );
    off_t offsetblocks = fileoffsetbytes/4096; assert( offsetblocks * 4096 == fileoffsetbytes );

    int fd = rdopenfd(fn, rdopenflags() );

    errno = 0;
    int rval = lseek(fd,fileoffsetbytes,SEEK_SET); // seek to the end of the file
    FAILERRNO;
    assert(rval!=-1);

    size_t bytesread = 0;
    if(is_aligned(x,4096)){
    	while(bytesread < length) bytesread += read(fd, &(x[bytesread]), length-bytesread);
    }
    else{
        fmt::print(stderr, "read_dio: File aligned, but memory not!\n"); assert(false);
        
    while(bytesread<length) {
        size_t bytes2read = min(alignedbytes, length-bytesread);
        errno = 0;
        size_t br = read(fd, alignedbuffer, bytes2read);
        FAILERRNO;
        assert(br == bytes2read);
        memcpy(&(x[bytesread]), alignedbuffer, bytes2read);
        bytesread += bytes2read;
    }
    }
    assert(bytesread==length);
    
    //fdatasync(fd);
    //posix_fadvise(fd, fileoffsetbytes, length, POSIX_FADV_DONTNEED);

    close(fd);
}

void ReadDirect::Blockingfread(const fs::path &fn, char *x, size_t length, off_t fileoffsetbytes) {
    if (fileoffsetbytes + (off_t) length > fs::file_size(fn)){
        fmt::print(stderr, "Trying to read offset {:d} of length {:d} from file {} of size {:d}\n",(intmax_t) fileoffsetbytes, (intmax_t) length, fn, (intmax_t) fs::file_size(fn));
        assert(0);
    }

    FILE *fp = fopen(fn.c_str(),"rb");
    assert(fp!=NULL);
    errno = 0;
    fseek(fp,fileoffsetbytes, SEEK_SET);
    FAILERRNO;
    errno = 0;
    size_t bytesread = fread( x ,  sizeof(char), length, fp);
    FAILERRNO;
    assert(bytesread == length);
    
    //int fd = fileno(fp);
    //fdatasync(fd);
    //posix_fadvise(fd, fileoffsetbytes, length, POSIX_FADV_DONTNEED);

    fclose(fp);
}

void ReadDirect::BlockingReadDirect(const fs::path &fn, char *x, size_t length, off_t fileoffsetbytes) {
    assert( fileoffsetbytes + (off_t) length <= fs::file_size(fn) ); // can't read off the end of the file

    size_t bytesleft = length;
    size_t memoryoffsetbytes = 0;

    size_t fob = fileoffsetbytes; // fob == local fileoffsetbytes
    size_t nonalignedfileoffset =  (4096 - fob%4096)%4096;
    if(nonalignedfileoffset || length < 4096) {
        size_t bytes2read = min(nonalignedfileoffset, length);
        Blockingfread( fn, x, bytes2read, fileoffsetbytes);
        bytesleft -= bytes2read;
        memoryoffsetbytes += bytes2read;
        fob += bytes2read;
    }
    
    if(bytesleft==0) return;
    
    // If we get here, we should have read to a 4096 boundary
    assert(fob%4096==0);
    
    if(bytesleft >= 4096)   {
        size_t alignedlength = bytesleft - bytesleft%4096;
        assert(alignedlength%4096==0);  
        if(alignedlength>0) {
            BlockingDirectReadAligned(fn, &(x[memoryoffsetbytes]), alignedlength, fob);
            memoryoffsetbytes += alignedlength;
            fob += alignedlength;
            bytesleft -= alignedlength;
        }
    }
    
    if(bytesleft==0) return;
    Blockingfread( fn, &(x[memoryoffsetbytes]), bytesleft, fob );
}

void ReadDirect::BlockingRead(const fs::path &fn, char *x, size_t length, off_t fileoffsetbytes) {
    ReadDirect::BlockingRead(fn, x, length, fileoffsetbytes, !allow_directio);
}

void ReadDirect::BlockingRead(const fs::path &fn, char *x, size_t length, off_t fileoffsetbytes, int no_dio) {
    if(no_dio || !allow_directio) 
        Blockingfread(fn,x,length,fileoffsetbytes);
    else
        BlockingReadDirect(fn,x,length,fileoffsetbytes);
}
#undef is_aligned