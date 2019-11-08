#include "io_interface.h"
#include "file.cpp"

// This is simple fopen commands.
// TODO: Not sure that this will work for files above 2 GB.
// TODO: This claims to be appending, but it appears to be just writing.


void IO_Initialize(char *logfn) { return; }    // Nothing to do
void IO_Terminate() { return; }


class DirectIO  {
public:

    DirectIO(void)  { } 
    ~DirectIO(void) {  }

    void directreadfd(char *x,  uint64 sizebytes, off_t fileoffset, uint64 memoryoffset, const char *fn) { 
    
        FILE *fp;
        fp = fopen(fn,"rb");
        assertf(fp!=NULL,"Can't open DIO file %s\n", fn);
        errno = 0;
        fseek(fp,fileoffset, SEEK_SET);
        FAILERRNO;
        errno = 0;
        uint64 bytesread = fread( x ,  sizeof(char), sizebytes, fp);
        FAILERRNO;
    
        assertf(bytesread == sizebytes, 
		"DIO file %s ended early.  Read %d bytes, expected %d\n", 
		fn, bytesread, sizebytes);
        fclose(fp);

    }
 
    void directwritefd( char *x, uint64 sizebytes, off_t fileoffset, uint64 memoryoffset, const char *fn) { 
        // write appending 

        FILE *fp;
        fp = fopen(fn,"ab"); // appending
        assertf(fp!=NULL,"Can't open DIO write file %s\n", fn);
        errno = 0;
        uint64 byteswritten = fwrite( x,  sizeof(char), sizebytes, fp );
        FAILERRNO;
        assertf(byteswritten == sizebytes, 
		"DIO file %s ended early.  Wrote %d bytes, expected %d.\n", 
		fn, byteswritten, sizebytes);
        fclose(fp);

    
    }
};

DirectIO DIO;


// This is all blocking, regardless of instruction!

void ReadFile(char *ram, uint64 sizebytes, int arenatype, int arenaslab,
	    const char *filename, off_t fileoffset, int blocking) {
    STDLOG(1,"Using IO_fopen library to read file %f\n", filename);
    
    BlockingIOReadTime.Start();
    
    // Read the file, wait to complete.
    DIO.directreadfd( ram, sizebytes, fileoffset, 
	    0,      // memory offset
	    filename);
    
    BlockingIOReadTime.Stop();
    BlockingIOReadBytes += sizebytes;

    IO_SetIOCompleted(arenatype, arenaslab);
    
    return;
}

void WriteFile(char *ram, uint64 sizebytes, int arenatype, int arenaslab, 
	    const char *filename, off_t fileoffset, int deleteafter, int blocking) {
    STDLOG(1,"Using IO_fopen library to write file %f\n", filename);
    
    BlockingIOWriteTime.Start();

    // Write the file
    DIO.directwritefd( ram, sizebytes, fileoffset, 
    	0, // memory offset
	filename);
    
    BlockingIOWriteTime.Stop();
    BlockingIOWriteBytes += sizebytes;
    
    if (deleteafter==IO_DELETE) IO_DeleteArena(arenatype, arenaslab);
    return;
}


