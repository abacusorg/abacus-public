#include "io_interface.h"
#include "file.cpp"


#include "read_dio.h"
#include "read_dio.cpp"
#include "write_dio.h"
#include "write_dio.cpp"

// This is all blocking, regardless of instruction!

void IO_Initialize(char *logfn) { return; }    // Nothing to do
void IO_Terminate() { return; }

void ReadFile(char *ram, uint64 sizebytes, int arenatype, int arenaslab,
	    const char *filename, off_t fileoffset, int blocking) {

    STDLOG(1,"Using IO_dio module to read file %f\n", filename);
    // Read the file, wait to complete.
    char fn[1024]; strcpy(fn, filename);
    int no_dio = !allow_directio_global;
    size_t diskbuffer = ((size_t) 4) << 20;  // 4 MB
    ReadDirect RD(no_dio, diskbuffer);
    
    char dir[1024];
    containing_dirname(filename, dir);
    
    BlockingIOReadTime[dir].Start();
    
    RD.BlockingRead( fn, ram, sizebytes, fileoffset);
    
    BlockingIOReadTime[dir].Stop();
    BlockingIOReadBytes[dir] += sizebytes;

    STDLOG(1,"Done reading file\n");
    IO_SetIOCompleted(arenatype, arenaslab);
    return;
}

void WriteFile(char *ram, uint64 sizebytes, int arenatype, int arenaslab, 
	    const char *filename, off_t fileoffset, int deleteafter, int blocking) {

    STDLOG(1,"Using IO_dio module to write file %f\n", filename);
    // Write the file
    char fn[1024]; strcpy(fn, filename);
    // When overwriting the state, we don't care if the file already exists
    // If we wanted to enforce that WriteFile can never overwrite files, we could delete the state file elsewhere
    //assertf(FileExists(fn)==0, "File %s already exists; not the intended use of WriteFile.\n", fn);
    assertf(fileoffset==0, "WriteFile fileoffest = %d.  Non-zero values not supported.\n", fileoffset);

    FILE * outfile = fopen(fn,"wb");
    assertf(outfile != NULL,"Touching file %s failed\n", fn);
    fclose(outfile);
    int no_dio = !allow_directio_global;
    size_t diskbuffer = ((size_t) 4) << 20;  // 4 MB
    WriteDirect WD(no_dio, diskbuffer);
    
    char dir[1024];
    containing_dirname(filename, dir);
    
    BlockingIOWriteTime[dir].Start();
    
    WD.BlockingAppend(fn, ram, sizebytes);
    
    BlockingIOWriteTime[dir].Stop();
    BlockingIOWriteBytes[dir] += sizebytes;

    STDLOG(1,"Done writing file\n");
    if (deleteafter==IO_DELETE) IO_DeleteArena(arenatype, arenaslab);
    return;
}


