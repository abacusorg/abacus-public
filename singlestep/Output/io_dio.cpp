#include "io_interface.h"
#include "file.cpp"


#include "read_dio.h"
#include "read_dio.cpp"
#include "write_dio.h"
#include "write_dio.cpp"

// This is all blocking, regardless of instruction!

void IO_Initialize(char *logfn) { return; }    // Nothing to do
void IO_Terminate() { return; }

void ReadFile(char *ram, uint64 sizebytes, int arena,
	    const char *filename, off_t fileoffset, int blocking) {

    STDLOG(1,"Using IO_dio module to read file %f\n", filename);
    // Read the file, wait to complete.
    char fn[1024]; strcpy(fn, filename);
    int ramdisk = io_ramdisk_global;
    int diskbuffer = 1024*512;
    ReadDirect RD(ramdisk, diskbuffer);
    
    BlockingIOReadTime.Start();
    
    RD.BlockingRead( fn, ram, sizebytes, fileoffset);
    
    BlockingIOReadTime.Stop();
    BlockingIOReadBytes += sizebytes;

    STDLOG(1,"Done reading file\n");
    IO_SetIOCompleted(arena);
    return;
}

void WriteFile(char *ram, uint64 sizebytes, int arena, 
	    const char *filename, off_t fileoffset, int deleteafter, int blocking) {

    STDLOG(1,"Using IO_dio module to write file %f\n", filename);
    // Write the file
    char fn[1024]; strcpy(fn, filename);
    assertf(FileExists(fn)==0, "File %s already exists; not the intended use of WriteFile.\n", fn);
    assertf(fileoffset==0, "WriteFile fileoffest = %d.  Non-zero values not supported.\n", fileoffset);

    //char cmd[1024];
    //sprintf(cmd,"touch %s", fn);

    //int rv = system(cmd);
    //assertf(rv!=-1, "Touching file %s failed\n", fn);

    FILE * outfile = fopen(fn,"a");
    assertf(outfile != NULL,"Touching file %s failed\n", fn);
    fclose(outfile);
    int ramdisk = io_ramdisk_global;
    int diskbuffer = 1024*512;
    WriteDirect WD(ramdisk, diskbuffer);
    
    BlockingIOWriteTime.Start();
    
    WD.BlockingAppend( fn, ram, sizebytes);
    
    BlockingIOWriteTime.Stop();
    BlockingIOWriteBytes += sizebytes;

    STDLOG(1,"Done writing file\n");
    if (deleteafter==IO_DELETE) IO_DeleteArena(arena);
    return;
}


