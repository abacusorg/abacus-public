// This is the structure that we pass through the pipe to the IO thread.


#ifndef IOREQUEST
#define IOREQUEST

#define IO_READ  1
#define IO_WRITE 2
#define IO_QUIT  3
#define IO_ERROR 4		// To signal trouble!

#include <libgen.h>

class ioacknowledge {
  // A trivial class to respond
  public:
    int     command; 		// use IO_READ, IO_WRITE, IO_QUIT
    int     arena;		// Which arena number this is, just to confirm
    ioacknowledge() { command = arena = -1; }
    ioacknowledge(int _command, int _arena) {
        command = _command; arena = _arena;
    }
    ~ioacknowledge() { }
};


class iorequest {
  public:
    char    *memory;		// Where the data is in memory
    uint64     sizebytes;		// How much data
    char    filename[1024];	// File name
    char    dir[1024];	// Directory name
    int     command; 		// use IO_READ, IO_WRITE, IO_QUIT
    int     arena;		// Which arena number this is
    off_t     fileoffset; 	// only used for reading
    int     deleteafterwriting; // use IO_DELETE, IO_KEEP
    int     blocking;		// use IO_BLOCKING, IO_NONBLOCKING

    void dumpior() {
        printf("IOR memory = %p ", memory);
        printf("sizebytes = %lld ", sizebytes);
        printf("filename = %s ", filename);
        printf("command = %d ", command);
        printf("arena = %d ", arena);
        printf("fileoffset = %lld ", fileoffset );
        printf("deleteafterwriting = %d ", deleteafterwriting);
        printf("blocking = %d\n", blocking );
    }

    iorequest() {
        memory = NULL;
        filename[0] = '\0';
        dir[0] = '\0';
        sizebytes = command = arena = fileoffset = deleteafterwriting = blocking = 0;
    }

    iorequest(
        char    *_memory,
        uint64     _sizebytes,
        const char    *_filename,
        int     _command,
        int     _arena,
        off_t     _fileoffset,
        int     _deleteafterwriting,
        int     _blocking) {
        
        memory = _memory;
        sizebytes = _sizebytes;
        strncpy(filename, _filename, 1024);
        
        // Get the directory of the file for logging purposes
        // Believe it or not, dirname modifies its argument
        char buffer[1024];
        strncpy(buffer,filename,1024);
        char *_dir = basename(dirname(buffer));
        strncpy(dir, _dir, 1024);
        int len = strlen(dir);  assertf(len < 1023, "Directory \"%s\" name too long!", dir);
        dir[len] = '/'; dir[len+1] = '\0';
        
        command = _command;
        arena = _arena;
        fileoffset = _fileoffset;
        deleteafterwriting = _deleteafterwriting;
        blocking = _blocking;
    }

    ~iorequest(void) { }
  
};

#endif //  IOREQUEST
