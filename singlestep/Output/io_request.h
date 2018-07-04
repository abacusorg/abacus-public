// This is the structure that we pass through the pipe to the IO thread.


#ifndef IOREQUEST
#define IOREQUEST

#define IO_READ  1
#define IO_WRITE 2
#define IO_QUIT  3
#define IO_ERROR 4		// To signal trouble!

#include "file.cpp"

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
    char    *memory = NULL;		// Where the data is in memory
    uint64     sizebytes = 0;		// How much data
    char    filename[1024] = "";	// File name
    char    dir[1024] = "";	// Directory name
    int     command = 0; 		// use IO_READ, IO_WRITE, IO_QUIT
    int     arena = 0;		// Which arena number this is
    off_t     fileoffset = 0; 	// only used for reading
    int     deleteafterwriting = 0; // use IO_DELETE, IO_KEEP
    int     blocking = 0;		// use IO_BLOCKING, IO_NONBLOCKING

    void dumpior() {
        printf("IOR memory = %p ", memory);
        printf("sizebytes = %lu ", sizebytes);
        printf("filename = %s ", filename);
        printf("command = %d ", command);
        printf("arena = %d ", arena);
        printf("fileoffset = %lu ", fileoffset );
        printf("deleteafterwriting = %d ", deleteafterwriting);
        printf("blocking = %d\n", blocking );
    }

    iorequest() {
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
        containing_dirname(filename, dir);
        
        command = _command;
        arena = _arena;
        fileoffset = _fileoffset;
        deleteafterwriting = _deleteafterwriting;
        blocking = _blocking;
    }

    ~iorequest(void) { }
  
};

#endif //  IOREQUEST
