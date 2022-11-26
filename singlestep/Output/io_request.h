// This is the structure that we pass through the pipe to the IO thread.


#ifndef IOREQUEST
#define IOREQUEST

enum iocommand { IO_READ = 1,
                    IO_WRITE,
                    IO_QUIT,
                    IO_ERROR
};

enum iomethod { IO_DIRECT,
                IO_FOPEN,
                IO_RAMDISK,
                IO_LIGHTCONE
};

#include "file.cpp"

class ioacknowledge {
  // A trivial class to respond
  public:
    int     command; 		// use IO_READ, IO_WRITE, IO_QUIT
    int     arenatype, arenaslab;		// Which arena number this is, just to confirm
    ioacknowledge() { command = arenatype = arenaslab = -1; }
    ioacknowledge(int _command, int _arenatype, int _arenaslab) {
        command = _command; arenatype = _arenatype, arenaslab = _arenaslab;
    }
    ~ioacknowledge() { }
};


class iorequest {
  public:
    char    *memory = NULL;		// Where the data is in memory
    uint64  sizebytes = 0;		// How much data
    char    filename[1024] = "";	// File name
    char    dir[1024] = "";	        // The name of the containing directory (i.e. not the full path)
    char    fulldir[1024] = "";
    char    justname[1024] = "";
    int     command = 0; 		// use IO_READ, IO_WRITE, IO_QUIT
    int     arenaslab = 0;		// Which slab number this is
    int     arenatype = 0;		// Which slab type this is
    off_t   fileoffset = 0; 	// only used for reading
    int     deleteafterwriting = 0; // use IO_DELETE, IO_KEEP
    int     blocking = 0;		// use IO_BLOCKING, IO_NONBLOCKING
    int     io_method = 0;  // use IO_DIRECT, IO_FOPEN
    int     do_checksum = 0;  // compute the crc
    FILE*   fp = nullptr; // File object for writing LightCones

    void dumpior() {
        printf("IOR memory = %p ", memory);
        printf("sizebytes = %lu ", sizebytes);
        printf("filename = %s ", filename);
        printf("command = %d ", command);
        printf("arenatype = %d ", arenatype);
        printf("arenaslab = %d ", arenaslab);
        printf("fileoffset = %lu ", fileoffset);
        printf("deleteafterwriting = %d ", deleteafterwriting);
        printf("blocking = %d\n", blocking);
    }

    iorequest() {
        memset(this, 0, sizeof(iorequest));   // Set to zero to appease valgrind
    }

    /* Construct an IO request
     * If fp is given, then that handle is used for IO with IO_LIGHTCONE
     * The filename should still be passed but will only be used for logging
     */
    iorequest(
        char    *_memory,
        uint64     _sizebytes,
        const char    *_filename,
        int     _command,
        int     _arenatype,
        int     _arenaslab,
        off_t     _fileoffset,
        int     _deleteafterwriting,
        int     _blocking,
        int     _do_checksum,
        FILE    *_fp = NULL) {

        memset(this, 0, sizeof(iorequest));   // Set to zero to appease valgrind

        memory = _memory;
        sizebytes = _sizebytes;
        strncpy(filename, _filename, 1024);

        if(_fp != NULL){
            io_method = IO_LIGHTCONE;
        } else {
            if(is_path_on_ramdisk(filename) || !allow_directio_global)
                io_method = IO_FOPEN;
            else
                io_method = IO_DIRECT;
        }

        // Get the directory of the file for logging purposes
        // Believe it or not, dirname modifies its argument
        containing_dirname(filename, dir);
        split_path(filename, fulldir, justname);

        command = _command;
        arenatype = _arenatype;
        arenaslab = _arenaslab;
        fileoffset = _fileoffset;
        deleteafterwriting = _deleteafterwriting;
        blocking = _blocking;
        do_checksum = _do_checksum;
        fp = _fp;
    }

    ~iorequest(void) { }

};

#endif //  IOREQUEST
