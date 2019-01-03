#define IOLOG(verbosity,...) { if (verbosity<=stdlog_threshold_global) { \
	LOG(iolog,__VA_ARGS__); iolog.flush(); } }

#define ioassertf(_mytest,...) do { \
    if (!(_mytest)) { \
            IOLOG(0,"Failed Assertion: %s\n", #_mytest); IOLOG(0,__VA_ARGS__); \
	    fprintf(stderr,"Failed Assertion: %s\n", #_mytest); \
	    fpprint(std::cerr, __VA_ARGS__); \
	    CrashIO(); \
    }} while(0)


#include <time.h>
#include "io_interface.h"
#include "io_request.h"
#include "ring.cpp"
#include "file.cpp"
#include "iolib.cpp"
#include <pthread.h>
#include <sched.h>
#include <sys/resource.h>
#include "threadaffinity.h"

class iothread {
public:
    iothread(char *logfn, int threadnum, int _io_core) {
        // Make the FIFO files
        // Note: we could move from ring buffers and FIFOs to a simple tbb:concurrent_bounded_queue,
        // as with the GPU module.  The current implementation is overkill since we don't need inter-process communication.
        pid_t pid = getpid();
        STDLOG(0,"Using pid %d in IO\n",pid);
        sprintf(IO_CMD_PIPE, "/tmp/iocmd_abacus.%d.%d", pid, threadnum);
        sprintf(IO_ACK_PIPE, "/tmp/ioack_abacus.%d.%d", pid, threadnum);
        remove_io_pipes();
        make_io_pipes();
        
        std::string _logfn = logfn;
        _logfn += "." + std::to_string(1ULL*threadnum);
        
        iolog.open(_logfn); 	// TODO: Probably need an error check

        int ramdisk = io_ramdisk_global;
        size_t diskbuffer = ((size_t) 128) << 10;  // 4 << 20 = 4 MB
        RD = new ReadDirect(ramdisk, diskbuffer);
        WD = new WriteDirect(ramdisk,diskbuffer);
        
        io_core = _io_core;

        /* 
        // Increase the CPU scheduling priority of the IO thread (probably has to be done with sudo)
        // Borrowed mostly from http://www.yonch.com/tech/82-linux-thread-priority
        // In practice, we haven't seen this speed up the IO, just make other parts of the code slower
        int res = 0;
        struct rlimit limit;
        res += getrlimit(RLIMIT_RTPRIO, &limit);
        assertf(res == 0, "Error %d getting resource limit\n", res);
        STDLOG(1,"Soft/hard priority limit: %d/%d\n", limit.rlim_cur, limit.rlim_max);

        pthread_attr_t tattr;
        int newprio = sched_get_priority_max(SCHED_FIFO);
        sched_param param;

        // initialized with default attributes
        res += pthread_attr_init(&tattr);

        // safe to get existing scheduling param
        res += pthread_attr_setinheritsched(&tattr, PTHREAD_EXPLICIT_SCHED);
        res += pthread_attr_setschedpolicy(&tattr, SCHED_FIFO);
        res += pthread_attr_getschedparam(&tattr, &param);

        // set the priority; others are unchanged
        STDLOG(1,"Changing IO thread priority from %d to %d\n", param.sched_priority, newprio);
        param.sched_priority = newprio;

        // setting the new scheduling param
        res += pthread_attr_setschedparam(&tattr, &param);*/

        // Launch io_thread() as a separate thread
        int res = 0;
        res += pthread_create(&io_pthread, NULL, iothread::start_thread, this);
        assertf(res == 0, "error %d starting io pthread!\n", res);
        STDLOG(0,"IO thread started!\n");

        // Open the pipes from the client side
        io_cmd = open(IO_CMD_PIPE, O_WRONLY);
        io_ack = open(IO_ACK_PIPE, O_RDONLY);
        STDLOG(0,"Done initializing IO\n");
    }
    
    static void *start_thread(void *iothread_obj){
        return ((iothread *) iothread_obj)->io_thread();
    }
    
    ~iothread() {
        STDLOG(0,"Terminating the IO thread\n");
        iorequest quitcmd; 
        quitcmd.command = IO_QUIT; 
        quitcmd.blocking = 1; 
        ssize_t ret = write(io_cmd, &quitcmd, sizeof(iorequest) ); 
        ioacknowledge quitcmdack; 
        ret = read(io_ack, &quitcmdack, sizeof(ioacknowledge) );  
        assertf(quitcmdack.command==IO_QUIT,
            "Error in IO acknowledgment\n"); 
        STDLOG(0,"Termination of IO thread confirmed\n");
        close(io_cmd);
        close(io_ack);
        remove_io_pipes();
        void * status =NULL;
        int rc = pthread_join(io_pthread, &status);
        delete RD;
        delete WD;
        assert(!rc);
        STDLOG(0,"Termination of IO complete\n");
        iolog.close();
    }
    
    void request(iorequest ior){
        ssize_t ret = write(io_cmd, &ior, sizeof(iorequest));
        if (ior.blocking) {
            wait_for_ioack(io_ack, ior.arenatype, ior.arenaslab);
            STDLOG(1,"Blocking IO returned\n");
        } else {
            STDLOG(2,"Non-blocking IO requested\n");
        }
        
    }
    
private:
    // FIFO file names
    char IO_CMD_PIPE[1024], IO_ACK_PIPE[1024];

    int fifo_cmd, fifo_ack;
    int io_cmd, io_ack;
    
    ReadDirect * RD;
    WriteDirect * WD;
    
    pthread_t io_pthread;
    int io_core = -1;
    
    std::ofstream iolog;
    
    void CrashIO() {
        IOLOG(0,"Crashing the IO thread; sending IO_ERROR to client!\n");
        ioacknowledge ioack(IO_ERROR,-1,-1);
        ssize_t ret = write(fifo_ack,&ioack, sizeof(ioacknowledge) );
        pthread_exit(NULL);
    }
    
    void wait_for_ioack(int io_ack, int arenatype, int arenaslab){
        ioacknowledge ioack;
        ssize_t ret = read(io_ack, &ioack, sizeof(ioacknowledge));
        assertf(ioack.arenatype == arenatype && ioack.arenaslab == arenaslab, "Error in IO acknowledgement for arena %d = %d, %d =%d\n", arenatype, ioack.arenatype, arenaslab, ioack.arenaslab);
    }

    // =====================================================================
    // The actual read/write commands

    void ReadIOR(iorequest *ior) {
        // Read the file, wait to complete.
        IOLOG(1,"Reading file %s\n", ior->filename);

        // Determine the ramdisk flag
        int ramdisk = -1;
        switch(ior->io_method){
            case IO_DIRECT:
                ramdisk = 0;
                break;
            case IO_FOPEN:
                ramdisk = 1;
                break;
            default:
                QUIT("Unknown IO method %d\n", ior->io_method);
        }

        RD->BlockingRead(ior->filename, ior->memory, ior->sizebytes, ior->fileoffset, ramdisk);

        IOLOG(1,"Done reading file\n");
        IO_SetIOCompleted(ior->arenatype, ior->arenaslab);
        return;
    }

    void WriteIOR(iorequest *ior) {
        IOLOG(1,"Writing file %s\n", ior->filename);
        // Write the file
        //ioassertf(FileExists(ior->filename)==0, 
        //	"File %s already exists; not the intended use of WriteFile.\n", ior->filename);
        //ioassertf(ior->fileoffset==0, 
        //	"WriteFile fileoffest = %d.  Non-zero values not supported.\n", ior->fileoffset);

        FILE * outfile = fopen(ior->filename,"wb");
        ioassertf(outfile != NULL,"Touching file %s failed\n", ior->filename);
        fclose(outfile);

        // Determine the ramdisk flag
        int ramdisk = -1;
        switch(ior->io_method){
            case IO_DIRECT:
                ramdisk = 0;
                break;
            case IO_FOPEN:
                ramdisk = 1;
                break;
            default:
                QUIT("Unknown IO method %d\n", ior->io_method);
        }

        WD->BlockingAppend(ior->filename, ior->memory, ior->sizebytes, ramdisk);

        IOLOG(1,"Done writing file\n");
        if (ior->deleteafterwriting==IO_DELETE) IO_DeleteArena(ior->arenatype, ior->arenaslab);
        return;
    }


    // ================================================================
    #define NINSTRUCTIONS 65536
    void *io_thread() {
        // This function runs as a stand-alone thread.
        // It must receive commands from the main program.
        // It inserts these commands into a set of ring buffers.
        // It then executes these commands in a priority order.

        if(io_core >= 0){
            set_core_affinity(io_core);
            STDLOG(0, "IO thread assigned to core %d\n", io_core);
        }
        else{
            STDLOG(0, "IO thread not bound to core\n");
        }

        /*{ // Related to the CPU scheduling experiment above
            // Double-check the thread priority
            int policy = 0, ret = 0;
            sched_param params;
            ret = pthread_getschedparam(io_pthread, &policy, &params);
            if (ret != 0){
               STDLOG(0, "Couldn't retrieve IO thread real-time scheduling params\n");
            } else {
                // Check the correct policy was applied
                if(policy != SCHED_FIFO) {
                    STDLOG(0,"IO thread scheduling is NOT SCHED_FIFO!\n");
                } else{
                    STDLOG(0, "IO thread schedule confirmed SCHED_FIFO\n");
                }
                // Print thread scheduling priority
                STDLOG(0,"IO thread scheduling priority is %d\n", params.sched_priority);
            }
        }*/

        IOLOG(0,"Opening IO pipes\n");
        fifo_cmd = open(IO_CMD_PIPE, O_RDONLY);
        fifo_ack = open(IO_ACK_PIPE, O_WRONLY);
        int highfd = fifo_cmd;
        fd_set set;

        // We maintain four instruction buffers, so as to sort by priority
        ringbuffer  read_blocking(NINSTRUCTIONS), read_nonblocking(NINSTRUCTIONS), 
                write_blocking(NINSTRUCTIONS), write_nonblocking(NINSTRUCTIONS);
        int quitflag = 0, wait_for_cmd = 0;
        IOLOG(0,"Starting IO loop\n");

        struct timeval timeout;
        timeout.tv_sec = 0;
        timeout.tv_usec = 1;

        while (1) {
            if (!quitflag) {
                // If quitflag!=0, skip this and quit if the buffers are empty.
                int fifo_not_empty = 1;
                do {
                    // Look at the fifo_cmd pipe to fill up the instruction buffers.
                    FD_ZERO(&set);
                    FD_SET(fifo_cmd,&set);
                    int ret = select(highfd+1,&set,NULL, NULL, &timeout);
                    assert(ret!=-1);
                    IOLOG(2,"Polling IO pipe: select() returned %d, wait_for_cmd %d\n", ret, wait_for_cmd);

                    // if wait_for_cmd==1 then we should wait for a cmd to appear.
                    // otherwise we would spin lock.  But if wait_for_cmd = 0, then 
                    // we are just trying to empty the cmd pipe without locking up.
                    if (wait_for_cmd == 1 || (ret>0 && FD_ISSET(fifo_cmd,&set))) {
                        IOLOG(2,"Reading from IO pipe\n");
                        // Read a command and place it in the buffers.
                        iorequest ior;
                        int nr = read(fifo_cmd, &ior, sizeof(iorequest));
                        assert(nr==sizeof(iorequest));
                        wait_for_cmd = 0;
                        // Put it in the buffer
                        if (ior.command==IO_READ) {
                        IOLOG(2,"Received IO read request: file = %s, arena type %d slab %d, blocking = %d\n", 
                            ior.filename, ior.arenatype, ior.arenaslab, ior.blocking);
                            if (ior.blocking==IO_BLOCKING) read_blocking.push(ior);
                            else read_nonblocking.push(ior);
                        } else if (ior.command==IO_WRITE) {
                        IOLOG(2,"Received IO write request: file = %s, arena type %d slab %d, blocking = %d\n", 
                            ior.filename, ior.arenatype, ior.arenaslab, ior.blocking);
                            if (ior.blocking==IO_BLOCKING) write_blocking.push(ior);
                            else write_nonblocking.push(ior);
                        } else if (ior.command==IO_QUIT) {
                        IOLOG(1,"Received QUIT command\n");
                            quitflag = 1; break;  // We're not taking any more commands!
                        } else assert(ior.command==IO_READ);  // Explode
                    } else fifo_not_empty = 0;
                } while(fifo_not_empty);
            }

            IOLOG(2,"Attempting to execute an IO command\n");

            // Do one instruction, chosen by priority
            if (write_blocking.isnotempty()) {
                IOLOG(2,"Starting blocking write\n"); 
                iorequest ior = write_blocking.pop(); 
                const char *dir = ior.dir;

                BlockingIOWriteTime[dir].Start();
                WriteIOR(&ior);
                BlockingIOWriteTime[dir].Stop();
                BlockingIOWriteBytes[dir] += ior.sizebytes;

                // Send an acknowledgement
                ioacknowledge ioack(IO_WRITE,ior.arenatype, ior.arenaslab);
                ssize_t ret = write(fifo_ack,&ioack, sizeof(ioacknowledge) );
                IOLOG(2,"IO_WRITE acknowledgement sent\n");
            } else if (read_blocking.isnotempty()) {
                IOLOG(2,"Starting blocking read\n"); 
                iorequest ior = read_blocking.pop();
                const char *dir = ior.dir;

                BlockingIOReadTime[dir].Start();
                ReadIOR(&ior);
                BlockingIOReadTime[dir].Stop();
                BlockingIOReadBytes[dir] += ior.sizebytes;

                // Send an acknowledgement
                ioacknowledge ioack(IO_READ,ior.arenatype, ior.arenaslab);
                ssize_t ret = write(fifo_ack,&ioack, sizeof(ioacknowledge) );
                IOLOG(2,"IO_READ acknowledgement sent\n");
            } else if (write_nonblocking.isnotempty()) {
                IOLOG(2,"Starting nonblocking write\n"); 
                iorequest ior = write_nonblocking.pop(); 
                const char *dir = ior.dir;

                NonBlockingIOWriteTime[dir].Start();
                WriteIOR(&ior);
                NonBlockingIOWriteTime[dir].Stop();
                NonBlockingIOWriteBytes[dir] += ior.sizebytes;
            } else if (read_nonblocking.isnotempty()) {
                IOLOG(2,"Starting nonblocking read\n");
                iorequest ior = read_nonblocking.pop();
                const char *dir = ior.dir;

                NonBlockingIOReadTime[dir].Start();
                ReadIOR(&ior);
                NonBlockingIOReadTime[dir].Stop();
                NonBlockingIOReadBytes[dir] += ior.sizebytes;
            } else if (quitflag) break;
            else wait_for_cmd = 1;	
            // We have no work to do, so we should read the cmd pipe to avoid spinlocking
        }

        // Terminate:  The IO is all done, so quit the pipes and remove the fifo files
        // Signal an acknowledgement.
        IOLOG(0,"Terminating the IO thread as planned.\n");
        ioacknowledge ioack(IO_QUIT,0,0);
        ssize_t ret = write(fifo_ack,&ioack, sizeof(ioacknowledge) );
        IOLOG(0,"IO_QUIT acknowledgement sent\n");
        close(fifo_cmd);
        close(fifo_ack);
        return NULL;		// Quit the thread!
    }

    void make_io_pipes(void) {
        STDLOG(0,"Making io pipes.\n");
        int ret_val;

        errno = 0;
        ret_val = mkfifo(IO_CMD_PIPE, 0666);
        assertf((ret_val!=-1)||(errno==EEXIST),
            "Error creating pipe IO_CMD_PIPE\n");

        errno = 0;
        ret_val = mkfifo(IO_ACK_PIPE, 0666);
        assertf((ret_val!=-1)||(errno==EEXIST),
            "Error creating pipe IO_ACK_PIPE\n");
    }

    void remove_io_pipes(void) {
        // Must remove old files
        STDLOG(0,"Deleting io pipe files\n");
        int ret = 0;
        if (FileExists(IO_ACK_PIPE)) {
            ret = remove(IO_ACK_PIPE);
        }
        assertf(ret == 0, "Error removing pipe IO_ACK_PIPE=\"%s\"\n", IO_ACK_PIPE);

        if (FileExists(IO_CMD_PIPE)) {
            ret = remove(IO_CMD_PIPE);
        }
        assertf(ret == 0, "Error removing pipe IO_CMD_PIPE=\"%s\"\n", IO_CMD_PIPE);
    }
};

// ================================================================ //

iothread **iothreads;
int niothreads;

void IO_Initialize(char *logfn) {
    // Count how many IO threads we need
    // IO thread numbers should be contiguous!
    niothreads = 1;
    for(int i = 0; i < P.nIODirs; i++)
        niothreads = max(niothreads, P.IODirThreads[i]);
    assertf(niothreads < MAX_IO_THREADS, "Too many io threads!\n");

    iothreads = new iothread*[niothreads];
    for(int i = 0; i < niothreads; i++){
        STDLOG(0,"Initializing IO thread %d\n", i + 1);
        iothreads[i] = new iothread(logfn, i + 1, P.IOCores[i]);
    }
}

void IO_Terminate() {
    for(int i = 0; i < niothreads; i++){
        STDLOG(0,"Terminating IO thread %d\n", i);
        delete iothreads[i];
    }

    delete[] iothreads;
}

// Return the ID of the IO thread that will handle this directory
// Threads are one-indexed
int GetIOThread(const char* dir){
    for(int i = 0; i < P.nIODirs; i++){
        if(strcmp(dir, P.IODirs[i]) == 0){
            return P.IODirThreads[i];
        }
    }
    
    return 1;
}


// Here are the actual interfaces for writing an arena
void ReadFile(char *ram, uint64 sizebytes, int arenatype, int arenaslab,
	    const char *filename, off_t fileoffset, int blocking) {
    
    STDLOG(2,"Using IO_thread module to read file %f, blocking %d\n", filename, blocking);
    iorequest ior(ram, sizebytes, filename, IO_READ, arenatype, arenaslab, fileoffset, 0, blocking);
    
    iothreads[GetIOThread(ior.dir) - 1]->request(ior);
}

void WriteFile(char *ram, uint64 sizebytes, int arenatype, int arenaslab, 
	    const char *filename, off_t fileoffset, int deleteafter, int blocking) {
    
    STDLOG(2,"Using IO_thread module to write file %f, blocking %d\n", filename, blocking);
    iorequest ior(ram, sizebytes, filename, IO_WRITE, arenatype, arenaslab, fileoffset, deleteafter, blocking );
    
    iothreads[GetIOThread(ior.dir) - 1]->request(ior);
}
