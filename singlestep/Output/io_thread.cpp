#include <time.h>
#include <pthread.h>
#include <sched.h>
#include <sys/resource.h>

#include "io_interface.h"
#include "io_request.h"
#include "ring.cpp"
#include "file.cpp"
#include "iolib.cpp"
#include "threadaffinity.h"

#define NINSTRUCTIONS 65536

#define IOLOG(verbosity,...) { if (verbosity<=stdlog_threshold_global) { \
    LOG(iolog,__VA_ARGS__); iolog.flush(); } }

#define ioassertf(_mytest,...) do { \
    if (!(_mytest)) { \
            IOLOG(0,"Failed Assertion: {:s}\n", #_mytest); IOLOG(0,__VA_ARGS__); \
        fmt::print(stderr,"Failed Assertion: {:s}\n", #_mytest); \
        fmt::print(stderr, __VA_ARGS__); \
        CrashIO(); \
    }} while(0)


class alignas(CACHE_LINE_SIZE) iothread {
public:
    iothread(fs::path logfn, int _threadnum, int _io_core)
        : io_core(_io_core),
          threadnum(_threadnum)
    {

        // Make the FIFO files
        // Note: we could move from ring buffers and FIFOs to a simple tbb:concurrent_bounded_queue,
        // as with the GPU module.  The current implementation is overkill since we don't need inter-process communication.
        pid_t pid = getpid();
        STDLOG(0,"Using pid {:d} in IO\n",pid);
        IO_CMD_PIPE = fmt::format("/tmp/iocmd_abacus.{:d}.{:d}", pid, threadnum);
        IO_ACK_PIPE = fmt::format("/tmp/ioack_abacus.{:d}.{:d}", pid, threadnum);
        remove_io_pipes();
        make_io_pipes();

        logfn += fmt::format(".{:d}", threadnum);

        iolog.open(logfn); 	// TODO: Probably need an error check

        int no_dio = !allow_directio_global;
        size_t diskbuffer = ((size_t) 128) << 10;  // 4 << 20 = 4 MB
        RD = new ReadDirect(no_dio, diskbuffer);
        WD = new WriteDirect(no_dio,diskbuffer);

        // Launch io_thread() as a separate thread
        int res = 0;
        res += pthread_create(&io_pthread, NULL, iothread::start_thread, this);
        assertf(res == 0, "error {:d} starting io pthread!\n", res);
        STDLOG(0,"IO thread started!\n");

        // Open the pipes from the client side
        io_cmd = open(IO_CMD_PIPE.c_str(), O_WRONLY);
        io_ack = open(IO_ACK_PIPE.c_str(), O_RDONLY);
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
        write(io_cmd, &quitcmd, sizeof(iorequest) );
        ioacknowledge quitcmdack;
        read(io_ack, &quitcmdack, sizeof(ioacknowledge) );
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

        // Report performance to the global structures
        ChecksumTime[threadnum] += checksum_timer.Elapsed();
        ChecksumBytes[threadnum] += checksum_bytes;
    }

    void request(iorequest ior){
        
        FifoWriteTimer.Start();
        write(io_cmd, &ior, sizeof(iorequest));
        FifoWriteTimer.Stop();
        if (ior.blocking) {
            wait_for_ioack(io_ack, ior.arenatype, ior.arenaslab);
            STDLOG(1,"Blocking IO returned\n");
        } else {
            STDLOG(3,"Non-blocking IO requested\n");
        }

    }

private:
    // FIFO file names
    fs::path IO_CMD_PIPE, IO_ACK_PIPE;

    int fifo_cmd, fifo_ack;
    int io_cmd, io_ack;

    ReadDirect * RD;
    WriteDirect * WD;

    pthread_t io_pthread;
    int io_core = -1;

    std::ofstream iolog;

    STimer checksum_timer;
    uint64 checksum_bytes = 0;

    int threadnum;

    void CrashIO() {
        IOLOG(0,"Crashing the IO thread; sending IO_ERROR to client!\n");
        ioacknowledge ioack(IO_ERROR,-1,-1);
        write(fifo_ack,&ioack, sizeof(ioacknowledge) );
        pthread_exit(NULL);
    }

    void wait_for_ioack(int io_ack, int arenatype, int arenaslab){
        ioacknowledge ioack;
        read(io_ack, &ioack, sizeof(ioacknowledge));
        assertf(ioack.arenatype == arenatype && ioack.arenaslab == arenaslab, "Error in IO acknowledgement for arena {:d} = {:d}, {:d} ={:d}\n", arenatype, ioack.arenatype, arenaslab, ioack.arenaslab);
    }

    // =====================================================================
    // The actual read/write commands

    void ReadIOR(iorequest *ior) {
        // Read the file, wait to complete.
        IOLOG(1,"Reading file {}\n", ior->filename);

        // Determine the ramdisk flag
        int no_dio = -1;
        switch(ior->io_method){
            case IO_DIRECT:
                no_dio = 0;
                break;
            case IO_FOPEN:
                no_dio = 1;
                break;
            default:
                QUIT("Unknown IO method {:d}\n", ior->io_method);
        }

        const char *dir = ior->dir;
        STimer *timer;
        if(ior->blocking){
            timer = &BlockingIOReadTime[dir];
            BlockingIOReadBytes[dir] += ior->sizebytes;
        }
        else{
            timer = &NonBlockingIOReadTime[dir];
            NonBlockingIOReadBytes[dir] += ior->sizebytes;
        }

        timer->Start();
        RD->BlockingRead(ior->filename, ior->memory, ior->sizebytes, ior->fileoffset, no_dio);
        timer->Stop();

        IOLOG(1,"Done reading file\n");
        IO_SetIOCompleted(ior->arenatype, ior->arenaslab);
        return;
    }

    void WriteIOR(iorequest *ior) {
        IOLOG(1,"Writing file {} using io_method {:d}\n", ior->filename, ior->io_method);
        // Write the file
        //ioassertf(FileExists(ior->filename)==0,
        //	"File {} already exists; not the intended use of WriteFile.\n", ior->filename);
        //ioassertf(ior->fileoffset==0,
        //	"WriteFile fileoffest = {:d}.  Non-zero values not supported.\n", ior->fileoffset);

        if (ior->io_method != IO_LIGHTCONE) {
            FILE * outfile = fopen(ior->filename,"wb");
            ioassertf(outfile != NULL,"Touching file {} failed\n", ior->filename);
            fclose(outfile);
        }

        // Determine the dio flag and if LightCone
        int no_dio = -1;
        switch(ior->io_method){
            case IO_DIRECT:
                no_dio = 0;
                break;
            case IO_FOPEN:
                no_dio = 1;
                break;
            case IO_LIGHTCONE:
                no_dio = 1;  // no-op, since we use the file pointer mechanism
                break;
            default:
                QUIT("Unknown IO method {:d}\n", ior->io_method);
        }

        if(ior->do_checksum){
            checksum_timer.Start();
            FileChecksums[ior->fulldir][ior->justname].ingest(ior->memory, ior->sizebytes);
            checksum_timer.Stop();
            checksum_bytes += ior->sizebytes;
        }

        const char *dir = ior->dir;
        STimer *timer;
        if(ior->blocking){
            timer = &BlockingIOWriteTime[dir];
            BlockingIOWriteBytes[dir] += ior->sizebytes;
        }
        else{
            timer = &NonBlockingIOWriteTime[dir];
            NonBlockingIOWriteBytes[dir] += ior->sizebytes;
        }

        // Use BlockingAppend for LightCones or non-LightCones
        if (ior->io_method==IO_LIGHTCONE)
        {
            timer->Start();
            WD->BlockingAppend(ior->fp, ior->memory, ior->sizebytes);
            timer->Stop();
        }
        else
        {
            timer->Start();
            WD->BlockingAppend(ior->filename, ior->memory, ior->sizebytes, no_dio);
            timer->Stop();
        }

        IOLOG(1,"Done writing file\n");
        if (ior->deleteafterwriting==IO_DELETE) IO_DeleteArena(ior->arenatype, ior->arenaslab);
        return;
    }


    // ================================================================
    void *io_thread() {
        // This function runs as a stand-alone thread.
        // It must receive commands from the main program.
        // It inserts these commands into a set of ring buffers.
        // It then executes these commands in a priority order.

        if(io_core >= 0){
            set_core_affinity(io_core);
            STDLOG(0, "IO thread assigned to core {:d}\n", io_core);
        }
        else{
            STDLOG(0, "IO thread not bound to core\n");
        }

        IOLOG(0,"Opening IO pipes\n");
        fifo_cmd = open(IO_CMD_PIPE.c_str(), O_RDONLY);
        fifo_ack = open(IO_ACK_PIPE.c_str(), O_WRONLY);
        int highfd = fifo_cmd;
        fd_set set;

        // We maintain four instruction buffers, so as to sort by priority
        // TODO: ringbuffers and IPC are total overkill unless we think we will ever
        // be siphoning off IO to another process (e.g. implementing our own ramdisk).
        // Otherwise, we could switch these to tbb::concurrent_queue, just like the GPU threads
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
                    IOLOG(3,"Polling IO pipe: select() returned {:d}, wait_for_cmd {:d}\n", ret, wait_for_cmd);

                    // if wait_for_cmd==1 then we should wait for a cmd to appear.
                    // otherwise we would spin lock.  But if wait_for_cmd = 0, then
                    // we are just trying to empty the cmd pipe without locking up.
                    if (wait_for_cmd == 1 || (ret>0 && FD_ISSET(fifo_cmd,&set))) {
                        IOLOG(3,"Reading from IO pipe\n");
                        // Read a command and place it in the buffers.
                        iorequest ior;
                        int nr = read(fifo_cmd, &ior, sizeof(iorequest));
                        assert(nr==sizeof(iorequest));
                        wait_for_cmd = 0;
                        // Put it in the buffer
                        if (ior.command==IO_READ) {
                        IOLOG(2,"Received IO read request: file = {}, arena type {:d} slab {:d}, blocking = {:d}\n",
                            ior.filename, ior.arenatype, ior.arenaslab, ior.blocking);
                            if (ior.blocking==IO_BLOCKING) read_blocking.push(ior);
                            else read_nonblocking.push(ior);
                        } else if (ior.command==IO_WRITE) {
                        IOLOG(2,"Received IO write request: file = {}, arena type {:d} slab {:d}, blocking = {:d}\n",
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

            IOLOG(3,"Attempting to execute an IO command\n");

            // Do one instruction, chosen by priority
            if (write_blocking.isnotempty()) {
                IOLOG(3,"Starting blocking write\n"); 
                iorequest ior = write_blocking.pop(); 
                
                WriteIOR(&ior);

                // Send an acknowledgement
                ioacknowledge ioack(IO_WRITE,ior.arenatype, ior.arenaslab);
                write(fifo_ack,&ioack, sizeof(ioacknowledge) );
                IOLOG(3,"IO_WRITE acknowledgement sent\n");
            } else if (read_blocking.isnotempty()) {

                IOLOG(3,"Starting blocking read\n"); 
                iorequest ior = read_blocking.pop();

                ReadIOR(&ior);

                // Send an acknowledgement
                ioacknowledge ioack(IO_READ,ior.arenatype, ior.arenaslab);
                write(fifo_ack,&ioack, sizeof(ioacknowledge) );
                IOLOG(3,"IO_READ acknowledgement sent\n");
            } else if (write_nonblocking.isnotempty()) {

                IOLOG(3,"Starting nonblocking write\n"); 
                iorequest ior = write_nonblocking.pop(); 

                WriteIOR(&ior);
            } else if (read_nonblocking.isnotempty()) {
                IOLOG(3,"Starting nonblocking read\n");
                iorequest ior = read_nonblocking.pop();

                ReadIOR(&ior);
            } else if (quitflag) break;
            else wait_for_cmd = 1;
            // We have no work to do, so we should read the cmd pipe to avoid spinlocking
        }

        // Terminate:  The IO is all done, so quit the pipes and remove the fifo files
        // Signal an acknowledgement.
        IOLOG(0,"Terminating the IO thread as planned.\n");
        ioacknowledge ioack(IO_QUIT,0,0);
        write(fifo_ack,&ioack, sizeof(ioacknowledge) );
        IOLOG(0,"IO_QUIT acknowledgement sent\n");
        close(fifo_cmd);
        close(fifo_ack);
        return NULL;		// Quit the thread!
    }

    void make_io_pipes(void) {
        STDLOG(0,"Making io pipes.\n");
        int ret_val;

        errno = 0;
        ret_val = mkfifo(IO_CMD_PIPE.c_str(), 0666);
        assertf((ret_val!=-1)||(errno==EEXIST),
            "Error creating pipe IO_CMD_PIPE\n");

        errno = 0;
        ret_val = mkfifo(IO_ACK_PIPE.c_str(), 0666);
        assertf((ret_val!=-1)||(errno==EEXIST),
            "Error creating pipe IO_ACK_PIPE\n");
    }

    void remove_io_pipes(void) {
        // Must remove old files
        STDLOG(0,"Deleting io pipe files\n");
        fs::remove(IO_ACK_PIPE);
        fs::remove(IO_CMD_PIPE);
    }
};

// ================================================================ //

#include <map>

iothread **iothreads;
size_t niothreads;

// Slab types, like light cones, may opt to append to a single file
// We keep an array of open file pointers for those types.
FILE **filepointers;
size_t nfilepointers;

void IO_Initialize(const fs::path &logfn, int NumTypes) {
    nfilepointers = NumTypes;
    filepointers = new FILE*[nfilepointers];
    for(size_t i = 0; i < nfilepointers; i++){
        filepointers[i] = NULL;
    }
    
    // Count how many IO threads we need
    // IO thread numbers should be contiguous!
    niothreads = 1;
    for(size_t i = 0; i < P.IODirThreads.size(); i++)
        niothreads = std::max(niothreads, (size_t) P.IODirThreads[i]);

    iothreads = new iothread*[niothreads];
    for(size_t i = 0; i < niothreads; i++){
        STDLOG(0,"Initializing IO thread {:d}\n", i + 1);
        int io_core = i < P.IOCores.size() ? P.IOCores[i] : -1;
        iothreads[i] = new iothread(logfn, i + 1, io_core);
    }
}

void IO_Terminate() {
    for(size_t i = 0; i < niothreads; i++){
        STDLOG(0,"Terminating IO thread {:d}\n", i);
        delete iothreads[i];
    }

    delete[] iothreads;

    for(size_t i = 0; i < nfilepointers; i++){
        if(filepointers[i] != NULL){
            int ret = fclose(filepointers[i]);
            assertf(ret == 0, "Error closing file pointer for type {:d}\n", i);
            filepointers[i] = NULL;
        }
    }
    delete[] filepointers;

    // Write the checksum files to their respective directories
    for(auto &diriter : FileChecksums){
        const fs::path dir = diriter.first;
        auto &_allcrc = diriter.second;

        // Sort this directory's checksums by filename
        std::map<std::string,CRC32> orderedcrc(_allcrc.begin(), _allcrc.end());

        // The format of this file should match that of GNU cksum so that one can do a simple diff
        // That's not to say that users should prefer cksum; it's slow!
        // We will instead provide a fast CRC32 utility.
        const fs::path checksumfn = dir / fmt::format("checksums{:s}.crc32", NodeString);
        FILE *fp = fopen(checksumfn.c_str(), "w");
        assertf(fp != NULL, "Failed to open file \"{}\"\n", checksumfn);

        for(auto &fileiter : orderedcrc){
            const std::string filename = fileiter.first;
            CRC32 crc = fileiter.second;
            fmt::print(fp, "{:d} {:d} {}\n", crc.finalize(), crc.size, filename);
        }
        fclose(fp);
    }
}

// Return the ID of the IO thread that will handle this directory
// Threads are one-indexed
int GetIOThread(const fs::path &dir){
    for(size_t i = 0; i < P.IODirThreads.size(); i++){
        if(dir == P.IODirs[i]){
            return P.IODirThreads[i];
        }
    }

    return 1;
}


// Here are the actual interfaces for writing an arena
void ReadFile(char *ram, uint64 sizebytes, int arenatype, int arenaslab, const fs::path &filename,
    off_t fileoffset, int blocking) {

    STDLOG(3,"Using IO_thread module to read file {}, blocking {:d}\n", filename, blocking);
    iorequest ior(ram, sizebytes, filename, IO_READ, arenatype, arenaslab, fileoffset, 0, blocking, 0);

    iothreads[GetIOThread(ior.dir) - 1]->request(ior);
}

void WriteFile(char *ram, uint64 sizebytes, int arenatype, int arenaslab, const fs::path &filename,
        off_t fileoffset, int deleteafter, int blocking, int do_checksum, int use_fp) {

    STDLOG(3,"Using IO_thread module to write file {}, blocking {:d}, use_fp {:d}\n", filename, blocking, use_fp);

    FILE *fp = NULL;
    if(use_fp){
        if(filepointers[arenatype] == NULL){
            filepointers[arenatype] = fopen(filename.c_str(), "wb");
            assertf(filepointers[arenatype] != NULL, "Failed to open file pointer for {}\n", filename);
        }

        fp = filepointers[arenatype];
    }

    iorequest ior(ram, sizebytes, filename, IO_WRITE, arenatype, arenaslab, fileoffset, deleteafter, blocking, do_checksum, fp);

    iothreads[GetIOThread(ior.dir) - 1]->request(ior);
}
