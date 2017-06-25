std::ofstream iolog;
#define IOLOG(verbosity,...) { if (verbosity<=stdlog_threshold_global) { \
	LOG(iolog,__VA_ARGS__); iolog.flush(); } }

#define ioassertf(_mytest,...) do { \
    if (!(_mytest)) { \
            IOLOG(0,"Failed Assertion: %s\n", #_mytest); IOLOG(0,__VA_ARGS__); \
	    fprintf(stderr,"Failed Assertion: %s\n", #_mytest); \
	    fpprint(std::cerr, __VA_ARGS__); \
	    CrashIO(); \
    }} while(0)


#include "io_interface.h"
#include "io_request.h"
#include "ring.cpp"
#include <time.h>

// FIFO file names
char IO_CMD_PIPE[1024], IO_ACK_PIPE[1024];

int fifo_cmd, fifo_ack; 

void CrashIO() {
    IOLOG(0,"Crashing the IO thread; sending IO_ERROR to client!\n");
    ioacknowledge ioack(IO_ERROR,-1);
    write(fifo_ack,&ioack, sizeof(ioacknowledge) );
    pthread_exit(NULL);    // Is this right?
}

// =====================================================================
// The actual read/write commands

#include "file.cpp"
#include "iolib.cpp"

ReadDirect * RD;
WriteDirect * WD;

void ReadIOR(iorequest *ior) {
    // Read the file, wait to complete.
    IOLOG(0,"Reading file %s\n", ior->filename);

    RD->BlockingRead( ior->filename, ior->memory, ior->sizebytes, ior->fileoffset);

    IOLOG(1,"Done reading file\n");
    IO_SetIOCompleted(ior->arena);
    return;
}

void WriteIOR(iorequest *ior) {
    IOLOG(0,"Writing file %s\n", ior->filename);
    // Write the file
    //ioassertf(FileExists(ior->filename)==0, 
    //	"File %s already exists; not the intended use of WriteFile.\n", ior->filename);
    //ioassertf(ior->fileoffset==0, 
    //	"WriteFile fileoffest = %d.  Non-zero values not supported.\n", ior->fileoffset);


    FILE * outfile = fopen(ior->filename,"wb");
    ioassertf(outfile != NULL,"Touching file %s failed\n", ior->filename);
    fclose(outfile);

    WD->BlockingAppend( ior->filename, ior->memory, ior->sizebytes);

    IOLOG(1,"Done writing file\n");
    if (ior->deleteafterwriting==IO_DELETE) IO_DeleteArena(ior->arena);
    return;
}


// ================================================================
#define NINSTRUCTIONS 65536
void *io_thread(void *data) {
    // This function runs as a stand-alone thread.
    // It must receive commands from the main program.
    // It inserts these commands into a set of ring buffers.
    // It then executes these commands in a priority order.

    if(data != NULL){
        const int assigned_core = *(int *) data;
        free(data);
        set_core_affinity(assigned_core);
        STDLOG(1, "IO thread assigned to core %d\n", assigned_core);
    }
    else{
        STDLOG(1, "IO thread not bound to core\n");
    }
    
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
		IOLOG(1,"Polling IO pipe: select() returned %d, wait_for_cmd %d\n", ret, wait_for_cmd);

		// if wait_for_cmd==1 then we should wait for a cmd to appear.
		// otherwise we would spin lock.  But if wait_for_cmd = 0, then 
		// we are just trying to empty the cmd pipe without locking up.
		if (wait_for_cmd == 1 || (ret>0 && FD_ISSET(fifo_cmd,&set))) {
		    IOLOG(1,"Reading from IO pipe\n");
		    // Read a command and place it in the buffers.
		    iorequest ior;
		    int nr = read(fifo_cmd, &ior, sizeof(iorequest));
		    assert(nr==sizeof(iorequest));
		    wait_for_cmd = 0;
		    // Put it in the buffer
		    if (ior.command==IO_READ) {
			IOLOG(1,"Received IO read request: file = %s, arena = %d, blocking = %d\n", 
			    ior.filename, ior.arena, ior.blocking);
		        if (ior.blocking==IO_BLOCKING) read_blocking.push(ior);
			    else read_nonblocking.push(ior);
		    } else if (ior.command==IO_WRITE) {
			IOLOG(1,"Received IO write request: file = %s, arena = %d, blocking = %d\n", 
			    ior.filename, ior.arena, ior.blocking);
		        if (ior.blocking==IO_BLOCKING) write_blocking.push(ior);
			    else write_nonblocking.push(ior);
		    } else if (ior.command==IO_QUIT) {
			IOLOG(1,"Received QUIT command\n");
		        quitflag = 1; break;  // We're not taking any more commands!
		    } else assert(ior.command==IO_READ);  // Explode
		} else fifo_not_empty = 0;
	    } while(fifo_not_empty);
	}

	IOLOG(1,"Attempting to execute an IO command\n");

	// Do one instruction, chosen by priority
	if (write_blocking.isnotempty()) {
	    IOLOG(1,"Starting blocking write\n"); 
	    iorequest ior = write_blocking.pop(); 
		
		BlockingIOWriteTime.Start();
	    WriteIOR(&ior);
		BlockingIOWriteTime.Stop();
		blocking_write_bytes += ior.sizebytes;
		
	    // Send an acknowledgement
	    ioacknowledge ioack(IO_WRITE,ior.arena);
	    write(fifo_ack,&ioack, sizeof(ioacknowledge) );
	    IOLOG(1,"IO_WRITE acknowledgement sent\n");
	} else if (read_blocking.isnotempty()) {
	    IOLOG(1,"Starting blocking read\n"); 
	    iorequest ior = read_blocking.pop();
		
		BlockingIOReadTime.Start();
	    ReadIOR(&ior);
		BlockingIOReadTime.Stop();
		blocking_read_bytes += ior.sizebytes;
		
	    // Send an acknowledgement
	    ioacknowledge ioack(IO_READ,ior.arena);
	    write(fifo_ack,&ioack, sizeof(ioacknowledge) );
	    IOLOG(1,"IO_READ acknowledgement sent\n");
	} else if (write_nonblocking.isnotempty()) {
	    IOLOG(1,"Starting nonblocking write\n"); 
	    iorequest ior = write_nonblocking.pop(); 
		
		NonBlockingIOWriteTime.Start();
	    WriteIOR(&ior);
		NonBlockingIOWriteTime.Stop();
		non_blocking_write_bytes += ior.sizebytes;
	} else if (read_nonblocking.isnotempty()) {
	    IOLOG(1,"Starting nonblocking read\n");
	    iorequest ior = read_nonblocking.pop();
		
		NonBlockingIOReadTime.Start();
	    ReadIOR(&ior);
		NonBlockingIOReadTime.Stop();
		non_blocking_read_bytes += ior.sizebytes;
	} else if (quitflag) break;
	else wait_for_cmd = 1;	
	    // We have no work to do, so we should read the cmd pipe to avoid spinlocking
    }

    // Terminate:  The IO is all done, so quit the pipes and remove the fifo files
    // Signal an acknowledgement.
    IOLOG(0,"Terminating the IO thread as planned.\n");
    ioacknowledge ioack(IO_QUIT,0);
    write(fifo_ack,&ioack, sizeof(ioacknowledge) );
    IOLOG(0,"IO_QUIT acknowledgement sent\n");
    close(fifo_cmd);
    close(fifo_ack);
    //pthread_exit(NULL);    // Is this right?
    return NULL;		// Quit the thread!
}


// ================================================================
// The commands below run in the client program (i.e. the compute thread).

void makeiopipes(void) {
    STDLOG(1,"Making io pipes.\n");
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
    STDLOG(1,"Deleting io pipe files\n");
    char cmd[1024];
    const char* rm_fmt = "rm -f %s";
    int ret = 0;
    if (FileExists(IO_ACK_PIPE)) {
        sprintf(cmd, rm_fmt, IO_ACK_PIPE);
        ret = system(cmd);
    }
    assertf(ret == 0, "Error removing pipe IO_ACK_PIPE=\"%s\"\n", IO_ACK_PIPE);
    
    if (FileExists(IO_CMD_PIPE)) {
        sprintf(cmd, rm_fmt, IO_CMD_PIPE);
        ret = system(cmd);
    }
    assertf(ret == 0, "Error removing pipe IO_CMD_PIPE=\"%s\"\n", IO_CMD_PIPE);
}

int io_cmd, io_ack;     // The pipe descriptor numbers
pthread_attr_t        io_attr1; // GLOBAL
pthread_t io_t1;

void IO_Initialize(char *logfn) {
    STDLOG(1,"Initializing IO\n");
    // Make the FIFO files
    pid_t pid = getpid();
    STDLOG(1,"Using pid %d in IO\n",pid);
    sprintf(IO_CMD_PIPE, "/tmp/iocmd_abacus.%d", pid);
    sprintf(IO_ACK_PIPE, "/tmp/ioack_abacus.%d", pid);
    remove_io_pipes();
    makeiopipes();

    iolog.open(logfn); 	// TODO: Probably need an error check


    int ramdisk = io_ramdisk_global;
    size_t diskbuffer = ((size_t) 1) << 29;
    RD = new ReadDirect(ramdisk, diskbuffer);
    WD = new WriteDirect(ramdisk,diskbuffer);
    
    // Determine core for IO thread
    int *io_core = (int*) malloc(sizeof(int));
    *io_core = 11;
    /*if(!free_cores.empty()){
        // look for a free core, starting from the last socket
        for(auto core_list
        
        // last free core on the last socket
        io_core = (int*) malloc(sizeof(int));
        *io_core = free_cores.back().back();
    }*/
    
    // Launch io_thread() as a separate thread
    int rc;
    rc = pthread_attr_init(&io_attr1);
    rc = pthread_attr_setdetachstate(&io_attr1, PTHREAD_CREATE_JOINABLE);
    int res = pthread_create(&io_t1, &io_attr1, io_thread, io_core);
    if(res) {
        printf("error %d\n", res);
    }
    STDLOG(1,"IO thread started!\n");

    // Open the pipes from the client side
    io_cmd = open(IO_CMD_PIPE, O_WRONLY);
    io_ack = open(IO_ACK_PIPE, O_RDONLY);
    STDLOG(1,"Done initializing IO\n");
    return;
}

void IO_Terminate() {
    STDLOG(1,"Terminating the IO thread\n");
    iorequest quitcmd; 
    quitcmd.command = IO_QUIT; 
    quitcmd.blocking = 1; 
    write(io_cmd, &quitcmd, sizeof(iorequest) ); 
    ioacknowledge quitcmdack; 
    read(io_ack, &quitcmdack, sizeof(ioacknowledge) );  
    assertf(quitcmdack.command==IO_QUIT,
    	"Error in IO acknowledgment\n"); 
    STDLOG(1,"Termination of IO thread confirmed\n");
    close(io_cmd);
    close(io_ack);
    remove_io_pipes();
    void * status =NULL;
    int rc = pthread_join(io_t1, &status);
    delete RD;
    delete WD;
    assert(!rc);
    STDLOG(1,"Termination of IO complete\n");
}

// Here are the actual interfaces for writing an arena

void ReadFile(char *ram, uint64 sizebytes, int arena,
	    const char *filename, off_t fileoffset, int blocking) {

    STDLOG(1,"Using IO_thread module to read file %f\n", filename);
    iorequest ior(ram, sizebytes, filename, IO_READ, arena, fileoffset, 0, blocking);
    write(io_cmd, &ior, sizeof(iorequest));
    
    if (blocking) {
        ioacknowledge ioack;
	read(io_ack, &ioack, sizeof(ioacknowledge));
	assertf(ioack.arena == arena, "Error in IO acknowledgement for arena %d = %d\n",
		arena, ioack.arena);
	STDLOG(1,"Blocking IO returned\n");
    } else {
        STDLOG(1,"Non-blocking IO requested\n");
    }
}

void WriteFile(char *ram, uint64 sizebytes, int arena, 
	    const char *filename, off_t fileoffset, int deleteafter, int blocking) {

    STDLOG(1,"Using IO_thread module to write file %f\n", filename);
    iorequest ior(ram, sizebytes, filename, IO_WRITE, arena, fileoffset, deleteafter, blocking );
    write(io_cmd, &ior, sizeof(iorequest));
    
    if (blocking) {
        ioacknowledge ioack;
	read(io_ack, &ioack, sizeof(ioacknowledge));
	assertf(ioack.arena == arena, "Error in IO acknowledgement for arena %d = %d\n",
		arena, ioack.arena);
	STDLOG(1,"Blocking IO returned\n");
    } else {
        STDLOG(1,"Non-blocking IO requested\n");
    }
}


