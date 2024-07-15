/* arenaallocator.cpp
 *
 * \brief Our low-level routine for handling memory.
 *
 * This allocates 
 * and deallocates arenas and provides a way to mark whether 
 * those arenas are ready for use.  Because I/O is being done 
 * by a separate thread, we have to be able to mark when I/O 
 * is completed.
 *
 * The ArenaAllocator class is exclusively accessed through
 * the SlabBuffer class.  SlabBuffer has an instance of
 * ArenaAllocator as a member variabled; i.e. SlabBuffer
 * no longer inherits from ArenaAllocator.
 * 
 * This file and class were formerly referred to as LinearBuffer.
*/


#ifndef INCLUDE_LB
#define INCLUDE_LB

#include "threadaffinity.h"

//#include "tbb/concurrent_queue.h"
#include "simple_concurrent_queue.cpp"

#include <mutex>
#include <atomic>

#define GUARDSIZE 4096
// Allocate this many bytes of guard space around the buffer
// This is used to check for overflowing writes which is useful
// but not strictly necessary

typedef struct {
    uint64 allocated_size, usable_size, max_usable_size, start_offset;
    int present;    // This will be 1 if the arena is available.  2 if allocated, but not available for normal work.
    int IsIOCompleted;
    int shared_mem;  // does this arena reside in shared memory?
    char *addr;
} arenainfo;

enum RamdiskArenaType { RAMDISK_NO,
                        RAMDISK_READSLAB,
                        RAMDISK_WRITESLAB,
                        RAMDISK_AUTO,
                        RAMDISK_AUTO_READSLAB,  // RAMDISK_READSLAB if on ramdisk else RAMDISK_NO
                        RAMDISK_AUTO_WRITESLAB  // RAMDISK_WRITESLAB if on ramdisk else RAMDISK_NO
};


class ArenaAllocator {
public:
    ArenaAllocator(int maximum_number_ids, int use_disposal_thread, int disposal_thread_core);
    ~ArenaAllocator(void);

    void report_peak(int force=0);
    void report_current();
    
    int64 total_allocation;
    std::atomic<int64> total_shm_allocation;  // atomic because the disposal thread may modify
    STimer ArenaMalloc;
    int numalloc, numreuse;
    int num_shm_alloc;
	
	
	float ArenaFree_elapsed;
    STimer DisposalThreadMunmap;

   // PTimer *ArenaFree;

    void Allocate(int id, uint64 s, int reuseid, int ramdisk = RAMDISK_NO, const char *ramdisk_fn=NULL);
    void DeAllocateArena(int id, int reuseID);
    void ResizeArena(int id, uint64 s);

    inline int IsArenaPresent(int id) { 
        // Returns 1 if the Arena is Allocated, 0 otherwise
        assert(id>=0 && id<maxids); 
        return (arena[id].present==1)?1:0; 
    }

    inline void SetIOCompletedArena(int id) { 
        // Declare the Arena completed with I/O
        assert(IsArenaPresent(id));
        arena[id].IsIOCompleted = 1; 
    }

    inline char *GetArenaPtr(int id) {
        // Returns a pointer to the usable area of the Arena
        assert(IsArenaPresent(id));
        return arena[id].addr+arena[id].start_offset; 
    }
    inline int IsIOCompleted(int id) { 
        // Returns 1 if the Arena has completed I/O and is ready for use.
        assert(IsArenaPresent(id));
        return arena[id].IsIOCompleted; 
    }

    uint64 ArenaSizeBytes(int id) { 
        // Return the usable size of the Arena
        assert(IsArenaPresent(id));
        return arena[id].usable_size; 
    }

    int ArenaRamdiskType(int id){
        return arena[id].shared_mem;
    }

    inline void MarkArenaUnavailable(int id) {
        /* This should be used sparingly.  It makes it so that calls
        to IsArenaPresent() will return False.  Pretty much the only 
        thing that will work after that is a call to DeAllocate() */
        /* One should *never* make the ReuseArena unavailable! */
        arena[id].present = 2;
        return;
    }


private:

    int maxids;
    arenainfo *arena;
    uint64 peak_alloc;
    uint64 peak_shm_alloc;
    int have_new_peak;
    std::mutex lb_mutex;
	std::mutex lb_freemutex;
	
    pthread_t disposal_thread;

    void DiscardArena(int id);

    void ResetArena(int id) {
        // This just nulls out the arena
            assert(id>=0 && id<maxids); 
        arena[id].present = 0;
        arena[id].IsIOCompleted = 0;
        arena[id].addr = NULL;
        arena[id].allocated_size = 0;
        arena[id].usable_size = 0;
        arena[id].max_usable_size = 0;
        arena[id].start_offset = 0;
        arena[id].shared_mem = 0;
    }

    void DumpArena(int id, int less) {
        assert(!arena[id].shared_mem);  // not implemented

        // Be wary of calling this from the IO thread; STDLOG is not thread-safe.
        char start[GUARDSIZE+1], end[GUARDSIZE+1];
        memcpy(start, arena[id].addr, GUARDSIZE); start[GUARDSIZE] = '\0';
        memcpy(end, arena[id].addr+GUARDSIZE+arena[id].usable_size, GUARDSIZE); end[GUARDSIZE] = '\0';
        if (less) { STDLOG(2,"Arena %d: %d %d %d %d %d %p\n", 
            id,
            arena[id].present,
            arena[id].IsIOCompleted,
            arena[id].allocated_size, 
            arena[id].usable_size, 
            arena[id].start_offset, 
            (void *) (arena[id].addr)); 
        } else STDLOG(2,"Arena %d: %d %d %d %d %d %p %s %s\n", 
            id,
            arena[id].present,
            arena[id].IsIOCompleted,
            arena[id].allocated_size, 
            arena[id].usable_size, 
            arena[id].start_offset, 
            (void *) (arena[id].addr), 
            start, end);
    }

    void SetGuardStart(int id) {
        assert(!arena[id].shared_mem);  // not supported

        char *c = arena[id].addr;
        for (int i=0; i<GUARDSIZE; i++) c[i] = 'a'+i%16;
    }
    /// Set a Guard at the end of the usable size
    void SetGuardEnd(int id) {
        assert(!arena[id].shared_mem);  // not supported

        //? char *c = arena[id].addr+arena[id].allocated_size-1;
        //? for (int i=0; i<size; i++, c--) *c = 'a'+i%16;
        char *c = arena[id].addr+GUARDSIZE+arena[id].usable_size;
        for (int i=0; i<GUARDSIZE; i++, c++) *c = 'a'+i%16;
    }
    int CheckGuardStart(int id) {
        assert(!arena[id].shared_mem);  // not supported

        // Return 0 if the guard is intact.
        // This routine can be invoked by the IO thread!
        // STDLOG(1,"Guard check id=%d, size=%d\n", id, size);
        char *c = arena[id].addr;
        for (int i=0; i<GUARDSIZE; i++) if (c[i] != 'a'+i%16) {
            STDLOG(0,"About to fail on id = %d, i=%d\n", id, i);
            return i+1;
        }
        return 0;
    }
    int CheckGuardEnd(int id) {
        assert(!arena[id].shared_mem);  // not supported

        // Return 0 if the guard is intact.
        // This routine can be invoked by the IO thread!
        // STDLOG(1,"Guard check id=%d, size=%d\n", id, size);
        //?char *c = arena[id].addr+arena[id].allocated_size-1;
        //?for (int i=0; i<size; i++, c--) if (*c != 'a'+i%16) {
        char *c = arena[id].addr+GUARDSIZE+arena[id].usable_size;
        for (int i=0; i<GUARDSIZE; i++, c++) if (*c != 'a'+i%16) {
            STDLOG(0,"About to fail on id = %d, i=%d\n", id, i);
            return i+1;
        }
        return 0;
    }

    uint64 ArenaFullSize(int id) { 
        // Return the maximum available size of the Arena
        assert(IsArenaPresent(id));
            return arena[id].max_usable_size;
    }


    // Fields related to the disposal thread
    int use_disposal_thread, disposal_thread_core;

    static void *start_thread(void *AA_obj){
        ArenaAllocator *_AA = (ArenaAllocator *) AA_obj;
        _AA->DisposalThreadLoop(_AA->disposal_thread_core);
        return NULL;
    }

    void DisposalThreadLoop(int core);

    struct disposal_item {
        void *addr;
        size_t size;
    };
    //tbb::concurrent_bounded_queue<struct disposal_item> disposal_queue;
    SimpleConcurrentQueue<struct disposal_item> disposal_queue;
};


ArenaAllocator::ArenaAllocator(int maximum_number_ids, int _use_disposal_thread, int _disposal_thread_core) {
    maxids = maximum_number_ids;
   // ArenaFree = new PTimer(maxids);
    arena = new arenainfo[maxids];
    total_allocation = 0;
    total_shm_allocation = 0;
    numalloc = numreuse = 0;
    num_shm_alloc = 0;
    ArenaFree_elapsed = 0.;
    peak_alloc = 0;
    peak_shm_alloc = 0;
    have_new_peak = 0;
    for(int i=0;i<maxids;i++)
        ResetArena(i);

    use_disposal_thread = _use_disposal_thread;
    disposal_thread_core = _disposal_thread_core;

    if(use_disposal_thread){
        assert(pthread_create(&disposal_thread, NULL, start_thread, this) == 0);
    }
}

ArenaAllocator::~ArenaAllocator(void) { 
    for(int i=0;i<maxids;i++) 
            if (arena[i].present) {
                STDLOG(0,"Arena %d was not deallocated before arena allocator destructor was called.\n", i);
                // This used to be a major concern, but it's not in the parallel code
                // fprintf(stderr,"Arena %d was not deallocated before arena allocator destructor was called.\n", i);
                assertf(arena[i].addr!=NULL, "Arena %d is present, but NULL\n", i);
                DiscardArena(i); 
            }
            else
                ResetArena(i);
    delete[] arena;
   // delete ArenaFree;

    if(use_disposal_thread){
        disposal_queue.push((struct disposal_item){NULL, 0});

        assert(pthread_join(disposal_thread, NULL) == 0);
        assertf(disposal_queue.empty(), "Disposal queue not empty!\n");
    }
}

void ArenaAllocator::report_peak(int force){
    if(force || have_new_peak){
        STDLOG(0,"Peak arena allocations (regular + shared) was %.3g GB\n", peak_alloc/1024./1024/1024);
        STDLOG(0,"Peak arena allocations (shared only) was %.3g GB\n", peak_shm_alloc/1024./1024/1024);
        STDLOG(0,"Executed %d arena fresh allocations, %d arena reuses\n", numalloc, numreuse);
        have_new_peak = 0;
    }
}

void ArenaAllocator::report_current(){
    STDLOG(1, "Current arena allocations (normal + shared) total %.3g GB\n", total_allocation/1024/1024/1024);
    STDLOG(1, "Current arena allocations (shared only) total %.3g GB\n", total_shm_allocation/1024/1024/1024);
}

#define LB_OVERSIZE 1.01

/// Allocating a new arena, reusing an old buffer if possible.
void ArenaAllocator::Allocate(int id, uint64 s, int reuseID, int ramdisk, const char *ramdisk_fn) {
    lb_mutex.lock();
    
    assertf(arena[id].present==0, "Error: Asking for Allocation of arena %d that already exists!\n", id);   // This is always a bad idea
    assertf(id < maxids, "Error: Asking for Allocation of arena %d, but maxids is %d\n", id, maxids);
    
    size_t ss;
    ss = sizeof(char);
    ss *=  s;

    std::string spath;

    STDLOG(4,"AA slab=%d, size %d, ramdisk %d\n", id, s, ramdisk);

    switch(ramdisk){
        case RAMDISK_NO:  // normal arena allocation

            // If no ramdisk allocation was requested, better not have received a path!
            assert(ramdisk_fn == NULL || strnlen(ramdisk_fn,1) == 0);

            // Let's check if the reuseID is available.
            // The reuseID is guaranteed to be unused.
            if (IsArenaPresent(reuseID) && ss<=arena[reuseID].max_usable_size) {
                // Yes, the arena allocation exists and is big enough
                arena[id] = arena[reuseID];
                ResetArena(reuseID);
                arena[id].usable_size = ss;
                arena[id].IsIOCompleted = 0;  
                numreuse++;
            } else {
                // No, we need to make a new one
                arena[id].usable_size = ss;
                arena[id].max_usable_size = ss*LB_OVERSIZE;
                arena[id].allocated_size = arena[id].max_usable_size + 2*GUARDSIZE;
                ArenaMalloc.Start();
                assertf(posix_memalign((void **) &(arena[id].addr), PAGE_SIZE, arena[id].allocated_size) == 0,
                        "posix_memalign failed trying to allocate %d bytes\n", arena[id].allocated_size);
                arena[id].present = 1;  
                ArenaMalloc.Stop();
                assert(arena[id].addr!=NULL);        // Crash if the malloc failed
                numalloc++;
                total_allocation += arena[id].allocated_size;
                arena[id].IsIOCompleted = 0;  
                arena[id].start_offset = GUARDSIZE;
                SetGuardStart(id);
            }

            // We need to reset the End Guard in either case.
            SetGuardEnd(id);
            break;

        case RAMDISK_WRITESLAB:  // create a new shared memory allocation
        case RAMDISK_READSLAB:  // attach an existing shared memory allocation

            {
            // If a ramdisk allocation was requested, must have received the path
            assert(ramdisk_fn != NULL);
            assert(strnlen(ramdisk_fn,1) > 0);
            STDLOG(2,"Mapping arena id %d from shared memory\n", id);

            // Shared memory arenas:
            // 1) do not have guard space
            // 2) are never "oversized"
            // 3) don't recycle old allocations

            arena[id].usable_size = ss;
            arena[id].max_usable_size = ss;
            arena[id].allocated_size = arena[id].max_usable_size;
            arena[id].shared_mem = ramdisk;
            arena[id].present = 1;

            ArenaMalloc.Start();  // just use the same timer so we notice if something's off the rails

            // Use open() to create the shm fd instead of shm_open()
            // They do the same thing, but open() takes a full path and allows subdirectories

            // Remove the shared memory handle if it exists
            // Technically, this is only necessary if we're overwriting the state
            // But then we need to know what slab types are state and what are convolution state
            // Our DIO library always does this and we've never had problems
            if(ramdisk == RAMDISK_WRITESLAB)
                unlink(ramdisk_fn);

            // We could think about letting the ArenaAllocator decide whether to make a new allocation
            // based on the (non)existence of a slab
            // But that would only work if
            // 1) ArenaAllocator knows about state vs convolution state slabs and which to overwrite
            // or 2) we always delete immediately after read (i.e. steps enforce overwriting)

            int shm_fd_flags = ramdisk == RAMDISK_WRITESLAB ? (O_RDWR | O_CREAT | O_EXCL) : O_RDWR;
            int fd = open(ramdisk_fn, shm_fd_flags, S_IRUSR | S_IWUSR);
            assertf(fd != -1, "Failed to open shared memory file at \"%s\"\n", ramdisk_fn);

            if(ramdisk == RAMDISK_WRITESLAB){
                // set the size
                int res = ftruncate(fd, ss);
                assertf(res == 0, "ftruncate on shared memory ramdisk_fn = %s to size = %d failed\n", ramdisk_fn, ss);
            } else {
                // check that the size is as expected
                struct stat shmstat;
                int res = fstat(fd, &shmstat);
                assertf(res == 0, "fstat on shared memory ramdisk_fn = %s\n", ramdisk_fn);
                assertf(shmstat.st_size == ss, "Found shared memory size %d; expected %d (ramdisk_fn = %s)\n", shmstat.st_size, ss, ramdisk_fn);
            }

            if(ss > 0){  // zero-length mmap is prohibited
                // map the shared memory fd to an address
                int mmap_flags = ramdisk == RAMDISK_WRITESLAB ? (PROT_READ | PROT_WRITE) : (PROT_READ | PROT_WRITE);  // the same
                arena[id].addr = (char *) mmap(NULL, ss, mmap_flags, MAP_SHARED, fd, 0);
                assertf((void *) arena[id].addr != MAP_FAILED, "mmap shared memory from fd = %d of size = %d failed\n", fd, ss);
                assertf(arena[id].addr != NULL, "mmap shared memory from fd = %d of size = %d failed\n", fd, ss);
            } else {
                arena[id].addr = NULL;
            }

            int res = close(fd);
            assertf(res == 0, "Failed to close fd %d\n", fd);

            num_shm_alloc++;

            ArenaMalloc.Stop();
            
            // Either the shared memory already exists or it was just allocated
            total_shm_allocation += arena[id].allocated_size;
            total_allocation += arena[id].allocated_size;

            arena[id].IsIOCompleted = 0;  
            arena[id].start_offset = 0;  //GUARDSIZE;
            break;
            }  // Need this to contain the variable initializations for this branch
        default:
            QUIT("Illegal value %d for ramdisk\n", ramdisk);
    }

    if (total_allocation > peak_alloc){
        peak_alloc = total_allocation;
        have_new_peak = 1;
    }
    if (total_shm_allocation > peak_shm_alloc){
        peak_shm_alloc = total_shm_allocation;
        have_new_peak = 1;
    }

    lb_mutex.unlock();
}

/// This discards an arena, freeing and resetting, no questions asked
void ArenaAllocator::DiscardArena(int id) {
	
	STimer ArenaFree; 
    //ArenaFree->Start(id);
	ArenaFree.Clear(); 
	ArenaFree.Start();

    if(arena[id].shared_mem){
        // This might free the underlying memory if the /dev/shm handle has been deleted
        // TODO: is there a way to check, just for logging/accounting purposes?

        if(arena[id].allocated_size > 0){
            if(use_disposal_thread){
                STDLOG(3, "Pushing %d to discard thread...\n", id);
                struct disposal_item di = {arena[id].addr, arena[id].allocated_size};
                disposal_queue.push(di);
                STDLOG(3, "Done pushing %d.\n", id);

                // If using munmap thread, it is responsible for decrementing total_shm_allocation
            } else {
                int res = munmap(arena[id].addr, arena[id].allocated_size);
                assertf(res == 0, "munmap failed\n");
                total_shm_allocation -= arena[id].allocated_size;  // only if not using thread
            }
            
        } else {
            assertf( arena[id].addr == NULL , 
            "Shared-memory arena %d of zero size has a non-null pointer %p!\n", id, arena[id].addr);
        }
    }
    else {
        assertf( arena[id].addr != NULL , 
            "Arena %d requested for deallocation, but it points to NULL\n", id); 

        // DumpArena(id,1);
        free(arena[id].addr);
    }

    arena[id].present = 0;
    total_allocation -= arena[id].allocated_size;
    //ArenaFree->Stop(id);
	
	ArenaFree.Stop();
	
	lb_freemutex.lock();
	ArenaFree_elapsed += ArenaFree.Elapsed(); 
	lb_freemutex.unlock();
	
    ResetArena(id);
}

void ArenaAllocator::DeAllocateArena(int id, int reuseID) {
    // This will deallocate and reset a given arena.
    // It is illegal to call this on an arena that hasn't been allocated.
    lb_mutex.lock();
    assert(id >= 0 && id < maxids); 
    assertf(arena[id].present>0, "Arena %d requested for deallocation, but it doesn't exist\n", id ); 
    arena[id].present = 1;   // If it was marked unavailable, fix that before it propagates

    // If we're asked to discard the reuse slab or a shared memory slab, just do it
    if (id == reuseID || arena[id].shared_mem) {
        DiscardArena(id);
        lb_mutex.unlock();
        return;
    }

    // Check that nothing has gone wrong
    assertf( CheckGuardStart(id)==0, "Arena %d failed its check of GuardStart\n", id);
    assertf( CheckGuardEnd(  id)==0, "Arena %d failed its check of GuardEnd\n", id);

    if (IsArenaPresent(reuseID)) {
        // We already have an arena available to reuse
        if (arena[id].max_usable_size<=arena[reuseID].max_usable_size) {
            // The reuse buffer is bigger, so discard the current one
            // TODO: Might be better to use <= here
            DiscardArena(id);
        } else {
            // Discard the reuse buffer and move the current storage over
            DiscardArena(reuseID);
            arena[reuseID] = arena[id];
            assertf(arena[reuseID].addr == arena[id].addr,
                "We set the pointers in arenas %d and %d equal, but it didn't work\n", reuseID, id);
            ResetArena(id);
        }
    } else {
        // We don't have a reuse buffer.  Move this one over.
        arena[reuseID] = arena[id];
        ResetArena(id);
    }

    assertf(IsArenaPresent(id)==0, "We tried to Deallocate arena %d, but it's still here.\n", id);
    lb_mutex.unlock();
}

void ArenaAllocator::ResizeArena(int id, uint64 s) {
    // We can resize an arena, but only up to the currently allocated size!
    // The primary use here is to shrink an arena to limit the output amount.
    // No memory is actually freed for other uses
    assert(IsArenaPresent(id));
    lb_mutex.lock();

    // TODO: I think this could be done with a combination of ftruncate and mremap()
    //  but would require saving the shared memory path
    assertf(!arena[id].shared_mem, "Arena %d allocated in shared memory; resizing currently not supported.  Tried to resize from %d to %d bytes.\n",
            id, arena[id].usable_size, s);

    assertf(s<=arena[id].max_usable_size, 
            "Can't Resize to a larger size! %d new %d old %d\n", 
            id, s, arena[id].max_usable_size);
    arena[id].usable_size = s;
    SetGuardEnd(id);
    lb_mutex.unlock();
    return;
}


void ArenaAllocator::DisposalThreadLoop(int core){
    struct disposal_item di;

    STDLOG(1, "Starting arena allocator disposal thread on core %d\n", core);

    if(core >= 0)
        set_core_affinity(core);

    int n = 0;
    while(true){
        disposal_queue.pop(di);  // blocking

        // NULL is a safe termination flag because munmap(NULL) is not safe!
        if(di.addr == NULL)
            break;

        if(di.size > 0){
            DisposalThreadMunmap.Start();
            int res = munmap(di.addr, di.size);
            DisposalThreadMunmap.Stop();
            total_shm_allocation -= di.size;
            STDLOG(2, "Just munmap()'d %d bytes. %d items left on queue\n", di.size, disposal_queue.size());
            assertf(res == 0, "munmap failed\n");
            n++;
        }
    }

    STDLOG(1, "Terminating arena allocator disposal thread.  Executed %d munmap()s.\n", n);
}

/**
 * This utility writes to the log the number of bytes allocated by malloc, and the heap size in use
 * by the program.  The difference is an estimate of the fragmentation.  Or it could simply
 * be due to under-use!  I'm not sure there's a way to tell.
 * It doesn't have anything to do with ArenaAllocator in particular; it just seemed this function
 * should go somewhere related to memory management.
 * Reference: https://gperftools.github.io/gperftools/tcmalloc.html
 */
#if defined(HAVE_LIBTCMALLOC_MINIMAL) || defined(HAVE_LIBTCMALLOC)
#include "gperftools/malloc_extension.h"
#endif

void ReportMemoryAllocatorStats(){
#if defined(HAVE_LIBTCMALLOC_MINIMAL) || defined(HAVE_LIBTCMALLOC)
    size_t bytes_allocated = 0;
    MallocExtension::instance()->GetNumericProperty("generic.current_allocated_bytes", &bytes_allocated);

    size_t bytes_total = 0;
    MallocExtension::instance()->GetNumericProperty("generic.heap_size", &bytes_total);

    STDLOG(2, "%.3g GiB held by mallocs; %.3g GiB held from system by allocator\n", bytes_allocated/1024./1024/1024, bytes_total/1024./1024/1024);
    STDLOG(2, "\t%.3g GiB (%.1f%%) held by allocator but not in use\n", (bytes_total - bytes_allocated)/1024./1024/1024, 100.*(bytes_total - bytes_allocated)/bytes_total);

    // This dumps a human-readable summary of the current allocator state to the log
    char logstr[2048];
    MallocExtension::instance()->GetStats(logstr, 2048);
    
    std::string line;
    std::istringstream logstream(logstr);
    // Split up the string on newlines so parsers don't break
    while(std::getline(logstream, line))
        STDLOG(2, "%s\n", line.c_str());

    // Malloc histogram
    // if you link against tcmalloc_minimal_debug (warning: slow!)
    // this will output a log2 histogram of mallocs by size
    int blocks;
    size_t total;
    int nHistBins = 64;
    int hist[nHistBins];
    MallocExtension::instance()->MallocMemoryStats(&blocks, &total, hist);
    if (blocks > 0){
        size_t startbytes = 0, endbytes = 1;
        const char *suffix = "  B";
        for(int i = 0; i < nHistBins; i++){
            if (i >= 11)
                suffix = "KiB";
            if (i >= 21)
                suffix = "MiB";
            if (i >= 31)
                suffix = "GiB";
            if (i >= 41){
                assertf(hist[i] == 0, "%d allocations found in bin %d.  Too big!\n", hist[i], i);
                continue;
            }

            
            if ((i % 10 == 1) && (i > 1)){
                startbytes >>= 10;
                endbytes >>= 10;
            }

            char rightalignnum[16];  // I think STDLOG has a bug; this is the workaround
            sprintf(rightalignnum, "%4zu", startbytes);
            STDLOG(3, "%s %s -- %4d %s: %d\n", rightalignnum, suffix, endbytes, suffix, hist[i]);
            startbytes = endbytes;
            endbytes <<= 1;
        }
    }
#endif
}

void ReleaseFreeMemoryToKernel(){
#if defined(HAVE_LIBTCMALLOC_MINIMAL) || defined(HAVE_LIBTCMALLOC)
    ReleaseFreeMemoryTime.Start();
    MallocExtension::instance()->ReleaseFreeMemory();
    ReleaseFreeMemoryTime.Stop();
#endif
}

#endif // INCLUDE_LB
