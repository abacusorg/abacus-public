/*** This is a set of classes to accumulate and store a cell-based listing
of objects in a Slab.  The challenge is that we don't know the lengths
as we're creating them, yet must support multi-threading.  The threading
is assumed to be strictly by Pencil, so we do buffering and indexing on
that basis.  We do not guarantee that the buffer is contiguous nor that 
the pencils are in order.  Each pencil is guaranteed to have contiguous
data, ordered by cell.

We limit the number of elements in a Pencil to fit in size int 
(note, this is the number of elements, not the number of bytes).

This is all templated.  

USAGE: 
Initialize SlabAccum<T> s 
Call s.setup(cpd, zwidth, maxsize);
maxsize should be some estimate of the size requirements per slab, 
but don't sweat it.

When starting a Pencil, call PencilAccum<T> *p = s.StartPencil(pencilnum).
This will detect the thread number and give you an object to append to.
Call p->append(val) to do so.
Call p->FinishCell() after each cell, and p->FinishPencil() after the pencil is done.

It is required that a Pencil be processed by a single thread and
that the cells be processed in order [0,cpd).

You shouldn't access values until FinishPencil is called.
Then s[j] will return the PencilAccum<T>, 
s[j][k] will return the CellPtr<T>,
and s[j][k][n] will return an individual element

2D: SlabAccum knows nothing about subslabs except that the z-dimension
has length `zwidth`. So one must use k indices relative to the lowest z cell,
which will usually be `node_z_start_ghost`.

s.get_slab_size() will return the total elements in the slab.
s.get_slab_bytes() will return the total bytes in the slab.
s[j]._size will return the number in a pencil.

s.copy_to_ptr(p) will copy the whole slab contents, in pencil order, to a 
contiguous buffer at p, which must be pre-allocated to size().  

s.build_pstart() builds an array s.pstart[j] that contains the starting offset for each pencil,
which is useful for indexing the contiguous buffer..

*/


#ifdef TEST
#include "test_driver_header.cpp"
STimer SlabAccumFree;
#endif

// -----------------------------------------------------------------

/// This is a lightweight indexing of the cells within a pencil
class CellAccum {
  public:
    int start;    // The starting index relative to the PencilAccum *data.
    int _size;	  // The size of the cell
    inline int size() { return _size; }
};

/// We give this simple class as a way to return the pointer to and size of
/// the cell-based data.
template <class T>
class CellPtr {
  public:
    T *start;
    int _size;
  
    inline CellPtr(T *data, CellAccum &c) { 
	start = data+c.start;
	_size = c._size;
    }
    ~CellPtr() { };

    inline int size() { return _size; }

    // Provide a [] to get to an individual element
    inline T operator[](int n) { return start[n]; }

    // Also provide a quick function to get a pointer to an individual element
    inline T *ptr(int n) { return start+n; }
};

// -----------------------------------------------------------------

/// This is just a light-weight vector, forming a linked list.
/// If we run out of space, then we make new space and copy [start,maxsize)
/// into the new space.
/// This buffer cannot exceed a 31-bit size.
template <class T>
class SlabAccumBuffer {
  private:
    T *data;
    SlabAccumBuffer<T> *previous;
    int size, maxsize, start;
    char padding[CACHE_LINE_SIZE-28];    // Just to make this cacheline safe

  public:
    SlabAccumBuffer() { 
	// We use a null constructor, so don't forget to use setup()
    	data = NULL; previous = NULL; 
	size = maxsize = start = 0;
    }
    void slabfree() {
        // We've been asked to free all of the lists.
	// We provide this destructor explicitly, so it can be called in a loop
    	if (data!=NULL) {
	    #ifdef TEST
	    fmt::print("Freeing memory at {:p}\n", data);
	    #endif
	    free(data); 
	    data = NULL;
	}
	if (previous!=NULL) {
	    // previous->slabfree();  // This is implied by the destructor
	    delete previous;   // We made it, so we need to get rid of it.
	    previous = NULL;
	}
    }
    ~SlabAccumBuffer() { 
        slabfree();
    }

    bool isNULL() { return data==NULL; }

    #define MAXSlabAccumBuffer 0x3fffffff
	// Cap at a billion elements, but we don't want 2*number to overflow int.
	// This is a huge amount, since this is for a pencil, not a slab.

    // Call setup to initialize the buffer
    void setup(int _maxsize) {
	// fmt::print("Allocating to size {:d}\n", _maxsize); fflush(NULL);
	size_t bytes = 4096;
	while (bytes<_maxsize*sizeof(T)) bytes <<= 1;
	    // Shift up to the next factor of 2 in bytes.
	    // Then figure out how many objects this is.
	maxsize = bytes/sizeof(T);
	if (maxsize>MAXSlabAccumBuffer) maxsize = MAXSlabAccumBuffer;
		// Save user from themselves
	int ret = posix_memalign((void **)&data, PAGE_SIZE, bytes); assert(ret==0);
	#ifdef TEST
	fmt::print("Allocated {:d} at position {:p}\n", maxsize, data);
	#endif
    }

    // If we already have data to load in, do it here.
    void load_from_ptr(void *p, int _size) {
	setup(_size);
	memcpy(data, p, sizeof(T)*_size);
	size = maxsize;
	start = 0;
    }

    // When we start work on a pencil, declare its starting point,
    // in case we need to move the data.
    void set_pencil_start() { start = size; }

    inline void append(T value) {
	if (size<maxsize) {
	    data[size++] = value;
	    return;
	} else {
	    // We've exhausted this buffer
	    // Save off the current buffer
	    // fmt::print("Exceeded a buffer: [{:d},{:d})\n", start, maxsize);
	    SlabAccumBuffer<T> *p = previous;
	    previous = new SlabAccumBuffer<T>;
	    previous->previous = p;
	    previous->data = data;
	    previous->size = start;
	    previous->maxsize = maxsize;
	    // Make a new buffer
	    if (2*start<maxsize) {
	    	// Apparently, this pencil is already >half the buffer size,
		// so let's make the next one bigger.
		if (maxsize>MAXSlabAccumBuffer/2) {
		    // Don't overflow an int
		    maxsize = MAXSlabAccumBuffer-1;
		    assert(size-start<maxsize);
		} else maxsize = 2*maxsize;
	    }
	    setup(maxsize);
	    // Copy the old data and set up the new buffer.
	    // fmt::print("Copying {:d} from old to new\n", size-start);
	    memcpy(data, previous->data+start, (size-start)*sizeof(T));
	    size = (size-start);
	    start = 0;
	    // Finally, append the new value
	    data[size++] = value;
	    return;
	}
    }
    
    // Need to provide a way to get to the start point, since it 
    // may have moved!
    inline T *get_pencil_start() { return data+start; }
    inline int get_pencil_size() { return size-start; }

    uint64 get_buffer_size() {
	// Return the total number of entries in all of the buffers
        uint64 total = size;
	SlabAccumBuffer<T> *p = previous;
	while (p!=NULL) {
	    total += p->size;
	    p = p->previous;
	}
	return total;
    }
};

// -----------------------------------------------------------------

/// The PencilAccum class is the way we want to interact with the work in
//  a pencil.
template <class T>
class PencilAccum {
    // These two variables are left in an ill-defined state after FinishPencil() is called.
    CellAccum *thiscell;     // Point to the current cell
    int thisstart;	// Start index of current cell

  public:
    SlabAccumBuffer<T> *buffer;    // Pointer to the buffer we're using
    CellAccum *cells;    // Points to an external allocation [0,zwidth)
    T *data;		// Eventually we point to the data, location of start index 0
    int _size;		// The number of elements in this pencil
    			// Important: this isn't set until the end!
			// use get_pencil_size() before FinishPencil;
			// use _size after

    // Null constructor, and nothing to destroy/free
    PencilAccum() {
    }
    ~PencilAccum() {
    }

    // Note that get_pencil_size() should only be used before FinishPencil
    inline int get_pencil_size() {
         return buffer->get_pencil_size();
    }

    inline int cellsize(int k) {
	// Get the size of a cell
        return cells[k].size();
    }
    inline T *cellstart(int k) {
	// Get a pointer to the start of the list for a cell
	return data+cells[k].start;
    }
    // Provide [] to get to the Cell-based values
    inline CellPtr<T> operator[](int j) { 
        return CellPtr<T>(data,cells[j]);
    }

    inline void StartPencil() {
	// Call this function at the start of the Pencil
        buffer->set_pencil_start();
	thiscell = cells;   // Start back at the beginning
	thisstart = buffer->get_pencil_size();
	return;
    }

    inline void append(T item) {
		// TODO: emplace semantics
	// Call this function to add one item to the current cell
	buffer->append(item);
	return;
    }
    
    inline void FinishCell() {
	// Call this function at the end of each Cell
	// Record this cell
        thiscell->start = thisstart;
	int nextstart = buffer->get_pencil_size();
	thiscell->_size = nextstart-thisstart;
	// Advance to the next one
	thisstart = nextstart;
	thiscell++;
	return;
    }

    inline void FinishPencil() {
	// The data may have moved, so we need to get the final location.
        data = buffer->get_pencil_start();
	_size = buffer->get_pencil_size();
	return;
    }
};

// --------------------------------------------------------------

/// This class allows us to accumulate objects of a particular type
/// in a cell and pencil-based manner.  It supports multi-threaded
/// work by spreading the threads over the pencils.  At the end, one
/// has a clean indexing scheme.

template <class T>
class SlabAccum {
    int cpd;                // number of y pencils
	int zwidth;             // number of z cells on this node
    CellAccum *cells;		// Will be allocated to [0,cpd*zwidth)

    // We have an accumulation buffer for each thread.
    int maxthreads;
    SlabAccumBuffer<T> *buffers;     // Vector [0,maxthreads)

  public:
    // We index these into Pencils
    PencilAccum<T> *pencils;	// Will be allocated [0,cpd)
    uint64 *pstart;   // Will be allocated [0,cpd)
	    // These indicate the index starts for each pencil, if one were to
	    // repack the pencils into a contiguous buffer.

    // Provide [] to get to the Pencil-based values
    inline PencilAccum<T> operator[](int j) const { return pencils[j]; }

    SlabAccum() {
        // Only have null constructor.
	// This keeps this object lightweight until one actually calls setup()
	pencils = NULL; cells = NULL; pstart = NULL; cpd = 0;
	buffers = NULL; maxthreads = 0;
    }

    // Prevent accidential duplicate ownership of buffers by deleting the copy operators
    SlabAccum(const SlabAccum&) = delete;
    SlabAccum& operator=(const SlabAccum&) = delete;

    void destroy() {
	SlabAccumFree.Start();
	if (cells!=NULL) free(cells); cells = NULL;
	if (pencils!=NULL) free(pencils); pencils = NULL;
	if (pstart!=NULL) free(pstart); pstart = NULL;
	if (buffers!=NULL) {
	    // In tcmalloc, we want to allocate and deallocate on the same thread
	    #pragma omp parallel for schedule(static,1)
	    for (int j=0; j<maxthreads; j++) {
		buffers[omp_get_thread_num()].slabfree();
	    }
	    for (int j=0; j<maxthreads; j++)
	    	assertf(buffers[j].isNULL(),
		    "We failed to free a thread-based malloc.\n");
	    delete[] buffers;
	}
	buffers = NULL;
	cpd = maxthreads = 0;
	SlabAccumFree.Stop();
    }
    ~SlabAccum() {
	destroy();
    }

    inline bool is_setup() const {
	// Returns true is setup has been called; false otherwise.
        return cpd>0;
    }

    /// The routine to initialize the SlabAccum
    /// 
    /// maxsize is not actually a firm maximum; 
    /// rather it is the user direction on the space to initially reserve.
    /// The arrays can grow.

    void setup(int _cpd, int _zwidth, int maxsize) {
	if (pencils==NULL) {
	    cpd = _cpd;
		zwidth = _zwidth;
	    int ret;
	    ret = posix_memalign((void **)&pencils, PAGE_SIZE, sizeof(PencilAccum<T>)*cpd); assert(ret==0);
	    ret = posix_memalign((void **)&cells, PAGE_SIZE, sizeof(CellAccum)*cpd*zwidth); assert(ret==0);
	    ret = posix_memalign((void **)&pstart, PAGE_SIZE, sizeof(uint64)*(cpd+1)); assert(ret==0);
	    for (int j=0; j<cpd; j++) pencils[j].cells = cells+j*zwidth;
	}
	if (buffers==NULL) {
	    // fmt::print("{:d}\n",(sizeof(SlabAccumBuffer<T>)));
	    assert(sizeof(SlabAccumBuffer<T>)==CACHE_LINE_SIZE);
	    // Adjust the SlabAccumBuffer padding, if this fails
	    maxthreads = omp_get_max_threads();
	    buffers = new SlabAccumBuffer<T>[maxthreads];
	    // Initialization of buffers, to reserve memory.

	    // maxsize is the user-supplied estimate per slab
	    // We then divide that across the threads, and then
	    // divide by another factor of 5 in the hopes of not having
	    // too much wasted space.  But the SlabBuffer then round up to the 
	    // nearest factor of 2 of bytes, in the hopes that malloc will get
	    // some re-use of these blocks. 
	    maxsize /= maxthreads;
	    maxsize *= 0.2;
	    maxsize = std::max(maxsize,1024);	// Don't waste time with small stuff
	    
	    // In tcmalloc, we want to allocate and deallocate on the same thread
	    #pragma omp parallel for schedule(static,1)
	    for (int j=0; j<maxthreads; j++) {
		buffers[omp_get_thread_num()].setup(maxsize);
	    }
	    // In the one test I ran, we did always have j==thread_num,
	    // but I chose not to trust this.
	    for (int j=0; j<maxthreads; j++)
	    	assertf(!buffers[j].isNULL(),
		    "We failed to assign a thread-based malloc.\n");
	}
    }

    /// Initialize work for this pencil, pnum=y
    /// Return a pointer to this PencilAccum
    PencilAccum<T> *StartPencil(int pnum) {
	int g = omp_get_thread_num();
	pencils[pnum].buffer = buffers+g;
	pencils[pnum].StartPencil();
	// fmt::print("Starting Pencil {:d} with Thread {:d}\n", pnum, g);
	return pencils+pnum;
    }

    uint64 get_slab_size() const {
        uint64 total = 0;
	for (int j=0; j<maxthreads; j++) total += buffers[j].get_buffer_size();
	return total;
    }

    size_t get_slab_bytes() const {
        return sizeof(T)*get_slab_size();
    }

    void build_pstart() {
	uint64 offset = 0;
	for (int j=0; j<cpd; j++) {
	    pstart[j] = offset;
	    offset += pencils[j]._size;
	}
	pstart[cpd] = offset;
    }

    /// Copy all of the data to a new buffer, making this contiguous 
    /// and in pencil order.  This does not allocate space; use size()
    /// to do that in advance.
	/// Also build the pstart[] array (so this will overwrite it
	/// if build_pstart() was called first, usually harmless).
    void copy_to_ptr(T *destination) {
		build_pstart();

	// This will only do anything if we're already in a parallel region
	#pragma omp taskloop
	for (int j=0; j<cpd; j++) {
		size_t offset = pstart[j];
	    memcpy(destination+offset, pencils[j].data, 
	    	sizeof(T)*pencils[j]._size);
	}
    }
    
    
    /// Same as copy_to_ptr, but accepts a function to do e.g. unit conversion
	template <class U, typename Func>
    void transform_to_ptr(U *destination, Func&& transform) {
		build_pstart();

		// #pragma omp parallel for schedule(dynamic)
        for (int j=0; j<cpd; j++) {
			size_t offset = pstart[j];
            for(int i = 0; i < pencils[j]._size; i++){
                destination[offset + i] = transform(pencils[j].data[i]);
            }
        }
    }


    /// Write all of the data to a file, making this contiguous 
    /// and in pencil order.  
    void dump_to_file(const fs::path &fname) const {
	FILE *fp = fopen(fname.c_str(), "wb");
	for (int j=0; j<cpd; j++) {
	    fwrite(pencils[j].data, sizeof(T), pencils[j]._size, fp);
	}
	fclose(fp);
    }

    /// Write the CellAccum of the data to a file, making this contiguous 
    /// and in pencil order.  
    void dump_cells_to_file(const fs::path &fname) const {
	FILE *fp = fopen(fname.c_str(), "w");
	for (int j=0; j<cpd; j++) {
	    for (int k=0; k<zwidth; k++) {
		fmt::print(fp, "{:d} {:d}:  {:d} {:d}\n",
			j,k, cells[j*zwidth+k].start, cells[j*zwidth+k].size());
	    }
	}
	fclose(fp);
    }

#ifndef TEST		// The test code doesn't provide SB
    void pack(int type, int slab) {
        // This will copy a SlabAccum into an unallocated arena
	// We need to allocate space for, and then fill, SB->(type,slab)
	// This does not delete the SlabAccum

	char *p;
	build_pstart();
	// Precompute the size needed and then allocate it
	uint64 size = 0;
	size+= sizeof(CellAccum)*cpd*zwidth;
	size+= sizeof(uint64)*(cpd+1);
	size+= sizeof(T)*pstart[cpd];
	SB->AllocateSpecificSize(type,slab,size);
	p = SB->GetSlabPtr(type,slab);

	memcpy(p, cells, sizeof(CellAccum)*cpd*zwidth);
	p+= sizeof(CellAccum)*cpd*zwidth;

	memcpy(p, pstart, sizeof(uint64)*(cpd+1));
	p+= sizeof(uint64)*(cpd+1);

	for (int j=0; j<cpd; j++) {
	    memcpy(p, pencils[j].data, sizeof(T)*pencils[j]._size);
	    p+= sizeof(T)*pencils[j]._size;
	}
	return;
    }

    void unpack(int type, int slab) {
	// This unpacks the given Arena into the given SlabAccum, which is 
	// assumed to be uninitialized.
	// This puts the SlabAccum into a static state, after Finish is called.
	// One cannot add any more data.

	// By creating buffers first, we avoid creating one per thread.
	maxthreads = 1;
	buffers = new SlabAccumBuffer<T>[maxthreads];
	setup(P.cpd, node_z_size_with_ghost, 0);
	    // Now cells[], pencils[], and pstart[] exist and pencils[].cells is filled.

	char *p = (char *) SB->GetSlabPtr(type,slab);  // Here's where our data is
	// Load the cells[] array
	memcpy(cells, p, sizeof(CellAccum)*cpd*zwidth);
	p+= sizeof(CellAccum)*cpd*zwidth;

	// Load the pstart[] array and the total size
	memcpy(pstart, p, sizeof(uint64)*(cpd+1));
	p+= sizeof(uint64)*(cpd+1);

	// Allocate the required buffer size
	assertf(pstart[cpd]<2e9, "This arena is too big to be unpacked into one SlabBuffer\n");
	buffers[0].load_from_ptr(p, pstart[cpd]);

	for (int j=0; j<cpd; j++) {
	    pencils[j].buffer = buffers;
	    pencils[j].data = buffers[0].get_pencil_start()+pstart[j];
	    pencils[j]._size = pstart[j+1]-pstart[j];
	}
	return;
    }
#endif

};    // End SlabAccum class

// --------------------------------------------------------------


#ifdef TEST

#define INT int

int main() {
    int cpd = 125;
    int nn = 100;
    SlabAccum<INT> s;
    s.setup(125, 1e5);

    #pragma omp parallel for schedule(static)
    for (int y=0; y<cpd; y++) {
	if (y==0) fmt::print("Running with {:d} threads\n", omp_get_num_threads());
	PencilAccum<INT> *p = s.StartPencil(y);
	for (int z=0; z<cpd; z++) {
	    for (int n=0;n<nn;n++) {
	        int val = n+nn*(z+y*cpd);
		p->append(val);
	    }
	    p->FinishCell();
	}
	p->FinishPencil();
    }

    int tot = 0;
    for (int y=0; y<cpd; y++) tot += s[y]._size;
    fmt::print("Found {:d} objects by buffer or {:d} objects by pencil.  Expected {:d}\n", 
    		(int)(s.get_slab_size()), tot, cpd*cpd*nn);

    for (int y=0; y<cpd; y+=20) {
	for (int z=0; z<cpd; z+=20) {
	    for (int n=0;n<nn;n+=20) {
	        int val = n+nn*(z+y*cpd);
		if (val-s[y][z][n]!=0) {
		    fmt::print("Error: {:d} {:d} {:d} = {:d}  vs  {:d}\n", y, z, n, val, s[y][z][n]);
		    fflush(NULL);
		    exit(1);
		}
	    }
	}
    }

    int *dest;
    dest = (int *) malloc(sizeof(INT)*s.get_slab_size()); 
    s.build_pstart();
    s.copy_to_ptr(dest);
    for (int y=0; y<cpd; y+=20) {
	for (int z=0; z<cpd; z+=20) {
	    for (int n=0;n<nn;n+=20) {
	        int val = n+nn*(z+y*cpd);
		if (val-dest[val]!=0) {
		    fmt::print("Error B: {:d} {:d} {:d} = {:d}  vs  {:d}\n", y, z, n, val, dest[val]);
		    fflush(NULL);
		    exit(1);
		}
	    }
	}
    }
    s.destroy();
    free(dest);
    fmt::print("Test successful!\n");
}

#endif 
