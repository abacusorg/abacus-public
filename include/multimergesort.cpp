/* MultiMergeSort -- Daniel Eisenstein, July 2017

The algorithm sorts cache-sized sub-lists and marks partitions in them.
Then it merges all of those partitions across the sub-lists, using a 
binary heap augmented with FIFO queues.  The queues are useful so that
we do more sequential memory access.

A big advantage in our Insert List problem is that the range of keys
is known in advance and the distribution of keys is expected to be 
approximately uniform.  So we can guess partition values that will
maintain rough load balancing.  In detail, we use uniform to guess at
the partitions, then linearly interpolate to improve that guess.
This should remove slow variations in the population and get the load
balancing pretty close.

The routine is invoked by creating an instance of the class MultiMergeSort
with a templated class.  I.e.,
      MultiMergeSort<MyType> mm;
      mm.sort(in, out, N, maxkey, sublistsize, fifosize);
where in,out are arrays of class MyType, with size N.
The *in array will be re-ordered (but not sorted).
The sorted output is written into *out.  Any input there is overwritten.

The class MyType is the objects to be sorted.  This must have a few properties:
We need to have a key() function, returning an unsigned int between 0 and maxkey.
key = MM_MAX_UINT = 0xffffffff is a special value and should not appear in the input list.
We need the < operation defined, returning based on key().
And there must be a set_max_key() function that sets the key to MM_MAX_UINT.
Formally, set_max_key() without an argument has to ensure that future key() evaluations 
will return MM_MAX_UINT.

sublistsize sets the size in bytes of the sublists that will be
sorted by single threads in the first pass.  We suspect that one
should choose this to be mildly smaller than the per-thread cache.
However, the partitioning will quickly yield cache-sized chunks,
so this may not be too stringent a requirement.  If 0 is entered,
this will default to 1 MB.

fifosize is the length of the fifo queues in the merge sort.  Ideally,
these would fit into a single threads cache, but this is not required.
There will be 1-1.5 queues per sublist, so one can estimate from
the size of the objects whether this will fit in cache.  We suspect
that values of order 8-16 are well chosen.  Smaller values will
defeat the sequential memory access opportunities.


Future possible work:  Could set maxkey (and minkey) adaptively after the first
sorting.  However, this has the disadvantage that the opportunity to do the first
bisection during the first sorting pass is lost.  So this is not interesting for
Abacus.

*/

#include <stdlib.h>
#include <cmath>
#include <algorithm>
#include <unistd.h>
#include <string.h>

#ifdef OMP
    #include <omp.h>
#else
    int get_omp_max_threads() { return 1; }
#endif

#include <sys/time.h>
#include "STimer.cc"
#include "PTimer.cc"
STimer MM_Sorting, MM_Merging, MM_RefinePartition, MM_Indexing;
PTimer MM_SortPart, MM_Bisect, MM_Popping;

#define MM_MAX_UINT 0xffffffff
#define MM_OMP_PLAN static

// For the heap, each node is a MergeFIFO object.

// pop() will pop the next object as well as the peek at the following key.

// Some lists are input lists (left=right=NULL).  When these are exhausted,
// they will pop() out objects with the maximum key value.  These will 
// propagate to other lists, but only at lowest priority.  Therefore,
// they will never pop() off of the head list unless one asks for more objects
// than are in the input list.  Of course, MergeTypes are never allowed to 
// have a valid key equal to this maximum value.  One can also detect the 
// exhaustion of the list by testing whether the following key is maximal.

// By propagating these maximum values, we don't need to explicitly signal 
// that a branch is dry.

template <class MergeType>
class MergeFIFO {
    unsigned int len;    // The max length of the buffer
    unsigned int ptr;    // The next item to read
    unsigned int leftkey, rightkey;   // The next keys in these lists
    MergeType *buf;     // Pointer to [0,len) of space.

    MergeFIFO<MergeType> *left, *right;    // Where we could refill the buffer
    	// If these are NULL, then it means we're at a leaf.

    inline void refill() {
        // Refill the buffer with as many points as available.
	// This routine is never called with an input list (left==NULL)
	for (unsigned int j=0; j<len; j++) {
	    if (leftkey<rightkey) left->pop(buf+j, leftkey);
	        else right->pop(buf+j, rightkey);
	}
        ptr = 0;
	return;
    }
    
  public:

    inline void pop(MergeType *out, unsigned int &nextkey) {
        // Return the next object, refilling the list if we can.
	// Return a peek at the following key
	// For an input list, return maximal keys at the end of the list.
	if (left!=NULL) {
	    // We're in a FIFO queue
	    *out = buf[ptr];
	    if (ptr<len-1) ptr++; else refill(); 
	    nextkey = buf[ptr].key();
	    return;
	} else {
	    // Else we're in an input list
	    if (ptr<len-1) {
		*out=buf[ptr++]; 
		nextkey = buf[ptr].key();
	    } else if (ptr==len-1) {
		*out=buf[ptr++]; 
		nextkey = MM_MAX_UINT;
	    } else {
		out->set_max_key(); nextkey = MM_MAX_UINT;
	    }
	    return;
	}
    }

    inline void pop(MergeType *out) {
	// When popping the head node, we may not need to look at the following key
        unsigned int tmp;
	pop(out, tmp);
	return;
    }

    void Setup(MergeFIFO<MergeType> *_left, MergeFIFO<MergeType> *_right, unsigned int _len) {
	// Given two input lists, set up the buffer.
	// The two input lists must already their buffers full.
	len = _len;
	buf = (MergeType *) malloc(len*sizeof(MergeType));
	left = _left;
	right = _right;
	// Need to initialize leftkey and rightkey, before the
	// first call to refill.  Note that the branch might be empty!
	if (left->len==0) leftkey = MM_MAX_UINT; else leftkey = left->buf[0].key();
	if (right->len==0) rightkey = MM_MAX_UINT; else rightkey = right->buf[0].key();
	// Now fill up the list
	refill();    
    }

    void Setup(MergeType *data, unsigned int _len) {
        // Set up a leaf that just points into the original data set
	buf = data;
	len = _len;
	ptr = 0;
	left = right = NULL;  // This signals that it is an input list, with no malloc.
	leftkey = rightkey = 0;   // These aren't used in this case, since refill() never called.
    }

    // The NULL constructor is not useful; must overwrite with Setup()
    MergeFIFO<MergeType>() { left = right = NULL; }

    ~MergeFIFO<MergeType>() { 
	// Free the buffer space.
	// Note that we are not freeing the input branches!
    	if (left!=NULL) free(buf); return; 
    }

};   // End MergeFIFO

// ====================================================


template <class MergeType> 
class MultiMergeSort {

private:

unsigned int bisection_search(MergeType *in, unsigned int size, unsigned int minvalue) {
    // Return the index of in[0..size) of the lowest element that is >= minvalue.
    // Return size if all of the elements are below minvalue
    unsigned int low = 0, high = size-1, mid;
    if (in[low].key()>=minvalue) return low;
    if (in[high].key()<minvalue) return high+1;
    // Otherwise, we search.
    while (high-low>1) {
        mid = (low+high)/2;
	if (in[mid].key()>=minvalue) high = mid; else low = mid;
    }
    return high;
}

void accumulate_partition_boundaries(unsigned int *index, 
		unsigned int partsize[], unsigned int pstart[]) {
    // Just sum over the sub-list to find the partition lengths and 
    // cumulation of them.
    for (int k=0; k<Nparts; k++) partsize[k] = 0;
    #pragma omp parallel for schedule(MM_OMP_PLAN)
    for (int k=0; k<Nparts; k++) 
	for (int j=0; j<Nsublist; j++) 
	    partsize[k] += index[j*Nparts+k+1]-index[j*Nparts+k];
    pstart[0]=0; 
    for (int k=0; k<Nparts; k++) pstart[k+1] = pstart[k]+partsize[k];
    return;
}

unsigned int interpolate(unsigned int minvalue[], unsigned int pstart[], unsigned int goal) {
    // Goal is the number of particles we're trying to find.
    // Look in the (minvalue,pstart) list to estimate where the goal would be.
    for (int i=0; i<Nparts; i++) {
	if (pstart[i+1]>=goal) {
	    // The goal is between pstart[i] and pstart[i+1].
	    // Also guaranteed that pstart[i+1]>pstart[i], else we would have already stopped
	    return minvalue[i]+int((float)(minvalue[i+1]-minvalue[i])*(goal-pstart[i])/(pstart[i+1]-pstart[i]));
	}
    }
}

public:

// We reveal these just in case someone wants to see what was chosen.
unsigned int Nsublist;    // Number of sublists for first-phase sorting
unsigned int Nparts;      // Number of partitions for second-phase merging


void mmsort(MergeType *a, MergeType *out, unsigned int N, unsigned int maxkey, unsigned int sublistsize, unsigned int FIFO_SIZE) {
    // Sort the list a[0..N-1].
    // Object MergeType should have a method key() to return an unsigned int to sort.
    // maxkey should be larger than all of the keys, but only just barely.
    // sublistsize is the L2 cache per core, in bytes, used to scope the work.

    // It turns out to be easier to write this as copying into a second copy
    // of the list. 

    // This routine is fairly liberal in allocating auxillary space for various
    // indices and queues.  Probably at the level of 10 MB.

    assert(maxkey!=MM_MAX_UINT);    // Can't use the MM_MAX_UINT value!
    if (FIFO_SIZE==0) FIFO_SIZE = 16;
    assert(FIFO_SIZE<257);   // Don't go crazy
    if (sublistsize<1e4) sublistsize = 1e6;
	// The user has probably entered MB or objects instead of
	// bytes; give them a reasonable modern number.

    if (N==0) return;    // No work to be done.
    if (N==1) {*out = *a; return;}  // That was quick too!

    // First step is to divide into N sub-lists, each to fit in sublistsize.
    // We will sort each one.

    MM_Indexing.Start();
    sublistsize /= sizeof(MergeType);   // Convert this to list units
    Nsublist = ceil((float)N/sublistsize);
    Nparts = omp_get_max_threads();

    // printf("For %d particles, planning %d sublists of %d particles. %d partitions.\n", 
    //      N, Nsublist, sublistsize, Nparts);

    // We aim to divide each sublist into Nparts parts, with partition values 
    // chosen to approximately balance the number of points in each partition when
    // summed over all sublists.
    // We have to find and store the starting point as index[SubList][Nparts]

    unsigned int minvalue[Nparts+1];  // The minimum key value for each partition.
    for (int k=0; k<Nparts; k++) minvalue[k] = ((float)maxkey/Nparts)*k;
    	// Here we assume that the keys are approximately homogenouesly distributed 
	// between 0 and maxkey.
    minvalue[Nparts] = maxkey;    // Just for interpolation

    unsigned int *index;
    assert(posix_memalign((void **)&index, 64, sizeof(unsigned int)*(Nsublist*Nparts+1))==0);
    index[Nsublist*Nparts] = N;
    for (int j=0; j<Nsublist; j++) index[j*Nparts] = j*sublistsize;
    // We will have to fill in the other values as we go.
    MM_Indexing.Stop();

    MM_Sorting.Start();
    #pragma omp parallel for schedule(MM_OMP_PLAN)
    for (int j=0; j<Nsublist; j++) {
        // Now we loop over each sublist.
	unsigned int start = index[j*Nparts];
	unsigned int size = index[(j+1)*Nparts]-start;

	// Sort, moving the result into the work list
	MM_SortPart.Start();
	std::sort(a+start, a+start+size);
	MM_SortPart.Stop();

	assert(a[start+size-1].key()<=maxkey);
	    // Just checking for input that would break the merge step

	// While this sublist may still be in cache, search for the breakpoints
	MM_Bisect.Start();
	for (int k=1; k<Nparts; k++) 
	    index[j*Nparts+k] = start + bisection_search(a+start, size, minvalue[k]);
	MM_Bisect.Stop();
    }
    MM_Sorting.Stop();

    // We now have Nsublist sorted sublists, each separated into Npart partitions,
    // with the bounding indices all saved.  

    MM_Indexing.Start();
    unsigned int partsize[Nparts];
	// This is the size of the partition, summed over sublists.
    unsigned int pstart[Nparts+1];
	// And here's the cumulation of these to find the starting point for each partition.
    accumulate_partition_boundaries(index, partsize, pstart);
    MM_Indexing.Stop();



    // Optional: if the partsize[] list reveals unbalanced sizes, then the next step will
    // not be well load-balanced.  One might prefer to set different minvalue[] targets and 
    // re-search for the index[] values.  

    // for (int k=0; k<Nparts; k++) 
        // printf("Partition %2d at key %7d has %7d particles\n", k, minvalue[k], partsize[k]);

    MM_RefinePartition.Start();

    // For each partition, we linearly interpolate the known points to try to 
    // estimate better minvalues.
    unsigned int newminvalue[Nparts];  
    newminvalue[0] = 0;
    #pragma omp parallel for schedule(MM_OMP_PLAN)
    for (int k=1; k<Nparts; k++) 
        newminvalue[k] = interpolate(minvalue, pstart, ((float)N/Nparts)*k);

    // Now we have our new partition points.  Repeat the bisection.
    #pragma omp parallel for schedule(MM_OMP_PLAN)
    for (int j=0; j<Nsublist; j++) {
	unsigned int start = index[j*Nparts];
	unsigned int size = index[(j+1)*Nparts]-start;
	for (int k=1; k<Nparts; k++) 
	    index[j*Nparts+k] = start + bisection_search(a+start, size, newminvalue[k]);
    }
    accumulate_partition_boundaries(index, partsize, pstart);

    MM_RefinePartition.Stop();

    // printf("Rebalancing the partition\");
    // for (int k=0; k<Nparts; k++) 
        // printf("Partition %02d at key %7d has %7d particles\n", k, newminvalue[k], partsize[k]);



    // Now we want to merge the Nsublist lists for each of the partitions, putting the result
    // back into the a[] list.
    
    MM_Merging.Start();
    #pragma omp parallel for schedule(MM_OMP_PLAN)
    for (int k=0; k<Nparts; k++) {
        MergeType *outpart = out+pstart[k]; 

	MergeFIFO<MergeType> *m = new MergeFIFO<MergeType>[3*Nsublist];;

	int j;
	for (j=0; j<Nsublist; j++)
	    m[j].Setup(a+index[j*Nparts+k], index[j*Nparts+k+1]-index[j*Nparts+k]);
	    // This loads the leaf-level input lists 
	    // with direct pointers to the partitions.
	    
	// Now we continue with a binary aggregation
	int n=Nsublist;
	while (n>1) {
	    int jj = j-n;    // Starting point of the last cycle
	    for (int i=0; i<n-1; i+=2, j++)
		m[j].Setup(m+jj+i, m+jj+i+1, FIFO_SIZE);
		// We define a new queue that will merge two of the previous set

		// If n was odd, then we will only have added n/2 (round down) merges.
		// The last of the last round is unmatched.  
		// So we carry forward (n/2) (round up) to the next round
	    n=ceil((float)n/2);
	}
	MergeFIFO<MergeType> *head = m+j-1;

	// Now we zip the lists together!
	MM_Popping.Start();
	for (unsigned int p=0; p<partsize[k]; p++) head->pop(outpart+p);
	MM_Popping.Stop();

	// Clean up
	delete[] m;
    }
    MM_Merging.Stop();

    // Clean up
    free(index);
    return;
}

};

#ifdef UNITTEST

// Here's a place-holder for the MergeType object:
// 
// We need to have a key() and set_key() function.
// And we need the < operation defined.
// key = MM_MAX_UINT is a special value, not users.
// set_key() without an argument has to ensure that key() will return MM_MAX_UINT

class MyMergeType {
  private:
    unsigned int k;    // A valid user key must be less than MM_MAX_UINT.
    float dat[11];

  public:
    inline unsigned int key() { return k; }
    inline void set_max_key() { k = MM_MAX_UINT; }
    bool operator< (const MyMergeType& b) const { return (k<b.k); }

    inline void set_key(unsigned int _k) { k = _k; }
};

// ================================

int main() {
    unsigned int N = 3e6;
    int Iter = 10;
    MyMergeType *il = new MyMergeType[N];
    MyMergeType *out = new MyMergeType[N];
    STimer time;
    MultiMergeSort<MyMergeType> mm;

    for (int i=0;i<Iter;i++) {
	#pragma omp parallel for schedule(static)
	for (int j=0;j<N;j++) il[j].set_key(floor(pow(drand48(),0.9)*1485*1485));
	    // Setting this up to be a little inhomogeneous, to see if the rebalancing is working.
	printf("."); fflush(NULL);
	time.Start();
	mm.mmsort(il, out, N, 1485*1485, 2e6, 16);
	time.Stop();
    }

    printf("\nTime to sort %d objects with %d threads, %d times: %f sec, %f Mobj/sec\n", N, omp_get_max_threads(), Iter, time.Elapsed(), Iter*N/1e6/time.Elapsed());
    printf("Wall-clock Time to Sort & Partition Sublists: %f\n", MM_Sorting.Elapsed());
    printf("    Sort Sublist (P): %f\n", MM_SortPart.Elapsed()/omp_get_max_threads());
    printf("    Bisect (P):       %f\n", MM_Bisect.Elapsed()/omp_get_max_threads());

    printf("Indexing:             %f\n", MM_Indexing.Elapsed());
    printf("Refine Partition:     %f\n", MM_RefinePartition.Elapsed());

    printf("Wall-clock Time to Merge the Partitions: %f\n", MM_Merging.Elapsed());
    printf("    Popping (P):      %f\n", MM_Popping.Elapsed()/omp_get_max_threads());

    for (int j=0;j<N-1;j++) if(out[j].key()>out[j+1].key())
	{ printf("Failed key(%d) = %d, key(%d) = %d\n", j, out[j].key(), j+1, out[j+1].key());
	return 1;
	}
    printf("List is verified to be sorted!\n");
    return 0;
}

#endif
