// multiappendlist.cpp

/* The class MultiAppendList is a multi-threaded append-in-place code.
When the process is completed, the objects will be contiguous, but 
there is no expectation that the objects will remain in order (that 
would be ill-defined in a multi-threaded context, but is not even 
guaranteed within a thread).

One allocates a big chunk of memory at the beginning.  Then the code
gives each thread a small segment (gap_size) of that 
list in which to append.  New segments are created invisibly when the
thread runs out of space.  One calls CollectGaps() when done; this 
copies the filled fragments of segments into the empty parts of segments
and creates a contiguous list.

The basic operation to append to the list is Push().
We also provide ShrinkMAL() and GrowMAL() in case one needs to adjust the
amount of space used by fiat.  Note that these require unsigned int's 
as input; there's not a check that these are >=0.  
ShrinkMAL() is the way to remove particles from the end of the list.
GrowMAL() would be an external way to add extra particles in bulk.

There is no provision to grow the originally allocated buffer!  But the
code will detect an overflow and throw an assertf().

The gap_size is a performance tuning parameter. Different MALs may prefer
different sizes. The Insert List is a particularly important MAL, and
its gap size can be set as P.InsertListGapElems.

We want a gap size big enough that we don't have threads
waiting for the mutex to clear.  But small enough that 
the gap collection at the end doesn't have big copies.
Gap sizes that are a multiple of the NUMA page size will
encourage locality, although our 40-byte ilstruct doesn't
divide power-of-two pages sizes evenly. One might try
the rather-heavy 5*PAGE_SIZE, though.

*/

#ifdef MALTEST
#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include "stdlog.cc"
typedef unsigned long long int uint64;
int omp_get_max_threads() { return 1;}
int omp_get_num_threads() { return 1;}
int omp_get_thread_num() { return 0;}
#endif

template <class T>
class MultiAppendList;

template <class T>
class MALgap {
    // The gap will be defined as [start,end)
    // As we fill in from the front, start will be incremented
    // nextid() will return a valid index number, where we can place a new 
    // particle.
    // size() and gapsize() return the used and unused sizes of this block.
  private:
    void make_next_gap(MultiAppendList<T> *mal);
  public:

    uint64 start;
    uint64 end;
    uint64 next;
    uint8_t pad[CACHE_LINE_SIZE-24];   // We want fill a cache line to avoid contention

    // When sorting, we sort by end, which is unique across gaps.
    // In principle, two gaps could have the same start
    bool operator< (const MALgap& b) const { return (end<b.end); }
  
    void reset_gap(uint64 length) { start = end = next = length; }
    MALgap() { reset_gap(0); }
    ~MALgap() {} ;


    inline uint64 nextid(MultiAppendList<T> *mal) {
        // Return the next index in line
        if (next==end) make_next_gap(mal);
        uint64 val = next;   // The index we will return
        next++;
        return val;
    }
    inline uint64 gapsize() { return end-next; }
    inline uint64 size() { return next-start; }
};  // End MALgap class


template <class T>
class MultiAppendList {
private:
    MALgap<T> *MALgaps;
    int Ngaps;

    void reset_gaps() {
        for (int j=0;j<Ngaps;j++) MALgaps[j].reset_gap(length);
    }

public: 
    T *list;
    uint64 length;
    uint64 longest;
    uint64 maxlist;
    const uint64 gap_size;

    MultiAppendList(uint64 maxlistsize, uint64 _gap_size)
            : gap_size(_gap_size) { 
        length = 0; 
        longest = 0;
        Ngaps = omp_get_max_threads();
        // we may try to grow the list by an extra block per thread
        maxlist = maxlistsize + gap_size*Ngaps;
        int ret = posix_memalign((void **) &list, PAGE_SIZE, sizeof(T) * maxlist);
        assertf(ret==0,"Failed to allocate MultiAppendList\n");
        MALgaps = new MALgap<T>[Ngaps];
        return;
    }
    ~MultiAppendList(void) { 
        free(list);
        delete[] MALgaps;
    }
    
    // Push to the next allowed position in this thread.
    inline void Push(T obj) {
        list[MALgaps[omp_get_thread_num()].nextid(this)] = obj;
        return;
    }

    inline void GrowMALGap(uint64 newlength) { 
        // One probably shouldn't use this version outside of the
        // MALgaps.
        assertf(newlength>=length, 
            "Illegal growing of MultiAppendList (length = {:d}, newlength = {:d})\n", length, newlength);
        assertf(newlength < maxlist, 
            "Illegal resizing of MultiAppendList (maxlist = {:d}, newlength = {:d}\n", maxlist, newlength);
        length = newlength; 
    }

    inline void ShrinkMAL(uint64 newlength) { 
        // Call this after CollectGaps()
        assertf(newlength<=length, 
            "Illegal shrinking of MultiAppendList (length = {:d}, newlength = {:d})\n", length, newlength);
        length = newlength; 
        reset_gaps();
    }
    
    inline void GrowMAL(uint64 newlength) { 
        // This is the version to use if one wants to add length
        // in an ad hoc manner.
        // Call this after CollectGaps()
        GrowMALGap(newlength);
        if (length>longest) longest = length;
        reset_gaps();
    }

    void CollectGaps() {
        // Given a list of gaps, we want to fill them in and reduce the 
        // size of the insert list.  This routine will change both the
        // order and the contents of the input vector; essentially,
        // the vector is inapplicable after use.

        // This routine is entirely serial.
        // In principle, it could be parallelized by recording a set of
        // planned moves, then splitting them up to balance the work 
        // across threads, then executing.  See parallel partitioning 
        // for an example.
        // However, the expected amount of work here seems modest,
        // of order NThreads * gap_size / 4 objects to copy.

        // std::sort( MALgaps, MALgaps+Ngaps );
        ips4o::sort(MALgaps, MALgaps+Ngaps);
            // Order the gaps by their end index
        int low = 0;          // low[next:end) is unfilled space
        int high = Ngaps-1;   // high[start:next) is filled space
        // We can assume that the highest gap is the end of the MAL.
        // But that doesn't always stay true, unless we do something about it.
        // Extend the filled part of each gap to the last unfilled point of the 
        // one below it.  
        for (int j=1;j<Ngaps;j++) MALgaps[j].start = MALgaps[j-1].end;

        while (high>low) {
            // When high==low, then there are no more gaps to fill
            if (MALgaps[low].gapsize() <= MALgaps[high].size()) {
                // We have enough particles to fill the gap
                uint64 copysize = MALgaps[low].gapsize();   // Number to move
                assert(MALgaps[high].next >= copysize);
                memcpy(list+MALgaps[low].next, 
                           list+MALgaps[high].next-copysize,
                       sizeof(T)*copysize);
                // And this means we're done with this low gap.
                MALgaps[high].next -= copysize;
                low++;
            } else {
                // The gap is more than the particles we have.
                uint64 copysize = MALgaps[high].size();   // Number to move
                assert(MALgaps[high].next >= copysize);
                memcpy(list+MALgaps[low].next, 
                           list+MALgaps[high].next-copysize,
                       sizeof(T)*copysize);
                // And this means we're done with the high set
                MALgaps[low].next += copysize;
                high--;
            }
        }
        ShrinkMAL(MALgaps[low].next);
	if (length>longest) longest = length;
        return;
    }

};

template <class T>
void MALgap<T>::make_next_gap(MultiAppendList<T> *MAL) {
    // Adjustments to the MAL list length can only happen one at a time

    #pragma omp atomic capture
    end = MAL->length += MAL->gap_size;

    next = start = end - MAL->gap_size;

    assertf(end <= MAL->maxlist,
        "Illegal resizing of MultiAppendList (maxlist = {:d}, newlength = {:d}\n", MAL->maxlist, end);
}
