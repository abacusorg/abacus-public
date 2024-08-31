/** \file The GroupLinkList is a pan-slab list of pairs of CellGroups,
where the pairs have at least one particle within the FOF distance.
Having found the CellGroups, we then form them into GlobalGroups by
traversing the graph defined by these GroupLinks.

*/

#ifdef TEST
#include "test_driver_header.cpp"
#endif

#ifndef STANDALONE_FOF
#include "tbb/parallel_sort.h"
// #include "parallel_stable_sort.h"
#endif

/** A LinkID is a unique reference to a CellGroup, using cell i,j,k
and then the group number.

We pack 12 bits of cell location (so CPD<4096!)
and then 20 bits of group number, to make a sortable object.
Use the sign bit to mark items for deletion

If the id is -1, then this link is unusable and marked for deletion.

2D: k is intended to be a relative value, rather than global.
This allows one use fewer bits for k.
*/

class LinkID {
  public:
    long long int id;

    inline LinkID() { return; };
    inline LinkID(int i, int j, int k, int n) {
        id = (((long long int)i*4096+(long long int)j)*4096+(long long int)k)*1048576+n;
    }
    inline integer3 localcell() const {
        integer3 c;
        c.x = (id&(uint64)0xfff00000000000)>>44;
        c.y = (id&(uint64)0xfff00000000)>>32;
        c.z = (id&(uint64)0xfff00000)>>20;
        return c;
    }
    inline int cellgroup() {
        return id&0xfffff;
    }
    inline void mark_for_deletion() {
        id = -1; return;
    }
    inline int slab() {
        return (id&(long long int)0xfff00000000000)>>44;
    }
    // For sorting the LinkIDs inside a completed GlobalGroup
    bool operator< (const LinkID& c) const { return (id<c.id); }
};

/** A GroupLink is simply a pair of LinkIDs.  
The class supplies a comparison operation for sorting.
*/

class GroupLink {
  public:
    LinkID a,b;
    inline GroupLink() { }
    inline GroupLink(LinkID _a, LinkID _b) { a = _a; b = _b; }
    bool operator< (const GroupLink& c) const { return (a.id<c.a.id); }
        // Sort on the first element

    inline void mark_for_deletion() { a.mark_for_deletion(); return; }
    inline bool is_marked_for_deletion() { return (a.id==-1);}
};

inline bool is_marked_for_deletion(GroupLink *gl, int val) {
    // This is the binary decision for whether a link will be swept
    // to the end of the list for deletion.
    // val is not used in this routine, but required by partitioning code.
    return gl->is_marked_for_deletion();
}



/** This is the list of GroupLinks.  

There is only one for the whole simulation, since links cross cells
and slabs.

When a link is added, it is actually added twice with the symmetric 
reflection.  During GlobalGroup finding, the pair is traversed in both
directions and each link is marked for deletion after it is considered.
A link that points back to an already used CellGroup is a dead-end,
as that group is already in the set for this GlobalGroup.

This class uses the MultiAppendList class, which allows multiple
threads to be appending into one buffer and then closes up the gaps
when done, using a minimum of re-copying.
*/

class GroupLinkList : public grid, public MultiAppendList<GroupLink> {
public:
    STimer GroupPartition, GroupSort;
    GroupLinkList(int cpd, uint64 maxsize) : 
    		grid(cpd), MultiAppendList<GroupLink>(maxsize, PAGE_SIZE/sizeof(GroupLink))  { 
    }
    ~GroupLinkList(void) { 
    }
    
    /// This pushes one pair onto the list,  reciprocally
    inline void DoublePush(LinkID _a, LinkID _b) {
	GroupLink gl(_a,_b);
	Push(gl);
	gl.a = _b;
	gl.b = _a;
	Push(gl);
	return;
    }

    /// This brings all of the links that are marked for deletion
    /// to the end of the list, and then shrinks the list to discard them. 
    void PartitionAndDiscard() {
	GroupPartition.Start();
	uint64 mid = ParallelPartition(list, length, is_marked_for_deletion, (int)0);
		// [0..mid) are to be retained.
		// [mid..length) are to be deleted.
	GroupPartition.Stop();
	ShrinkMAL(mid);
	return;
    }

    // This sorts the whole list of GroupLinks.  
    void Sort() {
	GroupSort.Start();
	#ifndef STANDALONE_FOF
    ips4o::parallel::sort(list, list+length);
	#else
	std::sort(list, list+length);
	#endif
	GroupSort.Stop();
	// Our other stuff doesn't work because these keys are uint64.
	// GroupLink *hlnew;
	// int ret = posix_memalign((void **)&glnew, CACHE_LINE_SIZE, sizeof(GroupLink)*(length));
	// assertf(ret==0, "GroupLinkList Sort() failed to allocate memory.\n");
	// MultiMergeSort<GroupLink> mm;
	// mm.mmsort(&(list[mid]), glnew, length, cpd*cpd, 2e6, 16);
    }

    void AsciiPrint() {
        for (uint64 j=0; j<length; j++) {
	    integer3 c1 = list[j].a.localcell();
	    integer3 c2 = list[j].b.localcell();
	    fmt::print("{:d} {:d} {:d} {:d} to {:d} {:d} {:d} {:d}\n",
	        c1.x, c1.y, c1.z, list[j].a.cellgroup(),
	        c2.x, c2.y, c2.z, list[j].b.cellgroup());
	}
	fmt::print("Dump of GLL length {:d}\n", length);
    }

    /** This searches the sorted GroupLinkList to 
    return the element in the list that is the first one greater or
    equal to the k=0,n=0 element for the cells i=slab,j.

    We assume the inputs are both in [0,cpd) range, already wrapped.
    LinkIDs are not unique nor contiguous.
    */
    GroupLink *Search(int slab, int j) {
	GroupLink ref;
	ref.a = LinkID(slab, j, 0, 0);
	if (length==0) return list;
	uint64 low = 0, high = length-1, mid;
	// fmt::print("Seeking {:d} {:d} = {:d}\n", slab, j, ref.a.id);
	// if (ref < list[low] || ref.a == list[low].a ) return list+low;
	if (!(list[low]<ref)) return list+low;    // Found at beginning
	if (list[high] < ref) return list+high+1;   // The end of the list
	while (high-low>1) {
	    mid = (low+high)/2;
	    if (list[mid]<ref) low = mid; else high = mid; 
	}
	// We end when the two are only 1 apart
	assertf(list[low]<ref && !(list[high]<ref),
		"Error in GLL index search; list not sorted.");
	return list+high;
    }
};   // End GroupLinkList



#ifdef TEST

int main() {
    GroupLinkList *GLL;
    GLL = new GroupLinkList(13,1e7);

    for (int iter = 0 ; iter<2; iter++) {
	fmt::print("Size of GLL: {:d}\n", (int)(GLL->length));
	#pragma omp parallel for schedule(static)
	for (int j=0; j<1000; j++) {
	    for (int k=0; k<1000; k++) {
		LinkID a(1,2,3,4), b(5,6,7,8);
		GLL->DoublePush(a,b);
	    }
	}
	fmt::print("Size of GLL: {:d}\n", (int)GLL->length);
	GLL->CollectGaps();
	fmt::print("Size of GLL: {:d}\n", (int)GLL->length);

	for (int j=0; j<1e6;j+=100) GLL->list[j].mark_for_deletion();
	GLL->PartitionAndDiscard();
	fmt::print("Size of GLL: {:d}\n", (int)GLL->length);
    }

    return 0;
}

#endif
