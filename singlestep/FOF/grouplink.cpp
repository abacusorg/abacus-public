#ifdef TEST
#include "test_driver_header.cpp"
#endif

#ifndef STANDALONE_FOF
#include "tbb/parallel_sort.h"
// #include "parallel_stable_sort.h"
#endif

class LinkID {
    // We pack 12 bytes of cell location (so CPD<4096!)
    // and then 20 bytes of group number, to make a sortable object.
    // Use the sign bit to mark items for deletion
  public:
    long long int id;

    inline LinkID() { return; };
    inline LinkID(int i, int j, int k, int n) {
        id = (((long long int)i*4096+(long long int)j)*4096+(long long int)k)*1048576+n;
    }
    inline integer3 cell() {
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
};

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

/* Not sure what's going on with this construction; vestigial from a parallel sort
struct GroupLinkSortOperator {
    inline bool operator() (const  GroupLink &pi, const GroupLink &pj ) const {
	return pi<pj;
    }
};
*/

inline bool is_marked_for_deletion(GroupLink *gl, int val) {
    // This is the binary decision for whether a link will be swept
    // to the end of the list for deletion.
    // val is not used in this routine, but required by partitioning code.
    return gl->is_marked_for_deletion();
}


class GroupLinkList : public grid, public MultiAppendList<GroupLink> {
public:
    STimer GroupPartition, GroupSort;
    GroupLinkList(int cpd, uint64 maxsize) : 
    		grid(cpd), MultiAppendList<GroupLink>(maxsize)  { 
    }
    ~GroupLinkList(void) { 
    }
    
    inline void DoublePush(LinkID _a, LinkID _b) {
	// We push reciprocally
	GroupLink gl(_a,_b);
	Push(gl);
	gl.a = _b;
	gl.b = _a;
	Push(gl);
	return;
    }

    void PartitionAndDiscard() {
	GroupPartition.Start();
	uint64 mid = ParallelPartition(list, length, (int)0, is_marked_for_deletion);
		// [0..mid) are to be retained.
		// [mid..length) are to be deleted.
	GroupPartition.Stop();
	ShrinkMAL(mid);
	return;
    }

    void Sort() {
	// We sort the whole list of GroupLinks.
	GroupSort.Start();
	// Put in a parallel sort
	#ifndef STANDALONE_FOF
	tbb::parallel_sort(list, list+length);
    // pss::parallel_stable_sort( list, list+length, GroupLinkSortOperator() );
	#else
	std::sort(list, list+length);
	#endif
	GroupSort.Stop();
	// Our other stuff doesn't work because these keys are uint64.
	// GroupLink *hlnew;
	// int ret = posix_memalign((void **)&glnew, 64, sizeof(GroupLink)*(length));
	// assertf(ret==0, "GroupLinkList Sort() failed to allocate memory.\n");
	// MultiMergeSort<GroupLink> mm;
	// mm.mmsort(&(list[mid]), glnew, length, cpd*cpd, 2e6, 16);
    }

    void AsciiPrint() {
        for (int j=0; j<length; j++) {
	    integer3 c1 = list[j].a.cell();
	    integer3 c2 = list[j].b.cell();
	    printf("%d %d %d %d to %d %d %d %d\n",
	        c1.x, c1.y, c1.z, list[j].a.cellgroup(),
	        c2.x, c2.y, c2.z, list[j].b.cellgroup());
	}
	printf("Dump of GLL length %lld\n", length);
    }

    GroupLink *Search(int slab, int j) {
        // We assume these are both in [0,cpd) range, already wrapped.
	// Return the element in the list that is the first one greater or
	// equal to the k=0,n=0 element.
	// LinkIDs are not unique nor contiguous
	// TODO: What if all elements are less than the goal value?
	// TODO: What if the list is empty?
	GroupLink ref;
	ref.a = LinkID(slab, j, 0, 0);
	if (length==0) return list;
	uint64 low = 0, high = length-1, mid;
	// printf("Seeking %d %d = %lld\n", slab, j, ref.a.id);
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
	printf("Size of GLL: %d\n", (int)(GLL->length));
	#pragma omp parallel for schedule(static)
	for (int j=0; j<1000; j++) {
	    for (int k=0; k<1000; k++) {
		LinkID a(1,2,3,4), b(5,6,7,8);
		GLL->DoublePush(a,b);
	    }
	}
	printf("Size of GLL: %d\n", (int)GLL->length);
	GLL->CollectGaps();
	printf("Size of GLL: %d\n", (int)GLL->length);

	for (int j=0; j<1e6;j+=100) GLL->list[j].mark_for_deletion();
	GLL->PartitionAndDiscard();
	printf("Size of GLL: %d\n", (int)GLL->length);
    }

    return 0;
}

#endif
