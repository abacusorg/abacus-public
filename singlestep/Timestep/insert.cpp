/* insert.cpp
 *
 *Class for the insert list.
 *
 *Particles are added to the end of the insert list, along with their
 *new cell number.  When requested, a partioning is done to place particles 
 *of a given slab at the end of the list; this slab is then sorted into
 *cell order.  We will then copy those particles out and delete them from
 *the insert list.
 *
 *Importantly, particle positions are adjusted when the particle is 
 *added to the insert list, so that the particle has a legal position
 *for its *new* cell.
 */

/*
 * FIXME: Theres a lot of vestigial and commented out code that should be cleaned up.
 *  We can rely on source control if it turns out later to be needed
 */

//for intel's fast parallel sort
// #include "tbb/parallel_sort.h"
// #include "parallel_stable_sort.h"
#include "parallel.partition.cpp"

#include "multimergesort.cpp"

#include <mutex>



class ilstruct {
  public:
    unsigned int k;   // The key based on y*cpd+z
    integer3 xyz;        // The new cell, wrapped
    posstruct pos;
    velstruct vel;
    auxstruct aux;

    // Some methods require for sorting
    inline unsigned int key() { return k; }
    inline void set_max_key() { k = MM_MAX_UINT; }
    bool operator< (const ilstruct& b) const { return (k<b.k); }
};   // end ilstruct


class ilgap {
    // The gap will be defined as [start,end)
    // As we fill in from the front, start will be incremented
    // nextid() will return a valid index number, where we can place a new 
    // particle.
    // size() and gapsize() return the used and unused sizes of this block.
  private:
    void make_next_gap();
  public:

    uint64 start;
    uint64 end;
    uint64 next;
    uint64 pad[5];   // We want fill a 64-byte cache line to avoid contention
    #define ILGAP_SIZE 512
    // We want something big enough that we don't have threads
    // waiting for the mutex to clear.  But small enough that 
    // the gap collection at the end doesn't have big copies.

    // When sorting, we sort by end, which is unique across gaps.
    // In principle, two gaps could have the same start
    bool operator< (const ilgap& b) const { return (end<b.end); }
  
    ilgap() { make_next_gap(); return; }
    ~ilgap() {} ;

    uint64 nextid() {
        // Return the next index in line
        uint64 val = next;   // The index we will return
        next++;
        if (next==end) make_next_gap();
            // Done with this gap; need a new allocation
        return val;
    }
    uint64 gapsize() { return end-next; }
    uint64 size() { return next-start; }
};  // End ilgap class


// The following routine is vestigial FIXME:Delete it.
struct GlobalSortOperator {
    inline bool operator() (const  ilstruct &pi, const ilstruct &pj ) const {
        // We must sort on pi.xyz.y*cpd+pi.xyz.z < pj.xyz.y*cpd+pj.xyz.z
        // Warning: This will get more complicated in the parallel code with 
        // the periodic wrap.
        if(pi.xyz.y - pj.xyz.y == 0)
            return pi.xyz.z < pj.xyz.z;
        return pi.xyz.y < pj.xyz.y;
        //return  pi.xyz.y*P.cpd + pi.xyz.z < pj.xyz.y*P.cpd + pj.xyz.z; 
    }
};

void ConfirmSorting(ilstruct *il, uint64 len) {
    for (uint64 j=0; j+1<len; j++) 
    assertf(il[j].key()<=il[j+1].key(), 
            "Insert list sorting failed: il[%d]=%d,  il[%d]=%d\n",
        j, il[j].key(), j+1, il[j+1].key());
    return;
}


// Partition function; returns true if `particle` belongs in `slab`
inline bool is_in_slab(ilstruct *particle, int slab){
    return particle->xyz.x == slab;
}


// WARNING: The insert list is written to accept uint64 sizes in its list.
// So be careful when using the indices you receive from it! 


class InsertList : public grid {
public: 
    ilstruct *il;
    uint64 length;
    uint64 maxil;

    uint64 n_sorted;

    InsertList(int cpd, uint64 maxilsize) : grid(cpd)  { 
        length = 0; 
        // we may try to grow the list by an extra block per thread
        maxil = maxilsize + ILGAP_SIZE*omp_get_max_threads();
        int ret = posix_memalign((void **) &il, 4096, sizeof(ilstruct) * maxilsize);
        assertf(ret==0,"Failed to allocate Insert List\n");
        n_sorted = 0;
    }
    ~InsertList(void) { 
        free(il);
    }

    // Push to the end of the list and grow
    inline void Push(posstruct  *pos, velstruct *vel, auxstruct *aux, integer3 xyz) {
        assertf(length<maxil+1, "Push overflows Insert List length\n");
        il[length].pos = *pos;
        il[length].vel = *vel;
        il[length].aux = *aux;
        il[length].xyz = xyz;
        il[length].k = xyz.y*cpd + xyz.z;
        length++;
    }
    
    // Push to a given location; don't grow
    inline void PushAt(posstruct  *pos, velstruct *vel, auxstruct *aux, integer3 xyz, uint64 i) {
        assertf(i<length, "Push overflows Insert List length\n");
        il[i].pos = *pos;
        il[i].vel = *vel;
        il[i].aux = *aux;
        il[i].xyz = xyz;
        il[i].k = xyz.y*cpd + xyz.z;
    }

    inline integer3 WrapToNewCell(posstruct *pos, int oldx, int oldy, int oldz) {
        // We have been given a particle that has drifted out of its cell.
        // Find the new cell, changing the position as appropriate.
        integer3 newcell = Position2Cell(*pos);
        // Important: this mapping does not apply a wrap.
        // Now we wrap, coordinating cell index and position
        while (newcell.x<0)    { newcell.x+=cpd; pos->x+=1.0; }
        while (newcell.x>=cpd) { newcell.x-=cpd; pos->x-=1.0; }
        while (newcell.y<0)    { newcell.y+=cpd; pos->y+=1.0; }
        while (newcell.y>=cpd) { newcell.y-=cpd; pos->y-=1.0; }
        while (newcell.z<0)    { newcell.z+=cpd; pos->z+=1.0; }
        while (newcell.z>=cpd) { newcell.z-=cpd; pos->z-=1.0; }
        return newcell;
    }

    inline integer3 LocalWrapToNewCell(posstruct *pos,
    int oldx, int oldy, int oldz) {
        // We have been given a particle that has drifted out of its cell.
        // Find the new cell, changing the position as appropriate.
        // This is for cell-referenced positions
        while (pos->x>halfinvcpd) {
            oldx+=1; pos->x-=invcpd;
        }
        while (pos->x<-halfinvcpd) {
            oldx-=1; pos->x+=invcpd;
        }
        while (pos->y>halfinvcpd) {
            oldy+=1; pos->y-=invcpd;
        }
        while (pos->y<-halfinvcpd) {
            oldy-=1; pos->y+=invcpd;
        }
        while (pos->z>halfinvcpd) {
            oldz+=1; pos->z-=invcpd;
        }
        while (pos->z<-halfinvcpd) {
            oldz-=1; pos->z+=invcpd;
        }
        return WrapCell(oldx, oldy, oldz);
    }

    inline void WrapAndPush(posstruct *pos, velstruct *vel, auxstruct *aux,
    int x, int y, int z, uint64 i) {
        integer3 newcell;
        
#ifdef GLOBALPOS
        newcell = WrapToNewCell(pos,x,y,z);
#else
        // Use LocalWrapToNewCell for cell-referenced positions
        newcell = LocalWrapToNewCell(pos,x,y,z);
#endif

        // Ensure that we are not trying to push a particle to a slab
        // that might already have finished
        int slab_distance = abs(x - newcell.x);
        if(slab_distance >= (P.cpd+1)/2){
            slab_distance = P.cpd - slab_distance;
        }
        if (slab_distance > FINISH_WAIT_RADIUS){
            posstruct p = *pos;
            velstruct v = *vel;
            auxstruct a = *aux;
            integer3 c = newcell;
            printf("Trying to push a particle to slab %d from slab %d.  This is larger than FINISH_WAIT_RADIUS = %d.\n",
                newcell.x, x, FINISH_WAIT_RADIUS);
            printf("pos: %e %e %e; vel: %e %e %e; aux: %llu; cell: %d %d %d\n", p.x, p.y, p.z, v.x, v.y, v.z, a.aux, c.x, c.y, c.z);
            assertf(slab_distance <= FINISH_WAIT_RADIUS,
                "Trying to push a particle to slab %d from slab %d.  This is larger than FINISH_WAIT_RADIUS = %d.",
                x, newcell.x, FINISH_WAIT_RADIUS);
        }
        
        //Push(pos,vel,aux,newcell);
        PushAt(pos,vel,aux,newcell,i);
    }

    void ShrinkIL(uint64 newlength) { 
        assertf(newlength<=length, 
        "Illegal shrinking of Insert List\n");
        length = newlength; 
    }
    
    void GrowIL(uint64 newlength) { 
        assertf(newlength>=length, 
        "Illegal growing of Insert List\n");
        assertf(newlength < maxil, 
        "Illegal resizing of Insert List\n");
        length = newlength; 
    }

    ilstruct * PartitionAndSort(int slab, uint64 *slablength);
    void CollectGaps(ilgap *gaps, int Ngaps);

    void DumpParticles(void);                // Debugging only
};


ilstruct *InsertList::PartitionAndSort(int slab, uint64 *_slablength) {
    // We've changed to an out-of-place sort.
    // This will return a pointer to NEWLY ALLOCATED memory that contains
    // the sorted particles for the slab.  These particles will be removed
    // from the InsertList.

    uint64 slablength = 0;
    FinishPartition.Start();

    uint64 mid = ParallelPartition(il, length, slab, is_in_slab);  // [0..mid-1] are not in slab, [mid..length-1] are in slab
    
    /*uint64 h = 0;
    uint64 mid = length;   // This will be one more than the last value

    while(h<mid) {
        if( il[h].xyz.x == slab ) {
            // If the particle at the end of the list is already in the right place, don't swap
            do {
                mid--;
            } while(mid > h && il[mid].xyz.x == slab);
            ilstruct tmp;
            tmp = il[mid];
            il[mid] = il[h];
            il[h] = tmp;
        }
        h++;
    }*/

    slablength = length - mid;

    FinishPartition.Stop();
    STDLOG(2, "Partition done; starting sort.\n");
    FinishSort.Start();

    ilstruct *ilnew;
    assert(posix_memalign((void **)&ilnew, 64, sizeof(ilstruct)*(slablength)) == 0);

    MultiMergeSort<ilstruct> mm;
    mm.mmsort(&(il[mid]), ilnew, slablength, cpd*cpd, 2e6, 16);
    for(uint64 i = 1 ; i < slablength; i++)
        if(ilnew[i-1].k > ilnew[i].k){
            printf("List of size %zu unsorted\n", slablength);
            FILE *f = fopen("ilparticles.dat", "wb");
            for(i = 0; i < slablength; i++)
                fwrite(&(il[mid + i].k), sizeof(unsigned int), 1, f);
            fclose(f);
            exit(0);
            break;
        }

    // only need memcpy if we want to return ilnew and free it
    //memcpy(ilnew, il + mid, sizeof(ilstruct)*slablength);
    //tbb::parallel_sort( &(il[mid]), &(il[length]), GlobalSortOperator() );
    //std::sort(ilnew, ilnew + slablength, GlobalSortOperator() );
    //pss::parallel_stable_sort( &(il[mid]), &(il[length]), GlobalSortOperator() );

    n_sorted += slablength;
    FinishSort.Stop();

    ConfirmSorting(ilnew, slablength);

    // Delete the slab particle from the insert list, since they're now in ilnew[]
    ShrinkIL(length - slablength);
    *_slablength = slablength;
    return ilnew;
}

void InsertList::DumpParticles(void) {
    for(uint64 i=0;i<length;i++) {
        posstruct p = il[i].pos;
        velstruct v = il[i].vel;
        auxstruct a = il[i].aux;
        integer3 c = il[i].xyz;
        printf("pos: %e %e %e; vel: %e %e %e; aux: %llu; cell: %d %d %d\n", p.x, p.y, p.z, v.x, v.y, v.z, a.aux, c.x, c.y, c.z);
    }
}


InsertList *IL;

std::mutex ilgap_mutex;

void ilgap::make_next_gap() {
    // Adjustments to the IL list length can only happen one at a time
    ilgap_mutex.lock();
    next = start = IL->length;
    IL->GrowIL(start+ILGAP_SIZE);
    end = IL->length;
    ilgap_mutex.unlock();
}

void InsertList::CollectGaps(ilgap *gaps, int Ngaps) {
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
    // of order NThreads * ILGAP_SIZE / 4 objects to copy.

    std::sort( gaps, gaps+Ngaps );
        // Order the gaps by their end index
    int low = 0;          // low[next:end) is unfilled space
    int high = Ngaps-1;   // high[start:next) is filled space
    // We can assume that the highest gap is the end of the IL.
    // But that doesn't always stay true, unless we do something about it.
    // Extend the filled part of each gap to the last unfilled point of the 
    // one below it.  
    for (int j=1;j<Ngaps;j++) gaps[j].start = gaps[j-1].end;

    while (high>low) {
        // When high==low, then there are no more gaps to fill
        if (gaps[low].gapsize() <= gaps[high].size()) {
            // We have enough particles to fill the gap
            uint64 copysize = gaps[low].gapsize();   // Number to move
            assert(gaps[high].next >= copysize);
            memcpy(IL->il+gaps[low].next, 
                       IL->il+gaps[high].next-copysize,
                   sizeof(ilstruct)*copysize);
            // And this means we're done with this low gap.
            gaps[high].next -= copysize;
            low++;
        } else {
            // The gap is more than the particles we have.
            uint64 copysize = gaps[high].size();   // Number to move
            assert(gaps[high].next >= copysize);
            memcpy(IL->il+gaps[low].next, 
                       IL->il+gaps[high].next-copysize,
                   sizeof(ilstruct)*copysize);
            // And this means we're done with the high set
            gaps[low].next += copysize;
            high--;
        }
    }
    IL->ShrinkIL(gaps[low].next);
    return;
}
