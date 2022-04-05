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

// TODO: Lots of VESTIGIAL code to remove.

// ======================== Insert List specific class ===================

// WARNING: The insert list is written to accept uint64 sizes in its list.
// So be careful not to re-size to 32bit when using the indices you receive from it! 

class ilstruct {
  public:
    unsigned int k;   // The key based on y*cpd+z
    int newslab;        // The new cell, wrapped
    posstruct pos;
    velstruct vel;
    auxstruct aux;

    // Some methods required for sorting
    inline unsigned int key() { return k; }
    inline void set_max_key() { k = MM_MAX_UINT; }
    bool operator< (const ilstruct& b) const { return (k<b.k); }

    // A method for dumping ASCII
    inline void print () {
        printf("pos: %e %e %e; vel: %e %e %e; aux: %lu; newslab: %d; key: %d\n", 
                pos.x, pos.y, pos.z, vel.x, vel.y, vel.z, 
                aux.aux, newslab, k);
        return;
    }

    // TOOD: not great to use P
    int celly(){
        return k/P.cpd;
    }

    int cellz(){
        return k%P.cpd;
    }
};   // end ilstruct

// Partition function; returns true if `particle` belongs in `slab`
inline bool is_in_slab(ilstruct *particle, int slab){
    return particle->newslab == slab;
}
void ConfirmSorting(ilstruct *il, uint64 len) {
    // Just an extra check that sorting is working
    for (uint64 j=0; j+1<len; j++) 
    assertf(il[j].key()<=il[j+1].key(), 
            "Insert list sorting failed: il[%d]=%d,  il[%d]=%d\n",
        j, il[j].key(), j+1, il[j+1].key());
    

    //STDLOG(1, "Sorting confirmed on insert list segment of length %d\n", len);
}

void ConfirmPartition(ilstruct *il, uint64 len, int slab) {
    // Just an extra check that partitioning is working
    int n_partitioned = 0;
    for (uint64 j=0; j<len; j++){
        if(il[j].newslab == slab)
            n_partitioned++;
        assertf((il[j].newslab == slab) == (n_partitioned > 0),
                "Insert list partitioning failed: il[%d].slab = %d; partition slab = %d; n_partitioned = %d\n",
                j, il[j].newslab, slab, n_partitioned);
    }

    STDLOG(2, "Partitioning confirmed on insert list of length %d for slab %d (found %d partitioned particles)\n", len, slab, n_partitioned);
}


class InsertList : public grid, public MultiAppendList<ilstruct> {
public: 
    uint64 n_sorted;

    STimer FinishPartition, FinishSort;

    InsertList(int cpd, uint64 maxilsize) : grid(cpd), MultiAppendList<ilstruct>(maxilsize)  { 
        n_sorted = 0;
        return;
    }
    ~InsertList(void) { 
    }

    // Push to the end of the list and grow
    inline void Push(posstruct  *pos, velstruct *vel, auxstruct *aux, integer3 xyz) {
        ilstruct il;
        il.pos = *pos;
        il.vel = *vel;
        il.aux = *aux;
        il.newslab = xyz.x;
        il.k = xyz.y*cpd + xyz.z;
        // TODO: emplace semantics
        MultiAppendList::Push(il);
    }

    inline integer3 WrapToNewCell(posstruct *pos, auxstruct *aux, int oldx, int oldy, int oldz) {
        // We have been given a particle that has drifted out of its cell.
        // Find the new cell, changing the position as appropriate.
        integer3 newcell = Position2Cell(*pos);
        // Important: this mapping does not apply a wrap.
        // Now we wrap, coordinating cell index and position
        while (newcell.x<0)    { newcell.x+=cpd; pos->x+=1.0; aux->clearLightCone(); }
        while (newcell.x>=cpd) { newcell.x-=cpd; pos->x-=1.0; aux->clearLightCone(); }
        while (newcell.y<0)    { newcell.y+=cpd; pos->y+=1.0; aux->clearLightCone(); }
        while (newcell.y>=cpd) { newcell.y-=cpd; pos->y-=1.0; aux->clearLightCone(); }
        while (newcell.z<0)    { newcell.z+=cpd; pos->z+=1.0; aux->clearLightCone(); }
        while (newcell.z>=cpd) { newcell.z-=cpd; pos->z-=1.0; aux->clearLightCone(); }
        return newcell;
    }

    inline integer3 LocalWrapToNewCell(posstruct *pos, auxstruct *aux,
        int oldx, int oldy, int oldz) {
        // We have been given a particle that has drifted out of its cell.
        // Find the new cell, changing the position as appropriate.
        // This is for cell-referenced positions
        assertf(fabs(pos->x)<0.5&&fabs(pos->y)<0.5&&fabs(pos->z)<0.5,
            "Particle has moved way too much: %f %f %f\n",
            pos->x, pos->y, pos->z);
        while (pos->x>halfinvcpd) {
            oldx+=1; pos->x-=invcpd;
            if (oldx==cpd) aux->clearLightCone();
        }
        while (pos->x<-halfinvcpd) {
            if (oldx==0) aux->clearLightCone();
            oldx-=1; pos->x+=invcpd;
        }
        while (pos->y>halfinvcpd) {
            oldy+=1; pos->y-=invcpd;
            if (oldy==cpd) aux->clearLightCone();
        }
        while (pos->y<-halfinvcpd) {
            if (oldy==0) aux->clearLightCone();
            oldy-=1; pos->y+=invcpd;
        }
        while (pos->z>halfinvcpd) {
            oldz+=1; pos->z-=invcpd;
            if (oldz==cpd) aux->clearLightCone();
        }
        while (pos->z<-halfinvcpd) {
            if (oldz==0) aux->clearLightCone();
            oldz-=1; pos->z+=invcpd;
        }
        return WrapCell(oldx, oldy, oldz);
    }

    inline void WrapAndPush(posstruct *pos, velstruct *vel, auxstruct *aux,
        int x, int y, int z) {
        integer3 newcell;
        
#ifdef GLOBALPOS
        newcell = WrapToNewCell(pos,x,y,z);
#else
        // Use LocalWrapToNewCell for cell-referenced positions
        newcell = LocalWrapToNewCell(pos,aux,x,y,z);
#endif

        /* REMOVING THIS CHECK, as it would be detected in another way,
           namely that the IL would not be empty at the end of the timestep.*/
        // Ensure that we are not trying to push a particle to a slab
        // that might already have finished
        int slab_distance = abs(x - newcell.x);
        if(slab_distance >= (P.cpd+1)/2)
            slab_distance = P.cpd - slab_distance;
        if (slab_distance > FINISH_WAIT_RADIUS){
            posstruct p = *pos;
            velstruct v = *vel;
            auxstruct a = *aux;
            integer3 c = newcell;
            printf("Trying to push a particle to slab %d from slab %d.  This is larger than FINISH_WAIT_RADIUS = %d.\n",
                newcell.x, x, FINISH_WAIT_RADIUS);
            printf("pos: %e %e %e; vel: %e %e %e; aux: %lu; cell: %d %d %d\n", p.x, p.y, p.z, v.x, v.y, v.z, a.aux, c.x, c.y, c.z);
            assertf(slab_distance <= FINISH_WAIT_RADIUS,
                "Trying to push a particle to slab %d from slab %d.  This is larger than FINISH_WAIT_RADIUS = %d.",
                x, newcell.x, FINISH_WAIT_RADIUS);
        }
        //*/   // DONE with REMOVED CHECK
        
        Push(pos,vel,aux,newcell);
    }

    ilstruct * PartitionAndSort(int slab, uint64 *slablength);

    void DumpParticles(void) {    // Debugging only
        for(uint64 i=0;i<length;i++) list[i].print();
        return;
    }
};

ilstruct *InsertList::PartitionAndSort(int slab, uint64 *_slablength) {
    // We've changed to an out-of-place sort.
    // This will return a pointer to NEWLY ALLOCATED memory that contains
    // the sorted particles for the slab.  These particles will be removed
    // from the InsertList.

    uint64 slablength = 0;
    FinishPartition.Start();

    // We've been pushing particles to spread-out locations in the MAL; now make them contiguous
    CollectGaps();

    uint64 mid = ParallelPartition(list, length, is_in_slab, slab);  // [0..mid-1] are not in slab, [mid..length-1] are in slab

    //ConfirmPartition(list, length, slab);
    
    // VESTIGIAL CODE, in case one doesn't trust the ParallelPartition code
    /*uint64 h = 0;
    uint64 mid = length;   // This will be one more than the last value

    while(h<mid) {
        if( list[h].xyz.x == slab ) {
            // If the particle at the end of the list is already in the right place, don't swap
            do {
                mid--;
            } while(mid > h && list[mid].xyz.x == slab);
            ilstruct tmp;
            tmp = list[mid];
            list[mid] = list[h];
            list[h] = tmp;
        }
        h++;
    }*/

    slablength = length - mid;

    FinishPartition.Stop();
    STDLOG(3, "Partition done, yielding %d particles; starting sort.\n", slablength);
    FinishSort.Start();

    ilstruct *ilnew;
    assert(posix_memalign((void **)&ilnew, CACHE_LINE_SIZE, sizeof(ilstruct)*(slablength)) == 0);

    MultiMergeSort<ilstruct> mm;
    mm.mmsort(&(list[mid]), ilnew, slablength, cpd*cpd, 2e6, 16);
    /* VESTIGIAL
    // The following appears to be a superfluous check as well
    for(uint64 i = 1 ; i < slablength; i++)
        if(ilnew[i-1].k > ilnew[i].k){
            printf("List of size %zu unsorted\n", slablength);
            FILE *f = fopen("ilparticles.dat", "wb");
            for(i = 0; i < slablength; i++)
                fwrite(&(list[mid + i].k), sizeof(unsigned int), 1, f);
            fclose(f);
            exit(0);
            break;
        }
    */

    // VESTIGIAL
    // only need memcpy if we want to return ilnew and free it
    // memcpy(ilnew, il + mid, sizeof(ilstruct)*slablength);
    // tbb::parallel_sort( &(il[mid]), &(il[length]), GlobalSortOperator() );
    // std::sort(ilnew, ilnew + slablength, GlobalSortOperator() );
    // pss::parallel_stable_sort( &(il[mid]), &(il[length]), GlobalSortOperator() );

    n_sorted += slablength;

    //ConfirmSorting(ilnew, slablength);

    // Delete the slab particle from the insert list, since they're now in ilnew[]
    ShrinkMAL(length - slablength);
    *_slablength = slablength;

    FinishSort.Stop();
    return ilnew;
}


InsertList *IL;
