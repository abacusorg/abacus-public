/* insert.cpp

Class for the insert list.

Particles are added to the end of the insert list, along with their
new cell number.  When requested, a partioning is done to place particles 
of a given slab at the end of the list; this slab is then sorted into
cell order.  We will then copy those particles out and delete them from
the insert list.

Importantly, particle positions are adjusted when the particle is 
added to the insert list, so that the particle has a legal position
for its *new* cell.
*/

//for intel's fast parallel sort
// #include "tbb/parallel_sort.h"
// #include "parallel_stable_sort.h"

#include "multimergesort.cpp"


class ilstruct {
  public:
    unsigned int k;   // The key based on y*cpd+z
    integer3 xyz;	// The new cell, wrapped
    posstruct pos;
    velstruct vel;
    auxstruct aux;

    // Some methods require for sorting
    inline unsigned int key() { return k; }
    inline void set_max_key() { k = MM_MAX_UINT; }
    bool operator< (const ilstruct& b) const { return (k<b.k); }
};   // end ilstruct



// The following routine is vestigial
struct GlobalSortOperator {
    inline bool operator() (const  ilstruct &pi, const ilstruct &pj ) const {
        // We must sort on pi.xyz.y*cpd+pi.xyz.z < pj.xyz.y*cpd+pj.xyz.z
        // Warning: This will get more complicated in the parallel code with 
        // the periodic wrap.
        // two conditionals instead of multiplication?
        return  ( (pi.xyz.y*P.cpd+pi.xyz.z
        <pj.xyz.y*P.cpd+pj.xyz.z) ) ; 
    }
};

void ConfirmSorting(ilstruct *il, uint64 len) {
    for (j=0; j+1<len; j++) 
    assertf(il[j].key()<=il[j+1].key(), 
    	"Insert list sorting failed: il[%d]=%d,  il[%d]=%d\n",
	j, il[j].key(), j+1, il[j+1].key());
    return;
}


// WARNING: The insert list is written to accept uint64 sizes in its list.
// So be careful when using the indices you receive from it! 


class InsertList : public grid {
public: 
    ilstruct *il;
    uint64 length;
    uint64 maxil;

    // Possibly these should be PTimer??
    STimer ILsort;
    uint64 n_sorted;

    InsertList(int cpd, uint64 maxilsize) : grid(cpd)  { 
        length = 0; 
        maxil = maxilsize;
        int ret = posix_memalign((void **) &il, 4096, sizeof(ilstruct) * maxilsize);
        assertf(ret==0,"Failed to allocate Insert List\n");
        n_sorted = 0;
    }
    ~InsertList(void) { 
        free(il);
    }

    // Push to the end of the list and grow
    inline void Push(posstruct  *pos, velstruct *vel, auxstruct *aux, integer3 xyz) {
        //assertf(length>=0, "Insert list has negative size!\n"); length is now unsigned so this is pointless
        assertf(length<maxil+1, "Push overflows Insert List length\n");
        il[length].pos = *pos;
        il[length].vel = *vel;
        il[length].aux = *aux;
        il[length].xyz = xyz;
	il[length].k = xyz.y*(P->cpd)+xyz.z;
        length++;
    }
    
    // Push to a given location; don't grow
    inline void PushAt(posstruct  *pos, velstruct *vel, auxstruct *aux, integer3 xyz, uint64 i) {
        assertf(i<length, "Push overflows Insert List length\n");
        il[i].pos = *pos;
        il[i].vel = *vel;
        il[i].aux = *aux;
        il[i].xyz = xyz;
	il[i].k = xyz.y*(P->cpd)+xyz.z;
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
            /*assertf(slab_distance <= FINISH_WAIT_RADIUS,
                "Trying to push a particle to slab %d from slab %d.  This is larger than FINISH_WAIT_RADIUS = %d.",
                x, newcell.x, FINISH_WAIT_RADIUS);*/
        }
        
        //Push(pos,vel,aux,newcell);
        PushAt(pos,vel,aux,newcell,i);
    }

    void ResetILlength(uint64 newlength) { 
        assertf(newlength<=length, 
        "Illegal resizing of Insert List\n");
        length = newlength; 
    }

    ilstruct * PartitionAndSort(int slab, uint64 *slablength);

    void DumpParticles(void);		// Debugging only
};




ilstruct * InsertList::PartitionAndSort(int slab, uint64 *slablength) {
    // We've changed to an out-of-place sort.
    // This will return a pointer to NEWLY ALLOCATED memory that contains
    // the sorted particles for the slab.  These particles will be removed
    // from the InsertList.

    FinishPartition.Start();

    uint64 h = 0;
    uint64 t = length;   // This will be one more than the last value

    while(h<t) {
        if( il[h].xyz.x == slab ) {
            // If the particle at the end of the list is already in the right place, don't swap
            do {
                t--;
            } while(t > h && il[t].xyz.x == slab);
            ilstruct tmp;
            tmp = il[t];
            il[t] = il[h];
            il[h] = tmp;
        }
        h++;
    }

    // invariant 0 .. t-1   are not in slab 
    // t .. length-1        are in slab

    *slabstart = t;
    *slablength = length - t;

    FinishPartition.Stop();
    FinishSort.Start();

    ilstruct *ilnew;
    assert(posix_memalign((void **)&ilnew, 64, sizeof(ilstruct)*(*slablength)==0);

    MultiMergeSort<ilstruct> mm;
    mm.mmsort(&(il[t]), ilnew, *slablength, P->cpd*P->cpd, 2e6, 16);

    // tbb::parallel_sort( &(il[t]), &(il[length]), GlobalSortOperator() );
    //std::sort( &(il[t]), &(il[length]), GlobalSortOperator() );
    //pss::parallel_stable_sort( &(il[t]), &(il[length]), GlobalSortOperator() );

    n_sorted += *slablength;
    FinishSort.Stop();

    // Delete the slab particle from the insert list, since they're now in ilnew[]
    IL->ResetILlength(IL->length - ilslablength);
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
