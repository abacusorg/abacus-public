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
//#include "tbb/parallel_sort.h"
#include "parallel_stable_sort.h"
typedef struct {
    posstruct pos;
    velstruct vel;
    auxstruct aux;
    integer3 xyz;	// The new cell, wrapped
} ilstruct;

struct GlobalSortOperator {
    inline bool operator() (const  ilstruct &pi, const ilstruct &pj ) const {
        // We must sort on pi.xyz.y*cpd+pi.xyz.z < pj.xyz.y*cpd+pj.xyz.z
        // Warning: This will get more complicated in the parallel code with 
        // the periodic wrap.
        return  ( (pi.xyz.y*P.cpd+pi.xyz.z
        <pj.xyz.y*P.cpd+pj.xyz.z) ) ; 
    }
};


// WARNING: The insert list is written to accept uint64 sizes in its list.
// So be careful when using the indices you receive from it! 


class InsertList : public grid {
public: 
    ilstruct *il;
    uint64 length;
    uint64 maxil;

    // Possibly these should be PTimer??
    STimer ILsort;
    double TotalLength;

    InsertList(int cpd, uint64 maxilsize) : grid(cpd)  { 
        length = 0; 
        maxil = maxilsize;
        il = (ilstruct *) malloc(sizeof(ilstruct) * maxilsize);
        assertf(il!=NULL,"Failed to allocate Insert List\n");
        TotalLength = 0;
    }
    ~InsertList(void) { 
        free(il);
    }

    inline void Push(posstruct  *pos, velstruct *vel, auxstruct *aux, integer3 xyz) {
        //assertf(length>=0, "Insert list has negative size!\n"); length is now unsigned so this is pointless
        assertf(length<maxil+1, "Push overflows Insert List length\n");
        il[length].pos = *pos;
        il[length].vel = *vel;
        il[length].aux = *aux;
        il[length].xyz = xyz;
        length++;
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

    void WrapAndPush(posstruct *pos, velstruct *vel, auxstruct *aux,
    int x, int y, int z) {
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
        
        Push(pos,vel,aux,newcell);
    }

    void ResetILlength(uint64 newlength) { 
        assertf(newlength>=0&&newlength<=length, 
        "Illegal resizing of Insert List\n");
        length = newlength; 
    }

    void PartitionAndSort(int slab, uint64 *slabstart, uint64 *slablength);

    void DumpParticles(void);		// Debugging only
};




void InsertList::PartitionAndSort(int slab, uint64 *slabstart, uint64 *slablength) {

    FinishPartition.Start();
    assert(length>=0);

    uint64 h = 0;
    uint64 t = length;   // This will be one more than the last value

    while(h<t) {
        if( il[h].xyz.x != slab ) {
            h++;
        }
        else {
            t--;
            if(t==h) break;
            ilstruct tmp;
            tmp = il[t];
            il[t] = il[h];
            il[h] = tmp;
        }
    }
    // At the end of this loop h==t

    //assert(t>=0); t is now unsigned so this is impossible
    assert(t<=length);

    // invariant 0 .. h - 1 (=t-1)   are not in slab 
    // t .. length-1        are in slab 

    *slabstart = t;
    *slablength = (length-1) - (t) + 1;

    FinishPartition.Stop();
    FinishSort.Start();

    //tbb::parallel_sort( &(il[t]), &(il[length]), GlobalSortOperator() );
    pss::parallel_stable_sort( &(il[t]), &(il[length]), GlobalSortOperator() );
    TotalLength += length;
    FinishSort.Stop();
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
