/** \file This file contains the routines to prepare plans for the 
Sink and Source pencils, and then to actually use them to load into
the GPU pinned memory.  The former action is conducted by the primary
threads, while the latter is usually performed by the GPU thread.

 */

#ifndef PENCIL_PLAN_H
#define PENCIL_PLAN_H

/// This copies a single Sink Pencil into the supplied location in pinned memory.
/// This requires that we convert from cell-centered coordinates to
/// a coordinate system centered on the middle cell of the pencil.

void SinkPencilPlan::copy_into_pinned_memory(List3<FLOAT> &pinpos, int start, int total, void *SinkPosSlab, int NearFieldRadius, uint64 Nslab) {
    // This version uses the slab XYZ ordering
    // and copies the X and Y pencils in one big go, instead of cell-by-cell

    // Copy cells contiguously into pinpos->X[start..start+total), Y[), Z[)
    // where total is the padded number of particles.
    // start is the offset from the beginning of the buffers
    int cumulative_number = 0;
    FLOAT dz;
    FLOAT cellsize = PP->invcpd;
    for (int c=0; c<NearFieldRadius*2+1; c++) {
        FLOAT *p = (FLOAT *)SinkPosSlab + cell[c].start + 2*Nslab;
        int N = cell[c].N;

        #ifdef GLOBAL_POS
        dz = cell[c].offset;
        #else
        dz = (c-NearFieldRadius)*cellsize;
        #endif
        FLOAT *d = pinpos.Z+start+cumulative_number;
        if (dz!=0) 
            // TODO: maybe there'd be a small gain in precomputing offsets and merging some of these loops
            #pragma simd assert
            for (int i=0; i<N; i++)
                d[i] = p[i]+dz;
        else
            memcpy(d, p, sizeof(FLOAT)*N);
        cumulative_number+=N;
    }
    assertf(cumulative_number==total, "Pencil contents doesn't match space supplied: %d vs %d\n", cumulative_number, total);

    // Now do xy in as few chunks as possible
    FLOAT *p1 = (FLOAT *)SinkPosSlab + cell[0].start;
    memcpy(pinpos.X+start, p1, sizeof(FLOAT)*contig);
    if(wrapcell >= 0){
        FLOAT *p2 = (FLOAT *)SinkPosSlab + cell[wrapcell].start;
        memcpy(pinpos.X+start+contig, p2, sizeof(FLOAT)*(cumulative_number-contig));
    }

    p1 = (FLOAT *)SinkPosSlab + cell[0].start + Nslab;
    memcpy(pinpos.Y+start, p1, sizeof(FLOAT)*contig);
    if(wrapcell >= 0){
        FLOAT *p2 = (FLOAT *)SinkPosSlab + cell[wrapcell].start + Nslab;
        memcpy(pinpos.Y+start+contig, p2, sizeof(FLOAT)*(cumulative_number-contig));
    }
}

/// This loads up the plan for a SinkPencil with pointers to the 
/// primary PosXYZ data, as well as the needed coordinate offsets.

int SinkPencilPlan::load(int x, int y, int z, int NearFieldRadius) {
    // Given the center cell index, load the cell information
    // Return the total number of particles in the cell (un-padded)
    contig = 0;
    wrapcell = -1;

    int total = 0;
    FLOAT cellsize = PP->invcpd;
    const int width = NearFieldRadius*2+1;
    for (int c=0; c<width; c++) {
        int zc = z+c-width/2;
        cellinfo *info = PP->CellInfo(x,y,zc);
        cell[c].start = info->startindex;
        cell[c].N = info->count;
        total += (int) cell[c].N;
        #ifdef GLOBAL_POS
            // Can use the z cell number to do this.
            cell[c].offset = (zc-PP->WrapSlab(zc))*cellsize;
        #endif

        // We could simply use the cell coords, but best to check cell continuity for sanity
        if(c > 0 && cell[c].start != cell[c-1].start + cell[c-1].N){
            assert(wrapcell == -1);
            wrapcell = c;
        } else if (wrapcell < 0) {
            contig += cell[c].N;
        }
    }
    return total;
}

void SinkPencilPlan::copy_from_pinned_memory(void *_pinacc, int start, 
    int total, void *SinkAccSlab, int sinkindex, int NearFieldRadius, uint64 Nslab) {
    // Copy pinacc, which holds the padded list 
    // accstruct[start..start+total)
    //  cells contiguously into pinpos->X[start..start+total), Y[), Z[)
    // where total is the padded number of particles.
    // start is the offset from the beginning of the buffers
    accstruct *pinacc = (accstruct *)_pinacc+start;
    int cumulative_number = 0;
    int k = sinkindex%P.cpd;    // Reduce this to the Z index for this pencil
    for(int c = 0; c < 2*NearFieldRadius + 1; c++)
        cumulative_number += cell[c].N;
    assertf(cumulative_number==total, "Pencil contents doesn't match space supplied");

    accstruct *p = (accstruct *)SinkAccSlab+cell[0].start;
    accstruct *pin = pinacc;
    // At the beginning of the row, we get to set the whole thing
    if(k == 0)
        memcpy(p, pin, sizeof(accstruct)*total);
    else {
        // We have two to three segments to process:
        // 1) Before the wrap [co-add]
        // 2) Possibly: after the wrap [co-add]
        // 3) and the last cell [set, unless past the wrap]
        int last_N = cell[2*NearFieldRadius].N;
        int coadd_contig = contig;
        if (wrapcell < 0)
            coadd_contig -= last_N;
        assert(coadd_contig >= 0);

        // first segment
        for (int t=0; t<coadd_contig; t++)
            p[t] += pin[t];
        pin += coadd_contig;

        if (wrapcell >= 0 && wrapcell < 2*NearFieldRadius){
            // second segment
            p = (accstruct *) SinkAccSlab + cell[wrapcell].start;
            coadd_contig = total - last_N - coadd_contig;
            assert(coadd_contig >= 0);
            for(int t=0; t < coadd_contig; t++)
                p[t] += pin[t];
            pin += coadd_contig;
        }
        // third segment
        p = (accstruct *) SinkAccSlab + cell[2*NearFieldRadius].start;
        if(wrapcell >= 0)
            for(int t=0; t < last_N; t++)
                p[t] += pin[t];
        else
            memcpy(p, pin, sizeof(accstruct)*last_N);
    }
}


/// This copies a single Source Pencil into the supplied location in pinned memory.
/// This requires that we convert from cell-centered coordinates to
/// a coordinate system centered on the middle cell of the pencil.

void SourcePencilPlan::copy_into_pinned_memory(List3<FLOAT> &pinpos, int start, int total, void **SourcePosSlab, int NearFieldRadius, uint64 *Nslab) {
    // Copy cells contiguously into pinpos->X[start..start+total), Y[), Z[)
    // where total is the padded number of particles.
    // start is the offset from the beginning of the buffers
    int cumulative_number = 0;
    FLOAT dx;
    FLOAT cellsize = PP->invcpd;
    int width = NearFieldRadius*2+1;
    for (int c=0; c<width; c++) {
	int N = cell[c].N;
    FLOAT *p = (FLOAT *)SourcePosSlab[c] + cell[c].start;
	#ifdef GLOBAL_POS
	    dx = cell[c].offset;
	#else 
            dx = (c-NearFieldRadius)*cellsize;
	#endif

        FLOAT *d = pinpos.X+start+cumulative_number;
        if (dx!=0) 
            #pragma simd assert
            for (int i=0; i<N; i++) d[i] = p[i]+dx;
        else memcpy(d, p, sizeof(FLOAT)*N);
        p+=Nslab[c]; memcpy(pinpos.Y+start+cumulative_number, p, sizeof(FLOAT)*N);
        p+=Nslab[c]; memcpy(pinpos.Z+start+cumulative_number, p, sizeof(FLOAT)*N);
        cumulative_number+=N;
    }
    assertf(cumulative_number==total, "Pencil contents doesn't match space supplied");
    return;
}

/// This loads up the plan for a SourcePencil with pointers to the 
/// primary PosXYZ data, as well as the needed coordinate offsets.

int SourcePencilPlan::load(int x, int y, int z, int NearFieldRadius) {
    // Given the center cell index, load the cell information
    // Return the total number of particles in the cell (un-padded)
    int total = 0;
    FLOAT cellsize = PP->invcpd;
    const int width = NearFieldRadius*2+1;
    for (int c=0; c<width; c++) {
        int xc = x+c-width/2;
	cellinfo *info = PP->CellInfo(xc,y,z);
	cell[c].start = info->startindex;
	cell[c].N = info->count;
        total += (int) cell[c].N;
        #ifdef GLOBAL_POS
            // Can use the x cell number to do this.
            cell[c].offset = (xc-PP->WrapSlab(xc))*cellsize;
        #endif
    }
    return total;
}

#endif
