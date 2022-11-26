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

void SinkPencilPlan::copy_into_pinned_memory(List3<FLOAT> &pinpos, int start, int total, FLOAT *SinkPosSlab, int NearFieldRadius, uint64 Nslab) {
    // This version uses the slab XYZ ordering
    // and copies the X and Y pencils in one big go, instead of cell-by-cell

    // Copy cells contiguously into pinpos->X[start..start+total), Y[), Z[)
    // where total is the padded number of particles.
    // start is the offset from the beginning of the buffers
    int nwritten = 0;
    FLOAT dz;
    FLOAT cellsize = CP->invcpd;
    int width = NearFieldRadius*2+1;
    for (int c=0; c<width; c++) {
        FLOAT *p = SinkPosSlab + cell[c].start + ghostoffset + 2*Nslab;
        int N = cell[c].N;

        #ifdef GLOBAL_POS
        dz = cell[c].offset;
        #else
        dz = (c-NearFieldRadius)*cellsize;
        #endif
        FLOAT *d = pinpos.Z+start+nwritten;
        if (dz!=0){
            // TODO: maybe there'd be a small gain in precomputing offsets and merging some of these loops
            #pragma omp simd
            for (int i=0; i<N; i++)
                d[i] = p[i]+dz;
        }
        else {
            memcpy(d, p, sizeof(FLOAT)*N);
        }

        /*// cell version of XY
        p = SinkPosSlab + cell[c].start + ghostoffset;
        d = pinpos.X + start + nwritten;
        memcpy(d, p, sizeof(FLOAT)*N);
        p = SinkPosSlab + cell[c].start + ghostoffset + Nslab;
        d = pinpos.Y + start + nwritten;
        memcpy(d, p, sizeof(FLOAT)*N);*/

        nwritten+=N;
    }
    assertf(nwritten==total, "Pencil contents doesn't match space supplied: %d vs %d\n",
        nwritten, total);

    // Now do xy in two chunks: before and after the wrap
    int N1 = 0, N2 = 0;
    if(wrapcell > 0){
        N1 = cell[wrapcell-1].start + cell[wrapcell-1].N - cell[0].start;
    } else {
        N1 = 0;
    }
    N2 = cell[width-1].start + cell[width-1].N - cell[wrapcell].start;
    assertf((N1 + N2) == total, "N1 (%d) + N2 (%d) != total (%d)\n", N1, N2, total);

    memcpy(pinpos.X+start,    SinkPosSlab+cell[0].start+ghostoffset,        sizeof(FLOAT)*N1);
    memcpy(pinpos.X+start+N1, SinkPosSlab+cell[wrapcell].start+ghostoffset, sizeof(FLOAT)*N2);

    memcpy(pinpos.Y+start,    SinkPosSlab+Nslab+cell[0].start+ghostoffset,        sizeof(FLOAT)*N1);
    memcpy(pinpos.Y+start+N1, SinkPosSlab+Nslab+cell[wrapcell].start+ghostoffset, sizeof(FLOAT)*N2);
}

/// This loads up the plan for a SinkPencil with pointers to the 
/// primary PosXYZ data, as well as the needed coordinate offsets.

int SinkPencilPlan::load(int x, int y, int z, int nfradius, int truncate) {
    // Given the center cell index, load the cell information
    // Return the total number of particles in the cell (un-padded)
    wrapcell = 0;

    cellinfo *startci = CP->CellInfo(x,y,node_z_start);
    ghostoffset = startci->startindex_with_ghost - startci->startindex;
    assert(ghostoffset >= 0);

    int total = 0;
    FLOAT cellsize = CP->invcpd;
    const int width = nfradius*2+1;
    for (int c=0; c<width; c++) {
        int zc = z + c - nfradius;

        int rightdist = CP->WrapSlab(zc - node_z_start);

        if( truncate && (rightdist >= node_z_size) ){
            // Edge of the domain, load a no-op cell
            cell[c].start = 0;
            cell[c].N = 0;
        } else {
            // Regular cell
            cellinfo *info = CP->CellInfo(x,y,zc);
            cell[c].start = info->startindex;
            cell[c].N = info->count;
        }

        // Identify if we cross the wrap, or leave the primary zone
        if(CP->WrapSlab(zc - (node_z_start + node_z_size)) == 0 ||
            CP->WrapSlab(zc) == node_z_start){
            assert(wrapcell == 0);
            wrapcell = c;
        }

        total += (int) cell[c].N;
        #ifdef GLOBAL_POS
            // Can use the z cell number to do this.
            cell[c].offset = (zc-CP->WrapSlab(zc))*cellsize;
        #endif
    }

    return total;
}

void SinkPencilPlan::copy_from_pinned_memory(void *_pinacc, int start, 
    int total, void *SinkAccSlab, int k, int NearFieldRadius, uint64 Nslab) {

    // Copy pinacc, which holds the padded list 
    // accstruct[start..start+total)
    //  cells contiguously into pinpos->X[start..start+total), Y[), Z[)
    // where total is the padded number of particles.
    // start is the offset from the beginning of the buffers
    accstruct *pinacc = (accstruct *)_pinacc+start;
    int cumulative_number = 0;
    int width = 2*NearFieldRadius + 1;
    for(int c = 0; c < width; c++)
        cumulative_number += cell[c].N;
    assertf(cumulative_number==total, "Pencil contents doesn't match space supplied");

    accstruct *p = (accstruct *)SinkAccSlab+cell[0].start;
    accstruct *pin = pinacc;
    // At the beginning of the row, we get to set the whole thing
    int nwritten = 0;
    if(k == 0){
        memcpy(p, pin, sizeof(accstruct)*total);
        nwritten += total;
    }
    else {
        int last_N;

        // count how many particles before the wrap or the last cell
        int coadd_contig;
        if(wrapcell > 0){
            coadd_contig = cell[wrapcell-1].start + cell[wrapcell-1].N - cell[0].start;
        } else {
            last_N = cell[2*NearFieldRadius].N;
            coadd_contig = total - last_N;
        }
        assert(coadd_contig >= 0);

        // up to the last cell, or the wrap: co-add
        for (int t=0; t<coadd_contig; t++)
            p[t] += pin[t];
        pin += coadd_contig;
        nwritten += coadd_contig;

        if (wrapcell > 0){
            // past the wrap: co-add
            p = (accstruct *) SinkAccSlab + cell[wrapcell].start;
            coadd_contig = total - coadd_contig;
            assert(coadd_contig >= 0);
            for(int t=0; t < coadd_contig; t++)
                p[t] += pin[t];
            nwritten += coadd_contig;
        } else {
            // last cell, not past the wrap: set
            p = (accstruct *) SinkAccSlab + cell[2*NearFieldRadius].start;
            memcpy(p, pin, sizeof(accstruct)*last_N);
            nwritten += last_N;
        }
    }

    assert(nwritten == total);

    // Check that all particles have finite accel
    /*for(int c = 0; c < 2*NearFieldRadius+1; c++){
        accstruct *p = (accstruct *) SinkAccSlab + cell[c].start;
        for(int i = 0; i < cell[c].N; i++){
            assertf(TOFLOAT3(p[i]).is_finite(), "p[%d of %d] in cell %d (wrapcell %d): %f %f %f\n", i, cell[c].N, c, wrapcell, p[i].x, p[i].y, p[i].z);
            //assertf(p[i].norm2() != 0., "p[%d of %d] in cell %d (wrapcell %d): %f %f %f\n", i, cell[c].N, c, wrapcell, p[i].x, p[i].y, p[i].z);
        }
    }*/
}


/// This copies a single Source Pencil into the supplied location in pinned memory.
/// This requires that we convert from cell-centered coordinates to
/// a coordinate system centered on the middle cell of the pencil.

void SourcePencilPlan::copy_into_pinned_memory(List3<FLOAT> &pinpos, int start, int total, FLOAT **SourcePosSlab, int NearFieldRadius, uint64 *Nslab) {
    // Copy cells contiguously into pinpos->X[start..start+total), Y[), Z[)
    // where total is the padded number of particles.
    // start is the offset from the beginning of the buffers
    int cumulative_number = 0;
    FLOAT dx;
    FLOAT cellsize = CP->invcpd;
    int width = NearFieldRadius*2+1;
    for (int c=0; c<width; c++) {
        int N = cell[c].N;
        FLOAT *p = SourcePosSlab[c] + cell[c].start;
        #ifdef GLOBAL_POS
            dx = cell[c].offset;
        #else 
            dx = (c-NearFieldRadius)*cellsize;
        #endif

        FLOAT *d = pinpos.X+start+cumulative_number;
        if (dx!=0) 
            #pragma omp simd
            for (int i=0; i<N; i++) d[i] = p[i]+dx;
        else memcpy(d, p, sizeof(FLOAT)*N);
        p+=Nslab[c]; memcpy(pinpos.Y+start+cumulative_number, p, sizeof(FLOAT)*N);
        p+=Nslab[c]; memcpy(pinpos.Z+start+cumulative_number, p, sizeof(FLOAT)*N);
        cumulative_number+=N;
    }
    assertf(cumulative_number==total, "Pencil contents doesn't match space supplied");
}

/// This loads up the plan for a SourcePencil with pointers to the 
/// primary PosXYZ data, as well as the needed coordinate offsets.

int SourcePencilPlan::load(int x, int y, int z, int NearFieldRadius) {
    // Given the center cell index, load the cell information
    // Return the total number of particles in the cell (un-padded)
    int total = 0;
    FLOAT cellsize = CP->invcpd;
    const int width = NearFieldRadius*2+1;
    for (int c=0; c<width; c++) {
        int xc = x+c-width/2;
        
        cellinfo *info = CP->CellInfo(xc,y,z);
        cell[c].start = info->startindex_with_ghost;
        cell[c].N = info->count;

        total += (int) cell[c].N;
        #ifdef GLOBAL_POS
            // Can use the x cell number to do this.
            cell[c].offset = (xc-CP->WrapSlab(xc))*cellsize;
        #endif
    }
    return total;
}

#endif
