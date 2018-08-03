/** \file This file contains the routines to prepare plans for the 
 * Sink and Source pencils, and then to actually use them to load into
 * the GPU pinned memory.  The former action is conducted by the primary
 * threads, while the latter is usually performed by the GPU thread.
 */

#ifndef PENCIL_PLAN_H
#define PENCIL_PLAN_H

/// This copies a single Sink Pencil into the supplied location in pinned memory.
/// This requires that we convert from cell-centered coordinates to
/// a coordinate system centered on the middle cell of the pencil.

void SinkPencilPlan::copy_into_pinned_memory(List3<FLOAT> &pinpos, int start, int total) {
    // Copy cells contiguously into pinpos->X[start..start+total), Y[), Z[)
    // where total is the padded number of particles.
    // start is the offset from the beginning of the buffers
    int cumulative_number = 0;
    FLOAT dz;
    for (int c=0; c<NFRADIUS*2+1; c++) {
        List3<FLOAT> p = cell[c].pos;
        int N = (int) p.N;
        dz = cell[c].offset;

        memcpy(pinpos.X+start+cumulative_number, p.X, sizeof(FLOAT)*N);
        memcpy(pinpos.Y+start+cumulative_number, p.Y, sizeof(FLOAT)*N);  
        FLOAT *d = pinpos.Z+start+cumulative_number;
        if (dz!=0) 
            #pragma simd assert
            for (int i=0; i<N; i++) d[i] = p.Z[i]+dz;
        else memcpy(d, p.Z, sizeof(FLOAT)*N);  
        cumulative_number+=N;
    }
    assertf(cumulative_number<=total, "Pencil contents exceed space supplied");
    return;
}

/// This loads up the plan for a SinkPencil with pointers to the 
/// primary PosXYZ data, as well as the needed coordinate offsets.

int SinkPencilPlan::load(int x, int y, int z) {
    // Given the center cell index, load the cell information
    // Return the total number of particles in the cell (un-padded)
    int total = 0;
    FLOAT cellsize = PP->invcpd;
    int width = NFRADIUS*2+1;
    for (int c=0; c<width; c++) {
        int zc = z+c-width/2;
        cell[c].pos = PP->PosXYZCell(x,y,zc);
        total += (int) cell[c].pos.N;
        #ifndef GLOBAL_POS
            // Local positions, just offset the cells
            cell[c].offset = (c-width/2)*cellsize;
        #else
            // Can use the z cell number to do this.
            cell[c].offset = (zc-PP->WrapSlab(zc))*cellsize;
        #endif
    }
    return total;
}

/// This copies a single Source Pencil into the supplied location in pinned memory.
/// This requires that we convert from cell-centered coordinates to
/// a coordinate system centered on the middle cell of the pencil.

void SourcePencilPlan::copy_into_pinned_memory(List3<FLOAT> &pinpos, int start, int total) {
    // Copy cells contiguously into pinpos->X[start..start+total), Y[), Z[)
    // where total is the padded number of particles.
    // start is the offset from the beginning of the buffers
    int cumulative_number = 0;
    FLOAT dx;
    int width = NFRADIUS*2+1;
    for (int c=0; c<width; c++) {
        List3<FLOAT> p = cell[c].pos;
        int N = (int) p.N;
        dx = cell[c].offset;

        FLOAT *d = pinpos.X+start+cumulative_number;
        if (dx!=0) 
            #pragma simd assert
            for (int i=0; i<N; i++) d[i] = p.X[i]+dx;
        else memcpy(d, p.X, sizeof(FLOAT)*N);
        memcpy(pinpos.Y+start+cumulative_number, p.Y, sizeof(FLOAT)*N);
        memcpy(pinpos.Z+start+cumulative_number, p.Z, sizeof(FLOAT)*N);
        cumulative_number+=N;
    }
    assertf(cumulative_number<=total, "Pencil contents exceed space supplied");
    return;
}

/// This loads up the plan for a SourcePencil with pointers to the 
/// primary PosXYZ data, as well as the needed coordinate offsets.

int SourcePencilPlan::load(int x, int y, int z) {
    // Given the center cell index, load the cell information
    // Return the total number of particles in the cell (un-padded)
    int total = 0;
    FLOAT cellsize = PP->invcpd;
    int width = NFRADIUS*2+1;
    for (int c=0; c<width; c++) {
        int xc = x+c-width/2;
        cell[c].pos = PP->PosXYZCell(xc,y,z);
        total += (int) cell[c].pos.N;
        #ifndef GLOBAL_POS
            // Local positions, just offset the cells
            cell[c].offset = (c-width/2)*cellsize;
        #else
            // Can use the x cell number to do this.
            cell[c].offset = (xc-PP->WrapSlab(xc))*cellsize;
        #endif
    }
    return total;
}

#endif
