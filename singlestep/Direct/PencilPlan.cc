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

void SinkPencilPlan::copy_into_pinned_memory(List3<FLOAT> &pinpos, int start, int total, void *SinkPosSlab) {
    // Copy cells contiguously into pinpos->X[start..start+total), Y[), Z[)
    // where total is the padded number of particles.
    // start is the offset from the beginning of the buffers
    int cumulative_number = 0;
    FLOAT dz;
    FLOAT cellsize = PP->invcpd;
    for (int c=0; c<NFRADIUS*2+1; c++) {
	// The pointer math has to be in posstruct; then cast
	FLOAT *p = (FLOAT *)((posstruct *)SinkPosSlab+cell[c].start);
	int N = cell[c].N;

	#ifdef GLOBAL_POS
        dz = cell[c].offset;
	#else
	dz = (c-NFRADIUS)*cellsize;
	#endif

        memcpy(pinpos.X+start+cumulative_number, p, sizeof(FLOAT)*N); p+=N;
        memcpy(pinpos.Y+start+cumulative_number, p, sizeof(FLOAT)*N); p+=N;
        FLOAT *d = pinpos.Z+start+cumulative_number;
        if (dz!=0) 
            #pragma simd assert
            for (int i=0; i<N; i++) d[i] = p[i]+dz;
        else memcpy(d, p, sizeof(FLOAT)*N);  
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
	cellinfo *info = PP->CellInfo(x,y,zc);
        cell[c].start = info->startindex;
        cell[c].N = info->count;
        //? cell[c].pos = PP->PosXYZCell(x,y,zc);
        total += (int) cell[c].N;
        #ifdef GLOBAL_POS
            // Can use the z cell number to do this.
            cell[c].offset = (zc-PP->WrapSlab(zc))*cellsize;
        #endif
    }
    return total;
}

/// Given the padded list of accelerations in pinned memory, 
/// copy into the unpadded cells in the given PartialAccSlab.  
/// This routine is called by the GPU code. 
void SinkPencilPlan::copy_from_pinned_memory(void *_pinacc, int start, 
	int total, void *SinkPartialAccSlab) {
    // Copy pinacc, which holds the padded list 
    // accstruct[start..start+total)
    //  cells contiguously into pinpos->X[start..start+total), Y[), Z[)
    // where total is the padded number of particles.
    // start is the offset from the beginning of the buffers
    accstruct *pinacc = (accstruct *)_pinacc;
    int cumulative_number = 0;
    for (int c=0; c<NFRADIUS*2+1; c++) {
	int N = cell[c].N;
	// The pointer math has to be in accstruct; then cast
	accstruct *p = (accstruct *)SinkPartialAccSlab+cell[c].start;

        memcpy(p, pinacc+start+cumulative_number, sizeof(accstruct)*N); 
        cumulative_number+=N;
    }
    assertf(cumulative_number<=total, "Pencil contents exceed space supplied");
    return;
}





/// This copies a single Source Pencil into the supplied location in pinned memory.
/// This requires that we convert from cell-centered coordinates to
/// a coordinate system centered on the middle cell of the pencil.

void SourcePencilPlan::copy_into_pinned_memory(List3<FLOAT> &pinpos, int start, int total, void **SourcePosSlab) {
    // Copy cells contiguously into pinpos->X[start..start+total), Y[), Z[)
    // where total is the padded number of particles.
    // start is the offset from the beginning of the buffers
    int cumulative_number = 0;
    FLOAT dx;
    FLOAT cellsize = PP->invcpd;
    int width = NFRADIUS*2+1;
    for (int c=0; c<width; c++) {
	FLOAT *p = (FLOAT *)((posstruct *)SourcePosSlab[c]+cell[c].start);
	int N = cell[c].N;
	#ifdef GLOBAL_POS
	    dx = cell[c].offset;
	#else 
            dx = (c-NFRADIUS)*cellsize;
	#endif

        FLOAT *d = pinpos.X+start+cumulative_number;
        if (dx!=0) 
            #pragma simd assert
            for (int i=0; i<N; i++) d[i] = p[i]+dx;
        else memcpy(d, p, sizeof(FLOAT)*N);
        p+=N; memcpy(pinpos.Y+start+cumulative_number, p, sizeof(FLOAT)*N);
        p+=N; memcpy(pinpos.Z+start+cumulative_number, p, sizeof(FLOAT)*N);
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
	cellinfo *info = PP->CellInfo(xc,y,z);
	cell[c].start = info->startindex;
	cell[c].N = info->count;

        //? cell[c].pos = PP->PosXYZCell(xc,y,z);
        total += (int) cell[c].N;
        #ifdef GLOBAL_POS
            // Can use the x cell number to do this.
            cell[c].offset = (xc-PP->WrapSlab(xc))*cellsize;
        #endif
    }
    return total;
}

#endif
