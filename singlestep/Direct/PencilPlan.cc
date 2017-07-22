#ifndef PENCIL_PLAN_H
#define PENCIL_PLAN_H

void SinkPencilPlan::copy_into_pinned_memory(List3<FLOAT> &pinpos, int start, int total) {
    // Copy cells contiguously into pinpos->X[start..start+total), Y[), Z[)
    // where total is the padded number of particles.
    // start is the offset from the beginning of the buffers
    int cumulative_number = 0;
    FLOAT dz;
    for (int c=0; c<NFRADIUS*2+1; c++) {
        FLOAT *p = cell[c].pos;
        int N = cell[c].num;
        dz = cell[c].offset;

        memcpy(pinpos.X+start+cumulative_number, p, sizeof(FLOAT)*N);
        p+=N; 
        memcpy(pinpos.Y+start+cumulative_number, p, sizeof(FLOAT)*N);  
        p+=N; 
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

int SinkPencilPlan::load(int x, int y, int z) {
    // Given the center cell index, load the cell information
    // Return the total number of particles in the cell (un-padded)
    int total = 0;
    FLOAT cellsize = PP->invcpd;
    int width = NFRADIUS*2+1;
    for (int c=0; c<width; c++) {
        int zc = z+c-width/2;
        cell[c].pos = PP->PosXYZCell(x,y,zc).X;
        cell[c].num = PP->NumberParticle(x,y,zc);
        total += cell[c].num;
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

void SourcePencilPlan::copy_into_pinned_memory(List3<FLOAT> &pinpos, int start, int total) {
    // Copy cells contiguously into pinpos->X[start..start+total), Y[), Z[)
    // where total is the padded number of particles.
    // start is the offset from the beginning of the buffers
    int cumulative_number = 0;
    FLOAT dx;
    int width = NFRADIUS*2+1;
    for (int c=0; c<width; c++) {
        FLOAT *p = cell[c].pos;
        int N = cell[c].num;
        dx = cell[c].offset;

        FLOAT *d = pinpos.X+start+cumulative_number;
        if (dx!=0) 
            #pragma simd assert
            for (int i=0; i<N; i++) d[i] = p[i]+dx;
            else memcpy(d, p, sizeof(FLOAT)*N);  
        p+=N; 
        memcpy(pinpos.Y+start+cumulative_number, p, sizeof(FLOAT)*N);  
        p+=N; 
        memcpy(pinpos.Z+start+cumulative_number, p, sizeof(FLOAT)*N);
        cumulative_number+=N;
    }
    assertf(cumulative_number<=total, "Pencil contents exceed space supplied");
    return;
}

int SourcePencilPlan::load(int x, int y, int z) {
    // Given the center cell index, load the cell information
    // Return the total number of particles in the cell (un-padded)
    int total = 0;
    FLOAT cellsize = PP->invcpd;
    int width = NFRADIUS*2+1;
    for (int c=0; c<width; c++) {
        int xc = x+c-width/2;
        cell[c].pos = PP->PosXYZCell(xc,y,z).X;
        cell[c].num = PP->NumberParticle(xc,y,z);
        total += cell[c].num;
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
