// TODO: This makes up a global variable NFWIDTH; need to change that.

class CellPencilPlan {
    FLOAT *pos;    
    // Cells are assumed to be x[0..N), y[0..N), z[0..N), contiguous,
    // with x[0] being the given position.
    int num;
    // Number of particles in the cell
    FLOAT offset;
    // The offset to be applied to x or z, relative to the center cell
}

class SinkPencilPlan {
  public:
    CellPencilPlan cell[2*NFWIDTH+1];
    // The cells are not assumed to be contiguous (e.g., periodic wraps)

    void copy_into_pinned_memory(List3<FLOAT> &pinpos, int start, int total) {
	// Copy cells contiguously into pinpos->X[start..start+total), Y[), Z[)
	// where total is the padded number of particles.
	// start is the offset from the beginning of the buffers
	int cumulative_number = 0;
	FLOAT dz;
	for (int c=0; c<2*NFWIDTH+1; c++) {
	    FLOAT *p = cell[c].pos;
	    int N = cell[c].num;
	    dz = cell[c].offset;

	    memcpy(pinpos.X+start+cumulative_number, p, sizeof(FLOAT)*N);
	    p+=N; 
	    memcpy(pinpos.Y+start+cumulative_number, p, sizeof(FLOAT)*N);  
	    p+=N; 
	    FLOAT *d = pinpos.Z+start+cumulative_number;
	    if (dz!=0) for (int i=0, i<N; i++) d[i] = p[i]+dz;
		else memcpy(d, p, sizeof(FLOAT)*N);  
	    cumulative_number+=N;
	}
	assertf(cumulative_number<=total, "Pencil contents exceed space supplied");
	return;
    }

    int load(int x, int y, int z) {
        // Given the center cell index, load the cell information
	// Return the total number of particles in the cell (un-padded)
	int total = 0;
	FLOAT cellsize = PP->invcpd;
	for (int c=0; c<= 2*NFWIDTH+1; c++) {
	    int zc = z+c-NFWIDTH;
	    cell[c].pos = (FLOAT *)PP->PosXYZCell(x,y,zc);
	    cell[c].num = PP->NumberParticle(x,y,zc);
	    total += cell[c].num;
	    #ifndef GLOBAL_POS
		// Local positions, just offset the cells
		cell[c].offset = (c-NFWIDTH)*cellsize;
	    #else
		// Can use the z cell number to do this.
		cell[c].offset = (zc-PP->WrapSlab(zc))*cellsize;
	    #endif
	}
	return total;
    }
};


class SourcePencilPlan {
    // TODO: The same as above, but with motion in the x direction
};

