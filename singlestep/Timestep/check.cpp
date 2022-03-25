/* check.cpp 
 * Future general validation code may be placed here.
 * Alternatively, this file could easily be refactored into the cellinfo class
 */

int are_cellinfo_legal(int slab, uint64 slabsize, uint64 slabsize_with_ghost) {
    // Return 1 if ok, return 0 if not.
    int cpd = CP->cpd;
    for (int j=0; j<cpd; j++)
	for (int k=0; k<cpd; k++) {
	    if (CP->CellInfo(slab,j,k)->legalvalue(slabsize, slabsize_with_ghost)==0) {
		assertf(0==9, "Failed for slab %d, j=%d, k=%d\n", slab, j, k);
	        return 0;
	    }
	}
    return 1;
}
