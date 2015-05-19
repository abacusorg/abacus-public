/* check.cpp 
Simple utilities  */

int are_cellinfo_legal(int slab, uint64 slabsize) {
    // Return 1 if ok, return 0 if not.
    int cpd = PP->cpd;
    for (int j=0; j<cpd; j++)
	for (int k=0; k<cpd; k++) {
	    if (PP->CellInfo(slab,j,k)->legalvalue(slabsize)==0) {
		assertf(0==9, "Failed for slab %d, j=%d, k=%d\n", slab, j, k);
	        return 0;
	    }
	}
    return 1;
}
