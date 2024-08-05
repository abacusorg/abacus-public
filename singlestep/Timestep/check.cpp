/* check.cpp 
 * Future general validation code may be placed here.
 * Alternatively, this file could easily be refactored into the cellinfo class
 */

int are_cellinfo_legal(int slab, uint64 slabsize, uint64 slabsize_with_ghost) {
    // Return 1 if ok, return 0 if not.
    int cpd = CP->cpd;
    NUMA_FOR(j,0,cpd, NO_CLAUSE, FALLBACK_STATIC){
		for (int k = node_z_start_ghost; k < node_z_start_ghost + node_z_size_with_ghost; k++) {
			
			int isghost = (k < node_z_start_ghost + GHOST_RADIUS) || (k >= node_z_start_ghost + node_z_size_with_ghost - GHOST_RADIUS);
			
			assertf(CP->CellInfo(slab,j,k)->legalvalue(slabsize, slabsize_with_ghost, isghost),
				"Failed for slab %d, j=%d, k=%d\n", slab, j, k);
		}
	}
	NUMA_FOR_END;
    return 1;
}
