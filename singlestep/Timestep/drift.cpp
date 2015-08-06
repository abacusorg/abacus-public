/* drift.cpp
This drifts the particles by a specified drift factor.  
It then checks to see if the particles have drifted outside of
their cell.  This is done in two steps: 
1) swap moving particles to the end of their cell list (multi-threaded)
2) add those particles to the insert list (single-threaded)
*/

void DriftCell(Cell c, FLOAT driftfactor) {
    int N = c.count();

	//#pragma simd
    for (int b = 0; b<N; b++) {
        // Drift the position
        c.pos[b] += c.vel[b] * driftfactor;
    }
}


void RebinCell(Cell c, int x, int y, int z) {
    // This should be ok regardless of box-center vs cell-center.
    posstruct cellcenter = PP->LocalCellCenter(x,y,z);
    int b = 0;
    int e = c.count();		

    while(b<e) {
        // Check if this particle is still in the cell
        posstruct residual = c.pos[b] - cellcenter;
        if (   fabs(residual.x)>PP->halfinvcpd
                || fabs(residual.y)>PP->halfinvcpd
                || fabs(residual.z)>PP->halfinvcpd) {
            // We will need to rebin this particle
            c.swap(b, e-1);
            e--;
        }
        else 
        b++;
    }
    
    c.ci->active = e;
}

void DriftAndCopy2InsertList(int slab, FLOAT driftfactor, 
void (*DriftCell)(Cell c, FLOAT driftfactor)) {
    int cpd = PP->cpd;
    DriftMoveRebin.Start();
    #pragma omp parallel for schedule(dynamic,1) 
    for(int y=0;y<cpd;y++)
        for(int z=0;z<cpd;z++) {
            // We'll do the drifting and rebinning separately because
            // sometimes we'll want special rules for drifting.
            DriftMove.Start();
            Cell c = PP->GetCell(slab ,y,z);
            (*DriftCell)(c,driftfactor);
            DriftMove.Stop();
            DriftRebin.Start();
            RebinCell(c, slab, y, z);
            DriftRebin.Stop();
        }
    DriftMoveRebin.Stop();

    // below has to be done serially

    DriftInsert.Start();
    uint64 ILbefore = IL->length;
    for(int y=0;y<cpd;y++) 
        for(int z=0;z<cpd;z++) {
            Cell c = PP->GetCell(slab,y,z);
            for(int q=c.active(); q<c.count(); q++) {
                IL->WrapAndPush( c.pos+q, c.vel+q, c.aux+q, slab,y,z);
            }
        }
    DriftInsert.Stop();
    STDLOG(1,"Drifting slab %d has rebinned %d particles.\n",
    slab, IL->length-ILbefore);
}
