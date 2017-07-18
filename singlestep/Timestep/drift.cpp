/* drift.cpp
This drifts the particles by a specified drift factor.  
It then checks to see if the particles have drifted outside of
their cell.  This is done in two steps: 
1) swap moving particles to the end of their cell list (multi-threaded)
2) add those particles to the insert list (single-threaded)
*/

void DriftCell(Cell c, FLOAT driftfactor) {
    int N = c.count();

    FLOAT3 *pos = (FLOAT3*) c.pos;
    FLOAT3 *vel = (FLOAT3*) c.vel;
    auxstruct *aux = c.aux;
    
    #pragma simd assert
    for (int b = 0; b<N; b++) {
        // Drift the position
        pos[b] += vel[b] * driftfactor;
        
        // set unused bit 48 to zero to force loading.
        // hopefully this makes rebinning faster?
        aux[b].aux &= ~((uint64) 1 << 48);
    }
}


int RebinCell(Cell c, int x, int y, int z) {
    // This should be ok regardless of box-center vs cell-center.
    #ifdef GLOBALPOS
    posstruct cellcenter = PP->LocalCellCenter(x,y,z);
    #endif
    int b = 0;
    int e = c.count();

    FLOAT halfinvcpd = PP->halfinvcpd;
    while(b<e) {
        // Check if this particle is still in the cell
        posstruct residual = c.pos[b];
        #ifdef GLOBALPOS
        residual -= cellcenter;
        #endif
        if (    std::abs(residual.x) > halfinvcpd
             || std::abs(residual.y) > halfinvcpd
             || std::abs(residual.z) > halfinvcpd) {
            // We will need to rebin this particle
            
            // It's possible the particle at the end of the list needs to be rebinned too
            // Keep it there if so (i.e don't swap it)
            posstruct eresidual;
            do {
                e--;
                eresidual = c.pos[e];
                // if (e <= b) break outer;  // faster to do here or in loop condition?
            } while(e > b && (   std::abs(eresidual.x) > halfinvcpd
                              || std::abs(eresidual.y) > halfinvcpd
                              || std::abs(eresidual.z) > halfinvcpd));
            c.swap(b, e);
        }
        // We know the particle we swapped with does not need rebinning,
        // so don't check it again
        b++;
    }
    
    c.ci->active = e;
    return c.count() - e;
}

void DriftAndCopy2InsertList(int slab, FLOAT driftfactor, 
void (*DriftCell)(Cell c, FLOAT driftfactor)) {
    STimer wc;
    PTimer move;
    PTimer rebin;
    
    int cpd = PP->cpd;
    
    // Number of rebinned particles in a given z-pencil
    uint32_t row_rebinned[cpd];  // false sharing hazard?
    uint64_t n_rebinned = 0;
    
    wc.Start();
    #pragma omp parallel for schedule(static)
    for(int y=0;y<cpd;y++){
        row_rebinned[y] = 0;
        for(int z=0;z<cpd;z++) {
            // We'll do the drifting and rebinning separately because
            // sometimes we'll want special rules for drifting.
            Cell c = PP->GetCell(slab ,y,z);
            move.Start();
            (*DriftCell)(c,driftfactor);
            move.Stop();
            rebin.Start();
            row_rebinned[y] += RebinCell(c, slab, y, z);
            rebin.Stop();
        }
    }
    wc.Stop();
    
    uint64_t row_offset[cpd];
    for(int y = 0; y < cpd; y++){
        row_offset[y] = n_rebinned;
        n_rebinned += row_rebinned[y];
    }
    
    double seq = move.Elapsed() + rebin.Elapsed();
    double f_move = move.Elapsed()/seq;
    double f_rebin = rebin.Elapsed()/seq;

    struct timeval seq_move = scale_timer(f_move, wc.get_timer() );
    struct timeval seq_rebin = scale_timer(f_rebin, wc.get_timer() );
    
    DriftMove.increment(seq_move);
    DriftRebin.increment(seq_rebin);

    DriftInsert.Start();
    
    // Grow the insert list by the total number of rebinned particles
    uint64 ILbefore = IL->length;
    IL->length += n_rebinned;
    
    // Now let each thread write into the insert list, starting at its row offset
    #pragma omp parallel for schedule(static)
    for(int y=0;y<cpd;y++){
        uint32_t this_row_rebinned = 0;
        uint64_t this_row_offset = row_offset[y];
        for(int z=0;z<cpd;z++) {
            Cell c = PP->GetCell(slab,y,z);
            int e = c.count();
            for(int q=c.active(); q<e; q++) {
                IL->WrapAndPush( c.pos+q, c.vel+q, c.aux+q, slab, y, z, ILbefore + this_row_offset + this_row_rebinned);
                this_row_rebinned++;
            }
        }
        assert(this_row_rebinned == row_rebinned[y]);
    }
    DriftInsert.Stop();
    STDLOG(1,"Drifting slab %d has rebinned %d particles.\n",
    slab, IL->length-ILbefore);
}
