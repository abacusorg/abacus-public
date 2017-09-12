/* drift.cpp
This drifts the particles by a specified drift factor.  
It then checks to see if the particles have drifted outside of
their cell.  If so, particle is moved to the insert list.
*/

void DriftCell(Cell &c, FLOAT driftfactor) {
    int N = c.count();

    FLOAT3 *pos = (FLOAT3*) c.pos;
    FLOAT3 *vel = (FLOAT3*) c.vel;
    
    #pragma simd assert
    for (int b = 0; b<N; b++) {
        // Drift the position
        pos[b] += vel[b] * driftfactor;
    }
}


int RebinCell(Cell &c, int x, int y, int z, ilgap *gap) {
    // This should be ok regardless of box-center vs cell-center.
    // When a particle is found to be out of the cell, then it is
    // moved to the insert list immediately.  A particle from the
    // end is then moved in.

    #ifdef GLOBALPOS
        posstruct cellcenter = PP->LocalCellCenter(x,y,z);
    #endif
    FLOAT halfinvcpd = PP->halfinvcpd;
    int b = 0;
    int e = c.count();

    while(b<e) {
        // List is not yet done
        // Check if this particle is still in the cell
        posstruct residual = c.pos[b];
        #ifdef GLOBALPOS
            residual -= cellcenter;
        #endif
        if (    std::abs(residual.x) > halfinvcpd
             || std::abs(residual.y) > halfinvcpd
             || std::abs(residual.z) > halfinvcpd) {
            // This particle is outside the cell and must be rebinned.
            // Push it onto the insert list at the location indicated 
            // by the ilgap controller.
            IL->WrapAndPush( c.pos+b, c.vel+b, c.aux+b, x, y, z, 
                    gap->nextid());
            // Now move the end particle into this location
            e--;    // Don't increment b, as we'll need to try again
            c.pos[b] = c.pos[e];
            c.vel[b] = c.vel[e];
            c.aux[b] = c.aux[e];
        } else b++;    // Particle was ok; keep moving along
    }
    c.ci->active = e;
    return c.count() - e;
}

void DriftAndCopy2InsertList(int slab, FLOAT driftfactor, 
            void (*DriftCell)(Cell &c, FLOAT driftfactor)) {
    // Drift an entire slab
    STimer wc;
    PTimer move;
    PTimer rebin;
    
    int cpd = PP->cpd;

    int nthreads = omp_get_max_threads();
    uint64 ILbefore = IL->length;
    ilgap *gaps = new ilgap[nthreads];
        // This reserves space on the insert list for each thread to use.
        // This class silently gets new space as needed.
    
    wc.Start();
    #pragma omp parallel for schedule(static)
    for(int y=0;y<cpd;y++){
        ilgap *thisgap = gaps + omp_get_thread_num();
            // The space pointer for this list
        for(int z=0;z<cpd;z++) {
            // We'll do the drifting and rebinning separately because
            // sometimes we'll want special rules for drifting.
            Cell c = PP->GetCell(slab ,y,z);
            move.Start();
            (*DriftCell)(c,driftfactor);
            move.Stop();
            rebin.Start();
            RebinCell(c, slab, y, z, thisgap);
            rebin.Stop();
        }
    }
    wc.Stop();

    DriftInsert.Start();
    IL->CollectGaps(gaps, nthreads);
    delete[] gaps;
    DriftInsert.Stop();
    STDLOG(1,"Drifting slab %d has rebinned %u particles.\n",
        slab, (uint32_t)(IL->length-ILbefore));
    
    // Compute timing by prorating the total wall-clock time 
    // by the time of the two major parts
    double seq = move.Elapsed() + rebin.Elapsed();
    double f_move = move.Elapsed()/seq;
    double f_rebin = rebin.Elapsed()/seq;

    struct timeval seq_move = scale_timer(f_move, wc.get_timer() );
    struct timeval seq_rebin = scale_timer(f_rebin, wc.get_timer() );
    
    DriftMove.increment(seq_move);
    DriftRebin.increment(seq_rebin);
    return;
}


