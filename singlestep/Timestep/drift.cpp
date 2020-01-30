/* drift.cpp
 * This drifts the particles by a specified drift factor.  
 * It then checks to see if the particles have drifted outside of
 * their cell.  If so, particle is moved to the insert list.
 */

inline void DriftCell(Cell &c, FLOAT driftfactor) {
    int N = c.count();

    FLOAT3 *pos = (FLOAT3*) c.pos;
    FLOAT3 *vel = (FLOAT3*) c.vel;
    
    #pragma simd assert
    for (int b = 0; b<N; b++) {
        // Drift the position
        pos[b] += vel[b] * driftfactor;
    }
}

// Unlike DriftCell, this operator drifts the whole slab
void DriftSlab(int slab, FLOAT driftfactor){
    FLOAT* pos = (FLOAT *) SB->GetSlabPtr(PosSlab, slab);
    FLOAT* vel = (FLOAT *) SB->GetSlabPtr(VelSlab, slab);

    uint64 N = 3*SS->size(slab);

    #pragma omp parallel for schedule(static)
    //#pragma simd assert
    for(uint64 i = 0; i < N; i++){
        pos[i] += vel[i]*driftfactor;
    }
}


inline int RebinCell(Cell &c, int x, int y, int z) {
    // This should be ok regardless of box-center vs cell-center.
    // When a particle is found to be out of the cell, then it is
    // moved to the insert list immediately.  A particle from the
    // end is then moved in.

    #ifdef GLOBALPOS
        posstruct cellcenter = CP->LocalCellCenter(x,y,z);
    #endif
    FLOAT halfinvcpd = CP->halfinvcpd;
    int b = 0;
    int e = c.count();

    while(b<e) {
        // List is not yet done
        // Check if this particle is still in the cell
        posstruct residual = c.pos[b];
        #ifdef GLOBALPOS
            residual -= cellcenter;
        #endif
        if ( std::abs(residual.x) > halfinvcpd ||
             std::abs(residual.y) > halfinvcpd ||
             std::abs(residual.z) > halfinvcpd   ) {
            // This particle is outside the cell and must be rebinned.
            // Push it onto the insert list at the next location 
            IL->WrapAndPush( c.pos+b, c.vel+b, c.aux+b, x, y, z); 
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

// An alternate, theoretically more efficient implementation
// So far, it appears to be very slightly slower in practice
// TODO: neither approach seems really vector friendly.  We could do one sweep to establish high/low, and another to rebin via bitmask
inline int RebinCell2(Cell &c, int x, int y, int z) {
    FLOAT halfinvcpd = CP->halfinvcpd;
    int b = 0;
    int count = c.count();
    int e = count-1;

    while(b <= e) {
        // Check if this particle is still in the cell
        posstruct residual = c.pos[b];

        if ( std::abs(residual.x) > halfinvcpd ||
             std::abs(residual.y) > halfinvcpd ||
             std::abs(residual.z) > halfinvcpd   ) {
            // This particle is outside the cell and must be rebinned.
            // Push it onto the insert list at the next location 
            IL->WrapAndPush( c.pos+b, c.vel+b, c.aux+b, x, y, z);

            // Now find a low particle, starting from the end
            while(e > b){
                posstruct eres = c.pos[e];
                if ( std::abs(eres.x) > halfinvcpd ||
                     std::abs(eres.y) > halfinvcpd ||
                     std::abs(eres.z) > halfinvcpd   )
                    IL->WrapAndPush( c.pos+e, c.vel+e, c.aux+e, x, y, z);
                else
                    break;
                e--;
            }
            if(e <= b)
                break;
            c.pos[b] = c.pos[e];
            c.vel[b] = c.vel[e];
            c.aux[b] = c.aux[e];
            e--;
        }
        b++;
    }
    c.ci->active = b;
    return c.count() - b;
}

void DriftAndCopy2InsertList(int slab, FLOAT driftfactor, 
            void (*DriftCell)(Cell &c, FLOAT driftfactor)) {
    // Drift an entire slab
    STimer wc;
    PTimer move;
    PTimer rebin;

    wc.Start();
    
    int cpd = CP->cpd;

    uint64 ILbefore = IL->length;

    NUMA_FOR(y,0,cpd)
        for(int z=0;z<cpd;z++) {
            // We'll do the drifting and rebinning separately because
            // sometimes we'll want special rules for drifting.
            Cell c;
            c = CP->GetCell(slab ,y,z);
            move.Start();
            (*DriftCell)(c,driftfactor);
            move.Stop();
            rebin.Start();
            RebinCell(c, slab, y, z);
            rebin.Stop();
        }
    }
    wc.Stop();
    
    STDLOG(3, "Before collecting gaps, IL has length %d\n", IL->length);

    STDLOG(2,"Drifting slab %d has rebinned %d particles (%d - %d).\n",
        slab, IL->length-ILbefore, IL->length, ILbefore);
    
    // Compute timing by prorating the total wall-clock time 
    // by the time of the two major parts
    double seq = move.Elapsed() + rebin.Elapsed();
    double f_move = move.Elapsed()/seq;
    double f_rebin = rebin.Elapsed()/seq;

    struct timespec seq_move = scale_timer(f_move, wc.get_timer() );
    struct timespec seq_rebin = scale_timer(f_rebin, wc.get_timer() );
    
    DriftMove.increment(seq_move);
    DriftRebin.increment(seq_rebin);
}

// Do the drift step on the whole slab at once, instead of cell-by-cell
// The efficiency gain from doing a slab sweep appears to outweigh the
// cache efficiency loss from not immediately rebinning a drifted cell
// But this may not always be the case, so let's leave both versions for now
void DriftSlabAndCopy2InsertList(int slab, FLOAT driftfactor, void (*DriftSlab)(int slab, FLOAT driftfactor)) {
    STimer move;
    STimer rebin;
    
    int cpd = CP->cpd;

    uint64 ILbefore = IL->length;
    
    move.Start();
    (*DriftSlab)(slab, driftfactor);
    move.Stop();

    rebin.Start();
    //#pragma omp parallel for schedule(static)
    //for(int y=0;y<cpd;y++){
    NUMA_FOR(y,0,cpd)
        for(int z=0;z<cpd;z++) {
            // We'll do the drifting and rebinning separately because
            // sometimes we'll want special rules for drifting.
            Cell c = CP->GetCell(slab ,y,z);
            RebinCell(c, slab, y, z);
        }
    }
    rebin.Stop();

    STDLOG(2,"Drifting slab %d has rebinned %d particles (%d - %d).\n",
        slab, IL->length-ILbefore, IL->length, ILbefore);
    
    DriftMove.increment(move.timer);
    DriftRebin.increment(rebin.timer);
}
