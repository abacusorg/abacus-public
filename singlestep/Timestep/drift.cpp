// Copyright 2012-2025 The Abacus Developers
// SPDX-License-Identifier: GPL-3.0-or-later

/* drift.cpp
 * This drifts the particles by a specified drift factor.  
 * It then checks to see if the particles have drifted outside of
 * their cell.  If so, particle is moved to the insert list.
 */

inline void DriftCell(Cell &c, FLOAT driftfactor) {
    int N = c.count();

    // __restrict__: tell C++ that despite these pointers
    // having the same type, they will never overlap.
    // #pragma omp simd appears to accomplish the same
    // thing in this case (probably more important
    // for function parameters).

    FLOAT3 * __restrict__ pos = (FLOAT3*) c.pos;
    FLOAT3 * __restrict__ vel = (FLOAT3*) c.vel;
    
    #pragma omp simd
    for (int b = 0; b<N; b++) {
        // Drift the position
        pos[b] += vel[b] * driftfactor;
    }
}

// Unlike DriftCell, this operator drifts a whole pencil at a time
void DriftPencil(int slab, int j, FLOAT driftfactor){
    FLOAT* __restrict__ pos = (FLOAT *) CP->PosCell(slab, j, node_z_start);
    FLOAT* __restrict__ vel = (FLOAT *) CP->VelCell(slab, j, node_z_start);

    // no ghost drift (ghost would use PencilLenWithGhost)
    uint64 Npen = 3*CP->PencilLen(slab, j);

    #pragma omp simd
    for(uint64 i = 0; i < Npen; i++){
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

/* Push an entire cell to the IL, marking all its particles as inactive.
 * This is used in the 2D code to copy boundary cells to the IL
 * so they can be sent to the neighbor.
 * Particles forced onto the IL this way might come back during the merge.
 */
void PushCellToIL(Cell &c, int x, int y, int z){
    uint64 count = c.ci->count;

    posstruct *p = c.pos;
    velstruct *v = c.vel;
    auxstruct *a = c.aux;
    // Note we're wrapping many particles that don't need wrapping
    // But we skipped the cell partition, so we don't know which is which
    for(uint64 i = 0; i < count; i++){
        IL->WrapAndPush(p+i, v+i, a+i, x, y, z);
    }

    // Mark all particles gone
    c.ci->active = 0;
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

    NUMA_FOR(y,0,cpd, NO_CLAUSE, FALLBACK_DYNAMIC){
        for(int z = node_z_start; z < node_z_start + node_z_size; z++) {
            // We'll do the drifting and rebinning separately because
            // sometimes we'll want special rules for drifting.
            Cell c = CP->GetCell(slab,y,z);
            move.Start();
            (*DriftCell)(c,driftfactor);
            move.Stop();
            rebin.Start();
            if( (z - node_z_start < MERGE_GHOST_RADIUS) || 
                ((node_z_start + node_z_size - z - 1) < MERGE_GHOST_RADIUS)){
                PushCellToIL(c, slab, y, z);  // near the ghost boundary; move whole cell to IL
            } else {
                RebinCell(c, slab, y, z);
            }
            rebin.Stop();
        }
    }
    NUMA_FOR_END;
    wc.Stop();

    STDLOG(2,"Drifting slab {:d} has rebinned {:d} particles ({:d} - {:d}).\n",
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

/* Do the drift step pencil-by-pencil, instead of cell-by-cell.
 * The efficiency gain of not stopping at each cell boundary
 * appears to outweigh the cache efficiency loss from not immediately
 * rebinning a drifted cell.
 * But this may not always be the case, so let's leave both versions for now
 */
void DriftPencilsAndCopy2InsertList(int slab, FLOAT driftfactor,
    void (*DriftPencil)(int slab, int y, FLOAT driftfactor)
    ) {
    STimer move;
    STimer rebin;
    
    int cpd = CP->cpd;

    uint64 ILbefore = IL->length;

    // We'll do the drifting and rebinning separately because
    // sometimes we'll want special rules for drifting.
    
    move.Start();
    NUMA_FOR(y,0,cpd, NO_CLAUSE, FALLBACK_DYNAMIC){
        (*DriftPencil)(slab, y, driftfactor);
    }
    NUMA_FOR_END;
    move.Stop();

    rebin.Start();
    
    NUMA_FOR(y,0,cpd, NO_CLAUSE, FALLBACK_DYNAMIC){
        // primary cells within MERGE_GHOST_RADIUS of the edge get pushed entirely to the IL
        for(int z = node_z_start; z < node_z_start + MERGE_GHOST_RADIUS; z++) {
            Cell c = CP->GetCell(slab,y,z);
            PushCellToIL(c, slab, y, z);
        }

        // then the middle
        for(int z = node_z_start + MERGE_GHOST_RADIUS; z < node_z_start + node_z_size - MERGE_GHOST_RADIUS; z++) {
            Cell c = CP->GetCell(slab,y,z);
            RebinCell(c, slab, y, z);
        }

        // then the primary cells at far edge
        for(int z = node_z_start + node_z_size - MERGE_GHOST_RADIUS; z < node_z_start + node_z_size; z++) {
            Cell c = CP->GetCell(slab,y,z);
            PushCellToIL(c, slab, y, z);
        }
    }
    NUMA_FOR_END;
    rebin.Stop();

    STDLOG(2,"Drifting slab {:d} has rebinned {:d} particles ({:d} - {:d}).\n",
        slab, IL->length-ILbefore, IL->length, ILbefore);
    
    DriftMove.increment(move.timer);
    DriftRebin.increment(rebin.timer);
}
