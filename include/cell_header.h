/* cell_header.h

This is a simple class to track information about cells for the append_arena.cpp
and pack14.cpp functions.
*/

#ifndef CELL_HEADER_CPP
#define CELL_HEADER_CPP

#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include "threevector.hh"

#define SHORT int16_t     // We want this to be 16-bits, must be signed.
// Might do a test to ensure that sizeof(char)=1.

class cell_header {
public:
    SHORT vscale, cpd;
    SHORT i,j,k;

    cell_header() { vscale = 0; cpd = 0; i = 0; j = 0; k = 0; }
    cell_header(SHORT ijk[3], SHORT _cpd, SHORT _vscale) { 
    	cpd = _cpd; vscale = _vscale; i = ijk[0]; j = ijk[1]; k = ijk[2];
    }
    cell_header(int ijk[3], int _cpd, int _vscale) { 
    	cpd = _cpd; vscale = _vscale; i = ijk[0]; j = ijk[1]; k = ijk[2];
    }
    ~cell_header(void) { }

    int islegal() { if (vscale>0) return 1; else return 0; }
    
    // Some simple unwrapping if one wants.
    void cellid(int cell[3]) { cell[0] = i; cell[1] = j; cell[2] = k; }

#ifdef __THREEVECTOR_CC__
    cell_header(integer3 ijk, int _cpd, int _vscale) { 
    	cpd = _cpd; vscale = _vscale; i = ijk.x; j = ijk.y; k = ijk.z;
    }
    integer3 cellid() {
        integer3 ijk(i,j,k); return ijk;
    }
#endif
};

#endif // CELL_HEADER_CPP
