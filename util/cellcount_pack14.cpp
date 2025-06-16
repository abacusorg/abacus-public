// Copyright 2012-2025 The Abacus Developers
// SPDX-License-Identifier: GPL-3.0-or-later

#define PROGNAME "cellcount_pack14"

// This simply outputs the counts in each cell.
// The cell number is provided in a comment line, but this does add a lot 
// of space.  
// Doing 'grep -v #' would remove those, which might then allow a 
// simple pipe to some kind of 2d plotting utility.

#include "util_header.h"
#include "util_main.c"

#include "../include/cell_header.h"
#include "../include/packN_storage.cpp"

typedef struct RVdouble {
	double pos[3];          // Global position, unit box
	double vel[3];          // Zspace, unit box
} RVdouble;

int range(double a) { a=fabs(a); return (a!=0&&(a<1e-6||a>1000.0))?1:0; }
int RVrange(RVdouble *rv) {
    return range(rv->pos[0]) +range(rv->pos[1]) +range(rv->pos[2])
    	+range(rv->vel[0]) +range(rv->vel[1]) +range(rv->vel[2]);
}

long long int count;
long long int cellcount;

// Nothing to do at the beginning or end
void start() { count = 0; cellcount = -1; return; }
void end() { 
    printf("%lld\n", cellcount);
    printf("# Total count: %lld\n", count); 
}

void print_data(FILE *fp) {
    cell_header current_cell;
    pack14 p;
    current_cell.vscale = 0;   // Make this illegal

    RVdouble rv;
    uint64_t id;

    while (fread(&p, sizeof(pack14), 1, fp)==1) {
	if (p.iscell()) {
	    current_cell = p.unpack_cell();
	    if (cellcount>=0) 
	    	printf("%lld\n", cellcount);
	    printf("# Cell %d %d %d\n", current_cell.i, current_cell.j, current_cell.k);
	    cellcount=0;
	} else {
	    assert(current_cell.islegal());
	    p.unpack(rv.pos, rv.vel, &id, current_cell);
	    cellcount++;
	    count++;
	}
    }
}
