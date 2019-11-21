// This takes the input and outputs binary RVfloat format.

#define PROGNAME "pack14_to_rvfloat"

#include "util_header.h"
#include "util_main.c"

#include "../include/cell_header.h"
#include "../include/packN_storage.cpp"

typedef struct RVfloat {
	float pos[3];          // Global position, unit box
	float vel[3];          // Zspace, unit box
} RVfloat;

int range(float a) { a=fabs(a); return (a!=0&&(a<1e-6||a>1000.0))?1:0; }
int RVrange(RVfloat *rv) {
    return range(rv->pos[0]) +range(rv->pos[1]) +range(rv->pos[2])
    	+range(rv->vel[0]) +range(rv->vel[1]) +range(rv->vel[2]);
}

void print_data(FILE *fp) {
    cell_header current_cell;
    pack14 p;
    current_cell.vscale = 0;   // Make this illegal

    RVfloat rv;
    uint64_t id;

    int count = 0;
    while (fread(&p, sizeof(pack14), 1, fp)==1) {
	if (p.iscell()) {
	    current_cell = p.unpack_cell();
	    // printf("# Cell %d %d %d\n", current_cell.i, current_cell.j, current_cell.k);
	} else {
	    assert(current_cell.islegal());
	    p.unpack(rv.pos, rv.vel, &id, current_cell);
	    fwrite(&rv, sizeof(RVfloat), 1, stdout);
	    count++;
	}
    }
}
