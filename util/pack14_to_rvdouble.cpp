// This takes the input and outputs binary RVdouble format.

#define PROGNAME "pack14_to_rvdouble"

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

void print_data(FILE *fp, const char *fn [[maybe_unused]]) {
    cell_header current_cell;
    pack14 p;
    current_cell.vscale = 0;   // Make this illegal

    RVdouble rv;
    uint64_t id;

    int count = 0;
    while (fread(&p, sizeof(pack14), 1, fp)==1) {
	if (p.iscell()) {
	    current_cell = p.unpack_cell();
	    // printf("# Cell %d %d %d\n", current_cell.i, current_cell.j, current_cell.k);
	} else {
	    assert(current_cell.islegal());
	    p.unpack(rv.pos, rv.vel, &id, current_cell);
	    fwrite(&rv, sizeof(RVdouble), 1, stdout);
	    count++;
	}
    }
}
