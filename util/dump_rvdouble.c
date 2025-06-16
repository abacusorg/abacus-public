/*
 * Copyright 2012-2025 The Abacus Developers
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#define PROGNAME "dump_rvdouble"

#include "util_header.h"
#include "util_main.c"


typedef struct RVdouble {
	double pos[3];          // Global position, unit box
	double vel[3];          // Zspace, unit box
} RVdouble;

int range(double a) { a=fabs(a); return (a!=0&&(a<1e-6||a>1000.0))?1:0; }
int RVrange(RVdouble *rv) {
    return range(rv->pos[0]) +range(rv->pos[1]) +range(rv->pos[2])
    	+range(rv->vel[0]) +range(rv->vel[1]) +range(rv->vel[2]);
}

void print_data(FILE *fp, const char *fn) {
	(void) fn;
    RVdouble rv;
    int count = 0;
    while (fread(&rv, sizeof(RVdouble), 1, fp)==1) {
        if (RVrange(&rv)) {
	    printf("%6d   %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e\n",
		count,
		rv.pos[0],
		rv.pos[1],
		rv.pos[2],
		rv.vel[0],
		rv.vel[1],
		rv.vel[2]
	    );
	} else {
	    printf("%6d   %10.7f %10.7f %10.7f %10.7f %10.7f %10.7f\n",
		count,
		rv.pos[0],
		rv.pos[1],
		rv.pos[2],
		rv.vel[0],
		rv.vel[1],
		rv.vel[2]
	    );
	}
	count++;
    }
}
