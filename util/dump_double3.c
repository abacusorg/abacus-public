#define PROGNAME "dump_double3"
#include "util_header.h"
#include "util_main.c"

typedef struct double3 {
	double pos[3];          
} double3;

int range(double a) { a=fabs(a); return (a!=0&&(a<1e-6||a>1000.0))?1:0; }
int RVrange(double3 *rv) {
    return range(rv->pos[0]) +range(rv->pos[1]) +range(rv->pos[2]);
}

void print_data(FILE *fp, const char *fn) {
    double3 rv;
    int count = 0;
    while (fread(&rv, sizeof(double3), 1, fp)==1) {
        if (RVrange(&rv)) {
	    printf("%6d   %10.3e %10.3e %10.3e\n",
		count,
		rv.pos[0],
		rv.pos[1],
		rv.pos[2]
	    );
	} else {
	    printf("%6d   %10.7f %10.7f %10.7f\n",
		count,
		rv.pos[0],
		rv.pos[1],
		rv.pos[2]
	    );
	}
	count++;
    }
}
