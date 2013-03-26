/* phdata -- Skip the ParseHeader header of a file and print the rest to stdout.
If no filename is given, read from stdin.  
If many file names are given, perform on all of them.
*/

#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <math.h>

typedef struct float3 {
	float pos[3];          
} float3;

int range(float a) { a=fabs(a); return (a!=0&&(a<1e-6||a>1000.0))?1:0; }
int RVrange(float3 *rv) {
    return range(rv->pos[0]) +range(rv->pos[1]) +range(rv->pos[2]);
}

void print_data(FILE *fp) {
    float3 rv;
    int count = 0;
    while (fread(&rv, sizeof(float3), 1, fp)==1) {
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

int main(int argc, char *argv[]) {
    if (argc==1) print_data(stdin);
    else {
	int f;
        for (f=1; f<argc; f++) {
	    FILE *fp;
	    fp = fopen(argv[f],"r");
	    if (fp==NULL) {
	        fprintf(stderr, "File %s not found or cannot be opened.\n", argv[f]);
		exit(1);
	    } else {
	        // fprintf(stderr, "Opened file %s.\n", argv[f]);
	    }
	    print_data(fp);
	    fclose(fp);
	}
    }
    return 0;
}
