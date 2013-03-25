/* dump_uint64 -- Convert binary to ASCII
*/

#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>

#define uint64 uint64_t

void print_data(FILE *fp) {
    uint64 rv;
    int count = 0;
    while (fread(&rv, sizeof(uint64), 1, fp)==1) {
	printf("%6d   0x%08lx %08lx  %15ld\n", count, rv>>32, rv&0xffff, rv);
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
