#define PROGNAME "dump_uint64"

#include "util_header.h"
#include "util_main.c"

#include <stdint.h>

#define uint64 uint64_t

void print_data(FILE *fp, const char *fn [[maybe_unused]]) {
    uint64 rv;
    int count = 0;
    while (fread(&rv, sizeof(uint64), 1, fp)==1) {
	printf("%6d   0x%08lx %08lx  %15ld\n", count, rv>>32, rv&0xffff, rv);
	count++;
    }
}
