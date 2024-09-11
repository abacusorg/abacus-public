/* phheader -- Print the ParseHeader header of a file to stdout.
This routine doesn't do any parsing -- just echos the text, 
stripping out the final demarker (^B\n).  We add a \n to the end of the file.
If no filename is given, read from stdin.  
If many file names are given, perform on all of them.
*/

#define PROGNAME "phheader"

#include "util_header.h"
#include "util_main.c"


void print_data(FILE *fp, const char *fn [[maybe_unused]]) {
    int c, clast;
    if ((clast=getc(fp))==EOF) return;
    while ((c=getc(fp))!=EOF) {
        if (clast==''&&c=='\n') { putchar('\n'); return; }
        putchar(clast); clast = c;
    }
    putchar(clast);   // If the demarker hasn't been found, we still have one more char
    putchar('\n');    // Add a \n as a formatting guard.
}
