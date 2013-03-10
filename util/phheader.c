/* phheader -- Print the ParseHeader header of a file to stdout.
This routine doesn't do any parsing -- just echos the text, 
stripping out the final demarker (^B^B).  We add a \n to the end of the file.
If no filename is given, read from stdin.  
If many file names are given, perform on all of them.
*/

#include <stdio.h>
#include <assert.h>
#include <stdlib.h>

void print_header(FILE *fp) {
    int c, clast;
    if ((clast=getc(fp))==EOF) return;
    while ((c=getc(fp))!=EOF) {
        if (clast==''&&c=='') return;
        putchar(clast); clast = c;
    }
    putchar(clast);   // If the demarker hasn't been found, we still have one more char
    putchar('\n');    // Add a \n as a formatting guard.
}

int main(int argc, char *argv[]) {
    if (argc==1) print_header(stdin);
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
	    print_header(fp);
	    fclose(fp);
	}
    }
    return 0;
}
