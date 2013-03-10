/* phdata -- Skip the ParseHeader header of a file and print the rest to stdout.
If no filename is given, read from stdin.  
If many file names are given, perform on all of them.
*/

#include <stdio.h>
#include <assert.h>
#include <stdlib.h>

void print_data(FILE *fp) {
    int c, clast;
    if ((clast=getc(fp))==EOF) return;
    while ((c=getc(fp))!=EOF) {
        if (clast==''&&c=='') break;
        clast = c;
    }
    // Now we're found the end of the header; continue.
    while ((c=getc(fp))!=EOF) putchar(c);   
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
