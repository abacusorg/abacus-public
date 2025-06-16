/*
 * Copyright 2012-2025 The Abacus Developers
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

// This is a simple code to drive all of the utilities.
// Generically, we take many file names and apply in sequence.
// If the first argument is -, then we apply to stdin.
// All output should go to stdout to encourage piping.

int main(int argc, char *argv[]) {
    if (argc==1) {
	fprintf(stderr,"No input given to %s.  Use - for stdin.\n", PROGNAME);
        return 1;
    }
    if (argc==2&&strcmp(argv[1],"-")==0) print_data(stdin);
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
