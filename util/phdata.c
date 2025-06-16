/*
 * Copyright 2012-2025 The Abacus Developers
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/* phdata -- Skip the ParseHeader header of a file and print the rest to stdout.
If no filename is given, read from stdin.  
If many file names are given, perform on all of them.
*/

#define PROGNAME "phdata"
#include "util_header.h"
#include "util_main.c"


void print_data(FILE *fp, const char *fn) {
    (void) fn;
    int c, clast;
    if ((clast=getc(fp))==EOF) return;
    while ((c=getc(fp))!=EOF) {
        if (clast==''&&c=='\n') break;
        clast = c;
    }
    // Now we're found the end of the header; continue.
    while ((c=getc(fp))!=EOF) putchar(c);   
}
