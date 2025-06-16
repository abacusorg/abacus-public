// Copyright 2012-2025 The Abacus Developers
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef INCLUDE_FACTORIAL
#define INCLUDE_FACTORIAL

double fact(int n) {
    uint64 f = 1;  // Do this as an int to get GCC vectorization
    for(int i=2;i<=n;i++) f *= i;
    return (double) f;
}

double fact2(int n) {
    uint64 f;
    
    f = 1;
    for(int i = n; i > 0; i -= 2)
        f *= i;
    
    return (double) f;
}



#endif // INCLUDE_FACTORIAL
