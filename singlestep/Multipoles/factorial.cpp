// Copyright 2012-2025 The Abacus Developers
// SPDX-License-Identifier: GPL-3.0-or-later

double fact(int n) {
    double f = 1.0;
    for(int i=2;i<=n;i++) f *= i;
    return f;
}

double fact2(int n) {
    double f;
    int i;

    f = 1.0;
    i = n;
    while(i>0) {
        f = f*i;
        i -= 2;
    }
    return f;
}



