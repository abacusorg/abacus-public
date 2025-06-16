/*
 * Copyright 2012-2025 The Abacus Developers
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

template <typename U>
    struct jpstruct {
        U x, y, z, m;
    };

    template <typename U, int size>
    struct ipstruct {
        U x[size], y[size], z[size], eps2[size];
    };

    template <typename U,int size>
    struct apstruct {
        U x[size], y[size], z[size], pot[size];
    };

