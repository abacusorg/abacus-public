/*
 * Copyright 2012-2025 The Abacus Developers
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/*

In the MPI parallel code, some of the MPI work falls way behind the theoretical
network speed on platforms like Summit.  The reason is unclear, but we worry
that when we get into this state we are querying the MPI progress thread via
MPI_Test so often that it may interrupt other MPI work.  MpiP profiling
showed that 1/3rd of our time spent in MPI calls goes to MPI_Test.  The
purpose of this class is to provide a global rate limiter on how often
we can call MPI_Test.

*/
#pragma once

#include <time.h>

#ifndef NSEC_PER_SEC
#define NSEC_PER_SEC 1000000000
#endif

class AbacusMPILimiter {
    struct timespec last;
    int64_t delta_ns;

public:

    // Constructor, with given number of milliseconds between calls
    AbacusMPILimiter(double delta_ms){
        delta_ns = (int64_t) (delta_ms*1e6);

        assert(clock_gettime(CLOCK_MONOTONIC, &last) == 0);
    }

    // Returns 1 if enough time has passed since the last time returning 1; else 0.
    int Try(){
        if(delta_ns <= 0)
            return 1;

        struct timespec now;
        assert(clock_gettime(CLOCK_MONOTONIC, &now) == 0);

        int64_t thisdelta;
        thisdelta = (int64_t)(now.tv_sec - last.tv_sec)*NSEC_PER_SEC + (now.tv_nsec - last.tv_nsec);

        if(thisdelta >= delta_ns){
            last = now;
            return 1;
        }

        return 0;
    }
};
