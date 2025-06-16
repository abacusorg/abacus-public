// Copyright 2012-2025 The Abacus Developers
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once
// This is a simple routine to provide a consistent subsampling, based on a PID hash.
#include "MurmurHash3.cpp"

#define SUBSAMPLE_SEED 0xBEEF


#define SUBSAMPLE_A 1
#define SUBSAMPLE_B 2

inline int is_subsample_particle(const int64_t pid, const double subsample_fracA, const double subsample_fracB){
    uint32_t hash = 0;
    MurmurHash3_x86_32(&pid, sizeof(pid), SUBSAMPLE_SEED, &hash);
    double prob = (double) hash / UINT32_MAX;

    if (prob < subsample_fracA || subsample_fracA == 1)
    	return SUBSAMPLE_A; 

    else if ( (prob >= subsample_fracA && prob < (subsample_fracA + subsample_fracB)) || subsample_fracB == 1)
    	return SUBSAMPLE_B;

    else return 0; 
}
