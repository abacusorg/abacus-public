#pragma once
// This is a simple routine to provide a consistent subsampling, based on a PID hash.
#include "MurmurHash3.cpp"

#define SUBSAMPLE_SEED 0xBEEF
inline int is_subsample_particle(const int64_t pid, const double subsample_frac){
    uint32_t hash = 0;
    MurmurHash3_x86_32(&pid, sizeof(pid), SUBSAMPLE_SEED, &hash);
    double prob = (double) hash / UINT32_MAX;
    return prob < subsample_frac;
}
