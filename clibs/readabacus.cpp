// Copyright 2012-2025 The Abacus Developers
// SPDX-License-Identifier: GPL-3.0-or-later

#include "cell_header.h"
#include "packN_storage.cpp"
#include "threevector.hh"
#include "particle_subsample.cpp"
#include <omp.h>

#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <fcntl.h>
#include <string.h>
#define __STDC_FORMAT_MACROS
#include <inttypes.h>
#include <stdint.h>
#include <errno.h>
#include "iolib.cpp"

// Unpack a buffer of packN data
// The particles are written into `out`
// number of read particles is returned
template <int N, typename T>
uint64_t unpack_packN(packN<N> *data, size_t datasize, int nthread, int zspace, double subsample_frac, ThreeVector<T> *posout, ThreeVector<T> *velout, uint64 *pidout){
    if(subsample_frac <= 0)
        return 0;

    // Is reading the whole file then parsing it more efficient than doing fread from disk?
    // This way certainly uses more memory.
    if(datasize%N != 0){ //ensure the file is sensible
        fmt::print(stderr, "[Error] Datasize {:d} of file not divisible by {:d}.  Is this a pack{:d} file?\n", datasize, N, N);
        exit(1);
    }
    uint64_t max_NP = datasize/N;  // upper limit on the number of particles
    
    if (datasize == 0){
        fmt::print("[Warning] empty pack{:d} buffer encountered\n", N);
        return 0;
    }

    int do_subsample = subsample_frac > 0 && subsample_frac < 1;  // do we need to bother hashing?
    if(do_subsample){
        fmt::print(stderr, "subsampling not implemented with parallel rewrite\n");
        exit(1);
    }

    int return_pos = posout != NULL;
    int return_vel = velout != NULL;
    int return_pid = pidout != NULL;
    
    /* For parallel processing, we do one pass over the data to count particles and distribute work
    */
    
    uint64 readstart[nthread+1];  // where to start reading
    readstart[0] = 0;
    readstart[nthread] = max_NP;
    uint64 n_per_thread = max_NP/nthread;  // last thread will have fewer real particles
    uint64 writestart[nthread];  // where to start writing
    writestart[0] = 0;
    uint64 np_real = 0;
    
    if (nthread > 1){
        uint64 n_thisthread = 0;
        int tid = 1;
        for(uint64 j = 0; j < max_NP; j++){
            packN<N> p = data[j];

            if (p.iscell()){
                // Can only start on a cell
                // Keep counting if we're on the last thread so we can get np_real as a sanity check
                if(n_thisthread >= n_per_thread && tid != nthread){
                    writestart[tid] = n_thisthread + writestart[tid-1];
                    readstart[tid] = j;
                    np_real += n_thisthread;
                    n_thisthread=0;
                    tid++;
                }
            } else {
                n_thisthread++;
            }
        }
        if (tid != nthread){
            // TODO: will fail for very small L0 files, for example.
            fmt::print(stderr, "Didn't find divisions for all threads! tid={:d}, nthread={:d}\n", tid, nthread);
            exit(1);  // technically not an error, but probably indicates a bug
        }
        np_real += n_thisthread;  // get the last thread's count
        readstart[nthread] = max_NP;  // force the last thread to do all the rest
    }
    
    #pragma omp parallel num_threads(nthread)
    {   
        int tid = omp_get_thread_num();
        uint64 jstart = readstart[tid];  // where to start reading
        uint64 jstop = readstart[tid+1];
        uint64_t i = writestart[tid];  // where to start writing
        
        if(!data[jstart].iscell()){
            fmt::print(stderr, "Thread not starting on a cell!\n");
            exit(1);
        }
        
        cell_header current_cell;
        current_cell.vscale = 0;   // Make this illegal
        
        for(uint64 j = jstart; j < jstop; j++){
            packN<N> p = data[j];

            if (p.iscell()) {
                current_cell = p.unpack_cell();
                if(!current_cell.islegal()){
                    fmt::print(stderr, "Illegal pack14 cell encountered.\n");
                    exit(1);
                }
            } else {
                ThreeVector<T> pos;
                ThreeVector<T> vel;
                uint64_t id;
                p.unpack(pos, vel, id, current_cell);

                if(do_subsample && is_subsample_particle(id, subsample_frac, 0.) != 0)
                    continue;

                if(return_pos){
                    posout[i] = pos;
                    if(zspace)
                        posout[i][2] += vel[2];
                }

                if(return_vel){
                    velout[i] = vel;
                }

                if(return_pid){
                    pidout[i] = id;
                }
                i++;
            }
        }
        
        if (nthread > 1){
            // All threads except last should have written up to writestart[tid+1]
            if(tid < nthread-1){
                if (i != writestart[tid+1]){
                    fmt::print(stderr, "Wrong thread count! i={:d} != writestart[tid+1]={:d} for tid={:d}\n", i, writestart[tid+1], tid);
                    exit(1);
                }
            }
            else{
                if(i != np_real){
                    fmt::print(stderr, "Last thread wrong thread count! i={:d} vs np_real={:d}\n", i, np_real);
                    exit(1);
                }
            }
        }
        else {
            // Didn't do a counting pass
            np_real = i;
        }
    }
    
    return np_real;
}

// Traditional C++ "using" aliases didn't seem to work in the extern "C" block
#define ALIAS(NAME,N,T) \
uint64_t unpack_##NAME(packN<N> *data, size_t datasize, int nthread, int zspace, double subsample_frac, ThreeVector<T> *posout, ThreeVector<T> *velout, uint64 *pidout){ \
    return unpack_packN<N,T>(data, datasize, nthread, zspace, subsample_frac, posout, velout, pidout); \
}

extern "C" {
    ALIAS(pack14, 14,double)
    ALIAS(pack14f,14,float)
    ALIAS(pack9, 9,double)
    ALIAS(pack9f,9,float)
}
