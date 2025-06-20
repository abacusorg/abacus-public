// Copyright 2012-2025 The Abacus Developers
// SPDX-License-Identifier: GPL-3.0-or-later

/* 
Compile with:
g++ -O2 -DPPTEST -march=native -fopenmp -lgomp -DOMP parallel.partition.cpp

This is a multi-threaded partition algorithm.  Given a list, a
boolean selection function, and a value, all selected objects
are moved to the end.

For Nt threads, divide the list into Nt roughly equal pieces, 
recording the beginning and ending point.

Then have each thread partition one piece, recording the midpoint.  
When done, the low objects for piece j will be in begin[j]..mid[j]-1 and
the high objects will be in mid[j]..begin[j+1]-1.  This partitioning
is the usual sweep from both ends, swapping out-of-place pairs.

Now in a serial section, we will develop a work plan for block swaps.
Each block swap will take a section of the high objects and move them
to a section of low objects.  The objects will not be reversed in 
their detailed order, but we will match up the first N high particles
to the last N low particles, etc.  

We first make a plan to create the blocks.  Then we do a second 
sweep to divide some blocks into multiple blocks and assign a balanced
amount of work to each thread.

Then we give the threads the job to carry out their lists of block swap.  

TODO: I'm unsure whether the boolean function that is passed in will
in fact be compiled in-line.  If not, that could really slow down the 
initial passes and we'd be better off hard-coding equality.
*/

// #define PPTEST     // Only if you want the unit-test (and verbose) executable

#ifdef PPTEST
#include <stdlib.h>
#include <vector>
#include <assert.h>
#include <stdio.h>
#include <omp.h>

#include <inttypes.h>
#ifndef uint64
#define uint64 uint64_t
#endif

#endif

class SwapBlock {
  public:
    uint64 highstart;
    uint64 lowstart;
    uint64 nswap;
    SwapBlock(SwapBlock *s) {
         highstart = s->highstart;
         lowstart = s->lowstart;
         nswap = s->nswap;
    }
    SwapBlock(uint64 _high, uint64 _low, uint64 _nswap) {
         highstart = _high;
         lowstart = _low;
         nswap = _nswap;
    }
    ~SwapBlock() { }

    template <typename ListType>
    void swap(ListType *list) {
        for (uint64 j=0; j<nswap; j++) std::swap(list[highstart+j], list[lowstart+j]);
    }
};

class PartitionBlock {
  public:
    uint64 begin, end;  // We work on the list from [begin, end-1]
            // There are end-begin points
    uint64 mid;     // The low points are [begin, mid-1], high are [mid, end-1]
            // There are end-mid high points, mid-begin low points
    std::vector<SwapBlock> threadplan;
        // This is the work for this thread

    PartitionBlock() { }
    ~PartitionBlock() { }

    uint64 numhigh() { return end-mid; }
    uint64 numlow() { return mid-begin; }

    template <typename ListType, typename... Args>
    void Partition(uint64 _begin, uint64 _end, ListType *list, 
            bool (is_high)(ListType *obj, Args... args), Args... args) {
        assert(_begin < _end);
        begin=_begin; end=_end;
        // Then do the partitioning, loading up mid
        uint64 a = begin, b = end-1;
        while (1) {
            while (!is_high(list+a,args...)&&a<b) a++;    // So list[a] is high
            while (is_high(list+b,args...)&&a<b) b--;    // So list[b] is low
            // If the boolean function call is slow and we want to code equality in-line.
            // Tests suggest the call is inlined, though.
            // while (list[a]!=slab&&a<b) a++;    // So list[a]==slab
            // while (list[b]==slab&&a<b) b--;    // So list[b]!=slab
            if (a==b) break;
            std::swap(list[a], list[b]);
        }
        // We have a==b when we exit, but we don't yet know if this value is low or high
        if (is_high(list+a,args...)) mid = a; else mid=a+1;
        return;
    }

    template <typename ListType>
    void swap(ListType *list) {
        for (uint64 p=0; p<threadplan.size(); p++) threadplan[p].swap(list);
        return;
    }
};


template <typename ListType, typename... Args>
uint64 ParallelPartition(ListType *list, uint64 N, bool (is_high)(ListType *obj, Args... args), Args... args) {
// Given a list of N objects and a slab
// partition the list, placing objects with value==slab at the end of the list.
// Return the starting index of the transition point.
// An arbitrary number of arguments may be passed to the is_high() partition function.
    
    if (N == 0)
        return 0;
    int nthread = omp_get_max_threads();
    #ifdef PPTEST
    fmt::print("Running on {:d} OMP threads\n", nthread);
    #endif

    if (2*N<(uint64) nthread) nthread = 1;    // If the list is too small, just do scalar
    std::vector<SwapBlock> plan;
    uint64 width = N/nthread+1;
    nthread = (int) ceil(1.*N/width);  // only use threads that will have work to do
    PartitionBlock *sublist = new PartitionBlock[nthread];

    // Split the list and partition the sub-lists
    #pragma omp parallel for schedule(static)
    for (int t=0; t<nthread; t++) 
        sublist[t].Partition(t*width, std::min((t+1)*width, N), list, is_high, args...);

    #ifdef PPTEST
    for (int t=0; t<nthread; t++) 
        fmt::print("Thread {:d}: begin={:d}, mid={:d}, end={:d}, {:d} low, {:d} high\n", t, 
                sublist[t].begin, sublist[t].mid, sublist[t].end, 
                sublist[t].numlow(), sublist[t].numhigh());
    #endif

    // Count the number of lows
    uint64 firsthigh = 0;
    for (int t=0; t<nthread; t++) firsthigh += sublist[t].numlow();

    // Now we need a work plan to do the large-scale swaps.
    // This is a lot of lines of code, but there's almost no work.
    int lowlist = nthread-1;    // The last sub-list
    uint64 lowactive = sublist[lowlist].numlow();  // The number that remain available
    uint64 lowend = sublist[lowlist].mid;  // The endpoint of the low list

    uint64 totalswap=0;     // The number of swaps we'll be doing

    for (int t=0; t<lowlist; t++) {
        // Need to find a home for the high particles in this sub-list.
        uint64 highactive = sublist[t].numhigh();    // The number that remain active
        uint64 highstart = sublist[t].mid;    // The first high object
        while (highactive>0 && t<lowlist) {
            // We have particles left to swap, and a potential place to do it
            uint64 nswap = std::min(lowactive, highactive);
            // We've found nswap particles that are mutually available
            // Record the trade plan
            if (nswap>0) {
                lowend -= nswap;   // We do this first, because we're marking the bottom index of both
                plan.push_back(SwapBlock(highstart, lowend, nswap));
                lowactive -= nswap;
                highactive -=nswap;
                highstart += nswap;
                totalswap += nswap;
            }
            if (lowactive==0) {
                // Go to the next sublist
                lowlist--;
                lowactive = sublist[lowlist].numlow();  
                lowend = sublist[lowlist].mid;  
            }
        }
        // When this is done, either we paired up all of the high points in sublist[t]
        // or we've run out of low points in all higher lists.
    }  
    // The first high point in the final list will be at highstart;

    #ifdef PPTEST
    for (uint64 j=0; j<plan.size(); j++) 
        fmt::print("Plan {:d}: high={:d}, low={:d}, nswap={:d}\n", 
                j, plan[j].highstart, plan[j].lowstart, plan[j].nswap);
    #endif                

    // We now want to spread the work more evenly over the threads.
    // We reuse the sublist class for this.
    uint64 perthread = totalswap/nthread+1;    // Seeking this much per thread
    #ifdef PPTEST
    fmt::print("Total second-stage swaps = {:d}, seeking {:d} per thread\n", totalswap, perthread);
    #endif
    int thread = 0;
    uint64 threadtot = 0;
    for (uint64 j=0; j<plan.size(); j++) {
        while (plan[j].nswap>0) {
            uint64 use = std::min(plan[j].nswap, perthread-threadtot);
            sublist[thread].threadplan.push_back(SwapBlock(plan[j].highstart, plan[j].lowstart, use));
            plan[j].highstart += use;
            plan[j].lowstart += use;
            plan[j].nswap -= use;
            threadtot += use;
            if (threadtot>=perthread && thread<nthread-1) {
                // We've exhausted this thread; move to the next one
                thread++;
                threadtot = 0;
            }
        }
    }

    #ifdef PPTEST
    for (int t=0; t<nthread; t++)
        for (uint64 j=0; j<sublist[t].threadplan.size(); j++) 
            fmt::print("Thread {:d}, Plan {:d}: high={:d}, low={:d}, nswap={:d}\n", 
                    t, j, sublist[t].threadplan[j].highstart, 
                    sublist[t].threadplan[j].lowstart, sublist[t].threadplan[j].nswap);
    #endif
    
    // Done with the serial planning.
    // Now let each thread do its swaps.
    #pragma omp parallel for schedule(static)
    for (int j=0; j<nthread; j++) sublist[j].swap(list);

    // All done!
    delete[] sublist;
    return firsthigh;
}


#ifdef PPTEST

#include <sys/time.h>
#include "STimer.cc"

inline bool int_equal(int *obj, int val) {
    if (*obj==val) return true;
    else return false;
}

int main() {
    uint64 N = 1001;
    uint64 max = 135;
    int *list;
    list = (int *)malloc(sizeof(int)*N);
    for (uint64 j=0; j<N; j++) list[j] = j%max;

    // for (uint64 j=0; j<N; j++) fmt::print("{:d} ", list[j]);
    // fmt::print("\n");

    STimer Part;
    Part.Start();

    int slab = 0;
    fmt::print("\nSlab Goal = {:d}\n", slab);
    uint64 mid = ParallelPartition(list, N, slab, int_equal);
    fmt::print("MidPoint = {:d}\n", mid);
    for (uint64 j=0;j<mid; j++) if(list[j]==slab) fmt::print("Fail: list[{:d}] = {:d}\n", j, list[j]);
    for (uint64 j=mid;j<N; j++) if(list[j]!=slab) fmt::print("Fail: list[{:d}] = {:d}\n", j, list[j]);
    for (uint64 j=0;j<N; j++) if(list[j]>max||list[j]<0) fmt::print("Fail: list[{:d}] = {:d}\n", j, list[j]);

    slab = 2;
    fmt::print("\nSlab Goal = {:d}\n", slab);
    mid = ParallelPartition(list, N, slab, int_equal);
    fmt::print("MidPoint = {:d}\n", mid);
    for (uint64 j=0;j<mid; j++) if(list[j]==slab) fmt::print("Fail: list[{:d}] = {:d}\n", j, list[j]);
    for (uint64 j=mid;j<N; j++) if(list[j]!=slab) fmt::print("Fail: list[{:d}] = {:d}\n", j, list[j]);
    for (uint64 j=0;j<N; j++) if(list[j]>max||list[j]<0) fmt::print("Fail: list[{:d}] = {:d}\n", j, list[j]);

    slab = 2;
    fmt::print("\nSlab Goal = {:d}\n", slab);
    mid = ParallelPartition(list, N, slab, int_equal);
    fmt::print("MidPoint = {:d}\n", mid);
    for (uint64 j=0;j<mid; j++) if(list[j]==slab) fmt::print("Fail: list[{:d}] = {:d}\n", j, list[j]);
    for (uint64 j=mid;j<N; j++) if(list[j]!=slab) fmt::print("Fail: list[{:d}] = {:d}\n", j, list[j]);
    for (uint64 j=0;j<N; j++) if(list[j]>max||list[j]<0) fmt::print("Fail: list[{:d}] = {:d}\n", j, list[j]);

    slab = 5;
    fmt::print("\nSlab Goal = {:d}\n", slab);
    mid = ParallelPartition(list, N, slab, int_equal);
    fmt::print("MidPoint = {:d}\n", mid);
    for (uint64 j=0;j<mid; j++) if(list[j]==slab) fmt::print("Fail: list[{:d}] = {:d}\n", j, list[j]);
    for (uint64 j=mid;j<N; j++) if(list[j]!=slab) fmt::print("Fail: list[{:d}] = {:d}\n", j, list[j]);
    for (uint64 j=0;j<N; j++) if(list[j]>max||list[j]<0) fmt::print("Fail: list[{:d}] = {:d}\n", j, list[j]);

    Part.Stop();
    fmt::print("Elapsed Time: {:f}\n", Part.Elapsed());

    // for (uint64 j=0; j<N; j++) fmt::print("{:d} ", list[j]);
    // fmt::print("\n");
    free(list);
    return 0;
}

#endif
