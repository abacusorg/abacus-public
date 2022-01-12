/*

On NUMA systems, having cores operate on memory within their NUMA node is important
for performance.  But OpenMP doesn't provide any mechanisms by which to assign threads
to loop iterations that take this into account.  It's easy to accomplish with static
scheduling, but then one sacrifies load balancing.

To work around this, we implement a custom OpenMP scheduler that schedules threads
dynamically within a NUMA node; we'll implement it as macro called NUMA_FOR that
replaces an ordinary "for" statement.

This file contains the macro defintion, a few global helper variables, and init
and destructor routines for the globals.

Usage
-----

Given the following code:

#pragma omp parallel for schedule(static)
for(int y = 0; y < cpd; y++){
    for(int z = 0; z < cpd; z++){
        ProcessCell(slab,y,z);
    }
}

Replace the pragma and outer for loop with NUMA_FOR, as follows:

NUMA_FOR(y,0,cpd)
    for(int z = 0; z < cpd; z++){
        ProcessCell(slab,y,z);
    }
}

N.B. There is no open brace after NUMA_FOR (it must open the loop itself
to set the iteration variable), but don't forget the usual close brace.

Remember that this is still ultimately dynamic scheduling.  So one needs to keep
in mind the scheduling overhead; don't use this if the loop iterations are ultra-fast.

This scheduler is unlikely to help with memory-bandwidth-bound loops, like the kick
or drift. It might help keep work within a NUMA node, but it probably wouldn't matter much.

Details
-------
We implement this using OpenMP tasks and atomic addition, so it should be fairly efficient.
We take care to cache-line pad the ints that count iterations by defining an "pint64" class
(aligned 64-bit int).  At the pencil level on summit, we haven't seen noticable scheduling
overhead over static.

This scheduler doesn't support any reduction clauses, partly because OpenMP tasks don't
yet natively support reductions.  So one has to implement any reductions manually, unfortunately.
The pint64 class may be useful for this.

The scheduling tries to assign iterations to nodes based on the fraction of threads in that
node.  The NUMA parts of the init function are pretty verbose, but boil down to counting
threads on nodes.

The iteration variable is of type int64_t; should be safe for most applications.

The original idea for implementing a custom scheduler came from: https://stackoverflow.com/a/30591616

*/

#ifndef __NUMA_FOR_H
#define __NUMA_FOR_H

#include "numa.h"

#ifdef ENABLE_NUMA_FOR

// This messy macro chain is used to generate unique variable names
// in case NUMA_FOR is used more than one time in a context
// From: https://stackoverflow.com/a/17624752
#define PP_CAT(a, b) PP_CAT_I(a, b)
#define PP_CAT_I(a, b) PP_CAT_II(~, a ## b)
#define PP_CAT_II(p, res) res
#define UNIQUE_NAME(base) PP_CAT(base, __LINE__)

// This is the primary macro definition
#define NUMA_FOR(I,START,END)\
pint64 UNIQUE_NAME(_next_iter)[_N_numa_nodes];\
int64_t UNIQUE_NAME(_last_iter)[_N_numa_nodes];\
\
UNIQUE_NAME(_next_iter)[0] = (int64_t) (START);\
UNIQUE_NAME(_last_iter)[_N_numa_nodes-1] = (int64_t) (END);\
for(int _nn = _N_numa_nodes-1; _nn > 0; _nn--){\
    UNIQUE_NAME(_next_iter)[_nn] = (int64_t) _numa_cumulative_nthread[_nn]*\
    (UNIQUE_NAME(_last_iter)[_N_numa_nodes-1] - UNIQUE_NAME(_next_iter)[0])/_numa_cumulative_nthread[_N_numa_nodes];\
    UNIQUE_NAME(_last_iter)[_nn-1] = UNIQUE_NAME(_next_iter)[_nn];\
}\
_Pragma("omp parallel")\
    for(int64_t I; ; ){\
        int _nn = _thread_numa_nodes[omp_get_thread_num()];\
        _Pragma("omp atomic capture")\
        I = UNIQUE_NAME(_next_iter)[_nn].i++;\
        if(I >= UNIQUE_NAME(_last_iter)[_nn])\
            break;\


int *_thread_numa_nodes;
int *_numa_cumulative_nthread;
int _N_numa_nodes;

#endif  // ENABLE_NUMA_FOR

void init_numa_for(int nthreads, int const *core_assignments){

#ifdef ENABLE_NUMA_FOR

    if(omp_get_proc_bind() == omp_proc_bind_false){
        STDLOG(1, "OMP_PROC_BIND = false; NUMA_FOR loops will use dynamic scheduling\n");

        // Set up a basic dynamic scheduler, ignoring NUMA
        _N_numa_nodes = 1;
        _thread_numa_nodes = new int[nthreads]();  // all 0
        _numa_cumulative_nthread = new int[_N_numa_nodes+1];
        _numa_cumulative_nthread[0] = 0;
        _numa_cumulative_nthread[1] = nthreads;

        return;
    }

    assertf(numa_available() != -1, "NUMA not available?\n");

    // First, determine which NUMA nodes have active CPUs and index them contiguously
    int max_node = numa_max_node();
    STDLOG(3, "Max NUMA node: %d\n", max_node);
    struct bitmask *bm = numa_allocate_cpumask();
    int numa_node_squashed_index[max_node+1];  // index of node, skipping empties
    _N_numa_nodes = 0;  // number of non-empty numa nodes
    for(int i = 0; i <= max_node; i++){
        if(!numa_bitmask_isbitset(numa_all_nodes_ptr,i)){
            numa_node_squashed_index[i] = -1;
            continue;
        }
        assertf(numa_node_to_cpus(i, bm) == 0, "NUMA detection failed (errno: %d)\n", errno);
        int cpus_on_node = numa_bitmask_weight(bm);
        if(cpus_on_node)
            numa_node_squashed_index[i] = _N_numa_nodes++;
        else
            numa_node_squashed_index[i] = -1;
        STDLOG(3,"%d cpus on NUMA node %d, squashed index %d\n", cpus_on_node, i, numa_node_squashed_index[i]);
    }
    numa_free_cpumask(bm);

    // Now init the thread->node mapping, using the squashed indices
    _thread_numa_nodes = new int[nthreads];  // numa node index (squashed) of OpenMP threads
    _numa_cumulative_nthread = new int[_N_numa_nodes+1];  // total of OpenMP threads before this node
    int *_nthread_per_numa_node = new int[_N_numa_nodes]();  // 0 init; number of OpenMP threads per NUMA node

    for(int i = 0; i < nthreads; i++){
        int nncpu = numa_node_of_cpu(core_assignments[i]);
        assertf(nncpu != -1, "NUMA failure!\n");
        _thread_numa_nodes[i] = numa_node_squashed_index[nncpu];
        assertf(_thread_numa_nodes[i] != -1, "NUMA failure!\n");
        _nthread_per_numa_node[_thread_numa_nodes[i]]++;
        STDLOG(3,"Thread %d on core %d on NUMA node %d\n",
            i, core_assignments[i], _thread_numa_nodes[i]);
    }

    for(int i = 0; i < _N_numa_nodes; i++)
        STDLOG(1,"NUMA_FOR configured with %d threads on NUMA node %d\n", _nthread_per_numa_node[i], i);

    _numa_cumulative_nthread[0] = 0;
    for(int i = 1; i < _N_numa_nodes+1; i++){
        _numa_cumulative_nthread[i] = _numa_cumulative_nthread[i-1] + _nthread_per_numa_node[i-1];
    }
    assert(_numa_cumulative_nthread[_N_numa_nodes] == nthreads);

    delete[] _nthread_per_numa_node;

#endif // ENABLE_NUMA_FOR
}


// Free any global variables
// Can't use NUMA_FOR after this
void finish_numa_for(){
#ifdef ENABLE_NUMA_FOR
    delete[] _thread_numa_nodes;
    delete[] _numa_cumulative_nthread;
#endif
}


#ifndef ENABLE_NUMA_FOR

// If disabled, fall back to dynamic scheduling

#define NUMA_FOR(I,START,END)\
_Pragma("omp parallel for schedule(dynamic)")\
for(int64_t I = START; I < END; I++){

#endif // ENABLE_NUMA_FOR

// Cache-line aligned variable template class
// Some typedefs for ints and doubles are below
template<typename T>
class alignas(CACHE_LINE_SIZE) padded{
public:
    T i;

    padded() = default;

    padded(const T &_i){
        i = _i;
    }

    padded& operator=(const T &rhs){
        i = rhs;
        return *this;
    }

    padded& operator++(int){
        i++;
        return *this;
    }

    template<typename U>
    padded& operator+=(const U& rhs){
        i += rhs;
        return *this;
    }

    template<typename U>
    friend U& operator+=(U& lhs, const padded& rhs) {
        lhs += rhs.i;
        return lhs;
    }

    template<typename U>
    friend bool operator>=(const U& lhs, const padded& rhs) {
        return lhs >= rhs.i;
    }

    template<typename U>
    bool operator!=(const U& rhs) const {
        return rhs != *this;
    }

    template<typename U>
    friend bool operator!=(const U& lhs, const padded& rhs) {
        return lhs != rhs.i;
    }

    friend padded operator+(const T& lhs, const padded& rhs) {
        return pint64(lhs + rhs.i);
    }

    friend padded operator+(const padded& lhs, const padded& rhs) {
        return padded(lhs.i + rhs.i);
    }

    operator int64_t() const {
        return i;
    }

    friend ostream& operator<<(ostream &lhs, const padded& rhs){
        lhs << rhs.i;
        return lhs;
    }
};

using pint64 = padded<int64_t>;
using pdouble = padded<double>;

#endif
