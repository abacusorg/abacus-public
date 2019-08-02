#include <omp.h>
#include <cstdio>
#include <cstdlib>
#include "PTimer.h"

PTimer::PTimer(void) : PTimer(omp_get_max_threads()) { }

PTimer::PTimer(int nthreads) { 
    nprocs = nthreads;
    assert(posix_memalign((void **) &timer, CACHE_LINE_SIZE, sizeof(padded_timespec)*nprocs) == 0);

    for(int g=0; g<nprocs; g++) {
        timer[g].tot.tv_sec = 0; timer[g].tot.tv_nsec = 0;
    }
    for(int g=0; g<nprocs; g++) timer[g].on = 0;
}

PTimer::~PTimer() {
    free(timer);
}

inline void PTimer::Start(void){
    Start(omp_get_thread_num());
}

inline void PTimer::Start(int thread_num) {
    int g = thread_num;
    
    assert(g < nprocs);  // If this fails omp_get_max_threads() may not be returning the global max # of threads
    /*if(timeron[g]){
        printf("Timer %d already on! nprocs = %d\n", g, nprocs);
    }*/
    
    assert(!timer[g].on);
    assert( clock_gettime( CLOCK_REALTIME, &(timer[g].tstart)) == 0 );
    timer[g].on = 1;
}

inline void PTimer::Stop(void){
    Stop(omp_get_thread_num());
}

inline void PTimer::Stop(int thread_num) {
    int g = thread_num;
    assert(timer[g].on);
    struct timespec dt;
    struct timespec *t = &(timer[g].tot);
    assert( clock_gettime( CLOCK_REALTIME, &dt) == 0 );
    t->tv_nsec += dt.tv_nsec-timer[g].tstart.tv_nsec;
    t->tv_sec += dt.tv_sec-timer[g].tstart.tv_sec;
    if (t->tv_nsec<0) { t->tv_nsec += 1000000000; t->tv_sec--; }
    else if (t->tv_nsec>1000000000) { t->tv_nsec -= 1000000000; t->tv_sec++; }
    timer[g].on = 0;
}

void PTimer::Clear(void) {
    for(int g=0; g<nprocs; g++)  assert(!timer[g].on);  
    for(int g=0; g<nprocs; g++) {
        timer[g].tot.tv_sec = 0; timer[g].tot.tv_nsec = 0;
    }
}

double PTimer::Elapsed(void) {
    for(int g=0; g<nprocs; g++)  assert(!timer[g].on);  
    double sum = 0;
    for(int g=0; g<nprocs; g++) sum += timer[g].tot.tv_sec + 1e-9*timer[g].tot.tv_nsec;
    return sum;
}
