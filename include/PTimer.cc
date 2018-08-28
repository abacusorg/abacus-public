#include <omp.h>
#include <cstdio>
#include <cstdlib>
#include "PTimer.h"

PTimer::PTimer(void) : PTimer(omp_get_max_threads()) { }

PTimer::PTimer(int nthreads) { 
    nprocs = nthreads;
    assert(posix_memalign((void **) &tuse, CACHE_LINE_SIZE, sizeof(padded_timeval)*nprocs) == 0);
    assert(posix_memalign((void **) &tstart, CACHE_LINE_SIZE, sizeof(padded_timeval)*nprocs) == 0);
    assert(posix_memalign((void **) &timer, CACHE_LINE_SIZE, sizeof(padded_timeval)*nprocs) == 0);
    assert(posix_memalign((void **) &timeron, CACHE_LINE_SIZE, sizeof(padded_int)*nprocs) == 0);

    for(int g=0; g<nprocs; g++) timerclear(&timer[g].t);
    for(int g=0; g<nprocs; g++) timeron[g].i = 0;
}

PTimer::~PTimer() {
    free(tuse);
    free(tstart);
    free(timer);
    free(timeron);
}

void PTimer::Start(void){
    Start(omp_get_thread_num());
}

void PTimer::Start(int thread_num) {
    int g = thread_num;
    
    assert(g < nprocs);  // If this fails omp_get_max_threads() may not be returning the global max # of threads
    /*if(timeron[g]){
        printf("Timer %d already on! nprocs = %d\n", g, nprocs);
    }*/
    
    assert(!timeron[g].i);
    assert( gettimeofday( &(tstart[g].t), (struct timezone *)NULL ) == 0 );
    timeron[g].i = 1;
}

void PTimer::Stop(void){
    Stop(omp_get_thread_num());
}

void PTimer::Stop(int thread_num) {
    int g = thread_num;
    
    assert(timeron[g].i);

    struct timeval dt;
    assert( gettimeofday( &(tuse[g].t), (struct timezone *)NULL ) == 0 );
    timersub(&(tuse[g].t), &(tstart[g].t), &dt);
    timeradd(&dt, &(timer[g].t), &(timer[g].t));
    timeron[g].i = 0;
}

void PTimer::Clear(void) {
    for(int g=0; g<nprocs; g++)  assert(!timeron[g].i);  
    for(int g=0; g<nprocs; g++) timerclear(&(timer[g].t));
}

double PTimer::Elapsed(void) {
    double sum = 0;
    for(int g=0; g<nprocs; g++) sum += timer[g].t.tv_sec + 1e-6*timer[g].t.tv_usec;
    return sum;
}

struct timeval PTimer::get_timer_seq(void) {
    struct timeval t;
    timerclear(&t);
    for(int g=0; g<nprocs; g++) timeradd(&(timer[g].t), &t, &t); 
    return t;
}
