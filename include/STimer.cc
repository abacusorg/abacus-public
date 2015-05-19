#include <cstdio>
#include <cstdlib>
#include "STimer.h"

struct timeval scale_timer(double s, struct timeval t) {
    struct timeval tp;
    tp.tv_sec = t.tv_sec * s;
    tp.tv_usec = t.tv_usec * s;
    return tp;
}

STimer::STimer(void) { 

    timerclear(&timer);
    timeron = 0;
}

STimer::~STimer() {
}

struct timeval STimer::get_timer(void) {
    return timer;
}

void STimer::increment(struct timeval dt) {
    timeradd(&dt, &(timer), &(timer));
}

void STimer::Start() {
    assert(!timeron); 
    assert( gettimeofday( &(tstart), (struct timezone *)NULL ) == 0 );
    timeron = 1;
}

void STimer::Stop(void) {
    assert( timeron );
    struct timeval dt;
    assert( gettimeofday( &(tuse), (struct timezone *)NULL ) == 0 );
    timersub(&(tuse), &(tstart), &dt);
    timeradd(&dt, &(timer), &(timer));
    timeron = 0;
}

void STimer::Clear(void) {
    assert(!timeron); 
    timerclear(&(timer));
}

double STimer::Elapsed(void) {
    return  timer.tv_sec + 1e-6*timer.tv_usec;
}
