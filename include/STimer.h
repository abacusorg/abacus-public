#ifndef STIMER
#define STIMER

#include <cassert>
#include <sys/time.h>

class STimer {
public:
    STimer();
    ~STimer();
    void Start(void);
    void Stop(void);
    double Elapsed(void);
    void Clear(void);
    void increment(struct timeval dt);
    struct timeval get_timer(void);
    int timeron;
    struct timeval tuse, tstart, timer;
};

struct timeval scale_timer(double s, struct timeval t);
#endif
