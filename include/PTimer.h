#ifndef PTIMER
#define PTIMER

#include <cassert>
#include <sys/time.h>

class PTimer {
public:
    PTimer(void); 
    ~PTimer();

    void Start(void);
    void Stop(void);
    double Elapsed(void);
    void Clear(void);
    struct timeval get_timer_seq(void);
    int nprocs;
    int *timeron;
    struct timeval *tuse, *tstart, *timer;
};

#endif
