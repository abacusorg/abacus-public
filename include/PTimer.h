#ifndef PTIMER
#define PTIMER

#include <cassert>
#include <sys/time.h>

class PTimer {
public:
    PTimer(void); 
    PTimer(int nthreads); 
    ~PTimer();

    void Start(void);
    void Stop(void);
    void Start(int thread_num);
    void Stop(int thread_num);
    double Elapsed(void);
    void Clear(void);
    struct timeval get_timer_seq(void);
    int nprocs;
    int *timeron;
    struct timeval *tuse, *tstart, *timer;
};

#endif
