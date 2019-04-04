#ifndef PTIMER
#define PTIMER

#include <cassert>
#include <sys/time.h>

class PTimer {
public:
    PTimer(void); 
    PTimer(int nthreads); 
    ~PTimer();

    // It really feels like there's a more elegant way to do this, but I haven't found it yet
    struct alignas(CACHE_LINE_SIZE) padded_timeval {
        struct timeval t;
    };

    struct alignas(CACHE_LINE_SIZE) padded_int {
        int i;
    };

    void Start(void);
    void Stop(void);
    void Start(int thread_num);
    void Stop(int thread_num);
    double Elapsed(void);
    void Clear(void);
    struct timeval get_timer_seq(void);
    int nprocs;
    padded_int *timeron;
    struct padded_timeval *tuse, *tstart, *timer;
};

#endif
