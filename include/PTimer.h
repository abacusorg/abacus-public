#ifndef PTIMER
#define PTIMER

#include <cassert>
#include <time.h>

class PTimer {
public:
    PTimer(void); 
    PTimer(int nthreads); 
    ~PTimer();

    // It really feels like there's a more elegant way to do this, but I haven't found it yet
    struct alignas(CACHE_LINE_SIZE) padded_timespec {
        struct timespec tstart;
        struct timespec tot;
        int on;
    };

    void Start(void);
    void Stop(void);
    void Start(int thread_num);
    void Stop(int thread_num);
    double Elapsed(void);
    void Clear(void);

    int nprocs;
    struct padded_timespec *timer;
};

#endif
