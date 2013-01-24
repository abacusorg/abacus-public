class PTimer {
public:
    PTimer(void); 
    ~PTimer();

    void Start(void);
    void Stop(void);

    double Elapsed(void);
    void Clear(void);

    int nprocs;

    int *timeron;

    struct timeval *tuse, *tstart, *timer;
};

PTimer::PTimer(void) { 

    nprocs = omp_get_num_procs();
    timeron = false;
    tuse = new timeval[nprocs];
    tstart = new timeval[nprocs];
    timer = new timeval[nprocs];
    timeron = new int[nprocs];

    for(int g=0; g<nprocs; g++) timerclear(&timer[g]);
    for(int g=0; g<nprocs; g++) timeron[g] = 0;
}

PTimer::~PTimer() {
    delete[] tuse;
    delete[] tstart;
    delete[] timer;
    delete[] timeron;
}

void PTimer::Start() {
    int g = omp_get_thread_num();

    assert(!timeron[g]);
    assert( gettimeofday( &(tstart[g]), (struct timezone *)NULL ) == 0 );
    timeron[g] = 1;
}

void PTimer::Stop(void) {
    int g = omp_get_thread_num();
    
    assert(timeron[g]);

    struct timeval dt;
    assert( gettimeofday( &(tuse[g]), (struct timezone *)NULL ) == 0 );
    timersub(&(tuse[g]), &(tstart[g]), &dt);
    timeradd(&dt, &(timer[g]), &(timer[g]));
    timeron[g] = 0;
}

void PTimer::Clear(void) {
    for(int g=0; g<nprocs; g++)  assert(!timeron[g]);  
    for(int g=0; g<nprocs; g++) timerclear(&(timer[g]));
}

double PTimer::Elapsed(void) {
    double sum = 0;
    for(int g=0; g<nprocs; g++) sum += timer[g].tv_sec + 1e-6*timer[g].tv_usec;
    return sum;
}
