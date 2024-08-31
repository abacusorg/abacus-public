/* dependency.cpp
 * 
 * Class to handle slab-based preconditions and actions.  We load it
 * up with precondition and action functions, then it will attempt
 * to execute the slabs in order as their preconditions become satisfied.
 * Slabs will never be executed more than once or out of order.
 * 
 * Also provides a timing wrapper to the actions.
 */

#ifndef INCLUDE_DEPENDENCY
#define INCLUDE_DEPENDENCY

#include "stdlog.cc"

enum SpinFlag { NOT_ENOUGH_RAM, WAITING_FOR_IO, WAITING_FOR_GPU, WAITING_FOR_MPI, NUM_SPIN_FLAGS };

class Dependency {
public:
    const char *name;

    STimer action_timer;
    STimer precon_timer;

    static STimer SpinTimer;

    virtual int Attempt() = 0;

    Dependency(const char* _name) : name(_name) { }
    virtual ~Dependency() { }
    
    double Elapsed(){
        return action_timer.Elapsed();
    }

    double ElapsedPrecon(){
        return precon_timer.Elapsed();
    }

    static double ElapsedSpinning(){
        return SpinTimer.Elapsed();
    }
};


class SlabDependency : public Dependency {
public:
    int cpd;
    int *_executed_status;        // An array to keep track of what we've done
    int number_of_slabs_executed;   // Total number of slabs so far
            // This is the number we report; it can differ from the raw number
        // because we may intentionally queue some for re-execution.
    int last_slab_executed;         // The last slab we did
    int raw_number_executed;        // This is the number of times the action()
            // has been run, which may differ 

    int64_t num_particles;
    int64_t num_particles_with_ghost;
    
    static int *spin_flags;
    static STimer *spin_timers;
    static STimer global_spin_timer2;
    static STimer global_precon_timer;

    // implementers must define
    virtual int precondition(int slab) = 0;
    virtual void action(int slab) = 0;

    SlabDependency(const char* _name, int _cpd, int _initialslab)
            : Dependency(_name),
              cpd(_cpd) {
        
        _executed_status = new int[cpd];
        for(int s=0;s<cpd;s++) _executed_status[s] = 0;

        number_of_slabs_executed = 0; 
        raw_number_executed = 0; 
                
        num_particles = 0;
        num_particles_with_ghost = 0;

        last_slab_executed = wrap(_initialslab-1);
    }

    ~SlabDependency(void) { 
        if(_executed_status != NULL)
            delete[] _executed_status;
    }
    
    virtual int done(int s) {
        return _executed_status[ wrap(s)];
    }
    int notdone(int s) {
        return !done(s);
    }
    /// We can compare alldone to cpd, but usually we want total_slabs_on_node
    int alldone(void) {
        return number_of_slabs_executed==cpd?1:0;
    }
    int alldone(int total) {
        return number_of_slabs_executed==total?1:0;
    }

    void do_action(int slab) {
        // We're taking an action!  Stop any spin timers
        for(int i = 0; i < NUM_SPIN_FLAGS; i++){
            spin_flags[i] = 0;
            if(spin_timers[i].timeron)
                spin_timers[i].Stop();
        }
        if(global_spin_timer2.timeron)
            global_spin_timer2.Stop();
        
        STDLOG_WITH_NAME(1, name, "Entering {:s} action for slab {:d}\n", name, slab);

        struct timespec stamp = SpinTimer.Stop();
        action_timer.StartFrom(stamp);
        
        action(slab);
        
        stamp = action_timer.Stop();
        SpinTimer.StartFrom(stamp);

        STDLOG_WITH_NAME(1, name, "Exited {:s} action for slab {:d}\n", name, slab);
        
        _executed_status[slab] = 1;
        last_slab_executed = slab;
        number_of_slabs_executed++;
        raw_number_executed++;

        num_particles += SS->size(slab);
        num_particles_with_ghost += SS->size_with_ghost(slab);

    }

    int wrap(int s) {
        while(s<0) s += cpd;
        while(s>=cpd) s -= cpd;
        return s;
    }

    virtual int Attempt(void) {

        // Take at most one action.
        int ws = wrap(last_slab_executed+1);
        if( notdone(ws) ){
            STimer wc;
            wc.Start();
            int precon = precondition(ws);
            wc.Stop();
            global_precon_timer.increment(wc.get_timer());
            precon_timer.increment(wc.get_timer());
            if(precon) do_action(ws);
            return 1;
        }
        return 0;
    }
    
    static void NotifySpinning(enum SpinFlag s){
        // If this is the first notification, we may not actually
        // be spinning yet.
        if(spin_flags[s]){
            // If it's the second time, start the timer.
            if(!spin_timers[s].timeron){
                spin_timers[s].Start();
                if(!global_spin_timer2.timeron)
                    global_spin_timer2.Start();
            }
        } else {
            spin_flags[s] = 1;
        }
    }

    void mark_to_repeat(int slab) {
        slab = wrap(slab);
        _executed_status[slab] = 0;
        number_of_slabs_executed--;
    }

    int return_done_range(int end) {
        // Returns begin, such that [begin,end) is done and begin-1 is not.
        // We don't check if end is done.
        // The check is wrapped, but the return value is not, so begin<=end.

        int begin = wrap(end-1);
        while (done(begin) && end-begin<cpd) {
            begin--;
        }
        return begin;
    }

    void force_done(int s) {
        // This overrides the actions and declares it done.
        // last_slab_executed is never updated, nor is any timing done.
        // This is intended to be used when we have executed the action()
        // on another parallel node and are moving the results over.
    
        _executed_status[wrap(s)] = 1;
        number_of_slabs_executed++;
    }
};

STimer Dependency::SpinTimer;

int *SlabDependency::spin_flags = new int[NUM_SPIN_FLAGS];
STimer *SlabDependency::spin_timers = new STimer[NUM_SPIN_FLAGS];
STimer SlabDependency::global_spin_timer2;
STimer SlabDependency::global_precon_timer;

#include "event_dependency.cpp"

#endif  // INCLUDE_DEPENDENCY
