/* dependency.cpp
 * 
 * Class to handle slab-based preconditions and actions.  We load it
 * up with precondition and action functions, then it will attempt
 * to execute the slabs in order as their preconditions become satisfied.
 * Slabs will never be executed more than once or out of order.
 * 
 * Also provides a timing wrapper to the actions.
 * 
 * TODO: DWF believes this to be a good starting point for implementing
 * multi-node parallelism. Each node can a limited set of slabs stored
 * locally, which it will execute normally. Slabs near the boundary
 * between nodes will have modified preconditions that check that
 * computation of the slab on the other node is complete *and* that
 * the slab has been pushed to where we can locally access it. The actions
 * for boundary slabs will also need to be modified to ensure the preconditions
 * will be met. 
 */

#ifndef INCLUDE_DEPENDENCY
#define INCLUDE_DEPENDENCY

enum SpinFlag { NOT_ENOUGH_RAM, WAITING_FOR_IO, WAITING_FOR_GPU, NUM_SPIN_FLAGS };

class Dependency : public STimer {
public:
    int cpd;
    int *_executed_status;	// An array to keep track of what we've done
    int number_of_slabs_executed;   // Total number of slabs so far
    int last_slab_executed;         // The last slab we did
    int (*precondition)(int slab);
    void (*action)(int slab);
    
    static int *spin_flags;
    static STimer *spin_timers;
    static STimer global_spin_timer;

    ~Dependency(void) { 
        if(_executed_status != NULL) 
            delete[] _executed_status; 

//        printf("%s : %e \n", _name, Elapsed() );
     }

    void instantiate(   int _cpd, int _initialslab, 
                        int (*_precondition)(int), void (*_action)(int) ) { 

        action       = _action;
        precondition = _precondition;
        cpd          = _cpd;

        if(_executed_status!=NULL) { // changing cpd
            delete[] _executed_status;
        }
        
        _executed_status = new int[cpd];
        for(int s=0;s<cpd;s++) _executed_status[s] = 0;

        number_of_slabs_executed = 0; 
        last_slab_executed = _initialslab-1;
    }
                    
    int done(int s) { return _executed_status[ wrap(s)];  }
    int notdone(int s) { return _executed_status[wrap(s)]?0:1;  }
    int alldone(void) { return number_of_slabs_executed==cpd?1:0; }

    void do_action(int slab) {
        // We're taking an action!  Stop any spin timers
        for(int i = 0; i < NUM_SPIN_FLAGS; i++){
            spin_flags[i] = 0;
            if(spin_timers[i].timeron)
                spin_timers[i].Stop();
        }
        if(global_spin_timer.timeron)
            global_spin_timer.Stop();
        
        Start();
        (*action)(slab);
        Stop();
	    _executed_status[slab] = 1;
	    last_slab_executed = slab;
	    number_of_slabs_executed++;
    }

    int wrap(int s) { while(s<0) s += cpd; while(s>=cpd) s -= cpd; return s; }

    void Attempt(void) {
        // Take at most one action.
        int ws = wrap(last_slab_executed+1);
        if( notdone(ws) && precondition(ws) ) do_action(ws);
    }
    
    static void NotifySpinning(enum SpinFlag s){
        assertf(s < NUM_SPIN_FLAGS, "Tried to use spin flag %d (only NUM_SPIN_FLAGS=%d exist)!", s, NUM_SPIN_FLAGS);
        // If this is the first notification, we may not actually
        // be spinning yet.
        if(spin_flags[s]){
            // If it's the second time, start the timer.
            if(!spin_timers[s].timeron){
                spin_timers[s].Start();
                if(!global_spin_timer.timeron)
                    global_spin_timer.Start();
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
	// The check is wrapped, but the return value is not.
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

int *Dependency::spin_flags = new int[NUM_SPIN_FLAGS];
STimer *Dependency::spin_timers = new STimer[NUM_SPIN_FLAGS];
STimer Dependency::global_spin_timer;

#endif  // INCLUDE_DEPENDENCY
