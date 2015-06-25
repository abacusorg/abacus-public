/* dependency.cpp

Class to handle slab-based preconditions and actions.  We load it
up with precondition and action functions, then it will attempt
to execute the slabs in order as their preconditions become satisfied.
Slabs will never be executed more than once or out of order.

Also provides a timing wrapper to the actions.

*/

#ifndef INCLUDE_DEPENDENCY
#define INCLUDE_DEPENDENCY

class Dependency : public STimer {
public:
    int cpd;
    int *_executed_status;	// An array to keep track of what we've done
    int number_of_slabs_executed;   // Total number of slabs so far
    int last_slab_executed;         // The last slab we did
    int (*precondition)(int slab);
    void (*action)(int slab);

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
        Start();
        (*action)(slab);
        Stop();
	    _executed_status[slab] = 1;
	    last_slab_executed = slab;
	    number_of_slabs_executed++;
    }

    int wrap(int s) { while(s<0) s += cpd; while(s>=cpd) s -= cpd; return s; }

    void Attempt(void) {
        //Start();
	// Take at most one action.
        if(number_of_slabs_executed==0) {
            // never executed anything before 
            int s = last_slab_executed+1;
            for(int slab=s;slab<s+cpd;slab++) {
                int ws = wrap(slab); 
                if(precondition(ws) && notdone(ws)) {
                    do_action(ws);
                    break;
                }
            }
        }
        else { 
            int ws = wrap(last_slab_executed+1);
            if( notdone(ws) && precondition(ws) ) do_action(ws);
        }
        //Stop();
    }
};

#endif  // INCLUDE_DEPENDENCY
