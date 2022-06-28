/**
 * event_dependency.cpp
 * 
 * An EventDependency is a Dependency that is not slab-oriented. That is,
 * we do not expect it to run CPD times. In fact, we want to check it every
 * iteration of the event loop in timestep.cpp.
 *
 * EventDependencies also don't use the precondition-action pattern. Instead,
 * they execute an action and return whether they did real work, or just
 * checked for work (i.e. spinning). If real work was done, then the resulting
 * time is obviously a combination of time checking for work and then doing work,
 * but the latter is likely to greatly outweigh the former, so it's a good
 * approximation to call it all action time.
 *
 */

class EventDependency : public Dependency {

    virtual int action() = 0;

public:
    EventDependency(const char *name) : Dependency(name) { }

    int Attempt(){
        
        struct timespec stamp = SpinTimer.Stop();
        
        STimer timer;
        timer.StartFrom(stamp);
        int didsomething = action();
        stamp = timer.Stop();
        
        if(didsomething) action_timer.increment(timer.get_timer());
        else SpinTimer.increment(timer.get_timer());
        
        SpinTimer.StartFrom(stamp);
        
        return didsomething;
    }
};
