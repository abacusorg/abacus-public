#ifndef INCLUDE_LOG
#define INCLUDE_LOG

#include <fstream>
#include <time.h>
#include <unistd.h>
#include "pprint.cc"
#include <time.h>
#include "STimer.h"
#include <mutex>

using namespace std;

// We will load this with the time at the beginning of the sim
struct timespec log_global_zero = { 0, 0 };

#define SET_LOG_REFERENCE_TIME do { \
	assert(clock_gettime(CLOCK_MONOTONIC, &log_global_zero) == 0); \
    } while (0)

std::mutex _log_mutex;

template<typename IO, typename... Args>
void _log(IO &out, const char *name, const char *s, Args... args) {
    std::string ss(s);
    _log_mutex.lock();
    if(s[0]=='+') {
        ss.erase(0,1);
        out << "+" << std::setw(38) << std::left << " ";
    }
    else {
	if (log_global_zero.tv_sec == 0) {
	    // First time, so initialize the global time reference
	    SET_LOG_REFERENCE_TIME;
	}
	struct timespec tnow, elapsed;
	assert(clock_gettime(CLOCK_MONOTONIC, &tnow) == 0);
	timespecsub(&tnow, &log_global_zero, &elapsed);
	double t = elapsed.tv_sec + 1e-9*elapsed.tv_nsec;
	out << std::right;
	fpprint(out, "%11.5f   ", t);
        out << std::setw(25) << std::left << name+std::string("()  ");

        // time_t tnow = time(NULL);
        // std::string time( ctime(&tnow) );
        // out << time.substr(0,time.length()-1) << "  " << std::setw(20) << std::left << name+std::string("()  ");
    }
    fpprint(out, ss.c_str(), args...);
    out.flush();
    _log_mutex.unlock();
}

#define LOG(stream, ...) _log(stream, __func__, __VA_ARGS__)
#define LOG_WITH_NAME(stream, name, ...) _log(stream, name, __VA_ARGS__)

#endif	// INCLUDE_LOG

/*
This routine provides a standardized logging function.
One should open a log with something like:

	std::ofstream out; out.open("mylog");

Then call with LOG(out,...) where the ... is printf-like string and list of
arguments.  The arguments can be C++ savvy; anything that C++ knows how to print
should echo in.

The time stamp and function name are automatically pre-pended.
If a string starts with +, then this is treated as a continuation of the
previous line.

 */
