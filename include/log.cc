#ifndef INCLUDE_LOG
#define INCLUDE_LOG

#include <fstream>
#include <time.h>
#include <unistd.h>
#include "pprint.cc"

using namespace std;

template<typename IO, typename... Args>
void _log(IO &out, const char *name, const char *s, Args... args) {
    std::string ss(s);
    if(s[0]=='+') {
        ss.erase(0,1);
        out << "+" << std::setw(45) << std::left << " ";
    }
    else {
        time_t tnow = time(NULL);
        std::string time( ctime(&tnow) );
        out << time.substr(0,time.length()-1) << "  " << std::setw(20) << std::left << name+std::string("()  ");
    }
    fpprint(out, ss.c_str(), args...);
}

#define LOG(stream, ...) _log(stream, __func__, __VA_ARGS__)

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
