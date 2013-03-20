#ifndef STDLOG_INCLUDE
#define STDLOG_INCLUDE

// Use: stdlog.open() and stdlog.close() at the top level code.
// STDLOG(...) uses pprint() so it should handle many C++ types.
// QUIT(...) writes to the log and aborts the program.
// assertf(test,...) does the test and writes to the log if failed.

// Please don't be bashful about writing to the Log!

// Please try to assure (via any of these mechanisms) that failures
// that could result from user errors or plausible dynamical events
// generate errors in the Log.  Try to limit bald assert() statements 
// to "these should never happen save by code bug" cases.

#include "log.cc"
std::ofstream stdlog;
#define STDLOG(...) { LOG(stdlog,__VA_ARGS__); stdlog.flush(); }

// Include a statement to put a global time stamp in the log.
#define STDLOG_TIMESTAMP do { \
	time_t tnow = time(NULL); \
	std::string time( ctime(&tnow) ); \
	STDLOG("Timestamp %s\n", time.substr(0,time.length()-1)); \
    } while (0)

#define QUIT(...) { STDLOG("Fatal error (QUIT)\n"); STDLOG(__VA_ARGS__); stdlog.flush(); \
        fprintf(stderr,"Fatal error (QUIT): "); \
	fpprint(std::cerr, __VA_ARGS__); \
	assert(0==98); }

// This is a form of assert, but it requires a printf style message.
// Be sure to terminate it with a \n.
// We encourage use of assertf and a clear explanation
// whenever user error could cause the assert fail.

#define assertf(_mytest,...) do { \
    if (!(_mytest)) { \
        STDLOG("Failed Assertion: %s\n", #_mytest); STDLOG(__VA_ARGS__); \
        fprintf(stderr,"Failed Assertion: %s\n", #_mytest); \
	fpprint(std::cerr, __VA_ARGS__); \
        assert(0==99); \
    }} while(0)

#endif // STDLOG_INCLUDE
