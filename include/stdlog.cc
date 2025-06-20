// Copyright 2012-2025 The Abacus Developers
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef STDLOG_INCLUDE
#define STDLOG_INCLUDE

// Use: stdlog.open() and stdlog.close() at the top level code.
// STDLOG(...) uses {fmt} so it should handle many C++ types.
// QUIT(...) writes to the log and aborts the program.
// assertf(test,...) does the test and writes to the log if failed.

// Please don't be bashful about writing to the Log!

// Please try to assure (via any of these mechanisms) that failures
// that could result from user errors or plausible dynamical events
// generate errors in the Log.  Try to limit bald assert() statements 
// to "these should never happen save by code bug" cases.

// STDLOG entries carry a verbosity flag.  Messages will be printed if the
// verbosity level is equal to or less than the global threshold value.
// verbosity=0 should be messages that we always want to print in production runs.
// verbosity=1 should be more detailed, e.g. listing the execution of major 
// 	portions of the code.
// verbosity=2 should be used for very detailed debugging.

#include "log.cc"
std::ofstream stdlog;
int stdlog_threshold_global = 1;

#define STDLOG(verbosity,...) do { if (verbosity<=stdlog_threshold_global) { \
	LOG(stdlog,__VA_ARGS__); } } while(0)

#define STDLOG_WITH_NAME(verbosity,name,...) do { if (verbosity<=stdlog_threshold_global) { \
	LOG_WITH_NAME(stdlog,name,__VA_ARGS__); } } while(0)

// C++-style variadic args don't work well across compilation units, so ensure a pure-string
// entry point for modules like the GPU directs
void stdlog_hook(int verbosity, const std::string &str){
    STDLOG(verbosity, str);
}

// Include a statement to put a global time stamp in the log.
#define STDLOG_TIMESTAMP do { \
	time_t tnow = time(NULL); \
	std::string time( ctime(&tnow) ); \
	STDLOG(0,"Timestamp {:s}\n", time.substr(0,time.length()-1)); \
    } while (0)

#define QUIT(...) do { fmt::print(stderr,"Fatal error (QUIT): "); \
    fmt::print(stderr, __VA_ARGS__); \
    STDLOG(0,"Fatal error (QUIT)\n"); STDLOG(0,__VA_ARGS__); \
	assert(0==98); } while(0)

#define WARNING(...) do { STDLOG(0,"WARNING:\n"); STDLOG(0,__VA_ARGS__); \
        fmt::print(stderr,"Warning: "); \
    fmt::print(stderr, __VA_ARGS__); } while(0)

// This is a form of assert, but it requires a printf style message.
// Be sure to terminate it with a \n.
// We encourage use of assertf and a clear explanation
// whenever user error could cause the assert fail.

#define assertf(_mytest,...) do { \
    if (UNLIKELY(!(_mytest))) { \
        STDLOG(0,"Failed Assertion: {:s}\n", #_mytest); STDLOG(1,__VA_ARGS__); \
        fmt::print(stderr,"Failed Assertion: {:s}\n", #_mytest); \
	fmt::print(stderr, __VA_ARGS__); \
        assert(0==99); \
    }} while(0)

#endif // STDLOG_INCLUDE
