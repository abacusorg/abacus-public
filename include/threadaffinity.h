#ifndef __THREADAFFINITY_H
#define __THREADAFFINITY_H

#include <sched.h>
#include <unistd.h>
#include <errno.h>

int set_core_affinity(int core_id) {
   int num_cores = sysconf(_SC_NPROCESSORS_ONLN);
   assertf(core_id < num_cores, "Tried to bind to core %d, but only %d cores detected\n", core_id, num_cores);

   cpu_set_t cpuset;
   CPU_ZERO(&cpuset);
   CPU_SET(core_id, &cpuset);

   pthread_t current_thread = pthread_self();    
   return pthread_setaffinity_np(current_thread, sizeof(cpu_set_t), &cpuset);
}

#endif