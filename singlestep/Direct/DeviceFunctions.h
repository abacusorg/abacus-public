/*
 * Copyright 2012-2025 The Abacus Developers
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifdef CUDADIRECT

void GetDeviceInfo(size_t *memGB, std::string &name);

int GetNGPU();
extern void GPUSetup(int cpd, uint64 MaxBufferSize,
		int numberGPUs, int bufferperdevice,
		const std::vector<int> &ThreadCoreStart, int NThreadCores,
		const std::vector<int> &GPUQueueAssignments,
		int *maxsinkblocks, int *maxsourceblocks, int UsePinnedGPUMemory);
void GPUReset();
void StartGPUBufferThreads();

void print_gpu_mem_usage();

void get_cuda_timers(void*);

// These are the sizes we allocate on the GPU/pinned memory
extern int MaxSinkBlocks, MaxSourceBlocks;
extern int MaxNSink, MaxNSource;
extern size_t MaxSinkSize, MaxSourceSize;

#else
// No GPU memory limits; could be RAM limits if we wanted
int MaxSinkBlocks = INT32_MAX, MaxSourceBlocks = INT32_MAX;
int MaxNSink = INT32_MAX, MaxNSource = INT32_MAX;
size_t MaxSinkSize = INT32_MAX, MaxSourceSize = INT32_MAX;
#endif
