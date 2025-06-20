/*
 * Copyright 2012-2025 The Abacus Developers
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef INCLUDE_IO_INTERFACE
#define INCLUDE_IO_INTERFACE

#include "checksums.h"

enum io_blocking_mode { IO_NONBLOCKING,
                        IO_BLOCKING };

enum io_deletion_mode { IO_KEEP,
                        IO_DELETE };

// Every implementation of ReadFile/WriteFile should record
// their performance in these vars
// The map keys are directory names
//#include <map>
#include "tbb/concurrent_unordered_map.h"
typedef tbb::concurrent_unordered_map<std::string, STimer> TimerMap;
typedef tbb::concurrent_unordered_map<std::string, uint64> SizeMap;


// ChecksumMap = {dirname: {filename:CRC32, ...}, ...}
typedef tbb::concurrent_unordered_map<std::string, tbb::concurrent_unordered_map<std::string, CRC32>> ChecksumMap;

TimerMap BlockingIOReadTime, BlockingIOWriteTime;
TimerMap NonBlockingIOReadTime, NonBlockingIOWriteTime;
SizeMap BlockingIOReadBytes, BlockingIOWriteBytes;
SizeMap NonBlockingIOReadBytes, NonBlockingIOWriteBytes;

ChecksumMap FileChecksums;
std::unordered_map<int,double> ChecksumTime;
std::unordered_map<int,uint64> ChecksumBytes;

int allow_directio_global = 1;    // Set to 0 if we're using a ramdisk, or a network file system

// Two quick functions so that the I/O routines don't need to know
// about the SB object.

void IO_SetIOCompleted(int arenatype, int arenaslab);
void IO_DeleteArena(int arenatype, int arenaslab);

void IO_Initialize(const fs::path &logfn);
void IO_Terminate();

void ReadFile(char *ram, uint64 sizebytes, int arenatype, int arenaslab,
    const fs::path &fn, off_t fileoffset, int blocking);
    // This prototype reads sizebytes into the location *ram, from file *fn
    // starting from an offset fileoffset.
    // If blocking is set, then don't return until it's done!
    // Otherwise, return immediately if the I/O module allows it.
    // If arena>=0, call SetIOCompleted

void WriteFile(char *ram, uint64 sizebytes, int arenatype, int arenaslab, const fs::path &fn,
    off_t fileoffset, int deleteafter, int blocking, int do_checksum, int use_fp);
    // This prototype writes sizebytes from the location *ram, to file *fn
    // starting from an offset fileoffset.
    // If blocking is set, then don't return until it's done!
    // Otherwise, return immediately if the I/O module allows it.
    // If arena>=0, consider whether to delete the arena.

#endif // INCLUDE_IO_INTERFACE
