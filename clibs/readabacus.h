/*
 * Copyright 2012-2025 The Abacus Developers
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef __READ_ABACUS_H
#define __READ_ABACUS_H

#include <stdint.h>
#include <stdlib.h>

// Two sets of functions: one that loads data into floats, the other into doubles.
#ifdef __cplusplus
extern "C" {
#endif
uint64_t read_pack14(char *fn, size_t offset, int ramdisk, int return_vel, int zspace, int return_pid, void *out);
uint64_t read_pack14f(char *fn, size_t offset, int ramdisk, int return_vel, int zspace, int return_pid, void *out);

#ifdef __cplusplus
}
#endif

#endif