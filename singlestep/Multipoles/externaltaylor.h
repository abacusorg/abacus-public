/*
 * Copyright 2012-2025 The Abacus Developers
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef EXTERNALTAYLOR
#define EXTERNALTAYLOR

#ifdef AVXMULTIPOLES
#ifndef __D4DECL__
#define __D4DECL__
typedef struct { double v[4]; } d4; 
#endif

extern "C" void DispatchTaylorAVXKernel(int order, d4 * px, d4 * py, d4 * pz, d4 * cx, d4 * cy, d4 * cz, d4 * tx, d4 * ty, d4 * tz, d4 * ax, d4 * ay, d4 * az ); 
#endif

#ifdef AVX512MULTIPOLES
#include "avx512_calls.h"

void DispatchTaylor512Kernel(int order, double *CT, FLOAT3 center, int n, FLOAT3 *xyz, FLOAT3 *acc);
#endif

#ifdef UNROLLEDMULTIPOLES
void DispatchTaylorUnrolledKernel(int order, FLOAT3 *p, int n, double3 center, double *CT, float3 *acc);
#endif

#ifdef VSXMULTIPOLES
void DispatchTaylorVSXKernel(int order, FLOAT3 *p, int n, double3 center, double *CT, float3 *acc);
#endif


#endif
