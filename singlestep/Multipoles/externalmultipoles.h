/*
 * Copyright 2012-2025 The Abacus Developers
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef __D4DECL__
#define __D4DECL__
typedef struct { double v[4]; } d4;
#endif

#ifndef EXTERNALMULTIPOLES
#define EXTERNALMULTIPOLES

void DispatchCartesian2Reduced(int order, double *cartesian, double *reduced);

#ifdef AVXMULTIPOLES
extern "C" void DispatchMultipoleAVXKernel(int order, d4 *ip1x, d4 *ip2x, d4 *ip1y, d4 *ip2y, d4 *ip1z, d4 *ip2z, d4 *cx, d4 *cy, d4 *cz, d4 *globalM, d4 *mass1, d4 *mass2); 
#endif

#ifdef AVX512MULTIPOLES
void DispatchMultipole512Kernel(int order, FLOAT3 *xyz, int n, FLOAT3 center, double *CM);
#endif

#ifdef UNROLLEDMULTIPOLES
void DispatchMultipoleUnrolledKernel(int order, FLOAT3 *particles, int n, double3 center, double *CM);
#endif

#ifdef VSXMULTIPOLES
void DispatchMultipoleVSXKernel(int order, FLOAT3 *p, int n, double3 center, double *CM);
#endif

#endif
