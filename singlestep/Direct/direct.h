/*
 * Copyright 2012-2025 The Abacus Developers
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef DIRECT
#define DIRECT

#ifdef AVXDIRECT
    #include <immintrin.h>
    #include "avxdirect.h"
    #include "avxdirectdouble.h"
    #include "avxdirectfloatNR.h"
#endif

class Direct {
public:

#ifdef AVXDIRECT
    AVXDirectDouble *directdouble;
    AVXDirectFloatNR *directfloatnr;
#endif

    Direct(void); 
    ~Direct(void);

    void Execute(FLOAT3 *sinks, FLOAT3 *sources, int nsinks, int nsources, FLOAT3 delta, FLOAT eps2, FLOAT3 *SA);
    void AVXExecute(FLOAT3 *sinks, FLOAT3 *sources, int nsinks, int nsources, FLOAT3 delta, FLOAT eps2, FLOAT3 *SA);
};

#endif
