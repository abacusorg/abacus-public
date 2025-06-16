// Copyright 2012-2025 The Abacus Developers
// SPDX-License-Identifier: GPL-3.0-or-later

// Generate FFT functions with various options.  See fft.h.

//aliased tsc window function from jeong (2010) chap. 7
#ifdef DOUBLEPRECISION
inline FLOAT window(FLOAT k){
    return 1 - squ(sin(k)) + 2./15. * squ(squ(sin(k)));
}
#else
inline FLOAT window(FLOAT k){
    return 1 - squ(sinf(k)) + 2.0f/15.0f * squ(squ(sinf(k)));
}
#endif

#define FFT_FUNC_NAME fftAndWindow
#define WINDOW
#include "fft.h"
#undef FFT_FUNC_NAME
#undef WINDOW

#define FFT_FUNC_NAME fftAndWindowInPlace
#define WINDOW
#define INPLACE
#include "fft.h"
#undef FFT_FUNC_NAME
#undef WINDOW
#undef INPLACE

#define FFT_FUNC_NAME fftAndWindowInPlace2D
#define WINDOW
#define INPLACE
#define TWO_D
#include "fft.h"
#undef FFT_FUNC_NAME
#undef WINDOW
#undef INPLACE
#undef TWO_D

#define FFT_FUNC_NAME fftAndWindow2D
#define WINDOW
#define TWO_D
#include "fft.h"
#undef FFT_FUNC_NAME
#undef WINDOW
#undef TWO_D

#define FFT_FUNC_NAME fft
#include "fft.h"
#undef FFT_FUNC_NAME

#define FFT_FUNC_NAME fftInPlace
#define INPLACE
#include "fft.h"
#undef FFT_FUNC_NAME
#undef INPLACE

#define FFT_FUNC_NAME fftInPlace2D
#define INPLACE
#define TWO_D
#include "fft.h"
#undef FFT_FUNC_NAME
#undef INPLACE
#undef TWO_D

#define FFT_FUNC_NAME fft2D
#define TWO_D
#include "fft.h"
#undef FFT_FUNC_NAME
#undef TWO_D
