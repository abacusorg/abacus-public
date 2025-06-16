// Copyright 2012-2025 The Abacus Developers
// SPDX-License-Identifier: GPL-3.0-or-later

// Generate inverse FFT functions with various options.  See ifft.h.

#define IFFT_FUNC_NAME ifftInPlace
#include "ifft.h"
#undef IFFT_FUNC_NAME

#define IFFT_FUNC_NAME ifftInPlace2D
#define TWO_D
#include "ifft.h"
#undef IFFT_FUNC_NAME
#undef TWO_D
