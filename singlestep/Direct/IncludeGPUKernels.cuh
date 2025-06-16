// Copyright 2012-2025 The Abacus Developers
// SPDX-License-Identifier: GPL-3.0-or-later

// By using the file twice, we create direct() and direct_density().

#define COMPUTE_FOF_DENSITY_SET
#include "DirectGPUKernels.cuh"

#undef COMPUTE_FOF_DENSITY_SET
#include "DirectGPUKernels.cuh"
