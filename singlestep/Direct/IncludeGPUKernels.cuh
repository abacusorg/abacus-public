// By using the file twice, we create direct() and direct_density().

#define COMPUTE_FOF_DENSITY_SET
#include "DirectGPUKernels.cuh"

#undef COMPUTE_FOF_DENSITY_SET
#include "DirectGPUKernels.cuh"
