// Generate TSC functions with various options.  See tsc.h.

inline int wrap(int result, int max){
	while(result <0) result += max;
	while(result >= max) result-=max;
	assert(result >=0);
	return result;
}

#define TSC_FUNC_NAME tsc_weighted
#define WEIGHTED
#include "tsc.h"
#undef TSC_FUNC_NAME
#undef WEIGHTED

#define TSC_FUNC_NAME tsc_weighted_rfft
#define WEIGHTED
#define RFFT
#include "tsc.h"
#undef TSC_FUNC_NAME
#undef WEIGHTED
#undef RFFT

#define TSC_FUNC_NAME tsc_weighted_rfft_2D
#define WEIGHTED
#define RFFT
#define TWO_D
#include "tsc.h"
#undef TSC_FUNC_NAME
#undef WEIGHTED
#undef RFFT
#undef TWO_D

#define TSC_FUNC_NAME tsc_weighted_2D
#define WEIGHTED
#define TWO_D
#include "tsc.h"
#undef TSC_FUNC_NAME
#undef WEIGHTED
#undef TWO_D

#define TSC_FUNC_NAME tsc
#include "tsc.h"
#undef TSC_FUNC_NAME

#define TSC_FUNC_NAME tsc_rfft
#define RFFT
#include "tsc.h"
#undef TSC_FUNC_NAME
#undef RFFT

#define TSC_FUNC_NAME tsc_rfft_2D
#define RFFT
#define TWO_D
#include "tsc.h"
#undef TSC_FUNC_NAME
#undef RFFT
#undef TWO_D

#define TSC_FUNC_NAME tsc_2D
#define TWO_D
#include "tsc.h"
#undef TSC_FUNC_NAME
#undef TWO_D
