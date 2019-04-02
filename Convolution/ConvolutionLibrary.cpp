#include "header.cpp"
#include "threevector.hh"

#include "STimer.h"

#include "file.h"

#include "factorial.h"

#include "basemultipoles.h"

#include "ConvolutionParametersStatistics.cpp"

#include "InCoreConvolution.cpp"

#ifdef PARALLEL
#include "ParallelConvolution.cpp"
#else
#include "OutofCoreConvolution.cpp"
#endif
