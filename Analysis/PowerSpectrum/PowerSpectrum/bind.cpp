// Copyright 2012-2025 The Abacus Developers
// SPDX-License-Identifier: GPL-3.0-or-later

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

// GSL
#include <gsl/gsl_sf_bessel.h>

// Parallel sorting for TSC
#include "tbb/parallel_sort.h"
#define SORT tbb::parallel_sort

//#include <algorithm>
//#define SORT std::sort

// TODO: rewrite this with a more specialized parallel numba partition algorithm

namespace py = pybind11;

// Sort of a list of particles, using the 'y' coordinate as the key
// This is used for dividing particles into y-chunks
class FLOAT3 {
    public:
    float r[3];
};
class FLOAT3w {
    public:
    float r[3];
    float w;
};

template <class T>
class ySortOperator {
    public:
    inline bool operator() (const T &pi, const T &pj ) const {
        return pi.r[1] < pj.r[1]; 
    }
};

void y_sorter(py::array_t<float, py::array::c_style | py::array::forcecast> particles){
    FLOAT3 *p = (FLOAT3 *) particles.mutable_unchecked<2>().mutable_data(0,0);
    uint64_t np = particles.shape(0);
    SORT(p, p + np, ySortOperator<FLOAT3>());
}

void y_sorter_weighted(py::array_t<float, py::array::c_style | py::array::forcecast> pw){
    FLOAT3w *p = (FLOAT3w *) pw.mutable_unchecked<2>().mutable_data(0,0);
    uint64_t np = pw.shape(0);
    SORT(p, p + np, ySortOperator<FLOAT3w>());
}

PYBIND11_MODULE(pslib, m) {
        m.def("y_sorter", &y_sorter);
        m.def("y_sorter_weighted", &y_sorter_weighted);
}
