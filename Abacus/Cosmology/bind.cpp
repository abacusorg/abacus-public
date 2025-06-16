// Copyright 2012-2025 The Abacus Developers
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmissing-field-initializers"
#include <pybind11/pybind11.h>
#pragma GCC diagnostic pop

#include "Cosmology.cpp"

namespace py = pybind11;

PYBIND11_MODULE(AbacusCosmo, m) {
    m.doc() = "The Abacus cosmology module";

    py::class_<MyCosmology>(m, "MyCosmology")
        .def(py::init<>())
        .def("w", &MyCosmology::w)
        .def("OmegaX", &MyCosmology::OmegaX)
        .def("dOmegaXdlna", &MyCosmology::dOmegaXdlna)
        .def_readwrite("Omega_m", &MyCosmology::Omega_m)
        .def_readwrite("Omega_smooth", &MyCosmology::Omega_smooth)
        .def_readwrite("Omega_K", &MyCosmology::Omega_K)
        .def_readwrite("Omega_DE", &MyCosmology::Omega_DE)
        .def_readwrite("H0", &MyCosmology::H0)
        .def_readwrite("w0", &MyCosmology::w0)
        .def_readwrite("wa", &MyCosmology::wa);

    py::class_<Cosmology>(m, "Cosmology")
        .def(py::init<double, MyCosmology&>())
        .def_readwrite("early", &Cosmology::early)
        .def_readwrite("today", &Cosmology::today)
        .def_readwrite("current", &Cosmology::current)
        .def_readwrite("next", &Cosmology::next)
        .def_readwrite("search", &Cosmology::search);

    py::class_<Epoch>(m, "Epoch")
        .def_readwrite("a", &Epoch::a)
        .def_readwrite("z", &Epoch::z)
        .def_readwrite("t", &Epoch::t)
        .def_readwrite("etaK", &Epoch::etaK)
        .def_readwrite("etaD", &Epoch::etaD)
        .def_readwrite("growth", &Epoch::growth)
        .def_readwrite("f_growth", &Epoch::f_growth)

        .def_readwrite("w", &Epoch::w)
        .def_readwrite("H", &Epoch::H)
        .def_readwrite("OmegaHat_X", &Epoch::OmegaHat_X)
        .def_readwrite("OmegaHat_m", &Epoch::OmegaHat_m)
        .def_readwrite("OmegaHat_K", &Epoch::OmegaHat_K);
}
