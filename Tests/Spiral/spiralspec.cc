// Copyright 2012-2025 The Abacus Developers
// SPDX-License-Identifier: GPL-3.0-or-later

// a second set of data
typedef struct {
    double3 kvec;
    double3 phase;
    double Across;
    fs::path datafile;
} spiralparam;

// and a second specialization
template <>
void ParseHeader::register_vars(spiralparam &P) {
    installvector("kvec", &(P.kvec.x), 3, 1, MUST_DEFINE);
    installvector("phase", &(P.phase.x), 3, 1, MUST_DEFINE);
    installscalar("Across", P.Across, MUST_DEFINE);
    installscalar("datafile", P.datafile, MUST_DEFINE);
}
