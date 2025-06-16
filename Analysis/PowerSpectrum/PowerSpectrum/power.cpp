// Copyright 2012-2025 The Abacus Developers
// SPDX-License-Identifier: GPL-3.0-or-later

/*
* power.cc
*    faster, more accurate power spectrum from input particles
*
*    Currently only incore. Out of core version is in dev
*
*  Created on: Sep 27, 2012
*      Author: dferrer
*  Updated: 2016, lgarrison
*/

#include "threevector.hh"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_sf_bessel.h>
#include <time.h>
#include <complex.h>
#include <fftw3.h>
#include <math.h>
#include <omp.h>
#include <fenv.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <fcntl.h>
#include <string.h>
#define __STDC_FORMAT_MACROS
#include <inttypes.h>

#include <stdint.h>
#include "iolib.cpp"

#ifdef DOUBLEPRECISION
#define FLOAT3 double3
#define FLOAT double
#define PI 3.141592653589793116
#define CABS std::abs
#else
#define FLOAT3 float3
#define FLOAT float
#define PI 3.141592653589793116f
#define CABS std::abs
#endif

#include "readabacus.h"

#define squ(x) ((x)*(x))
#define arr(a,x,y,z) a[gridN1D*gridN1D *(x) + gridN1D*(y) + (z)]
#define arr_2D(a,x,y) a[gridN1D*(x) + (y)]
#define arrR(a,x,y,z) a[gridN1D*2*(gridN1D/2+1) *(x) + 2*(gridN1D/2+1)*(y) + (z)]
#define arrR_2D(a,x,y) a[2*(gridN1D/2+1)*(x) + (y)]

#define den(x,y,z) arr(density,x,y,z)
#define den_2D(x,y) arr_2D(density,x,y)
#define denR(x,y,z) arrR(density,x,y,z)
#define denR_2D(x,y) arrR_2D(density,x,y)

extern "C" {
void set_num_threads(int nthreads){
    if(nthreads <= 0){
        nthreads = omp_get_num_procs();
    }
    omp_set_num_threads(nthreads);
}
    
    
#include "tsc.cpp"

// Rotate z-hat to align with rotate_to
// The domain should be zero-centered
// The resulting particles are trimmed if they go outside the "safe" zone
// The "safe" zone is a cube of edge length boxsize/sqrt(3)
// We probably could save more particles that we currently do in some cases
uint64_t rotate_positions(FLOAT3* positions, FLOAT* weights, uint64_t NP, double boxsize, double *rotate_to){
    // Normalize the new zhat
    double3 new_zhat(rotate_to[0], rotate_to[1], rotate_to[2]);
    new_zhat /= new_zhat.norm();
    double3 zhat(0.,0.,1.);
    double3 rotation_axis = zhat.cross(new_zhat);
    
    // We just normalized the vectors
    // so double check that round-off error isn't putting us outside [-1,1]
    double dot = zhat.dot(new_zhat);
    if(dot > 1.)
        dot = 1.;
    else if(dot < -1.)
        dot = -1;
    double rotation_angle = acos(dot);
    
    // Given an angle and axis, use Rodrigues' rotation formula to rotate
    double cos_ra = cos(rotation_angle);
    double sin_ra = sin(rotation_angle);
    double one_minus_cos_ra = 1. - cos_ra;
    uint64_t NP_keep = 0;
    for(uint64_t i = 0; i < NP; i++){
        FLOAT3 new_pos = positions[i]*cos_ra + rotation_axis.cross(positions[i])*sin_ra + rotation_axis*rotation_axis.dot(positions[i])*one_minus_cos_ra;
        // Discard particles outside the inner "safe" cube
        if(fabs(new_pos.x) > boxsize/(2*sqrt(3.)) || fabs(new_pos.y) > boxsize/(2*sqrt(3.)) || fabs(new_pos.z) > boxsize/(2*sqrt(3.)))
            continue;
        positions[NP_keep] = new_pos;
        if(weights != NULL)
            weights[NP_keep] = weights[NP];
        NP_keep++;
    }
    
    return NP_keep;
}

    
void do_tsc(FLOAT3 * positions, FLOAT * weights, FLOAT * density, uint64_t NP, uint64_t gridN1D, double boxsize, double *rotate_to, int projected, int rfft){
/*
 * Wrapper for all in-memory triangle-shaped cloud binning functions.
 *
 * Parameters
 * ----------
 * FLOAT3 *positions:
 *      The locations of the particles
 * FLOAT *weights:
 *      The binning weights.  No weights if NULL.
 * FLOAT *density:
 *      The output array for the binned density field
 * uint64_t NP:
 *      The length of the positions array
 * uint64_t gridN1D:
 *      The side length of the density grid
 * double boxsize:
 *      The physical domain side length.
 *      Positions are assumed to already be a zero-centered box of size `boxsize`.
 * double *rotate_to:
 *      A three-vector to which to align the z-axis of the simulation box
 * int projected:
 *      Whether to project the z-axis of the box and do 2D binning
 * int rfft:
 *      Whether the the output density array has shape suitable for in-place rfft
 */
    
    if(rotate_to != NULL){
        // This will cull particles
        NP = rotate_positions(positions, weights, NP, boxsize, rotate_to);
        // The binning functions will take this and appropriately scale the safe region
        // to the full grid domain
        boxsize /= sqrt(3.);
    }
    
    if(weights)
        if(projected)
            if (rfft)
                tsc_weighted_rfft_2D(positions, weights, density, NP, gridN1D, boxsize);
            else
                tsc_weighted_2D(positions, weights, density, NP, gridN1D, boxsize);
        else
            if (rfft)
                tsc_weighted_rfft(positions, weights, density, NP, gridN1D, boxsize);
            else
                tsc_weighted(positions, weights, density, NP, gridN1D, boxsize);
    else
        if(projected)
            if (rfft)
                tsc_rfft_2D(positions, density, NP, gridN1D, boxsize);
            else
                tsc_2D(positions, density, NP, gridN1D, boxsize);
        else
            if (rfft)
                tsc_rfft(positions, density, NP, gridN1D, boxsize);
            else
                tsc(positions, density, NP, gridN1D, boxsize);
}
    
#include "fft.cpp"
#include "ifft.cpp"
    
void shell_fft(FLOAT *kmag, FLOAT rmin, FLOAT rmax, FLOAT *out, uint64_t n){
    #pragma omp parallel for schedule(dynamic)
    for(uint64_t i = 0; i < n; i++){
        out[i] = rmax*rmax*gsl_sf_bessel_j1(rmax*kmag[i]) - rmin*rmin*gsl_sf_bessel_j1(rmin*kmag[i]);
    }
}

void annulus_fft(FLOAT *kmag, FLOAT smin, FLOAT smax, FLOAT *out, uint64_t n){
    #pragma omp parallel for schedule(dynamic)
    for(uint64_t i = 0; i < n; i++){
        out[i] = smax*gsl_sf_bessel_J1(smax*kmag[i]) - smin*gsl_sf_bessel_J1(smin*kmag[i]);
    }
}

}

