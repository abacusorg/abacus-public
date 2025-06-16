/*
 * Copyright 2012-2025 The Abacus Developers
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
#include <assert.h>
#include <string.h>
#include <sys/stat.h>
#include "io_packN.h"
#include "io_util.h"
#include "../universal_constants.h"
#include "../check_syscalls.h"
#include "../config_vars.h"
#include "../config.h"
#include "../particle.h"
#include <inttypes.h>

// Check that the packN file has a header
// and that the data is a multiple of N bytes
// Return the size of the data portion of the file.
void packN_validate_file(int packN, char *filename, uint64_t *header_size, uint64_t *data_size) {

    struct stat st;
    stat(filename, &st);
    off_t file_size = st.st_size;

    *header_size = packN_get_header_size(filename);
    *data_size = file_size - *header_size;
    
    if(*data_size % packN != 0){
        fprintf(stderr, "[Error] The size (%lu) of the pack%d data in %s is not a multiple of %d bytes\n",
            *data_size, packN, filename, packN);
        exit(1);
    }

    if (*header_size <= 0){
        fprintf(stderr, "[Error] No packN header in %s!\n", filename);
        exit(1);
    }
    if(*data_size <= 0){
        fprintf(stderr, "[Error] The size (%lu) of the packN data in %s is too small\n", *data_size, filename);
        exit(1);
    }
}

void packN_parse_header_info(packNheader *header){
    // Set the corresponding global variables
    h0 = header->H0 / 100;
    BOX_SIZE = header->BoxSizeHMpc;
    FORCE_RES = header->SofteningLength / (header->hMpc ? 1. : h0);  // Convert to h^-1 units if necessary
    TOTAL_PARTICLES = header->NP / (DOWNSAMPLE*DOWNSAMPLE*DOWNSAMPLE);
    if(TOTAL_PARTICLES * (DOWNSAMPLE*DOWNSAMPLE*DOWNSAMPLE) != header->NP){
        fprintf(stderr, "[Error] Downsample-per-dim = %" PRId64 " does not divide NP = %lld evenly.\n", DOWNSAMPLE, header->NP);
        exit(1);
    }
    if(IGNORE_PARTICLE_IDS && DOWNSAMPLE != 1){
        fprintf(stderr, "[Error] Downsampling needs PIDs!\n");
    }

    Om = header->Omega_M;
    Ol = header->Omega_DE;
    SCALE_NOW = header->ScaleFactor;
    FORCE_RES_PHYS_MAX = FORCE_RES*SCALE_NOW;
    PARTICLE_MASS = header->ParticleMassHMsun * (DOWNSAMPLE*DOWNSAMPLE*DOWNSAMPLE);
    W0 = header->w0;
    WA = header->wa;

    header->ppd_int = (uint64_t)round(header->ppd);
    
    AVG_PARTICLE_SPACING = cbrt(PARTICLE_MASS / (Om*CRITICAL_DENSITY));

    if(header->MaxPID > 0)
        MAX_PID = header->MaxPID;
    else if(header->packN == 14)
        MAX_PID = header->NP-1;
    else if(header->packN == 9)
        MAX_PID = (header->ppd_int-1) | ((header->ppd_int-1) << 16) | ((header->ppd_int-1) << 32);
    
    // Do a bit of validation
    if(header->packN == 14){
        if ((strncmp(header->OutputFormat, "Packed", 6) != 0)
            && (strncmp(header->OutputFormat, "Pack14", 6) != 0)){
            fprintf(stderr, "[Error] Header indicates that the OutputFormat is \"%s\" instead of \"Packed\" or \"Pack14\"", header->OutputFormat);
            exit(1);
        }
    } else if (header->packN == 9){
        if (strncmp(header->OutputFormat, "Pack9", 5) != 0){
            fprintf(stderr, "[Error] Header indicates that the OutputFormat is \"%s\" instead of \"Pack9\"", header->OutputFormat);
            exit(1);
        }
    }

    
    if (fabs(Om + Ol - 1.0) > 1e-5) {
        fprintf(stderr, "[Error] Halo Finder Not Currently Configured to Run on Cosmologies with Curvature. (Omega_Matter = %f, Omega_Lambda = %f!)\n", Om, Ol);
        exit(1);
    }
}
