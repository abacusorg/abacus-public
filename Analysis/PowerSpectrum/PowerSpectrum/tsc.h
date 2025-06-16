/*
 * Copyright 2012-2025 The Abacus Developers
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/* This is a macro template file meant to be loaded via #include "tsc.h"
 *  Must define:
 *  TSC_FUNC_NAME:
 *       The function name
 *  FLOAT, FLOAT3:
 *       The data types
 *  den, denR, den_2D, denR_2D:
 *       Ways of accessing the density array
 *
 *  Can define:
 *  WEIGHTED:
 *       will add FLOAT *weights as the second argument and will weight points accordingly
 *  RFFT:
 *      Assumes the last array dimension has length 2(gridN1d/2 + 1)
 *  TWO_D:
 *       2D binning.  Default is 3D.
 */

void TSC_FUNC_NAME(FLOAT3 * positions,
#ifdef WEIGHTED
                   FLOAT * weights,
#endif
                   FLOAT * density, uint64_t NP, uint64_t gridN1D, double boxsize){
    uint64_t n;
    for(n = 0; n < NP; n++){
        double px = (positions[n].x/boxsize + .5) * gridN1D;
        double py = (positions[n].y/boxsize + .5) * gridN1D;
#ifndef TWO_D
        double pz = (positions[n].z/boxsize + .5) * gridN1D;
#endif

        //round to nearest cell center (we offset the grid .5 so we can use floor instead of round)
        int ix = floor(px+.5);
        int iy = floor(py+.5);
#ifndef TWO_D
        int iz = floor(pz+.5);
#endif

        //calculate distance to cell center
        double dx = ix - px;
        double dy = iy - py;
#ifndef TWO_D
        double dz = iz - pz;
#endif

        //find the tsc weights for each dimension
        double wx   = .75 -     squ(dx);
        double wxm1 = .5 * squ(.5 + dx);
        double wxp1 = .5 * squ(.5 - dx);
        double wy   = .75 -     squ(dy);
        double wym1 = .5 * squ(.5 + dy);
        double wyp1 = .5 * squ(.5 - dy);
#ifndef TWO_D
        double wz   = .75 -     squ(dz);
        double wzm1 = .5 * squ(.5 + dz);
        double wzp1 = .5 * squ(.5 - dz);
#endif

        //find the wrapped x,y,z grid locations of the points we need to change
        int ixm1 =wrap(ix-1,gridN1D);
        int ixw = wrap(ix  ,gridN1D);
        int ixp1 =wrap(ix+1,gridN1D);
        int iym1 =wrap(iy-1,gridN1D);
        int iyw = wrap(iy  ,gridN1D);
        int iyp1 =wrap(iy+1,gridN1D);
#ifndef TWO_D
        int izm1 =wrap(iz-1,gridN1D);
        int izw = wrap(iz  ,gridN1D);
        int izp1 =wrap(iz+1,gridN1D);
#endif
        
#ifdef WEIGHTED
        double C = weights[n];
        #define WEIGHT_FAC * C
#else
        #define WEIGHT_FAC
#endif
        
#ifdef TWO_D        
    #ifdef RFFT
        #define DEN denR_2D
    #else
        #define DEN den_2D
    #endif
        
        // change the 9 cells the cloud touches
        DEN(ixm1,iym1) += wxm1*wym1 WEIGHT_FAC;
        DEN(ixw, iym1) += wx  *wym1 WEIGHT_FAC;
        DEN(ixp1,iym1) += wxp1*wym1 WEIGHT_FAC;

        DEN(ixm1,iyw ) += wxm1*wy   WEIGHT_FAC;
        DEN(ixw, iyw ) += wx  *wy   WEIGHT_FAC;
        DEN(ixp1,iyw ) += wxp1*wy   WEIGHT_FAC;

        DEN(ixm1,iyp1) += wxm1*wyp1 WEIGHT_FAC;
        DEN(ixw, iyp1) += wx  *wyp1 WEIGHT_FAC;
        DEN(ixp1,iyp1) += wxp1*wyp1 WEIGHT_FAC;        
#else
    #ifdef RFFT
        #define DEN denR
    #else
        #define DEN den
    #endif

        //change the 27 cells that the cloud touches
        DEN(ixm1,iym1,izm1) += wxm1*wym1*wzm1 WEIGHT_FAC;
        DEN(ixw, iym1,izm1) += wx  *wym1*wzm1 WEIGHT_FAC;
        DEN(ixp1,iym1,izm1) += wxp1*wym1*wzm1 WEIGHT_FAC;

        DEN(ixm1,iyw ,izm1) += wxm1*wy  *wzm1 WEIGHT_FAC;
        DEN(ixw, iyw ,izm1) += wx  *wy  *wzm1 WEIGHT_FAC;
        DEN(ixp1,iyw ,izm1) += wxp1*wy  *wzm1 WEIGHT_FAC;

        DEN(ixm1,iyp1,izm1) += wxm1*wyp1*wzm1 WEIGHT_FAC;
        DEN(ixw, iyp1,izm1) += wx  *wyp1*wzm1 WEIGHT_FAC;
        DEN(ixp1,iyp1,izm1) += wxp1*wyp1*wzm1 WEIGHT_FAC;
        
        DEN(ixm1,iym1,izw ) += wxm1*wym1*wz   WEIGHT_FAC;
        DEN(ixw, iym1,izw ) += wx  *wym1*wz   WEIGHT_FAC;
        DEN(ixp1,iym1,izw ) += wxp1*wym1*wz   WEIGHT_FAC;

        DEN(ixm1,iyw ,izw ) += wxm1*wy  *wz   WEIGHT_FAC;
        DEN(ixw, iyw ,izw ) += wx  *wy  *wz   WEIGHT_FAC;
        DEN(ixp1,iyw ,izw ) += wxp1*wy  *wz   WEIGHT_FAC;

        DEN(ixm1,iyp1,izw ) += wxm1*wyp1*wz   WEIGHT_FAC;
        DEN(ixw, iyp1,izw ) += wx  *wyp1*wz   WEIGHT_FAC;
        DEN(ixp1,iyp1,izw ) += wxp1*wyp1*wz   WEIGHT_FAC;
        
        DEN(ixm1,iym1,izp1) += wxm1*wym1*wzp1 WEIGHT_FAC;
        DEN(ixw, iym1,izp1) += wx  *wym1*wzp1 WEIGHT_FAC;
        DEN(ixp1,iym1,izp1) += wxp1*wym1*wzp1 WEIGHT_FAC;

        DEN(ixm1,iyw ,izp1) += wxm1*wy  *wzp1 WEIGHT_FAC;
        DEN(ixw, iyw ,izp1) += wx  *wy  *wzp1 WEIGHT_FAC;
        DEN(ixp1,iyw ,izp1) += wxp1*wy  *wzp1 WEIGHT_FAC;

        DEN(ixm1,iyp1,izp1) += wxm1*wyp1*wzp1 WEIGHT_FAC;
        DEN(ixw, iyp1,izp1) += wx  *wyp1*wzp1 WEIGHT_FAC;
        DEN(ixp1,iyp1,izp1) += wxp1*wyp1*wzp1 WEIGHT_FAC;
#endif
    }
}

#undef WEIGHT_FAC
#undef WZ
#undef IZW
#undef DEN
