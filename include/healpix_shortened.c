/* -----------------------------------------------------------------------------
 *
 *  Copyright (C) 1997-2016 Krzysztof M. Gorski, Eric Hivon, Martin Reinecke,
 *                          Benjamin D. Wandelt, Anthony J. Banday,
 *                          Matthias Bartelmann,
 *                          Reza Ansari & Kenneth M. Ganga
 *
 *
 *  This file is part of HEALPix.
 *
 *  HEALPix is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  HEALPix is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with HEALPix; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 *  For more information about HEALPix see http://healpix.sourceforge.net
 *
 *---------------------------------------------------------------------------*/


/* Abacus only needs one healpix routine, and we don't want to deal with 
   cfitsio, so we're redacting the file down to bare minimum.  This was 
   from Healpix v3.60.0, downloaded January 10, 2020. */ 

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>

// #include "chealpix.h"

static const double twothird=2.0/3.0;
static const double pi=3.141592653589793238462643383279502884197;
static const double twopi=6.283185307179586476925286766559005768394;
static const double halfpi=1.570796326794896619231321691639751442099;
static const double inv_halfpi=0.6366197723675813430755350534900574;

/*! Returns the remainder of the division \a v1/v2.
    The result is non-negative.
    \a v1 can be positive or negative; \a v2 must be positive. */
static double fmodulo (double v1, double v2)
  {
  if (v1>=0)
    return (v1<v2) ? v1 : fmod(v1,v2);
  double tmp=fmod(v1,v2)+v2;
  return (tmp==v2) ? 0. : tmp;
  }

/* utab[m] = (short)(
      (m&0x1 )       | ((m&0x2 ) << 1) | ((m&0x4 ) << 2) | ((m&0x8 ) << 3)
    | ((m&0x10) << 4) | ((m&0x20) << 5) | ((m&0x40) << 6) | ((m&0x80) << 7)); */
static const short utab[]={
#define Z(a) 0x##a##0, 0x##a##1, 0x##a##4, 0x##a##5
#define Y(a) Z(a##0), Z(a##1), Z(a##4), Z(a##5)
#define X(a) Y(a##0), Y(a##1), Y(a##4), Y(a##5)
X(0),X(1),X(4),X(5)
#undef X
#undef Y
#undef Z
};

/* 64bit functions */

#ifndef __BMI2__
static int64_t spread_bits64 (int v)
  {
  return  (int64_t)(utab[ v     &0xff])
       | ((int64_t)(utab[(v>> 8)&0xff])<<16)
       | ((int64_t)(utab[(v>>16)&0xff])<<32)
       | ((int64_t)(utab[(v>>24)&0xff])<<48);
  }

static int64_t xyf2nest64 (int64_t nside, int ix, int iy, int face_num)
  {
  return (face_num*nside*nside) + spread_bits64(ix) + (spread_bits64(iy)<<1);
  }

#else

static int64_t xyf2nest64 (int64_t nside, int ix, int iy, int face_num)
  {
  return (face_num*nside*nside)
    + _pdep_u64(ix,0x5555555555555555ull) + _pdep_u64(iy,0xaaaaaaaaaaaaaaaaull);
  }

#endif

static int64_t ang2pix_nest_z_phi64 (int64_t nside_, double z, double s,
  double phi)
{
  double za = fabs(z);
  double tt = fmodulo(phi,twopi) * inv_halfpi; /* in [0,4) */
  int face_num, ix, iy;

  if (za<=twothird) /* Equatorial region */
    {
    double temp1 = nside_*(0.5+tt);
    double temp2 = nside_*(z*0.75);
    int64_t jp = (int64_t)(temp1-temp2); /* index of  ascending edge line */
    int64_t jm = (int64_t)(temp1+temp2); /* index of descending edge line */
    int64_t ifp = jp/nside_;  /* in {0,4} */
    int64_t ifm = jm/nside_;
    face_num = (ifp==ifm) ? (ifp|4) : ((ifp<ifm) ? ifp : (ifm+8));

    ix = jm & (nside_-1);
    iy = nside_ - (jp & (nside_-1)) - 1;
    }
  else /* polar region, za > 2/3 */
    {
    int ntt = (int)tt, jp, jm;
    double tp, tmp;
    if (ntt>=4) ntt=3;
    tp = tt-ntt;
    if (s>-2.)
      tmp = nside_*s/sqrt((1.+za)/3.);
    else
      tmp = nside_*sqrt(3*(1-za));

    jp = (int64_t)(tp*tmp); /* increasing edge line index */
    jm = (int64_t)((1.0-tp)*tmp); /* decreasing edge line index */
    if (jp>=nside_) jp = nside_-1; /* for points too close to the boundary */
    if (jm>=nside_) jm = nside_-1;
    if (z >= 0)
      {
      face_num = ntt;  /* in {0,3} */
      ix = nside_ - jm - 1;
      iy = nside_ - jp - 1;
      }
    else
      {
      face_num = ntt + 8; /* in {8,11} */
      ix =  jp;
      iy =  jm;
      }
    }

  return xyf2nest64(nside_,ix,iy,face_num);
}

void vec2pix_nest64(int64_t nside, const double *vec, int64_t *ipix)
{
  double vlen=sqrt(vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2]);
  double cth = vec[2]/vlen;
  double sth=(fabs(cth)>0.99) ? sqrt(vec[0]*vec[0]+vec[1]*vec[1])/vlen : -5;
  *ipix=ang2pix_nest_z_phi64 (nside,cth,sth,atan2(vec[1],vec[0]));
}

