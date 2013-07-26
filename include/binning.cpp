/*
 * binning.cpp
 * Seperates out the triangular shaped cloud binning code from the power spectrum code. Other binning methods can also go in here
 *  Created on: Jul 20, 2013
 *      Author: dferrer
 */
#ifndef BINNING_CPP
#define BINNING_CPP
#define arr(a,x,y,z) a[gridN1D*gridN1D *x + gridN1D*y + z]
#define den(x,y,z) arr(density,x,y,z)
#define squ(x) ((x)*(x))

#include <stdatomic.h>

inline int wrap(int input, int max){
	int result = input % max;
	while(result <0){
		result += max;
	}
	return result;
}

void tsc(FLOAT3 * positions,FLOAT3 cc, FLOAT * density, long long int NP, int gridN1D,FLOAT boxsize){
	long long int n;
	
#ifdef PARALLELBIN
	#pragma omp parallel for schedule(dynamic,128)
	#endif
	for(n = 0; n < NP; n++){
		FLOAT px = (positions[n].x+cc.x)/boxsize * gridN1D;
		FLOAT py = (positions[n].y+cc.y)/boxsize * gridN1D;
		FLOAT pz = (positions[n].z+cc.z)/boxsize * gridN1D;

		//round to nearest cell center (we offset the grid .5 so we can use floor instead of round)
#ifdef DOUBLEPRECISION
		int ix = floor(px+.5);
		int iy = floor(py+.5);
		int iz = floor(pz+.5);
#else
		int ix = floor(px+.5f);
		int iy = floor(py+.5f);
		int iz = floor(pz+.5f);
#endif
		//calculate distance to cell center
		FLOAT dx = ix - px;
		FLOAT dy = iy - py;
		FLOAT dz = iz - pz;

		//find the tsc weights for each dimension
#ifdef DOUBLEPRECISION
		FLOAT wx =  .75 -      squ(dx);
		FLOAT wxm1 = .5 * squ(.5 + dx);
		FLOAT wxp1 = .5 * squ(.5 - dx);
		FLOAT wy =  .75 -      squ(dy);
		FLOAT wym1 = .5 * squ(.5 + dy);
		FLOAT wyp1 = .5 * squ(.5 - dy);
		FLOAT wz =  .75 -      squ(dz);
		FLOAT wzm1 = .5 * squ(.5 + dz);
		FLOAT wzp1 = .5 * squ(.5 - dz);
#else
		FLOAT wx =  .75f -       squ(dx);
		FLOAT wxm1 = .5f * squ(.5f + dx);
		FLOAT wxp1 = .5f * squ(.5f - dx);
		FLOAT wy =  .75f -       squ(dy);
		FLOAT wym1 = .5f * squ(.5f + dy);
		FLOAT wyp1 = .5f * squ(.5f - dy);
		FLOAT wz =  .75f -       squ(dz);
		FLOAT wzm1 = .5f * squ(.5f + dz);
		FLOAT wzp1 = .5f * squ(.5f - dz);
#endif

		//find the wrapped x,y,z grid locations of the points we need to change
		int ixm1 =wrap(ix-1,gridN1D);
		int iym1 =wrap(iy-1,gridN1D);
		int izm1 =wrap(iz-1,gridN1D);
		int ixw = wrap(ix,gridN1D);
		int iyw = wrap(iy,gridN1D);
		int izw = wrap(iz,gridN1D);
		int ixp1 =wrap(ix+1,gridN1D);
		int iyp1 =wrap(iy+1,gridN1D);
		int izp1 =wrap(iz+1,gridN1D);

		//change the 27 cells that the cloud touches
		atomic_fetch_add( den(ixm1,iym1,izm1) , wxm1*wym1*wzm1);
		atomic_fetch_add( den(ixw, iym1,izm1) , wx  *wym1*wzm1);
		atomic_fetch_add( den(ixp1,iym1,izm1) , wxp1*wym1*wzm1);

#pragma omp atomic
		den(ixm1,iyw ,izm1) += wxm1*wy  *wzm1;
#pragma omp atomic
		den(ixw, iyw ,izm1) += wx  *wy  *wzm1;
#pragma omp atomic
		den(ixp1,iyw ,izm1) += wxp1*wy  *wzm1;

#pragma omp atomic
		den(ixm1,iyp1,izm1) += wxm1*wyp1*wzm1;
#pragma omp atomic
		den(ixw, iyp1,izm1) += wx  *wyp1*wzm1;
#pragma omp atomic
		den(ixp1,iyp1,izm1) += wxp1*wyp1*wzm1;

#pragma omp atomic
		den(ixm1,iym1,izw ) += wxm1*wym1*wz  ;
#pragma omp atomic
		den(ixw, iym1,izw ) += wx  *wym1*wz  ;
#pragma omp atomic
		den(ixp1,iym1,izw ) += wxp1*wym1*wz  ;

#pragma omp atomic
		den(ixm1,iyw ,izw ) += wxm1*wy  *wz  ;
#pragma omp atomic
		den(ixw, iyw ,izw ) += wx  *wy  *wz  ;
#pragma omp atomic
		den(ixp1,iyw ,izw ) += wxp1*wy  *wz  ;

#pragma omp atomic
		den(ixm1,iyp1,izw ) += wxm1*wyp1*wz  ;
#pragma omp atomic
		den(ixw, iyp1,izw ) += wx  *wyp1*wz  ;
#pragma omp atomic
		den(ixp1,iyp1,izw ) += wxp1*wyp1*wz  ;

#pragma omp atomic
		den(ixm1,iym1,izp1) += wxm1*wym1*wzp1;
#pragma omp atomic
		den(ixw, iym1,izp1) += wx  *wym1*wzp1;
#pragma omp atomic
		den(ixp1,iym1,izp1) += wxp1*wym1*wzp1;

#pragma omp atomic
		den(ixm1,iyw ,izp1) += wxm1*wy  *wzp1;
#pragma omp atomic
		den(ixw, iyw ,izp1) += wx  *wy  *wzp1;
#pragma omp atomic
		den(ixp1,iyw ,izp1) += wxp1*wy  *wzp1;

#pragma omp atomic
		den(ixm1,iyp1,izp1) += wxm1*wyp1*wzp1;
#pragma omp atomic
		den(ixw, iyp1,izp1) += wx  *wyp1*wzp1;
#pragma omp atomic
		den(ixp1,iyp1,izp1) += wxp1*wyp1*wzp1;
	}
}
#endif
