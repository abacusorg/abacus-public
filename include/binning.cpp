/*
 * binning.cpp
 * Seperates out the triangular shaped cloud binning code from the power spectrum code. Other binning methods can also go in here
 *  Created on: Jul 20, 2013
 *      Author: dferrer
 */
#define arr(a,x,y,z) a[gridN1D*gridN1D *x + gridN1D*y + z]
#define den(x,y,z) arr(density,x,y,z)
#define squ(x) ((x)*(x))

void tsc(FLOAT3 * positions,double3 cc, FLOAT * density, long long int NP, int gridN1D,double boxsize){
	long long int n;
	for(n = 0; n < NP; n++){
		double px = (positions[n].x+cc)/boxsize * gridN1D;
		double py = (positions[n].y+cc)/boxsize * gridN1D;
		double pz = (positions[n].z+cc)/boxsize * gridN1D;

		//round to nearest cell center (we offset the grid .5 so we can use floor instead of round)
		int ix = floor(px+.5);
		int iy = floor(py+.5);
		int iz = floor(pz+.5);

		//calculate distance to cell center
		double dx = ix - px;
		double dy = iy - py;
		double dz = iz - pz;

		//find the tsc weights for each dimension
		double wx = .75 -       squ(dx);
		double wxm1 = .5 * squ(.5 + dx);
		double wxp1 = .5 * squ(.5 - dx);
		double wy = .75 -       squ(dy);
		double wym1 = .5 * squ(.5 + dy);
		double wyp1 = .5 * squ(.5 - dy);
		double wz = .75 -       squ(dz);
		double wzm1 = .5 * squ(.5 + dz);
		double wzp1 = .5 * squ(.5 - dz);

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
		den(ixm1,iym1,izm1) += wxm1*wym1*wzm1;
		den(ixw, iym1,izm1) += wx  *wym1*wzm1;
		den(ixp1,iym1,izm1) += wxp1*wym1*wzm1;

		den(ixm1,iyw ,izm1) += wxm1*wy  *wzm1;
		den(ixw, iyw ,izm1) += wx  *wy  *wzm1;
		den(ixp1,iyw ,izm1) += wxp1*wy  *wzm1;

		den(ixm1,iyp1,izm1) += wxm1*wyp1*wzm1;
		den(ixw, iyp1,izm1) += wx  *wyp1*wzm1;
		den(ixp1,iyp1,izm1) += wxp1*wyp1*wzm1;

		den(ixm1,iym1,izw ) += wxm1*wym1*wz  ;
		den(ixw, iym1,izw ) += wx  *wym1*wz  ;
		den(ixp1,iym1,izw ) += wxp1*wym1*wz  ;

		den(ixm1,iyw ,izw ) += wxm1*wy  *wz  ;
		den(ixw, iyw ,izw ) += wx  *wy  *wz  ;
		den(ixp1,iyw ,izw ) += wxp1*wy  *wz  ;

		den(ixm1,iyp1,izw ) += wxm1*wyp1*wz  ;
		den(ixw, iyp1,izw ) += wx  *wyp1*wz  ;
		den(ixp1,iyp1,izw ) += wxp1*wyp1*wz  ;

		den(ixm1,iym1,izp1) += wxm1*wym1*wzp1;
		den(ixw, iym1,izp1) += wx  *wym1*wzp1;
		den(ixp1,iym1,izp1) += wxp1*wym1*wzp1;

		den(ixm1,iyw ,izp1) += wxm1*wy  *wzp1;
		den(ixw, iyw ,izp1) += wx  *wy  *wzp1;
		den(ixp1,iyw ,izp1) += wxp1*wy  *wzp1;

		den(ixm1,iyp1,izp1) += wxm1*wyp1*wzp1;
		den(ixw, iyp1,izp1) += wx  *wyp1*wzp1;
		den(ixp1,iyp1,izp1) += wxp1*wyp1*wzp1;
	}
}
