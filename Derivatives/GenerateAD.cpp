// Copyright 2012-2025 The Abacus Developers
// SPDX-License-Identifier: GPL-3.0-or-later

#include "../include/header.cpp"
#include "quad_double.cpp"
#include "../include/threevector.hh"
#include <omp.h>

typedef ThreeVector<qd_real> qd_real3;

#define ODIM (32+1)

int _rmap[33*33];
int TL;

#define rmap(i,j,k) _rmap[ (i)*33 + j ]

void makeEvenrmap(void) {
    int m = 0;
    for(int i=0;i<=32;i++)
        for(int j=0;j<=32;j++) {
            _rmap[ i*33 + j ] = m;
            m++;
        }

    TL = m;
}

// read in the analytic far derivatives
qd_real AD[ODIM][ODIM];
void GetAD(int maxorder) {
    assert(maxorder<ODIM);
    
    for(int a=0;a<ODIM;a+=2) for(int b=0;b<ODIM;b+=2) AD[a][b] = 0;

    std::ifstream fp("Order32"); //reads in infinite lattice derivs from file, Order32?
    assert(fp.is_open());
    std::string num;
    for(int i=0; i<79; i++) {
        int a, b;
        fp >> a >> b >> num;
        AD[a][b] = qd_real(num);
    }
    fp.close();

    for(int a=0;a<ODIM;a+=2) 
        for(int b=0;b<a;b+=2) 
            if( a+b>2 && a+b<ODIM) 
                AD[a][b] = AD[b][a];
}

//Use the recursion relationships discussed in Section 2.2 to construct the derivatives tensor D^(ABC)({x,y,z}). Will be used in order to calculate the sum over the infinite lattice \mathcal{D}^(ABC) = \sum_{n} D^(ABC) (n).
void QD_ab0(qd_real3 p, qd_real *d) {
    int i,j;
    qd_real x,y,z,Rm1,Rm2,Rm3,Rm5,Rm7;

    int order = 32;

    x = p.x;
    y = p.y;
    z = p.z;

    Rm1 = "1";
    Rm1 = Rm1/sqrt(x*x+y*y+z*z);
    Rm2 = Rm1 * Rm1;
    Rm3 = Rm1 * Rm2;
    Rm5 = Rm3 * Rm2;
	
    d[rmap(0,0,0)] = Rm1;
    if(order>=1) {
        d[rmap(1,0,0)] = -x*Rm3;
        d[rmap(0,1,0)] = -y*Rm3;
    }
    if(order>=2) {
        d[rmap(1,1,0)] = 3*x*y*Rm5;
    }

    for(i=2;i<=order;i++) d[rmap(i,0,0)] = -Rm2*( (2*(i)-1)*x*d[rmap(i-1,0,0)] + (i-1)*(i-1)*d[rmap(i-2,0,0)] );
    for(i=2;i<=order-1;i++) d[rmap(i,1,0)] = -Rm2*( (2*(i)-1)*x*d[rmap(i-1,1,0)] + (i-1)*(i-1)*d[rmap(i-2,1,0)] + 2*y*d[rmap(i,0,0)] );

    for(j=2;j<=order;j++) d[rmap(0,j,0)] = -Rm2*( (2*(j)-1)*y*d[rmap(0,j-1,0)] + (j-1)*(j-1)*d[rmap(0,j-2,0)] );
    for(j=2;j<=order-1;j++) d[rmap(1,j,0)] = -Rm2*( 2*x*d[rmap(0,j,0)] + (2*(j)-1)*y*d[rmap(1,j-1,0)] + (j-1)*(j-1)*d[rmap(1,j-2,0)] );

    for(i=2;i<=order;i++) 
        for(j=2;j<=order-i;j++)
            d[rmap(i,j,0)] = -Rm2*( (2*(i)-1)*x*d[rmap(i-1,j,0)] + (i-1)*(i-1)*d[rmap(i-2,j,0)] + 2*(j)*y*d[rmap(i,j-1,0)] + (j)*(j-1)*d[rmap(i,j-2,0)] );
}


void doradii(int radius, qd_real *dd) {
    qd_real rd[8192]; 

    for(int m=0;m<2048;m++) dd[m] = 0;

    for(int l1=-radius;l1<=radius;l1++) //why not first octant only? --> because these will be used to subtract off the nearest neighbors from infinite lattice sum. 
        for(int l2=-radius;l2<=radius;l2++)
            for(int l3=-radius;l3<=radius;l3++) {
                qd_real x2 = l1*l1 + l2*l2 + l3*l3;
                if(x2>"0" && x2<= qd_real(radius*radius) ) {
                    qd_real3 r;
                    r.x = l1; r.y = l2; r.z = l3;
                    QD_ab0(r,rd);
                    for(int m=0;m<2048;m++) dd[m] += rd[m];
                }
            }

}


// make the far derivatives -- the analytic derivatives minus the lattice sum of derivatives out to innerradius
qd_real AnalyticDerivatives[ODIM][ODIM][ODIM];
void MakeAnalyticDerivatives(int maxorder, int innerradius) {
    assert(maxorder < ODIM);

    GetAD(maxorder); //fills up AD[a][b] from binary file Order32

    for(int a=0;a<ODIM;a++) 
        for(int b=0;b<ODIM;b++) 
            for(int c=0;c<ODIM;c++) 
                AnalyticDerivatives[a][b][c] = 0;

    for(int a=0;a<=maxorder;a+=2)
        for(int b=0;b<=maxorder-a;b+=2) 
                if( a+b>2) {
                    AnalyticDerivatives[a][b][0] = AD[a][b]; 
                }

    qd_real rd2[8192];

    doradii(innerradius, rd2);

    for(int a=0;a<=maxorder;a+=2)
        for(int b=0;b<=maxorder-a;b+=2)
                if( a+b>2) {
					fmt::print("{: 2d} {: 2d} {: 2f} {: 2f} {: 2f} \n", a, b, double(AnalyticDerivatives[a][b][0]), double(rd2[ rmap(a,b,0)] ),   double(AnalyticDerivatives[a][b][0] - rd2[ rmap(a,b,0)] )   );
                    AnalyticDerivatives[a][b][0] -= rd2[ rmap(a,b,0) ];
                }

}


int main(void) {
   unsigned int old_cw;
    fpu_fix_start(&old_cw);

    cout.precision(qd_real::_ndigits);
    int nd = qd_real::_ndigits + 8;

    int order = 32;

    makeEvenrmap();

	for(int innerradius=1; innerradius<=9; innerradius+=1) {
		if(innerradius == 9){innerradius = 16;};
		
        fmt::print("making innerradius = {:d}\n", innerradius);

        MakeAnalyticDerivatives(order, innerradius);

        fs::path fname = fmt::format("AD{:02d}_{:03d}.dat",order,innerradius);
        fmt::print("fname = {}\n", fname);
        ofstream fout(fname,ios::out);
        fout.precision(qd_real::_ndigits);
        
        fout << order << "  " << innerradius << endl;
        for(int a=0; a<=order; a+=2) {
            for(int b=a; b<=order-a; b+=2) {
                fout << setw(3) << a << "  " << setw(3) << b << "  ";
                fout << setw(nd) << AnalyticDerivatives[a][b][0] << endl;
            }
        }
        fout.close();
    }

    fpu_fix_end(&old_cw);
}
