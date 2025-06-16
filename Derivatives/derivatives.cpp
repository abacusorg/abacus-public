// Copyright 2012-2025 The Abacus Developers
// SPDX-License-Identifier: GPL-3.0-or-later

class Derivatives : public basemultipoles {
public:
    Derivatives(int order) : basemultipoles(order) { }

    void ReducedDerivatives(double3 p, double *d);
    void Reduced2Complete(double *rd, double *cd);
    //double factorial(int n) { double f = 1; for(int i=2;i<=n;i++) f*= i; return f; } //only used by DerivativeCiprianiSilvi. Consider removing.
    //double DerivativeCiprianiSilvi(double3 r, int a, int b, int c); //this function appears out of us. no calls to it within abacus. consider removing. 
}; 

//Use the recursion relationships discussed in Section 2.2 to construct the derivatives tensor D^(ABC)({x,y,z}). Will be used in order to calculate the sum over the infinite lattice \mathcal{D}^(ABC) = \sum_{n} D^(ABC) (n).
void Derivatives::ReducedDerivatives(double3 p, double *d) {
    int i,j;
    double x,y,z,Rm1,Rm2,Rm3,Rm5,Rm7;
	
    x = p.x;
    y = p.y;
    z = p.z;

	//Eqn. 2.25
    Rm1 = sqrt(1.0/(x*x+y*y+z*z));
    Rm2 = Rm1 * Rm1;
    Rm3 = Rm1 * Rm2;
    Rm5 = Rm3 * Rm2;
    Rm7 = Rm5 * Rm2;

    d[rmap(0,0,0)] = Rm1;
    if(order>=1) {
        d[rmap(1,0,0)] = -x*Rm3;
        d[rmap(0,1,0)] = -y*Rm3;
        d[rmap(0,0,1)] = -z*Rm3;
    }
    if(order>=2) {
        d[rmap(1,1,0)] = 3*x*y*Rm5;
        d[rmap(1,0,1)] = 3*x*z*Rm5;
        d[rmap(0,1,1)] = 3*y*z*Rm5;
    }
    if(order>=3) {
        d[rmap(1,1,1)] = -15*x*y*z*Rm7;
    }
	//

	//Eqn. 2.26
    for(i=2;i<=order;i++) d[rmap(i,0,0)] = -Rm2*( (2*(i)-1)*x*d[rmap(i-1,0,0)] + (i-1)*(i-1)*d[rmap(i-2,0,0)] );
    for(i=2;i<=order-1;i++) d[rmap(i,1,0)] = -Rm2*( (2*(i)-1)*x*d[rmap(i-1,1,0)] + (i-1)*(i-1)*d[rmap(i-2,1,0)] + 2*y*d[rmap(i,0,0)] );
    for(i=2;i<=order-1;i++) d[rmap(i,0,1)] = -Rm2*( (2*(i)-1)*x*d[rmap(i-1,0,1)] + (i-1)*(i-1)*d[rmap(i-2,0,1)] + 2*z*d[rmap(i,0,0)] );
    for(i=2;i<=order-2;i++) d[rmap(i,1,1)] = -Rm2*( (2*(i)-1)*x*d[rmap(i-1,1,1)] + (i-1)*(i-1)*d[rmap(i-2,1,1)] + 2*y*d[rmap(i,0,1)] + 2*z*d[rmap(i,1,0)] );

    for(j=2;j<=order;j++) d[rmap(0,j,0)] = -Rm2*( (2*(j)-1)*y*d[rmap(0,j-1,0)] + (j-1)*(j-1)*d[rmap(0,j-2,0)] );
    for(j=2;j<=order-1;j++) d[rmap(0,j,1)] = -Rm2*( (2*(j)-1)*y*d[rmap(0,j-1,1)] + (j-1)*(j-1)*d[rmap(0,j-2,1)] + 2*z*d[rmap(0,j,0)] );
    for(j=2;j<=order-1;j++) d[rmap(1,j,0)] = -Rm2*( 2*x*d[rmap(0,j,0)] + (2*(j)-1)*y*d[rmap(1,j-1,0)] + (j-1)*(j-1)*d[rmap(1,j-2,0)] );
    for(j=2;j<=order-2;j++) d[rmap(1,j,1)] = -Rm2*( 2*x*d[rmap(0,j,1)] + (2*(j)-1)*y*d[rmap(1,j-1,1)] + (j-1)*(j-1)*d[rmap(1,j-2,1)] + 2*z*d[rmap(1,j,0)] );
	//
	
	//Eqn. 2.27
    for(i=2;i<=order;i++)
        for(j=2;j<=order-i;j++)
            d[rmap(i,j,0)] = -Rm2*( (2*(i)-1)*x*d[rmap(i-1,j,0)] + (i-1)*(i-1)*d[rmap(i-2,j,0)] + 2*(j)*y*d[rmap(i,j-1,0)] + (j)*(j-1)*d[rmap(i,j-2,0)] );

    for(i=2;i<=order-1;i++)
        for(j=2;j<=order-1-i;j++)
            d[rmap(i,j,1)] = -Rm2*( (2*(i)-1)*x*d[rmap(i-1,j,1)] + (i-1)*(i-1)*d[rmap(i-2,j,1)] + 2*(j)*y*d[rmap(i,j-1,1)] + (j)*(j-1)*d[rmap(i,j-2,1)] + 2*z*d[rmap(i,j,0)] );
	//
}

//Use the trace-free property of the derivatives tensor to reconstruct remaining terms. (Eqn. 2.24).
void Derivatives::Reduced2Complete(double *rd, double *cd) {
    for(int k=0;k<=1;k++) 
        for(int i=0;i<=order-k;i++) 
            for(int j=0;j<=order-k-i;j++) 
                cd[cmap(i,j,k)] = rd[rmap(i,j,k)];

    for(int k=2;k<=order;k++)
        for(int i=0;i<=order-k;i++)
            for(int j=0;j<=order-k-i;j++)
                cd[ cmap(i,j,k) ] = - cd[ cmap(i+2,j,k-2) ] - cd[ cmap(i,j+2,k-2) ];
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//This function appears to be outdated. Grep-ing the abacus directory doesn't return any calls to this. Consider deleting? 
//double Derivatives::DerivativeCiprianiSilvi(double3 r, int a, int b, int c) {
//    double s,nr;
//    s = 0;

//    nr = r.norm();

//    for(int l=0;l<=a/2;l++) {
//        for(int m=0;m<=b/2;m++) {
//            for(int n=0;n<=c/2;n++) {
//                s += pow(-1.0,l+m+n) *
//                     factorial( 2*(a+b+c) - 2*(l+m+n)  )/
//                     (factorial(l) * factorial(m) * factorial(n)) /
//                     ( factorial(a-2*l)  * factorial(b-2*m) * factorial(c-2*n) ) /
//                     factorial(a+b+c-l-m-n) *
//                     pow(r.x/nr, a-2*l) *
//                     pow(r.y/nr, b-2*m) *
//                     pow(r.z/nr, c-2*n);
//            }
//        }
//    }

//    return s/pow(nr,a+b+c+1)/pow(2.0,a+b+c)*( factorial(a) * factorial(b) *factorial(c) *  pow(-1.0,a+b+c) ); 
//}
