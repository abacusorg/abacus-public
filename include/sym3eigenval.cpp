#ifdef SYMTEST 
#include <math.h>
#include <stdio.h>
#include <algorithm>
#include <stdlib.h>
#include <assert.h>
#endif 

void FindEigenvalues(double vxx, double vxy, double vxz, 
		double vyy, double vyz, double vzz, float sigmav[3]) {
    // Given a 3x3 symmetric matrix, find the eigenvalues.
    // https://en.wikipedia.org/wiki/Eigenvalue_algorithm#3.C3.973_matrices
    // TODO: Do we really need this to be in double precision?
    double p1 = vxy*vxy+vxz*vxz+vyz*vyz;
    // printf("%f %f %f\n", vxx, vxy, vxz);
    // printf("%f %f %f\n", vxy, vyy, vyz);
    // printf("%f %f %f\n\n", vxz, vyz, vzz);
    if (p1==0) {   // Diagonal matrix
	sigmav[0] = vxx; sigmav[1] = vyy; sigmav[2] = vzz;
	// Put the eigenvalues in descending order
	if (sigmav[1]<sigmav[2]) std::swap(sigmav[1], sigmav[2]);
	if (sigmav[0]<sigmav[1]) std::swap(sigmav[0], sigmav[1]);
	if (sigmav[1]<sigmav[2]) std::swap(sigmav[1], sigmav[2]);
    } else {
	double q = (vxx+vyy+vzz)/3.0;
	vxx -= q; vyy -= q; vzz -= q;
	double p = sqrt((vxx*vxx+vyy*vyy+vzz*vzz+2.0*p1)/6.0);
	// Our matrix is now A' = A-q*I.  
	// We're supposed to compute B = (A-q*I)/p and then r=det(B)/2
	// Instead we'll compute det(A') and then divide by 2*p**3.
	double r = vxx*vyy*vzz + 2.0*vxy*vxz*vyz - vxx*vyz*vyz - vyy*vxz*vxz - vzz*vxy*vxy;
	r /= 2.0*p*p*p;
	if (r > 1 && r < 1.00000001) r = 1.;
	else if (r < -1 && r > -1.00000001) r = -1.;
	assert(r>=-1 && r<=1);
	double phi = acos(r)/3.0;
	sigmav[0] = q + 2.0*p*cos(phi);
	sigmav[2] = q + 2.0*p*cos(phi+2.0/3.0*M_PI);
	sigmav[1] = 3.0*q - sigmav[0] - sigmav[2];
	// These eigenvalues are already guaranteed to be in descending order.
    }

    // FP imprecision may result in slightly negative eigenvalue
    if(sigmav[0] < 0) sigmav[0] = 0.;
    if(sigmav[1] < 0) sigmav[1] = 0.;
    if(sigmav[2] < 0) sigmav[2] = 0.;
    return;
}

#ifdef SYMTEST 
int main(int argc, char *argv[]) {
    // Enter vxx, vyy, vzz, vxy, vxz, vyz
    double in[6]; 
    float out[3];
    assert(argc==7);
    for (int j=0;j<6; j++) in[j] = atof(argv[j+1]);
    FindEigenvalues(
    	in[0], in[3], in[4], in[1], in[5], in[2], out);
    printf("%f %f %f\n", out[0], out[1], out[2]);
    return 0;
}

/* 

./a.out 1 1 1 0.1 0 0
./a.out 1 1 1 0 0.1 0
./a.out 1 1 1 0 0 0.1

./a.out 1 3 5 0 0 0.00001

./a.out 1 3 5 1 2 3
7.872983 1.000000 0.127017
Confirmed against Wolfram

*/

#endif 
