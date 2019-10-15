#ifdef SYMTEST 
#include <math.h>
#include <stdio.h>
#include <algorithm>
#include <stdlib.h>
#include <assert.h>
#endif 



#ifdef OLD_CODE
// TODO!  THIS IS CURRENTLY BROKEN
/* I think that we have a conceptual problem: we need to subtract the 
first basis vector from the rest.  Need to do Gauss-Jordan w/pivoting?
Maybe we should be trying to construct orthogonal basis vectors first?
*/

// This will project away one vector, but only if it is non-zero to start with.
// And only if the dot product is substantial
#define NullProject(a,x,y,z) norm=x*x+y*y+z*z; dot=(a[0]*x+a[1]*y+a[2]*z); \
    if (norm>scale&&dot*dot>norm*tol) { dot/=norm; a[0]-=dot*x; a[1]-=dot*y; a[2]-=dot*z;} \
    printf("%f %f %f %f %f\n", a[0],a[1], a[2], norm, dot); \
    while(0)

/// We are given A-lambda*I and we need to find the null-space vector,
/// to be returned in major[3].
/// We will default to the axis given by fallback.
/// We use scale as the typical scale of the eigenvalues, so that
/// we can judge when we've gotten toward round-off problems.
void FindNullSpace(double fvxx, double fvyy, double fvzz, 
    double vxy, double vxz, double vyz, float scale, float *major, int fallback)
{
    // We're going to proceed by setting up a vector and then 
    // projecting out the 3 columns we've been given.
    double axis[3], dot, norm;
    double tol = 1e-7;
    scale = tol*tol*scale*scale;
    int j=0; 

    do {
        axis[0] = axis[1] = axis[2] = 0.0; 
        int k=fallback+j; if (k>2) k-=3; axis[k] = 1.0;
        NullProject(axis, fvxx, vxy, vxz);
        NullProject(axis,  vxy,fvyy, vyz);
        NullProject(axis,  vxz, vyz,fvzz);
        norm = axis[0]*axis[0]+axis[1]*axis[1]+axis[2]*axis[2];
        j++;
    } while (norm<scale && j<3);
    // It will be very rare to fail this check, but if it happens it means that
    // the initial direction has fallen purely in the non-null plane.
    // Repeat with another choice.

    if (norm<scale) {
        /* If we got here, it's because we were handed columns that spanned
        3-dim and hence didn't have a null space.  But that would happen
        if the eigenvalues were imperfect. */
        major[0] = major[1] = major[2] = 0; major[fallback] = 1.0;
    } else {
        norm = 1/sqrt(norm);
        major[0] = axis[0]*norm;
        major[1] = axis[1]*norm;
        major[2] = axis[2]*norm;
    }
    return;
}
        
    

#endif



void FindNullSpace(double fvxx, double fvyy, double fvzz, 
    double vxy, double vxz, double vyz, float scale, float *major, int fallback)
{
    // We're seeking the null space of the matrix.
    // This is code is very brute force:
    // It is extremely likely that the cross product between 
    // two of the columns is the null vector.  If not, then 
    // we cascade through other options.

    // We rescale the matrix so that we can assess when norms are close to zero
    // and to avoid over/under-flow.  [square norms are quartic in the code units, which can be larger/smaller than expected.]
    scale = 1.0/scale;
    fvxx *= scale;
    fvyy *= scale;
    fvzz *= scale;
     vxy *= scale;
     vxz *= scale;
     vyz *= scale;
    // printf("%f %f %f\n", fvxx, vxy, vxz);
    // printf("%f %f %f\n", vxy, fvyy, vyz);
    // printf("%f %f %f\n", vxz, vyz, fvzz);
    // printf("Cross 1 by 2\n");
    major[0] = vxy*vyz - vxz*fvyy;
    major[1] = vxz*vxy - fvxx*vyz;
    major[2] = fvxx*fvyy - vxy*vxy;
    float axis[3];
    float norm = major[0]*major[0]+major[1]*major[1]+major[2]*major[2];
    if (fabs(norm)<1e-15) {
        // We failed to find a subspace.  Try the first x third column
        // printf("Cross 1 by 3\n");
        major[0] = vxy*fvzz - vxz*vyz;
        major[1] = vxz*vxz - fvxx*fvzz;
        major[2] = fvxx*vyz - vxz*vxy;
        norm = major[0]*major[0]+major[1]*major[1]+major[2]*major[2];
    }
    if (fabs(norm)<1e-15) {
        // We failed to find a subspace.  Try the second x third column
        // printf("Cross 2 by 3\n");
        major[0] = fvyy*fvzz - vyz*vyz;
        major[1] = vyz*vxz - vxy*fvzz;
        major[2] = vxy*vyz - vxz*fvyy;
        norm = major[0]*major[0]+major[1]*major[1]+major[2]*major[2];
    }
    if (fabs(norm)<1e-15) {
        // Still failed; this indicates a degenerate eigenvalue.
        // All columns must be parallel.  However, some may be zero.
        // Pick the largest one
        axis[0] = fvxx*fvxx + vxy*vxy + vxz*vxz;
        axis[1] = fvyy*fvyy + vxy*vxy + vyz*vyz;
        axis[2] = fvzz*fvzz + vxz*vxz + vyz*vyz;
        if (axis[0]>=axis[1] && axis[0]>=axis[2]) {
            axis[0] = fvxx; axis[1] = vxy; axis[2] = vxz;
        } else {
            if (axis[1]>=axis[2]) {
                axis[0] = vxy; axis[1] = fvyy; axis[2] = vyz;
            } else {
                axis[0] = vxz; axis[1] = vyz; axis[2] = fvzz;
            }
        }
        // Now make up an axis orthogonal to the biggest column
        // Use the cross with (1,0,0)
        // printf("Cross with x: %e %e %e\n", axis[0], axis[1], axis[2]);
        major[0] = 0.0; major[1] = axis[2]; major[2] = -axis[1];
        norm = major[1]*major[1]+major[2]*major[2];
    }
    if (fabs(norm)<1e-15) {
        // Still failed; use a cross with (0,1,0)
        // printf("Cross with y\n");
        major[0] = -axis[2]; major[1] = 0.0; major[2] = axis[0];
        norm = major[0]*major[0]+major[2]*major[2];
    } 
    if (fabs(norm)<1e-15) {
        // Still failed; use a cross with (0,0,1)
        // printf("Cross with z\n");
        major[0] = axis[1]; major[1] = -axis[0]; major[2] = 0.0;
        norm = major[0]*major[0]+major[1]*major[1];
    } 
    if (fabs(norm)<1e-15) {
        // We have a triple degeneracy, punt.
        // printf("Punt\n");
        major[0] = 0.0; major[1] = 0.0; major[2] = 0.0; norm = 1.0;
        major[fallback] = 1.0;
    }
    // printf("Norm %e\n", norm);
    norm = sqrt(norm);
    major[0] /= norm; major[1] /= norm; major[2] /= norm;
    return;
}

void FindEigenvalues(double vxx, double vxy, double vxz, 
		double vyy, double vyz, double vzz, float sigmav[3],
        float *major, float *minor) {
    // Given a 3x3 symmetric matrix, find the eigenvalues.
    // https://en.wikipedia.org/wiki/Eigenvalue_algorithm#3.C3.973_matrices
    // Fill sigmav[3] with eigenvalues in descending order
    // If major[3] and minor[3] are given, fill with these eigenvectors
    // TODO: Do we really need this to be in double precision?
    double p;
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
        double fvxx = vxx-q; double fvyy = vyy-q; double fvzz = vzz-q;
        p = sqrt((fvxx*fvxx+fvyy*fvyy+fvzz*fvzz+2.0*p1)/6.0);
        // Our matrix is now A' = A-q*I.  
        // We're supposed to compute B = (A-q*I)/p and then r=det(B)/2
        // Instead we'll compute det(A') and then divide by 2*p**3.
        double r = fvxx*fvyy*fvzz + 2.0*vxy*vxz*vyz - fvxx*vyz*vyz - fvyy*vxz*vxz - fvzz*vxy*vxy;
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

    // If we're requested to compute eigenvectors, do that
    if (major!=NULL) {
        // We do this by finding the null space of (A-lam*I)
        float scale = fabs(sigmav[0])+fabs(sigmav[2]);
        FindNullSpace(vxx-sigmav[0], vyy-sigmav[0], vzz-sigmav[0],
            vxy, vxz, vyz, scale, major, 0);
    }
    if (minor!=NULL) {
        float scale = fabs(sigmav[0])+fabs(sigmav[2]);
        FindNullSpace(vxx-sigmav[2], vyy-sigmav[2], vzz-sigmav[2],
            vxy, vxz, vyz, scale, minor, 2);
    }
    return;
}

#ifdef SYMTEST 
int main(int argc, char *argv[]) {
    // Enter vxx, vyy, vzz, vxy, vxz, vyz
    double in[6]; 
    float out[3];
    float major[3], minor[3];
    assert(argc==7);
    for (int j=0;j<6; j++) in[j] = atof(argv[j+1]);
    FindEigenvalues(
    	in[0], in[3], in[4], in[1], in[5], in[2], out, major, minor);
    printf("%f %f %f   %f %f %f    %f %f %f\n", out[0], out[1], out[2], 
        major[0], major[1], major[2],
        minor[0], minor[1], minor[2]);
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
