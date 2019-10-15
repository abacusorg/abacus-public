//#define SYMTEST
#ifdef SYMTEST 
#include <math.h>
#include <stdio.h>
#include <algorithm>
#include <stdlib.h>
#include <assert.h>
#include <float.h>
#endif 


#ifdef BROKEN_CODE
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
        
    

//SOURCE: http://www.mymathlib.com/c_source/matrices/eigen/jacobi_cyclic_method.c  
void FindEigensystem(double vxx, double vxy, double vxz, 
        double vyy, double vyz, double vzz, double eigenvalues[], double * eigenvectors)
{
   int n = 3; // num dimensions
   double A[9]; 

    A[0] = vxx; 
    A[1] = vxy; 
    A[2] = vxz; 
    A[3] = vxy; 
    A[4] = vyy;
    A[5] = vyz;
    A[6] = vxz;
    A[7] = vyz;
    A[8] = vzz; 

   int row, i, j, k, m;
   double *pAk, *pAm, *p_r; double *p_e;
   double threshold_norm;
   double threshold;
   double tan_phi, sin_phi, cos_phi, tan2_phi, sin2_phi, cos2_phi;
   double sin_2phi, cos_2phi, cot_2phi;
   double dum1;
   double dum2;
   double dum3;
   double r;
   double max;

          // Initialize the eigenvalues to the identity matrix.

   for (p_e = eigenvectors, i = 0; i < n; i++)
      for (j = 0; j < n; p_e++, j++)
         if (i == j) *p_e = 1.0; else *p_e = 0.0;
  
            // Calculate the threshold and threshold_norm.
 
   for (threshold = 0.0, pAk = A, i = 0; i < ( n - 1 ); pAk += n, i++) 
      for (j = i + 1; j < n; j++) threshold += *(pAk + j) * *(pAk + j);
   threshold = sqrt(threshold + threshold);
   threshold_norm = threshold * DBL_EPSILON;
   max = threshold + 1.0;
   while (threshold > threshold_norm) {
      threshold /= 10.0;
      if (max < threshold) continue;
      max = 0.0;
      for (pAk = A, k = 0; k < (n-1); pAk += n, k++) {
         for (pAm = pAk + n, m = k + 1; m < n; pAm += n, m++) {
            if ( fabs(*(pAk + m)) < threshold ) continue;

                 // Calculate the sin and cos of the rotation angle which
                 // annihilates A[k][m].

            cot_2phi = 0.5 * ( *(pAk + k) - *(pAm + m) ) / *(pAk + m);
            dum1 = sqrt( cot_2phi * cot_2phi + 1.0);
            if (cot_2phi < 0.0) dum1 = -dum1;
            tan_phi = -cot_2phi + dum1;
            tan2_phi = tan_phi * tan_phi;
            sin2_phi = tan2_phi / (1.0 + tan2_phi);
            cos2_phi = 1.0 - sin2_phi;
            sin_phi = sqrt(sin2_phi);
            if (tan_phi < 0.0) sin_phi = - sin_phi;
            cos_phi = sqrt(cos2_phi); 
            sin_2phi = 2.0 * sin_phi * cos_phi;
            cos_2phi = cos2_phi - sin2_phi;

                     // Rotate columns k and m for both the matrix A 
                     //     and the matrix of eigenvectors.

            p_r = A;
            dum1 = *(pAk + k);
            dum2 = *(pAm + m);
            dum3 = *(pAk + m);
            *(pAk + k) = dum1 * cos2_phi + dum2 * sin2_phi + dum3 * sin_2phi;
            *(pAm + m) = dum1 * sin2_phi + dum2 * cos2_phi - dum3 * sin_2phi;
            *(pAk + m) = 0.0;
            *(pAm + k) = 0.0;
            for (i = 0; i < n; p_r += n, i++) {
               if ( (i == k) || (i == m) ) continue;
               if ( i < k ) dum1 = *(p_r + k); else dum1 = *(pAk + i);
               if ( i < m ) dum2 = *(p_r + m); else dum2 = *(pAm + i);
               dum3 = dum1 * cos_phi + dum2 * sin_phi;
               if ( i < k ) *(p_r + k) = dum3; else *(pAk + i) = dum3;
               dum3 = - dum1 * sin_phi + dum2 * cos_phi;
               if ( i < m ) *(p_r + m) = dum3; else *(pAm + i) = dum3;
            }
            for (p_e = eigenvectors, i = 0; i < n; p_e += n, i++) {
               dum1 = *(p_e + k);
               dum2 = *(p_e + m);
               *(p_e + k) = dum1 * cos_phi + dum2 * sin_phi;
               *(p_e + m) = - dum1 * sin_phi + dum2 * cos_phi;
            }
         }
         for (i = 0; i < n; i++)
            if ( i == k ) continue;
            else if ( max < fabs(*(pAk + i))) max = fabs(*(pAk + i));
      }
   }
   for (pAk = A, k = 0; k < n; pAk += n, k++) eigenvalues[n-k-1] = (float) *(pAk + k); 

    // canned code spits these out in a different order than we want. 
    // someday I'll come back and do this properly, but for now
    // bless me father for I'm about to sin. 

    double e00, e01, e02, e10, e12, e20, e21, e22; 

    e00 = eigenvectors[2]; 
    e01 = eigenvectors[5];
    e02 = eigenvectors[8];
    e10 = eigenvectors[1];
    e12 = eigenvectors[7]; 
    e20 = eigenvectors[0]; 
    e21 = eigenvectors[1];
    e22 = eigenvectors[2]; 

    eigenvectors[0] = e00; 
    eigenvectors[1] = e01; 
    eigenvectors[2] = e02; 
    eigenvectors[3] = e10; 
    eigenvectors[5] = e12;
    eigenvectors[6] = e20; 
    eigenvectors[7] = e21;
    eigenvectors[8] = e22; 
}


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
    FindEigenvalues( in[0], in[3], in[4], in[1], in[5], in[2], out, major, minor);
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

1 1 2
1 3 3
2 3 5

./a.out 1 3 5 1 2 3
7.872983 1.000000 0.127017
Confirmed against Wolfram

w/ jacobi-cyclic addition: (FindEigensystem)
7.872983 1.000000 0.127017 (old method)
7.872983 1.000000 0.127017 (FindEigensystem)
     0.306461 0.543844 0.781227
     -0.408248 0.816497 -0.408248
     0.859893 -0.408248 0.306461

Wolfram gives eigenvalues:
{7.87298, 1., 0.127017}

and eigenvectors:
{{-0.306461, -0.543844, -0.781227}, 
{0.408248, -0.816497, 0.408248}, 
{-0.859893, -0.193823, 0.472247}}

*/

#endif 
