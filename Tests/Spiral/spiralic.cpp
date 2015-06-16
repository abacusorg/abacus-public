#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include "threevector.hh"
using namespace std;

#define TWOPI (2*M_PI)

#define MIN(a,b) ( (a)<(b) ? (a):(b) )

class Spiral {
public:

    Spiral(long long int N1D, double3 kvec, double3 phase); 
    ~Spiral();

    void Create(double3 *pp, double3 *vv, int *id, long long int np, double Ainitial, double Across);
    void Project(double3 *pp, double3 *vv, int *id, long long int np);

    void WrapPosition(double3 &p);
    double3 cross(double3 v1, double3 v2);
    double3 RotatePointAboutLine(double3 p, double3 n, double angle);

    long long int N1D, N;
    double3 kvec, phase;
};

Spiral::Spiral(long long int N1D, double3 kvec, double3 phase) : N1D(N1D), kvec(kvec), phase(phase) {
    N = N1D*N1D*N1D;
}

Spiral::~Spiral() {
}

// Create a one-d wave of perturbations with a given wave vector and phase.
// Positions are in dimensionless units and "velocity" is in specific congugate momentum
// (dx/dt) / a^2, in dimensionless units.
void Spiral::Create(double3 *pp, double3 *vv, int *id, long long int _n, double Ainitial, double Across) {

    assert(N==_n);

    printf("   N1D: %lld\n", N1D);
    printf("     N: %lld\n", N);
    printf("  kvec: % e  % e  % e\n", kvec.x, kvec.y, kvec.z);
    printf(" phase: % e  % e  % e\n", phase.x, phase.y, phase.z);
    printf(" Ainit: % e\n", Ainitial);
    printf("Across: % e\n", Across);
    
    kvec *= 2*M_PI; // we specify k-vector w/o the 2 pi for ease of writing...

    double       DofA = Ainitial;           // linearized perturnbation growth factor for Omega_M=1, K=0
    double    DdotofA = pow(Ainitial,1.5);  // time derivative of growth factor for Omega_M=1, K=0
    double  DofAcross = Across;             // growth factor at desired time of 1st crossing

    double       knrm = kvec.norm();
    double3 direction = kvec/knrm;
    double Amplitude = 1.0/(DofAcross*knrm);

    double      dx = 1.0/((double) N1D);  // spacing of unperturbed positions q
    int p = 0;
    for(int i=-N1D/2; i<N1D/2; i++) {
        for(int j=-N1D/2; j<N1D/2; j++) {
            for(int k=-N1D/2; k<N1D/2; k++) {
                
                double3 q = dx*double3(i,j,k); // unperturbed coordinate
                double kdotq = kvec.dot(q);
                double3 perturbation;
                perturbation.x = direction.x*Amplitude*sin(kdotq+phase.x);
                perturbation.y = direction.y*Amplitude*sin(kdotq+phase.y);
                perturbation.z = direction.z*Amplitude*sin(kdotq+phase.z);
                
                double3 pos = q + DofA*perturbation; WrapPosition(pos);
                double3 vel = DofA*perturbation;
                pp[p] = pos;
                vv[p] = vel;
                id[p] = p;
                p++;
            }
        }
    }
    assert(p==N);
}

void Spiral::WrapPosition(double3 &p) {
    if(p.x>=0.5) p.x -= 1.0; if(p.x<-0.5) p.x += 1.0;
    if(p.y>=0.5) p.y -= 1.0; if(p.y<-0.5) p.y += 1.0;
    if(p.z>=0.5) p.z -= 1.0; if(p.z<-0.5) p.z += 1.0;
}

double3 Spiral::cross(double3 v1, double3 v2) {
    double3 res;
    res.x = v1.y*v2.z - v1.z*v2.y;
    res.y = v1.z*v2.x - v1.x*v2.z;
    res.z = v1.x*v2.y - v1.y*v2.x;
    return res;
}

double3 Spiral::RotatePointAboutLine(double3 p, double3 n, double angle) {
    double norm = n.norm();
    n /= norm;
    double ca = cos(angle); double sa = sin(angle);
    double pdotn = p.dot(n);
    double3 pcrossn = p.cross(n);
    double3 pp;
    pp = ca*p + (1.0-ca)*pdotn*n + sa*pcrossn;
    return pp;
}

void Spiral::Project(double3 *pp, double3 *vv, int *id, long long int np) {

    // we want sol'm along kvec rotated parallel to xhat
    double3 xhat; xhat.x = 1; xhat.y = 0; xhat.z = 0;
    double xnorm = xhat.norm(); xhat /=xnorm;
    double3 khat; khat = kvec;
    double knorm = khat.norm(); khat /= knorm;
    double3 rotaxis = xhat.cross(kvec);
    double angle = acos( xhat.dot(khat) );
    double X = 1.0/knorm;

//    fprintf(fp,"#  Rotating solution %f degrees about ", angle*180/M_PI);
//    fprintf(fp," (% f, % f, % f)\n", rotaxis.x, rotaxis.y, rotaxis.z);
    printf("Rotating solution %f degrees about ", angle*180/M_PI);
    printf("(% f, % f, % f)\n", rotaxis.x, rotaxis.y, rotaxis.z);

    double3 pmin = double3(1e30);
    for(int p=0; p<N; p++) {

        double3 pos;
        pos = pp[p];
        if(fabs(angle)>0.0) pos = RotatePointAboutLine(pos, rotaxis, angle);

        // merge the multiple spirals due to wrapping
        if(pos.x >= X/2)  pos.x -= X; if(pos.x < -X/2) pos.x += X;
        if(pos.y >= X/2)  pos.y -= X; if(pos.y < -X/2) pos.y += X;
        if(pos.z >= X/2)  pos.z -= X; if(pos.z < -X/2) pos.z += X;

        double3 vel = vv[p];
        if(fabs(angle)>0.0) vel = RotatePointAboutLine(vel, rotaxis, angle);
        pp[p] = pos;
        vv[p] = vel;

        pmin.x = MIN(pmin.x,pos.x);
        pmin.y = MIN(pmin.y,pos.y);
        pmin.z = MIN(pmin.z,pos.z);
    }

    // put corner at -0.5 as per usual
    for(int p=0; p<N; p++) {
        pp[p] += - pmin - double3(0.5);
    }
}

void GenSpiral( long long int n1d, double ainitial, double across, double3 *pp, double3 *vv, int *id,double3 kvec,double3 phase) {

    long long int np = n1d*n1d*n1d;
    assert(np>0);



    Spiral SP(n1d, kvec, phase);
    SP.Create(pp, vv, id, np, ainitial, across);
}

void writespiral(char *fn, double3 * pos, double3 * vel, long long int np){
	FILE *outfile = fopen(fn,"wb");
	for(long long int i = 0; i < np; i++){
		fwrite(pos+i,1,sizeof(double3),outfile);
		fwrite(vel+i,1,sizeof(double3),outfile);
	}
	fclose(outfile);
}

int main(int argc, char **argv){

	if (argc != 11) {
		printf("Usage: makespiralic <n1d> <a initial> <across> <kvec x y z> <phase x y z> <output filename> ");
		return 1;
	}
	long long int n1d = atoi(argv[1]);
	double3 kvec,phase;
	double ainitial = atof(argv[2]);
	double across = atof(argv[3]);
	kvec.x = atof(argv[4]);
	kvec.y = atof(argv[5]);
	kvec.z = atof(argv[6]);
	phase.x = atof(argv[7]);
	phase.y = atof(argv[8]);
	phase.z = atof(argv[9]);
	long long int np = n1d*n1d*n1d;
	double3 * pos = new double3[np];
	double3 * vel = new double3[np];
	int *id = new int[np];

	GenSpiral(n1d, ainitial, across,pos,vel,id,kvec,phase);
	writespiral(argv[10],pos,vel,np);

	return 0;
}
