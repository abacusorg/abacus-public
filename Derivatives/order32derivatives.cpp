class Order32Derivatives;
#include "../include/PTimer.cc"
//double Runtime = 0.0;
//double numLoops = 0.0;
int OPTIMIZE = 1;

class Order32Derivatives {
public:
    Order32Derivatives(int inner);
    void SumNearDerivatives(double3 r, int radius, double *rd, int order);
    void InfiniteDerivativeSummation(double3 r, double *rd, int order);
    void Derivative(double3 r, double *rdnear, double *rdfar, int order);

private:
    qd_real AnalyticDerivatives[33][33][33];
    int inner;
};


#define ODIM 33


Order32Derivatives::Order32Derivatives(int inner) : inner(inner)  {
    for(int a=0;a<ODIM;a++)
       for(int b=0;b<ODIM;b++)
            for(int c=0;c<ODIM;c++)
                AnalyticDerivatives[a][b][c] = 0;

    int tmporder, tmpinneradius;
    char fname[1024];
    sprintf(fname,"AD32_%03d.dat",inner);
    printf("reading derivatives from %s\n", fname);
    ifstream fin(fname,ios::in);
    if(!fin.is_open()) {
        fprintf(stderr,"couldn't open \"%s\"\n",fname);
        exit(1);
    }
    fin >> tmporder >> tmpinneradius;
    if(tmporder<32 || tmpinneradius!=inner) {
        fprintf(stderr,"file data does not agree with arguments\n");
        exit(1);
    }
    char num[1024];
    for(int a=0; a<=32; a+=2) {
        for(int b=a; b<=32-a; b+=2) {
            int atmp, btmp;
            fin >> atmp >> btmp >> num;
//                   printf("READINGING %3d %3d  %s\n", atmp, btmp, num);
            assert(a==atmp && b==btmp);
            AnalyticDerivatives[a][b][0] = qd_real(num);
            AnalyticDerivatives[b][a][0] = qd_real(num);
        }
    }
    fin.close();

    for(int c=2;c<=32;c+=2) {
        for(int a=0;a<=32-c;a++) {
            for(int b=0;b<=32-a-c;b++) {
                AnalyticDerivatives[a][b][c] =
                    -AnalyticDerivatives[a+2][b][c-2]
                    -AnalyticDerivatives[a][b+2][c-2];
            }
        }
    }

}


#define FORALL_RADII(n,r)  for(n.x=-r;n.x<=r;n.x++) for(n.y=-r;n.y<=r;n.y++) for(n.z=-r;n.z<=r;n.z++)
void Order32Derivatives::SumNearDerivatives(double3 r, int radius, double *rd, int order) {
    integer3 n;


    basemultipoles bm(order);
    int rml = (order+1)*(order+1);

    double d[rml];

    Derivatives RD(order);


    for(int m=0;m<rml;m++) rd[m] = 0;

    FORALL_RADII(n, radius) {
        if( (n.norm2() > 0) && (n.norm2() <= radius*radius) ) {
            double3 x = r + n;
            RD.ReducedDerivatives(x,d);
            for(int m=0;m<rml;m++) rd[m] += d[m];
        }
    }
}


void Order32Derivatives::InfiniteDerivativeSummation(double3 r, double *rd, int order) {
	if( !OPTIMIZE ){
		 ////////////////////////////////////////////////////////////////////////
		//NO OPTIMIZATION (factorial computed each loop, multipole bound = 32)
		qd_real fmem[33];  
		fmem[0] = 1;
		fmem[1] = 1;
		for(int i=2;i<=32;i++) { 
			fmem[i] = fmem[i-1] * i;
		}		
	
	    qd_real rx[33],ry[33],rz[33]; 
	
	    rx[0] = 1; rx[1] = to_qd_real(r.x);
	    ry[0] = 1; ry[1] = to_qd_real(r.y);
	    rz[0] = 1; rz[1] = to_qd_real(r.z);
	
	    for(int i=2;i<=32;i++) { 
			rx[i] = rx[i-1] * rx[1] ; 
			ry[i] = ry[i-1] * ry[1] ;
			rz[i] = rz[i-1] * rz[1] ;
	    }
	
	    int i,j,k;
	    int a,b,c;
    
	    int m = 0;
	    FORALL_REDUCED_MULTIPOLES_BOUND(a,b,c,order) {
	        qd_real s = 0;
	        FORALL_COMPLETE_MULTIPOLES_BOUND(i,j,k,32) { 
	            if( a+i+b+j+c+k <= 32 )  {
	                if( ((a+i)%2==0) && ((b+j)%2==0) && ((c+k)%2==0) ) {
	                    qd_real x = rx[i] * ry[j] * rz[k] / (fmem[i] * fmem[j] * fmem[k]) * AnalyticDerivatives[a+i][b+j][c+k];
	                    s += x;
	                }
	            }
	        }
	        rd[m] = to_double(s);
	        m++;        
	    }
		////////////////////////////////////////////////////////////////////////
	}
	
	else{
		////////////////////////////////////////////////////////////////////////
		//OPTIMIZED VERSION (factorial fixed, multipole bound = 24)
	
	    qd_real rx[9],ry[9],rz[9]; 

	    rx[0] = 1; rx[1] = to_qd_real(r.x);
	    ry[0] = 1; ry[1] = to_qd_real(r.y);
	    rz[0] = 1; rz[1] = to_qd_real(r.z);
	
	    for(int i=2;i<=8;i++) {
			rx[i] = rx[i-1] * rx[1] /i; 
			ry[i] = ry[i-1] * ry[1] /i;
			rz[i] = rz[i-1] * rz[1] /i;
	    }
	
	    int i,j,k;
	    int a,b,c;
    
	    int m = 0;
	    FORALL_REDUCED_MULTIPOLES_BOUND(a,b,c,order) {
	        qd_real s = 0;
	        FORALL_COMPLETE_MULTIPOLES_BOUND(i,j,k,8) { 
	            if( a+i+b+j+c+k <= 8 )  {//change 32 to 24 // precompute kmax = 32 - a-i-b-j-c for innermost loop and remove if statement //--> k< kmax 
	                if( ((a+i)%2==0) && ((b+j)%2==0) && ((c+k)%2==0) ) {
						qd_real x = rx[i] * ry[j] * rz[k]  * AnalyticDerivatives[a+i][b+j][c+k];
	                    s += x;
	                }
	            }
	        }
	        rd[m] = to_double(s);
	        m++;        
	    }
		////////////////////////////////////////////////////////////////////////
	}

}

void Order32Derivatives::Derivative(double3 r, double *rdfar, double *rdnear, int order) {
	//PTimer TIMER;
	//TIMER.Clear();
	//TIMER.Start();
    SumNearDerivatives( r, inner , rdnear, order);
	//TIMER.Stop();
	//double timeElapsed = TIMER.Elapsed();
	//printf("%f\n", timeElapsed); 
	//Runtime = Runtime + timeElapsed;
	//numLoops = numLoops + 1.0;
	//printf("numLoops = %f, timeElapsed = %f, average runtime = %f\n", numLoops, Runtime, Runtime/numLoops);
	//TIMER.Clear();
	
	//PTimer TIMER;
	//TIMER.Clear();
	//TIMER.Start();
	
    InfiniteDerivativeSummation(r, rdfar, order);
	
	//TIMER.Stop();
	//double timeElapsed = TIMER.Elapsed();
	//printf("%f\n", timeElapsed); 
	//Runtime = Runtime + timeElapsed;
	//numLoops = numLoops + 1.0;
	//printf("numLoops = %f, timeElapsed = %f, average runtime = %f\n", numLoops, Runtime, Runtime/numLoops);
	//TIMER.Clear();
}
