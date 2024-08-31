class Order32Derivatives;

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
#define MAXLOOPORDER 32 //Recommended value for runtime/accuracy balance: 24. 

Order32Derivatives::Order32Derivatives(int inner) : inner(inner)  {
    for(int a=0;a<ODIM;a++)
       for(int b=0;b<ODIM;b++)
            for(int c=0;c<ODIM;c++)
                AnalyticDerivatives[a][b][c] = 0;

    int tmporder, tmpinneradius;
    fs::path fname = fmt::format("AD32_{:03d}.dat",inner);
    fmt::print("reading derivatives from {}\n", fname);
    ifstream fin(fname,ios::in);
    if(!fin.is_open()) {
        fmt::print(stderr,"couldn't open \"{}\"\n",fname);
        exit(1);
    }
    fin >> tmporder >> tmpinneradius;
    if(tmporder<32 || tmpinneradius!=inner) {
        fmt::print(stderr,"file data does not agree with arguments\n");
        exit(1);
    }
    char num[1024];
    for(int a=0; a<=32; a+=2) {
        for(int b=a; b<=32-a; b+=2) {
            int atmp, btmp;
            fin >> atmp >> btmp >> num;
            assert(a==atmp && b==btmp);
            AnalyticDerivatives[a][b][0] = qd_real(num);  //read in K = 1 derivatives tensor \sum_{\mathcal{B}\left(\mathbf{n}, L_{\mathrm{outer}}, \infty\right)} D^{(a+d)(b+e)(c+f)}(\mathbf{n}) from AD32_innerRadius.dat file. We need this for Eqn. 4.8. Note that AD32_innerRadius.dat contains only the scalars we need to reconstruct the derivatives tensor. 
            AnalyticDerivatives[b][a][0] = qd_real(num);
        }
    }
    fin.close();

    for(int c=2;c<=32;c+=2) {
        for(int a=0;a<=32-c;a++) {
            for(int b=0;b<=32-a-c;b++) {
                AnalyticDerivatives[a][b][c] = //Reconstruct the rest of the derivatives with third index != 0,1 using the trace-free relationship, Eqn. 2.24.
                    -AnalyticDerivatives[a+2][b][c-2]
                    -AnalyticDerivatives[a][b+2][c-2];
            }
        }
    }

}

//Calculate the inner far field region's derivatives tensor \sum\limits_{\mathcal{B}\left(\mathbf{n}, 1, L_{\mathrm{outer}}\right)} D^{ABC}(\mathbf{n} + \mathbf{c}_{jkl}) explicitly by lopping over all lattice vectors \mathbf{n} that fall between L_{inner} (inner_radius) and L_{outer} (far_radius)
#define FORALL_RADII(n,r)  for(n.x=-r;n.x<=r;n.x++) for(n.y=-r;n.y<=r;n.y++) for(n.z=-r;n.z<=r;n.z++)
void Order32Derivatives::SumNearDerivatives(double3 r, int radius, double *rd, int order) {
    integer3 n;


    basemultipoles bm(order);
    int rml = (order+1)*(order+1);

    double d[rml];

    Derivatives RD(order);

    for(int m=0;m<rml;m++) rd[m] = 0;

    FORALL_RADII(n, radius) { //for some radius int radius, do all values of lattice vector n that gives norm(n) < radius
        if( (n.norm2() > 0) && (n.norm2() <= radius*radius) ) { //this makes sure we're outside the home box and that we're within L_{outer} (far_radius) 
            double3 x = r + n;
            RD.ReducedDerivatives(x,d); //Use the recursion and trace free relationships discussed in Section 2.2. 
            for(int m=0;m<rml;m++) rd[m] += d[m];
        }
    }
}

//Compute the outer far field's contribution to the total derivatives tensor by using the K=1 toy problem trick described in Section 4. 
void Order32Derivatives::InfiniteDerivativeSummation(double3 r, double *rd, int order) {
	int maxLoopOrder = MAXLOOPORDER; 
	
    qd_real rx[maxLoopOrder+1],ry[maxLoopOrder+1],rz[maxLoopOrder+1]; 

    rx[0] = 1; rx[1] = to_qd_real(r.x);
    ry[0] = 1; ry[1] = to_qd_real(r.y);
    rz[0] = 1; rz[1] = to_qd_real(r.z);
		
    for(int i=2;i<=maxLoopOrder;i++) {
		rx[i] = rx[i-1] * rx[1] /i; //Referring to Eqn. 4.8, this equals \mathbf{c}_{jkl,x}^a/a!
		ry[i] = ry[i-1] * ry[1] /i; // 									 \mathbf{c}_{jkl,y}^b/b!
		rz[i] = rz[i-1] * rz[1] /i; // 									 \mathbf{c}_{jkl,z}^c/c!
    }

    int a,b,c; //Eqn. 4.8's a, b, c.
    int d,e,f; //Eqn. 4.8's d, e, f. 

    int m = 0;
	
	for(int f = 0; f <= 1; f++){
		for(int d = 0; d <= order - f; d++){
			for(int e = 0; e <= order - f - d; e++) //for a fixed d, e, f, compute \sum_{\mathcal{B}\left(\mathbf{n}, L_{\mathrm{outer}}, \infty\right)}D^{def}(\mathbf{n} + \mathbf{c}_{jkl}) by looping over a, b, c to some maximum loop order p (which can be higher than the rest of the calculation's order, if desired). 
			{
				qd_real s = 0;
				
				for(int a = d%2; a <= maxLoopOrder; a+=2)
				{
					for(int b = e%2; b <= maxLoopOrder-a; b+=2)
					{
						int cmax = maxLoopOrder - d - e - f - a - b;
						for(int c = f%2;c <= cmax; c+=2)
						{
							qd_real x = rx[a] * ry[b] * rz[c]  * AnalyticDerivatives[a+d][b+e][c+f]; //Eqn. 4.8
							s += x; 
						 }
					}	
				}
						
				rd[m] = to_double(s); 
				m++;  
						
			}
		}
	}
}

//Compute the far field derivatives -- specifically, the first two terms of Equation 4.1.
void Order32Derivatives::Derivative(double3 r, double *rdfar, double *rdnear, int order) {
	
	//Compute inner far field's contribution to the derivatives tensor by explicitly looping over all \mathbf{c}_{jkl} between L_{inner} and L_{outer} (but outside of the home box -- that portion of the inner far field was dealt with in CreateDerivatives.cpp's FormDerivatives) and using the trace-free and recursion relations discussed in Section 2.2 
    SumNearDerivatives( r, inner , rdnear, order);
	//Compute the outer far field's contribution to the total derivatives tensor by using the K=1 toy problem trick described in Section 4 and given by Eqn. 4.8. 
    InfiniteDerivativeSummation(r, rdfar, order);

}
