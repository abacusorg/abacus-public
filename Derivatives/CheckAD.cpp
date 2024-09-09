#include "../include/header.cpp"
#include "quad_double.cpp"
#include "../include/threevector.hh"

typedef ThreeVector<qd_real> qd_real3;


#define N 4


#define ODIM (32+1)

#define FORALLREPLICAS(n,r) \
            for(n.x=-r;n.x<=r;n.x++) \
                for(n.y=-r;n.y<=r;n.y++) \
                    for(n.z=-r;n.z<=r;n.z++) 


qd_real fact(int n) {
    assert(n>=0);
    qd_real f=1;
    for(int i=2;i<=n;i++) f*=i;
    return f;
}

// (-1)^n as an integer
int powofmone(int n) {
    int s=1;
    for(int i=0; i<n; i++) s *= -1;
    return s;
}

// Read in derivatives from file
// (Made with MakeDerivatives.cpp)
qd_real AnalyticDerivatives[ODIM][ODIM][ODIM];
void MakeAnalyticDerivatives(int maxorder, int innerradius) {
    assert(maxorder < ODIM);

    for(int a=0;a<ODIM;a++) 
        for(int b=0;b<ODIM;b++) 
            for(int c=0;c<ODIM;c++) 
                AnalyticDerivatives[a][b][c] = 0;

    int tmporder, tmpinneradius;
    fs::path fname = fmt::format("AD32_{:03d}.dat", innerradius);
    fmt::print("reading derivatives from {:s}\n", fname);
    ifstream fin(fname,ios::in);
    if(!fin.is_open()) {
        fmt::print(stderr,"couldn't open \"{:s}\"\n",fname);
        exit(1);
    }
    fin >> tmporder >> tmpinneradius;
    if(tmporder<maxorder || tmpinneradius!=innerradius) {
        fmt::print(stderr,"file data does not agree with arguments\n");
        exit(1);
    }
    std::string num;
    for(int a=0; a<=32; a+=2) {
        for(int b=a; b<=32-a; b+=2) {
            int atmp, btmp;
            fin >> atmp >> btmp >> num;
//                   fmt::print("READINGING {:3d} {:3d}  {:s}\n", atmp, btmp, num);
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

    for(int c=2;c<=32;c+=2) {
        for(int a=0;a<=32-c;a++) {
            for(int b=0;b<=32-a-c;b++) {
                if( (AnalyticDerivatives[a][b][c]  > qd_real("0")) || (AnalyticDerivatives[a][b][c] < qd_real("0")) ) { 
                    // fmt::print("{:d} {:d} {:d}   ", a,b,c);
                    // fmt::print("{}\n", AnalyticDerivatives[a][b][c]);
                }
            }
        }
    }


    fmt::print("Done making derivatives \n");
}

qd_real PI = qd_real::_pi;
qd_real SPI = sqrt(qd_real::_pi);
qd_real PI2 = 2*qd_real::_pi;

// Ewald must go out to a radius of 8 for 63 digits
qd_real3 ewaldacc( qd_real3 ri,  qd_real3 rj) {
    qd_real r2,nr,h2;
    qd_real3 r,dr;
    integer3 n,h;
    qd_real3 ewald;
    ewald.zero();
    r=ri-rj;
    FORALLREPLICAS(n,8) {
        dr=r-n;
        r2 = dr.norm2(); nr = sqrt(r2);
        ewald += (erfc(SPI*nr) + 2*nr * exp(-PI*r2))/(r2*nr) * dr;
    }
    FORALLREPLICAS(h,8) {
        h2 = h.x*h.x + h.y*h.y + h.z*h.z;
        if(h2 > qd_real("0") ) ewald  += 2/h2 * exp(-PI*h2) * sin(PI2*h.dot(r)) * h ;
    }
    return ewald;
}

qd_real3 pos[N];
qd_real mass[N];

qd_real M[ODIM][ODIM][ODIM];

void MakePoints(void) {
#if 0
    srand48(234234234);

    for(int i=0;i<N;i++) {
        pos[i].x = to_qd_real(drand48()-0.5);
        pos[i].y = to_qd_real(drand48()-0.5);
        pos[i].z = to_qd_real(drand48()-0.5);
        mass[i] = to_qd_real(drand48());
    }
#else
    assert(N==4);
    pos[0].x =  qd_real(1)/3;
    pos[0].y =  qd_real(1)/3;
    pos[0].z =  qd_real(1)/3;
    mass[0] = 1;
    pos[1].x = -qd_real(1)/3;
    pos[1].y =  qd_real(1)/3;
    pos[1].z =  qd_real(1)/3;
    mass[1] = 1;
    pos[2].x = -qd_real(1)/3;
    pos[2].y = -qd_real(1)/3;
    pos[2].z =  qd_real(1)/3;
    mass[2] = 1;
    pos[3].x = -qd_real(1)/3;
    pos[3].y = -qd_real(1)/3;
    pos[3].z = -qd_real(1)/3;
    mass[3] = 1;
#endif

}

void MakeMultipolesAcc(int order) {
    
    for(int i=0;i<=order;i++) 
        for(int j=0;j<=order-i;j++) 
            for(int k=0;k<=order-i-j;k++) 
                M[i][j][k] = 0;

    for(int i=0;i<=order;i++) 
        for(int j=0;j<=order-i;j++) 
            for(int k=0;k<=order-i-j;k++) 
                for(int p=0;p<N;p++) 
                    if(p!=100) {
                        qd_real3 q;
                        q=pos[p];
                        qd_real term = mass[p] * pow(q.x,i) * pow(q.y,j) * pow(q.z,k);
                        M[i][j][k] += term;
                    }

    for(int i=0;i<=order;i++) 
        for(int j=0;j<=order-i;j++) 
            for(int k=0;k<=order-i-j;k++) 
                M[i][j][k] /= (fact(i) * fact(j) * fact(k)); 

}



#define FORALL_RADII(n,r)  for(n.x=-r;n.x<=r;n.x++) for(n.y=-r;n.y<=r;n.y++) for(n.z=-r;n.z<=r;n.z++) 

#define FOR(a,b,c) for(a=b;a<=c;a++)
#define FORALL_COMPLETE_MULTIPOLES_BOUND(a,b,c,u)  FOR(a,0,u) FOR(b,0,u-a) FOR(c,0,u-a-b)





qd_real T[ODIM][ODIM][ODIM];

void MakeTaylor(int order) {

    for(int k=0;k<=1;k++)     
        for(int i=0;i<=order-k;i++) 
            for(int j=0;j<=order-k-i;j++)  {
                qd_real s = 0;
                for(int a=0;a<=order-(i+j+k);a++) 
                    for(int b=0;b<=order-(i+j+k)-a;b++) 
                        for(int c=0;c<=order-(i+j+k)-a-b;c++) {
                            s += M[a][b][c] *
                                 AnalyticDerivatives[a+i][b+j][c+k];
                        }
                 T[i][j][k] = s;
            }

    for(int k=2;k<=order;k++)  
        for(int i=0;i<=order-k;i++)
            for(int j=0;j<=order-k-i;j++)  
                T[i][j][k] = -T[i+2][j][k-2] -T[i][j+2][k-2];

    for(int i=0;i<=order;i++)
        for(int j=0;j<=order-i;j++)  
            for(int k=0;k<=order-i-j;k++) 
                T[i][j][k] *= qd_real(powofmone(i+j+k)) /
                    (fact(i)*fact(j)*fact(k));
}

qd_real3 EvaluateTaylorAcc(int order, qd_real3 rp) {
    qd_real fi,fij,fijk;

    qd_real3 a;

    a.zero(); //intializes three-vector to zeros (see /include/threevector.hh)

    fi = 1;
    for(int i=0;i<=order-1;i++) {
        fij = fi;
        for(int j=0;j<=order-1-i;j++) {
            fijk = fij;
            for(int k=0;k<=order-1-i-j;k++) {
                a.x += qd_real(i+1) * T[i+1][j][k] * fijk;
                a.y += qd_real(j+1) * T[i][j+1][k] * fijk;
                a.z += qd_real(k+1) * T[i][j][k+1] * fijk;
                fijk *= rp.z;
            }
            fij *= rp.y;
        }
        fi *= rp.x;
    }

    return a;
}


qd_real3 acc( qd_real3 r) {
    qd_real ir = 1/r.norm();
    qd_real ir3 = ir*ir*ir;
    qd_real3 a;
    a=ir3*r;
    return a;
}


qd_real3 latticeacc(int radius,   qd_real3 r) {
    integer3 n;
    qd_real3 a;
    
    a.zero(); //intializes three-vector to zeros (see /include/threevector.hh)
    FORALL_RADII(n,radius) {
        if(n.norm2() <=radius*radius) {
            qd_real3 x;
            x=r+n;
            a+=acc(x);
        }
    }
    return a;
}

void TestAcc(int order, int innerradius) {
    MakeAnalyticDerivatives(order, innerradius);
    MakePoints();

    

    fmt::print("Making Multipoles \n");
    MakeMultipolesAcc(order);
    fmt::print("Done Making Multipoles \n");

    fmt::print("Making Taylor \n");
    MakeTaylor(order);
    fmt::print("End Making Taylor\n");

    
    for(int j=0;j<N;j++) { 

        qd_real3 ewald;
        ewald.zero();
        for(int i=0;i<N;i++) {
            if(i!=j) {
                ewald+=mass[i]*ewaldacc(pos[i],pos[j]); //true acc
            }
        }

        qd_real3 l;
        l.zero();
        for(int i=0;i<N;i++) {
            if(i!=j) {
                qd_real3 r;
                r=pos[i]-pos[j];
                l+=mass[i]*latticeacc(innerradius,r); //near acc
            }
        }
        
        qd_real3 red;
        red.zero();
        
        for(int i=0;i<N;i++) {
            if(i!=j) {
                qd_real3 r;
                r=pos[i]-pos[j];
                red+= qd_real(4)/qd_real(3)*PI*mass[i]*r;     //redlack-grindlay term
            }
        }
        
        qd_real3 t;
        t = EvaluateTaylorAcc(order, pos[j]); //far field acc
        
        qd_real3 ans;
        ans=l-red+t; //near acc - redlack grindlay + far field acceleration
        qd_real3 err;
        err=ans-ewald; //abs error

        qd_real3 re;
        re=err/ewald.norm(); //relative error

#if 0        
        fmt::print("ans: \n");
        fmt::print("  {}  {}  {}\n", ans.x, ans.y, ans.z);
        fmt::print("ewald: \n");
        fmt::print("  {}  {}  {}\n", ewald.x, ewald.y, ewald.z);
        fmt::print("rel err: \n");
        fmt::print("  {}  {}  {}\n", re.x, re.y, re.z);
        fmt::print("\n");
#else
        fmt::print("relerr = {: e}  {: e}  {: e}\n",
               (double)re.x, (double)re.y, (double)re.z);
#endif
    }
    fmt::print("\n");
}


int main(void) {
   unsigned int old_cw;
    fpu_fix_start(&old_cw);

    cout.precision(qd_real::_ndigits);
//     int nd = qd_real::_ndigits + 8;

    for(int innerradius=1; innerradius<=9; innerradius+=1) {
		if(innerradius==9){innerradius = 16;};
        fmt::print("Acc:\n");
        TestAcc(32,innerradius);
    }

    fpu_fix_end(&old_cw);
}
