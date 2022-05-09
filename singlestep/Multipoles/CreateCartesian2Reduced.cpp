#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "threevector.hh"
#include "basemultipoles.h"
#include "basemultipoles.cpp"


#define LLI long long int 

#define FOR(a,b,c) for(a=b;a<=c;a++) 
#define FORALL_COMPLETE_MULTIPOLES_BOUND(a,b,c,u) \
            FOR(a,0,u) FOR(b,0,u-a) FOR(c,0,u-a-b)
#define FORALL_REDUCED_MULTIPOLES_BOUND(a,b,c,u) \
            FOR(c,0,1) FOR(a,0,u-c) FOR(b,0,u-c-a) 
#define TRACEFREERECURSION(a,b,c,u) \
            FOR(c,2,u) FOR(a,0,u-c) FOR(b,0,u-c-a)
#define FORALL_REDUCED_MULTIPOLES_ORDER(a,b,c,p) \
            FORALL_REDUCED_MULTIPOLES_BOUND(a,b,c,p) if(a+b+c==p) 


inline int lmap(int a, int b, int c) {
    return a*(MAXORDER+1)*(MAXORDER+1)+b*(MAXORDER+1)+c;
}

int    _cmap[(MAXORDER+1)*(MAXORDER+1)*(MAXORDER+1)];
int    _rmap[(MAXORDER+1)*(MAXORDER+1)*(MAXORDER+1)];

typedef struct { int a,b,c; } tuple3;

tuple3 _invcmap[(MAXORDER+1)*(MAXORDER+1)*(MAXORDER+1)];
tuple3 _invrmap[(MAXORDER+1)*(MAXORDER+1)];

inline int cmap(int a, int b, int c) { return _cmap[lmap(a,b,c)]; }
inline int rmap(int a, int b, int c) { return _rmap[lmap(a,b,c)]; }

tuple3 invcmap(int n) { return _invcmap[n]; }
tuple3 invrmap(int n) { return _invrmap[n]; }

int completemultipolelength;
int reducedmultipolelength;

int ncomplete(int order) { return ((order+1)*(order+2))/2; }
int nreduced(int order) { return (2*order+1); }

void fillmultipoleindices(int order) {

    completemultipolelength = (((order+1)*(order+2)*(order+3))/6);
    reducedmultipolelength = (order+1)*(order+1);

    int m = 0;
    int a,b,c;
    FORALL_COMPLETE_MULTIPOLES_BOUND(a,b,c,order) {
        _cmap[ lmap(a,b,c) ] = m;
        _invcmap[ m ].a = a;
        _invcmap[ m ].b = b;
        _invcmap[ m ].c = c;
        m++;
    }

    m = 0;
    FORALL_REDUCED_MULTIPOLES_BOUND(a,b,c,order) {
        _rmap[ lmap(a,b,c) ] = m;
        _invrmap[m].a = a;
        _invrmap[m].b = b;
        _invrmap[m].c = c;
        m++;
    }
}


#define ULLI unsigned long long int 
ULLI factorial(int n) {
    static ULLI  _factorial[MAXORDER+1] = { 0 };
    assert(n>=0);
    assert(n<=MAXORDER);
    if(_factorial[n]==0) {
        ULLI f=1;
        for(int i=2;i<=n;i++) f*=i;
        _factorial[n] = f;
    }
    return _factorial[n];
}

ULLI doublefactorial(int n) {
    static ULLI  _doublefactorial[2*MAXORDER+1] = { 0 };
    assert(n>=0);
    assert(n<=2*MAXORDER);
    if(_doublefactorial[n]==0) {
        ULLI df = 1;
        for(int i=n;i>1;i-=2) df *= i;
        assert(df>=1);
        _doublefactorial[n] = df;
    }
    return _doublefactorial[n];
}

ULLI zeta(int n, int m) {
    static ULLI _zeta[MAXORDER+1][MAXORDER+1] = { {0} };
    assert(n>=0); 
    assert(n<=MAXORDER);
    assert(m>=0);
    assert(m<=MAXORDER);
    assert(n>=m);
    assert(n>=2*m);

    if(_zeta[n][m]==0) {
        ULLI s = 0;
        s =  factorial(n);
        s /= doublefactorial(2*m); 
            assert(s*doublefactorial(2*m)==factorial(n));
        s /= factorial(n-2*m); 
            assert(s * factorial(n-2*m) * doublefactorial(2*m) == factorial(n));
        assert(s>=1);
        _zeta[n][m] = s;
    }
    return _zeta[n][m];
}

ULLI combination(int n, int m) {
    static ULLI _combination[MAXORDER+1][MAXORDER+1] = { {0} };
    assert(n>=0);
    assert(n<=MAXORDER);
    assert(m>=0);
    assert(m<=MAXORDER);
    assert(n>=m);
    if(_combination[n][m]==0) {
        ULLI s;
        s = factorial(n);
        s /= factorial(m); 
            assert(s*factorial(m)==factorial(n));
        s /= factorial(n-m); 
            assert(s*factorial(m)*factorial(n-m)==factorial(n));
        assert(s>=1);
        _combination[n][m] = s;
    }
    return _combination[n][m];
}

struct Polynomial {
    int termcount;
    LLI C[1023];

    Polynomial(void) {
        termcount = 0;
        for(int i=0;i<1023;i++) C[i] = 0;
    }
};

LLI ff[1024];
Polynomial rp[1024];


Polynomial ReducedPolynomial(int a, int b, int c) {

    assert(a>=0);
    assert(a<=MAXORDER);
    assert(b>=0);
    assert(b<=MAXORDER);
    assert(c>=0);
    assert(c<=1);
    assert(a+b+c<=MAXORDER);

    Polynomial p;

    for(int l=0;l<=a/2;l++) {
        for(int m=0;m<=b/2;m++) {
            for(int n=0;n<=c/2;n++) {
                LLI abcsign = 1; if( (a+b+c)%2 == 1) abcsign = -1;
                LLI lmnsign = 1; if( (l+m+n)%2 == 1) lmnsign = -1;
                LLI f = abcsign * lmnsign *
                            doublefactorial( 2*(a+b+c) - 2*(l+m+n) - 1 ) *
                            zeta(a,l) *
                            zeta(b,m) *
                            zeta(c,n);

                for(int i=0;i<=l+m+n;i++) {
                    for(int j=0;j<=l+m+n-i;j++) {
                        LLI f2 = f * combination(l+m+n,i) * combination(l+m+n-i,j);

                        int xp = 2*i + a-2*l;
                        int yp = 2*j + b-2*m;
                        int zp = 2*l + 2*m - 2*i - 2*j  + c;
                        assert(xp>=0);
                        assert(yp>=0);
                        assert(zp>=0);

                        int q = cmap(xp,yp,zp);
                        p.C[q] += f2;
                    }
                }
            }
        }
    }

    int termcount = 0;
    for(int i=0;i<completemultipolelength;i++) if(p.C[i]!=0) termcount++;
    p.termcount = termcount;
    return p;
}

Polynomial ConstMultiply(Polynomial pp, LLI f) {
    Polynomial p;

    for(int q = 0;q<completemultipolelength;q++) {
        if( pp.C[q] != 0) {
            p.C[q] = pp.C[q] * f;
        }
    }
    
    p.termcount = pp.termcount;
    return p;
}

void DividePoly( Polynomial &p, LLI x) {
    for(int i=0;i<completemultipolelength;i++) {
        LLI q = p.C[i];
        if(q!=0) {
            LLI qq;
            qq = q/x;
            assert(qq*x == q );
            p.C[i] = qq;
        }
    }
}

int DoesDivide( Polynomial p, LLI x ) {
    for(int i=0;i<completemultipolelength;i++) {
        LLI q = p.C[i];
        if(q!=0) {
            LLI qq;
            qq = q/x;
            if(qq*x != q ) return 0;
        }
    }
    return 1;
}

void MaxDividePoly( Polynomial &p, LLI *factor ) {
    LLI f = 1;
    LLI ff;
    
    ff=3;
    while( DoesDivide( p, ff ) ) {
        f = f * ff;
        DividePoly(p,ff);
    }

    ff=5;
     while( DoesDivide( p, ff ) ) {
        f = f * ff;
        DividePoly(p,ff);
    }

    ff=7;
     while( DoesDivide( p, ff ) ) {
        f = f * ff;
        DividePoly(p,ff);
    }

    ff=11;
     while( DoesDivide( p, ff ) ) {
        f = f * ff;
        DividePoly(p,ff);
    }

    ff=13;
     while( DoesDivide( p, ff ) ) {
        f = f * ff;
        DividePoly(p,ff);
    }

    *factor = f;
}


void DumpCartesian2Reduced(FILE *fp, int order, int a, int b, int c) {
    assert(c==0 || c==1);
    fprintf(fp,"double Cartesian2Reduced%d_%d_%d_%d(double *cm) { \n",order,a,b,c);
    fprintf(fp,"double s = 0; \n");
    for(int i=0;i<completemultipolelength;i++) {
        LLI p = rp[rmap(a,b,c)].C[i];
        if(p!=0) {
            fprintf(fp,"s += %lld * cm[%d];\n", p, i);
        }
    }
    fprintf(fp,"s *= %lld;\n", ff[rmap(a,b,c)]);
    fprintf(fp,"return s;\n");
    fprintf(fp,"}\n");
}

void CreateC2R(int order) {
    assert( order >= 0);
    assert( order <= MAXORDER );

    fillmultipoleindices(order);

    FILE *fp;   
    char fn[1024];
    sprintf(fn,"Cartesian2Reduced%d.cpp",order);
    fp = fopen(fn,"w");
    assert(fp!=NULL);

    int a,b,c;
    FORALL_REDUCED_MULTIPOLES_BOUND(a,b,c,order) {
        if(a+b+c>0) {
            rp[ rmap(a,b,c) ] = ReducedPolynomial(a,b,c);
            MaxDividePoly( rp[ rmap(a,b,c) ] , &ff[rmap(a,b,c)]); 
        }
    }

    fprintf(fp,"double Cartesian2Reduced%d_0_0_0(double *cm) { return cm[0]; }\n",order);
    FORALL_REDUCED_MULTIPOLES_BOUND(a,b,c,order) {
        if(a+b+c>0) {
            DumpCartesian2Reduced(fp,order,a,b,c);
        }
    }

    fprintf(fp, "\ntemplate <int Order>\n"
        "void Cartesian2Reduced(double *cm, double *rm);\n\n");

    fprintf(fp,"template <>\nvoid Cartesian2Reduced<%d>(double *cm, double *rm) { \n", order);
    FORALL_REDUCED_MULTIPOLES_BOUND(a,b,c,order) {
        int rm = rmap(a,b,c);
        fprintf(fp," rm[%d] = Cartesian2Reduced%d_%d_%d_%d(cm); \n", rm, order,a,b,c);
    }
    fprintf(fp,"}\n");

    fclose(fp);
}

void CreateDispatchFunction(){
    FILE *fp;
    char fn[1024];
    sprintf(fn, "Cartesian2ReducedDispatch.cpp");
    fp = fopen(fn, "w");
    assert(fp != NULL);

    fprintf(fp,
        "#include <cstdlib>\n"
        "#include <cstdio>\n"
        "template <int Order>\n"
        "void Cartesian2Reduced(double *cm, double *rm);\n\n"
        "void DispatchCartesian2Reduced(int order, double *cartesian, double *reduced) {\n"
        "   switch(order){\n"
    );

    for(int i = 0; i <= MAXORDER; i++){
        fprintf(fp, "        case %d: Cartesian2Reduced<%d>(cartesian, reduced); break;\n", i, i);
    }

    fprintf(fp, "\n        default: fprintf(stderr, \"Error: unknown order in dispatch\\n\"); exit(1); break;\n"
        "    }\n}\n");
    assert(fclose(fp) == 0);
}

int main(void) {

    CreateDispatchFunction();

    for(int order=0;order<=MAXORDER;order++)
        CreateC2R(order); 

}
