#ifndef BASEMULTIPOLES
#define BASEMULTIPOLES


#define FOR(a,b,c) for(a=b;a<=c;a++) 

#define FORALL_COMPLETE_MULTIPOLES_BOUND(a,b,c,u) \
            FOR(a,0,u) FOR(b,0,u-a) FOR(c,0,u-a-b)

#define FORALL_REDUCED_MULTIPOLES_BOUND(a,b,c,u) \
            FOR(c,0,1) FOR(a,0,u-c) FOR(b,0,u-c-a) 

#define TRACEFREERECURSION(a,b,c,u) \
            FOR(c,2,u) FOR(a,0,u-c) FOR(b,0,u-c-a)

#define MAXORDER 16

class basemultipoles { 
public:
    basemultipoles(int order);
    virtual ~basemultipoles(void) { } 

    int cmap(int a, int b, int c);
    int rmap(int a, int b, int c);

    int _cmap[(MAXORDER+1)*(MAXORDER+1)*(MAXORDER+1)];
    integer3 _invcmap[(MAXORDER+1)*(MAXORDER+1)*(MAXORDER+1)];

    int _rmap[(MAXORDER+1)*(MAXORDER+1)*(MAXORDER+1)];
    integer3 _invrmap[(MAXORDER+1)*(MAXORDER+1)*(MAXORDER+1)];

    integer3 invcmap(int m);
    integer3 invrmap(int m);

    int reducedmultipolelength;
    int completemultipolelength;
    int rml, cml;
    int order;

    inline int lmap(int a, int b, int c);
}; 

#endif
