#include "header.cpp"
#include "threevector.hh"

#include "basemultipoles.h"

inline int basemultipoles::lmap(int a, int b, int c) {
    return a*(MAXORDER+1)*(MAXORDER+1)+b*(MAXORDER+1)+c;
}

basemultipoles::basemultipoles(int order) : order(order) {
    assert(order <= MAXORDER);
    completemultipolelength = (((order+1)*(order+2)*(order+3))/6);
    reducedmultipolelength = (order+1)*(order+1);
    cml = completemultipolelength;
    rml =  reducedmultipolelength;

    int m = 0;
    for(int a=0;a<=order;a++)
        for(int b=0;b<=order-a;b++)
            for(int c=0;c<=order-a-b;c++)  {
                _cmap[ lmap(a,b,c) ] = m;
                _invcmap[m].x = a;
                _invcmap[m].y = b;
                _invcmap[m].z = c;
                m++;
            }
    assert(m==completemultipolelength);

    m = 0;
    for(int c=0;c<=1;c++) 
        for(int a=0;a<=order-c;a++)
            for(int b=0;b<=order-c-a;b++) {
                _rmap[ lmap(a,b,c) ] = m;
                _invrmap[m].x = a;
                _invrmap[m].y = b;
                _invrmap[m].z = c;
                m++;
            }
    assert(m==reducedmultipolelength);

}

int basemultipoles::cmap(int a, int b, int c) {
#ifdef DEBUG
    assert(a>=0);
    assert(b>=0);
    assert(c>=0);
    assert(a<=order); 
    assert(b<=order);
    assert(c<=order);
    assert( a + b + c <= order);
#endif
    
    return _cmap[ lmap(a,b,c) ];
}

integer3 basemultipoles::invcmap(int m) {
#ifdef DEBUG
    assert(m >= 0);
    assert(m < completemultipolelength);
#endif

    assert(_invcmap[m].x >= 0 && _invcmap[m].x <= order );
    assert(_invcmap[m].y >= 0 && _invcmap[m].y <= order );
    assert(_invcmap[m].z >= 0 && _invcmap[m].z <= order );

    return _invcmap[m];
}

int basemultipoles::rmap(int a, int b, int c) {
#ifdef DEBUG
    assert(a>=0);
    assert(b>=0);
    assert(c>=0);
    assert(a<=order);
    assert(b<=order);
    assert(c<=1);
    assert( a + b + c <= order);
#endif

    return _rmap[ lmap(a,b,c) ];
}

integer3 basemultipoles::invrmap(int m) {
#ifdef DEBUG
    assert(m >= 0);
    assert(m < reducedmultipolelength);
#endif

    return _invrmap[m];
}

