#include "config.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include <fmt/format.h>

#define FOR(a,b,c) for(a=b;a<=c;a++) 

int main(void) {
    FILE *fp;
    
    fp = fopen("CMAVX.c","w"); assert(fp!=NULL);

#ifdef AVXMULTIPOLES
    int i,j,k,ORDER,maxm,m;

    fmt::print(fp,"#include \"avxsseabrev.h\"  \n"
                "#include <stdio.h>\n"
                "#include <stdlib.h>\n");

    fmt::print(fp," typedef struct {{ double v[4]; }} d4; \n");

    fmt::print(fp,"#define X      \"%ymm0\" \n");
    fmt::print(fp,"#define Y      \"%ymm1\" \n");
    fmt::print(fp,"#define Z      \"%ymm2\" \n");

    fmt::print(fp,"#define XI     \"%ymm3\" \n");
    fmt::print(fp,"#define XIJ    \"%ymm4\" \n");
    fmt::print(fp,"#define XIJK   \"%ymm5\" \n");

    fmt::print(fp,"#define CX     \"%ymm4\" \n");
    fmt::print(fp,"#define CY     \"%ymm5\" \n");
    fmt::print(fp,"#define CZ     \"%ymm12\" \n");

    fmt::print(fp,"#define X2      \"%ymm7\" \n");
    fmt::print(fp,"#define Y2      \"%ymm8\" \n");
    fmt::print(fp,"#define Z2      \"%ymm9\" \n");

    fmt::print(fp,"#define XI2     \"%ymm11\" \n");
    fmt::print(fp,"#define XIJ2    \"%ymm12\" \n");
    fmt::print(fp,"#define XIJK2   \"%ymm13\"  \n");

    fmt::print(fp,"#define P      \"%ymm15\" \n");

    //fmt::print(fp,"extern \"C\"{{\n");

    
    for(ORDER=1;ORDER<=MAXORDER;ORDER++) {
        maxm = ((ORDER+1)*(ORDER+2)*(ORDER+3))/6;
        m = 0;

        fmt::print(fp, "void MultipoleKernel{:d}(d4 *ip1x, d4 *ip2x, d4 *ip1y, d4 *ip2y, d4 *ip1z, d4 *ip2z, d4 *cx, d4 *cy, d4 *cz, d4 *globalM, d4 *masses1, d4 *masses2) {{\n", ORDER); 

        fmt::print(fp, " VLOADPD(*cx->v, CX); VLOADPD(*cy->v, CY); VLOADPD(*cz->v, CZ); \n");
        fmt::print(fp, " VLOADPD(*masses1->v,XI); VLOADPD(*masses2->v,XI2); \n");

        fmt::print(fp, " VLOADPD(*ip1x->v, X); VLOADPD(*ip1y->v, Y); VLOADPD(*ip1z->v, Z); \n");
        fmt::print(fp, " VLOADPD(*ip2x->v, X2); VLOADPD(*ip2y->v, Y2); VLOADPD(*ip2z->v, Z2);  \n");

        fmt::print(fp, " VSUBPD(CX,X,X); VSUBPD(CX,X2,X2); \n");
        fmt::print(fp, " VSUBPD(CY,Y,Y); VSUBPD(CY,Y2,Y2); \n");
        fmt::print(fp, " VSUBPD(CZ,Z,Z); VSUBPD(CZ,Z2,Z2); \n");

        fmt::print(fp, "VLOADPD(*globalM->v,P);\n");

        FOR(i,0,ORDER) {
            fmt::print(fp, "VMOVAPD(XI,XIJ);\n");
            fmt::print(fp, "VMOVAPD(XI2,XIJ2);\n");
        
            FOR(j,0,ORDER-i) {
                fmt::print(fp, "VMOVAPD(XIJ,XIJK);\n");
                fmt::print(fp, "VMOVAPD(XIJ2,XIJK2);\n");

                FOR(k,0,ORDER-i-j) { 
    
                    fmt::print(fp, "VADDPD(XIJK,P,P);\n");
                    fmt::print(fp, "VADDPD(XIJK2,P,P);\n");

                    fmt::print(fp, "VSTORPD(P,*globalM->v);\n");
                    m++; 
                    if(m!=maxm) {
                        fmt::print(fp, "globalM++;\n");
                        fmt::print(fp, "VLOADPD(*globalM->v,P);\n");
                    }

                    if(k<ORDER-i-j) fmt::print(fp, "VMULPD(Z,XIJK,XIJK);\n");
                    if(k<ORDER-i-j) fmt::print(fp, "VMULPD(Z2,XIJK2,XIJK2);\n");
                }
                if(j<ORDER-i) fmt::print(fp, "VMULPD(Y,XIJ,XIJ);\n");
                if(j<ORDER-i) fmt::print(fp, "VMULPD(Y2,XIJ2,XIJ2);\n");
            }
            if(i<ORDER) fmt::print(fp, "VMULPD(X,XI,XI);\n");
            if(i<ORDER) fmt::print(fp, "VMULPD(X2,XI2,XI2);\n");
        }
        fmt::print(fp, "}}\n");
    }

    //fmt::print(fp,"}}\n");

    fmt::print(fp,"#undef X       \n");
    fmt::print(fp,"#undef Y       \n");
    fmt::print(fp,"#undef Z       \n");
    fmt::print(fp,"#undef XI      \n");
    fmt::print(fp,"#undef XIJ     \n");
    fmt::print(fp,"#undef XIJK    \n");

    fmt::print(fp,"#undef X2       \n");
    fmt::print(fp,"#undef Y2       \n");
    fmt::print(fp,"#undef Z2       \n");
    fmt::print(fp,"#undef XI2      \n");
    fmt::print(fp,"#undef XIJ2     \n");
    fmt::print(fp,"#undef XIJK2    \n");

    fmt::print(fp,"#undef P        \n");

    fmt::print(fp,"\nvoid DispatchMultipoleAVXKernel(int order, d4 *ip1x, d4 *ip2x, d4 *ip1y, d4 *ip2y, d4 *ip1z, d4 *ip2z, d4 *cx, d4 *cy, d4 *cz, d4 *globalM, d4 *masses1, d4 *masses2){{\n"
                "   switch(order) {{\n");

    for(i = 1; i <= MAXORDER; i++){
        fmt::print(fp, "       case {:d}: MultipoleKernel{:d}(ip1x, ip2x, ip1y, ip2y, ip1z, ip2z, cx, cy, cz, globalM, masses1, masses2); break;\n", i, i);
    }
    fmt::print(fp, "       default: fprintf(stderr, \"Unknown order in multipole dispatch\\n\"); exit(1); break;\n"
                "   }}\n}}");

#endif
    fclose(fp);

    return 0;
}
