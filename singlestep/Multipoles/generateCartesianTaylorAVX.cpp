// Copyright 2012-2025 The Abacus Developers
// SPDX-License-Identifier: GPL-3.0-or-later

#include "config.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include <fmt/format.h>

int ORDER;

#define FOR(a,b,c) for(a=b;a<=c;a++) 

int main(void) {
    FILE *fp;

    fp = fopen("ETAVX.c","w");
    assert(fp!=NULL);

#ifdef AVXMULTIPOLES
    int i,j,k;

    fmt::print(fp,"#include \"avxsseabrev.h\" \n"
                "#include <stdio.h>\n"
                "#include <stdlib.h>\n");
    fmt::print(fp,"#include \"assert.h\" \n");

    fmt::print(fp,"#define AX      \"%ymm0\" \n");
    fmt::print(fp,"#define AY      \"%ymm1\" \n");
    fmt::print(fp,"#define AZ      \"%ymm2\" \n");
    fmt::print(fp,"#define FIJK    \"%ymm3\" \n");
    fmt::print(fp,"#define FIJ     \"%ymm4\" \n");
    fmt::print(fp,"#define FI      \"%ymm5\" \n");
    fmt::print(fp,"#define X       \"%ymm6\" \n");
    fmt::print(fp,"#define Y       \"%ymm7\" \n");
    fmt::print(fp,"#define Z       \"%ymm8\" \n");
    fmt::print(fp,"#define TX      \"%ymm9\" \n");
    fmt::print(fp,"#define TY      \"%ymm10\" \n");
    fmt::print(fp,"#define TZ      \"%ymm11\" \n");

    fmt::print(fp," typedef struct {{ double v[4]; }} d4; \n");

    //fmt::print(fp,"extern \"C\"{{\n");


    for(ORDER=1;ORDER<=MAXORDER;ORDER++) {

        fmt::print(fp,"void TaylorKernel{:d}( d4 *px, d4 *py, d4 *pz,  \n", ORDER);
        fmt::print(fp,"                    d4 *cx, d4 *cy, d4 *cz,  \n");
        fmt::print(fp,"                    d4 *tx, d4 *ty, d4 *tz,  \n");
        fmt::print(fp,"                    d4 *ax, d4 *ay, d4 *az  )  {{ \n"); 

        // assumption ax = 1,1,1,1 to start off the FI 

          fmt::print(fp,"\tVLOADPD(*ax->v, FI); \n");
        
          fmt::print(fp,"\tVLOADPD(*cx->v, AX); \n");
          fmt::print(fp,"\tVLOADPD(*cy->v, AY); \n");
          fmt::print(fp,"\tVLOADPD(*cz->v, AZ); \n");
        
         fmt::print(fp,"\tVLOADPD(*px->v, X); \n");
         fmt::print(fp,"\tVLOADPD(*py->v, Y); \n");
         fmt::print(fp,"\tVLOADPD(*pz->v, Z); \n");
 
         fmt::print(fp,"\tVSUBPD(AX,X,X); \n");
         fmt::print(fp,"\tVSUBPD(AY,Y,Y); \n");
         fmt::print(fp,"\tVSUBPD(AZ,Z,Z); \n");

        fmt::print(fp,"\tVXORPS(AX,AX,AX);  \n");
        fmt::print(fp,"\tVXORPS(AY,AY,AY);  \n");
        fmt::print(fp,"\tVXORPS(AZ,AZ,AZ);  \n");
        
        FOR(i,0,ORDER-1) {
            fmt::print(fp,"\t\tVMOVAPD(FI,FIJ);  \n");
            
            FOR(j,0,ORDER-1-i) {
                fmt::print(fp,"\t\t\tVMOVAPD(FIJ,FIJK); \n");
                
                FOR(k,0,ORDER-1-i-j) {
                    fmt::print(fp,"\t\t\t\t VLOADPD(*tx->v,TX); \n\n");
                    fmt::print(fp,"\t\t\t\t tx++;            \n");
                    
                    fmt::print(fp,"\t\t\t\t VLOADPD(*ty->v,TY); \n\n");
                    fmt::print(fp,"\t\t\t\t ty++;            \n");
                    
                    fmt::print(fp,"\t\t\t\t VLOADPD(*tz->v,TZ); \n\n");
                    fmt::print(fp,"\t\t\t\t tz++;            \n");
                                        
                    fmt::print(fp,"\t\t\t\t VMULPD(FIJK, TX, TX ); \n");
                    fmt::print(fp,"\t\t\t\t VADDPD(TX, AX, AX); \n");
                    
                    fmt::print(fp,"\t\t\t\t VMULPD(FIJK, TY, TY ); \n");
                    fmt::print(fp,"\t\t\t\t VADDPD(TY, AY, AY); \n");
                    
                    fmt::print(fp,"\t\t\t\t VMULPD(FIJK, TZ, TZ); \n");
                    fmt::print(fp,"\t\t\t\t VADDPD(TZ, AZ, AZ); \n");
                    
                    if(k<ORDER-1-i-j) fmt::print(fp,"\t\t\t\t VMULPD(Z,FIJK, FIJK); \n");
                }
                if(j<ORDER-1-i) fmt::print(fp,"\t\t\t VMULPD(Y,FIJ,FIJ); \n");
            }
            if(i<ORDER-1) fmt::print(fp,"\t\t VMULPD(X,FI,FI); \n");
        }
        fmt::print(fp,"\t VSTORPD(AX, *ax->v); VSTORPD(AY, *ay->v); VSTORPD(AZ, *az->v); \n");
        fmt::print(fp,"}}\n");
    }
    //fmt::print(fp,"}}\n");
    
    fmt::print(fp,"#undef AX      \n");
    fmt::print(fp,"#undef AY      \n");
    fmt::print(fp,"#undef AZ      \n");
    fmt::print(fp,"#undef FIJK    \n");
    fmt::print(fp,"#undef FIJ     \n");
    fmt::print(fp,"#undef FI      \n");
    fmt::print(fp,"#undef X       \n");
    fmt::print(fp,"#undef Y       \n");
    fmt::print(fp,"#undef Z       \n");
    fmt::print(fp,"#undef Q0      \n");
    fmt::print(fp,"#undef Q1      \n");
    fmt::print(fp,"#undef Q2      \n");

    fmt::print(fp,"\nvoid DispatchTaylorAVXKernel(int order, d4 *px, d4 *py, d4 *pz, \n"
                "                    d4 *cx, d4 *cy, d4 *cz, \n"
                "                    d4 *tx, d4 *ty, d4 *tz, \n"
                "                    d4 *ax, d4 *ay, d4 *az){{\n"
                "   switch(order) {{\n");

    for(i = 1; i <= MAXORDER; i++){
        fmt::print(fp, "       case {:d}: TaylorKernel{:d}(px, py, pz, cx, cy, cz, tx, ty, tz, ax, ay, az); break;\n", i, i);
    }
    fmt::print(fp, "       default: fprintf(stderr, \"Unknown order in Taylor dispatch\\n\"); exit(1); break;\n"
                "   }}\n}}");

#endif

    fclose(fp);

    return 0;
}
