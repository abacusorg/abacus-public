#ifndef EXTERNALTAYLOR
#define EXTERNALTAYLOR

#ifndef __D4DECL__
#define __D4DECL__
typedef struct { double v[4]; } d4; 
#endif

extern "C" void TaylorKernel0( d4 * px, d4 * py, d4 * pz, d4 * cx, d4 * cy, d4 * cz, d4 * tx, d4 * ty, d4 * tz, d4 * ax, d4 * ay, d4 * az ); 
extern "C" void TaylorKernel1( d4 * px, d4 * py, d4 * pz, d4 * cx, d4 * cy, d4 * cz, d4 * tx, d4 * ty, d4 * tz, d4 * ax, d4 * ay, d4 * az ); 
extern "C" void TaylorKernel2( d4 * px, d4 * py, d4 * pz, d4 * cx, d4 * cy, d4 * cz, d4 * tx, d4 * ty, d4 * tz, d4 * ax, d4 * ay, d4 * az ); 
extern "C" void TaylorKernel3( d4 * px, d4 * py, d4 * pz, d4 * cx, d4 * cy, d4 * cz, d4 * tx, d4 * ty, d4 * tz, d4 * ax, d4 * ay, d4 * az ); 
extern "C" void TaylorKernel4( d4 * px, d4 * py, d4 * pz, d4 * cx, d4 * cy, d4 * cz, d4 * tx, d4 * ty, d4 * tz, d4 * ax, d4 * ay, d4 * az ); 
extern "C" void TaylorKernel5( d4 * px, d4 * py, d4 * pz, d4 * cx, d4 * cy, d4 * cz, d4 * tx, d4 * ty, d4 * tz, d4 * ax, d4 * ay, d4 * az ); 
extern "C" void TaylorKernel6( d4 * px, d4 * py, d4 * pz, d4 * cx, d4 * cy, d4 * cz, d4 * tx, d4 * ty, d4 * tz, d4 * ax, d4 * ay, d4 * az ); 
extern "C" void TaylorKernel7( d4 * px, d4 * py, d4 * pz, d4 * cx, d4 * cy, d4 * cz, d4 * tx, d4 * ty, d4 * tz, d4 * ax, d4 * ay, d4 * az ); 
extern "C" void TaylorKernel8( d4 * px, d4 * py, d4 * pz, d4 * cx, d4 * cy, d4 * cz, d4 * tx, d4 * ty, d4 * tz, d4 * ax, d4 * ay, d4 * az ); 
extern "C" void TaylorKernel9( d4 * px, d4 * py, d4 * pz, d4 * cx, d4 * cy, d4 * cz, d4 * tx, d4 * ty, d4 * tz, d4 * ax, d4 * ay, d4 * az ); 
extern "C" void TaylorKernel10( d4 * px, d4 * py, d4 * pz, d4 * cx, d4 * cy, d4 * cz, d4 * tx, d4 * ty, d4 * tz, d4 * ax, d4 * ay, d4 * az ); 
extern "C" void TaylorKernel11( d4 * px, d4 * py, d4 * pz, d4 * cx, d4 * cy, d4 * cz, d4 * tx, d4 * ty, d4 * tz, d4 * ax, d4 * ay, d4 * az ); 
extern "C" void TaylorKernel12( d4 * px, d4 * py, d4 * pz, d4 * cx, d4 * cy, d4 * cz, d4 * tx, d4 * ty, d4 * tz, d4 * ax, d4 * ay, d4 * az ); 
extern "C" void TaylorKernel13( d4 * px, d4 * py, d4 * pz, d4 * cx, d4 * cy, d4 * cz, d4 * tx, d4 * ty, d4 * tz, d4 * ax, d4 * ay, d4 * az ); 
extern "C" void TaylorKernel14( d4 * px, d4 * py, d4 * pz, d4 * cx, d4 * cy, d4 * cz, d4 * tx, d4 * ty, d4 * tz, d4 * ax, d4 * ay, d4 * az ); 
extern "C" void TaylorKernel15( d4 * px, d4 * py, d4 * pz, d4 * cx, d4 * cy, d4 * cz, d4 * tx, d4 * ty, d4 * tz, d4 * ax, d4 * ay, d4 * az ); 
extern "C" void TaylorKernel16( d4 * px, d4 * py, d4 * pz, d4 * cx, d4 * cy, d4 * cz, d4 * tx, d4 * ty, d4 * tz, d4 * ax, d4 * ay, d4 * az ); 

#ifdef AVX512MULTIPOLES
#include "avx512_calls.h"

extern "C" void Taylor512Kernel0(AVX512_DOUBLES &px, AVX512_DOUBLES &py, AVX512_DOUBLES &pz, AVX512_DOUBLES &cx, AVX512_DOUBLES &cy, AVX512_DOUBLES &cz,
                                 AVX512_DOUBLES *tx, AVX512_DOUBLES *ty, AVX512_DOUBLES *tz, AVX512_DOUBLES &ax, AVX512_DOUBLES &ay, AVX512_DOUBLES &az ); 
extern "C" void Taylor512Kernel1(AVX512_DOUBLES &px, AVX512_DOUBLES &py, AVX512_DOUBLES &pz, AVX512_DOUBLES &cx, AVX512_DOUBLES &cy, AVX512_DOUBLES &cz,
                                 AVX512_DOUBLES *tx, AVX512_DOUBLES *ty, AVX512_DOUBLES *tz, AVX512_DOUBLES &ax, AVX512_DOUBLES &ay, AVX512_DOUBLES &az );
extern "C" void Taylor512Kernel2(AVX512_DOUBLES &px, AVX512_DOUBLES &py, AVX512_DOUBLES &pz, AVX512_DOUBLES &cx, AVX512_DOUBLES &cy, AVX512_DOUBLES &cz,
                                 AVX512_DOUBLES *tx, AVX512_DOUBLES *ty, AVX512_DOUBLES *tz, AVX512_DOUBLES &ax, AVX512_DOUBLES &ay, AVX512_DOUBLES &az );
extern "C" void Taylor512Kernel3(AVX512_DOUBLES &px, AVX512_DOUBLES &py, AVX512_DOUBLES &pz, AVX512_DOUBLES &cx, AVX512_DOUBLES &cy, AVX512_DOUBLES &cz,
                                 AVX512_DOUBLES *tx, AVX512_DOUBLES *ty, AVX512_DOUBLES *tz, AVX512_DOUBLES &ax, AVX512_DOUBLES &ay, AVX512_DOUBLES &az );
extern "C" void Taylor512Kernel4(AVX512_DOUBLES &px, AVX512_DOUBLES &py, AVX512_DOUBLES &pz, AVX512_DOUBLES &cx, AVX512_DOUBLES &cy, AVX512_DOUBLES &cz,
                                 AVX512_DOUBLES *tx, AVX512_DOUBLES *ty, AVX512_DOUBLES *tz, AVX512_DOUBLES &ax, AVX512_DOUBLES &ay, AVX512_DOUBLES &az );
extern "C" void Taylor512Kernel5(AVX512_DOUBLES &px, AVX512_DOUBLES &py, AVX512_DOUBLES &pz, AVX512_DOUBLES &cx, AVX512_DOUBLES &cy, AVX512_DOUBLES &cz,
                                 AVX512_DOUBLES *tx, AVX512_DOUBLES *ty, AVX512_DOUBLES *tz, AVX512_DOUBLES &ax, AVX512_DOUBLES &ay, AVX512_DOUBLES &az );
extern "C" void Taylor512Kernel6(AVX512_DOUBLES &px, AVX512_DOUBLES &py, AVX512_DOUBLES &pz, AVX512_DOUBLES &cx, AVX512_DOUBLES &cy, AVX512_DOUBLES &cz,
                                 AVX512_DOUBLES *tx, AVX512_DOUBLES *ty, AVX512_DOUBLES *tz, AVX512_DOUBLES &ax, AVX512_DOUBLES &ay, AVX512_DOUBLES &az );
extern "C" void Taylor512Kernel7(AVX512_DOUBLES &px, AVX512_DOUBLES &py, AVX512_DOUBLES &pz, AVX512_DOUBLES &cx, AVX512_DOUBLES &cy, AVX512_DOUBLES &cz,
                                 AVX512_DOUBLES *tx, AVX512_DOUBLES *ty, AVX512_DOUBLES *tz, AVX512_DOUBLES &ax, AVX512_DOUBLES &ay, AVX512_DOUBLES &az );
extern "C" void Taylor512Kernel8(AVX512_DOUBLES &px, AVX512_DOUBLES &py, AVX512_DOUBLES &pz, AVX512_DOUBLES &cx, AVX512_DOUBLES &cy, AVX512_DOUBLES &cz,
                                 AVX512_DOUBLES *tx, AVX512_DOUBLES *ty, AVX512_DOUBLES *tz, AVX512_DOUBLES &ax, AVX512_DOUBLES &ay, AVX512_DOUBLES &az );
extern "C" void Taylor512Kernel9(AVX512_DOUBLES &px, AVX512_DOUBLES &py, AVX512_DOUBLES &pz, AVX512_DOUBLES &cx, AVX512_DOUBLES &cy, AVX512_DOUBLES &cz,
                                 AVX512_DOUBLES *tx, AVX512_DOUBLES *ty, AVX512_DOUBLES *tz, AVX512_DOUBLES &ax, AVX512_DOUBLES &ay, AVX512_DOUBLES &az );
extern "C" void Taylor512Kernel10(AVX512_DOUBLES &px, AVX512_DOUBLES &py, AVX512_DOUBLES &pz, AVX512_DOUBLES &cx, AVX512_DOUBLES &cy, AVX512_DOUBLES &cz,
                                 AVX512_DOUBLES *tx, AVX512_DOUBLES *ty, AVX512_DOUBLES *tz, AVX512_DOUBLES &ax, AVX512_DOUBLES &ay, AVX512_DOUBLES &az );
extern "C" void Taylor512Kernel11(AVX512_DOUBLES &px, AVX512_DOUBLES &py, AVX512_DOUBLES &pz, AVX512_DOUBLES &cx, AVX512_DOUBLES &cy, AVX512_DOUBLES &cz,
                                 AVX512_DOUBLES *tx, AVX512_DOUBLES *ty, AVX512_DOUBLES *tz, AVX512_DOUBLES &ax, AVX512_DOUBLES &ay, AVX512_DOUBLES &az );
extern "C" void Taylor512Kernel12(AVX512_DOUBLES &px, AVX512_DOUBLES &py, AVX512_DOUBLES &pz, AVX512_DOUBLES &cx, AVX512_DOUBLES &cy, AVX512_DOUBLES &cz,
                                 AVX512_DOUBLES *tx, AVX512_DOUBLES *ty, AVX512_DOUBLES *tz, AVX512_DOUBLES &ax, AVX512_DOUBLES &ay, AVX512_DOUBLES &az );
extern "C" void Taylor512Kernel13(AVX512_DOUBLES &px, AVX512_DOUBLES &py, AVX512_DOUBLES &pz, AVX512_DOUBLES &cx, AVX512_DOUBLES &cy, AVX512_DOUBLES &cz,
                                 AVX512_DOUBLES *tx, AVX512_DOUBLES *ty, AVX512_DOUBLES *tz, AVX512_DOUBLES &ax, AVX512_DOUBLES &ay, AVX512_DOUBLES &az );
extern "C" void Taylor512Kernel14(AVX512_DOUBLES &px, AVX512_DOUBLES &py, AVX512_DOUBLES &pz, AVX512_DOUBLES &cx, AVX512_DOUBLES &cy, AVX512_DOUBLES &cz,
                                 AVX512_DOUBLES *tx, AVX512_DOUBLES *ty, AVX512_DOUBLES *tz, AVX512_DOUBLES &ax, AVX512_DOUBLES &ay, AVX512_DOUBLES &az );
extern "C" void Taylor512Kernel15(AVX512_DOUBLES &px, AVX512_DOUBLES &py, AVX512_DOUBLES &pz, AVX512_DOUBLES &cx, AVX512_DOUBLES &cy, AVX512_DOUBLES &cz,
                                 AVX512_DOUBLES *tx, AVX512_DOUBLES *ty, AVX512_DOUBLES *tz, AVX512_DOUBLES &ax, AVX512_DOUBLES &ay, AVX512_DOUBLES &az );
extern "C" void Taylor512Kernel16(AVX512_DOUBLES &px, AVX512_DOUBLES &py, AVX512_DOUBLES &pz, AVX512_DOUBLES &cx, AVX512_DOUBLES &cy, AVX512_DOUBLES &cz,
                                 AVX512_DOUBLES *tx, AVX512_DOUBLES *ty, AVX512_DOUBLES *tz, AVX512_DOUBLES &ax, AVX512_DOUBLES &ay, AVX512_DOUBLES &az );
#endif


#endif
