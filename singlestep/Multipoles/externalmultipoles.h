#ifndef __D4DECL__
#define __D4DECL__
typedef struct { double v[4]; } d4;
#endif

#ifndef EXTERNALMULTIPOLES
#define EXTERNALMULTIPOLES

 extern "C" void MultipoleKernel1(   d4 *ip1x, d4 *ip2x, d4 *ip1y, d4 *ip2y, d4 *ip1z, d4 *ip2z, d4 *cx, d4 *cy, d4 *cz, d4 *globalM, d4 *mass1, d4 *mass2); 
 extern "C" void MultipoleKernel2(   d4 *ip1x, d4 *ip2x, d4 *ip1y, d4 *ip2y, d4 *ip1z, d4 *ip2z, d4 *cx, d4 *cy, d4 *cz, d4 *globalM, d4 *mass1, d4 *mass2); 
 extern "C" void MultipoleKernel3(   d4 *ip1x, d4 *ip2x, d4 *ip1y, d4 *ip2y, d4 *ip1z, d4 *ip2z, d4 *cx, d4 *cy, d4 *cz, d4 *globalM, d4 *mass1, d4 *mass2); 
 extern "C" void MultipoleKernel4(   d4 *ip1x, d4 *ip2x, d4 *ip1y, d4 *ip2y, d4 *ip1z, d4 *ip2z, d4 *cx, d4 *cy, d4 *cz, d4 *globalM, d4 *mass1, d4 *mass2); 
 extern "C" void MultipoleKernel5(   d4 *ip1x, d4 *ip2x, d4 *ip1y, d4 *ip2y, d4 *ip1z, d4 *ip2z, d4 *cx, d4 *cy, d4 *cz, d4 *globalM, d4 *mass1, d4 *mass2); 
 extern "C" void MultipoleKernel6(   d4 *ip1x, d4 *ip2x, d4 *ip1y, d4 *ip2y, d4 *ip1z, d4 *ip2z, d4 *cx, d4 *cy, d4 *cz, d4 *globalM, d4 *mass1, d4 *mass2); 
 extern "C" void MultipoleKernel7(   d4 *ip1x, d4 *ip2x, d4 *ip1y, d4 *ip2y, d4 *ip1z, d4 *ip2z, d4 *cx, d4 *cy, d4 *cz, d4 *globalM, d4 *mass1, d4 *mass2); 
 extern "C" void MultipoleKernel8(   d4 *ip1x, d4 *ip2x, d4 *ip1y, d4 *ip2y, d4 *ip1z, d4 *ip2z, d4 *cx, d4 *cy, d4 *cz, d4 *globalM, d4 *mass1, d4 *mass2); 
 extern "C" void MultipoleKernel9(   d4 *ip1x, d4 *ip2x, d4 *ip1y, d4 *ip2y, d4 *ip1z, d4 *ip2z, d4 *cx, d4 *cy, d4 *cz, d4 *globalM, d4 *mass1, d4 *mass2); 
 extern "C" void MultipoleKernel10(   d4 *ip1x, d4 *ip2x, d4 *ip1y, d4 *ip2y, d4 *ip1z, d4 *ip2z, d4 *cx, d4 *cy, d4 *cz, d4 *globalM, d4 *mass1, d4 *mass2); 
 extern "C" void MultipoleKernel11(   d4 *ip1x, d4 *ip2x, d4 *ip1y, d4 *ip2y, d4 *ip1z, d4 *ip2z, d4 *cx, d4 *cy, d4 *cz, d4 *globalM, d4 *mass1, d4 *mass2); 
 extern "C" void MultipoleKernel12(   d4 *ip1x, d4 *ip2x, d4 *ip1y, d4 *ip2y, d4 *ip1z, d4 *ip2z, d4 *cx, d4 *cy, d4 *cz, d4 *globalM, d4 *mass1, d4 *mass2); 
 extern "C" void MultipoleKernel13(   d4 *ip1x, d4 *ip2x, d4 *ip1y, d4 *ip2y, d4 *ip1z, d4 *ip2z, d4 *cx, d4 *cy, d4 *cz, d4 *globalM, d4 *mass1, d4 *mass2); 
 extern "C" void MultipoleKernel14(   d4 *ip1x, d4 *ip2x, d4 *ip1y, d4 *ip2y, d4 *ip1z, d4 *ip2z, d4 *cx, d4 *cy, d4 *cz, d4 *globalM, d4 *mass1, d4 *mass2); 
 extern "C" void MultipoleKernel15(   d4 *ip1x, d4 *ip2x, d4 *ip1y, d4 *ip2y, d4 *ip1z, d4 *ip2z, d4 *cx, d4 *cy, d4 *cz, d4 *globalM, d4 *mass1, d4 *mass2); 
 extern "C" void MultipoleKernel16(   d4 *ip1x, d4 *ip2x, d4 *ip1y, d4 *ip2y, d4 *ip1z, d4 *ip2z, d4 *cx, d4 *cy, d4 *cz, d4 *globalM, d4 *mass1, d4 *mass2); 

#ifdef AVX512MULTIPOLES
 extern "C" void Multipole512Kernel1( AVX512_DOUBLES &px, AVX512_DOUBLES &py, AVX512_DOUBLES &pz, AVX512_DOUBLES &cx, AVX512_DOUBLES &cy, AVX512_DOUBLES &cz, AVX512_DOUBLES *CM );
 extern "C" void Multipole512Kernel2( AVX512_DOUBLES &px, AVX512_DOUBLES &py, AVX512_DOUBLES &pz, AVX512_DOUBLES &cx, AVX512_DOUBLES &cy, AVX512_DOUBLES &cz, AVX512_DOUBLES *CM );
 extern "C" void Multipole512Kernel3( AVX512_DOUBLES &px, AVX512_DOUBLES &py, AVX512_DOUBLES &pz, AVX512_DOUBLES &cx, AVX512_DOUBLES &cy, AVX512_DOUBLES &cz, AVX512_DOUBLES *CM );
 extern "C" void Multipole512Kernel4( AVX512_DOUBLES &px, AVX512_DOUBLES &py, AVX512_DOUBLES &pz, AVX512_DOUBLES &cx, AVX512_DOUBLES &cy, AVX512_DOUBLES &cz, AVX512_DOUBLES *CM );
 extern "C" void Multipole512Kernel5( AVX512_DOUBLES &px, AVX512_DOUBLES &py, AVX512_DOUBLES &pz, AVX512_DOUBLES &cx, AVX512_DOUBLES &cy, AVX512_DOUBLES &cz, AVX512_DOUBLES *CM );
 extern "C" void Multipole512Kernel6( AVX512_DOUBLES &px, AVX512_DOUBLES &py, AVX512_DOUBLES &pz, AVX512_DOUBLES &cx, AVX512_DOUBLES &cy, AVX512_DOUBLES &cz, AVX512_DOUBLES *CM );
 extern "C" void Multipole512Kernel7( AVX512_DOUBLES &px, AVX512_DOUBLES &py, AVX512_DOUBLES &pz, AVX512_DOUBLES &cx, AVX512_DOUBLES &cy, AVX512_DOUBLES &cz, AVX512_DOUBLES *CM );
 extern "C" void Multipole512Kernel8( AVX512_DOUBLES &px, AVX512_DOUBLES &py, AVX512_DOUBLES &pz, AVX512_DOUBLES &cx, AVX512_DOUBLES &cy, AVX512_DOUBLES &cz, AVX512_DOUBLES *CM );
 extern "C" void Multipole512Kernel9( AVX512_DOUBLES &px, AVX512_DOUBLES &py, AVX512_DOUBLES &pz, AVX512_DOUBLES &cx, AVX512_DOUBLES &cy, AVX512_DOUBLES &cz, AVX512_DOUBLES *CM );
 extern "C" void Multipole512Kernel10( AVX512_DOUBLES &px, AVX512_DOUBLES &py, AVX512_DOUBLES &pz, AVX512_DOUBLES &cx, AVX512_DOUBLES &cy, AVX512_DOUBLES &cz, AVX512_DOUBLES *CM );
 extern "C" void Multipole512Kernel11( AVX512_DOUBLES &px, AVX512_DOUBLES &py, AVX512_DOUBLES &pz, AVX512_DOUBLES &cx, AVX512_DOUBLES &cy, AVX512_DOUBLES &cz, AVX512_DOUBLES *CM );
 extern "C" void Multipole512Kernel12( AVX512_DOUBLES &px, AVX512_DOUBLES &py, AVX512_DOUBLES &pz, AVX512_DOUBLES &cx, AVX512_DOUBLES &cy, AVX512_DOUBLES &cz, AVX512_DOUBLES *CM );
 extern "C" void Multipole512Kernel13( AVX512_DOUBLES &px, AVX512_DOUBLES &py, AVX512_DOUBLES &pz, AVX512_DOUBLES &cx, AVX512_DOUBLES &cy, AVX512_DOUBLES &cz, AVX512_DOUBLES *CM );
 extern "C" void Multipole512Kernel14( AVX512_DOUBLES &px, AVX512_DOUBLES &py, AVX512_DOUBLES &pz, AVX512_DOUBLES &cx, AVX512_DOUBLES &cy, AVX512_DOUBLES &cz, AVX512_DOUBLES *CM );
 extern "C" void Multipole512Kernel15( AVX512_DOUBLES &px, AVX512_DOUBLES &py, AVX512_DOUBLES &pz, AVX512_DOUBLES &cx, AVX512_DOUBLES &cy, AVX512_DOUBLES &cz, AVX512_DOUBLES *CM );
 extern "C" void Multipole512Kernel16( AVX512_DOUBLES &px, AVX512_DOUBLES &py, AVX512_DOUBLES &pz, AVX512_DOUBLES &cx, AVX512_DOUBLES &cy, AVX512_DOUBLES &cz, AVX512_DOUBLES *CM );
#endif

#endif
