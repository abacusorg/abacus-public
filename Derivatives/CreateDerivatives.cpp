#include <stdlib.h>
#include <fstream>

// Set precision of calculation. If QUAD_DOUBLE = 1, we will use
// quad double precision (quad_double.cpp). If QUAD_DOUBLE = 0,
// we will use double precision. 
#ifndef QUAD_DOUBLE
#define QUAD_DOUBLE 0
#endif

#if QUAD_DOUBLE
        #include "quad_double.cpp"
#else  
        typedef double qd_real;
        #define to_double(x) x;
        #define to_qd_real(x) x;
        #define qd_real(x) std::stod(x);
#endif
        

#include "header.cpp"
#include "threevector.hh"
#include "fftw3.h"

typedef ThreeVector<qd_real> qd_real3;

#include "../Multipoles/basemultipoles.cpp"
#include "derivatives.cpp"
#include "order32derivatives.cpp"
#include "../include/STimer.cc"


#define ULLI unsigned long long int 

ULLI CPD;
ULLI CPD3;
ULLI CPDHALF;

ULLI CompressedMultipoleLengthXY;
ULLI CompressedMultipoleLengthXYZ;

ULLI linearindex(int x, int y, int z) { return x*CPD*CPD + y*CPD + z; }

ULLI linearFFTW(int i, int j, int k) {
    int indx_i = i<0?CPD+i:i;
    int indx_j = j<0?CPD+j:j;
    int indx_k = k<0?CPD+k:k;

    ULLI l = indx_i*CPD*CPD + indx_j*CPD + indx_k;
    assert(l>=0);
    assert(l<CPD*CPD*CPD);
    return l;
}

#define RINDEXY(x,y)    (((x)*(x+1))/2 + (y))
#define RINDEXYZ(x,y,z) (((z)*(1+CPD)*(3+CPD))/8 + ((x)*((x)+1))/2 + y)
#define sign(i)         (i>0?1:-1) 

//STimer myTimer; //if timing performance, uncomment. 

//calculate the derivatives tensor for every cell vector, \mathcal{D}^{ABC}_{jkl}. (Eqn. 4.1)
void FormDerivatives(int inner_radius, int order, int far_radius, int slabnumber) {
                
    int nthread = omp_get_max_threads();
    fmt::print("Initializing CreateDerivatives with {:d} OMP threads.\n", nthread);        
    double *rdfar[nthread];
    double *rdnear[nthread];
    double *FD[nthread];

    for(int i = 0; i < nthread; i++){
        rdfar[i] = new double[1024];
        rdnear[i] = new double[1024];
        FD[i] = new double[1024];
    }

    Order32Derivatives OD(far_radius);

    Derivatives RD(order);

    int rml = (order+1)*(order+1);

    ULLI fdsize = sizeof(double) * rml*CompressedMultipoleLengthXY;
    fmt::print("Allocating {:d} GB for slab derivs\n", (int) (fdsize/(1<<30)) );
    double *FDSlab = (double *) malloc(fdsize);
    assert(FDSlab != NULL);

        //We are now going to create the far field derivatives tensor for our specified number of cells per dimension (CPD, or K in the write-up) and order p. 
        //See Section 4: 'Pre-computing the far field derivatives tensor' in the Abacus Force Solver Documentation to follow along. The notation used in the comments corresponds to the notation in the documentation. Code comments in parantheses relate the code parameters to the write-up's notation.  
                
        //For every cell vector \mathbf{c}_{jkl}, compute \mathcal{D}^{ABC}_{jkl} for a given K (CPD), L_{inner} (inner radius), L_{outer} (outer radius), and p (order). 
    for(int k=0;k<=CPDHALF;k++) { //loop over z direction within the home box. Each z-slice through the simulation box will be called a slab and we will store the computed derivatives tensor for each slab in slab_deriv_CPD_order_innerRadius_farRadius_0__k files. 
        if (slabnumber>=0 && slabnumber!=k) continue;
            // If we were asked to do only one slab, skip each iteration of this loop until you reach the desired slab within the home box. 

        // Check if the desired file already exists.  If so, skip the computation in this loop.
        fs::path fn = fmt::format("slabderiv_{:d}_{:d}_{:d}_{:d}_{:d}__{:d}",CPD,order,inner_radius,far_radius, 0, k);
        FILE *fp;
        fp = fopen(fn.c_str(),"rb");
        if(fp!=NULL) {
            // Yes, the file already exists!  Don't repeat the work, so skip this iteration of the loop..
            fclose(fp);
            continue;  
        }
        

        // If the desired file doesn't exist, make it. 
        #pragma omp parallel for schedule(dynamic,1)
        for(int i=0;i<=CPDHALF;i++) { //loop over x direction within the home box. 
            for(int j=0;j<=i;j++) { //loop over y direction within the home box -- we only need to do one triangle of one quadrant because of the symmetries of the derivatives tensor, as discussed below Equation 3.9. 
                int p = omp_get_thread_num();

                for(int m=0;m<rml;m++) FD[p][m] = 0;

                double3 r;
                r.x = i; r.x = r.x/CPD; //find the components of the cell vector \mathbf{c}_{jkl} (in code notation, \mathbf{c}_{ijk} = { i/CPD, j/CPD, k/CPD} ). 
                r.y = j; r.y = r.y/CPD;
                r.z = k; r.z = r.z/CPD;
                                                                        
                if( (abs(i)>inner_radius) || (abs(j)>inner_radius) || (abs(k)>inner_radius) ) //check if the given cell vector \mathbf{c}_{jkl} falls outside of L_{inner}. If not, we will use equation 4.1a to calculate the derivatives tensor. If it does, we will use equation 4.1b: 
                                {
                    RD.ReducedDerivatives(r,&(FD[p][0])); //Calculate D^{ABC}(\mathbf{c}_{jkl}) explicitly using the trace-free and recursion relations discussed in Section 2.2 and coded in the file `derivatives.cpp'. 
                }
                                

                OD.Derivative( r, &(rdnear[p][0]), &(rdfar[p][0]),  order ); //Call order32derivatives.cpp to calculate the first two terms in Equation 4.1 -- 
                                //1)  Contribution to the derivatives tensor from the outer far field region (L_{outer} < n): \sum\limits_{ \mathcal{B}\left(\mathbf{n}, L_{\mathrm{outer}}, \infty\right)}D^{ABC}(\mathbf{n} +  \mathbf{c}_{jkl}) 
                                //2)  Contribution to the derivatives tensor from the inner far field region (L_{inner} < n <= L_{outer}): \sum\limits_{\mathcal{B}\left(\mathbf{n}, 1, L_{\mathrm{outer}}\right)} D^{ABC}(\mathbf{n} + \mathbf{c}_{jkl}) 
                                
                                //Sum together all contributions to the derivatives tensor (Eqn. 4.1) and store. 
                                for(int m=0;m<rml;m++) {
                    ULLI idx = m*CompressedMultipoleLengthXY + RINDEXY(i,j);
                    FDSlab[idx] = FD[p][m] + rdnear[p][m] + rdfar[p][m]; //combine all of the above
                }
            }
        }
                
        // Store FDSlab (containing derivatives tensor for the given slab) to disk. 
        fp = fopen(fn.c_str(),"wb");
        assert(fp!=NULL);
        fwrite(&(FDSlab[0]), sizeof(double), rml*CompressedMultipoleLengthXY, fp); 
        fclose(fp);
    }
    free(FDSlab);

    for(int i = 0; i < nthread; i++){
        delete[] rdfar[i];
        delete[] rdnear[i];
        delete[] FD[i];
    }

    return;
}

//calculate total far derivatives tensor for entire simulation box by merging components for every cell. {\mathcal{D}^{ABC}_{jkl}} = \mathcal{D}^{ABC}.
void MergeDerivatives(int inner_radius, int order, int far_radius, double *FarDerivatives) {
    int rml = (order+1)*(order+1);
    ULLI fdsize = sizeof(double) * rml*CompressedMultipoleLengthXY;
    fmt::print("Allocating {:d} GB for far derivs\n", (int) (fdsize/(1<<30)) );
    double *FDSlab = (double *) malloc(fdsize);
    assert(FDSlab != NULL);
        
        //Merge all slab's derivatives tensors into a single farderivatives file. This will be used to calculate \hat{\mathcal{D}}^{ABC}: the Fourier transform of the derivatives tensor required to calculate the Taylor coefficients for any given sink cell by convolving the multipole moments and the derivatives tensor. 

    for(int k=0;k<=CPDHALF;k++) {
                // Fetch FDSlab from disk
                fs::path fn = fmt::format("slabderiv_{:d}_{:d}_{:d}_{:d}_{:d}__{:d}",CPD,order,inner_radius,far_radius, 0, k);
                FILE *fp;
                fp = fopen(fn.c_str(),"rb");
                assert(fp!=NULL);
                ULLI sizeread = fread(&(FDSlab[0]), sizeof(double), rml*CompressedMultipoleLengthXY, fp); 
                assert(sizeread == rml*CompressedMultipoleLengthXY);
                fclose(fp);
                
                #pragma omp parallel for schedule(dynamic,1) 
                for(int m=0;m<rml;m++) {
                    for(int i=0;i<=CPDHALF;i++) {
                        for(int j=0;j<=i;j++) {
                            ULLI idxy = m*CompressedMultipoleLengthXY + RINDEXY(i,j);
                            ULLI idx = m*CompressedMultipoleLengthXYZ + RINDEXYZ(i,j,k);
                            FarDerivatives[idx] = FDSlab[idxy];
                        }
                    }
                }
                // FarDerivatives is in the order [m][k][i][j], but subject to the 
                // constraing that j<=i and that i,j,k are all in the first octant.
    }
    free(FDSlab);
    return;
}

//create empty fourier files to store the fourier transforms of the far derivatives in, to be calculated in Part2. 
void CreateFourierFiles(int order, int inner_radius, int far_radius) {
    for(int z=0;z<(CPD+1)/2;z++) {
        FILE *fp;
        fs::path fn = fmt::format("fourierspace_{:d}_{:d}_{:d}_{:d}_{:d}",CPD,order,inner_radius,far_radius,z);
        int cpd = CPD;
        fp = fopen(fn.c_str(),"w");
        assert(fp!=NULL);
        fclose(fp);
    }
}









// //calculate Fast Fourier transform of far derivatives tensor and store.
// void Part2(int order, int inner_radius, int far_radius) {
//
//
//         STimer part2timer;
//         part2timer.Start();
//
//     basemultipoles bm(order);
//     int cpd = CPD;
//
//     ULLI fdsize = sizeof(double) * CompressedMultipoleLengthXYZ;
//     fmt::print("Allocating {:d} GB\n", (int) (2*fdsize/(1<<30)) );
//     double *FarDerivatives_ab = (double *) malloc(fdsize);
//     double *FarDerivatives_ba = (double *) malloc(fdsize);
//     assert( FarDerivatives_ab != NULL );
//     assert( FarDerivatives_ba != NULL );
//
//     ULLI tdsize = sizeof(double)*CPD3;
//     fmt::print("Allocating {:d} GB\n", (int) (2*tdsize/(1<<30)) );
//     double *tdprime = (double *) malloc(tdsize);
//     double *td = (double *) malloc(tdsize);
//     assert(tdprime!=NULL);
//     assert(td!=NULL);
//
//     ULLI tmpDsize = sizeof(Complex)*CPD3;
//     fmt::print("Allocating {:d} GB\n", (int) (tmpDsize/(1<<30)) );
//     Complex *tmpD   = (Complex *) malloc(tmpDsize);
//     assert(tmpD!=NULL);
//
//     fftw_plan plan_forward_1d;
//      Complex in_1d[CPD];
//     Complex out_1d[CPD];
//
//     fftw_plan plan_forward_1d_r2c;
//     double   in_r2c[CPD];
//     Complex out_r2c[CPD];
//
//     //plan_forward_1d  =  fftw_plan_dft_1d( CPD, (fftw_complex *) &(in_1d[0]), (fftw_complex *) &(out_1d[0]), FFTW_FORWARD, FFTW_PATIENT);
//     //plan_forward_1d_r2c = fftw_plan_dft_r2c_1d(CPD,  &(in_r2c[0]), (fftw_complex *) &(out_r2c[0]), FFTW_PATIENT);
//
//
//
//         //NAM: commented this out for large CPD run overnight.
//         //int wisdomExists = fftw_import_wisdom_from_filename("Part2.wisdom");
//         //if (!wisdomExists)
//         //     fmt::print("No wisdom file exists!\n");
//
//
//     plan_forward_1d  =  fftw_plan_dft_1d( CPD, (fftw_complex *) &(in_1d[0]), (fftw_complex *) &(out_1d[0]), FFTW_FORWARD, FFTW_PATIENT);
//     plan_forward_1d_r2c = fftw_plan_dft_r2c_1d(CPD,  &(in_r2c[0]), (fftw_complex *) &(out_r2c[0]), FFTW_PATIENT);
//
//         //if(!wisdomExists)
//         //        fmt::print("Exporting wisdom to file == {:d}\n", fftw_export_wisdom_to_filename("Part2.wisdom"));
//         //
//         //NAM: end commented out section.
//
//
//
//
//
//     ULLI sizeread;
//
//
//         //myTimer.Start(); //NAM
//
//
//
//
//
//
//
//         part2timer.Stop();
//         fmt::print("About to begin multipoles lap. Time elapsed = {:f}\n", part2timer.Elapsed()); fflush(NULL);
//         part2timer.Start();
//
//     int a,b,c;
//     FORALL_REDUCED_MULTIPOLES_BOUND(a,b,c,order) {
//
//         int mab = bm.rmap(a,b,c);
//         int mba = bm.rmap(b,a,c);
//
//         FILE *fpfar;
//         fs::path fpfar_fn = "farderivatives";
//
//                 //Open and unpack farderivatives file.
//         fpfar = fopen(fpfar_fn.c_str(),"rb");
//         assert(fpfar!=NULL);
//         fseek(fpfar, mab*CompressedMultipoleLengthXYZ*sizeof(double)  , SEEK_SET );
//         sizeread = fread(&(FarDerivatives_ab[0]), sizeof(double), CompressedMultipoleLengthXYZ, fpfar);
//         assert(sizeread == CompressedMultipoleLengthXYZ);
//         fclose(fpfar);
//
//         fpfar = fopen(fpfar_fn.c_str(),"rb");
//         assert(fpfar!=NULL);
//         fseek(fpfar, mba*CompressedMultipoleLengthXYZ*sizeof(double) , SEEK_SET );
//         sizeread = fread(&(FarDerivatives_ba[0]), sizeof(double), CompressedMultipoleLengthXYZ, fpfar);
//         assert(sizeread == CompressedMultipoleLengthXYZ);
//         fclose(fpfar);
//
//
//         for(int k=0;k<=CPDHALF;k++)
//             for(int i=0;i<=CPDHALF;i++)
//                 for(int j=0;j<=CPDHALF;j++)  {
//                     ULLI l = linearindex( i+CPDHALF, j+CPDHALF, k+CPDHALF );
//                     if(j<=i) td[l] = FarDerivatives_ab[RINDEXYZ(i,j,k)];
//                     else     td[l] = FarDerivatives_ba[RINDEXYZ(j,i,k)];
//                 }
//
//         for(int k=0;k<=CPDHALF;k++)
//             for(int i=0;i<=CPDHALF;i++)
//                 for(int j=0;j<=CPDHALF;j++)
//                     for(int sx = -1; sx<=1; sx+=2)
//                         for(int sy=-1; sy<=1; sy+=2)
//                             for(int sz=-1; sz<=1; sz+=2) {
//                                 if( sx + sy + sz < 3 ) {
//
//                                     int s = 1;
//                                     if(a%2==1) s *= sign(sx);
//                                     if(b%2==1) s *= sign(sy);
//                                     if(c%2==1) s *= sign(sz);
//
//                                     int xi = sx*i + CPDHALF;
//                                     int yj = sy*j + CPDHALF;
//                                     int zk = sz*k + CPDHALF;
//
//                                     assert(xi>=0); assert(xi<CPD);
//                                     assert(yj>=0); assert(yj<CPD);
//                                     assert(zk>=0); assert(zk<CPD);
//
//                                     ULLI  l = linearindex( i+CPDHALF, j+CPDHALF, k+CPDHALF );
//                                     ULLI ll = linearindex( xi, yj, zk );
//
//                                     if(s>0)
//                                         td[ll] =  td[l];
//                                     else
//                                         td[ll] = -td[l];
//
//                                 }
//                             }
//
//         // now write back into the global array
//
//         int m = mab;
//
//         for(int i=0;i<CPD;i++)
//             for(int j=0;j<CPD;j++)
//                 for(int k=0;k<CPD;k++) {
//                     ULLI l = linearFFTW(i-CPDHALF,j-CPDHALF,k-CPDHALF);
//                     ULLI ll = linearindex(i,j,k);
//                     tdprime[l] = td[ll];
//                 }
//
//         double *tmpreal = tdprime;
//
//
//                 part2timer.Stop();
//                 fmt::print("Starting yz fftw {:d} {:d} {:d}. Time elapsed = {:f}\n", a,b,c,part2timer.Elapsed()); fflush(NULL);
//                 part2timer.Start();
//
//         for(int x=0;x<CPD;x++) {
//             for(int y=0;y<CPD;y++) {
//                 for(int z=0;z<CPD;z++) in_r2c[z] = tmpreal[x*CPD*CPD + y*CPD + z];
//                 fftw_execute(plan_forward_1d_r2c);
//                 for(int z=0;z<(CPD+1)/2;z++) tmpD[ x*CPD*(CPD+1)/2 + y*(CPD+1)/2 + z  ] = out_r2c[z];
//             }
//
//             for(int z=0;z<(CPD+1)/2;z++) {
//                 for(int y=0;y<CPD;y++) in_1d[y] = tmpD[ x*CPD*(CPD+1)/2 + y*(CPD+1)/2 + z ];
//                  fftw_execute(plan_forward_1d);
//                 for(int y=0;y<CPD;y++) tmpD[ x*CPD*(CPD+1)/2 + y*(CPD+1)/2 + z ] = out_1d[y];
//             }
//         }
//
//
//                 part2timer.Stop();
//                 fmt::print("Finished yz fftw {:d} {:d} {:d}. Time elapsed = {:f}\n",a,b,c, part2timer.Elapsed()); fflush(NULL);
//                 part2timer.Start();
//
//         for(int z=0;z<(CPD+1)/2;z++)
//             for(int y=0;y<CPD;y++) {
//                 for(int x=0;x<CPD;x++) in_1d[x] = tmpD[ x*CPD*(CPD+1)/2 + y*(CPD+1)/2 + z ];
//                  fftw_execute(plan_forward_1d);
//                 for(int x=0;x<CPD;x++) tmpD[  x*CPD*(CPD+1)/2 + y*(CPD+1)/2 + z ] = out_1d[x];
//             }
//
//
//         for(int z=0;z<(CPD+1)/2;z++) {
//             if( ((a+b+c)%2) == 0 )
//                 for(int x=0;x<(CPD+1)/2;x++)
//                     for(int y=0;y<=x;y++)
//                         tmpreal[ RINDEXY(x,y) ] = real(conj(tmpD[x*CPD*(CPD+1)/2 + y*(CPD+1)/2 + z]));
//             else
//                 for(int x=0;x<(CPD+1)/2;x++)
//                     for(int y=0;y<=x;y++)
//                         tmpreal[ RINDEXY(x,y) ] = imag(conj(tmpD[x*CPD*(CPD+1)/2 + y*(CPD+1)/2 + z]));
//
//             FILE *fp;
//             fs::path fn = fmt::format("fourierspace_{:d}_{:d}_{:d}_{:d}_{:d}",cpd,order,inner_radius,far_radius,z);
//             fp = fopen(fn.c_str(),"r+b");
//             assert(fp!=NULL);
//             fseek(fp, m*CompressedMultipoleLengthXY*sizeof(double), SEEK_SET );
//             fwrite( &(tmpreal[0]), sizeof(double), CompressedMultipoleLengthXY, fp);
//             fclose(fp);
//         }
//
//
//                 part2timer.Stop();
//                 fmt::print("Finished writing {:d} {:d} {:d}. Time elapsed = {:f}\n",a,b,c, part2timer.Elapsed()); fflush(NULL);
//                 part2timer.Start();
//     }
//
// }














//calculate Fast Fourier transform of far derivatives tensor and store. 
void Part2(int order, int inner_radius, int far_radius, int MultipoleStart) {


        STimer part2timer;
        part2timer.Start();

    fftw_init_threads();


    basemultipoles bm(order);
    int cpd = CPD;

    ULLI fdsize = sizeof(double) * CompressedMultipoleLengthXYZ;
    fmt::print("Allocating {:d} GB for FarDerivatives_ab and _ba in Part2\n", (int) (2*fdsize/(1<<30)) );
    double *FarDerivatives_ab = (double *) malloc(fdsize);
    double *FarDerivatives_ba = (double *) malloc(fdsize);
    assert( FarDerivatives_ab != NULL );
    assert( FarDerivatives_ba != NULL );

    ULLI tdsize = sizeof(double)*CPD3;
    fmt::print("Allocating {:d} GB for tdprime in Part2\n", (int) (2*tdsize/(1<<30)) );
    double *tdprime = (double *) malloc(tdsize);
    double *td = (double *) malloc(tdsize);
    assert(tdprime!=NULL);
    assert(td!=NULL);

    // ULLI tmpDsize = sizeof(Complex)*CPD3;
    ULLI tmpDsize = sizeof(Complex)*CPD*CPD*(CPD+1)/2;
    fmt::print("Allocating {:d} GB for tmpD in Part2\n", (int) (tmpDsize/(1<<30)) );
    Complex *tmpD   = (Complex *) malloc(tmpDsize);
    assert(tmpD!=NULL);

    int rml = (order+1)*(order+1);
    int MultipoleBuffer = 9;     // We're going to save up this many multipoles before writing out.
    int MultipoleNumber[rml];
    ULLI diskbuffersize = sizeof(double)*MultipoleBuffer*(CPD+1)/2*CompressedMultipoleLengthXY;
        // This is enough to hold all the z slabs for these multipoles
    fmt::print("Allocating {:d} GB for diskbuffer in Part2, enough for {:d} multipoles, {:d}\n", (int) (diskbuffersize/(1<<30)), MultipoleBuffer, diskbuffersize/sizeof(double));
    double *diskbuffer   = (double *) malloc(diskbuffersize);
    assert(diskbuffer!=NULL);
    int MultipoleCount = 0;

    if (MultipoleStart!=0) {
        if (MultipoleStart%MultipoleBuffer>0) {
             // We only deal with entire buffers, round down
             fmt::print("Apparent size of fourier file was {:d}, but didn't divide buffer size {:d} evenly.  Adjusting\n", MultipoleStart, MultipoleBuffer);
             MultipoleStart = (MultipoleStart/MultipoleBuffer)*MultipoleBuffer;
        }
    }
    fmt::print("We will start at multipole {:d}\n", MultipoleStart);


    // The X-Y FFTs are just a simple 2D C->C FFT
    fftw_plan plan_forward_2d;
    Complex *in_2d;
    Complex *out_2d;

    // We're going to setup to do CPD FFTs at once, each of size CPD
    ULLI CPDpad = (CPD/16+1)*16;    // Just want to pad out to separate these FFTs.
    // We are supposed to put the [CPD][CPD] info into [CPD][CPDpad] space.  FFT is on the rightmost index.
    fftw_plan plan_forward_1d_r2c;
    double   *in_r2c;
    Complex *out_r2c;

    in_2d = new Complex[CPD*CPD];
    out_2d = new Complex[CPD*CPD];
    in_r2c = new double[CPDpad*CPD];
    out_r2c = new Complex[CPDpad*CPD];
    // This is overkill: only needed half the elements, i.e., CPDpad/2; don't need to economize


    // Old 1-d code
    //plan_forward_1d  =  fftw_plan_dft_1d( CPD, (fftw_complex *) &(in_1d[0]), (fftw_complex *) &(out_1d[0]), FFTW_FORWARD, FFTW_PATIENT);
    //plan_forward_1d_r2c = fftw_plan_dft_r2c_1d(CPD,  &(in_r2c[0]), (fftw_complex *) &(out_r2c[0]), FFTW_PATIENT);

        fmt::print("Planning fftw with omp\n");
        fftw_plan_with_nthreads(omp_get_max_threads());
        fmt::print("Done planning with {:d} threads\n", omp_get_max_threads());


    // This is the plan to do the 2d CPD*CPD complex-to-complex XY FFTs
    plan_forward_2d  =  fftw_plan_dft_2d( CPD, CPD,
        (fftw_complex *) &(in_2d[0]),
        (fftw_complex *) &(out_2d[0]),
        FFTW_FORWARD, FFTW_PATIENT);

    // This is the plan to do CPD real-to-complex Z FFTs, each of size CPD, with a stride of CPDpad
    int fftw_n[] = { (int) CPD };    // This is how long the FFTs are
    int howmany = CPD;                // How many FFTs we're doing
    int idist, odist; idist = odist = CPDpad;   // This is how the FFTs are separated in memory
    int istride, ostride; istride = ostride = 1;
    // iembed and oembed are just NULL
    plan_forward_1d_r2c  =  fftw_plan_many_dft_r2c( 1, fftw_n, howmany,
        &(in_r2c[0]), NULL, istride, idist,
        (fftw_complex *) &(out_r2c[0]), NULL, ostride, odist,
        FFTW_PATIENT);

        fmt::print("Plans created.\n");




    ULLI sizeread;


        part2timer.Stop();
        fmt::print("About to begin multipoles lap. Time elapsed = {:f}\n", part2timer.Elapsed()); fflush(NULL);
        part2timer.Start();

    int a,b,c;
    FORALL_REDUCED_MULTIPOLES_BOUND(a,b,c,order) {

        int mab = bm.rmap(a,b,c);
        int mba = bm.rmap(b,a,c);
        if (MultipoleCount<MultipoleStart) {
            MultipoleCount++;
            continue;
        }

        FILE *fpfar;
        fs::path fpfar_fn = fmt::format("farderivatives_{:d}_{:d}_{:d}_{:d}",CPD,order,inner_radius,far_radius);
        //Open and unpack farderivatives file.
        fpfar = fopen(fpfar_fn.c_str(),"rb");
        assert(fpfar!=NULL);
        fseek(fpfar, mab*CompressedMultipoleLengthXYZ*sizeof(double)  , SEEK_SET );
        sizeread = fread(&(FarDerivatives_ab[0]), sizeof(double), CompressedMultipoleLengthXYZ, fpfar);
        assert(sizeread == CompressedMultipoleLengthXYZ);
        fclose(fpfar);

        fpfar = fopen(fpfar_fn.c_str(),"rb");
        assert(fpfar!=NULL);
        fseek(fpfar, mba*CompressedMultipoleLengthXYZ*sizeof(double) , SEEK_SET );
        sizeread = fread(&(FarDerivatives_ba[0]), sizeof(double), CompressedMultipoleLengthXYZ, fpfar);
        assert(sizeread == CompressedMultipoleLengthXYZ);
        fclose(fpfar);
                
                part2timer.Stop();
                // fmt::print("Done reading farderivs {:d} {:d} {:d}. Time elapsed = {:f}\n", a,b,c,part2timer.Elapsed()); fflush(NULL);
                part2timer.Start();

        // We've loaded in [k][i][j] with j<=i for one [m], only in first octant.
        // First, we reflect [i][j] to fill the first octant.


        // fmt::print("A: {:d} {:d} {:d} {:d}: \n", a,b,c, order); fflush(NULL);
                
        #pragma omp parallel for schedule(static)
        for(int k=0;k<=CPDHALF;k++)
            for(int i=0;i<=CPDHALF;i++)
                for(int j=0;j<=CPDHALF;j++)  {
                    ULLI l = linearindex( i+CPDHALF, j+CPDHALF, k+CPDHALF );
                    if(j<=i) td[l] = FarDerivatives_ab[RINDEXYZ(i,j,k)];
                    else     td[l] = FarDerivatives_ba[RINDEXYZ(j,i,k)];
                }
                                
                                
        fmt::print("."); fflush(NULL);

        // Now we fill the other octants

        #pragma omp parallel for schedule(static)
        for(int i=0;i<=CPDHALF;i++)
            for(int j=0;j<=CPDHALF;j++)
                for(int sx = -1; sx<=1; sx+=2)
                    for(int sy=-1; sy<=1; sy+=2)
                        for(int sz=-1; sz<=1; sz+=2) {
                            if( sx + sy + sz < 3 ) {

                                int s = 1;
                                if(a%2==1) s *= sign(sx);
                                if(b%2==1) s *= sign(sy);
                                if(c%2==1) s *= sign(sz);

                                int xi = sx*i + CPDHALF;
                                int yj = sy*j + CPDHALF;

                                // assert(xi>=0); assert(xi<CPD);
                                // assert(yj>=0); assert(yj<CPD);

                                // These are the starting indices for k=0
                                ULLI  l = linearindex( i+CPDHALF, j+CPDHALF, CPDHALF );
                                ULLI ll = linearindex( xi, yj, CPDHALF );

                                // If sz=+1, we're copying [CPDHALF,2*CPDHALF] onto itself.
                                // If sz=-1, we're reflecting [CPDHALF,2*CPDHALF] into [CPDHALF,0]
                                if(s>0) {
                                    if (sz>0) for(int k=0;k<=CPDHALF;k++) td[ll++] =  td[l++];
                                         else for(int k=0;k<=CPDHALF;k++) td[ll--] =  td[l++];
                                } else {
                                    if (sz>0) for(int k=0;k<=CPDHALF;k++) td[ll++] = -td[l++];
                                         else for(int k=0;k<=CPDHALF;k++) td[ll--] = -td[l++];
                                }

                                // for(int k=0;k<=CPDHALF;k++,l++,ll+=sz) {
                                    // int zk = sz*k + CPDHALF;
                                    // assert(zk>=0); assert(zk<CPD);
                                    // ULLI  l = linearindex( i+CPDHALF, j+CPDHALF, k+CPDHALF );
                                    // ULLI ll = linearindex( xi, yj, zk );

                                    // We are mapping from z = k+CPDHALF to zk = CPDHALF+k*sz

                                    // if(s>0)
                                      //   td[ll] =  td[l];
                                    // else
                                        // td[ll] = -td[l];
                                // }

                                /* We actually want to put everything in a tdprime spot that is 
                                shifted down by -CPDHALF in each coordinate, and then periodic wrapped.

                                If we had started with td already having things in the correct octant,
                                then these reflections would have 
                                */
                            } // else td[l] = td[l];
                        }

        // now write back into the global array
                                                        
        part2timer.Stop();
        // fmt::print("Done with td population {:d} {:d} {:d}. Time elapsed = {:f}\n", a,b,c,part2timer.Elapsed()); fflush(NULL);
        part2timer.Start();

        int m = mab;

        // Do a CPDHALF 3-d cyclic translation of the values.
        // DJE TODO: With care, this could be done in place, which would save a big array

                #pragma omp parallel for schedule(static)
        for(int i=0;i<CPD;i++)
            for(int j=0;j<CPD;j++) {
                ULLI l = linearFFTW(i-CPDHALF,j-CPDHALF,-CPDHALF);
                ULLI ll = linearindex(i,j,0);
                for(int k=0;k<CPDHALF;k++) tdprime[l++] = td[ll++];
                l -= CPD;    // Need to reset to the beginning of the row
                // ULLI l = linearFFTW(i-CPDHALF,j-CPDHALF,0);
                // ULLI ll = linearindex(i,j,CPDHALF);
                for(int k=CPDHALF;k<CPD;k++) tdprime[l++] = td[ll++];
            }

                // for(int k=0;k<CPD;k++) {
                    // ULLI l = linearFFTW(i-CPDHALF,j-CPDHALF,k-CPDHALF);
                    // ULLI ll = linearindex(i,j,k);
                    // tdprime[l] = td[ll];
                // }


        double *tmpreal = tdprime;

                part2timer.Stop();
                // fmt::print("Starting yz fftw {:d} {:d} {:d}. Time elapsed = {:f}\n", a,b,c,part2timer.Elapsed()); fflush(NULL);
                part2timer.Start();

        for(int x=0;x<CPD;x++) {
                        #pragma omp parallel for schedule(static)
            for(int y=0;y<CPD;y++) {
                for(int z=0;z<CPD;z++) {
                                        in_r2c[y*CPDpad+z] = tmpreal[x*CPD*CPD + y*CPD + z];
                                }
                        }
                    fftw_execute(plan_forward_1d_r2c);   // Z FFT, done R->C
                        #pragma omp parallel for schedule(static)
            for(int y=0;y<CPD;y++) {
                for(int z=0;z<(CPD+1)/2;z++) {
                                        tmpD[ x*CPD*(CPD+1)/2 + y*(CPD+1)/2 + z ] = out_r2c[y*CPDpad+z];
                                }
                        }
                }

                part2timer.Stop();
                // fmt::print("Finished z fftw {:d}. Time elapsed = {:f}\n",m, part2timer.Elapsed()); fflush(NULL);
                part2timer.Start();


        // DJE TODO: The above moves the data from tmpreal=tdprime (which could have been td)
        // into tmpD.  If td had gotten the CPD+1 padding, this could have been in place.

        // Now we have a purely complex problem, but we need to do the X-Y FFTs
        // on the (CPD+1)/2 separate KZ problems.
        // Since we're outputing different files for each KZ, we're going to handle each
        // one at a time -- and write out the result directly.
        // This is important, because we have to do a nasty transpose -- fetching one KZ entry
        // per skewer, and we don't want to do that more than we have to.

        for(int z=0;z<(CPD+1)/2;z++) {
                #pragma omp parallel for schedule(static)
            for(int x=0;x<CPD;x++)
                for(int y=0;y<CPD;y++) in_2d[x*CPD+y] = tmpD[ x*CPD*(CPD+1)/2 + y*(CPD+1)/2 + z ];
                // This gather requires a ton of random access!
            fftw_execute(plan_forward_2d);    // X-Y FFT, done C->C
                    // Now the results are in out_2d[x*CPD+y]


/* OLD CODE
            // Put pragma here
            for(int z=0;z<(CPD+1)/2;z++)
        for(int x=0;x<CPD;x++) {
                for(int y=0;y<CPD;y++) tmpD[ x*CPD*(CPD+1)/2 + y*(CPD+1)/2 + z ] = out_1d[z*CPDpad+y];
        }

        for(int z=0;z<(CPD+1)/2;z++) {
            // Put pragma here
            for(int y=0;y<CPD;y++)
                for(int x=0;x<CPD;x++) in_1d[x] = tmpD[ x*CPD*(CPD+1)/2 + y*(CPD+1)/2 + z ];
            fftw_execute(plan_forward_1d);  // X FFT done, C->C
            // Put pragma here
            for(int y=0;y<CPD;y++)
                for(int x=0;x<CPD;x++) tmpD[  x*CPD*(CPD+1)/2 + y*(CPD+1)/2 + z ] = out_1d[x];
        }

        for(int z=0;z<(CPD+1)/2;z++) 
END OLD CODE */


            // We're going to store the results here, so that we can gather up some multipoles
            // before writing out.
            ULLI offset = ((MultipoleCount%MultipoleBuffer)*(CPD+1)/2 + z)
                                *CompressedMultipoleLengthXY;
            double *buf = diskbuffer+offset;
            // if (z==0) {
            //     fmt::print("Saving Multipole {:d} (position {:d} in file) from buffer {:d}; offset {:d}\n", MultipoleCount, m, MultipoleCount%MultipoleBuffer, offset); fflush(NULL);
            // }

            if( ((a+b+c)%2) == 0 )
                                #pragma omp parallel for schedule(static)
                for(int x=0;x<(CPD+1)/2;x++)
                    for(int y=0;y<=x;y++) 
                        buf[ RINDEXY(x,y) ] = real(conj(out_2d[x*CPD+y]));
                        // tmpreal[ RINDEXY(x,y) ] = real(conj(tmpD[x*CPD*(CPD+1)/2 + y*(CPD+1)/2 + z]));
            else
                                #pragma omp parallel for schedule(static)
                for(int x=0;x<(CPD+1)/2;x++)
                    for(int y=0;y<=x;y++)
                        buf[ RINDEXY(x,y) ] = imag(conj(out_2d[x*CPD+y]));
                        // tmpreal[ RINDEXY(x,y) ] = imag(conj(tmpD[x*CPD*(CPD+1)/2 + y*(CPD+1)/2 + z]));
        }

                part2timer.Stop();
                // fmt::print("Finished xy FFT {:d}. Time elapsed = {:f}\n",MultipoleCount, part2timer.Elapsed()); fflush(NULL);
                part2timer.Start();

        // Now we're done with this multipole; consider whether to write out
        MultipoleNumber[MultipoleCount] = m;
        MultipoleCount++;

        if (MultipoleCount%MultipoleBuffer>0 && MultipoleCount<rml) continue;  // Go on to next multipole

        int nMult = (MultipoleCount-1)%MultipoleBuffer+1;
                // fmt::print("Preparing to write {:d} multipoles\n",nMult); fflush(NULL);
        // This is the number of multipoles we have saved up.
        // If we ended at mod 0, then we must have a full set.

        // Here we're doing highly threaded output.  Will the disk like this QueueDepth?
        // Optional pragma!
        // #pragma omp parallel for schedule(static)
        for(int z=0;z<(CPD+1)/2;z++) {
            // We loop over files
            FILE *fp;
            fs::path fn = fmt::format("fourierspace_{:d}_{:d}_{:d}_{:d}_{:d}",cpd,order,inner_radius,far_radius,z);
            fp = fopen(fn.c_str(),"r+b");
            assert(fp!=NULL);
            for (int j=0; j<nMult; j++) {
                int thisMult = MultipoleCount-nMult+j;
                // if (z==0) {
                //     fmt::print("Writing Multipole {:d} (position {:d} in file)from buffer {:d}\n", thisMult, MultipoleNumber[thisMult], thisMult%MultipoleBuffer); fflush(NULL);
                // }
                fseek(fp, MultipoleNumber[thisMult]*CompressedMultipoleLengthXY*sizeof(double), SEEK_SET );
                double *buf = diskbuffer+((thisMult%MultipoleBuffer)*(CPD+1)/2 + z)
                                *CompressedMultipoleLengthXY;
                fwrite( &(buf[0]), sizeof(double), CompressedMultipoleLengthXY, fp);
                // fseek(fp, m*CompressedMultipoleLengthXY*sizeof(double), SEEK_SET );
                // fwrite( &(tmpreal[0]), sizeof(double), CompressedMultipoleLengthXY, fp);
            }
            fclose(fp);
        }

                part2timer.Stop();
                // fmt::print("Finished writing {:d} {:d} {:d}. Time elapsed = {:f}\n",a,b,c, part2timer.Elapsed()); fflush(NULL);
                part2timer.Start();

    }


    delete[] in_2d;
    delete[] out_2d;
    delete[] in_r2c;
    delete[] out_r2c;
}


int main(int argc, char **argv) {
        
#if QUAD_DOUBLE
    fmt::print("QUAD_DOUBLE = 1. Using 256-bit precision.\n");
#else
    fmt::print("QUAD_DOUBLE = 0. Using 64-bit precision.\n");
#endif
    //myTimer.Start();
    
    //unpack user inputs.
    if( argc!=5 && argc!=6 ) {
        fmt::print("Usage: {:s} CPD ORDER INNERRADIUS FARRADIUS [slab]\n", argv[0]);
        fmt::print("Slab number is optional\n");
        exit(1);
    }

    CPD = atoi(argv[1]);
    int cpd = CPD;
    CPD3 = CPD*CPD*CPD;
    fmt::print("CPD = {:d} \n", cpd);
    assert(CPD%2==1);

    CPDHALF = (CPD-1)/2;

    CompressedMultipoleLengthXY  = ((1+CPD)*(3+CPD))/8;
    CompressedMultipoleLengthXYZ = ((1+CPD)*(1+CPD)*(3+CPD))/16;

    int order = atoi(argv[2]);
    fmt::print("order = {:d} \n", order );

    int inner_radius = atoi(argv[3]);
    assert(inner_radius <= (CPD-1)/2 );
    fmt::print("inner_radius = {:d} \n", inner_radius );

    int far_radius = atoi(argv[4]);
    assert( (far_radius >= 1 && far_radius <= 8) || (far_radius == 16));
    fmt::print("far_radius = {:d} \n", far_radius );

    int slabnumber = -1;
    if (argc==6) {
        slabnumber = atoi(argv[5]);
        assert( slabnumber>=0 && slabnumber<=CPDHALF);
        fmt::print("Doing slab {:d} only.", slabnumber);
    }

    int rml = (order+1)*(order+1);

    //check if the derivatives tensor we're about to calculate already stored on disk? If so, don't repeat the work! 
    fs::path fn = fmt::format("fourierspace_{:d}_{:d}_{:d}_{:d}_{:d}",cpd,order,inner_radius,far_radius, (cpd+1)/2-1);
    // This is the last file

    int MultipoleStart;
    FILE *fp = fopen(fn.c_str(),"rb");
    fmt::print(stderr, "Trying to find derivativesfile={} on disk\n", fn);
    if(fp==NULL) {
        // Derivatives don't exist
        MultipoleStart = 0;
    } else {
        // Derivatives files do exist
        MultipoleStart = floor(fs::file_size(fn)/(CompressedMultipoleLengthXY)/sizeof(double));
        // This uses the file size to determine how many multipoles were completed in the last z file
        if (MultipoleStart==rml) {
            fmt::print("Derivatives already present \n");
            exit(0);
        }
        fmt::print("Derivative files were started, but appear to contain only {:d} of {:d} multipoles\n",
            MultipoleStart, rml);
        // If the fourier files are corrupted, delete them to start over
    }
        
    //myTimer.Stop();
    //fmt::print("Time spent in set up = {:f} \n", myTimer.Elapsed());
    //myTimer.Clear();
    //myTimer.Start();
    
    //For more detailed comments, see individual functions, above. 
    //calculate the derivatives tensor for every cell vector, \mathcal{D}^{ABC}_{jkl}. (Eqn. 4.1). 
    FormDerivatives(inner_radius, order, far_radius, slabnumber);
    if (slabnumber>=0) exit(0);
    //myTimer.Stop();
    //double FormDerivsRuntime = myTimer.Elapsed();
    //fmt::print("Time spent in FormDerivatives = {:f} \n", FormDerivsRuntime);
    //myTimer.Clear();
    //myTimer.Start();
    //double mergeDerivsRuntime;
        
        
    
    FILE *fpfar;
    fs::path fpfar_fn = fmt::format("farderivatives_{:d}_{:d}_{:d}_{:d}",CPD,order,inner_radius,far_radius);
    fpfar = fopen(fpfar_fn.c_str(),"rb");
        
    if (fpfar == NULL){    
        // Next we merge the individual files, effectively doing a big transpose
        ULLI fdsize = sizeof(double) * rml*CompressedMultipoleLengthXYZ;
                
        //calculate total far derivatives tensor for entire simulation box by merging components for every cell. {\mathcal{D}^{ABC}_{jkl}} = \mathcal{D}^{ABC}.
        fmt::print("Allocating {:d} GB before MergeDerivatives\n", (int) (fdsize/(1<<30)) );
        double *FarDerivatives = (double *) malloc(fdsize);
        assert(FarDerivatives != NULL);
        MergeDerivatives(inner_radius, order, far_radius, FarDerivatives);
        //myTimer.Stop();
        //mergeDerivsRuntime = myTimer.Elapsed();
        //fmt::print("Time spent in Merge Derivatives = {:f} \n", mergeDerivsRuntime);
        //myTimer.Clear();
        //myTimer.Start();
        
        // write out FarDerivatives file
        fpfar = fopen(fpfar_fn.c_str(),"wb");
        assert(fpfar!=NULL);
        fwrite(&(FarDerivatives[0]), sizeof(double), rml*CompressedMultipoleLengthXYZ, fpfar); 
        free(FarDerivatives);
    }
    fclose(fpfar);


    //create empty fourier files to store the fourier transforms of the far derivatives in, to be calculated in Part2. 
    if (MultipoleStart==0) CreateFourierFiles(order,inner_radius, far_radius);
    //myTimer.Stop();
    //double CFFRuntime = myTimer.Elapsed();
    //fmt::print("Time spent in CreateFourierFiles = {:f} \n", CFFRuntime);
    //myTimer.Clear();
    //myTimer.Start();
    
    //calculate Fast Fourier transforms of far derivatives tensor and store. 
    Part2(order, inner_radius, far_radius, MultipoleStart);
    //myTimer.Stop();
    //double Part2Runtime = myTimer.Elapsed();
    //fmt::print("Time spent in Part2 (not including fftw plan)= {:f} \n", Part2Runtime);
    //myTimer.Clear();
    //myTimer.Start();
    
    
    
    //myTimer.Stop();
    //double totalRuntime = myTimer.Elapsed();
    //fmt::print("Total time elapsed: {} \n", totalRuntime);
    
    //FILE * outFile;
    //outFile = fopen("/home/nam/AbacusProject/timing/alan_CreateDerivatives_LoopOptimization_MultipoleOrder_vs_FR_timing.txt", "a");
    //fmt::print(outFile, "{:d} {:d} {:d} {:d} {:d} {:f} {:f} {:f} {:f}\n", 24, CPD, order, inner_radius, far_radius, FormDerivsRuntime, mergeDerivsRuntime, CFFRuntime, Part2Runtime);
    //fclose (outFile); 
    
    fmt::print("\nRUN COMPLETE: double precision = {:d}, InfDerivSum maxLoopOrder = {:d}, CPD = {:d}, order = {:d}, innerRadius = {:d}, farRadius = {:d}.\n",
        1-QUAD_DOUBLE, MAXLOOPORDER, cpd,order,inner_radius,far_radius);
}
