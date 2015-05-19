#include "stdheader.cpp"
#include "threevector.hh"
#include "fftw3.h"
#define MAXORDER 16
#include "quad_double.cpp"

typedef ThreeVector<qd_real> qd_real3;

#include "../Multipoles/basemultipoles.cpp"
#include "derivatives.cpp"
#include "order32derivatives.cpp"

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

    int l = indx_i*CPD*CPD + indx_j*CPD + indx_k;
    assert(l>=0);
    assert(l<CPD*CPD*CPD);
    return l;
}

#define RINDEXY(x,y)    (((x)*(x+1))/2 + (y))
#define RINDEXYZ(x,y,z) (((z)*(1+CPD)*(3+CPD))/8 + ((x)*((x)+1))/2 + y)
#define sign(i)         (i>0?1:-1) 


void FormDerivatives(int inner_radius, int order, int far_radius, int slabnumber) {
    double rdfar[32][1024];
    double rdnear[32][1024];
    double FD[32][1024];

    Order32Derivatives OD(far_radius);

    Derivatives RD(order);

    int rml = (order+1)*(order+1);

    ULLI fdsize = sizeof(double) * rml*CompressedMultipoleLengthXY;
    printf("Allocating %d GB\n", (int) (fdsize/(1<<30)) );
    double *FDSlab = (double *) malloc(fdsize);

    for(int k=0;k<=CPDHALF;k++) {
	if (slabnumber>=0 && slabnumber!=k) continue;
	    // We were asked to do only one slab

	// Perhaps this file already exists?  If so, we'll skip.
	char fn[1024];
	sprintf(fn,"slabderiv_%llu_%d_%d_%d_%d__%d",CPD,order,inner_radius,far_radius, 0, k);
	FILE *fp;
	fp = fopen(fn,"rb");
	if(fp!=NULL) {
	    // Yes, we've found them!  Don't repeat the work.
	    fclose(fp);
	    continue;  
	}

	// Not on disk, so we have to make them.
	#pragma omp parallel for schedule(dynamic,1) 
        for(int i=0;i<=CPDHALF;i++) {
            for(int j=0;j<=i;j++) {
                int p = omp_get_thread_num();

                for(int m=0;m<rml;m++) FD[p][m] = 0;

                double3 r;
                r.x = i; r.x = r.x/CPD;
                r.y = j; r.y = r.y/CPD;
                r.z = k; r.z = r.z/CPD;

                if( (abs(i)>inner_radius) || (abs(j)>inner_radius) || (abs(k)>inner_radius) ) {
                    RD.ReducedDerivatives(r,&(FD[p][0]));
                }

                OD.Derivative( r, &(rdnear[p][0]), &(rdfar[p][0]),  order );

                for(int m=0;m<rml;m++) {
                    ULLI idx = m*CompressedMultipoleLengthXY + RINDEXY(i,j);
                    FDSlab[idx] = FD[p][m] + rdnear[p][m] + rdfar[p][m];
                }
            }
        }
	// Now we store FDSlab to disk
	fp = fopen(fn,"wb");
	assert(fp!=NULL);
	fwrite(&(FDSlab[0]), sizeof(double), rml*CompressedMultipoleLengthXY, fp); 
	fclose(fp);
    }
    free(FDSlab);
    return;
}


void MergeDerivatives(int inner_radius, int order, int far_radius, double *FarDerivatives) {
    int rml = (order+1)*(order+1);
    ULLI fdsize = sizeof(double) * rml*CompressedMultipoleLengthXY;
    printf("Allocating %d GB\n", (int) (fdsize/(1<<30)) );
    double *FDSlab = (double *) malloc(fdsize);

    for(int k=0;k<=CPDHALF;k++) {
	// Fetch FDSlab from disk
	char fn[1024];
	sprintf(fn,"slabderiv_%llu_%d_%d_%d_%d__%d",CPD,order,inner_radius,far_radius, 0, k);
	FILE *fp;
	fp = fopen(fn,"rb");
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
    }
    free(FDSlab);
    return;
}



void Part2(int order, int inner_radius, int far_radius) { 
    
    basemultipoles bm(order);
    int cpd = CPD;

    ULLI fdsize = sizeof(double) * CompressedMultipoleLengthXYZ;
    printf("Allocating %d GB\n", (int) (2*fdsize/(1<<30)) );
    double *FarDerivatives_ab = (double *) malloc(fdsize);
    double *FarDerivatives_ba = (double *) malloc(fdsize);
    assert( FarDerivatives_ab != NULL );
    assert( FarDerivatives_ba != NULL );

    ULLI tdsize = sizeof(double)*CPD3;
    printf("Allocating %d GB\n", (int) (2*tdsize/(1<<30)) );
    double *tdprime = (double *) malloc(tdsize);
    double *td = (double *) malloc(tdsize);
    assert(tdprime!=NULL);
    assert(td!=NULL);

    ULLI tmpDsize = sizeof(Complex)*CPD3;
    printf("Allocating %d GB\n", (int) (tmpDsize/(1<<30)) );
    Complex *tmpD   = (Complex *) malloc(tmpDsize);
    assert(tmpD!=NULL);

    fftw_plan plan_forward_1d;
     Complex in_1d[CPD];
    Complex out_1d[CPD];

    fftw_plan plan_forward_1d_r2c;
    double   in_r2c[CPD];
    Complex out_r2c[CPD];

     plan_forward_1d  =  fftw_plan_dft_1d( CPD, (fftw_complex *) &(in_1d[0]), (fftw_complex *) &(out_1d[0]), FFTW_FORWARD, FFTW_MEASURE);
    plan_forward_1d_r2c = fftw_plan_dft_r2c_1d(CPD,  &(in_r2c[0]), (fftw_complex *) &(out_r2c[0]), FFTW_MEASURE);


    ULLI sizeread;

    int a,b,c;
    FORALL_REDUCED_MULTIPOLES_BOUND(a,b,c,order) {
        
        int mab = bm.rmap(a,b,c);
        int mba = bm.rmap(b,a,c);

        FILE *fpfar;
        char fpfar_fn[1024];
        sprintf(fpfar_fn,"farderivatives");

        fpfar = fopen(fpfar_fn,"rb");
        assert(fpfar!=NULL);
        fseek(fpfar, mab*CompressedMultipoleLengthXYZ*sizeof(double)  , SEEK_SET );
        sizeread = fread(&(FarDerivatives_ab[0]), sizeof(double), CompressedMultipoleLengthXYZ, fpfar);
	assert(sizeread == CompressedMultipoleLengthXYZ);
        fclose(fpfar);

        fpfar = fopen(fpfar_fn,"rb");
        assert(fpfar!=NULL);
        fseek(fpfar, mba*CompressedMultipoleLengthXYZ*sizeof(double) , SEEK_SET );
        sizeread = fread(&(FarDerivatives_ba[0]), sizeof(double), CompressedMultipoleLengthXYZ, fpfar);
	assert(sizeread == CompressedMultipoleLengthXYZ);
        fclose(fpfar);


        for(int k=0;k<=CPDHALF;k++)
            for(int i=0;i<=CPDHALF;i++)
                for(int j=0;j<=CPDHALF;j++)  {
                    ULLI l = linearindex( i+CPDHALF, j+CPDHALF, k+CPDHALF );
                    if(j<=i) td[l] = FarDerivatives_ab[RINDEXYZ(i,j,k)];
                    else     td[l] = FarDerivatives_ba[RINDEXYZ(j,i,k)];
                }

        for(int k=0;k<=CPDHALF;k++)
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
                                    int zk = sz*k + CPDHALF;

                                    assert(xi>=0); assert(xi<CPD);
                                    assert(yj>=0); assert(yj<CPD);
                                    assert(zk>=0); assert(zk<CPD);

                                    ULLI  l = linearindex( i+CPDHALF, j+CPDHALF, k+CPDHALF );
                                    ULLI ll = linearindex( xi, yj, zk );

                                    if(s>0) 
                                        td[ll] =  td[l];
                                    else
                                        td[ll] = -td[l];

                                }
                            }

        // now write back into the global array 
        
        int m = mab;

        for(int i=0;i<CPD;i++)
            for(int j=0;j<CPD;j++)
                for(int k=0;k<CPD;k++) {
                    ULLI l = linearFFTW(i-CPDHALF,j-CPDHALF,k-CPDHALF);
                    ULLI ll = linearindex(i,j,k);
                    tdprime[l] = td[ll];
                }

        double *tmpreal = tdprime;

        for(int x=0;x<CPD;x++) {
            for(int y=0;y<CPD;y++) {
                for(int z=0;z<CPD;z++) in_r2c[z] = tmpreal[x*CPD*CPD + y*CPD + z];
                fftw_execute(plan_forward_1d_r2c);
                for(int z=0;z<(CPD+1)/2;z++) tmpD[ x*CPD*(CPD+1)/2 + y*(CPD+1)/2 + z  ] = out_r2c[z];
            }

            for(int z=0;z<(CPD+1)/2;z++) {
                for(int y=0;y<CPD;y++) in_1d[y] = tmpD[ x*CPD*(CPD+1)/2 + y*(CPD+1)/2 + z ];
                 fftw_execute(plan_forward_1d); 
                for(int y=0;y<CPD;y++) tmpD[ x*CPD*(CPD+1)/2 + y*(CPD+1)/2 + z ] = out_1d[y];
            }
        }

        for(int z=0;z<(CPD+1)/2;z++)
            for(int y=0;y<CPD;y++) {
                for(int x=0;x<CPD;x++) in_1d[x] = tmpD[ x*CPD*(CPD+1)/2 + y*(CPD+1)/2 + z ];
                 fftw_execute(plan_forward_1d);
                for(int x=0;x<CPD;x++) tmpD[  x*CPD*(CPD+1)/2 + y*(CPD+1)/2 + z ] = out_1d[x];
            }


        for(int z=0;z<(CPD+1)/2;z++) {
            if( ((a+b+c)%2) == 0 )
                for(int x=0;x<(CPD+1)/2;x++) 
                    for(int y=0;y<=x;y++) 
                        tmpreal[ RINDEXY(x,y) ] = real(conj(tmpD[x*CPD*(CPD+1)/2 + y*(CPD+1)/2 + z]));
            else
                for(int x=0;x<(CPD+1)/2;x++) 
                    for(int y=0;y<=x;y++) 
                        tmpreal[ RINDEXY(x,y) ] = imag(conj(tmpD[x*CPD*(CPD+1)/2 + y*(CPD+1)/2 + z]));

            FILE *fp;
            char fn[1024];
            sprintf(fn,"fourierspace_%d_%d_%d_%d_%d",cpd,order,inner_radius,far_radius,z);
            fp = fopen(fn,"r+b");
            assert(fp!=NULL);
            fseek(fp, m*CompressedMultipoleLengthXY*sizeof(double), SEEK_SET );
            fwrite( &(tmpreal[0]), sizeof(double), CompressedMultipoleLengthXY, fp);
            fclose(fp);
        }
    }

}


void CreateFourierFiles(int order, int inner_radius, int far_radius) {
    for(int z=0;z<(CPD+1)/2;z++) {
        FILE *fp;
        char fn[1024];
        int cpd = CPD;
        sprintf(fn,"fourierspace_%d_%d_%d_%d_%d",cpd,order,inner_radius,far_radius,z);
        fp = fopen(fn,"w");
        assert(fp!=NULL);
        fclose(fp);
    }
}



int main(int argc, char **argv) {

    if( argc!=5 && argc!=6 ) {
        printf("Usage: CreateDerivatives CPD ORDER INNERRADIUS FARRADIUS <slab>\n");
        printf("Slab number is optional\n");
        exit(1);
    }

    CPD = atoi(argv[1]);
    int cpd = CPD;
    CPD3 = CPD*CPD*CPD;
    printf("CPD = %d \n", cpd);
    assert(CPD%2==1);

    CPDHALF = (CPD-1)/2;

    CompressedMultipoleLengthXY  = ((1+CPD)*(3+CPD))/8;
    CompressedMultipoleLengthXYZ = ((1+CPD)*(1+CPD)*(3+CPD))/16;

    int order = atoi(argv[2]);
    printf("order = %d \n", order );

    int inner_radius = atoi(argv[3]);
    assert(inner_radius <= (CPD-1)/2 );
    printf("inner_radius = %d \n", inner_radius );

    int far_radius = atoi(argv[4]);
    assert( (far_radius==8) || (far_radius==16) );
    printf("far_radius = %d \n", far_radius );

    int slabnumber = -1;
    if (argc==6) {
	slabnumber = atoi(argv[5]);
	assert( slabnumber>=0 && slabnumber<=CPDHALF);
	printf("Doing slab %d only.", slabnumber);
    }

    int rml = (order+1)*(order+1);

    char fn[1024];
    sprintf(fn,"fourierspace_%d_%d_%d_%d_%d",cpd,order,inner_radius,far_radius, 0);

    FILE *fp;
    fp = fopen(fn,"rb");
    fprintf(stderr, "Trying to find derivativesfile=%s on disk\n", fn);
    if(fp!=NULL) {
        printf("Derivatives already present \n");
        exit(0);
    }
    FormDerivatives(inner_radius, order, far_radius, slabnumber);
    if (slabnumber>=0) exit(0);
    
    FILE *fpfar;
    char fpfar_fn[1024];

    sprintf(fpfar_fn,"farderivatives");
    fpfar = fopen(fpfar_fn,"rb");
    if (fpfar == NULL){
    //fclose(fpfar);
    
	// We were only asked to do one slab.
    // Next we merge the individual files, effectively doing a big transpose
    ULLI fdsize = sizeof(double) * rml*CompressedMultipoleLengthXYZ;
    printf("Allocating %d GB\n", (int) (fdsize/(1<<30)) );
    double *FarDerivatives = (double *) malloc(fdsize);
    MergeDerivatives(inner_radius, order, far_radius, FarDerivatives);

    // write out FarDerivatives file
    fpfar = fopen(fpfar_fn,"wb");
    assert(fpfar!=NULL);
    fwrite(&(FarDerivatives[0]), sizeof(double), rml*CompressedMultipoleLengthXYZ, fpfar); 
    free(FarDerivatives);
    }
    fclose(fpfar);

    // Now do the FFTs based on seeking on the disk file
    CreateFourierFiles(order,inner_radius, far_radius);
    Part2(order, inner_radius, far_radius);
}
