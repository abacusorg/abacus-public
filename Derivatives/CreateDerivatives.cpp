#include "header.cpp"
#include "threevector.hh"
#include "fftw3.h"
#define MAXORDER 16
#include "quad_double.cpp"

typedef ThreeVector<qd_real> qd_real3;

#include "../Multipoles/basemultipoles.cpp"
#include "derivatives.cpp"
#include "order32derivatives.cpp"

int CPD;
int CPDHALF;

int CompressedMultipoleLengthXY;
int CompressedMultipoleLengthXYZ;

double *FarDerivatives;
double *tdprime;
double *td;

double *tmpreal;
Complex *tmpD;


int linearindex(int x, int y, int z) { return x*CPD*CPD + y*CPD + z; }

int linearFFTW(int i, int j, int k) {
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

void Part2(int order); 

void FormDerivatives(int inner_radius, int order, int far_radius) {
    double rdfar[32][1024];
    double rdnear[32][1024];
    double FD[32][1024];

    Order32Derivatives OD(far_radius);

    Derivatives RD(order);

    int reducedmultipolelength = (order+1)*(order+1);

    for(int k=0;k<=CPDHALF;k++) {
        for(int i=0;i<=CPDHALF;i++) {
            #pragma omp parallel for schedule(dynamic,1)
            for(int j=0;j<=i;j++) {
                int p = omp_get_thread_num();

                for(int m=0;m<reducedmultipolelength;m++) FD[p][m] = 0;

                double3 r;
                r.x = i; r.x = r.x/CPD;
                r.y = j; r.y = r.y/CPD;
                r.z = k; r.z = r.z/CPD;

                if( (abs(i)>inner_radius) || (abs(j)>inner_radius) || (abs(k)>inner_radius) ) {
                    RD.ReducedDerivatives(r,&(FD[p][0]));
                }

                OD.Derivative( r, &(rdnear[p][0]), &(rdfar[p][0]),  order );

                for(int m=0;m<reducedmultipolelength;m++) 
                    FarDerivatives[m*CompressedMultipoleLengthXYZ + RINDEXYZ(i,j,k) ] = FD[p][m] + rdnear[p][m] + rdfar[p][m];
            }
        }
    }

    Part2(order);

}

void Part2(int order) { 
    
    basemultipoles bm(order);

    int a,b,c;
    FORALL_REDUCED_MULTIPOLES_BOUND(a,b,c,order) {

        int m = bm.rmap(a,b,c);

        for(int k=0;k<=CPDHALF;k++)
            for(int i=0;i<=CPDHALF;i++)
                for(int j=0;j<=CPDHALF;j++)  {
                    int l = linearindex( i+CPDHALF, j+CPDHALF, k+CPDHALF );
                    if(j<=i) td[l] = FarDerivatives[             m*CompressedMultipoleLengthXYZ + RINDEXYZ(i,j,k)];
                    else     td[l] = FarDerivatives[bm.rmap(b,a,c)*CompressedMultipoleLengthXYZ + RINDEXYZ(j,i,k)];
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

                                    int  l = linearindex( i+CPDHALF, j+CPDHALF, k+CPDHALF );
                                    int ll = linearindex( xi, yj, zk );

                                    if(s>0) 
                                        td[ll] =  td[l];
                                    else
                                        td[ll] = -td[l];

                                }
                            }

        // now write back into the global array 

        for(int i=0;i<CPD;i++)
            for(int j=0;j<CPD;j++)
                for(int k=0;k<CPD;k++) {
                    int l = linearFFTW(i-CPDHALF,j-CPDHALF,k-CPDHALF);
                    int ll = linearindex(i,j,k);
                    tdprime[l] = td[ll];
                }

        FILE *fp;
        char fn[1024];
        sprintf(fn,"realspace%d",m);
        fp = fopen(fn,"w");
        assert(fp!=NULL);
        fclose(fp);
        fp = fopen(fn,"r+b");
        assert(fp!=NULL);
        fwrite( &(tdprime[0]), sizeof(double), CPD*CPD*CPD, fp);
        fclose(fp);
    }
}

void FormFFTD(char *fn, int order, int inner_radius, int far_radius) {

//     basemultipoles bm(order);

    fftw_plan plan_forward_1d;
     Complex in_1d[CPD];
    Complex out_1d[CPD];

    fftw_plan plan_forward_1d_r2c;
    double   in_r2c[CPD];
    Complex out_r2c[CPD];

     plan_forward_1d  =  fftw_plan_dft_1d( CPD, (fftw_complex *) &(in_1d[0]), (fftw_complex *) &(out_1d[0]), FFTW_FORWARD, FFTW_MEASURE);
    plan_forward_1d_r2c = fftw_plan_dft_r2c_1d(CPD,  &(in_r2c[0]), (fftw_complex *) &(out_r2c[0]), FFTW_MEASURE);

    for(int z=0;z<(CPD+1)/2;z++) {
        FILE *fp;
        char fn[1024];
        sprintf(fn,"fourierspace_%d_%d_%d_%d_%d",CPD,order,inner_radius,far_radius,z);
        fp = fopen(fn,"w");
        assert(fp!=NULL);
        fclose(fp);
    }

    int a,b,c;
    int m = 0;
    FORALL_REDUCED_MULTIPOLES_BOUND(a,b,c,order) {
        // int m = bm.rmap(a,b,c);

        FILE *fp;
        char fn[1024];
        sprintf(fn,"realspace%d",m);
        fp = fopen(fn,"rb");
        assert(fp!=NULL);
        int dummy = fread( &(tmpreal[0]), sizeof(double), CPD*CPD*CPD, fp);
        assert(dummy == CPD*CPD*CPD);
        fclose(fp);

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
            sprintf(fn,"fourierspace_%d_%d_%d_%d_%d",CPD,order,inner_radius,far_radius,z);
            fp = fopen(fn,"r+b");
            assert(fp!=NULL);
            fseek(fp, m*CompressedMultipoleLengthXY*sizeof(double), SEEK_SET );
            fwrite( &(tmpreal[0]), sizeof(double), CompressedMultipoleLengthXY, fp);
            fclose(fp);
        }

        m++;
    }
}


int main(int argc, char **argv) {

    if( argc!=5 ) {
        printf("Usage: CreateDerivatives CPD ORDER INNERRADIUS FARRADIUS\n");
        exit(1);
    }

    CPD = atoi(argv[1]);
    printf("CPD = %d \n", CPD);
    assert(CPD%2==1);

    CPDHALF = (CPD-1)/2;

    CompressedMultipoleLengthXY  = ((1+CPD)*(3+CPD))/8;
    CompressedMultipoleLengthXYZ = ((1+CPD)*(1+CPD)*(3*CPD))/16;

    int order = atoi(argv[2]);
    printf("order = %d \n", order );

    int inner_radius = atoi(argv[3]);
    assert(inner_radius <= (CPD-1)/2 );
    printf("inner_radius = %d \n", inner_radius );

    int far_radius = atoi(argv[4]);
    assert( (far_radius==8) || (far_radius==16) );
    printf("far_radius = %d \n", far_radius );

    int rml = (order+1)*(order+1);

    tmpreal = new double[CPD*CPD*CPD];
    tmpD    = new Complex[CPD*CPD*CPD];

    tdprime = new double[CPD*CPD*CPD];
         td = new double[CPD*CPD*CPD];   

    char fn[1024];
    sprintf(fn,"fourierspace_%d_%d_%d_%d_%d",CPD,order,inner_radius,far_radius, 0);

    FILE *fp;
    fp = fopen(fn,"rb");
    fprintf(stderr, "Trying to find derivativesfile=%s on disk\n", fn);
    if(fp!=NULL) {
        printf("Derivatives already present \n");
        exit(1);
    }
    else {
        FarDerivatives = new double[ rml*CompressedMultipoleLengthXYZ ];
        FormDerivatives(inner_radius, order, far_radius);
        FormFFTD(fn,order, inner_radius, far_radius);

        for(int i=0;i<rml;i++) {
            char syscmd[1024];
            sprintf(syscmd,"rm realspace%d", i);
            int dummy = system(syscmd);
            assert(dummy!=-1);
        }
    }

}
