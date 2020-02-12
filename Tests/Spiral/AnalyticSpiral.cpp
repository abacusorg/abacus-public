#include <fftw3.h>
#include <math.h>
#include <cassert>
#include <stdio.h>
#include <stdlib.h>
class AnalyticSpiral {
public:
    AnalyticSpiral(float Ainitial, float Across, float Astop, float Astep, int grid1d, float fsmooth);
    ~AnalyticSpiral(void);
    void PMintegrate( float Aexpn, float Astep);

    int blitzNX;
    int blitzNG1X;

    double *x;
    double *px;
    
    float Ainitial;
    float Across;
    float Astep;
    float Astop; 
    float fclustered;

    double *g;
    double *Green;

    fftw_complex *rho;
    fftw_complex *fftrho;
    fftw_complex *phi;


    fftw_plan plan_forward;
    fftw_plan plan_backward;

}; 

AnalyticSpiral::AnalyticSpiral( float _Ainitial, float _Across, float _Astop, float _Astep, int _grid1d, float _fsmooth ) {
    Ainitial = _Ainitial;
    Across   = _Across;
    Astop    = _Astop;
    Astep    = _Astep;
    fclustered  = 1.0-_fsmooth;

    blitzNX = _grid1d;
    blitzNG1X = blitzNX - 1;

    rho = new fftw_complex[blitzNX];
    fftrho = new fftw_complex[blitzNX];
    phi = new fftw_complex[blitzNX];
         
    
    x = new double[blitzNX];
    px = new double[blitzNX];
    
    g = new double[blitzNX];
    Green = new double[blitzNX];


    plan_forward = fftw_plan_dft_1d(blitzNX,&(rho[0]), &(fftrho[0]), FFTW_FORWARD, FFTW_ESTIMATE );
    plan_backward = fftw_plan_dft_1d(blitzNX, &(fftrho[0]), &(phi[0]), FFTW_BACKWARD, FFTW_ESTIMATE );


    double kw,kxq,ampx,da,da32;

    kw=2.0*M_PI/blitzNX;


    da = Ainitial - 0.5*Astep;
    da32 = pow(da,1.5);
    double f_growth = (sqrt(1.0+24.0*fclustered)-1.0)/4.0;

    int i;
    for(i=0;i<blitzNX;i++)   {
        kxq=kw*i;
        ampx=1.0/(Across*kw);
        x[i]=i+Ainitial*ampx*sin(kxq);
        if(x[i]<0.0) x[i]+=blitzNX;
        if(x[i]>=blitzNX) x[i]-=blitzNX;
        px[i]=da32*ampx*sin(kxq)*f_growth;
    }

    double Aexpn;

    Aexpn = Ainitial;

    while(Aexpn<Astop) {
      PMintegrate(Aexpn, Astep);
      if(Aexpn+Astep>Astop) Astep=Astop-Aexpn+0.00001;
      Aexpn+=Astep;
    }

      for(i=0;i<blitzNX;i++) { 
        x[i] /= blitzNX;
        px[i] /= blitzNX;
    }


}

AnalyticSpiral::~AnalyticSpiral(void) {
    fftw_destroy_plan(plan_forward );
    fftw_destroy_plan(plan_backward);
}


void AnalyticSpiral::PMintegrate( float Aexpn, float Astep ) {

    int i,ip1,im1,p,l; 
    double xc,tx,dx,cc;
    float ahalf,faexpn,fahalf,ax;

    double wx,gg;

    wx = 2 * M_PI/blitzNX;

    for(i=0;i<blitzNX;i++) {
        l = i;

        if(l==0) Green[i]= 0;
        else {
            if( l > blitzNX/2 ) l = -(blitzNX-l);
            gg = -0.5 * ( cos(wx*l) - 1  );
            Green[i] =  1.0/gg;
        }
    }


    for (i=0;i<blitzNX;i++)  {
        rho[i][0]=0.0;
        rho[i][1]=0;
    }

    for (p=0;p<blitzNX;p++) {
        i=(int) floor(x[p]); xc=(float)i; dx=x[p]-xc; tx=1.0-dx;
        ip1 = (i+1)&blitzNG1X;
        rho[i  ][0] +=tx;
        rho[ip1][0] +=dx;
    }

    fftw_execute(plan_forward);

    for(i=0;i<blitzNX;i++) {
        fftrho[i][0] *= (Green[i]);
        fftrho[i][1] *= (Green[i]);
    }
    
    fftw_execute(plan_backward);

   
    for(i=0;i<blitzNX;i++)  {
        ip1 = (i+1)&blitzNG1X;
        im1 = (i-1)&blitzNG1X;
        g[i] = -(phi[ip1][0]-phi[im1][0]);
    }

    // Aexpn is the current scale factor.
    // For Omega=1 and H0=1, the kickfactor is da/sqrt(a).
    // The driftfactor is da/a**1.5.
    // The code below implements this, accounting for the factor of Aexpn
    // in the 'cc' normalization of the accelerations!
    // The resulting px is the canonical momentum, not the comoving velocity!

    cc = -(3.0*1.0/8.0/Aexpn)/(float)(blitzNX)*fclustered;

    ahalf=Aexpn+0.5*Astep;
    faexpn=sqrt(Aexpn)*Astep;
    fahalf=sqrt(ahalf)*Astep/(ahalf*ahalf);

    for(p=0;p<blitzNX;p++) {
        i=(int) floor(x[p]); 
        xc=(float)i; 
        dx=x[p]-xc; 
        tx=1.0-dx;

        ip1 = (i+1)&blitzNG1X;

        ax = (tx*g[i] + dx*g[ip1])*0.5*cc;
        px[p]+=faexpn*ax;
        x[p]+=fahalf*px[p];

        if(x[p]<0.0) x[p]+=blitzNX;
        if(x[p]>=blitzNX) x[p]-=blitzNX;
    }
}


int main(int argc, char **argv) {

    if(argc!=5) {
        printf("usage: AnalyticSpiral <ainitial> <across> <afinal> <fsmooth>\n");
        exit(1);
    }

    double ainitial   = atof(argv[1]);
    double across   = atof(argv[2]);
    double afinal    = atof(argv[3]);
    double fsmooth    = atof(argv[4]);
    assert(fsmooth>=0 && fsmooth<1.0);

    AnalyticSpiral AS(ainitial, across, afinal, 0.0001, 8192, fsmooth);

    FILE *fp;
    fp = fopen("analytic","w");
    assert(fp!=NULL);
    for(int p=0;p<8192;p++) {
        fprintf(fp,"%e %e \n", AS.x[p]-0.5, AS.px[p] );
    }
    fclose(fp);
}
