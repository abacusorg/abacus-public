#include "redlack.h"

#ifndef REDLACK_CPP
#define REDLACK_CPP
int Redlack::wrapcellindex(int x) { 
    if(x>=cpd) x-=cpd; if(x<0) x+=cpd; assert(x>=0 && x<cpd); 
    return x; 
}

void Redlack::ClearVariables(void) {
    for(int i=0;i<cpd;i++) { MassSlabX[i]=0; MassSlabY[i]=0; MassSlabZ[i]=0; }
    for(int i=0;i<3*MAXCPD;i++) { redlack[i]=0; }
    globaldipole = double3(0);
    MassTotal = 0.0;
}

Redlack::Redlack(int _cpd) {
    cpd = _cpd;
    invcpd = 1.0/((double) cpd);
    invcpd3 = invcpd*invcpd*invcpd;
    cpdhalf = (cpd-1)/2;
    cpdp1half = (cpd+1)/2;
    fourthirdspi = 4.0/3.0*M_PI;
    double rc = ((double) cpdhalf)/((double) cpd);
    redlackconstant=double3(rc,rc,rc);
    ClearVariables();
}

// If we're multi-node, then we need to add up the MassSlab vectors
void Redlack::GatherRedlack() {
    #ifdef PARALLEL
    MPI_REDUCE_TO_ZERO(MassSlabX, cpd, MPI_DOUBLE, MPI_SUM);
    MPI_REDUCE_TO_ZERO(MassSlabY, cpd, MPI_DOUBLE, MPI_SUM);
    MPI_REDUCE_TO_ZERO(MassSlabZ, cpd, MPI_DOUBLE, MPI_SUM);
    MPI_REDUCE_TO_ZERO((double *)&globaldipole, 3, MPI_DOUBLE, MPI_SUM);
    #endif
    return;
}

void Redlack::convolve_one_RG_direction(double *MassSlab, double *redlack) {
    redlack[0] = 0;   
    for(int ix=-cpdhalf;ix<=cpdhalf;ix++) 
        redlack[0] += MassSlab[wrapcellindex(ix)]*ix;
    redlack[0] /= cpd;
    for(int ix=1; ix<cpd; ix++) 
        redlack[ix] =  redlack[ix-1] + MassSlab[wrapcellindex(ix-cpd/2-1)]-MassTotal/cpd;
    // Subtracting by MassTotal/cpd shifts to the center of the domain.
}

void Redlack::ComputeRedlack(void) {
    MakeRedlack.Start();
    MassTotal = 0.0;
    for (int ix=0; ix<cpd; ix++) MassTotal += MassSlabX[ix];
    // For unit particle mass, MassTotal = NP

    convolve_one_RG_direction(MassSlabX, &(redlack[0]));
    convolve_one_RG_direction(MassSlabY, &(redlack[cpd]));
    convolve_one_RG_direction(MassSlabZ, &(redlack[2*cpd]));

    MakeRedlack.Stop();
}

void Redlack::ApplyRedlack( int slab, 
                            FLOAT3 *slabacceleration, FLOAT3 *slabpos, 
                            int *count, int *offset, FLOAT3 *cc, int np){
	TaylorRedlack.Start();

    #pragma omp parallel for schedule(static)
    for(int y=0;y<cpd;y++){
        for(int z=0;z<node_z_size;z++) {
            int i = y*node_z_size + z;
            double3 dip2( redlack[slab], redlack[cpd+y], redlack[2*cpd + (z + node_z_start)] );
            double3 dc = globaldipole + dip2; //- (redlackconstant+cc[i])*np;
            
            #pragma omp simd
            for(int q = offset[i]; q < offset[i] + count[i]; q++) {
                double3 red = fourthirdspi * ( dc - (slabpos[q]-cc[i])*MassTotal );
                slabacceleration[q] += red;
            }
        }
    }

    TaylorRedlack.Stop();
}

void Redlack::WriteOutAuxiallaryVariables(char *writedirectory) {
    FILE *fp;
    char fn[1024];

    sprintf(fn,"%s/redlack", writedirectory);
    fp = fopen(fn,"wb");
    assert(fp!=NULL);
    fwrite(&(redlack[0]), sizeof(double), MAXCPD*3, fp);
    fclose(fp);

    sprintf(fn,"%s/globaldipole", writedirectory);
    fp = fopen(fn,"wb");
    assert(fp!=NULL);
    double3 gd = -globaldipole;
    fwrite(&(gd), sizeof(double3), 1, fp);
    fwrite(&MassTotal, sizeof(double), 1, fp);
    fclose(fp);
}

void Redlack::ReadInAuxiallaryVariables(char *readdirectory) {
    FILE *fp;
    char fn[1024];
    int rv;

    sprintf(fn,"%s/redlack", readdirectory);
    fp = fopen(fn,"rb");
    assert(fp!=NULL);
    rv = fread(redlack, sizeof(double), MAXCPD*3, fp);
    assert(rv==MAXCPD*3);
    fclose(fp);

    sprintf(fn,"%s/globaldipole", readdirectory);
    fp = fopen(fn,"rb");
    assert(fp!=NULL);
    rv = fread(&globaldipole, sizeof(double3), 1, fp);
    rv = rv && fread(&MassTotal, sizeof(double), 1, fp);
    assert(rv==1);
    fclose(fp);
}
#endif
