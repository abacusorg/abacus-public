#ifndef REDLACK_H
#define REDLACK_H

class Redlack {
public:
    Redlack(int cpd);
    ~Redlack(void) {}

    double MassSlabX[MAXCPD];
    double MassSlabY[MAXCPD];
    double MassSlabZ[MAXCPD];

    double redlack[3*MAXCPD];
    double3 globaldipole;

    void ApplyRedlack(int slab, FLOAT3 *slabacceleration, FLOAT3 *slabpos,
                        int *count, int *offset, int *ghost_offsets,
                        FLOAT3 *cc, int np);
    void ClearVariables(void);

    void GatherRedlack(void);
    void ComputeRedlack(void);
    void convolve_one_RG_direction(double *MassSlab, double *redlack);

    void ReadInAuxiallaryVariables(char *readdirectory);
    void WriteOutAuxiallaryVariables(char *writedirectory);

    int wrapcellindex(int x); 

    int cpdhalf;
    int cpdp1half;
    int cpd;
    double invcpd;
    double invcpd3;
    double fourthirdspi;
    double3 redlackconstant;
    double MassTotal;
    STimer MakeRedlack;
    STimer TaylorRedlack;
}; 

#endif
