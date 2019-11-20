#ifndef _IO_PACKN_H_
#define _IO_PACKN_H_
#include <stdint.h>
#include "../particle.h"

typedef struct {
    double SofteningLength;
    double ppd;
    uint64_t ppd_int;
    long long int NP;
    double H0;
    double Omega_M;
    double Omega_DE;
    char OutputFormat[1024];
    double ScaleFactor;
    double BoxSizeHMpc;
    double ParticleMassHMsun;
    double VelZSpace_to_kms;
    double w0;
    double wa;
    int hMpc;
    long long int MaxPID;
    int packN;
} packNheader;

int packN_get_header_size(char* filename);
void packN_validate_file(int packN, char *filename, uint64_t *header_size, uint64_t *data_size);
void packN_parse_header_info(packNheader *header);

void load_particles_packN(int packN, char *filename, struct particle **p, int64_t *num_p);

#endif
