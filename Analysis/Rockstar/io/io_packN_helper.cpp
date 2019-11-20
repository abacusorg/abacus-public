#include "ParseHeader.hh"
#include <stdlib.h>
#include <cassert>
#include <cstring>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <fcntl.h>
#include <stdint.h>
#include <errno.h>
#include "file.h"
#include "file.cpp"
#include "read_dio.cpp"
#include "cell_header.h"
#include "packN_storage.cpp"

extern "C" {
#include "io_packN.h"
#include "../check_syscalls.h"
#include "../config_vars.h"
#include "../bounds.h"
}

// Define the header fields that we are interested in
template <>
void ParseHeader::register_vars(packNheader &P) {
    // SofteningLength is duplicated in the .par and state headers
    // We use a stride of 0 here to trick ParseHeader into overwriting SofteningLength
    // with the last-occurring value
    installvector("SofteningLength", &(P.SofteningLength), 2, 0, MUST_DEFINE);
    installscalar("ppd", P.ppd, MUST_DEFINE);
    installscalar("NP", P.NP, MUST_DEFINE);
    installscalar("H0", P.H0, MUST_DEFINE);
    installscalar("Omega_M", P.Omega_M, MUST_DEFINE);
    installscalar("Omega_DE", P.Omega_DE, MUST_DEFINE);
    installscalar("OutputFormat", P.OutputFormat, MUST_DEFINE);
    installscalar("ScaleFactor", P.ScaleFactor, MUST_DEFINE);
    installscalar("BoxSizeHMpc", P.BoxSizeHMpc, MUST_DEFINE);
    installscalar("ParticleMassHMsun", P.ParticleMassHMsun, MUST_DEFINE);
    installscalar("VelZSpace_to_kms", P.VelZSpace_to_kms, MUST_DEFINE);
    installscalar("w0", P.w0, MUST_DEFINE);
    installscalar("wa", P.wa, MUST_DEFINE);
    installscalar("hMpc", P.hMpc, MUST_DEFINE);

    P.MaxPID = -1;
    installscalar("MaxPID", P.MaxPID, DONT_CARE);
}

extern "C" int packN_get_header_size(char* filename){
    HeaderStream hstream(filename);
    hstream.SkipHeader();
    hstream.Close();
    
    int header_size = hstream.bufferlength;
    
    return header_size;
}

// Loads the header from "filename" into "header",
void packN_load_header(char* filename, packNheader *header){
    ParseHeader PH;
    PH.register_vars(*header);
    HeaderStream hstream(filename);
    PH.ReadHeader(hstream);
    hstream.Close();
}

// Use direct IO to read the data portion of the file into a buffer
char *packN_load_data(int packN, char* filename, uint64_t header_size, uint64_t data_size, char **pid_buffer, size_t *pid_data_size){
    char *buffer = (char *) check_posix_memalign(4096, data_size);
    int ramdisk = 1;  // filesystems like El Gato /rsgrps seem to strongly prefer fread
	ReadDirect rd(ramdisk,1024*1024);
	rd.BlockingRead(filename, (char*)buffer, data_size, header_size);

    if(packN == 9 && !IGNORE_PARTICLE_IDS){
        // Pack9 PIDs are in another file
        size_t fnlen = strnlen(filename,1024);
        assert(strcmp(filename +(fnlen-4), ".dat") == 0);

        char pidfn[fnlen+6];  // + "_pids" and null byte
        strncpy(pidfn, filename, fnlen-4);
        strcpy(pidfn+fnlen-4, "_pids.dat");

        // We can do this here because there's no header
        struct stat pid_stat;
        int ret = stat(pidfn, &pid_stat);
        assert(ret == 0);
        *pid_data_size = pid_stat.st_size;
        *pid_buffer = (char*)check_posix_memalign(4096, *pid_data_size);

        rd.BlockingRead(pidfn, *pid_buffer, *pid_data_size, 0);
    }
    
    return buffer;
}

template<int N>
uint64_t _load_particles_packN(packNheader &header, FILE* buffer_file, FILE* pid_buffer_file, struct particle *p)
{   
    // Now read the particles from the buffer
    cell_header current_cell;
    packN<N> pN;
    current_cell.vscale = 0;   // Make this illegal

    uint64_t id;
    uint64_t idx, idy, idz; // particle lattice indices
    uint64_t ppd = (uint64_t)(header.ppd+.5);
    uint64_t ds = DOWNSAMPLE;
    uint64_t i = 0;
    while (fread(&pN, sizeof(packN<N>), 1, buffer_file)==1) {
        if (pN.iscell()) {
            current_cell = pN.unpack_cell();
        } else {
            assert(current_cell.islegal());
            float *pos = p[i].pos;
            float *vel = p[i].pos+3;
            pN.unpack(pos, vel, &id, current_cell);
            
            // Do unit conversions
            for(int j = 0; j < 3; j++){
                pos[j] *= header.BoxSizeHMpc;
                vel[j] *= header.VelZSpace_to_kms;
            }
            wrap_into_box(pos);  // wrap to [0,BOXSIZE)
            
            if(N == 9){
                if(!IGNORE_PARTICLE_IDS) {
                    // pack9 PIDs come in a separate file
                    assert(fread(&id, sizeof(id), 1, pid_buffer_file)==1);
                } else {
                    id = 0;  // will be changed later to a simple counter
                }
            }

            p[i].id = id;

            assert(p[i].id >= 0 && p[i].id <= MAX_PID);

            idx = id / (ppd*ppd);
            idy = (id % (ppd*ppd)) / ppd;
            idz = (id % (ppd*ppd)) % ppd;
            
            // Only advance the pointers if we are keeping this particle
            if (IGNORE_PARTICLE_IDS ||
                (idx%ds == 0 && idy%ds == 0 && idz%ds == 0)){
                i++;
            }

            // Particle is now loaded; do one last sanity check!
            for(int j = 0; j < 3; j++){
                assert(std::isfinite(pos[j]));
                assert(std::isfinite(vel[j]));
            }
        }
    }
    
    return i;
}


extern "C" void load_particles_packN(int packN, char *filename, struct particle **p, int64_t *num_p)
{
    //printf("Loading packN particles from %s\n", filename);
    packNheader header;
    uint64_t header_size, data_size;

    header.packN = packN;
    
    // File sanity checks; also, get the header and data sizes
    packN_validate_file(packN, filename, &header_size, &data_size);
    //printf("Validated file.\n");
    
    // Populates the header struct
    packN_load_header(filename, &header);
    //printf("Loaded header.\n");
    
    // Take the header info and load parameters into global vars
    packN_parse_header_info(&header);
    //printf("Parsed header.\n");

    // Don't know how many particles to allocate (some could be cells), so use the upper limit
    // Doesn't matter if we run over a bit.  At the next realloc, it will get resized appropriately.
    // Probably want to change this to do one of the following
    // -read the particles into a new buffer
    // -do a counting pass of the particles
    // -or do one big malloc instead of lots of reallocs
    uint64_t max_new_particles = data_size / packN;
    *p = (struct particle *)check_realloc(*p, ((*num_p)+max_new_particles)*sizeof(struct particle), "Allocating particles.");
    
    // Now read the data section into the buffer
    char *pid_buffer = NULL;
    size_t pid_buffer_len = 0;
    char *buffer = packN_load_data(packN, filename, header_size, data_size, &pid_buffer, &pid_buffer_len);
    //printf("Loaded data.\n");
    FILE *buffer_file = fmemopen(buffer,data_size,"rb");
    FILE *pid_buffer_file;
    
    uint64_t new_particles;
    if(packN == 14)
        new_particles = _load_particles_packN<14>(header, buffer_file, NULL, *p + *num_p);
    else if (packN == 9){
        if(!IGNORE_PARTICLE_IDS)
            pid_buffer_file = fmemopen(pid_buffer, pid_buffer_len,"rb");
        new_particles = _load_particles_packN<9>(header, buffer_file, pid_buffer_file, *p + *num_p);
    } else {
        fprintf(stderr, "[Error] Unknown packN %d\n", packN);
        exit(1);
    }

    check_fclose(buffer_file);
    free(buffer);

    if(packN == 9 && !IGNORE_PARTICLE_IDS){
        check_fclose(pid_buffer_file);
        free(pid_buffer);
    }
    
    // Read more particles than we alloc'd
    if(new_particles > max_new_particles){
        fprintf(stderr, "[Error] Read %lu particles from %s; expected %lu at most.\n", new_particles, filename, max_new_particles);
        exit(1);
    }
    
    // Read fewer particles than we alloc'd
    //if(new_particles < max_new_particles){
    //    *p = (struct particle *)check_realloc(*p, ((*num_p)+new_particles)*sizeof(struct particle), "Shrinking particle allocation.");
    //}

    *num_p += new_particles;
    //printf("Done loading %lu new particles\n", new_particles);
}

