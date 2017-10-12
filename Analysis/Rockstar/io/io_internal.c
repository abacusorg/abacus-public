#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "meta_io.h"
#include "../config_vars.h"
#include "../rockstar.h"
#include "../groupies.h"
#include "../check_syscalls.h"
#include "../version.h"
#include "../halo.h"

void *output_buffer = NULL;
int64_t buffered = 0;

void fill_binary_header(struct binary_output_header *bh,
			int64_t snap, int64_t chunk) {
  bh->magic = ROCKSTAR_MAGIC;
  bh->snap = snap;
  bh->chunk = chunk;
  bh->scale = SCALE_NOW;
  bh->Om = Om;
  bh->Ol = Ol;
  bh->h0 = h0;
  bh->box_size = BOX_SIZE;
  bh->particle_mass = PARTICLE_MASS;  
  snprintf(bh->rockstar_version, VERSION_MAX_SIZE, "%s", ROCKSTAR_VERSION);
  bh->format_revision = HALO_FORMAT_REVISION;
}


void read_binary_header_config(struct binary_output_header *bh) {
  SCALE_NOW = bh->scale;
  Om = bh->Om;
  Ol = bh->Ol;
  h0 = bh->h0;
  BOX_SIZE = bh->box_size;
  PARTICLE_MASS = bh->particle_mass;
}

#include "MurmurHash3.cpp"
#define SUBSAMPLE_SEED 0xBEEF
inline int is_subsample_particle(const int64_t pid, const double subsample_frac){
    uint32_t hash = 0;
    MurmurHash3_x86_32(&pid, sizeof(pid), SUBSAMPLE_SEED, &hash);
    double prob = (double) hash / UINT32_MAX;
    return prob < subsample_frac;
}

void output_particles_internal(int64_t snap, int64_t chunk, double fraction) {
  char buffer[1024];
  struct binary_output_header bh;
  FILE *output;
  int64_t i, j;
  
  memset(&bh, 0, sizeof(struct binary_output_header));
  fill_binary_header(&bh, snap, chunk);
  bh.num_particles = num_p;
  bh.particle_type = PARTICLE_TYPE_FULL;
  get_output_filename(buffer, 1024, snap, chunk, "rbin");
  if (fraction>0 && fraction < 1)
    get_output_filename(buffer, 1024, snap, chunk, "rbin_wl");

  output = check_fopen(buffer, "wb");
  check_fwrite(&bh, sizeof(struct binary_output_header), 1, output);
  if (fraction==0 || fraction == 1) {
    check_fwrite(p, sizeof(struct particle), num_p, output);
  } else {
    for (i=0,j=0; i<num_p; i++) {
      if (is_subsample_particle(p[i].id, fraction)) {
        check_fwrite(p+i, sizeof(struct particle), 1, output);
        j++;
      }
    }
    bh.num_particles = j;
    bh.particle_mass /= fraction;
    rewind(output);
    check_fwrite(&bh, sizeof(struct binary_output_header), 1, output);
  }
  fclose(output);
}


void load_particles_internal(char *filename, struct particle **part, int64_t *num_part) {
  FILE *input;
  struct binary_output_header bh;
  
  input = check_fopen(filename, "rb");
  check_fread(&bh, sizeof(struct binary_output_header), 1, input);
  assert(bh.magic == ROCKSTAR_MAGIC);
  assert(bh.particle_type == PARTICLE_TYPE_FULL);
  assert(bh.num_particles >= 0);
  read_binary_header_config(&bh);
  check_realloc_s(*part, sizeof(struct particle), bh.num_particles);
  check_fread(*part, sizeof(struct particle), bh.num_particles, input);
  *num_part = bh.num_particles;
  fclose(input);
}


inline void _clear_buffer(FILE *output) {
  if(buffered == 0)
      return;
  check_fwrite(output_buffer, buffered, 1, output);
  buffered = 0;
}

inline void _append_to_buffer(void *src, int64_t size, FILE *output) {
  if (buffered + size > OUTPUT_BUFFER_SIZE) _clear_buffer(output);
  if (size > OUTPUT_BUFFER_SIZE) { check_fwrite(src, size, 1, output);
      return;
  }
  memcpy(((char *)output_buffer)+buffered, src, size);
  buffered+=size;
}

// Used to update the binary files with the SO (BGC2) masses
void output_binary(int64_t id_offset, int64_t snap, int64_t chunk, float *bounds, int64_t output_particles) {
    assert(output_particles == 0);
    
  float max[3]={0}, min[3]={0};
  char buffer[1024];
  int64_t i,j, id=0;
  FILE *output;
  struct binary_output_header bheader;

  if (!output_buffer) check_realloc_s(output_buffer, 1, OUTPUT_BUFFER_SIZE);

  get_output_filename(buffer, 1024, snap, chunk, "bin");
  output = check_fopen(buffer, "r+b");

  load_binary_header(snap, chunk, &bheader);
  _append_to_buffer(&bheader, sizeof(struct binary_output_header), output);

  //Output Halos
  if (num_halos) { 
    memcpy(min, halos[0].pos, sizeof(float)*3);
    memcpy(max, halos[0].pos, sizeof(float)*3);
  }
  for (i=0; i<num_halos; i++) {
    if (!_should_print(halos+i, bounds)) continue;
    for (j=0; j<3; j++) {
      if (min[j]>halos[i].pos[j]) min[j] = halos[i].pos[j];
      if (max[j]<halos[i].pos[j]) max[j] = halos[i].pos[j];
    }
    _append_to_buffer(halos + i, sizeof(struct halo), output);
    id++;
  }

  _clear_buffer(output); 
  fclose(output);
}

void output_binary_subsample(int64_t id_offset, int64_t snap, int64_t chunk, float *bounds, int64_t output_particles) {
    assert(output_particles == 1);
    
  float max[3]={0}, min[3]={0};
  struct halo tmp;
  char fnbuffer[1024];
  int64_t i,j, id=0;
  FILE *output, *ss_output;
  struct binary_output_header bheader;

  if (!output_buffer) check_realloc_s(output_buffer, 1, OUTPUT_BUFFER_SIZE);

  get_output_filename(fnbuffer, 1024, snap, chunk, "bin");
  output = check_fopen(fnbuffer, "wb");

  get_particle_output_filename(fnbuffer, 1024, snap, chunk, "bin");
  ss_output = check_fopen(fnbuffer, "wb");

  memset(&bheader, 0, sizeof(struct binary_output_header));
  _append_to_buffer(&bheader, sizeof(struct binary_output_header), output);
  _clear_buffer(output);

  //Output Halos
  if (num_halos) {
    memcpy(min, halos[0].pos, sizeof(float)*3);
    memcpy(max, halos[0].pos, sizeof(float)*3);
  }
  for (i=0; i<num_halos; i++) {
    if (!_should_print(halos+i, bounds)) continue;
    for (j=0; j<3; j++) {
      if (min[j]>halos[i].pos[j]) min[j] = halos[i].pos[j];
      if (max[j]<halos[i].pos[j]) max[j] = halos[i].pos[j];
    }
    halos[i].subsamp_start = bheader.num_particles;
    halos[i].parent_id = -100;  // signal not yet valid
      
    // Output subsample particles for this halo
    // Modifies the subsamp fields of the original halo structs so they are preserved when overwritten by the SO outputs
    halos[i].subsamp_len = 0;
    for (j = halos[i].p_start; j < halos[i].p_start + halos[i].num_p; j++){
        if(is_subsample_particle(p[j].id, SUBSAMPLE_FRAC)){
            _append_to_buffer(p + j, sizeof(p[0]), ss_output);
            halos[i].subsamp_len++;
        }
    }
    bheader.num_particles += halos[i].subsamp_len;
    _clear_buffer(ss_output);
     
    //tmp = halos[i];
    //tmp.id = id+id_offset;
    halos[i].id = id+id_offset;
    
    _append_to_buffer(halos + i, sizeof(struct halo), output);
    bheader.num_halos++;
    _clear_buffer(output);
    
    id++;
  }
  
  fclose(ss_output);

  //Output header
  fill_binary_header(&bheader, snap, chunk);
  bheader.particle_type = PARTICLE_TYPE_FULL;
  if (bounds) memcpy(bheader.bounds, bounds, sizeof(float)*6);
  else { memcpy(bheader.bounds, min, sizeof(float)*3);
    memcpy(&(bheader.bounds[3]), max, sizeof(float)*3);
  }
  rewind(output);
  check_fwrite(&bheader, sizeof(struct binary_output_header), 1, output);
  fclose(output);
}

void load_binary_header(int64_t snap, int64_t chunk, 
			struct binary_output_header *bheader) {
  char buffer[1024];
  FILE *input;
  get_output_filename(buffer, 1024, snap, chunk, "bin");
  input = check_fopen(buffer, "rb");
  check_fread(bheader, sizeof(struct binary_output_header), 1, input);
  assert(bheader->magic == ROCKSTAR_MAGIC);
  fclose(input);
}


void load_binary_halos(int64_t snap, int64_t chunk, 
      struct binary_output_header *bheader, struct halo **halos,
		       int64_t **part_ids, int64_t coalesced)
{
  char buffer[1024];
  FILE *input;
  int64_t i, j;

  if (!coalesced) get_output_filename(buffer, 1024, snap, chunk, "bin");
  else get_output_filename(buffer, 1024, snap, chunk, "coalesced.bin");
  input = check_fopen(buffer, "rb");
  check_fread(bheader, sizeof(struct binary_output_header), 1, input);
  assert(bheader->magic == ROCKSTAR_MAGIC);
  assert(bheader->num_halos >= 0);
  assert(bheader->num_particles >= 0);
  check_realloc_s(*halos, sizeof(struct halo), bheader->num_halos);
  i = check_fread(*halos, sizeof(struct halo), bheader->num_halos, input);
  
  if (i!=bheader->num_halos) {
    fprintf(stderr, "[Error] Truncated input file %s!\n", buffer);
    exit(1);
  }
    
  if(part_ids != NULL){
      check_realloc_s(*part_ids, sizeof(int64_t), bheader->num_particles);
      j = check_fread(*part_ids, sizeof(int64_t), bheader->num_particles, input);
      if (j!=bheader->num_particles) {
        fprintf(stderr, "[Error] Truncated input file %s!\n", buffer);
        exit(1);
      }
  }
    
  fclose(input);

  read_binary_header_config(bheader);
}

void load_binary_halos_subsample(int64_t snap, int64_t chunk, struct halo *halos)
{
  char fnbuffer[1024];
  FILE *input;
  int64_t i;
  struct binary_output_header bheader;

  get_output_filename(fnbuffer, 1024, snap, chunk, "bin");
  input = check_fopen(fnbuffer, "rb");
  check_fread(&bheader, sizeof(struct binary_output_header), 1, input);
  assert(bheader.magic == ROCKSTAR_MAGIC);
  assert(bheader.num_halos >= 0);
  assert(bheader.num_particles >= 0);
  i = check_fread(halos, sizeof(struct halo), bheader.num_halos, input);
  
  if (i!=bheader.num_halos) {
    fprintf(stderr, "[Error] Truncated input file %s!\n", fnbuffer);
    exit(1);
  }
    
  fclose(input);

  //read_binary_header_config(bheader);
}
