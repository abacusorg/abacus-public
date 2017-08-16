#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <assert.h>
#include "../config_vars.h"
#include "../config.h"
#include "../rockstar.h"
#include "../io/meta_io.h"
#include "../check_syscalls.h"

struct halo *ss_halos;

#define FAST3TREE_TYPE struct halo
#define FAST3TREE_PREFIX SUBSAMPLE
#include "../fast3tree.c"

#define GROUP_LIST ss_halos
// r is based on the main mass definition (nominally rvir)
#define RADIUS r
// The halo radii are defined in kpc/h, while the halo separations are Mpc/h
#define RADIUS_CONVERSION 1e-3
#define parent parent_id
#include "../parents.c"
#undef parent

int sort_halos_by_id(const void *a, const void *b) {
  const struct halo *c = a;
  const struct halo *d = b;
  if (c->id < d->id) return -1;
  if (c->id > d->id) return 1;
  return 0;
}

// A utility function that loads halo records from all blocks
// then calculates parent/sub relationships
// We have to load all halos simultaneously because halos in one block
// can be subs of halos in another
void find_parents_binary(int64_t snap){
    int64_t i, total_groups = 0;
    int64_t nhalo_per_file[NUM_WRITERS];
    int64_t *first_ids = NULL;
    
    first_ids = (int64_t*) malloc(sizeof(int64_t)*NUM_WRITERS);
    
    //printf("Opening files and counting halos...\n");
    // Open all tables to get the total number of halos
    for (i=0; i<NUM_WRITERS; i++) {
        struct binary_output_header header;
        load_binary_header(snap, i, &header);
        nhalo_per_file[i] = header.num_halos;
        total_groups += header.num_halos;
    }
    
    ss_halos = (struct halo*) malloc(sizeof(struct halo)*total_groups);
    
    //printf("Loading halos...\n");
    // Load all halos
    int64_t num_groups_read = 0;
    for (i=0; i<NUM_WRITERS; i++) {
        if (nhalo_per_file[i] <= 0) continue;
        load_binary_halos_subsample(snap, i, ss_halos + num_groups_read);
        first_ids[i] = ss_halos[num_groups_read].id;
        num_groups_read += nhalo_per_file[i];
    }

    //printf("Finding parents...\n");
    find_parents(total_groups);
    
    // Finally, check that all halos have been processed
    for(i=0; i < total_groups; i++){
        if(ss_halos[i].parent_id == -100){
            //ss_halos[i].parent_id = -1;
            fprintf(stderr, "[Error] Uninitialized parent_id " PRId64 " found at halo " PRId64 ".\n", ss_halos[i].parent_id, i);
            exit(1);
        }
    }
    
    // Have to re-sort the halos because find_parents scrambles them,
    // possibly across files!
    qsort(ss_halos, total_groups, sizeof(struct halo), sort_halos_by_id);
    
    //printf("Writing back halos\n");
    // Now, write back out the halos
    int64_t num_groups_written = 0;
    for (i=0; i<NUM_WRITERS; i++) {
        if (nhalo_per_file[i] <= 0) continue;
        
        // Check that we've found the correct block of halos for this file
        // We assume the halo IDs match the indices in the ss_halos array
        assert(ss_halos[first_ids[i]].id == first_ids[i]);
        
        // Write the new records
        char fn[1024];
        FILE *of;
        get_output_filename(fn, 1024, snap, i, "bin");
        of = fopen(fn, "r+b");
        check_fseeko(of, sizeof(struct binary_output_header), SEEK_SET);
        check_fwrite(ss_halos + first_ids[i], sizeof(struct halo), nhalo_per_file[i], of);
        
        num_groups_written += nhalo_per_file[i];
        
        // Close the file
        fclose(of);
    }
    printf("\n");
    assert(num_groups_written == total_groups);
    
    free(first_ids);
    free(ss_halos);
}

int main(int argc, char **argv)
{
  int64_t i, snap=0, did_config = 0;
  if (argc < 3) {
    printf("Usage: %s -c config\n", argv[0]);
    exit(1);
  }

  for (i=1; i<argc-1; i++) {
    if (!strcmp("-c", argv[i])) { do_config(argv[i+1]); i++; did_config=1; }
    if (!strcmp("-s", argv[i])) { snap = atoi(argv[i+1]); i++; }
  }
  if (!did_config){
      fprintf(stderr, "[Error] No config!\n");
      exit(1);
  }
      
  if (strlen(SNAPSHOT_NAMES)) 
    read_input_names(SNAPSHOT_NAMES, &snapnames, &NUM_SNAPS);
  if (strlen(BLOCK_NAMES))
    read_input_names(BLOCK_NAMES, &blocknames, &NUM_BLOCKS);

  printf("Identifying parent/sub halos...\n");
  find_parents_binary(snap);
  return 0;
}
