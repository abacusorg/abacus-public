/** This code describes a file that will direct the nodes as to 
which slabs to work on.  The file is simply NNode numbers, which 
list the first slab of each node.  A node needs that number as well
as the next one (cyclicly wrapped).

When finishing the timestep, we need to gather up the first finished
slabs and construct a new file.
**/

/// Read the global file to set two global variables for this node:
/// first_slab_on_node and total_slab_on_node
void ReadNodeSlabs() {
    #ifndef PARALLEL
        first_slab_on_node = 0; total_slabs_on_node = P.cpd;
        first_slab_finished = -1;   // Just a silly value
        return;
    #else
        int neighbor = (MPI_rank+1)%MPI_size;
        char fname[1024];
        int value, last_slab;
        // TODO: This needs to be the Global Read State
        sprintf(fname, "%s/nodeslabs", P.ReadStateDirectory);
        FILE *fp;
        fp = fopen(fname,"r");
        // assertf(fp!=NULL, "Couldn't find nodeslabs file %s\n", fname);
        if (fp==NULL) {
            // We couldn't find a file, so let's make up something
            first_slab_on_node = floor((float)P.cpd*MPI_rank/MPI_size);
            last_slab = floor((float)P.cpd*(MPI_rank+1)/MPI_size);
        } else {
            for (int j=0; j<MPI_size; j++) {
                int nread = fscanf(fp, "%d", &value);
                assertf(nread==1, "Couldn't read entry %j from NodeSlabs file\n", j);
                if (j==MPI_rank) first_slab_on_node = value;
                if (j==neighbor) last_slab = value;
            }
            fclose(fp);
        }
        total_slabs_on_node = last_slab - first_slab_on_node;
        if (total_slabs_on_node<0) total_slabs_on_node += P.cpd;
        STDLOG(1,"Read NodeSlab file: will do %d slabs from [%d,%d)\n",
            total_slabs_on_node, first_slab_on_node, last_slab);
    #endif
    return;
}

/// Gather the information from all of the nodes to node 0 and have
/// it write the global NodeSlabs file.
void WriteNodeSlabs() {
    #ifndef PARALLEL
        // Note that first_slab_finished is only set in PARALLEL
        return;
    #else
        if (MPI_size==1) return;   // We don't want to use this file if we're serial
        int first[MPI_size];
        for (int j=0; j<MPI_size; j++) first[j]=0;
        first[MPI_rank] = first_slab_finished;
        // MPI: Send first_slab_finished to fill in this element in the vector
        // We just do this as a boring summation.
        MPI_REDUCE_TO_ZERO(first, MPI_size, MPI_INT, MPI_SUM);

        if (MPI_rank==0) {
            char fname[1024];
            // TODO: This needs to be the Global Write State
            sprintf(fname, "%s/nodeslabs", P.WriteStateDirectory);
            FILE *fp;
            fp = fopen(fname,"w");
            assertf(fp!=NULL, "Couldn't create nodeslabs file %s\n", fname);
            for (int j=0; j<MPI_size; j++) {
                fprintf(fp, "%d\n", first[j]);
            }
            fclose(fp);
            STDLOG(1, "Wrote the NodeSlab file to %s\n", fname);
        }
    #endif
}
