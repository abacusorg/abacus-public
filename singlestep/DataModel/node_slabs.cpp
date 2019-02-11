/** This code describes a file that will direct the nodes as to 
which slabs to work on.  The file is simply NNode numbers, which 
list the first slab of each node.  A node needs that number as well
as the next one (cyclicly wrapped).

When finishing the timestep, we need to gather up the first finished
slabs and construct a new file.
**/

/// Read the global file to set two global variables for this node:
/// first_slab_on_node and total_slab_on_node
void ReadNodeSlabs(int get_all_nodes = 0, int * first_slabs_all = NULL, int * total_slabs_all = NULL) {
    #ifndef PARALLEL
        first_slab_on_node = 0; total_slabs_on_node = P.cpd;
        first_slab_finished = -1;   // Just a silly value
        return;
    #else
        int neighbor = (MPI_rank+1)%MPI_size;
        char fname[1024];
        int value, last_slab;
		
		int *last_slabs = new int[MPI_size];
			

        sprintf(fname, "%s/nodeslabs", P.ReadStateDirectory); //NAM DE TODO have convolution look at MultipoleDirectory for node x domain. Check 0th step --> what comes first, singlestep or convolve? 
        FILE *fp;
        fp = fopen(fname,"r");
        // assertf(fp!=NULL, "Couldn't find nodeslabs file %s\n", fname);
        if (fp==NULL) {
			int offset = 17;
            // We couldn't find a file, so let's make up something. +3 is to test periodic wrapping! 
            first_slab_on_node = (int)(floor((float)P.cpd*MPI_rank/MPI_size) + offset)%P.cpd;
            last_slab = (int)(floor((float)P.cpd*(MPI_rank+1)/MPI_size) + offset)%P.cpd;

			if (get_all_nodes){
				for (int j=0; j<MPI_size; j++) {
					first_slabs_all[j] = (int)(floor((float)P.cpd*j/MPI_size) + offset)%P.cpd;
					total_slabs_all[j] = floor((float)P.cpd*(j+1)/MPI_size) - floor((float)P.cpd*j/MPI_size);			
																					
				}				
			}
			
        } else {		
			for (int j=0; j<MPI_size; j++) {

				
                int nread = fscanf(fp, "%d", &value);
                assertf(nread==1, "Couldn't read entry %j from NodeSlabs file\n", j);
                if (j==MPI_rank) first_slab_on_node = value;
                if (j==neighbor) last_slab = value;
				
				if (get_all_nodes) {

					first_slabs_all[j] = value;

					if (j>0) last_slabs[(j-1)%MPI_size] = value;
					else if (j==0) last_slabs[MPI_size-1] = value; 						
				}
						
            }
			
			
			if (get_all_nodes) {
				for (int j=0; j<MPI_size; j++) {
					total_slabs_all[j] = last_slabs[j] - first_slabs_all[j]  ; 
			        if (total_slabs_all[j]<0) total_slabs_all[j] += P.cpd;
			        STDLOG(1,"Read NodeSlab file: node %d will do %d slabs starting with %d, last %d\n", j, total_slabs_all[j], first_slabs_all[j], last_slabs[j]);
				
				}
			}
            fclose(fp);
        }
		
		total_slabs_on_node = last_slab - first_slab_on_node;

		
        if (total_slabs_on_node<0) total_slabs_on_node += P.cpd;
        STDLOG(1,"Read NodeSlab file: will do %d slabs from [%d,%d)\n",
            total_slabs_on_node, first_slab_on_node, last_slab);
		if (get_all_nodes) {	
			assert(total_slabs_all[MPI_rank] == total_slabs_on_node) ;
			assert(first_slabs_all[MPI_rank] == first_slab_on_node) ;
		}	
			
		delete last_slabs; 	
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
			char mname[1024];
			
            // TODO: This needs to be the Global Write State
            sprintf(fname, "%s/nodeslabs", P.WriteStateDirectory); //NAM DE TODO consider putting another copy in multipole directory so convolve can look at it. 
			sprintf(mname, "%s/nodeslabs", P.MultipoleDirectory);
			
            FILE *fp, *fm; 
            fp = fopen(fname,"w"); fm = fopen(mname, "w");
            assertf(fp!=NULL, "Couldn't create nodeslabs file %s\n", fname);
            assertf(fm!=NULL, "Couldn't create nodeslabs file %s\n", mname);
			
            for (int j=0; j<MPI_size; j++) {
                fprintf(fp, "%d\n", first[j]);
                fprintf(fm, "%d\n", first[j]);
				
            }
			
			
            STDLOG(1, "Wrote the NodeSlab file to %s and %s\n", fname, mname);
			
			
        }
    #endif
}
