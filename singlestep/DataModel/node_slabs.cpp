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
        // In the 2D decomposition, NodeSlabs only relates to the x-decomposition.
        // So one should use MPI_rank_x and MPI_size_x in this function

        int neighbor = (MPI_rank_x+1)%MPI_size_x;
        int value, last_slab;
		
		int *last_slabs = new int[MPI_size_x];
			

        //NAM DE TODO have convolution look at MultipoleDirectory for node x domain. Check 0th step --> what comes first, singlestep or convolve? 
        FILE *fp = fopen((P.ReadStateDirectory / "nodeslabs").c_str(),"r");
        // assertf(fp!=NULL, "Couldn't find nodeslabs file {}\n", fname);
        if (fp==NULL) {

			int offset = 17;
            // We couldn't find a file, so let's make up something. +3 is to test periodic wrapping! 
            first_slab_on_node = (int)(floor((float)P.cpd*MPI_rank_x/MPI_size_x) + offset)%P.cpd;
            last_slab = (int)(floor((float)P.cpd*(MPI_rank_x+1)/MPI_size_x) + offset)%P.cpd;

			if (get_all_nodes){
				for (int j=0; j<MPI_size_x; j++) {
					first_slabs_all[j] = (int)(floor((float)P.cpd*j/MPI_size_x) + offset)%P.cpd;
					total_slabs_all[j] = floor((float)P.cpd*(j+1)/MPI_size_x) - floor((float)P.cpd*j/MPI_size_x);			
																					
				}				
			}
			
        } else {		
			for (int j=0; j<MPI_size_x; j++) {

				
                int nread = fscanf(fp, "%d", &value);
                assertf(nread==1, "Couldn't read entry {:d} from NodeSlabs file\n", j);
                if (j==MPI_rank_x) first_slab_on_node = value;
                if (j==neighbor) last_slab = value;
				
				if (get_all_nodes) {

					first_slabs_all[j] = value;

					if (j>0) last_slabs[(j-1)%MPI_size_x] = value;
					else if (j==0) last_slabs[MPI_size_x-1] = value; 						
				}
						
            }
			
			
			if (get_all_nodes) {
				for (int j=0; j<MPI_size_x; j++) {
					total_slabs_all[j] = last_slabs[j] - first_slabs_all[j]  ; 
			        if (total_slabs_all[j]<0) total_slabs_all[j] += P.cpd;
			        STDLOG(1,"Read NodeSlab file: node {:d} will do {:d} slabs starting with {:d}, last {:d}\n", j, total_slabs_all[j], first_slabs_all[j], last_slabs[j]);
				
				}
			}
            fclose(fp);
        }
		
		total_slabs_on_node = last_slab - first_slab_on_node;

		
        if (total_slabs_on_node<0) total_slabs_on_node += P.cpd;
        STDLOG(1,"Read NodeSlab file: will do {:d} slabs from [{:d},{:d})\n",
            total_slabs_on_node, first_slab_on_node, last_slab);
		if (get_all_nodes) {	
			assert(total_slabs_all[MPI_rank_x] == total_slabs_on_node) ;
			assert(first_slabs_all[MPI_rank_x] == first_slab_on_node) ;
		}	
			
		delete[] last_slabs; 	
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
        if (MPI_size_x==1) return;   // We don't want to use this file if we're serial
        if (MPI_rank_z > 0) return;  // Only one row needs to compute nodeslabs

        int first[MPI_size_x];
        for (int j=0; j<MPI_size_x; j++) first[j]=0;
        first[MPI_rank_x] = first_slab_finished;
        // MPI: Send first_slab_finished to fill in this element in the vector
        // We just do this as a boring summation.
        MPI_Reduce(MPI_rank_x != 0 ? first : MPI_IN_PLACE, first, MPI_size_x, MPI_INT, MPI_SUM, 0, comm_1d_x);

        if (MPI_rank_x==0) {
			
            //NAM DE TODO consider putting another copy in multipole directory so convolve can look at it. 
            fs::path fname = P.WriteStateDirectory / "nodeslabs";
            fs::path mname = P.MultipoleDirectory / "nodeslabs";
			
            FILE *fp, *fm; 
            fp = fopen(fname.c_str(),"w"); fm = fopen(mname.c_str(), "w");
            assertf(fp!=NULL, "Couldn't create nodeslabs file {}\n", fname);
            assertf(fm!=NULL, "Couldn't create nodeslabs file {}\n", mname);
			
            for (int j=0; j<MPI_size_x; j++) {
                fmt::print(fp, "{:d}\n", first[j]);
                fmt::print(fm, "{:d}\n", first[j]);
				
            }
			
			
            STDLOG(1, "Wrote the NodeSlab file to {} and {}\n", fname, mname);
			
			
        }
    #endif
}
