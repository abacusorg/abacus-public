// This is the top-level code for group-finding.
// We include this file in the program.


// TODO: Are the SlabAccum maxsize set sensibly?

#include "fof_sublist.cpp"
	// Code to do FOF on a set (e.g., a cell)

#include "slab_accum.cpp"
	// Code to establish templated slab-based storage of flexible size 
	// that is cell indexed and multi-threaded by pencil

#include "grouplink.cpp"
	// Code to store the list of links between CellGroups.
	// This is GCC->GLL

#include "cellgroup.cpp"
	// Code to define and create CellGroups.



class GroupFindingControl {
    // This is the control class for all group finding.
  public:
    int cpd;		
    FOFloat linking_length;    // In code units
    FOFloat boundary;	// The distance from the origin that indicates that something
    		// is within linking_length of the edge
    int GroupRadius;    // Cell radius for group finding
    Float invcpd;
    uint64 np; // the number of particles in the simulation; used to compute buffers

    // Parameters for finding subgroups
    FOFloat linking_length_level1;
    FOFloat linking_length_level2;
    int minhalosize;	// The minimum size for a level 1 halo to be outputted

    uint64 pPtot, fPtot, fGtot, CGtot, GGtot, Ltot;

    SlabAccum<CellGroup> *cellgroups;   // Allocated [0,cpd), one for each slab
    int *cellgroups_status;     // Allocated [0,cpd) for each slab. 
    	// 0 means CellGroups not found.
    	// 1 means CellGroups found but slab not closed.
    	// 2 means slab closed (so CellGroups may be deleted).
    GroupLinkList *GLL;

    STimer CellGroupTime, CreateFaceTime, FindLinkTime, 
    	SortIndexLinks, FindGlobalGroupTime, GatherGroups, ScatterGroups;

    GroupFindingControl(FOFloat _linking_length, int _cpd, FOFloat _invcpd, int _GroupRadius, uint64 _np) {
	cpd = _cpd; 
	linking_length = _linking_length;
	invcpd = _invcpd;
	boundary = (invcpd/2.0-linking_length);
	np = _np;
	GroupRadius = _GroupRadius;
	cellgroups = new SlabAccum<CellGroup>[cpd];
	cellgroups_status = new int[cpd];
	for (int j=0;j<cpd;j++) cellgroups_status[j] = 0;
	GLL = new GroupLinkList(cpd, np/cpd*linking_length/_invcpd*3*15);    
	// This is a MultiAppendList, so the buffer cannot grow. 
	// It will be storing links between cells.  Most particles do 
	// not generate such a link; on average 3*linking_length/cellsize do. 
	// But we will have up to 10 slabs worth in here.
	// So we need to be a bit generous.
        pPtot = fPtot = fGtot = CGtot = GGtot = Ltot = 0;
	return;
    }

    ~GroupFindingControl() { 
        delete[] cellgroups;
        delete[] cellgroups_status;
	delete GLL;
    }
    int WrapSlab(int s) {
	while (s<0) s+=cpd;
	while (s>=cpd) s-=cpd;
	return s;
    }
    integer3 WrapCell(int i, int j, int k) {
        return integer3(WrapSlab(i), WrapSlab(j), WrapSlab(k));
    }
    
    void ConstructCellGroups(int slab);
    void DestroyCellGroups(int slab);


    void report() {
	 // TODO: Eventually these should be STDLOG
	 printf("Found %lld cell groups (including boundary singlets)\n", CGtot);
	 printf("Used %lld pseudoParticles, %lld faceParticles, %lld faceGroups\n",
	     pPtot, fPtot, fGtot);
	 printf("Found %lld links between groups.\n", Ltot);
	 printf("Found %lld global groups\n", GGtot);

	 printf("Timings: \n");
	 printf("Finding Cell Groups:    %f sec\n",
			CellGroupTime.Elapsed());
	 printf("Creating Faces:         %f sec\n",
			CreateFaceTime.Elapsed());
	 printf("Finding Group Links:    %f sec\n",
			FindLinkTime.Elapsed());
	 printf("Sort & Index Links:     %f sec\n",
			SortIndexLinks.Elapsed());
	 printf("Find Global Groups:     %f sec\n",
			FindGlobalGroupTime.Elapsed());
	 printf("Gather Group Particles: %f sec\n",
			GatherGroups.Elapsed());
	 printf("Scatter Group Particles: %f sec\n",
			ScatterGroups.Elapsed());
    }
};

// We keep the code for Constructing CellGroups in here, because
// it is where the cellgroups[] are created and destroyed.

void GroupFindingControl::ConstructCellGroups(int slab) {
    // Construct the Cell Groups for this slab
    CellGroupTime.Start();
    slab = WrapSlab(slab);
    cellgroups[slab].setup(cpd, 1024);     // TODO: Pick a better maxsize!
    FOFcell doFOF[omp_get_max_threads()];
    #pragma omp parallel for schedule(static)
    for (int g=0; g<omp_get_num_threads(); g++) 
    	doFOF[g].setup(linking_length, boundary);

    #pragma omp parallel for schedule(dynamic,1)
    for (int j=0; j<cpd; j++) {
	float *aligned;
	int ret = posix_memalign((void **)&(aligned), 64, 8*sizeof(float)); assert(ret==0);
	int g = omp_get_thread_num();
        PencilAccum<CellGroup> *cg = cellgroups[slab].StartPencil(j);
        for (int k=0; k<cpd; k++) {
	    // Find CellGroups in (slab,j,k).  Append results to cg.
	    Cell c = PP->GetCell(slab, j, k);
	    doFOF[g].findgroups(c.pos, c.vel, c.aux, c.count());
	    // printf("Cell %d %d %d: Found %d cell groups with %d particles, plus %d boundary singlets\n", slab,j,k,doFOF[g].ngroups, doFOF[g].nmultiplets, doFOF[g].nsinglet_boundary-doFOF[g].nmultiplets);
	    for (int gr=0; gr<doFOF[g].ngroups; gr++) {
		CellGroup tmp(doFOF[g].groups[gr], boundary, aligned);
		// printf("Group %d: %d %d, %f %f %f to %f %f %f, %f %02x\n",
		// 	gr, tmp.start, tmp.size(), 
		// 	aligned[0], aligned[1], aligned[2],
		// 	aligned[4], aligned[5], aligned[6], boundary, tmp.n >> 25);
		cg->append(tmp);
	    }
	    // Also need to look at the singlets!
	    for (int p=doFOF[g].nmultiplets; p<doFOF[g].nsinglet_boundary; p++) {
		CellGroup tmp(p, c.pos[p], boundary);
		// printf("Singl %d: %f %f %f %lld %02x\n",
		   // p, c.pos[p].x, c.pos[p].y, c.pos[p].z, c.aux[p].val, tmp.n >>25);
	        cg->append(tmp);
	    }
	    // for (int p=0; p<doFOF[g].np; p++) {
	        // if (c.aux[p].val == 1516) {
		    // printf("Singl %d: %f %f %f %lld     %d\n",
		       // p, c.pos[p].x, c.pos[p].y, c.pos[p].z, c.aux[p].val,
		       // doFOF[g].nsinglet_boundary);
		// }
	    // }
	    cg->FinishCell();
	}
	cg->FinishPencil();
	free(aligned);
    }
    for (int g=0; g<omp_get_num_threads(); g++) 
    	doFOF[g].destroy();
    uint64 tot = cellgroups[slab].get_slab_size();
    CGtot += tot;
    cellgroups_status[slab] = 1;
    STDLOG(1,"Found %lld cell groups in slab %d\n", tot, slab);
    CellGroupTime.Stop();
    return;
}

void GroupFindingControl::DestroyCellGroups(int slab) {
    slab = WrapSlab(slab);
    cellgroups[slab].destroy();
    cellgroups_status[slab] = 2;
    return;
}

GroupFindingControl *GFC;	
	// We have one global instance of this, initialized in proepi


#include "findgrouplinks.cpp"
	// Code to search between pairs of cells and find the linked groups,
	// which get added to GLL.

#include "globalgroup.cpp"
	// Code to traverse the links and find the GlobalGroups as 
	// sets of CellGroups (stored by their LinkIDs).

#include "halostat.cpp"
	// Code to compute L1 halo properties

// ===================== Global Groups: Find & Process ===============


class GlobalGroupSlab {
  public:
    SlabAccum<GlobalGroup> globalgroups;
    SlabAccum<LinkID> globalgrouplist;
    
    int slab;    // This slab number
    posstruct *pos;  // Big vectors of all of the pos/vel/aux for these global groups
    velstruct *vel;
    auxstruct *aux;
    uint64 np;

    GlobalGroupSlab(int _slab) {
        pos = NULL; vel = NULL; aux = NULL; np = 0;
	slab = _slab;
    }
    void destroy() {
	// TODO: If we use arenas, change this
        if (pos!=NULL) free(pos); pos = NULL;
        if (vel!=NULL) free(vel); vel = NULL;
        if (aux!=NULL) free(aux); aux = NULL;
	np = 0;
    }
    ~GlobalGroupSlab() { destroy(); }
    void setup(uint64 _np) {
	destroy();
	np = _np;
	// TODO: May eventually prefer these to be arenas.
	int ret;
        ret = posix_memalign((void **)&pos, 4096, sizeof(posstruct)*np); assert(ret==0);
        ret = posix_memalign((void **)&vel, 4096, sizeof(velstruct)*np); assert(ret==0);
        ret = posix_memalign((void **)&aux, 4096, sizeof(auxstruct)*np); assert(ret==0);
    }

    void GatherGlobalGroups() {
	// Copy all of the particles into a single list
	GFC->GatherGroups.Start();
	int diam = 2*GFC->GroupRadius+1;
	assert(GFC->cpd>=4*GFC->GroupRadius+1);
	// TODO: This registers the periodic wrap using the cells.
	// However, this will misbehave if CPD is smaller than the group diameter,
	// because the LinkIDs have already been wrapped.
	uint64 total_particles = 0;
	uint64 *pstart = new uint64[GFC->cpd];
	for (int j=0; j<GFC->cpd; j++) {
	    pstart[j] = total_particles;    // So we have the starting indices
	    uint64 this_pencil = 0;
	    // TODO: check this syntax
	    #pragma omp parallel for schedule(static) reduction(+:this_pencil)
	    for (int k=0; k<GFC->cpd; k++)
		for (int n=0; n<globalgroups[j][k].size(); n++) 
		    this_pencil += globalgroups[j][k][n].np;
	    total_particles += this_pencil;
	}
	// Now we can allocate these buffers
	setup(total_particles);

	// Now copy the particles into these structures
	#pragma omp parallel for schedule(static)
	for (int j=0; j<GFC->cpd; j++)
	    for (int k=0; k<GFC->cpd; k++)
		for (int n=0; n<globalgroups[j][k].size(); n++) {
		    // Process globalgroups[j][k][n]
		    // Compute where we'll put the particles, and update this starting point
		    globalgroups[j][k].ptr(n)->start += pstart[j];
		    uint64 start = globalgroups[j][k][n].start;

		    LinkID *cglink = globalgrouplist.pencils[j].data
		    			+globalgroups[j][k][n].cellgroupstart;
		    	// This is where we'll find the CG LinkIDs for this GG
		    integer3 firstcell(slab,j,k);
		    for (int c=0; c<globalgroups[j][k][n].ncellgroups; c++, cglink++) {
		        // Loop over CellGroups
			integer3 cellijk = cglink->cell();
			CellGroup *cg = LinkToCellGroup(*cglink);
			// printf("%d %d %d %d in %d %d %d n=%d\n", 
			    // j, k, n, c, cellijk.x, cellijk.y, cellijk.z, cg->size());
			// CellGroup cg = GFC->cellgroups[cellijk.x][cellijk.y][cellijk.z][cglink->cellgroup()];
			Cell cell = PP->GetCell(cellijk);
			// Copy the particles, including changing positions to 
			// the cell-centered coord of the first cell.  
			// Note periodic wrap.
			memcpy(vel+start, cell.vel+cg->start, sizeof(velstruct)*cg->size());
			memcpy(aux+start, cell.aux+cg->start, sizeof(auxstruct)*cg->size());
			cellijk -= firstcell;
			if (cellijk.x> diam) cellijk.x-=GFC->cpd;
			if (cellijk.x<-diam) cellijk.x+=GFC->cpd;
			if (cellijk.y> diam) cellijk.y-=GFC->cpd;
			if (cellijk.y<-diam) cellijk.y+=GFC->cpd;
			if (cellijk.z> diam) cellijk.z-=GFC->cpd;
			if (cellijk.z<-diam) cellijk.z+=GFC->cpd;
			posstruct offset = GFC->invcpd*(cellijk);
			// printf("Using offset %f %f %f\n", offset.x, offset.y, offset.z);
			for (int p=0; p<cg->size(); p++) pos[start+p] = offset+cell.pos[cg->start+p];
			start += cg->size();
		    } // End loop over cellgroups in this global group
		} // End loop over globalgroups in a cell
	// End loop over cells
	delete[] pstart;
	GFC->GatherGroups.Stop();
	return;
    }

    void ScatterGlobalGroups() {
        // Write the information from pos,vel,aux back into the original Slabs
	int diam = 2*GFC->GroupRadius+1;
	GFC->ScatterGroups.Start();

	#pragma omp parallel for schedule(static)
	for (int j=0; j<GFC->cpd; j++)
	    for (int k=0; k<GFC->cpd; k++)
		for (int n=0; n<globalgroups[j][k].size(); n++) {
		    // Process globalgroups[j][k][n]
		    // Recall where the particles start
		    uint64 start = globalgroups[j][k][n].start;

		    LinkID *cglink = globalgrouplist.pencils[j].data
		    			+globalgroups[j][k][n].cellgroupstart;
		    	// This is where we'll find the CG LinkIDs for this GG
		    integer3 firstcell(slab,j,k);
		    for (int c=0; c<globalgroups[j][k][n].ncellgroups; c++, cglink++) {
		        // Loop over CellGroups
			integer3 cellijk = cglink->cell();
			CellGroup *cg = LinkToCellGroup(*cglink);
			Cell cell = PP->GetCell(cellijk);
			// Copy the particles, including changing positions to 
			// the cell-centered coord of the first cell.  
			// Note periodic wrap.
			memcpy(cell.vel+cg->start, vel+start, sizeof(velstruct)*cg->size());
			memcpy(cell.aux+cg->start, aux+start, sizeof(auxstruct)*cg->size());
			cellijk -= firstcell;
			if (cellijk.x> diam) cellijk.x-=GFC->cpd;
			if (cellijk.x<-diam) cellijk.x+=GFC->cpd;
			if (cellijk.y> diam) cellijk.y-=GFC->cpd;
			if (cellijk.y<-diam) cellijk.y+=GFC->cpd;
			if (cellijk.z> diam) cellijk.z-=GFC->cpd;
			if (cellijk.z<-diam) cellijk.z+=GFC->cpd;
			posstruct offset = GFC->invcpd*(cellijk);
			// printf("Using offset %f %f %f\n", offset.x, offset.y, offset.z);
			for (int p=0; p<cg->size(); p++) 
			     cell.pos[cg->start+p] = pos[start+p] - offset;
			start += cg->size();
		    } // End loop over cellgroups in this global group
		} // End loop over globalgroups in a cell
	// End loop over cells
	GFC->ScatterGroups.Stop();
	return;
    }


    // TODO: Need to alter GFC to initialize the params

    void FindSubGroups() {
	// Process each group, looking for L1 and L2 subgroups.
	FOFcell FOFlevel1[omp_get_max_threads()], FOFlevel2[omp_get_max_threads()];
	posstruct **L1pos;
	velstruct **L1vel;
	auxstruct **L1aux;
	#pragma omp parallel for schedule(static)
	for (int g=0; g<omp_get_num_threads(); g++) {
	    FOFlevel1[g].setup(GFC->linking_length_level1, 1e10);
	    FOFlevel2[g].setup(GFC->linking_length_level2, 1e10);
	    L1pos[g] = (posstruct *)malloc(sizeof(posstruct)*1e6);
	    L1vel[g] = (velstruct *)malloc(sizeof(velstruct)*1e6);
	    L1aux[g] = (auxstruct *)malloc(sizeof(auxstruct)*1e6);
	}
	
	#pragma omp parallel for schedule(static)
	for (int j=0; j<GFC->cpd; j++) {
	    int g = omp_get_thread_num();
	    for (int k=0; k<GFC->cpd; k++) {
		for (int n=0; n<globalgroups[j][k].size(); n++) {
		    if (globalgroups[j][k][n].np<GFC->minhalosize) continue;
		    // We have a large-enough global group to process
		    posstruct *grouppos = pos+globalgroups[j][k][n].start;
		    velstruct *groupvel = vel+globalgroups[j][k][n].start;
		    auxstruct *groupaux = aux+globalgroups[j][k][n].start;
		    int groupn = globalgroups[j][k][n].np;
		    FOFlevel1[g].findgroups(grouppos, NULL, NULL, groupn);
		    // Now we've found the L1 groups
		    for (int a=0; a<FOFlevel1[g].ngroups; a++) {
			int size = FOFlevel1[g].groups[a].n;
		        if (size<GFC->minhalosize) continue;
			// The group is big enough.
			FOFparticle *start = FOFlevel1[g].p+FOFlevel1[g].groups[a].start;
			// Particle indices are in start[0,size).index()
			// which deref the grouppos, groupvel, groupaux list.

			// We now need to find the L2 subgroups.
			// Make a temporary list of the particles in the L1 group
			for (int b=0; b<size; b++) {
			    L1pos[g][b] = grouppos[start[b].index()];
			    L1vel[g][b] = groupvel[start[b].index()];
			    L1aux[g][b] = groupaux[start[b].index()];
			}
			FOFlevel2[g].findgroups(L1pos[g], NULL, NULL, size);

			// TODO: Merger trees require tagging the L2 particles 
			// in the original aux array.  This can be done:
			// The L2 index() gives the position in the L1 array,
			// and that index() gets back to aux.

			ComputeStats(size, L1pos[g], L1vel[g], L1aux[g], FOFlevel2[g], slab, j, k);
		    } // Done with this L1 halo
		} // Done with this group
	    }
	}

	// TODO: destroy the FOFcells, if needed
    }

};


void SimpleOutput(int slab, GlobalGroupSlab &GGS) {
    // This just writes two sets of ASCII files to see the outputs,
    // but it also helps to document how the global group information is stored.
    char fname[200];
    sprintf(fname, "/tmp/out.group.%03d", slab);
    FILE *fp = fopen(fname,"w");
    for (int j=0; j<GFC->cpd; j++)
	for (int k=0; k<GFC->cpd; k++)
	    for (int n=0; n<GGS.globalgroups[j][k].size(); n++)
		fprintf(fp, "%d %d %d %d %d\n", (int)GGS.globalgroups[j][k][n].start,
				GGS.globalgroups[j][k][n].np, GGS.slab, j, k);
    fclose(fp);
    // Remember that these positions are relative to the first-cell position,
    // which is why we include that cell ijk in the first outputs
    /*
    sprintf(fname, "/tmp/out.pos.%03d", slab);
    fp = fopen(fname,"w");
    for (int p=0; p<GGS.np; p++)
	fprintf(fp, "%f %f %f %d\n", GGS.pos[p].x, GGS.pos[p].y, GGS.pos[p].z, (int)GGS.aux[p].pid());
    fclose(fp);
    */
}

void FindAndProcessGlobalGroups(int slab) {
    slab = GFC->WrapSlab(slab);
    GlobalGroupSlab GGS(slab);
    GGS.globalgroups.setup(GFC->cpd, 1024);   // TODO: Correct size to start?
    GGS.globalgrouplist.setup(GFC->cpd, 1024);   // TODO: Correct size to start?
    CreateGlobalGroups(GGS.slab, GGS.globalgroups, GGS.globalgrouplist);
    STDLOG(1,"Closed global groups in slab %d, finding %lld groups involving %lld cell groups\n", slab, GGS.globalgroups.size(), GGS.globalgrouplist.size());
    GFC->GGtot += GGS.globalgroups.get_slab_size();

    // Now process and output each one....
    // Start by gathering all of the particles into a contiguous set.
    GGS.GatherGlobalGroups();
    STDLOG(1,"Gathered %lld particles from global groups in slab %d\n", GGS.np, slab);
    // The GGS.globalgroups[j][k][n] now reference these as [start,start+np)

    // TODO: Lots to do here
    // For now, maybe we should just output the group multiplicity and the PIDs,
    // as a way of demonstrating that we have something.
    SimpleOutput(slab, GGS);

    GGS.ScatterGlobalGroups();
}


/* Having found the global groups, we should:

Gather a list of the particles, including the cell-offset positions.

Compute sub-group finding on that list, perhaps just returning a list of 
particle indices.

Tag the tagable particles in the L2 groups.
Compute statistics on the L1 sub-groups (do we need to make a permuted list
for this, or just dereference the indices?).
Output the L1 group information.
Output the tagable particles in the groups.

Execute microstepping on the groups.

Copy those pos/vel back into the original slabs, removing cell-offset positions.

*/







#ifdef DONT_COMPILE


// ======================= Code for the top-level driver ===============
// TODO: Move this to timestep.cpp and add to the setup and looping.
// TODO: Need to add timing and logging.

Dependency MakeCellGroups;
int MakeCellGroupsPrecondition(int slab) {
    // Need the new velocities because we're going to rearrange particles
    if( Kick.notdone(slab) ) return 0;
    
    // Also need the auxs, because we're going to re-order
    if( !LBW->IOCompleted( AuxSlab, slab ) ) {
        Dependency::NotifySpinning(WAITING_FOR_IO);
        return 0;
    }
    return 1;
    // TODO: Is it ok that the AccSlab is not being reordered now?
}

void MakeCellGroupAction(int slab) {
    GFC->CreateCellGroups(slab);
    return;
}

Dependency FindGroupLinks;
int FindGroupLinksPrecondition(int slab) {
    // We want to find all links between this slab and the one just behind
    for (int j=-1; j<=0; j++) 
	if (MakeCellGroups.notdone(slab+j)) return 0;
    return 1;
}

void FindGroupLinksAction(int slab) {
    // Find links between slab and slab-1
    FindGroupLinks(slab);
    return;
}
 
Dependency DoGlobalGroups;
int DoGlobalGroupsPrecondition(int slab) {
    // We're going to close all CellGroups in this slab.
    // GlobalGroups can span 2*GroupRadius+1.
    // But even though we usually encounter a CellGroup in its minimum slab,
    // we could be anywhere in the first instance.  So we have to query a big range.
    for (int j=-2*GFC->GroupRadius; j<=2*GFC->GroupRadius; j++) 
	if (FindGroupLinks.notdone(slab+j)) return 0;
    return 1;
}
void DoGlobalGroupsAction(int slab) {
    FindAndProcessGlobalGroups(slab);
    // TODO: LightCone Outputs might have to go here, since we need the GlobalGroups

    // At this point, all CellGroups in this slab have been closed and used
    // so we can delete them.
    GFC->DestroyCellGroups(slab);
    return;
}

#endif 
