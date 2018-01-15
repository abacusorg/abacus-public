// This is the top-level code for group-finding.
// We include this file in the program.


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

#include "multiplicity_stats.cpp"


class GlobalGroupSlab;


class GroupFindingControl {
    // This is the control class for all group finding.
  public:
    int cpd;		
    FOFloat linking_length;    // In code units
    FOFloat boundary;	// The distance from the origin that indicates that something
    		// is within linking_length of the edge
    int GroupRadius;    // Cell radius for group finding
    FOFloat invcpd;
    uint64 np; // the number of particles in the simulation; used to compute buffers
    uint64 particles_per_pencil;	// Typical number of particles per pencil

    // Parameters for finding subgroups
    FOFloat linking_length_level1;
    FOFloat linking_length_level2;
    int minhalosize;	// The minimum size for a level 1 halo to be outputted

    uint64 pPtot, fPtot, fGtot, CGtot, GGtot, Ltot;
    int largest_GG;

    MultiplicityStats L0stats, L1stats;

    SlabAccum<CellGroup> *cellgroups;   // Allocated [0,cpd), one for each slab
    int *cellgroups_status;     // Allocated [0,cpd) for each slab. 
    	// 0 means CellGroups not found.
    	// 1 means CellGroups found but slab not closed.
    	// 2 means slab closed (so CellGroups may be deleted).
    GroupLinkList *GLL;

    GlobalGroupSlab *globalslabs;    // Allocated [0,cpd), one for each slab

    STimer CellGroupTime, CreateFaceTime, FindLinkTime, 
    	SortLinks, IndexLinks, FindGlobalGroupTime, IndexGroups, 
	GatherGroups, ScatterGroups, ProcessLevel1;
    PTimer L1FOF, L2FOF, L1Tot;
    PTimer IndexLinksSearch, IndexLinksIndex;

    void setupGGS();

    GroupFindingControl(FOFloat _linking_length, 
    	FOFloat _linking_length_level1, FOFloat _linking_length_level2,
    int _cpd, FOFloat _invcpd, int _GroupRadius, int _minhalosize, uint64 _np) {
	cpd = _cpd; 
	linking_length = _linking_length;
	linking_length_level1 = _linking_length_level1;
	linking_length_level2 = _linking_length_level2;
	minhalosize = _minhalosize;
	invcpd = _invcpd;
	boundary = (invcpd/2.0-linking_length);
	np = _np;
	particles_per_pencil = np/cpd/cpd;
	GroupRadius = _GroupRadius;
	cellgroups = new SlabAccum<CellGroup>[cpd];
	cellgroups_status = new int[cpd];
	for (int j=0;j<cpd;j++) cellgroups_status[j] = 0;
	GLL = new GroupLinkList(cpd, np/cpd*linking_length/_invcpd*3*15);    
	setupGGS();
	// This is a MultiAppendList, so the buffer cannot grow. 
	// It will be storing links between cells.  Most particles do 
	// not generate such a link; on average 3*linking_length/cellsize do. 
	// But we will have up to 10 slabs worth in here.
	// So we need to be a bit generous.
        pPtot = fPtot = fGtot = CGtot = GGtot = Ltot = 0;
	largest_GG = 0;
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
	 STDLOG(0,"Found %d cell groups (including boundary singlets)\n", CGtot);
	 STDLOG(0,"Used %d pseudoParticles, %d faceParticles, %d faceGroups\n",
	     pPtot, fPtot, fGtot);
	 STDLOG(0,"Found %d links between groups.\n", Ltot);
	 STDLOG(0,"Found %d global groups\n", GGtot);
	 STDLOG(0,"Longest GroupLink list was %d, compared to %d allocation\n", GLL->longest, GLL->maxlist);
	 STDLOG(0,"Largest Global Group has %d particles\n", largest_GG);

	 STDLOG(0,"L0 group multiplicity distribution:\n");
	 L0stats.report_multiplicities();

	 STDLOG(0,"L1 group multiplicity distribution:\n");
	 L1stats.report_multiplicities();

	 float total_time = CellGroupTime.Elapsed()+
			CreateFaceTime.Elapsed()+
			FindLinkTime.Elapsed()+
			SortLinks.Elapsed()+
			IndexLinks.Elapsed()+
			FindGlobalGroupTime.Elapsed()+
			IndexGroups.Elapsed()+
			GatherGroups.Elapsed()+
			ProcessLevel1.Elapsed()+
			ScatterGroups.Elapsed();
	 STDLOG(0,"\nTimings: \n");
	 #define RFORMAT(a) a.Elapsed(), a.Elapsed()/total_time*100.0
	 STDLOG(0,"Finding Cell Groups:     %8.4f sec (%5.2f%%)\n",
			RFORMAT(CellGroupTime));
	 STDLOG(0,"Creating Faces:          %8.4f sec (%5.2f%%)\n",
			RFORMAT(CreateFaceTime));
	 STDLOG(0,"Finding Group Links:     %8.4f sec (%5.2f%%)\n",
			RFORMAT(FindLinkTime));
	 STDLOG(0,"Sort Links:              %8.4f sec (%5.2f%%)\n",
			RFORMAT(SortLinks));
	 STDLOG(0,"Index Links:             %8.4f sec (%5.2f%%)\n",
			RFORMAT(IndexLinks));
	 // printf("     Searching:               %8.4f sec\n", IndexLinksSearch.Elapsed());
	 STDLOG(0,"Indexing (P):                %8.4f sec\n", IndexLinksIndex.Elapsed());
	 STDLOG(0,"Find Global Groups:      %8.4f sec (%5.2f%%)\n",
			RFORMAT(FindGlobalGroupTime));
	 STDLOG(0,"Index Global Groups:     %8.4f sec (%5.2f%%)\n",
			RFORMAT(IndexGroups));
	 STDLOG(0,"Gather Group Particles:  %8.4f sec (%5.2f%%)\n",
			RFORMAT(GatherGroups));
	 STDLOG(0,"Level 1 & 2 Processing:  %8.4f sec (%5.2f%%)\n",
			RFORMAT(ProcessLevel1));
	 STDLOG(0,"Level 1 FOF:                  %8.4f sec (%5.2f%%)\n",
			RFORMAT(L1FOF));
	 STDLOG(0,"Level 2 FOF:                  %8.4f sec (%5.2f%%)\n",
			RFORMAT(L2FOF));
	 STDLOG(0,"Level 1 Total:                %8.4f sec (%5.2f%%)\n",
			RFORMAT(L1Tot));
	 STDLOG(0,"Scatter Group Particles: %8.4f sec (%5.2f%%)\n",
			RFORMAT(ScatterGroups));
	 STDLOG(0,"Total Booked Time:       %8.4f sec (%5.2f Mp/sec)\n", total_time, np/total_time*1e-6);
	 #undef RFORMAT
    }
};

// We keep the code for Constructing CellGroups in here, because
// it is where the cellgroups[] are created and destroyed.

void GroupFindingControl::ConstructCellGroups(int slab) {
    // Construct the Cell Groups for this slab
    CellGroupTime.Start();
    slab = WrapSlab(slab);
    cellgroups[slab].setup(cpd, particles_per_pencil);     
    FOFcell doFOF[omp_get_max_threads()];
    #pragma omp parallel for schedule(static)
    for (int g=0; g<omp_get_max_threads(); g++) 
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
	    doFOF[g].findgroups(c.pos, c.vel, c.aux, c.acc, c.count());
	    // printf("Cell %d %d %d: Found %d cell groups with %d particles, plus %d boundary singlets\n", slab,j,k,doFOF[g].ngroups, doFOF[g].nmultiplets, doFOF[g].nsinglet_boundary-doFOF[g].nmultiplets);
	    // We need to clear the L0 & L1 bits for this timestep
	    for (int p=0; p<c.count(); p++) c.aux[p].reset_L01_bits();
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
    #pragma omp parallel for schedule(static)
    for (int g=0; g<omp_get_max_threads(); g++) 
    	doFOF[g].destroy();
    uint64 tot = cellgroups[slab].get_slab_size();
    CGtot += tot;
    cellgroups_status[slab] = 1;
    STDLOG(1,"Found %d cell groups in slab %d\n", tot, slab);
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


// ===================== Output Field Particles ===============

#include "halostat.cpp"
	// Code to compute L1 halo properties

uint64 GatherTaggableFieldParticles(int slab, RVfloat *pv, TaggedPID *pid) {
    // Gather all of the taggable particles that aren't in L1 groups into two vectors,
    // converting to global positions.
    // Space must be allocated beforehand.
    // Returns the number of elements used.
    // Warning: This must be called after ScatterGlobalGroupsAux()
    // and before ScatterGlobalGroups()
    slab = GFC->WrapSlab(slab);
    uint64 nfield = 0;
    for (int j=0; j<GFC->cpd; j++)
	for (int k=0; k<GFC->cpd; k++) {
	    // Loop over cells
	    posstruct offset = PP->CellCenter(slab, j, k);
	    Cell c = PP->GetCell(slab, j, k);
	    for (int p=0; p<c.count(); p++)
		if (c.aux[p].is_taggable() && !c.aux[p].is_L1()) {
		    // We found a taggable field particle
		    posstruct r = c.pos[p] + offset;
		    velstruct v = c.vel[p];
		    pv[nfield] = RVfloat(r.x, r.y, r.z, v.x, v.y, v.z);
		    pid[nfield] = c.aux[p].pid();
		    nfield++;
		}
	}
    return nfield;
}

// ===================== Global Groups ===============

#include "globalgroup.cpp"
	// Code to traverse the links and find the GlobalGroups as 
	// sets of CellGroups (stored by their LinkIDs).

void GroupFindingControl::setupGGS() {
    globalslabs = new GlobalGroupSlab[cpd];
}

void FindAndProcessGlobalGroups(int slab) {
    slab = GFC->WrapSlab(slab);
    GlobalGroupSlab *GGS = GFC->globalslabs+slab;
    GGS->setup(slab);
    GGS->CreateGlobalGroups();
    STDLOG(1,"Closed global groups in slab %d, finding %d groups involving %d cell groups\n", slab, GGS->globalgroups.get_slab_size(), GGS->globalgrouplist.get_slab_size());
    GFC->GGtot += GGS->globalgroups.get_slab_size();

    // Now process and output each one....
    // Start by gathering all of the particles into a contiguous set.
    GGS->GatherGlobalGroups();
    STDLOG(1,"Gathered %d particles from global groups in slab %d\n", GGS->np, slab);
    GFC->largest_GG = std::max(GFC->largest_GG, GGS->largest_group);
    // TODO: This largest_GG work is now superceded by MultiplicityHalos
    // The GGS->globalgroups[j][k][n] now reference these as [start,start+np)
	
	// Check if, by going from ReadState to WriteState, we are crossing a L1Output_dlna checkpoint
	// Also always output if we're doing a TimeSlice output
	// Due to C's mod behavior, this assumes log(a) < 0
	int do_output = fmod(log(WriteState.ScaleFactor), P.L1Output_dlna) < fmod(log(ReadState.ScaleFactor), P.L1Output_dlna)
					|| log(WriteState.ScaleFactor) - log(ReadState.ScaleFactor) >= P.L1Output_dlna
					|| ReadState.DoTimeSliceOutput;
	if(do_output)
		GGS->FindSubGroups();
    GGS->ScatterGlobalGroupsAux();

    #ifdef ASCII_TEST_OUTPUT
    GGS->SimpleOutput();
    #endif
	if(do_output)
		GGS->HaloOutput();

    GGS->ScatterGlobalGroups();
    GGS->destroy();
    return;
}




