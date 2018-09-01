/** \file This is the top-level code for group-finding.
We include this file in the program, and it includes the rest.
*/

// Use GLOG instead of STDLOG to write to the lastrun.groupstats file
#define GLOG(verbosity,...) { if (verbosity<=stdlog_threshold_global) { \
	LOG(*grouplog,__VA_ARGS__); grouplog->flush(); } }


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
class MicrostepControl;


/** This is the control class for all group finding.

It contains slab-based lists of the CellGroups and GlobalGroups,
as well as the GroupLinkList.  It also contains stats and timings
and generates a detailed log report.

*/

class GroupFindingControl {
  public:
    int cpd;		
    FOFloat linking_length;    ///< In code units
    FOFloat boundary;	///< The distance from the origin that indicates that something
    		///< is within linking_length of the edge
    int GroupRadius;    ///< Cell radius for group finding
    FOFloat invcpd;
    uint64 np; ///< the number of particles in the simulation; used to compute buffers
    uint64 particles_per_pencil;	///< Typical number of particles per pencil
    uint64 particles_per_slab;	///< Typical number of particles per slab

    // Parameters for finding subgroups
    FOFloat linking_length_level1;
    FOFloat linking_length_level2;
    FOFloat SOdensity1;	///< The density threshold for SO L1
    FOFloat SOdensity2;	///< The density threshold for SO L2
    int minhalosize;	///< The minimum size for a level 1 halo to be outputted

    uint64 pPtot, fPtot, fGtot, CGtot, GGtot, Ltot, CGactive;
    int largest_GG;
    FLOAT maxFOFdensity;
    double meanFOFdensity;

    MultiplicityStats L0stats, L1stats;
    long long numdists1, numdists2;	
    	///< The total number of distances computed
    long long numsorts1, numsorts2;	
    	///< The total number of sorting elements 
    long long numcenters1, numcenters2;	
    	///< The total number of centers considered

    SlabAccum<CellGroup> *cellgroups;   // Allocated [0,cpd), one for each slab
    int *cellgroups_status;     // Allocated [0,cpd) for each slab. 
    	// 0 means CellGroups not found.
    	// 1 means CellGroups found but slab not closed.
    	// 2 means slab closed (so CellGroups may be deleted).
    GroupLinkList *GLL;

    GlobalGroupSlab **globalslabs;    // Allocated [0,cpd), one for each slab
    MicrostepControl **microstepcontrol;

    STimer CellGroupTime, CreateFaceTime, FindLinkTime, 
    	SortLinks, IndexLinks, FindGlobalGroupTime, IndexGroups, 
	GatherGroups, ScatterAux, ScatterGroups, ProcessLevel1, OutputLevel1;
    PTimer L1FOF, L2FOF, L1Tot;
    PTimer IndexLinksSearch, IndexLinksIndex;
	
	std::ofstream *grouplog;

    // number of timeslice output particles written to disk
    uint64 n_L0_output = 0;

    void setupGGS();

    /// The constructor, where we load the key parameters
    ///
    /// The 2nd and 3rd parameters give the L1 and L2 parameters,
    /// either for FOF linking length or SO overdensity, depending on 
    /// that global definition.
    /// FOF lengths should be comoving lengths in code units (unit box)
    GroupFindingControl(FOFloat _linking_length, 
    	FOFloat _linking_length_level1, FOFloat _linking_length_level2,
    int _cpd, int _GroupRadius, int _minhalosize, uint64 _np) {
#ifdef STANDALONE_FOF
    grouplog = &stdlog;
#else
    std::string glogfn = P.LogDirectory;
    glogfn += "/lastrun.groupstats";
    grouplog = new std::ofstream();
    grouplog->open(glogfn);
#endif

	char onoff[5];
	#ifdef AVX_FOF
	    sprintf(onoff, "on");
	#else
	    sprintf(onoff, "off");
	#endif
	STDLOG(0,"Group finding sizeof(FOFloat)=%d, sizeof(FLOAT)=%d, AVX_FOF is %s\n", sizeof(FOFloat), sizeof(FLOAT), onoff);

	cpd = _cpd; 
	linking_length = _linking_length;
	#ifdef SPHERICAL_OVERDENSITY
	    SOdensity1 = _level1;
	    SOdensity2 = _level2;
	    STDLOG(0,"Planning for L1/2 group finding with SO: %f and %f\n", 
	    	SOdensity1, SOdensity2);
	#else
	    linking_length_level1 = _level1;
	    linking_length_level2 = _level2;
	    STDLOG(0,"Planning for L1/2 group finding with FOF: %f and %f\n", 
	    	linking_length_level1, linking_length_level2);
	#endif
	minhalosize = _minhalosize;
	invcpd = 1. / (double) _cpd;
	boundary = (invcpd/2.0-linking_length);
	np = _np;
	particles_per_pencil = np/cpd/cpd;
	particles_per_slab = np/cpd;
	GroupRadius = _GroupRadius;
	cellgroups = new SlabAccum<CellGroup>[cpd];
	cellgroups_status = new int[cpd];
	for (int j=0;j<cpd;j++) cellgroups_status[j] = 0;
	GLL = new GroupLinkList(cpd, np/cpd*linking_length/invcpd*3*15);    
    STDLOG(1,"Allocated %.2f GB for GroupLinkList\n", sizeof(GroupLink)*GLL->maxlist/1024./1024./1024.)
	
	setupGGS();
	// This is a MultiAppendList, so the buffer cannot grow. 
	// It will be storing links between cells.  Most particles do 
	// not generate such a link; on average 3*linking_length/cellsize do. 
	// But we will have up to 10 slabs worth in here.
	// So we need to be a bit generous.
        pPtot = fPtot = fGtot = CGtot = GGtot = Ltot = 0;
	CGactive = 0;
	maxFOFdensity = 0.0;
	largest_GG = 0;
	numdists1 = numsorts1 = numcenters1 = 0;
	numdists2 = numsorts2 = numcenters2 = 0;
	return;
    }

    ~GroupFindingControl();

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
	
    /// This generates the log report
    void report() {
	 GLOG(0,"Considered %f G particles as active\n", CGactive/1e9);
	 // The FOFdensities are weighted by b^2-r^2.  When integrated,
	 // that yields a mass at unit density of 
   	 // (2/15)*4*PI*b^5*np
	 float FOFunitdensity = P.np*4.0*M_PI*2.0/15.0*pow(WriteState.DensityKernelRad2,2.5)+1e-30;
	 GLOG(0,"Maximum reported density = %f (%e in code units)\n", maxFOFdensity/FOFunitdensity, maxFOFdensity);
	 meanFOFdensity /= P.np-WriteState.DensityKernelRad2;
	 GLOG(0,"Mean reported non-self density = %f (%e in code units)\n", meanFOFdensity/FOFunitdensity, meanFOFdensity);
	 GLOG(0,"Found %f G cell groups (including boundary singlets)\n", CGtot/1e9);
	 GLOG(0,"Used %f G pseudoParticles, %f G faceParticles, %f G faceGroups\n",
	     pPtot/1e9, fPtot/1e9, fGtot/1e9);
	 GLOG(0,"Found %f M links between groups.\n", Ltot/1e6);
	 GLOG(0,"Found %f M global groups\n", GGtot/1e6);
	 GLOG(0,"Longest GroupLink list was %f M, compared to %f M allocation (%f MB)\n", GLL->longest/1e6, GLL->maxlist/1e6, GLL->maxlist/1024/1024*sizeof(GroupLink));
	 GLOG(0,"Largest Global Group has %d particles\n", largest_GG);

	 GLOG(0,"L0 group multiplicity distribution:\n");
	 L0stats.report_multiplicities(grouplog);

	 GLOG(0,"L1 & L2 groups min size = %d\n", minhalosize);
	 GLOG(0,"L1 groups required %f G distances, %f G sorts, %f G centers\n", numdists1/1e9, numsorts1/1e9, numcenters1/1e9);
	 GLOG(0,"L2 groups required %f G distances, %f G sorts, %f G centers\n", numdists2/1e9, numsorts2/1e9, numcenters2/1e9);
	 GLOG(0,"L1 group multiplicity distribution:\n");
	 L1stats.report_multiplicities(grouplog);

	 float total_time = CellGroupTime.Elapsed()+
			CreateFaceTime.Elapsed()+
			FindLinkTime.Elapsed()+
			SortLinks.Elapsed()+
			IndexLinks.Elapsed()+
			FindGlobalGroupTime.Elapsed()+
			IndexGroups.Elapsed()+
			GatherGroups.Elapsed()+
			ProcessLevel1.Elapsed()+
            ScatterAux.Elapsed()+
			ScatterGroups.Elapsed();
	 GLOG(0,"Timings: \n");
	 #define RFORMAT(a) a.Elapsed(), a.Elapsed()/total_time*100.0
	 GLOG(0,"Finding Cell Groups:     %8.4f sec (%5.2f%%)\n",
			RFORMAT(CellGroupTime));
	 GLOG(0,"Creating Faces:          %8.4f sec (%5.2f%%)\n",
			RFORMAT(CreateFaceTime));
	 GLOG(0,"Finding Group Links:     %8.4f sec (%5.2f%%)\n",
			RFORMAT(FindLinkTime));
	 GLOG(0,"Sort Links:              %8.4f sec (%5.2f%%)\n",
			RFORMAT(SortLinks));
	 GLOG(0,"Index Links:             %8.4f sec (%5.2f%%)\n",
			RFORMAT(IndexLinks));
	 // printf("     Searching:               %8.4f sec\n", IndexLinksSearch.Elapsed());
	 GLOG(0,"Indexing (P):                %8.4f sec\n", IndexLinksIndex.Elapsed());
	 GLOG(0,"Find Global Groups:      %8.4f sec (%5.2f%%)\n",
			RFORMAT(FindGlobalGroupTime));
	 GLOG(0,"Index Global Groups:     %8.4f sec (%5.2f%%)\n",
			RFORMAT(IndexGroups));
	 GLOG(0,"Gather Group Particles:  %8.4f sec (%5.2f%%)\n",
			RFORMAT(GatherGroups));
	 GLOG(0,"Level 1 & 2 Processing:  %8.4f sec (%5.2f%%)\n",
			RFORMAT(ProcessLevel1));
	 GLOG(0,"Level 1 FOF (P):               %8.4f sec (%5.2f%%)\n",
			RFORMAT(L1FOF));
	 GLOG(0,"Level 2 FOF (P):               %8.4f sec (%5.2f%%)\n",
			RFORMAT(L2FOF));
	 GLOG(0,"Level 1 Total (P):             %8.4f sec (%5.2f%%)\n",
			RFORMAT(L1Tot));
	 GLOG(0,"Level 1 Output:          %8.4f sec (%5.2f%%)\n",
            RFORMAT(OutputLevel1));
	 GLOG(0,"Scatter Aux:             %8.4f sec (%5.2f%%)\n",
            RFORMAT(ScatterAux));
	 GLOG(0,"Scatter Group Particles: %8.4f sec (%5.2f%%)\n",
			RFORMAT(ScatterGroups));
	 GLOG(0,"Total Booked Time:       %8.4f sec (%5.2f Mp/sec)\n", total_time, np/total_time*1e-6);
	 #undef RFORMAT
    }
};

/// We keep the code for Constructing CellGroups in here, because
/// it is where the cellgroups[] are created and destroyed.

void GroupFindingControl::ConstructCellGroups(int slab) {
    // Construct the Cell Groups for this slab
    CellGroupTime.Start();
    slab = WrapSlab(slab);
    cellgroups[slab].setup(cpd, particles_per_slab/20);     
    	// Guessing that the number of groups is 20-fold less than particles
    FOFcell doFOF[omp_get_max_threads()];
    #pragma omp parallel for schedule(static,1)
    for (int g=0; g<omp_get_max_threads(); g++) 
    	doFOF[g].setup(linking_length, boundary);

    uint64 _CGactive = 0; 
    FLOAT _maxFOFdensity = 0.0;
    double _meanFOFdensity = 0.0;
    FLOAT DensityKernelRad2 = WriteState.DensityKernelRad2;
    #pragma omp parallel for schedule(dynamic,1) reduction(+:_CGactive) reduction(max:_maxFOFdensity) reduction(+:_meanFOFdensity)
    for (int j=0; j<cpd; j++) {
	int g = omp_get_thread_num();
        PencilAccum<CellGroup> *cg = cellgroups[slab].StartPencil(j);
        for (int k=0; k<cpd; k++) {
	    // Find CellGroups in (slab,j,k).  Append results to cg.
	    Cell c = PP->GetCell(slab, j, k);

	    int active_particles = c.count();
	    #ifdef COMPUTE_FOF_DENSITY
	    if (DensityKernelRad2>0) {
	        // The FOF-scale density is in acc.w.  
		// All zeros cannot be in groups, partition them to the end
		for (int p=0; p<active_particles; p++) {
		    _meanFOFdensity += c.acc[p].w;
			// This will be the mean over all particles, not just
			// the active ones
		    if (c.acc[p].w>DensityKernelRad2) {
			// Active particle; retain and accumulate stats
			_maxFOFdensity=std::max(_maxFOFdensity, c.acc[p].w);
		    } else {
		        // We found an inactive particle; swap to end.
			active_particles--;
			std::swap(c.pos[p], c.pos[active_particles]);
			std::swap(c.vel[p], c.vel[active_particles]);
			std::swap(c.aux[p], c.aux[active_particles]);
			std::swap(c.acc[p], c.acc[active_particles]);
			p--; // Need to try this location again
		   }
		}
	    }
	    // By limiting the particles being considered, we automatically
	    // will exclude these particles from being in the boundary
	    // singlet set.  So they will not be in the CellGroups and 
	    // hence never considered further.
	    #endif
	    _CGactive += active_particles;

	    doFOF[g].findgroups(c.pos, c.vel, c.aux, c.acc, active_particles);
	    // We need to clear the L0 & L1 bits for this timestep
	    for (int p=0; p<c.count(); p++) c.aux[p].reset_L01_bits();
	    for (int gr=0; gr<doFOF[g].ngroups; gr++) {
		CellGroup tmp(doFOF[g].groups[gr], boundary);
		cg->append(tmp);
	    }
	    // Also need to look at the singlets!
	    for (int p=doFOF[g].nmultiplets; p<doFOF[g].nsinglet_boundary; p++) {
		CellGroup tmp(p, c.pos[p], boundary);
	        cg->append(tmp);
	    }
	    cg->FinishCell();
	}
	cg->FinishPencil();
    }
    // Best if we destroy on the same thread, for tcmalloc
    #pragma omp parallel for schedule(static,1)
    for (int g=0; g<omp_get_max_threads(); g++) 
    	doFOF[g].destroy();
    uint64 tot = cellgroups[slab].get_slab_size();
    CGtot += tot;
    cellgroups_status[slab] = 1;
    CGactive += _CGactive;
    meanFOFdensity += _meanFOFdensity;
    maxFOFdensity = std::max(maxFOFdensity, _maxFOFdensity);
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

#ifdef STANDALONE_FOF
GroupFindingControl *GFC = NULL;	
	// We have one global instance of this, initialized in proepi
#endif

#include "findgrouplinks.cpp"
	// Code to search between pairs of cells and find the linked groups,
	// which get added to GLL.


// ===================== Output Field Particles ===============

#include "halostat.cpp"
	// Code to compute L1 halo properties

/** Gather all of the taggable particles that aren't in L1 groups into
two vectors, converting to global positions.  

Space must be allocated beforehand.  Returns the number of elements used.
Warning: This must be called after ScatterGlobalGroupsAux() and before
ScatterGlobalGroups()
*/

uint64 GatherTaggableFieldParticles(int slab, RVfloat *pv, TaggedPID *pid, FLOAT unkickfactor) {
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
            if(c.acc != NULL)
                v -= unkickfactor*TOFLOAT3(c.acc[p]);
		    pv[nfield] = RVfloat(r.x, r.y, r.z, v.x, v.y, v.z);
		    pid[nfield] = c.aux[p].pid();
		    nfield++;
		}
	}
    return nfield;
}

// ===================== Global Groups ===============

#include "globalgroup.cpp"

/// Code to traverse the links and find the GlobalGroups as 
/// sets of CellGroups (stored by their LinkIDs).

void GroupFindingControl::setupGGS() {
    globalslabs = new GlobalGroupSlab*[cpd];
    microstepcontrol = new MicrostepControl*[cpd];
    for(int i = 0; i < cpd; i++){
        globalslabs[i] = NULL;
        microstepcontrol[i] = NULL;
    }
}

GroupFindingControl::~GroupFindingControl() { 
    delete[] cellgroups;
    for (int j=0;j<cpd;j++) assert(cellgroups_status[j] == 2);
    delete[] cellgroups_status;
    delete GLL;
    #ifndef STANDALONE_FOF
    grouplog->close();
    delete grouplog;
    #endif

    delete[] globalslabs;
}

/** The driver routine to find the GlobalGroups in a slab
and then find and output the subgroups.
*/

void FindAndProcessGlobalGroups(int slab) {
    slab = GFC->WrapSlab(slab);
    assert(GFC->globalslabs[slab] == NULL);
    GlobalGroupSlab *GGS = new GlobalGroupSlab;
    GFC->globalslabs[slab] = GGS;
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
	
	int do_output;
#ifndef STANDALONE_FOF
    do_output = ReadState.DoGroupFindingOutput;
#else
	do_output = 1;
#endif
	if(do_output)
		GGS->FindSubGroups();
    GGS->ScatterGlobalGroupsAux();

    #ifdef ASCII_TEST_OUTPUT
    GGS->SimpleOutput();
    #endif
	if(do_output)
		GGS->HaloOutput();

#ifndef STANDALONE_FOF
    // We have a split time slice output model where non-L0 particles and L0 particles go in separate files
    if(ReadState.DoTimeSliceOutput){
        FLOAT unkickfactor = WriteState.FirstHalfEtaKick;
        STDLOG(1,"Outputting L0 group particles in slab %d with unkick factor %f\n", slab, unkickfactor);
        GFC->n_L0_output += GGS->L0TimeSliceOutput(unkickfactor);
    }
#endif
}

void FinishGlobalGroups(int slab){
	slab = GFC->WrapSlab(slab);
    GlobalGroupSlab *GGS = GFC->globalslabs[slab];

	// pos,vel have been updated in the group-local particle copies by microstepping
    // now push these updates to the original slabs
    GGS->ScatterGlobalGroups();
    delete GGS;
    GFC->globalslabs[slab] = NULL;
}
