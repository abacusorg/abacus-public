/** \file This is the top-level code for group-finding.
We include this file in the program, and it includes the rest.
*/

// Use GLOG instead of STDLOG to write to the lastrun.groupstats file
// #define GLOG(verbosity,...) { if (verbosity<=stdlog_threshold_global) { \
//         LOG(*grouplog,__VA_ARGS__); grouplog->flush(); } }

#define GLOG(verbosity,...) { if (verbosity<=stdlog_threshold_global) fmt::print(reportfp,__VA_ARGS__); }


#include "fof_sublist.cpp"
        // Code to do FOF on a set (e.g., a cell)

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
    int zstart;  ///< z domain includes ghosts
    int zend;
    int zwidth;

    FOFloat linking_length;    ///< In code units
    FOFloat boundary;        ///< The distance from the origin that indicates that something
                    ///< is within linking_length of the edge
    int GroupRadius;    ///< Cell radius for group finding
    FOFloat invcpd;
    uint64 np; ///< the number of particles in the simulation; used to compute buffers
    uint64 particles_per_pencil;        ///< Typical number of particles per pencil
    uint64 particles_per_slab;        ///< Typical number of particles per slab

    // Parameters for finding subgroups
    FOFloat linking_length_level1;
    FOFloat linking_length_level2;
    FOFloat SOdensity1;        ///< The density threshold for SO L1
    FOFloat SOdensity2;        ///< The density threshold for SO L2
    int minhalosize;        ///< The minimum size for a level 1 halo to be outputted

    // Globals for SO finding
    FOFloat FOFhalfcell;     // The half cell size
    FOFloat SOpartition;  // The radial binning

    uint64 pPtot, fPtot, fGtot, CGtot, GGtot, Ltot, CGactive;
    int largest_GG;
    FLOAT maxFOFdensity;
    double meanFOFdensity;

    MultiplicityStats L0stats, L1stats;
    // BOT
    MultiplicityStats L2stats;
  
    long long numdists1, numdists2;        
            ///< The total number of distances computed
    long long numsorts1, numsorts2;        
            ///< The total number of sorting elements 
    long long numcenters1, numcenters2;        
            ///< The total number of centers considered
    long long numcg1, numcg2;        
            ///< The total number of cell groups used
    long long numgroups1, numgroups2;        
            ///< The total number of groups found
    int max_group_diameter;

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
        DeferGroups, ClearDefer,
        GatherGroups, ScatterAux, ScatterGroups, ProcessLevel1, OutputLevel1;
    PTimer L1FOF, L2FOF, L1Tot;
    PTimer IndexLinksSearch, IndexLinksIndex;
        
        // std::ofstream *grouplog;

    // number of timeslice output particles written to disk
    uint64 n_L0_output = 0;

    int use_aux_dens; ///< Whether acc has the density, or it's already in the aux

    void setupGGS();

    /// The constructor, where we load the key parameters
    ///
    /// The 2nd and 3rd parameters give the L1 and L2 parameters,
    /// either for FOF linking length or SO overdensity, depending on 
    /// that global definition.
    /// FOF lengths should be comoving lengths in code units (unit box)
    GroupFindingControl(FOFloat _linking_length, 
            FOFloat _level1, FOFloat _level2,
        int _cpd, int _zstart, int _zwidth,
        int _GroupRadius, int _minhalosize, uint64 _np, int _use_aux_dens) {
        /*
    #ifdef STANDALONE_FOF
        grouplog = &stdlog;
    #else
        fs::path glogfn = P.LogDirectory;
        glogfn += "/lastrun"; glogfn += NodeString; glogfn += ".groupstats";
        grouplog = new std::ofstream();
        grouplog->open(glogfn);
    #endif
    */
    


        char onoff[5];
        #ifdef AVXFOF
            sprintf(onoff, "on");
        #else
            sprintf(onoff, "off");
        #endif
        STDLOG(0,"Group finding sizeof(FOFloat)={:d}, sizeof(FLOAT)={:d}, AVXFOF is {:s}\n", sizeof(FOFloat), sizeof(FLOAT), onoff);

        cpd = _cpd; 
        zstart = _zstart;
        zwidth = _zwidth;
        zend = zstart + zwidth;
        
        linking_length = _linking_length;
        #ifdef SPHERICAL_OVERDENSITY
            SOdensity1 = _level1;
            SOdensity2 = _level2;
            STDLOG(0,"Planning for L1/2 group finding with SO: {:f} and {:f}\n", 
                SOdensity1, SOdensity2);
        #else
            linking_length_level1 = _level1;
            linking_length_level2 = _level2;
            STDLOG(0,"Planning for L1/2 group finding with FOF: {:f} and {:f}\n", 
                linking_length_level1, linking_length_level2);
        #endif
        FOFhalfcell = FOF_RESCALE/2.0*CP->invcpd;     // The half cell size
        SOpartition = FOFhalfcell*2.0*sqrt(3.0)/3.0;  // The radial binning
        STDLOG(1,"SO parameters: FOFhalfcell = {:f}, SOpartition = {:f}\n", FOFhalfcell, SOpartition);

        minhalosize = _minhalosize;
        invcpd = 1. / (double) _cpd;
        boundary = (invcpd/2.0-linking_length);
        np = _np;
        particles_per_pencil = np*zwidth/cpd/cpd/cpd;
        particles_per_slab = np*zwidth/cpd/cpd;
        GroupRadius = _GroupRadius;
        cellgroups = new SlabAccum<CellGroup>[cpd];
        cellgroups_status = new int[cpd];
        for (int j=0;j<cpd;j++) cellgroups_status[j] = 0;
        GLL = new GroupLinkList(cpd, np/cpd*linking_length/invcpd*3*15);  // TODO: does not appear to rescale for 1D or 2D
        STDLOG(1,"Allocated {:.2f} GB for GroupLinkList\n", sizeof(GroupLink)*GLL->maxlist/1024./1024./1024.);
        
        setupGGS();
        // This is a MultiAppendList, so the buffer cannot grow. 
        // It will be storing links between cells.  Most particles do 
        // not generate such a link; on average 3*linking_length/cellsize do. 
        // But we will have up to 10 slabs worth in here.
        // So we need to be a bit generous.
        pPtot = fPtot = fGtot = CGtot = GGtot = Ltot = 0;
        CGactive = 0;
        maxFOFdensity = 0.0;
        meanFOFdensity = 0.;
        largest_GG = 0;
        numdists1 = numsorts1 = numcenters1 = numcg1 = numgroups1 = 0;
        numdists2 = numsorts2 = numcenters2 = numcg2 = numgroups2 = 0;
        max_group_diameter = 0;

        use_aux_dens = _use_aux_dens;
        STDLOG(1,"Using {:s} densities in group finding\n", use_aux_dens ? "aux" : "acc");

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
    void report(FILE *reportfp) {
         GLOG(0,"Considered {:f} G particles as active\n", CGactive/1e9);
         // The FOFdensities are weighted by b^2-r^2.  When integrated,
         // that yields a mass at unit density of 
            // (2/15)*4*PI*b^5*np
         GLOG(0,"Maximum reported density = {:f} ({:e} in code units)\n", maxFOFdensity/WriteState.FOFunitdensity, maxFOFdensity);
         meanFOFdensity /= P.np;
         GLOG(0,"Mean reported non-self density = {:f} ({:e} in code units)\n", meanFOFdensity/WriteState.FOFunitdensity, meanFOFdensity);
         GLOG(0,"Found {:f} G cell groups (including boundary singlets)\n", CGtot/1e9);
         GLOG(0,"Used {:f} G pseudoParticles, {:f} G faceParticles, {:f} G faceGroups\n",
             pPtot/1e9, fPtot/1e9, fGtot/1e9);
         GLOG(0,"Found {:f} M links between groups.\n", Ltot/1e6);
         GLOG(0,"Found {:f} M global groups\n", GGtot/1e6);
         GLOG(0,"Longest GroupLink list was {:f} M, compared to {:f} M allocation ({:f} MB)\n", GLL->longest/1e6, GLL->maxlist/1e6, GLL->maxlist/1024./1024.*sizeof(GroupLink));
         GLOG(0,"Widest L0 Diameter reached {:d} slabs from the first\n", max_group_diameter);

#ifdef PARALLEL
         MPI_REDUCE_TO_ZERO(&max_group_diameter, 1, MPI_INT, MPI_MAX);
#endif
         WriteState.MaxGroupDiameter = max_group_diameter; 

         GLOG(0,"Largest Global Group has {:d} particles\n", largest_GG);
         WriteState.MaxL0GroupSize = largest_GG; 

         GLOG(0,"L0 group multiplicity distribution:\n");
         L0stats.report_multiplicities(reportfp);

         GLOG(0,"L1 & L2 groups min size = {:d}\n", minhalosize);
         GLOG(0,"L1 groups required {:f} G distances, {:f} G sorts, {:f} G centers, {:f} G cg\n", numdists1/1e9, numsorts1/1e9, numcenters1/1e9, numcg1/1e9);
         GLOG(0,"L2 groups required {:f} G distances, {:f} G sorts, {:f} G centers, {:f} G cg\n", numdists2/1e9, numsorts2/1e9, numcenters2/1e9, numcg2/1e9);
         GLOG(0,"L1 group multiplicity distribution:\n");
         GLOG(0,"Total number of L1 groups considered {:f} M\n", numgroups1/1e6);
         L1stats.report_multiplicities(reportfp);
         GLOG(0,"L2 group multiplicity distribution:\n");
         GLOG(0,"Total number of L2 groups considered {:f} M\n", numgroups2/1e6);
     L2stats.report_multiplicities(reportfp);
     
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
     float total_cycle;
         GLOG(0,"Timings: \n");
         #define RFORMAT(a) a.Elapsed(), a.Elapsed()/total_time*100.0
         #define CFORMAT(a) a.Elapsed(), a.Elapsed()/total_cycle*100.0
         GLOG(0,"Finding Cell Groups:     {:8.4f} sec ({:5.2f}%)\n",
                        RFORMAT(CellGroupTime));
         GLOG(0,"Creating Faces:          {:8.4f} sec ({:5.2f}%)\n",
                        RFORMAT(CreateFaceTime));
         GLOG(0,"Finding Group Links:     {:8.4f} sec ({:5.2f}%)\n",
                        RFORMAT(FindLinkTime));
         GLOG(0,"Sort Links:              {:8.4f} sec ({:5.2f}%)\n",
                        RFORMAT(SortLinks));
         GLOG(0,"Index Links:             {:8.4f} sec ({:5.2f}%)\n",
                        RFORMAT(IndexLinks));
         // fmt::print("     Searching:               {:8.4f} sec\n", IndexLinksSearch.Elapsed());
         GLOG(0,"Indexing (P):                {:8.4g} cyc\n", IndexLinksIndex.Elapsed());
         GLOG(0,"Defer Groups:            {:8.4f} sec ({:5.2f}%)\n",
                        RFORMAT(DeferGroups));
         GLOG(0,"Find Global Groups:      {:8.4f} sec ({:5.2f}%)\n",
                        RFORMAT(FindGlobalGroupTime));
         GLOG(0,"Clear Deferrals:         {:8.4f} sec ({:5.2f}%)\n",
                        RFORMAT(ClearDefer));
         GLOG(0,"Index Global Groups:     {:8.4f} sec ({:5.2f}%)\n",
                        RFORMAT(IndexGroups));
         GLOG(0,"Gather Group Particles:  {:8.4f} sec ({:5.2f}%)\n",
                        RFORMAT(GatherGroups));
         GLOG(0,"Level 1 & 2 Processing:  {:8.4f} sec ({:5.2f}%)\n",
                        RFORMAT(ProcessLevel1));
     total_cycle = L1Tot.Elapsed();
         GLOG(0,"Level 1 FOF (P):               {:8.4g} cyc ({:5.2f}%)\n",
                        CFORMAT(L1FOF));
         GLOG(0,"Level 2 FOF (P):               {:8.4g} cyc ({:5.2f}%)\n",
                        CFORMAT(L2FOF));
         GLOG(0,"Level 1 Total (P):             {:8.4g} cyc ({:5.2f}%)\n",
                        CFORMAT(L1Tot));
         GLOG(0,"Level 1 Output:          {:8.4f} sec ({:5.2f}%)\n",
            RFORMAT(OutputLevel1));
         GLOG(0,"Scatter Aux:             {:8.4f} sec ({:5.2f}%)\n",
            RFORMAT(ScatterAux));
         GLOG(0,"Scatter Group Particles: {:8.4f} sec ({:5.2f}%)\n",
                        RFORMAT(ScatterGroups));
         GLOG(0,"Total Booked Time:       {:8.4f} sec ({:5.2f} Mp/sec)\n", total_time, np/total_time*1e-6);
         #undef RFORMAT
    }
};

/// We keep the code for Constructing CellGroups in here, because
/// it is where the cellgroups[] are created and destroyed.

void GroupFindingControl::ConstructCellGroups(int slab) {
    // Construct the Cell Groups for this slab
    CellGroupTime.Start();
    slab = WrapSlab(slab);
    cellgroups[slab].setup(cpd, zwidth, particles_per_slab/20);
            // Guessing that the number of groups is 20-fold less than particles
    int nthread = omp_get_max_threads();
    FOFcell doFOF[nthread];
    #pragma omp parallel for schedule(static,1)
    for (int g=0; g<nthread; g++) 
            doFOF[g].setup(linking_length, boundary);

    FLOAT DensityKernelRad2 = WriteState.DensityKernelRad2;
    FLOAT L0DensityThreshold = WriteState.L0DensityThreshold;

    if (L0DensityThreshold>0) {
        // We've been asked to compare the kernel density to a particular
        // density (unit of cosmic mean) to make a particle L0 eligible.
        // Remember that this density does include the self-count.
        // The FOFdensities are weighted by b^2-r^2.  When integrated,
        // that yields a mass at unit density of 
        // (2/15)*4*PI*b^5*np
        L0DensityThreshold *= WriteState.FOFunitdensity;  // Now in code units
    } else {
        L0DensityThreshold = 0.0;  
        // Was DensityKernelRad2 but now the self-count has been subtracted.
    }

    NUMA_FOR(j,0,cpd, reduction(+:CGactive,meanFOFdensity) reduction(max:maxFOFdensity), FALLBACK_DYNAMIC){
        int g = omp_get_thread_num();
        PencilAccum<CellGroup> *cg = cellgroups[slab].StartPencil(j);
        for (int k=zstart; k<zend; k++) {  // global z
            // Find CellGroups in (slab,j,k).  Append results to cg.
            Cell c = CP->GetCell(slab, j, k);

            int active_particles = c.count();
            #ifdef COMPUTE_FOF_DENSITY
            if (DensityKernelRad2>0) {
                // The FOF-scale density is in acc.w.
                // acc is NULL in 2D group finding steps
                // All zeros cannot be in groups, partition them to the end
                for (int p=0; p<active_particles; p++) {
                    FLOAT dens;
                    if(use_aux_dens){
                        dens = c.aux[p].get_density();
                    } else {
                        dens = c.acc[p].w;
                        c.aux[p].set_compressed_density(dens);
                    }
                    meanFOFdensity += dens;
                        // This will be the mean over all particles, not just
                        // the active ones

                    if (dens>L0DensityThreshold) {
                        // Active particle; retain and accumulate stats
                        maxFOFdensity = std::max(maxFOFdensity, dens);
                    } else {
                        // We found an inactive particle; swap to end.
                        active_particles--;
                        std::swap(c.pos[p], c.pos[active_particles]);
                        std::swap(c.vel[p], c.vel[active_particles]);
                        std::swap(c.aux[p], c.aux[active_particles]);
                        if(c.acc) std::swap(c.acc[p], c.acc[active_particles]);
                        p--; // Need to try this location again
                   }
                }
            }
            // By limiting the particles being considered, we automatically
            // will exclude these particles from being in the boundary
            // singlet set.  So they will not be in the CellGroups and 
            // hence never considered further.
            #endif
            CGactive += active_particles;

            doFOF[g].findgroups(c.pos, c.vel, c.aux, c.acc, active_particles);

            // We need to clear the L0 & L1 bits for this timestep
        // This has been moved to the merge, so the bits never get written out
            // for (int p=0; p<c.count(); p++) c.aux[p].reset_L01_bits();

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
    NUMA_FOR_END;

    // Best if we destroy on the same thread, for tcmalloc
    #pragma omp parallel for schedule(static,1)
    for (int g=0; g<omp_get_max_threads(); g++) 
            doFOF[g].destroy();
    uint64 tot = cellgroups[slab].get_slab_size();
    CGtot += tot;
    cellgroups_status[slab] = 1;
    STDLOG(2,"Found {:d} cell groups in slab {:d}\n", tot, slab);
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


#include "sofind.cpp"

// ===================== Output Field Particles ===============

#include "halostat.cpp"
    // Code to compute L1 halo properties

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
    for (int j=0;j<cpd;j++) assert(cellgroups_status[j] != 1);
    delete[] cellgroups_status;
    delete GLL;
    /*
    #ifndef STANDALONE_FOF
    grouplog->close();
    delete grouplog;
    #endif
    */

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
    STDLOG(2,"Closed global groups in slab {:d}, finding {:d} groups involving {:d} cell groups\n",
        slab, GGS->globalgroups.get_slab_size(), GGS->globalgrouplist.get_slab_size());
    GFC->GGtot += GGS->globalgroups.get_slab_size();

    // Now process and output each one....
    // Start by gathering all of the particles into a contiguous set.
    GGS->GatherGlobalGroups();
    STDLOG(1,"Gathered {:d} particles from global groups in slab {:d}\n", GGS->np, slab);
    GFC->largest_GG = std::max(GFC->largest_GG, GGS->largest_group);
    // TODO: This largest_GG work is now superceded by MultiplicityHalos
    // The GGS->globalgroups[j][k][n] now reference these as [start,start+np)
        
    // ReadState.DoGroupFindingOutput is decided in InitGroupFinding()
    if(ReadState.DoGroupFindingOutput)
            GGS->FindSubGroups();
    GGS->ScatterGlobalGroupsAux();

    // Output the information about the Global Groups
    #ifdef ASCII_TEST_OUTPUT
        GGS->SimpleOutput();
    #endif
        if(ReadState.DoGroupFindingOutput)
            GGS->HaloOutput();

#ifndef STANDALONE_FOF
    // We have a split time slice output model where non-L0 particles and L0 particles go in separate files
    if(ReadState.DoTimeSliceOutput){
        FLOAT unkickfactor = WriteState.FirstHalfEtaKick;
        STDLOG(1,"Outputting L0 group particles in slab {:d} with unkick factor {:f}\n", slab, unkickfactor);
        GFC->n_L0_output += GGS->L0TimeSliceOutput(unkickfactor);
    }
#endif
}

void FinishGlobalGroups(int slab){
    slab = GFC->WrapSlab(slab);
    GlobalGroupSlab *GGS = GFC->globalslabs[slab];

    // pos,vel have been updated in the group-local particle copies by microstepping
    // now push these updates to the original slabs
    if(GFC->microstepcontrol[slab]!=NULL) {
        GGS->ScatterGlobalGroups();
    }
    delete GGS;
    GFC->globalslabs[slab] = NULL;
}
