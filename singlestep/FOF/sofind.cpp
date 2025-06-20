// Copyright 2012-2025 The Abacus Developers
// SPDX-License-Identifier: GPL-3.0-or-later

/** \file Code to perform spherical overdensity finding on provided particle
sets, in particular the L0 groups.

The algorithms here are purely serial; it is expected that multithreading
will be applied to operate on different L0 groups with different instances
of SOfind.

We are provided with the density of each particle, using a 1-r^2/b^2
kernel where b is the FOF L0 linking length.  We use this to find the
densest particles.  For a uniform density, this kernel has a volume of
(2/5)(4*pi/3)b^3, which is (1.17b)^3.  Given that we use b around 0.20,
that's about 0.25 times the interparticle spacing.  Hence, a uniform
overdensity of 64 would yield about 1 count.  SO(180) is about 3 counts.
For large FOF groups, we expect the percolation density to be around 80,
so 1-2 counts.

If we think that uniform-density percolation happens when the average
counts in a sphere is 1, then we would find N = 2/5 average count in
such a circumstance.  Perhaps this estimate of percolation is missing
a constant of order unity; could consult the More++ paper for that.

Alternatively, if we have a group with 20 particles inside of R_180,
then using (4*pi/3) 180 R_180^3 = M_180 implies a radius R_180 of 
about 0.30 times the interparticle spacing.  If we now approximate
the group density profile as a singular isothermal sphere, then we
have M(<r) = M_180 (R/R_180) and rho(R) = (M_180/R_180) (1/4*pi*R^2).
We can then integrate to find the number of particles in the central
density, N_0 = \int_0^b 4*pi r^2 dr rho (1-r^2/b^2) = (M_180/R_180) (2/3b).

Hence, if we have M_180 = 20 and R_180 = 0.3, with b=0.2 we expect N_0 =
20*(4/9) = 9.  So our innermost density is around 9 for a small group.
Larger groups will have N_0 scaling as M_200^(2/3), so it will get quite
large.  That said, the SIS 1/r^2 density singularity is over-optimistic.
But the key thing is that we don't need to search down to the smallest
density particles.

These numbers all get modified by b^2 in code units to get to the code
densities.  We'll prefer to quote as a cosmic density.

For now, we tie the required minimum FOF density to be the same as the 
SO overdensity.  This is reasonably conservative, as for a big group,
one expects the local (FOF-scale) density at the edge to be ~3 times
less than the SO overdensity of that edge.  But it is also conservative
in the other direction, because for a peaked distribution even in a 
small group, the local density at the center should be higher, as
described above.  Still, this may be a tuning choice.

So far, we are using a greedy algorithm.  We find the densest particle, find
the distance to all others and sort in increasing order.  Then work from
the outside to find the largest radius that has an enclosed overdensity
above the threshold.  That defines the first group.  Then repeat with
the next densest particle, excluding those particles already in a 
group.

We are not implementing any unbinding here.

*/

#define SO_CACHE 1024

#define SO_UNASSIGNED 999999    // A huge number, more than we have groups

class alignas(16) SOcellgroup {
  public:
    FOFparticle cellcenter; // The center of the cell in GlobalGroup coords 
    //FOFloat d2;         // The minimum distance to the current halo nucleus position
    int firstbin;       // The radial bin this cell first appears in.
    int start[5];       // The starting point of particles in the L0 list
    int active;         // Has a cell been touched for this halo center
    FOFloat d2_furthest;  // furthest distance of particle in cell from current nucleus
    
    // We will allow 4 partitions of this list: 
    // [ start[j], start[j+1] ) for j=0,1,2,3.
    // So start[4] is start[0]+np_in_cellgroup
    
    // Null constructors & destructors
    SOcellgroup() = default;
    ~SOcellgroup() = default;

    void load(int start_first, int start_next, FOFparticle _center) {
        cellcenter = _center;
        start[0] = start_first;
        start[1] = -1;
        start[4] = start_next;
        active = 0;
        d2_furthest = 0.;
    }

    /// Return the number of particles in a given bin; input must be [0,4)
    inline int binsize(int bin) {
        // assert(bin>=0 && bin<4);   // TODO: Can hope to remove this
        return start[bin+1]-start[bin];
    }
    
    /// Compute the distance to the nearest point in the cell, from the supplied halocenter
    void compute_d2(FOFparticle *halocen) {
        FOFloat distx = fabs(cellcenter.x-halocen->x)-GFC->FOFhalfcell;
        FOFloat disty = fabs(cellcenter.y-halocen->y)-GFC->FOFhalfcell;
        FOFloat distz = fabs(cellcenter.z-halocen->z)-GFC->FOFhalfcell;
        if (distx < 0.) {
            distx = 0;
        }
        if (disty < 0.) {
            disty = 0.;
        }
        if (distz < 0.) {
            distz = 0.;
        }
        // The minimum distance to the current halo nucleus position
        FOFloat d2 = distx*distx+disty*disty+distz*distz;
        firstbin = floor(sqrt(d2)/GFC->SOpartition);

        return;
    }
};


/** This class contains the workspace to do SO finding on a 
given contiguous set of particles.

The input particle sets are not permuted, and no ability to 
permute is supported.   The FOFparticle listing is in the SO
group order, and one uses the indexes therein to access the 
original data for those particles.
*/

class SOcell {
  public:
    FOFparticle *p;     ///< The particles
    FOFloat *density;   ///< The densities
    FOFloat *min_inv_den;  ///< The minimum inverse enclosed densities (i.e., highest enclosed density)

    int *halo_index;       ///< List of halo indices for each particle
    FOFloat *d2buffer;    ///< A buffer of distances
    FOFloat *d2_bin;      ///< A buffer of sorted distances
    int np;             ///< The number of particles in this group
    int maxsize;        ///< The maximum number of particles
    
    FOFloat threshold;  ///< The density threshold we are applying
    FOFloat xthreshold;
    FOFloat FOFunitdensity; //  cosmic mean in FOF units
    FOFloat mag_roche;    /// Condition for a satellite halo to eat up particles belonging to a larger halo
    FOFloat min_central;   ///< The minimum FOF-scale density to require for a central particle.
    FOFloat *twothirds;    ///< We compare x^(3/2) to integers so much that we'll cache integers^(2/3)
    FOFloat min_radius2;    ///< We require that R_Delta be at least 0.5*KernelRadius
    FOFloat Rdensmax2;    ///< The scale over which a particle must have the highest kernal density of all neighbors in order to be called a maximum.

    FOFloat alpha_eligible2;  ///< How to scale R_Delta when assigning eligibility.
        // TODO: Note that we must have Rdensmax2 be smaller than this radius. 

    FOFgroup *groups;   ///< The list of found groups
    int ngroups;        ///< The number of groups

    FOFloat *halo_radius2; ///< The distance to the threshold squared for each halo center B.H.
    int *center_particle;  ///< The index of the particle used in each group.

    SOcellgroup *socg;  ///< a list of the cell groups
    unsigned int *cellindex;     ///< a list for each particle in pos of which cell it belongs to
    FOFloat *d2_active; ///< A buffer of distances to the particles in active cells
    integer3 refcell;   ///< the cell index triple for the first particle, -128
  
    int maxcg;          ///< The maximum number of cellgroups
    int ncg;            ///< The active number of cellgroups
    
    FOFTimer Total;
    FOFTimer Copy, Sweep, Distance, Search, Sort; 
    long long numdists;  ///< Total number of distances we compute
    long long numsorts;  ///< Total number of sorting elements
    long long numcenters; ///< Total number of SO centers considered
    long long numcg;     ///< Total number of cell groups considered
    long long numgroups; ///< Total number of groups defined (even if smaller than minmass)

    int use_aux_dens;  ///< Densities from acc or aux?

    char pad[CACHE_LINE_SIZE];    // Just to ensure an array of these always fall on
        // a different cache line

    inline void reset(int _size) {
        int ret;
        np = ngroups = 0;
        if (_size+16<maxsize) return;    // Do nothing if we have enough space
        maxsize = _size+16;   // Oversize to limit the number of re-allocations
        assertf(maxsize<800e6, "Maxsize {:d} is too large\n", maxsize);
            // Very large sizes would overflow the way we store ints inside floats in FOFparticles
        
        if (p!=NULL) free(p);
        ret = posix_memalign((void **)&p, CACHE_LINE_SIZE, sizeof(FOFparticle)*maxsize);  assert(ret == 0);
        memset(p, 0, sizeof(FOFparticle)*maxsize);

        if (d2buffer!=NULL) free(d2buffer);
        ret = posix_memalign((void **)&d2buffer, CACHE_LINE_SIZE, sizeof(FOFloat)*maxsize);  assert(ret == 0);

        if (d2_bin!=NULL) free(d2_bin);
        ret = posix_memalign((void **)&d2_bin, CACHE_LINE_SIZE, sizeof(FOFloat)*maxsize);  assert(ret == 0);

        if (d2_active!=NULL) free(d2_active);
        ret = posix_memalign((void **)&d2_active, CACHE_LINE_SIZE, sizeof(FOFloat)*maxsize);  assert(ret == 0);

        if (groups!=NULL) free(groups);
        ret = posix_memalign((void **)&groups, CACHE_LINE_SIZE, sizeof(FOFgroup)*maxsize);  assert(ret == 0);

        if (density!=NULL) free(density);
        ret = posix_memalign((void **)&density, CACHE_LINE_SIZE, sizeof(FOFloat)*maxsize);  assert(ret == 0);

        if (min_inv_den!=NULL) free(min_inv_den);
        ret = posix_memalign((void **)&min_inv_den, CACHE_LINE_SIZE, sizeof(FOFloat)*maxsize);  assert(ret == 0);

        if (halo_radius2!=NULL) free(halo_radius2); 
        ret = posix_memalign((void **)&halo_radius2, CACHE_LINE_SIZE, sizeof(FOFloat)*maxsize);  assert(ret == 0);

        if (center_particle!=NULL) free(center_particle); 
        ret = posix_memalign((void **)&center_particle, CACHE_LINE_SIZE, sizeof(int)*maxsize);  assert(ret == 0);

        if (cellindex!=NULL) free(cellindex);
        ret = posix_memalign((void **)&cellindex, CACHE_LINE_SIZE, sizeof(unsigned int)*maxsize);  assert(ret == 0);
    
        if (halo_index!=NULL) free(halo_index);
        ret = posix_memalign((void **)&halo_index, CACHE_LINE_SIZE, sizeof(int)*maxsize);  assert(ret == 0);
    }

    void setup_socg(int size) {
        maxcg = size;
        ncg = 0;
        if (socg!=NULL) free(socg);
        int ret = posix_memalign((void **)&socg, CACHE_LINE_SIZE, sizeof(SOcellgroup)*maxcg);  assert(ret == 0);
    }
        
    void resize_socg() {
        SOcellgroup *newlist;
        int ret = posix_memalign((void **)&newlist, CACHE_LINE_SIZE, sizeof(SOcellgroup)*maxcg*2);  assert(ret == 0);
        memcpy(newlist, socg, sizeof(SOcellgroup)*maxcg);
        maxcg *= 2;   // Double the previous
        free(socg);
        socg = newlist;
        return;
    }

    SOcell() {
        p = NULL;
        d2buffer = NULL;
        d2_bin = NULL;
        d2_active = NULL;
        groups = NULL;
        density = NULL;
        min_inv_den = NULL;
        halo_radius2 = NULL; 
        center_particle = NULL;
        halo_index = NULL;

        socg = NULL;
        cellindex = NULL;
    
        maxsize = -1;
        numdists = 0;
        numsorts = 0;
        numcenters = 0;
        numcg = 0;
        numgroups = 0;
        reset(32768);
        setup_socg(2048);


        mag_roche = P.SO_RocheCoeff;    /// Condition for a satellite halo to eat up particles belonging to a larger halo

        // We need to convert the R_Delta^2
        min_radius2 = 0.25*WriteState.DensityKernelRad2*FOF_RESCALE*FOF_RESCALE;

        // We adopt the DensityKernelRadius as the local density maximum criteria
        Rdensmax2 = WriteState.DensityKernelRad2*FOF_RESCALE*FOF_RESCALE;
        assertf(Rdensmax2<GFC->SOpartition*GFC->SOpartition, "SO Local Density is bigger than SOpartition\n");   // This violates the algorithm below, which assumes that the first bin contains the local density criteria.

        alpha_eligible2 = P.SO_alpha_eligible;
        alpha_eligible2 *= alpha_eligible2;

        /*
        if (omp_get_thread_num()==0) 
            STDLOG(1,"Setting up SO with mag_roche= {:f} min_radius= {:f} Rdensmax2= {:f} alpha_eligible= {:f}\n", 
                mag_roche, sqrt(min_radius2)/FOF_RESCALE, sqrt(Rdensmax2)/FOF_RESCALE, sqrt(alpha_eligible2));
        */

        int ret = posix_memalign((void **)&twothirds, 64, sizeof(FOFloat)*(SO_CACHE+2));  assert(ret == 0);
        for (int j=0; j<SO_CACHE+2; j++) twothirds[j] = pow(j,2.0/3.0);
        // We will have terrible bugs if these aren't true!
    #ifdef AVXFOF
        assertf(sizeof(FOFparticle)==16, "FOFparticle is illegal size!");
    #endif
        assertf(sizeof(FOFgroup)%16==0, "FOFgroup is illegal size!");

        use_aux_dens = GFC->use_aux_dens;

        return;
    }

    void destroy() {
        if (p!=NULL) free(p); p = NULL;
        if (d2buffer!=NULL) free(d2buffer); d2buffer = NULL;
        if (d2_bin!=NULL) free(d2_bin); d2_bin = NULL;
        if (d2_active!=NULL) free(d2_active); d2_active = NULL;
        if (groups!=NULL) free(groups); groups = NULL;
        if (density!=NULL) free(density); density = NULL;
        if (min_inv_den!=NULL) free(min_inv_den); min_inv_den = NULL;
        if (halo_radius2!=NULL) free(halo_radius2); halo_radius2 = NULL; 
        if (center_particle!=NULL) free(center_particle); center_particle = NULL; 
        if (halo_index!=NULL) free(halo_index); halo_index = NULL;

        if (socg!=NULL) free(socg); socg = NULL;
        if (cellindex!=NULL) free(cellindex); cellindex = NULL;
    }
    ~SOcell() {
        destroy();
        free(twothirds);
    }

    // Co-add the timers in x to this one
    void coadd_timers(SOcell &x) {
        Total.increment(x.Total.get_timer());
        Copy.increment(x.Copy.get_timer());
        Sweep.increment(x.Sweep.get_timer());
        Distance.increment(x.Distance.get_timer());
        Search.increment(x.Search.get_timer());
        Sort.increment(x.Sort.get_timer());
    }

    /// These are density thresholds in interparticle units,
    /// where the cosmic mean is unity.
    ///
    /// threshold refers to the spherical overdensity we'll require.
    /// min_central refers to the minimum FOFscale density that we require
    /// to make a particle eligible to be the seed of a SO halo.
    /// One must keep min_central to be comfortably larger than the density 
    /// at the edge of larger halos, which is threshold/3 for SIS.  
    /// We input _min_central_mass, which is used to compute this.
    void setup (FOFloat _threshold, FOFloat min_central_mass) {
        // We're comparing the threshold to overdensity density:
        // mass/(4*PI/3)/r^3.  So that means r^3*threshold*4*PI/3 < count
        // in interparticle length units, which are PPD times the code units.
        // Define T = (threshold*4*PI/3)**(1/3), so that
        // we just compare (r*T)**3 to count.  The positions are in 
        // FOF_RESCALE units, so we have to multiply by PPD/FOF_RESCALE
        // to get to interparticle length units.
        threshold = _threshold;
        xthreshold = pow(threshold*4.0*M_PI/3.0*P.np, 1.0/3.0)/FOF_RESCALE;
        // But we have distance squared, so we'll square this, so that
        // multiplying this by d2 gives (r*T)**2, and ^(3/2) of that 
        // is compared to count.
        xthreshold *= xthreshold;

        // Cosmic unit density yields a count in our FOFscale densities of this:
        // The (b^2-r^2) weighted density for the cosmic mean is (2/5)(4*PI/3) b^5.
        // And that density for a SIS of a given Delta and M_Delta is 
        // (2/3)(4*PI/3) Delta R_Delta^2 b^3, so that in cosmic units it's
        // (5/3) (1/b^2) (M_Delta^2 Delta/48 pi^2)^(1/3).
        // Keeping this larger than C*Delta requires M_Delta > Delta b^3 sqrt(48 pi^2 C^3/125)
        // which is about 2 Delta b^3 C^1.5.  This is why we scale the min_central_mass
        // with Delta.
        // This is further documented at the top of the file.
        FOFunitdensity = WriteState.FOFunitdensity;
        FOFloat M_D = min_central_mass;
        FOFloat sigma3 = M_D*sqrt(threshold*P.np/(48*M_PI*M_PI));
        // Density for a SIS with mass M_D
        min_central = 5./WriteState.DensityKernelRad2*pow(sigma3,2./3); 
        min_central /= (P.np);
        min_central *= FOFunitdensity;
        return;
    }

// =================================================================
// These are still class functions, but we're moving the indentation out by 1
    
/// Consider the list d2_active[start,last) and partition it on the value r2_part
/// Return the index of the first value above r2_part 
int partition_only(FOFloat r2_part, int start, int last) {
    if (start==last) {
        return last;   // could happen if there are no particles in one of the partitions
    }
    int s = start;
    while (s<last && d2_active[s]<=r2_part) s++;
    // Advance start to the first high spot
    while (s<last) {
        last--; // Consider next element in the upper list
        if (d2_active[last]<=r2_part) {
            // We've found an out of place one; flip it with start
            std::swap(p[s],p[last]);
            std::swap(density[s],density[last]);
            std::swap(d2_active[s],d2_active[last]);
            std::swap(halo_index[s],halo_index[last]);
            std::swap(min_inv_den[s],min_inv_den[last]);
            s++;
            while (s<last && d2_active[s]<=r2_part) s++;
            // Advance s to the next high spot
        }
    }
    // We're done when start = last.  last is the first high particle
    return last;
}



/// Sets up a given cellgroup for use.  This involves computing the
/// distance to all particles and then partitioning them into 4 shells
/// on the global registration around the already declared halo center.
void partition_cellgroup(SOcellgroup *cg, FOFparticle *center) {
    // Compute the distances from all particles in group cg to 
    // supplied halo center.
    // But the AVX compute_d2() function creates misalignment, so for
    // now we copy the results into the properly indexed array d2_active.
    
    // d2_active has all dist2 to the particles in bins [0,r]
    // (with same indexing as positions *p) to be used later in the loop 
    // from which partition_cellgroup is called;

    int len = cg->start[4]-cg->start[0];
    Search.Stop();
    Distance.Start();
    FOFloat *d2 = compute_d2(center, p+cg->start[0], len, d2buffer, numdists);
    Distance.Stop();
    Search.Start();
    for (int j=0; j<len; j++) {
        d2_active[cg->start[0]+j] = d2[j];
        if (d2[j] > cg->d2_furthest) cg->d2_furthest = d2[j];
    }

    // Partition into 4 radial parts
    FOFloat r2[5];
    for (int j=1; j<4; j++) {
        r2[j] = GFC->SOpartition*(cg->firstbin+j); r2[j] *= r2[j];
    }
    // These now contain the boundaries (r2[0] is irrelevant)
    cg->start[2] = partition_only(r2[2],cg->start[0],cg->start[4]);
    cg->start[1] = partition_only(r2[1],cg->start[0],cg->start[2]);
    cg->start[3] = partition_only(r2[3],cg->start[2],cg->start[4]);
}

  
/// Searches for the density crossing in this shell, assuming a mass interior
/// to it.  Returns -1 if not found; else returns square distance of threshold.
/// This code should only be called if the outer edge of the shell, combined
/// with the mass interior to the shell, would fall below the threshold.
FOFloat partial_search(int len, int mass, FOFloat shell_max_rad2, int &size_thresh, FOFloat &inv_enc_den) {
    // number of particles within threshold in that partition
    size_thresh = 0;
    FOFloat x = 0;

    if (len==0) {
        // It is rare, but this could get called on an empty shell.
        // We'll return the answer for the outer edge of the shell.
        // shell_max_rad2 needs to be supplied in the same units as d2_bin[].
        x = shell_max_rad2*xthreshold;
        inv_enc_den = x*sqrt(x)/((size_thresh+mass)*threshold);
        return shell_max_rad2;
    }
    
    // Sort the distances in increasing order
    // std::sort(d2_bin, d2_bin+len);
    ips4o::sort(d2_bin, d2_bin+len);
    numsorts += len;
    // Now sweep in from the center to find the threshold
    for (int j=0; j<len; j++) {
        x = d2_bin[j]*xthreshold;
        size_thresh = j+1; // we want the rightmost on the left side of the density threshold //TODO: ASK
        if (x*sqrt(x)>(size_thresh+mass) && d2_bin[j]>=min_radius2) {
            break;
            // This particle is below the overdensity threshold
        }
    }
    // record result
    // assertf(size_thresh+mass>0, "Found a zero mass interior to a SO shell, len = {:d}", len);
    inv_enc_den = (x*sqrt(x))/((size_thresh+mass)*threshold);
    
    if (size_thresh==len) {
        // Didn't find a threshold crossing, so return a signal of this.
        return -1.0;
    }
    return (d2_bin[size_thresh-1]);
}
    
// ======================  Search SOcellgroup Threshold  ========================
    
/// We will search outward from the given halo center, using the cellgroups
/// and shells to order the search.  Returns the square radius of the chosen 
/// density threshold, the mass interior to that, as well as 
/// the inverse enclosed density implied.  The latter is supposed to be 
/// very close to the inverse of the density threshold, but it might differ
/// if we ran out of particles before getting down to that density.

FOFloat search_socg_thresh(FOFparticle halocenter, FOFloat halocentraldensity, int &mass, FOFloat &inv_enc_den, int &is_not_density_maximum) {
    mass = 0;
    FOFloat FOFr2 = 0;
    FOFloat x = 0;
    
    int size_bin;
    int size_thresh;
    FOFloat d2_thresh = 0.;
    int size_partition;
      
    // Compute the distance to all of the cellgroup centers,
    // and find the furthest one.
    Search.Stop();
    Distance.Start();
    int furthest_firstbin = -1;
    for (int i = 0; i<ncg; i++) {
        socg[i].active = 0;   // Reset this flag
        socg[i].compute_d2(&halocenter);
        // uses minimum distance to halo center and gives us firstbin
        if (socg[i].firstbin > furthest_firstbin) {
            furthest_firstbin = socg[i].firstbin;
        }
    }
    Distance.Stop();
    Search.Start();

    // If we find no density threshold crossing, we'll return furthest edge -- not ideal
    // Proceed outward through the shells
    for (int r = 0; r<furthest_firstbin+4; r++) {
        // Compute the outer boundary of this shell
        FOFr2 = GFC->SOpartition*(r+1);
        FOFr2 *= FOFr2;
        x = xthreshold*FOFr2; 

        // Partition the newly touched cells 
        for (int i = 0; i<ncg; i++) {
    
            if (socg[i].firstbin == r) {
                socg[i].active = 1;
                // Get start[1] thru start[3] for every new SOgroupcell
                // Compute the d2 in particle order and then use this
                // for partition -- stored in d2_active when fn is called
                partition_cellgroup(socg+i, &halocenter);
                if (r==0) {
                    // We want to check if this particle is a density maximum.
                    // This means that it must be denser than all other particles
                    // within Rdensmax.
                    // We know that all particles higher in density must already
                    // be ineligible, else they would have been the center!
                    for (int j=socg[i].start[0]; j<socg[i].start[1]; j++) {
                        if (d2_active[j]<Rdensmax2) {
                            // It's close enough to test
                            if (density[j]>halocentraldensity)
                                is_not_density_maximum = 1;  
                            else halo_index[j] = -abs(halo_index[j]);  // Make ineligible
                        }
                    }
                }
            }
        }

        if (is_not_density_maximum) {
            // We have found that this center is not acceptable.  Sound the alarm.
            return 0.0;
        }

        // Is there enough mass within to skip or not enough to look
        if (x*sqrt(x) < mass) { 
            // Enough mass, so crossing cannot be inside and we just add that mass
            for (int i = 0; i<ncg; i++) {
                int bin = r-socg[i].firstbin;
                if (bin>=0 && bin<4) {
                    // Number of particles in this cell in this partition
                    mass += socg[i].binsize(bin);
                }
            }
        }
        else { 
            // Not enough mass, so there could be a density crossing.
            // We need to gather the distances to all the particles in this radial bin
            size_bin = 0;   // Number of particles within this radial bin
            for (int i = 0; i<ncg; i++) {
                // Is this SOcellgroup crossed by the radial bin?
                int bin = r-socg[i].firstbin;
                if (bin>=0 && bin<4) {
                    size_partition = socg[i].binsize(bin);
                    // Number of particles for that radial bin in that SOcellgroup
                    for (int j=0; j<size_partition; j++) {
                        d2_bin[size_bin+j] = d2_active[socg[i].start[bin]+j];
                    }
                    size_bin += size_partition;
                }
            }
            // d2_bin is used in partial_search and modified through
            // swapping and sorting; d2_active has all FOF particles in bins 0 through
            // r partitioned and should remain intact till done with the halo center.

            // Search for density threshold in list, given previous mass.
            d2_thresh = partial_search(size_bin, mass, FOFr2, size_thresh, inv_enc_den);
            if (d2_thresh > 0.0) {
                // If something was found, record it
                mass += size_thresh;
                return d2_thresh;
            }
            else {
                // False alarm: add all mass in that radial bin
                mass += size_bin;
            }
        }
    }
    // If nothing was found return the largest distance to a particle encountered 
    if (d2_thresh <= 0.0) {
        FOFloat d2_max = 0.;
        for (int i = 0; i<ncg; i++) {            
            if (socg[i].d2_furthest > d2_max) d2_max = socg[i].d2_furthest;
        }
      
        x = xthreshold*d2_max;
        inv_enc_den = (x*sqrt(x))/(mass*threshold);
        return d2_max;
    }
    
    QUIT("search_socg_thresh should have returned by now.\n");
    // Record inverse density and threshold radius

    return 0;  // just to silence the compiler
}

  
// ======================  Partition ===========================

/// Partition halos on halo with index i in the [start,last) region.
/// Also find the index of the densest particle in the high-separation
/// region.
///
/// We also return 'size', which should be the first index of the high
/// separation region, relative to start.  This may differ from 
/// what was found externally if there is some conspiracy about 
/// floating point equality comparisons of d2.  
void partition_on_index(int halo_i, int start, int last, int &size) {
    if (start==last) {
        size = 0;
        return;    
    }

    FOFloat maxdens=-1.0;
    int densest=-1;
    int s = start;
    while (s<last && halo_index[s]==halo_i) {
        // Advance start to the first high spot
        if (density[s]>maxdens) { maxdens=density[s]; densest=s; }
        s++;
    }
    while (s<last) {
        // s is now pointed to a high entry; we need to decrement last to find a low entry to swap.
        last--; // Consider next low element in the last part of the list
        if (halo_index[last]==halo_i) {
            // We've found an out of place one; flip it with start
            std::swap(p[s],p[last]);
            std::swap(density[s],density[last]);
            std::swap(halo_index[s],halo_index[last]);
            while (s<last && halo_index[s]==halo_i) {
                if (density[s]>maxdens) { maxdens=density[s]; densest=s; }
                s++;
                // Advance s to the next high spot
            }
        }
    }
    // We're done when start = last.  last is the first high particle.
    size = last-start;
    // Put the densest particle at the starting point.
    if (size>1) {
        std::swap(p[densest],p[start]);
        std::swap(density[densest],density[start]);
        std::swap(halo_index[densest],halo_index[start]);
    }
    return;
}




// ======================  SO Algorithm  ===========================
int greedySO() {
    int start = 0;  // Index of the first particle of the current group
    int densest = -1;
    FOFloat maxdens = -1.0;

    halo_radius2[0] = 0.;  /// halo[0] is the unassigned fluff, counting starts at 1 B.H.
    int count = 1;        /// The number of the group
    
    Sweep.Start();
    // Loop over all particles in the L0 group
    for (int j=0; j<np; j++) {
        // Initialize the minimum inverse enclosed density and halo indices arrays
        min_inv_den[j] = 1.e30;
        halo_index[j] = SO_UNASSIGNED;
        // An unassigned halo has halo_index=SO_UNASSIGNED; 
        // first halo has 1; second 2, third 3, etc.
        // Negative values indicate that a particle is ineligible to be a center
        // Looks for densest particle in L0
        if (density[j]>maxdens) {
            maxdens=density[j];
            densest = j;
        }
    }
    Sweep.Stop();

    // First halo center is the densest particle in the group
    start = densest;
    
    while (start>=0) {

        // If density is insufficient, do not form a new halo.
        // But we always try at least one central.
        if (density[start]<min_central && count>1) break;
        numcenters++;
        center_particle[count] = p[start].index();

        Search.Start();
        FOFloat inv_enc_den;
        int mass;
        // Call search function and obtain inverse density and threshold distance (d2SO), 
        // but also mass within the threshold.
        // Threshold distance d2SO has the property that all particles within it
        // satisfy the overdensity condition
        int is_not_density_maximum=0;
        FOFloat d2SO = search_socg_thresh(p[start],density[start],mass,inv_enc_den,is_not_density_maximum);
        halo_radius2[count] = d2SO;
        // Threshold distance d2SO has the property that particles
        // found in it can never be eligible centers
        Search.Stop();

        int densest = -1; // index of the next densest eligible particle
        FOFloat maxdens = -1.0; // density of the next densest eligible particle

        if (is_not_density_maximum) {
            // We've discovered that this center was a bad choice -- 
            // not a local maximum in the density field.
            // Need to search again.
            Sweep.Start();
            for (int j=0; j<np; j++) {
                if (halo_index[j]>=0 && density[j]>maxdens) {
                    maxdens=density[j];
                    densest = j;
                }
            }
            start = densest;
            Sweep.Stop();
            continue;
        }
        // Otherwise, we have a good center, and we need to process all the
        // particles inside of d2SO (and find the next candidate center).
        
        // Inverse density used when interpolating
        FOFloat inv_d;
        // Inverse product of the enclosed density and threshold distance squared
        FOFloat inv_d2del = inv_enc_den/(d2SO);

        FOFloat Religible2 = alpha_eligible2*d2SO;

        Sweep.Start();

        // Here we record the particles which preliminarily are attributed to the current
        // halo center and look for the next densest particle.
        // d2_active has distances to the halo center of all particles in bins 0 through
        // r. The value of p, density and d2_active for a given particle 
        // within the threshold has the same array index.
        // Loop over all cells, ncg, and split them in two -- active and inactive cells,
        // where the active ones have their minimum distance to halo center within d2SO.
        for (int i = 0; i<ncg; i++) {
            if (socg[i].active == 0) {
                // This cell wasn't used, but it might contain the next halo center.
                // Check all the particles to see if they are eligible and densest.
                for (int j=socg[i].start[0]; j<socg[i].start[4]; j++) {
                    if (density[j]>maxdens && halo_index[j]>=0) {
                        maxdens=density[j];
                        densest = j;
                    }
                }
            } else {
                // This cell was used.
                // Make it inactive to reset the array for the next forming halo
                // socg[i].active = 0; // We're now doing this elsewhere.
                // Loop over all particles in that cell
                for (int j=socg[i].start[0]; j<socg[i].start[4]; j++) {
                    if (d2_active[j]>d2SO) {
                        // This particle is outside of the halo.
                        // Just check if it's the next center
                        if (density[j]>maxdens && halo_index[j]>=0) {
                            maxdens=density[j];
                            densest = j;
                        }
                        continue;  // Go on to the next particle
                    }
                    // This particle is within the density threshold and might
                    // be part of this halo.

                    // Further, it is ineligible to be a halo center;
                    // this is signalled by a negative halo_index.
                    // halo_index[j] = (halo_index[j]==0)?(-1):-abs(halo_index[j]);
                        // The -1 here will get overwritten just below
                    if (d2_active[j]<=Religible2)
                        halo_index[j] = -abs(halo_index[j]);

                    // Interpolate to get the density of the given particle 
                    // Dens(r)=Dens_SO d2SO/r^2
                    inv_d = d2_active[j]*(inv_d2del);
                    // If this is the highest enclosed density yet encountered
                    // for this particle by a sufficient factor mag_roche, 
                    // adopt this halo as its current choice and update max density.
                    if (min_inv_den[j] > mag_roche*inv_d) {
                        halo_index[j] = (halo_index[j]<0)?(-count):count;
                        // Preserve the sign as the indicator of eligibility.
                        min_inv_den[j] = inv_d;
                    }
                    // N.B., in the current code, I think we can only make
                    // negative halo_index[] values.
                }  
            }
        }

        Sweep.Stop();
        count++;
        start = densest;
    }
    return count;
}

/// Having assigned a halo ID to all particles, we now consolidate these into 
/// contiguous sets.  This breaks the SOcellgroup sorting.  The densest particle
/// is put as the first particle in each group.
void partition_halos(int count) {
    
    Sort.Start();
    for (int j=0; j<np; j++) {
      halo_index[j] = abs(halo_index[j]);
    }
    
    int size = 0;
    int start = 0;
    int halo_ind;
    for (int i=0; i<count; i++) {
        // TODO: Document why we're putting the unassigned particles in the second slot
        if (i == 0) halo_ind = 1; // start with first subhalo                
        else if (i == 1) halo_ind = SO_UNASSIGNED; // deal with the unassigned fluff  
        else halo_ind = i; // then rest

        // Partition so as to move left all particles in halo halo_ind.
        // This also moves the densest particle to the start of its list.
        partition_on_index(halo_ind, start, np, size);
	
        // Mark the group.  
        if (size > 0) {
            if (halo_ind < SO_UNASSIGNED) groups[ngroups++] = FOFgroup(start,size,halo_radius2[halo_ind],center_particle[halo_ind]); 
            start += size;
        }
    }
    Sort.Stop();
}

// ================ Routines for cell indexes and cell groups ================

/// We want our indices to be more local, so let's get the values for one cell.
/// And then we subtract 128, so that the delta(cell) is around 128 +- few.
inline void set_reference_cell(posstruct &p) {
    refcell.x = floor((p.x+CP->halfinvcpd)*CP->cpd)-128;
    refcell.y = floor((p.y+CP->halfinvcpd)*CP->cpd)-128;
    refcell.z = floor((p.z+CP->halfinvcpd)*CP->cpd)-128;
}

/// This provides a simple parsing of the positions back into uniquely
/// numbered cell indices.  
// NOTE: This assumes that particles occupy [-halfinvcpd,+halfinvcpd) in cells
inline unsigned int compute_cellindex(posstruct &p) {
    unsigned int i = floor((p.x+CP->halfinvcpd)*CP->cpd)-refcell.x;
    unsigned int j = floor((p.y+CP->halfinvcpd)*CP->cpd)-refcell.y;
    unsigned int k = floor((p.z+CP->halfinvcpd)*CP->cpd)-refcell.z;
    // assertf(i>=0&&i<256, "Bad cell index i={:d}", i);
    // assertf(j>=0&&j<256, "Bad cell index j={:d}", j);
    // assertf(k>=0&&k<256, "Bad cell index k={:d}", k);
    return (i<<16)|(j<<8)|k;
}

/// Given the cellindex number, return the cell center
inline FOFparticle compute_cellcenter(unsigned int cellidx) {
    int k = static_cast<int>(cellidx&0xff);
    int j = static_cast<int>((cellidx&0xff00)>>8);
    int i = static_cast<int>((cellidx&0xff0000)>>16);
    // assertf(i>=0&&i<256, "Bad cell index i={:d}", i);
    // assertf(j>=0&&j<256, "Bad cell index j={:d}", j);
    // assertf(k>=0&&k<256, "Bad cell index k={:d}", k);

    return {{CP->invcpd*(i+refcell.x),
             CP->invcpd*(j+refcell.y),
             CP->invcpd*(k+refcell.z)}, 0};
}
  
/// Given the cellindex[] array, we want to scan through
void load_socg() {
    unsigned int lastidx = 0x0f000000;   // An impossible index
    int laststart = -1;
    ncg = 0;

    for (int j=0; j<np; j++) {
        if (cellindex[j] == lastidx) continue;
            // else we've found the start of a new group
        if (laststart==-1) {
            // It's the first group, so just start it.
            lastidx = cellindex[j]; laststart = j; continue;
        }
        // else record the old group and start anew
        socg[ncg++].load(laststart, j, compute_cellcenter(lastidx));
        if (ncg>=maxcg) resize_socg();   
                // This guards against overrunning the socg memory
        lastidx = cellindex[j]; laststart = j; 
    }
    // Always have one group left at the end
    socg[ncg++].load(laststart, np, compute_cellcenter(lastidx));
    numcg += ncg;
    return;
}


// =================== Main calling routine =================================

/// This takes in all of the positions, velocities, and aux for
/// a cell of size n.

/// It copies pos[0,n) into the FOFparticle array.
/// vel and aux are not used, but acc[0,n) is required, as it carries
/// the FOF-scale densities that we use.
/// Then perform the SO.
/// Return the number of multiplet groups.
/// pos[], etc. will be unchanged.

int findgroups(posstruct *pos, velstruct *vel [[maybe_unused]], auxstruct *aux, FLOAT3p1 *acc, int n) {
    Total.Start();
    reset(n);
    np = n;

    // Load the particles
    Copy.Start();
    set_reference_cell(pos[0]);
    for (int j=0; j<np; j++) {
        p[j] = FOFparticle(pos[j],j);
        if(use_aux_dens) density[j] = aux[j].get_density();
        else density[j] = acc[j].w;
        cellindex[j] = compute_cellindex(pos[j]);
    }
    // Remove: Copy.Stop() was here

    // Load the cellgroup list from the indices
    load_socg();
    Copy.Stop();

    // TODO question: After this point, do we ever use the particle cellindex again?
    // If not, then let's not permute it. BTH I don't think it's being permuted
    // TODO: In fact, it seems we could have re-used halo_index[].
    // BTH Am in that case going to use halo_index here instead of cellindex --> should it be indbuffer
    // I think we agreed that the particle index is in fact all that is needed
    // if we want to put the L1 particles back into cellgroup order.

    int count_halos = greedySO();
    partition_halos(count_halos);
    
    // Sort the particles in each group into their original ordering,
    // as this restores the cell groups for L2.
    // TODO: Didn't add any code to skip this for L2
    // Note: I left the first particle unsorted, as it is planned to be the densest one.  This is a tiny inefficiency for L2: one extra group.


    Copy.Start();
    for (int g=0; g<ngroups; g++) {
        // std::sort(p+groups[g].start+1, p+groups[g].start+groups[g].n);
        ips4o::sort(p+groups[g].start+1, p+groups[g].start+groups[g].n);
    }
    numgroups += ngroups;
    Copy.Stop();

    Total.Stop();
    return ngroups;
}

};    // End of SOCell class


#ifdef UNUSED_CODE
// ======================  SO Hist Search  ===========================
#define SO_HIST_SEARCH 16     ///< The size to switch to the histogram based method
// Timings suggest this is a very small cross-over point.

/// We are given a list of d2[0,len).  We want to find the largest value
/// of d2 such that the count of all values <= that value satisfies
/// the overdensity criterium.  We will accelerate this with a 
/// variant of a radix sort: we don't actually need to sort the whole
/// list, but just find this one intercept.
FOFloat hist_search(FOFloat *d2use, int len) {
    #define SO_BINS 16
    FOFloat norm = SO_BINS*xthreshold*pow(len,-2.0/3.0);
    FOFloat inorm = xthreshold/norm;
    // d2*xthreshold = x, which is to be compared x^(3/2) vs count
    // d2*xthreshold*SO_BINS/len^(2/3) = d2*norm = our integer binning
    // We want bin*inorm = x, so inorm = xthreshold/norm
    // bin*inorm/xthreshold = d2 = bin/norm

    int low = 0;        /// The low-edge of the lowest bin we'll search
    int hi = SO_BINS;   /// The high-edge of the highest bin we'll search
    int size = 0;       /// The number of particles interior to the search range

    if (len > SO_HIST_SEARCH) {
        // Histogram the particles to determine low, high, and size
        // of a more efficient range.
        int bins[SO_BINS];
        for (int j=0; j<SO_BINS; j++) bins[j]=0;
        for (int j=0; j<len; j++) {
            int b = floor(d2use[j]*norm);
            if (b<SO_BINS) bins[b]++;
        }
        // Next, cumulate the bins so that bins[b] contains
        // all particles in bins <= b, i.e., all particles less than b+1
        for (int j=0; j<SO_BINS; j++) size += bins[j];

        // Now we've counted the particles in these bins.
        // If the overdensity criteria is satisfied at the outer
        // boundary of the bin, then we know those particles are in.

        // If the criteria could be satisfied by placing the mass in the
        // bin at the inner boundary, then the boundary could be in that bin.
        // Unfortunately, there is guarantee that such bins are consecutive.

        // If the criteria can't be satisifed even placing all of the mass
        // exterior to that bin at the inner boundary, then we know the boundary
        // is not exterior to that point.  So let's bound the region
        // before committing resources to a sort.
        // We're going to consider whether bin 'low' is all inside
        for (low = SO_BINS-1;low>=0;low--) {
            FOFloat x = (low)*inorm;   // The inner boundary
            if (hi==SO_BINS && x*sqrt(x)<=size) hi = low;   
            // Putting all of this bin's mass at the inner edge,
            // we could possibly get above the density threshold.  
            // So the boundary is interior to (hi+1)*inorm.
            // Can only set this the first (outermost) time
            x = (low+1)*inorm;   // The outer boundary
            if (x*sqrt(x)<=size) break;
            // Even putting all of this bin's mass at the outer edge, 
            // this bin and its interior exceeds the threshold.
            size -= bins[low];
        }
        low++; hi++; 
        // Now low*inorm is the lower bound on the boundary
        // We need to search the range [low,hi]*inorm
        // size contains all of the mass interior to low.
    }

    // We're going to load up particles in range [low,hi]*inorm into 
    // the sorting array starting at size+1
    int start = size;    
    for (int j=0; j<len; j++) {
        FOFloat b = d2use[j]*norm;
        if (b>=low && b<hi) d2_bin[++size] = d2use[j];
        // Need the >= here because if low==0, we need the d2=b=0 self-count
    }
    // Now the particles to be considered are in (start,size]
    // Sort this smaller segment
    std::sort(d2_bin+start+1, d2_bin+size+1);
    numsorts += size-start;

    // Now sweep in from the end to find the threshold
    int max = std::max(start,SO_CACHE);
    for (; size>max; size--) {
        FOFloat x = d2_bin[size]*xthreshold;
        if (x*sqrt(x)<size) goto FoundSORadius;
        // This particle exceeds the overdensity threshold
    }
    for (; size>start; size--) {
        FOFloat x = d2_bin[size]*xthreshold;
        if (x<twothirds[size]) goto FoundSORadius;
        // This particle exceeds the overdensity threshold
        // We have cached the size^(2/3) values, just to avoid
        // the square roots
    }

    FoundSORadius:
    if (size==start) return low/norm;
        // We don't have the particle from the next lower bin, 
        // so just return the boundary
    return d2_bin[size];
}
#endif
