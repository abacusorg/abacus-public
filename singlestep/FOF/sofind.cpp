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
    SOcellgroup() { } 
    ~SOcellgroup() { }

    void load(int start_first, int start_next, FOFparticle _center) {
        cellcenter = _center;
        start[0] = start_first;
        start[1] = -1;
        start[4] = start_next;
        active = 0;
        d2_furthest = 0.;
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
    FOFloat *min_inv_den;  ///< The min densities

    int *halo_inds;       ///< List of halo indices for each particle
    FOFloat *d2buffer;    ///< A buffer of distances
    FOFloat *d2_bin;      ///< A buffer of sorted distances
    int np;             ///< The number of particles in this group
    int maxsize;        ///< The maximum number of particles
    
    FOFloat threshold;  ///< The density threshold we are applying
    FOFloat xthreshold;
    FOFloat FOFunitdensity; //  cosmic mean in FOF units
    FOFloat mag_roche = 2.;    /// Condition for a satellite halo to eat up particles belonging to a larger halo
    FOFloat min_central;   ///< The minimum FOF-scale density to require for a central particle.
    FOFloat *twothirds;    ///< We compare x^(3/2) to integers so much that we'll cache integers^(2/3)

    FOFgroup *groups;   ///< The list of found groups
    int ngroups;        ///< The number of groups

    FOFloat *halo_thresh2; ///< The distance to the threshold squared for each halo center B.H.

    SOcellgroup *socg;  ///< a list of the cell groups
    int *cellindex;     ///< a list for each particle in pos of which cell it belongs to
    FOFloat *d2_active; ///< A buffer of distances to the particles in active cells
    integer3 refcell;   ///< the cell index triple for the first particle, -128
  
    int maxcg;          ///< The maximum number of cellgroups
    int ncg;            ///< The active number of cellgroups
    
    FOFTimer Total;
    FOFTimer Copy, Sweep, Distance, Search; 
    long long numdists;  ///< Total number of distances we compute
    long long numsorts;  ///< Total number of sorting elements
    long long numcenters; ///< Total number of SO centers considered


    char pad[CACHE_LINE_SIZE];    // Just to ensure an array of these always fall on
        // a different cache line

    inline void reset(int _size) {
        int ret;
        np = ngroups = 0;
        if (_size+16<maxsize) return;    // Do nothing if we have enough space
        maxsize = _size+16;   // Oversize to limit the number of re-allocations
        assertf(maxsize<8e6, "Maxsize %d is too large\n", maxsize);
            // This is required because FOFparticle is storing the
            // indices in a float.  We could in principle work around
            // this, by using additional bits in a negative exponent.
            // But current Abacus plans do not produce cells this big.
        
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

        if (halo_thresh2!=NULL) free(halo_thresh2); //B.H.
        ret = posix_memalign((void **)&halo_thresh2, CACHE_LINE_SIZE, sizeof(FOFloat)*maxsize);  assert(ret == 0);

        if (cellindex!=NULL) free(cellindex);
        ret = posix_memalign((void **)&cellindex, CACHE_LINE_SIZE, sizeof(int)*maxsize);  assert(ret == 0);
    
        if (halo_inds!=NULL) free(halo_inds);
        ret = posix_memalign((void **)&halo_inds, CACHE_LINE_SIZE, sizeof(int)*maxsize);  assert(ret == 0);
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
        halo_thresh2 = NULL; // B.H.
        halo_inds = NULL;

        socg = NULL;
        cellindex = NULL;
    
        maxsize = -1;
        numdists = 0;
        numsorts = 0;
        numcenters = 0;
        reset(32768);
        setup_socg(2048);

        int ret = posix_memalign((void **)&twothirds, 64, sizeof(FOFloat)*(SO_CACHE+2));  assert(ret == 0);
        for (int j=0; j<SO_CACHE+2; j++) twothirds[j] = pow(j,2.0/3.0);
        // We will have terrible bugs if these aren't true!
    #ifdef AVXFOF
        assertf(sizeof(FOFparticle)==16, "FOFparticle is illegal size!");
    #endif
        assertf(sizeof(FOFgroup)%16==0, "FOFgroup is illegal size!");
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
        if (halo_thresh2!=NULL) free(halo_thresh2); halo_thresh2 = NULL; // B.H.
        if (halo_inds!=NULL) free(halo_inds); halo_inds = NULL;

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
    }

    /// These are density thresholds in interparticle units,
    /// where the cosmic mean is unity.
    ///
    /// threshold refers to the spherical overdensity we'll require.
    /// min_central refers to the minimum FOFscale density that we require
    /// to make a particle eligible to be the seed of a SO halo.
    /// One probably wants min_central = threshold/3 or so.
    void setup (FOFloat _threshold, FOFloat _min_central) {
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
        // TODO: Document better
        FOFunitdensity = P.np*4.0*M_PI*2.0/15.0*pow(WriteState.DensityKernelRad2,2.5);
        FOFloat M_D = 35.;
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
            std::swap(halo_inds[s],halo_inds[last]);
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
    FOFloat *d2 = compute_d2(center, p+cg->start[0], len, d2buffer, numdists);
    for (int j=0; j<len; j++) {
        d2_active[cg->start[0]+j] = d2[j];
        if (d2[j] > cg->d2_furthest) cg->d2_furthest = d2[j];
    }

    // Partition into 4 radial parts
    FOFloat r2[5];
    for (int j=1; j<=4; j++) {
        r2[j] = GFC->SOpartition*(cg->firstbin+j); r2[j] *= r2[j];
    }
    // These now contain the boundaries (r2[0] is irrelevant)
    cg->start[2] = partition_only(r2[2],cg->start[0],cg->start[4]);
    cg->start[1] = partition_only(r2[1],cg->start[0],cg->start[2]);
    cg->start[3] = partition_only(r2[3],cg->start[2],cg->start[4]);
}

  
/// Searches for the density crossing in this shell, assuming a mass interior
/// to it.  Returns -1 if not found; else returns square distance of threshold.
FOFloat partial_search(int len, int mass, int &size_thresh, FOFloat &inv_enc_den) {
    // number of particles within threshold in that partition
    size_thresh = 0;
    
    // Sort the distances in increasing order
    std::sort(d2_bin, d2_bin+len); 
    numsorts += len;
    // Now sweep in from the center to find the threshold
    FOFloat x;
    for (int j=0; j<len; j++) {
        x = d2_bin[j]*xthreshold;
        size_thresh = j+1; // we want the rightmost on the left side of the density threshold //TODO: ASK
        if (x*sqrt(x)>(size_thresh+mass) && d2_bin[j]>=.25*WriteState.DensityKernelRad2) {
            break;
            // This particle is below the overdensity threshold
        }
    }
    // record result
    assertf(size_thresh+mass>0, "Found a zero mass interior to a SO shell, len = %d", len);
    inv_enc_den = (x*sqrt(x))/((size_thresh+mass)*threshold);
    
    if (size_thresh==len) {
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

FOFloat search_socg_thresh(FOFparticle *halocenter, int &mass, FOFloat &inv_enc_den) {
    mass = 0;
    FOFloat FOFr2 = 0;
    FOFloat x = 0;
    
    int size_bin;
    int size_thresh;
    FOFloat d2_thresh;
    FOFloat d2_max = 0.;
    int size_partition;
      
    // Compute the distance to all of the cellgroup centers,
    // and find the furthest one.
    int furthest_firstbin = -1;
    for (int i = 0; i<ncg; i++) {
        socg[i].compute_d2(halocenter);
        // uses minimum distance to halo center and gives us firstbin
        if (socg[i].firstbin > furthest_firstbin) {
            furthest_firstbin = socg[i].firstbin;
        }
    }

    // If we find no density threshold crossing, we'll return furthest edge -- not ideal
    FOFloat r2found; // TODO: REMOVE
    //r2found = (furthest_firstbin)*GFC->SOpartition;
    //r2found *= r2found;
    //inv_enc_den = 1./threshold;//1.e30; //TESTING no defined inv
    
    // Proceed outward through the shells
    for (int r = 0; r<furthest_firstbin+4; r++) {
        // Compute the inner boundary of this shell
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
                partition_cellgroup(socg+i, halocenter);
            }
        }

        // Is there enough mass within to skip or not enough to look
        if (x*sqrt(x) < mass) { 
            // Enough mass, so crossing cannot be inside and we just add that mass
            for (int i = 0; i<ncg; i++) {
                if (socg[i].firstbin >= r-3 && socg[i].firstbin <= r) {
                    // Number of particles in this cell in this partition
                    mass += socg[i].start[r-socg[i].firstbin+1]-socg[i].start[r-socg[i].firstbin];
                }
            }
        }
        else { 
            // Not enough mass, so there could be a density crossing.
            // Number of all particles that are within this radial bin
            size_bin = 0;
            // for every SOcellgroup crossed by the radial bin
            for (int i = 0; i<ncg; i++) {
                if (socg[i].firstbin >= r-3 && socg[i].firstbin <= r) {
                    // Number of particles for that radial bin in that SOcellgroup
                    size_partition = socg[i].start[r-socg[i].firstbin+1]-socg[i].start[r-socg[i].firstbin];
                    for (int j = 0; j<size_partition; j++) {
                        d2_bin[size_bin+j] = d2_active[socg[i].start[r-socg[i].firstbin]+j];
                    }
                    size_bin+=size_partition;
                }
            }
            // d2_bin is used in partial_search and modified through
            // swapping and sorting; d2_active has all FOF particles in bins 0 through
            // r partitioned and should remain intact till done with the halo center
            // However, we need an array d2_bin
            // to save only the particle distances in r
            
            // Search for density threshold in list, given previous mass.
            Distance.Start();
            d2_thresh = partial_search(size_bin, mass, size_thresh, inv_enc_den);
        
            Distance.Stop();
            if (d2_thresh > 0.0) {
                // If something was found, record it
                mass += size_thresh;

                //r2found = d2_thresh;
                return d2_thresh;
            }
            else {
                // False alarm: add all mass in that radial bin
                mass += size_bin;

                //continue;
            }
        }
    }
    // If nothing was found return the largest distance to a particle encountered 
    if (d2_thresh <= 0.0) {
        for (int i = 0; i<ncg; i++) {            
            if (socg[i].d2_furthest > d2_max) d2_max = socg[i].d2_furthest;
        }
      
        x = xthreshold*d2_max;
        inv_enc_den = (x*sqrt(x))/(mass*threshold);
        return d2_max;
    }
    printf("it should never get here\n");
    // Record inverse density and threshold radius
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
int partition_and_index(int *halos, int halo_i, int start, int last, int &size) {
  
    FOFloat maxdens=-1.0;
    int densest=-1;
    if (start==last) {
        size = 0;
        return -1;    // Nothing to do
    }

    int halo_ip; // next halo
    if (halo_i == 0) halo_ip = halo_i+2; // unassigned particle is after 1st halo
    else if (halo_i == 1) halo_ip = halo_i-1; // first halo
    else halo_ip = halo_i+1;
    int s = start;
    while (s<last && halos[s]==halo_i) s++;
    // Advance start to the first high spot
    while (s<last) {
        // BTH Is it the wording of next element in upper list confusing? It should just mean
        // on the other side of the array -- otherwise it does look for the last low spot
        last--; // Consider next low element in the last part of the list
        if (halos[last]==halo_i) {
            // We've found an out of place one; flip it with start
            std::swap(p[s],p[last]);
            std::swap(density[s],density[last]);
            std::swap(halos[s],halos[last]);
            s++;
            while (s<last && halos[s]==halo_i) s++;
            // Advance s to the next high spot
        }
        if (density[last]>maxdens && halos[last]==halo_ip) {
          maxdens = density[last];
          densest = last;
        }
        // last points at a high particle.  Check if it's the densest
    }
    // We're done when start = last.  last is the first high particle
    size = last-start;
    return densest;
}




// ======================  SO Algorithm  ===========================
int greedySO() {
    int start = 0;  // Index of the first particle of the current group
    int densest = -1;
    FOFloat maxdens = -1.0;

    halo_thresh2[0] = 0.;  /// halo[0] is the unassigned fluff, counting starts at 1 B.H.
    int count = 1;        /// The number of the group
    
    Sweep.Start();
    // Loop over all particles in the L0 group
    for (int j=0; j<np; j++) {
        // Initialize the minimum inverse enclosed density and halo indices arrays
        min_inv_den[j] = 1.e30;
        halo_inds[j] = 0;
        // An unassigned halo has halo_inds=0; first halo has 1; second 2, third 3, etc.
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

        Search.Start();
        FOFloat inv_enc_den;
        int mass;
        // Call search function and obtain inverse density and threshold distance (d2SO), 
        // but also mass within the threshold.
        // Threshold distance d2SO has the property that all particles within it
        // satisfy the overdensity condition
        FOFloat d2SO = search_socg_thresh(p+start,mass,inv_enc_den);
        halo_thresh2[count] = d2SO;
        // Threshold distance d2SO has the property that particles
        // found in it can never be eligible centers
        
        Search.Stop();
        
        // Inverse density used when interpolating
        FOFloat inv_d;
        // Inverse product of the enclosed density and threshold distance squared
        FOFloat inv_d2del = inv_enc_den/(d2SO);

        int densest = -1; // index of the next densest eligible particle
        FOFloat maxdens = -1.0; // density of the next densest eligible particle

        Sweep.Start();

        // Here we record the particles which preliminarily are attributed to the current
        // halo center and look for the next densest particle.
        // d2_active has distances to the halo center of all particles in bins 0 through
        // r. The value of p, density and d2_active for a given particle 
        // within the threshold has the same array index.
        // Loop over all cells, ncg, and split them in two -- active and inactive cells,
        // where the active ones have their minimum distance to halo center within d2SO.
        for (int i = 0; i<ncg; i++) {
            if (socg[i].active == 1) {
                // If a cell is active, make it inactive so that you reset the array
                // for the next forming halo
                socg[i].active = 0;
		// Loop over all particles in that cell
                for (int j=socg[i].start[0]; j<socg[i].start[4]; j++) {
            // for those within the threshold distance, set their halo_inds to a negative number to
		    // make them ineligible to be halo centers
                    if (d2_active[j]<d2SO) {
		        // If this particle has not yet been assigned to a halo, set it to -1
		        // as a placeholder (will be changed within loop, see ***). If already assigned,
		        // make sure the halo_index remains negative (i.e. particle ineligible
		        // to ever be a halo center)
                        halo_inds[j] = (halo_inds[j]==0)?(-1):-abs(halo_inds[j]);
                    }
                    // For those outside the threshold distance which are
		    // still active (i.e. halo_inds is non-negative) and eligible, i.e. satisfies the
		    // local density to max enclosed density ratio criterion,
		    // check whether any can be the next densest particle.
                    else if (density[j]>maxdens && halo_inds[j]>=0) {
                        maxdens=density[j];
                        densest = j;
			// If the particle is outside the density threshold, we don't need to look
			// further into it as it cannot be assigned to this halo
                        if (d2_active[j]>d2SO) continue;
                    }
		    // These are the particles in this cell which are not eligible to be halo center
		    // and also outside the threshold distance. They cannot be assigned to this halo.
                    else if (d2_active[j]>d2SO)  continue;

		    // This snippet is for ALL the particles within the density threshold.
                    // Interpolate to get the density of the given particle Dens(r)=Dens_SO d2SO/r^2
                    inv_d = d2_active[j]*(inv_d2del);
                    // If j is the densest particle seen so far, i.e. its enclosed density
		    // with respect to this halo center is the largest yet, mark this as its halo
                    if (min_inv_den[j] > inv_d) {
                        // Update the max dens for that particle
                      if (min_inv_den[j] > mag_roche*inv_d) {
                            halo_inds[j] = (halo_inds[j]<0)?(-count):count;
                      }
                        min_inv_den[j] = inv_d;
                        // If this particle has already been marked as ineligible (i.e. has negative
			// halo_inds), preserve the sign and just change its halo assignment (*** notice
			// that the halo_inds in this case will definitely change within this loop)
			// The halo index gets updated only if the enclosed density of the particle
            // with respect to the newcomer is  mag_roche times its largest enclosed density so far
                    }
                }
            }
            // Loop over the remaining particles in that cell outside the threshold distance 
            // which are eligible to be halo centers to check whether any of them can be the next densest
            else {
                for (int j=socg[i].start[0]; j<socg[i].start[4]; j++) {
                    if (density[j]>maxdens && halo_inds[j]>=0) {
                        maxdens=density[j];
                        densest = j;
                    }
                }
            }
        }

        Sweep.Stop();
        count++;

        start = densest;

    }
    return count;
}

void partition_halos(int count) {
    
    for (int j=0; j<np; j++) {
      halo_inds[j] = abs(halo_inds[j]);
    }
    
    
    Sweep.Start();
    int size = 0;
    int start = 0;
    int next_densest;
    int halo_ind;
    for (int i=0; i<count; i++) {
        if (i == 0) halo_ind = 1; // start with first subhalo                
        else if (i == 1) halo_ind = 0; // deal with the unassigned fluff     
        else halo_ind = i; // then rest
      
        // rearrange so that you first sort to the left all particles in halo 0,
        // then all in halo 1, ... i, ... count;
        // and finding the densest particle of those to the right of i while doing so
        next_densest = partition_and_index(halo_inds, halo_ind, start, np, size);
        
        // TODO: it would be mildly more cache friendly to move the 
        // sorting by particle ID number here.
        
        // Mark the group.  Note that the densest particle may not be first.
        // TODO: Will write more documentation 
        if (next_densest < 0) {
            numcenters++;
            if (halo_ind > 0) groups[ngroups++] = FOFgroup(start,size,halo_thresh2[i]); // B.H.
            start += size;
            continue;
        }
        if (size > 0) {
            numcenters++;
            if (halo_ind > 0) groups[ngroups++] = FOFgroup(start,size,halo_thresh2[i]); // B.H.
            start += size;
        }
      
        // Swap so that remaining particles start with the densest.
        std::swap(p[start], p[next_densest]);
        std::swap(density[start], density[next_densest]);
        std::swap(halo_inds[start], halo_inds[next_densest]);
    }
    Sweep.Stop();
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
inline int compute_cellindex(posstruct &p) {
    int i = floor((p.x+CP->halfinvcpd)*CP->cpd)-refcell.x;
    int j = floor((p.y+CP->halfinvcpd)*CP->cpd)-refcell.y;
    int k = floor((p.z+CP->halfinvcpd)*CP->cpd)-refcell.z;
    assertf(i>=0&&i<256, "Bad cell index i=%d", i);
    assertf(j>=0&&j<256, "Bad cell index j=%d", j);
    assertf(k>=0&&k<256, "Bad cell index k=%d", k);
    return (i<<16)|(j<<8)|k;
}

/// Given the cellindex number, return the cell center
inline FOFparticle compute_cellcenter(int cellidx) {
    int k = (cellidx&0xff);
    int j = (cellidx&0xff00)>>8;
    int i = (cellidx&0xff0000)>>16;
    assertf(i>=0&&i<256, "Bad cell index i=%d", i);
    assertf(j>=0&&j<256, "Bad cell index j=%d", j);
    assertf(k>=0&&k<256, "Bad cell index k=%d", k);
    posstruct p;
    p.z = CP->invcpd*(k+refcell.z);
    p.y = CP->invcpd*(j+refcell.y);
    p.x = CP->invcpd*(i+refcell.x);
    return FOFparticle(p,0);
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

int findgroups(posstruct *pos, velstruct *vel, auxstruct *aux, FLOAT3p1 *acc, int n) {
    Total.Start();
    reset(n);
    np = n;

    // Load the particles
    Copy.Start();
    set_reference_cell(pos[0]);
    for (int j=0; j<np; j++) {
        p[j] = FOFparticle(pos[j],j);
        density[j] = acc[j].w; 
        cellindex[j] = compute_cellindex(pos[j]);
    }
    // Remove: Copy.Stop() was here

    // Load the cellgroup list from the indices
    load_socg();
    Copy.Stop();

    // TODO question: After this point, do we ever use the particle cellindex again?
    // If not, then let's not permute it. BTH I don't think it's being permuted
    // TODO: In fact, it seems we could have re-used halo_inds[].
    // BTH Am in that case going to use halo_inds here instead of cellindex --> should it be indbuffer
    // I think we agreed that the particle index is in fact all that is needed
    // if we want to put the L1 particles back into cellgroup order.

    int count_halos = greedySO();
    partition_halos(count_halos);
    
    // Sort the particles in each group into their original ordering,
    // as this restores the cell groups for L2.
    // TODO: Didn't add any code to skip this for L2
    // Note: I left the first particle unsorted, as it is planned to be the densest one.  This is a tiny inefficiency for L2: one extra group.


    for (int g=0; g<ngroups; g++) {
        std::sort(p+groups[g].start+1, p+groups[g].start+groups[g].n);
    }

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
