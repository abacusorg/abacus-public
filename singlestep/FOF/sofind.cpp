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

// B.H. new

class alignas(16) SOcellgroup {
  public:
    FOFparticle center; // The center of the cell in GlobalGroup coords 
    FOFLOAT d2;         // The minimum distance to the current halo nucleus position
    int firstbin;       // The radial bin this cell first appears in.
    int start[5];       // The starting point of particles in the L0 list
    // We will allow 4 partitions of this list: 
    // [ start[j], start[j+1] ) for j=0,1,2,3.
    // So start[4] is start[0]+np_in_cellgroup

    // constants, maybe better location for these
    FOFloat FOFhalfcell = FOF_RESCALE/2.0*CP->invcpd;     // The half cell size 
    FOFloat SOpartition = FOFhalfcell*2.0*sqrt(3.0)/3.0;  // The radial binning
    
    // Null constructors & destructors
    SOcellgroup() { } 
    ~SOcellgroup() { }

    void load(int start_first, int start_next, FOFparticle _center) {
        center = _center;
        start[0] = start_first;
        start[1] = -1;
        start[4] = start_next; 
    }

    /// Compute the distance to the nearest point in the cell, from the supplied halocenter
    void compute_d2(FOFparticle *halocen) {
        dx = fabs(center.x-halocen->x)-FOFhalfcell;
        dy = fabs(center.y-halocen->y)-FOFhalfcell;
        dz = fabs(center.z-halocen->z)-FOFhalfcell;
        if (dx<0) dx = 0.0;
        if (dy<0) dy = 0.0;
        if (dz<0) dz = 0.0;

        d2 = dx*dx+dy*dy+dz*dz;
        firstbin = floor(sqrt(d2)/SOpartition);
        return;
    }

    /// TODO: Provide a sorting command, based on d2
    // B.H. I think this just means define a sorting function
    // DJE: if we define a < operator, then the standard sort will just work.
    // But we may opt not to sort this list at all.
    // B.H. note to me: refrain from this for now
};


// B.H. end new

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
    FOFloat *density;	///< The densities
    FOFloat *min_inv_den_part;	///< The min densities

    int *halo_part;	///< The min densities
    FOFloat *d2buffer;    ///< A buffer of distances
    FOFloat *d2sort;      ///< A buffer of sorted distances
    int np;             ///< The number of particles in this group
    int maxsize;        ///< The maximum number of particles
    
    FOFloat threshold;  ///< The density threshold we are applying
    FOFloat xthreshold;
    FOFloat mag_loc = 2.;    /// Condition for creating a center
    FOFloat inner_rad2 = .8*.8; /// What is the inner radius of Delta prime with respect to Delta
    FOFloat min_central;   ///< The minimum FOF-scale density to require for a central particle.
    FOFloat *twothirds;	  ///< We compare x^(3/2) to integers so much that we'll cache integers^(2/3)

    FOFgroup *groups;   ///< The list of found groups
    int ngroups;        ///< The number of groups


    // B.H. new
    SOcellgroup *socg; ///< a list of the cell groups
    int *cellindex; ///< a list for each particle in pos of which cell it belongs to
    FOFloat *d2comb;      ///< A buffer of combined distances
  
    int maxcg;          ///< The maximum number of cellgroups
    int ncg;            ///< The active number of cellgroups
    // B.H. end new
    
    FOFTimer Total;
    FOFTimer Copy, Sweep, Distance, Search; 
    long long numdists;	///< Total number of distances we compute
    long long numsorts;	///< Total number of sorting elements
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
        // printf("Allocating FOFCell to maxsize = %d\n", maxsize);
        if (p!=NULL) free(p);
        ret = posix_memalign((void **)&p, CACHE_LINE_SIZE, sizeof(FOFparticle)*maxsize);  assert(ret == 0);
        memset(p, 0, sizeof(FOFparticle)*maxsize);
        if (d2buffer!=NULL) free(d2buffer);
        ret = posix_memalign((void **)&d2buffer, CACHE_LINE_SIZE, sizeof(FOFloat)*maxsize);  assert(ret == 0);
        if (d2sort!=NULL) free(d2sort);
        ret = posix_memalign((void **)&d2sort, CACHE_LINE_SIZE, sizeof(FOFloat)*maxsize);  assert(ret == 0);
	if (d2comb!=NULL) free(d2comb);
	ret = posix_memalign((void **)&d2comb, CACHE_LINE_SIZE, sizeof(FOFloat)*maxsize);  assert(ret == 0);
        if (groups!=NULL) free(groups);
        ret = posix_memalign((void **)&groups, CACHE_LINE_SIZE, sizeof(FOFgroup)*maxsize);  assert(ret == 0);
        if (density!=NULL) free(density);
        ret = posix_memalign((void **)&density, CACHE_LINE_SIZE, sizeof(FOFloat)*maxsize);  assert(ret == 0);
        if (min_inv_den_part!=NULL) free(min_inv_den_part);
        ret = posix_memalign((void **)&min_inv_den_part, CACHE_LINE_SIZE, sizeof(FOFloat)*maxsize);  assert(ret == 0);

        // B.H. new
        if (cellindex!=NULL) free(cellindex);
        ret = posix_memalign((void **)&cellindex, CACHE_LINE_SIZE, sizeof(int)*maxsize);  assert(ret == 0);
        // B.H. end new
	
        if (halo_part!=NULL) free(halo_part);
        ret = posix_memalign((void **)&halo_part, CACHE_LINE_SIZE, sizeof(int)*maxsize);  assert(ret == 0);
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
        d2sort = NULL;
	d2comb = NULL;
        groups = NULL;
        density = NULL;
        min_inv_den_part = NULL;
        halo_part = NULL;

        // B.H. new 
        socg = NULL;
        cellindex = NULL;
        // B.H. end new
	
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
        if (d2sort!=NULL) free(d2sort); d2sort = NULL;
	if (d2comb!=NULL) free(d2comb); d2comb = NULL;
        if (groups!=NULL) free(groups); groups = NULL;
        if (density!=NULL) free(density); density = NULL;
        if (min_inv_den_part!=NULL) free(min_inv_den_part); min_inv_den_part = NULL;
        if (halo_part!=NULL) free(halo_part); halo_part = NULL;

        // B.H. new 
        if (socg!=NULL) free(socg); socg = NULL;
        if (cellindex!=NULL) free(cellindex); cellindex = NULL;
        // B.H. end new
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
        // tuk REMOVE OR DOCUMENT b.h.
        // Cosmic unit density yields a count in our FOFscale densities of this:
        FOFloat FOFunitdensity = P.np*4.0*M_PI*2.0/15.0*pow(WriteState.DensityKernelRad2,2.5)+1e-30;
        FOFloat M_D = 30.;
        FOFloat alpha_safety = 1.;
        FOFloat sigma3 = M_D*sqrt(threshold*P.np/(48*M_PI*M_PI));
        min_central = alpha_safety*5./WriteState.DensityKernelRad2*pow(sigma3,2./3); 
        min_central /= (P.np);
        min_central *= FOFunitdensity;
        return;
    }

    
    int partition_only(FOFloat *d2_part, FOFloat r2_part, 
    	int start, int last) {

	if (start==last) {
	    return -1;    // Nothing to do, but this should never happen
	}
	int s = start;
        while (s<last && d2_part[s]<=r2_part) s++;
		// Advance start to the first high spot
	while (s<last) {
	    last--; // Consider next element in the upper list
	    if (d2_part[last]<=r2_part) {
		// We've found an out of place one; flip it with start
		std::swap(p[s],p[last]);
		std::swap(density[s],density[last]);
		std::swap(d2_part[s],d2_part[last]);
		std::swap(halo_part[s],halo_part[last]); // B.H. depending on when you call this
		s++;
		while (s<last && d2_part[s]<=r2_part) s++;
		// Advance s to the next high spot
	    }
	}
	// We're done when start = last.  last is the first high particle
	return last;
    }
    
    void partition_cellgroup(SOcellgroup *cg, FOFparticle &center) {
        // Compute the distances from all particles in group cg to 
        // supplied halo center.
        // TODO: Ideally these would go into the d2[] list with the 
        // same indexing as the p list.  
        // TODO for DJE: But we need to watch for misalignment
        // coming from the AVX compute_d2() function.  May need to alter that code.
        //compute_d2(center, p+cg->start[0], cg->start[4]-cg->start[0], d2buffer, numdists);
        // B.H. might be wrong but attempt need to def remove the other one
        FOFloat *d2 = compute_d2(center, p+cg->start[0], cg->start[4]-cg->start[0], d2buffer, numdists);

	
        // Partition into 4 radial parts
        FOFloat r2[4];
        for (int j=1; j<4; j++) {
            r2[j] = SOpartition*(cg->firstbin+j); r2[j] *= r2[j];
        }
        // These now contain the boundaries (r[0] is irrelevant)
	cg.start[2] = partition_only(d2,r2[2],p+cg->start[0],p+cg->start[4]);
	cg.start[1] = partition_only(d2,r2[1],p+cg->start[0],p+cg->start[2]);
	cg.start[3] = partition_only(d2,r2[3],p+cg->start[2],p+cg->start[4]);
    }


    // B.H. new

    FOFloat partial_search(FOFloat *d2use, int len, int m, int &partial, FOFloat &inv_enc_den) {
	int size = 0;	// We're going to unit-index d2sort[]
	for (int j=0; j<len; j++) 
	  d2sort[++size] = d2use[j]; // maybe optimize tuk
	// Now we've filled d2sort in the range [1,size] with 
	// all the points.
	
	// Sort the distances in increasing order
	std::sort(d2sort+1, d2sort+size+1); // d2sort starts at 1, not 0
	numsorts += size;
	
	// Now sweep in from the center to find the threshold
	partial = 1;
	FOFloat x;
	for (int j=1; j<=len; j++) {
	  x = d2sort[j]*xthreshold;
	  partial = j; // we want the rightmost on the left side of the density threshold
	  if (x*sqrt(x)>(j+m) && d2sort[partial]>=.5*WriteState.DensityKernelRad2) {
	    break;
	    // This particle exceeds the overdensity threshold
	  }
	}
	
	inv_enc_den = (x*sqrt(x))/((partial+m)*threshold);
	if (size==len) return -1;
	return (d2sort[partial]); 
    }


    // ======================  Fill SOcellgroup array  ===========================
    // for a given halo center fill an array of SOcell groups
    void fill_socg(FOFparticle hcenter, SOcellgroup *socg, int *cellindex) {
        // B.H. new 
        // please do more elegantly, Boryana
        // DJE: I wrote a version of this below, but it's rather similar
        int count_cg = 0;
        int start = 0;
        socg[count_cg]->center = hcenter; // center of halo check syntax in fof_sublink
        socg[count_cg]->start_first = start; // where first particle starts for this SOcellgroup
        count_cg++;
        // Then scan through list, looking for changes in cellindex to trigger new SOcellgroup
        for (int j=1; j<np; j++) {
            if (cellindex[j]-cellindex[j-1] != 0) {
                start = j;
                socg[count_cg-1]->start_next = start; // where next SOcellgroup starts

                socg[count_cg]->center = hcenter; // center of halo check syntax
                socg[count_cg]->start_first = start; // where first particle starts for this SOcellgroup
                count_cg++;
            }
        }
    }
    


    // ======================  Search SOcellgroup Threshold  ===========================
    
    FOFloat search_socg_thresh(FOFparticle hcenter, SOcellgroup *socg, FOFloat &inv_enc_den) {
        int mass = 0;
	FOFloat FOF_CONV = FOF_RESCALE*r*CP->invcpd; // FOF_RESCALE/cpd gives same units as pcle pos 
	FOFloat FOFr2 = 0;
	FOFloat x = 0;

	// Find furthest first bin
	int furthest_firstbin = -1;
	for (int i = 0; i<ncg; i++) {
	    if (socg[i]->firstbin > furthest_firstbin) {
	        furthest_firstbin = socg[i]->firstbin;
	    }
	}

	// If nothing is found, return furthest edge -- not ideal
	FOFloat r2found = (furthest_firstbin+4)*FOF_CONV; 
        FOFloat inv_enc_den = 0.;

	// Consider each bin r outwardly 
	for (int r=0; r<furthest_firstbin+4; r++) {
	    FOFr2 = FOF_CONV*r;
       	    FOFr2 *= FOFr2;
	    x = xthreshold*FOFr2;
	    
	    // Partition the newly touched cells 
	    for (int i = 0; i<ncg; i++) {
	        if (socg[i]->firstbin == r) {
		    // Get start[1] thru start[3] for every new SOgroupcell
		    partition_cellgroup(socg[i], hcenter); 
		}
	    }

	    // Is there enough mass within to skip or not enough to look
	    if (x*sqrt(x) < mass) { 
 	        // Enough mass, so crossing cannot be inside and we just add that mass
	        for (int i = 0; i<ncg; i++) {
	            if (socg[i]->firstbin >= r-3 && socg[i]->firstbin <= r) {
		        // Number of particles in this cell in this partition
		        mass += start[r-socg[i]->firstbin+1]-start[r-socg[i]->firstbin];
		    }
		}
	    }
	    
	    // Not enough mass, so could have crossing
	    else { 
	        // Number of all particles that are within this radial bin
	        int size_comb = 0;
		// number of particles within threshold in that partition
		int size_partial = 0;
		
		// for every SOcellgroup crossed by the radial bin
		for (int i = 0; i<ncg; i++) {
	            if (socg[i]->firstbin >= r-3 && socg[i]->firstbin <= r) {
		      
		        // Number of particles for that radial bin in that SOcellgroup
 		        int size_partition = (socg[i]->start[r-firstbin+1]-socg[i]->start[r-firstbin]);
			// Compute d2 for them
			// B.H. I should be able to use what I computed already in partition but not sure how 
			FOFloat *d2_partition = compute_d2(center, p+socg[i]->start[r-firstbin], size_partition, d2buffer, numdists);
			
			//Copy d2 of particles into combined list.
			for (int j = 0; j<size_partition; j++) {
			    d2comb[size_comb++] = d2_partition[j];
			}
		    }
		}
		
		// Search for density threshold in list, given previous mass.
		FOFloat d2_thresh = partial_search(d2comb, size_partition, mass, size_partial, inv_enc_den);

		// If something was found, record it
	        if (d2_thresh > 0) {
		    mass += size_partial;
		    r2found = d2_thresh;
		    break;  
		}
		// False alarm -- add all mass in that radial bin
		else {
		  mass += size_comb;
		}
	    }
	}
	// Record inverse density and threshold radius
	inv_enc_den = (x*sqrt(x))/(mass*threshold);
	return r2found;
    }   
    // B.H. new end


  
    // ======================  Partition ===========================

    /// Partition on d2buffer<=d2SO in the [start,last) region.
    /// Also find the index of the densest particle in the high-separation
    /// region.
    /// 
    /// Need to supply a pointer to d2 that match [start,last)
    /// indexing of the p and density list
    ///
    /// We also return 'size', which should be the first index of the high
    /// separation region, relative to start.  This may differ from 
    /// what was found externally if there is some conspiracy about 
    /// floating point equality comparisons of d2.
    int partition_and_index(FOFloat *d2, FOFloat d2SO, 
    	int start, int last, int &size) {

	FOFloat maxdens=-1.0;
	int densest=-1;
	if (start==last) {
	    return -1;    // Nothing to do, but this should never happen
	}
	int s = start;
        while (s<last && d2[s]<=d2SO) s++;
		// Advance start to the first high spot
	while (s<last) {
	    last--; // Consider next element in the upper list
	    if (d2[last]<=d2SO) {
		// We've found an out of place one; flip it with start
		std::swap(p[s],p[last]);
		std::swap(density[s],density[last]);
		std::swap(d2[s],d2[last]);
		s++;
		while (s<last && d2[s]<=d2SO) s++;
		// Advance s to the next high spot
	    }
	    // last points at a high particle.  Check if it's the densest
	    if (density[last]>maxdens) {
		maxdens = density[last]; densest = last;
	    }
	}
	// We're done when start = last.  last is the first high particle
	size = last-start;
	return densest;
    }
  
    int partition_alt_and_index(int *halos, int halo_i, 
    	int start, int last, int &size) {
        // similar to the original partition function
        int N_p = last;
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
	    last--; // Consider next element in the upper list
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
	int size = 0;	    /// The number of particles interior to the search range

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
	    if (b>=low && b<hi) d2sort[++size] = d2use[j];
	    // Need the >= here because if low==0, we need the d2=b=0 self-count
	}
	// Now the particles to be considered are in (start,size]
	// Sort this smaller segment
	std::sort(d2sort+start+1, d2sort+size+1);
	numsorts += size-start;

	// Now sweep in from the end to find the threshold
	int max = std::max(start,SO_CACHE);
	for (; size>max; size--) {
	    FOFloat x = d2sort[size]*xthreshold;
	    if (x*sqrt(x)<size) goto FoundSORadius;
		// This particle exceeds the overdensity threshold
	}
	for (; size>start; size--) {
	    FOFloat x = d2sort[size]*xthreshold;
	    if (x<twothirds[size]) goto FoundSORadius;
		// This particle exceeds the overdensity threshold
		// We have cached the size^(2/3) values, just to avoid
		// the square roots
	}

	FoundSORadius:
	if (size==start) return low/norm;
	    // We don't have the particle from the next lower bin, 
	    // so just return the boundary
	return d2sort[size];
    }

    // ---------------------------------------------------
    /// This is the original simple algorithm, retained for comparison testing

    FOFloat original_search_out(FOFloat *d2use, int len) {
	// We know that the group cannot be bigger than N=len,
	// and therefore that any x>len^{2/3} cannot matter.
	FOFloat x0 = pow(len,2.0/3.0)/xthreshold;
	int size = 0;	// We're going to unit-index d2sort[]


	for (int j=0; j<len; j++) 
	    if (d2use[j]<=x0) 
		d2sort[++size] = d2use[j];
	// Now we've filled d2sort in the range [1,size] with only 
	// the points that matter.
	// Sort the distances in increasing order
	std::sort(d2sort+1, d2sort+size+1);
	numsorts += size;
	// Now sweep in from the end to find the threshold
	for (; size>1; size--) {
	  FOFloat x = d2sort[size]*xthreshold;	  
	  if (x*sqrt(x)<size)   break;
	  // This particle exceeds the overdensity threshold
	}
	return d2sort[size]; 
    }



    FOFloat original_search(FOFloat *d2use, int len, FOFloat &inv_enc_den) {
	int size = 0;	// We're going to unit-index d2sort[]
	for (int j=0; j<len; j++) 
	  d2sort[++size] = d2use[j]; // maybe optimize tuk
	// Now we've filled d2sort in the range [1,size] with 
	// all the points.
	// Sort the distances in increasing order
	std::sort(d2sort+1, d2sort+size+1); // d2sort starts at 1, not 0
	numsorts += size;
	// Now sweep in from the center to find the threshold
	size = 1;
	FOFloat x;
	for (int j=1; j<=len; j++) {
	  x = d2sort[j]*xthreshold;
	  size = j; // we want the rightmost on the left side of the density threshold
	  if (x*sqrt(x)>(j) && d2sort[size]>=.5*WriteState.DensityKernelRad2) {
	    break;
	    // This particle exceeds the overdensity threshold
	  }
	}
	inv_enc_den = (x*sqrt(x))/(size*threshold);
	//if (size==len) printf("IT HAPPENS\n");//return d2sort[size]*1.34; // if too big radius, size_pr is also outside

	
	//if (size<=5) return d2sort[5];
	//if (size==1) {
	//return 0.5*(d2sort[size]+d2sort[size+1]); // no division by 0 if threshold
	// is right after the first particle
	//}
	return (d2sort[size]); 
    }
  
    // ======================  SO Algorithm  ===========================
    void greedySO() {
	int start = 0;	// Index of the first particle of the current group
	int densest = -1;
	FOFloat maxdens = -1.0;

	// TESTING
	FOFloat FOFunitdensity = P.np*4.0*M_PI*2.0/15.0*pow(WriteState.DensityKernelRad2,2.5);
	
	int count = 1;        /// The number of the group
	
	Sweep.Start();
	for (int j=0; j<np; j++) {
	  // initializing the arrays
	  min_inv_den_part[j] = 1.e30;
	  halo_part[j] = 0; // unassigned is 0; first is 1; second, third is 2, 3, etc.
	  // looking for densest pcle
	  if (density[j]>maxdens) {
	    maxdens=density[j];
	    densest = j;
	  }
	}
	Sweep.Stop();
	start = densest;	

	while (start>=0) {
	// Find the densest particle, move it to the front
	  // No eligible central is left.  But we always 
	  // try at least one central.

	    // TESTING
	    //if (count>50) break; // cap on number of halos per FOF group	    
	    Distance.Start();
	    // compute the distance to all particles
	    FOFloat *d2use = compute_d2(p+start, p, np, d2buffer, numdists);
	    Distance.Stop();
	    // d2use points to element start, and our interest is [start,np)
	    // In this next part, we need to arrive at a value of 
	    // d2SO with the property that all of the particles with 
	    // d2 <= d2SO are the largest body that satisfies 
	    // the overdensity condition.
	    Search.Start();

	    // B.H. new
	    // First fill the socg array with all cells in which pcles of this FOF reside
	    // and distances d2 to this nucleus (i.e. call fill_socg(p[start],np)).
	    // Instead of computing distances to all particles as done in compute_d2 
	    // call compute_d2(p[start]) for each socg entry 
	    // to get the distances to all cells
	    // (so move function here and rename to compute_cell_d2 or keep in other class but
	    // need to check syntax). Finally substitute original_search with
	    // search_socg_thresh(p[start], socg) and voila
	    // B.H. end new
	    
	    FOFloat inv_enc_den;
	    FOFloat d2SO = original_search(d2use,np,inv_enc_den);
	    
	    FOFloat d2SO_pr = d2SO*inner_rad2;
	    Search.Stop();
	    FOFloat inv_d;	    // inverse density
	    //FOFloat inv_d2del = 1./(d2SO*threshold); 	    // product of density and distance threshold
	    FOFloat inv_d2del = inv_enc_den/(d2SO); 	    // product of density and distance threshold
	    int size;

	    int densest = -1;
	    FOFloat maxdens = -1.0;

	    Sweep.Start();
	    for  (int j=0; j<np; j++) {
	      // for those within prime, set to inactive and count
	      if (d2use[j]<d2SO_pr) {
		//if (halo_part[j] == 0) halo_part[j] = -1;
		//else halo_part[j] |= 0x80000000;
		halo_part[j] = -(abs(halo_part[j])+(halo_part[j]==0));
	      }
	      // look for the next densest particle which is still active
	      // check the second condition tuk
	      // TESTING
	      //else if (density[j]>maxdens && density[j]*min_inv_den_part[j]>mag_loc && halo_part[j]>=0) {
	      else if (density[j]>maxdens && density[j]*min_inv_den_part[j]/FOFunitdensity>mag_loc && halo_part[j]>=0) {
		  maxdens=density[j];
		  densest = j;
		  if (d2use[j]>d2SO) continue; // in this way the next line of code is only exec if pcle j is within d2SO
		}
	      else if (d2use[j]>d2SO) continue; // in this way the next line of code is only exec if pcle j is within d2SO
		
	      // interpolate to get the density of the particle
	      inv_d = d2use[j]*(inv_d2del);
	      // if j is the densest particle seen so far, mark this as its halo
	      if (min_inv_den_part[j] > inv_d) {
		min_inv_den_part[j] = inv_d;
		//halo_part[j] = ((halo_part[j] & 0x80000000) | count);
		halo_part[j] = ((halo_part[j] >= 0)-(halo_part[j] < 0))*count;
	      }
	    }
	    Sweep.Stop();

	    

	    // if you have no particles in the innermost sphere, give up on this halo
	    count++; 
	    start = densest;
	}


	/*
	// for testing 
	// except we now have negative numbers, too
	int size;
	start = 0;
	
	for (int i=1; i<count; i++) {
	    size = 0;
	    //printf("i = %6i\n",i);
	    for (int j=0; j<np; j++) {
	      if ((i == halo_part[j]) || (i == -halo_part[j])) size++;
	    }
	    //printf("count, size_group = %6i, %6i\n",i,size);
	    if (size > 0) {
	      numcenters++;
	      groups[ngroups++] = FOFgroup(start,size);
	      start += size;
	    }
	}
	size = 0;
	*/

	for (int j=0; j<np; j++) {
	  halo_part[j] = abs(halo_part[j]);
	}
	
	Sweep.Start();
	int size = 0;
	start = 0;
	int next_densest;
	int halo_ind;
	for (int i=0; i<count; i++) {
	  if (i == 0) halo_ind = 1; // start with first subhalo                
          else if (i == 1) halo_ind = 0; // deal with the unassigned fluff     
          else halo_ind = i; // rest
	  //for (int i=0; i<count; i++) {
	  //if (i == 0) halo_ind = 1; // start with first subhalo
	  //else if (i == 1) halo_ind = 0; // deal with the unassigned fluff
	  //else halo_ind = i; // rest
	  
	  // rearrange so that you first sort to the left all particles in halo 0,
	  // then all in halo 1, ... i, ... count;
	  // and finding the densest particle of those to the right of i while doing so
	  next_densest = partition_alt_and_index(halo_part, halo_ind, start, np, size);
	  //printf("count, size_group, next dens, np = %6i, %6i, %6i, %6i\n",halo_ind,size,next_densest, np);
	  // Mark the group.  Note that the densest particle may not be first.
	  if (next_densest < 0) {
	    numcenters++;
	    if (halo_ind > 0) groups[ngroups++] = FOFgroup(start,size);
	    start += size;
	    continue;
	  }
	  if (size > 0) {
	    numcenters++;
	    if (halo_ind > 0) groups[ngroups++] = FOFgroup(start,size);
	    start += size;
	  }
	  
	  // swap so that remaining particles start with the densest
	  std::swap(p[start], p[next_densest]);
	  std::swap(density[start], density[next_densest]);
	  std::swap(halo_part[start], halo_part[next_densest]);
	}
	Sweep.Stop();
	
    }
    }
    }

// ---------------- Routines for cell indexes and cell groups -------------

    // This provides a simple parsing of the positions back into uniquely
    // numbered cell indices.  
    // NOTE: This assumes that particles occupy [-halfinvcpd,+halfinvcpd) in cells
    inline int compute_cellindex(posstruct &p) {
        int i = floor((p.x+CP->halfinvcpd)*CP->cpd)+128;
        int j = floor((p.y+CP->halfinvcpd)*CP->cpd)+128;
        int k = floor((p.z+CP->halfinvcpd)*CP->cpd)+128;
        return (i<<16)|(j<<8)|k;
    }

    /// Given the cellindex number, return the cell center
    inline FOFparticle compute_cellcenter(int cellidx) {
        int k = (cellidx&0xff);
        int j = (cellidx&0xff00)>>8;
        int i = (cellidx&0xff0000)>>16;
        posstruct p;
        p.z = CP->invcpd*(k-128);
        p.y = CP->invcpd*(j-128);
        p.x = CP->invcpd*(i-128);
        return FOFparticle(p,0);
    }

    /// Given the cellindex[] array, we want to scan through
    void load_socg() {
        int lastidx = 0x0f000000;   // An impossible index
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


// ------------------- Main calling routine ---------------------------------

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
        for (int j=0; j<np; j++) {
            p[j] = FOFparticle(pos[j],j);
            density[j] = acc[j].w;   // TODO: Check syntax
            cellindex[j] = compute_cellindex(pos[j]);
        }
        Copy.Stop();

        // Load the cellgroup list from the indices
        load_socg();

	
// B.H. say which cell each particle in SOcell belongs to and then loop over each
// member of cellindex which should have as many members as the whole pos array
// to find i guess the start[0] aka first particle belonging to this cell and possibly
// total number of dudes in that cell who are contiguous members


/* TODO:

    Loop over SOcellgroup's to compute_d2() for this halo center
    // B.H. aka find the minimum distance to each cell in the SOcellgroup array
    // aka distances should be same size as socellgroup array
*/

    // TODO question: After this point, do we ever use the particle cellindex again?
    // If not, then let's not permute it.
    // I think we agreed that the particle index is in fact all that is needed
    // if we want to put the L1 particles back into cellgroup order.

        greedySO();

    // Optional: If L1, then we could sort the FOFparticles of each outgoing group into particle index order, as that will re-order into cellgroups and speed up L2.  
        // TODO: Didn't add any code to skip this for L2
        // Note: I left the first particle unsorted, as it is planned to be the densest one.  This is a tiny inefficiency for L2: one extra group.
        for (int g=0; g<ngroups; g++) {
            std::sort(p+groups[g].start+1, p+groups[g].start+groups[g].n);
        }

        Total.Stop();
        return ngroups;
    }

};

