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
    int *act_part;	///< The min densities
    int *halo_part;	///< The min densities
    FOFloat *d2buffer;    ///< A buffer of distances
    FOFloat *d2sort;      ///< A buffer of sorted distances
    int np;             ///< The number of particles in this group
    int maxsize;        ///< The maximum number of particles
    
    FOFloat threshold;  ///< The density threshold we are applying
    FOFloat xthreshold;
    // B.H.
    FOFloat threshold_pr;  ///< The prime density threshold we are applying
    FOFloat xthreshold_pr;
    FOFloat magi = 2.2;
    FOFloat mag_loc = 1.;
    FOFloat min_central;   ///< The minimum FOF-scale density to require for a central particle.
    FOFloat *twothirds;	  ///< We compare x^(3/2) to integers so much that we'll cache integers^(2/3)

    FOFgroup *groups;   ///< The list of found groups
    int ngroups;        ///< The number of groups
    
    
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
        if (groups!=NULL) free(groups);
        ret = posix_memalign((void **)&groups, CACHE_LINE_SIZE, sizeof(FOFgroup)*maxsize);  assert(ret == 0);
        if (density!=NULL) free(density);
        ret = posix_memalign((void **)&density, CACHE_LINE_SIZE, sizeof(FOFloat)*maxsize);  assert(ret == 0);
	if (min_inv_den_part!=NULL) free(min_inv_den_part);
        ret = posix_memalign((void **)&min_inv_den_part, CACHE_LINE_SIZE, sizeof(FOFloat)*maxsize);  assert(ret == 0);
	if (act_part!=NULL) free(act_part);
        ret = posix_memalign((void **)&act_part, CACHE_LINE_SIZE, sizeof(int)*maxsize);  assert(ret == 0);
	if (halo_part!=NULL) free(halo_part);
        ret = posix_memalign((void **)&halo_part, CACHE_LINE_SIZE, sizeof(int)*maxsize);  assert(ret == 0);
    }

    SOcell() {
        p = NULL;
        d2buffer = NULL;
        d2sort = NULL;
        groups = NULL;
        density = NULL;
	min_inv_den_part = NULL;
	halo_part = NULL;
	act_part = NULL;
	maxsize = -1;
        numdists = 0;
        numsorts = 0;
        numcenters = 0;
        reset(32768);
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
        if (groups!=NULL) free(groups); groups = NULL;
        if (density!=NULL) free(density); density = NULL;
	if (min_inv_den_part!=NULL) free(min_inv_den_part); min_inv_den_part = NULL;
	if (act_part!=NULL) free(act_part); act_part = NULL;
	if (halo_part!=NULL) free(halo_part); halo_part = NULL;
	
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
	// B.H.
	threshold_pr = threshold*magi;
	xthreshold = pow(threshold*4.0*M_PI/3.0*P.np, 1.0/3.0)/FOF_RESCALE;
	xthreshold_pr = pow(threshold_pr*4.0*M_PI/3.0*P.np, 1.0/3.0)/FOF_RESCALE;
	// But we have distance squared, so we'll square this, so that
	// multiplying this by d2 gives (r*T)**2, and ^(3/2) of that 
	// is compared to count.
	xthreshold *= xthreshold;
	xthreshold_pr *= xthreshold_pr;

	// Cosmic unit density yields a count in our FOFscale densities of this:
	FOFloat FOFunitdensity = P.np*4.0*M_PI*2.0/15.0*pow(WriteState.DensityKernelRad2,2.5)+1e-30;
	FOFloat M_D = 30.;
	FOFloat alpha_safety = 1.;
	FOFloat sigma3 = M_D*sqrt(threshold*P.np/(48*M_PI*M_PI));
	min_central = alpha_safety*5./WriteState.DensityKernelRad2*pow(sigma3,2./3); 
	if (P.np > 0) min_central /= (P.np);
	min_central *= FOFunitdensity;
	
	
	
        return;
    }

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
	    return -1;    // Nothing to do
	}
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
	    if (density[last]>maxdens && halos[last]==halo_i+1) {
		  maxdens = density[last];
		  densest = last;
	    }
	    // last points at a high particle.  Check if it's the densest
	    
	}
	// We're done when start = last.  last is the first high particle
	size = last-start;
	
	return densest;
    }

    int partition_only(FOFloat *d2, FOFloat d2SO, 
    	int start, int last, int &size) {

	if (start==last) return -1;    // Nothing to do, but this should never happen

	int s = start;
        while (s<last && d2[s]<=d2SO) s++;
		// Advance start to the first high spot
	while (s<last) {
	    last--; // Consider next element in the upper list
	    if (d2[last]<=d2SO) {
		// We've found an out of place one; flip it with start
		std::swap(p[s],p[last]);
		
		std::swap(d2[s],d2[last]);
		s++;
		while (s<last && d2[s]<=d2SO) s++;
		// Advance s to the next high spot
	    }

	}
	// We're done when start = last.  last is the first high particle
	size = last-start;
	return 0;
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


    FOFloat original_search(FOFloat *d2use, int len) {
	int size = 0;	// We're going to unit-index d2sort[]
	for (int j=0; j<len; j++) 
	  d2sort[++size] = d2use[j];
	// Now we've filled d2sort in the range [1,size] with 
	// all the points.
	// Sort the distances in increasing order
	std::sort(d2sort+1, d2sort+size+1); // d2sort starts at 1, not 0
	numsorts += size;
	// Now sweep in from the center to find the threshold
	size = 1; 
	for (int j=1; j<=len; j++) {
	  FOFloat x = d2sort[j]*xthreshold;
	  size = j; // we want the rightmost on the left side of the density threshold
	  if (x*sqrt(x)>(j)) {
	    break;
	    // This particle exceeds the overdensity threshold
	  }
	}

	if (size==len) {
	  size = len;
	  d2sort[size] *= 1.34; // to ensure that if too big size_pr is also outside
	}
	if (size==1) {
	  return 0.5*(d2sort[size]+d2sort[size+1]); //i.e. don't want division by 0 if threshold
	  // is right after the first particle
	}
	return (d2sort[size]); 
    }
  
    // ======================  SO Algorithm  ===========================
    void greedySO() {
	int start = 0;	// Index of the first particle of the current group
	int densest = -1;

	FOFloat maxdens = -1.0;
	
	// B.H. delta prime determines the cores and delta the edges
	int count = 0;        /// B.H.The number of groups
	
	int ind_delta; /// which particle crosses the target density
	int ind_delta_pr; /// which particle crosses the target density prime
	
	Sweep.Start();
	for (int j=0; j<np; j++) {
	  // initializing the three arrays
	  min_inv_den_part[j] = 1.e30;
	  act_part[j] = 1;
	  halo_part[j] = np+1;
	  // looking for densest pcle
	  if (density[j]>maxdens) {
	    maxdens=density[j];
	    densest = j;
	  }
	}
	Sweep.Stop();
	start = densest;	
	
	while (start<np) {
	// Find the densest particle, move it to the front
	    if (start>0 && density[densest]<min_central) break;
	    	// No eligible central is left.  But we always 
		// try at least one central.
	    	
	    //if (count>50) break; // cap on number of halos per FOF group
	    
	    // for testing purposes
	    //if (np>256) break;
	    //if (np<9900) break;
	    
	    Distance.Start();
	    // B.H. compute the distance to all particles
	    FOFloat *d2use = compute_d2(p+start, p, np, d2buffer, numdists);
	    Distance.Stop();
	    // d2use points to element start, and our interest is [start,np)
	    // In this next part, we need to arrive at a value of 
	    // d2SO with the property that all of the particles with 
	    // d2 <= d2SO are the largest body that satisfies 
	    // the overdensity condition.
	    Search.Start();

	    
	    FOFloat d2SO = original_search(d2use, np)+1.e-2;
	    
	    FOFloat d2SO_pr = d2SO*.9;//*.9;//*0.75;// used to be .5 but a bit weird

	    // inverse density
	    FOFloat inv_d;
	   
	    FOFloat d2del = (d2SO*threshold); 
	    

	    int densest = -1;
	    FOFloat maxdens = -1.0;

	    Sweep.Start();
	    int size = 0;
	    int size_pr = 0;
	    int sums = 0;


	    for  (int j=0; j<np; j++) {
	      // for those within prime, set to inactive and count
	      if (d2use[j]<d2SO_pr) {
		act_part[j] = 0;
		size_pr++;
		size++;
	      }
	      // look for the next densest particle which is still active and
	      // only separating into two else statements
	      // so as to keep track of size of particles within Delta
	     
	      else if (d2use[j]<d2SO) {	      
		size++;
		if (density[j]>maxdens && density[j]*min_inv_den_part[j]>mag_loc && act_part[j]==1) {
		  maxdens=density[j];
		  densest = j;
		}
	      }
	      else {
		if (density[j]>maxdens && density[j]*min_inv_den_part[j]>mag_loc && act_part[j]==1) {
		  maxdens=density[j];
		  densest = j;
		}
		continue; // in this way the next line of code is only exec if pcle j is within d2SO
	      }
	      sums ++;
	      // interpolate to get the density of all the particles
	     
	      inv_d = d2use[j]/(d2del);
	     
	      
	      // if j is the densest particle seen so far, mark this as its halo
	      if (min_inv_den_part[j] > inv_d) {
		min_inv_den_part[j] = inv_d;
		halo_part[j] = count;
	      }
	    }
	    
	    // if you have no particles in the innermost sphere, give up on this halo
	    if (size < 1) {
	      break;
	    }
	    
	    Search.Stop();
	    
	    
	    
	    Sweep.Stop();

	    start = densest;

	    count++; 
	    if (start==-1) {
	      break;
	    }
	}


	/*
	// for testing  
	int size;
	start = 0;
       
	for (int i=0; i<count; i++) {
	  size = 0;
	  printf("i = %6i\n",i);
	  for (int j=0; j<np; j++) {
	    if (i == halo_part[j])) {
	      size++;
	    }
	  }
	  printf("count, size_group = %6i, %6i\n",i,size);
	  numcenters++; // new
	  //groups[ngroups++] = FOFgroup(start,size); // new
	  start+=size; // new
	}
	size = 0;
	*/

	
	int size = 0;
	start = 0;
	int next_densest;
	for (int i=0; i<count; i++) {
	  // rearrange so that you first sort to the left all particles in halo 0,
	  // then all in halo 1, ... i, ... count;
	  // and finding the densest particle of those to the right of i while doing so
	  next_densest = partition_alt_and_index(halo_part, i, start, np, size);
	  
	  // Mark the group.  Note that the densest particle may not be first.
	  if (next_densest < 0) {
	    numcenters++;
	    groups[ngroups++] = FOFgroup(start,size);
	    continue;
	  }
	  numcenters++;
	  groups[ngroups++] = FOFgroup(start,size);
	  start += size;

	  // swap so that remaining particles start with the densest
	  std::swap(p[start], p[next_densest]);
	  std::swap(density[start], density[next_densest]);
	  std::swap(halo_part[start], halo_part[next_densest]);

	}
	
    }


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
        }
        Copy.Stop();

	greedySO();

	Total.Stop();
	return ngroups;
    }
};


