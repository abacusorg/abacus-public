/** Code to perform spherical overdensity finding on provided particle
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

This will use a greedy algorithm.  We find the densest particle, find
the distance to all others and sort in increasing order.  Then work from
the outside to find the largest radius that has an enclosed overdensity
above the threshold.  That defines the first group.  Then repeat with
the next densest particle, excluding those particles already in a 
group.

We are not implementing any unbinding here.

*/

class SOcell {
  public:
    FOFparticle *p;     ///< The particles
    float *density;	///< The densities
    float *d2buffer;    ///< A buffer of distances
    float *d2sort;      ///< A buffer of sorted distances
    int np;             ///< The number of particles in this group
    int maxsize;        ///< The maximum number of particles
    
    FOFloat threshold;  ///< The density threshold we are applying
    FOFloat xthreshold;
    FOFloat min_central;   ///< The minimum FOF-scale density to require for a central particle.

    FOFgroup *groups;   ///< The list of found groups
    int ngroups;        ///< The number of groups
    
    FOFTimer time_total;
    FOFTimer time_copy;
    long long numdists;

    char pad[64];    // Just to ensure an array of these always fall on
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
        ret = posix_memalign((void **)&p, 64, sizeof(FOFparticle)*maxsize);  assert(ret == 0);
        memset(p, 0, sizeof(FOFparticle)*maxsize);
        if (d2buffer!=NULL) free(d2buffer);
        ret = posix_memalign((void **)&d2buffer, 64, sizeof(float)*maxsize);  assert(ret == 0);
        if (d2sort!=NULL) free(d2sort);
        ret = posix_memalign((void **)&d2sort, 64, sizeof(float)*maxsize);  assert(ret == 0);
        if (groups!=NULL) free(groups);
        ret = posix_memalign((void **)&groups, 64, sizeof(FOFgroup)*maxsize);  assert(ret == 0);
        if (density!=NULL) free(density);
        ret = posix_memalign((void **)&density, 64, sizeof(float)*maxsize);  assert(ret == 0);
    }

    SOcell() {
        p = NULL;
        d2buffer = NULL;
        d2sort = NULL;
        groups = NULL;
        density = NULL;
        maxsize = -1;
        numdists = 0;
        reset(32768);
        // We will have terrible bugs if these aren't true!
        assertf(sizeof(FOFparticle)==16, "FOFparticle is illegal size!");
        assertf(sizeof(FOFgroup)%16==0, "FOFgroup is illegal size!");
        return;
    }

    void destroy() {
        if (p!=NULL) free(p); p = NULL;
        if (d2buffer!=NULL) free(d2buffer); d2buffer = NULL;
        if (d2sort!=NULL) free(d2sort); d2sort = NULL;
        if (groups!=NULL) free(groups); groups = NULL;
        if (density!=NULL) free(density); density = NULL;
    }
    ~SOcell() {
        destroy();
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
	// mass/(4*PI/3)/r^3.  So that means r^3*threshold*3/4/PI = count
	// Define T = (threshold*3/4/PI)**(1/3) in code units, so that
	// we just compare (r*T)**3 to count.  The positions are in 
	// FOF_RESCALE units, so the interparticle spacing is P.invcpd*FOF_RESCALE
	threshold = _threshold;
	xthreshold = pow(threshold*3.0/4.0/M_PI, 1.0/3.0)*P.cpd/FOF_RESCALE;

	// Cosmic unit density yields a count in our FOFscale densities of this:
	float FOFunitdensity = P.np*4.0*M_PI*2.0/15.0*pow(WriteState.DensityKernelRad2,2.5)+1e-30;
	min_central = _min_central*FOFunitdensity;
        return;
    }

    /// Partition on d2buffer<=d2SO in the [start,last) region.
    /// Also find the index of the densest particle in the high-separation
    /// region.
    /// 
    /// Need to supply a pointer to d2 that match [start,last)
    /// indexing of the p and density list
    ///
    /// We also return 'end', which should be the first index of the high
    /// separation region, relative to start.  This may differ from 
    /// what was found externally if there is some conspiracy about 
    /// floating point equality comparisons of d2.
    int partition_and_index(float *d2,
    	float d2SO, int start, int last, int &end) {

	float maxdens=-1.0;
	int densest=-1;
	if (start==last) return -1;    // Nothing to do, but this should never happen

        while (start<last && d2[start]<=d2SO) start++;
		// Advance start to the first high spot
	while (start<last) {
	    last--; // Consider next element in the upper list
	    if (d2buffer[last]<=d2SO) {
		// We've found an out of place one; flip it with start
		std::swap(p[start],p[last]);
		std::swap(density[start],density[last]);
		std::swap(d2[start],d2[last]);
		start++;
		while (start<last && d2[start]<=d2SO) start++;
		// Advance start to the next high spot
	    }
	    // last points at a high particle.  Check if it's the densest
	    if (density[last]>maxdens) {
		maxdens = density[last]; densest = last;
	    }
	}
	// We're done when start = last.  last is the first high particle
	end = last-start;
	return densest;
    }


    void greedySO() {
	int start = 0;	// Index of the first particle of the current group
	int densest = -1; 
	float maxdens = -1.0;
	for (int j=0; j<np; j++) 
	    if (density[j]>maxdens) {
	    	maxdens=density[j]; densest = j; 
	    }

	while (start<np) {
	// Find the densest particle, move it to the front
	    if (start>0 && density[densest]<min_central) break;
	    	// No eligible central is left.  But we always 
		// try at least one central.
	    std::swap(p[start], p[densest]);
	    std::swap(density[start], density[densest]);
	    // Compute the distances to all of the unassigned particles
	    int len = np-start-1;    // TODO: Maybe not -1?
	    float *d2use = compute_d2(p+start, p+start, len, d2buffer, numdists);
	    // d2use points to element start
	    // Sort the distances in increasing order
	    memcpy(d2sort+1, d2use, sizeof(float)*len);
	    // Change notation so that start particle is index 1,
	    // last actual particle in the list is index len+1 = np-start
	    std::sort(d2sort+2, d2sort+len+2);
	    d2sort[1] = 0.0;
	    // Now sweep in from the end to find the threshold
	    int end = len+1;   // Start with the last particle
	    for (; end>1; end--) {
		FOFloat x = d2sort[end]*xthreshold;
		if (x*sqrt(x)<end) break;
		    // This particle exceeds the overdensity threshold
	    }
	    // Now the group runs [1,end] relative to start-1.
	    // which means [start,start+end) in the original list
	    FOFloat d2SO = d2sort[end];	
	    // This is the d2 value that marks the output.
	    // Partition the original lists based on this boundary.
	    // I expect it is easier to do this partition than to 
	    // sort the p & density lists along with d2; we don't
	    // need the particles to be radial sorted.
	    // Use this opportunity to return the index of the densest
	    // particle in the exterior set.
	    densest = partition_and_index(d2use-start, d2SO, start, np, end);
	    // Mark the group.  Note that the densest particle is first.
	    groups[ngroups++] = FOFgroup(start,end);
	    start += end;
	}
    }

    int findgroups(posstruct *pos, velstruct *vel, auxstruct *aux, accstruct *acc, int n) {
        // This takes in all of the positions, velocities, and aux for
        // a cell of size n.

        // It copies pos[0,n) into the FOFparticle array.
	// vel and aux are not used, but acc[0,n) is required, as it carries
	// the FOF-scale densities that we use.
        // Then perform the SO.
        // Return the number of multiplet groups.
        // pos[], etc. will be unchanged.

        time_total.Start();
        reset(n);
        np = n;

        // Load the particles
        time_copy.Start();
        for (int j=0; j<np; j++) {
            p[j] = FOFparticle(pos[j],j);
	    density[j] = acc[j].w;   // TODO: Check syntax
        }
        time_copy.Stop();

	greedySO();

	time_total.Stop();
	return ngroups;
    }
};


