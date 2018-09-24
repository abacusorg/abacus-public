/** This provides a simple statistical summary of the group sizes,
particularly generating a histogram.  
*/

#define MS_NBIN 24
#define MS_MIN 16
class MultiplicityStats {
  public:
    uint64 ngroups;
    uint64 largest;
    uint64 tot, tot2;
    uint64 pad[4];
    uint64 count[MS_NBIN], sumn[MS_NBIN], sumn2[MS_NBIN];
	// Count will be log2 bins.  np=1 is in bin 0.
	// np = 2-3 in bin 1, 4-7 in bin 2, etc

    MultiplicityStats() {
	ngroups = largest = tot = tot2 = 0;
	for (int j=0; j<MS_NBIN; j++) count[j] = sumn[j] = sumn2[j] = 0;
    }
    void add(MultiplicityStats &b) {
        // Add b to this one
	ngroups += b.ngroups;
	largest = std::max(largest, b.largest);
	tot += b.tot;
	tot2 += b.tot2;
	for (int j=0; j<MS_NBIN; j++) {
	    count[j] += b.count[j];
	    sumn[j] += b.sumn[j];
	    sumn2[j] += b.sumn2[j];
	}
    }

    /// Given a group multiplicity, add it to the stats
    void push(uint64 np) {
	if (np>largest) largest=np;
	ngroups++;
	tot += np;
	tot2 += np*np;
	if (np<MS_MIN) return;
	int j=0;
	int n = np>>1;
	while (n!=0) { j++; n = n>>1; }
	assertf(j<MS_NBIN, "Group is larger than maximum multiplicity stats range!");
	count[j]++;
	sumn[j] += np;
	sumn2[j] += (uint64) np*np;
	return;
    }
	
    /// This generates a report for the log
    void report_multiplicities(std::ofstream *grouplog) {
        GLOG(0,"Total number of groups %f M\n", ngroups/1e6);
	int j, m, nbin;
        GLOG(0,"Groups contain %f M particles\n", tot/1e6);
	GLOG(0,"Average group has %f particles and %f pairs\n", 
		(float)tot/ngroups, (float)tot2/ngroups);
        GLOG(0,"Largest Group contains %u particles\n", largest);
	for (nbin=MS_NBIN-1; nbin>=0 && count[nbin]==0; nbin--);
	    // nbin is now the number of the highest non-empty bin
	GLOG(2,"Max bin is %d\n", nbin);
	for (j=0,m=1; j<=nbin; j++, m*=2)
	    if (count[j]>0) 
		GLOG(0,"%7d -- %7d: %8u groups, %6.3f%% of particles, %6.3f%% of pairs\n",
		    (m<MS_MIN?MS_MIN:m), m*2-1, count[j],
		    100.0*sumn[j]/tot, 100.0*sumn2[j]/tot2);
    }
};
