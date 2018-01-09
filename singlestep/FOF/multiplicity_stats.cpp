#define MS_NBIN 24
#define MS_MIN 32
class MultiplicityStats {
  public:
    uint64 ngroups;
    uint64 largest;
    uint64 totn, totn2;
    uint64 pad[4];
    uint64 count[MS_NBIN], sumn[MS_NBIN], sumn2[MS_NBIN];
	// Count will be log2 bins.  np=1 is in bin 0.
	// np = 2-3 in bin 1, 4-7 in bin 2, etc

    MultiplicityStats() {
	ngroups = largest = totn = totn2 = 0;
	for (int j=0; j<MS_NBIN; j++) count[j] = sumn[j] = sumn2[j] = 0;
    }
    void add(MultiplicityStats &b) {
        // Add b to this one
	ngroups += b.ngroups;
	largest = std::max(largest, b.largest);
	totn += b.totn;
	totn2 += b.totn2;
	for (int j=0; j<MS_NBIN; j++) {
	    count[j] += b.count[j];
	    sumn[j] += b.sumn[j];
	    sumn2[j] += b.sumn2[j];
	}
    }

    void push(int np) {
        // Given a group multiplicity, add it to the list
	if (np>largest) largest=np;
	ngroups++;
	totn += np;
	totn2 = += np*np;
	if (np<MS_MIN) return;
	int j=0;
	int n = np>>1;
	while (n!=0) { j++; n = n>>1; }
	assert(j<MS_NBIN);
	count[j]++;
	sumn[j] += np;
	sumn2[j] += np*np;
	return;
    }

    void report() {
        printf("Total number of groups %llu\n", ngroups);
	int j, m, nbin;
        printf("Groups contain %llu particles\n", tot);
	printf("Average group has %f particles and %f pairs\n", 
		(float)tot/ngroups, (float)tot2/ngroups);
        printf("Largest Group contains %llu particles\n", largest);
	for (nbin=MS_NBIN-1; count[nbin]==0 && nbin>=0; nbin--);
	    // nbin is now the number of the highest non-empty bin
	printf("Max bin is %d\n", nbin);
	nbin = MS_NBIN-1;
	for (j=0,m=1; j<=nbin; j++, m*=2)
	    printf("%5d -- %5d: %7llu groups, %6.3f%% of particles, %6.3f%% of pairs\n",
	    	m, m*2-1, count[j],
		100.0*sumn[j]/tot, 100.0*sumn2[j]/tot2);
    }
};
