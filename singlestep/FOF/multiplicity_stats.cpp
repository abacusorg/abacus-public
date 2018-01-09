#define MS_NBIN 24
class MultiplicityStats {
  public:
    uint64 ngroups;
    uint64 largest;
    uint64 count[MS_NBIN], sumn[MS_NBIN], sumn2[MS_NBIN];
	// Count will be log2 bins.  np=1 is in bin 0.
	// np = 2-3 in bin 1, 4-7 in bin 2, etc

    MultiplicityStats() {
	ngroups = largest = 0;
	for (int j=0; j<MS_NBIN; j++) count[j] = sumn[j] = sumn2[j] = 0;
    }
    void add(MultiplicityStats &b) {
        // Add b to this one
	ngroups += b.ngroups;
	largest = std::max(largest, b.largest);
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
	int j=0;
	np = np>>1;
	while (np!=0) { j++; np = np>>1; }
	count[j]++;
	sumn[j] += np;
	sumn2[j] += np*np;
	return;
    }

    void report() {
        printf("Total number of groups %llu\n", ngroups);
	uint64 tot = 0, tot2 = 0;
	int j, m, nbin;
	for (j=0; j<MS_NBIN; j++) { 
	    tot += sumn[j];
	    tot2 += sumn2[j];
	}
        printf("Particles in groups %llu\n", tot);
	printf("Average group has %f particles and %f pairs\n", 
		(float)tot/ngroups, (float)tot2/ngroups);
	for (nbin==MS_NBIN-1; nbin>=0; nbin--) 
	    if (count[nbin]>0) break;
	    // nbin is now the number of the highest non-empty bin
	for (j=0,m=1; j<=nbin; j++, m*=2)
	    printf("%5d -- %5d: %7llu groups, %6.3f%% of particles, %6.3f%% of pairs\n",
	    	m, m*2-1, count[j],
		100.0*sumn[j]/tot, 100.0*sumn2[j]/tot2);
    }
};
