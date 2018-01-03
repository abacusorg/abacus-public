#define N_LARGEST_SUBHALOS 3

class HaloStat {
  public:
    uint64_t id;	// A unique halo number.
    uint64_t npstart;	// Where to start counting in the particle output
    uint64_t npout;	// Number of particles written out.
    uint32_t N;	// The number of particles in this halo
    uint32_t subhalo_N[N_LARGEST_SUBHALOS];

    float x[3];      // Center of mass position
    float v[3];      // Center of mass velocity
    float sigmav[3];  // sqrt(Eigenvalues) of the velocity tensor
    float r25, r50, r75, r90;   // Radii of this percentage of mass
    float vcirc_max, rvcirc_max;  // max velocity and radius thereof

    // The largest subhalo center of mass
    float subhalo_x[3];   
    float subhalo_v[3];
    // The profile properties computed from that center point
    float subhalo_sigmav[3];
    float subhalo_r25, subhalo_r50, subhalo_r75, subhalo_r90;
    float subhalo_vcirc_max, subhalo_rvcirc_max;
};

class TaggedPID {
  public:
    // Cram the PID down into 5 bytes
    unsigned char v[5];
    uint64_t pid() {
        return ((((uint64_t)v[4]*256+(uint64_t)v[3])*256+(uint64_t)v[2])*256+(uint64_t)v[1])*256+(uint64_t)v[0];
    }
    TaggedPID(uint64_t p) {
	uint64_t q = p/256;
	v[0] = (unsigned_char)(p-q*256);
	p = q;
	q = p/256;
	v[1] = (unsigned_char)(p-q*256);
	p = q;
	q = p/256;
	v[2] = (unsigned_char)(p-q*256);
	p = q;
	q = p/256;
	v[3] = (unsigned_char)(p-q*256);
	p = q;
	q = p/256;
	v[4] = (unsigned_char)(p-q*256);
	assert(q==0);
	return;
    }
}
