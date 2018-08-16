/** We compute properties of each L1 group and store them in a HaloStat
object.
*/

#ifndef __HALOSTAT_HH
#define __HALOSTAT_HH

#define N_LARGEST_SUBHALOS 3

class HaloStat {
  public:
    uint64_t id;	///< A unique halo number.
    uint64_t npstart;	///< Where to start counting in the particle output
    uint64_t npout;	///< Number of taggable particles pos/vel/aux written out.

    uint64_t taggedstart;   ///< Where to start counting in the tagged particle output
    uint64_t ntagged;	    ///< Number of tagged particle PIDs written out.

    uint32_t N;	///< The number of particles in this halo
    uint32_t subhalo_N[N_LARGEST_SUBHALOS];   ///< The number of particles in the largest L2 subhalos
    uint32_t L0_N;    ///< The number of particles in the L0 parent group

    float x[3];      ///< Center of mass position
    float v[3];      ///< Center of mass velocity
    float sigmav[3];  ///< sqrt(Eigenvalues) of the velocity tensor
    float r25, r50, r75, r90;   ///< Radii of this percentage of mass
    // float r100;
    float vcirc_max, rvcirc_max;  ///< max velocity and radius thereof
    // float xcen[3];	///< Center of the SO sphere
    // float densitycen;  ///< FOF-scale density of the SO central particle

    // The largest (most massive) subhalo center of mass
    float subhalo_x[3];   ///< Center of mass pos of the largest L2 subhalo
    float subhalo_v[3];   ///< Center of mass vel of the largest L2 subhalo
    // The profile properties computed from that center point
    float subhalo_sigmav[3]; ///< sqrt(Eigenvalues) of the velocity tensor of the L1 halo from L2 center
    float subhalo_r25, subhalo_r50, subhalo_r75, subhalo_r90;
    	///< Radii of this percentage of mass, relative to L2 center
    // float subhalo_r100;
    float subhalo_vcirc_max, subhalo_rvcirc_max;   ///< max circular velocity and radius thereof, relative to L2 center
    // float subhalo_xcen[3];	///< Center of the SO sphere for largest L2 subhalo
    // float subhalo_densitycen;  ///< FOF-scale density of the central particle of the largest L2 subhalo
};


class RVfloat {
  public:
    float pos[3];
    float vel[3];
    RVfloat(float px, float py, float pz, float vx, float vy, float vz) {
    	pos[0] = px; pos[1] = py; pos[2] = pz;
    	vel[0] = vx; vel[1] = vy; vel[2] = vz;
    }
};

class TaggedPID {
  public:
    uint64_t _pid;
    uint64_t pid() { return _pid; }
    TaggedPID(uint64_t p) { _pid = p; }
};

class RVfloatPID {
  public:
    uint64_t pid;
    float pos[3];
    float vel[3];
    RVfloatPID(uint64_t _pid, float px, float py, float pz, float vx, float vy, float vz) {
        pid = _pid;
        pos[0] = px; pos[1] = py; pos[2] = pz;
        vel[0] = vx; vel[1] = vy; vel[2] = vz;
    }
};

/*
class TaggedPID {
  public:
    // Cram the PID down into 5 bytes
    unsigned char v[5];
    uint64_t pid() {
        return ((((uint64_t)v[4]*256+(uint64_t)v[3])*256+(uint64_t)v[2])*256+(uint64_t)v[1])*256+(uint64_t)v[0];
    }
    TaggedPID(uint64_t p) {
	uint64_t q = p/256;
	v[0] = (unsigned char)(p-q*256);
	p = q;
	q = p/256;
	v[1] = (unsigned char)(p-q*256);
	p = q;
	q = p/256;
	v[2] = (unsigned char)(p-q*256);
	p = q;
	q = p/256;
	v[3] = (unsigned char)(p-q*256);
	p = q;
	q = p/256;
	v[4] = (unsigned char)(p-q*256);
	assert(q==0);
	return;
    }
};
*/

#endif
