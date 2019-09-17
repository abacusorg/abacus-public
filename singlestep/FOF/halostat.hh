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
    //float r25, r50, r75, r90;   ///< Radii of this percentage of mass
    float r100; ///<Radius of 100% of mass 
	int16_t r10, r25, r50, r67, r75, r90; ///<Expressed as ratios of r100, and scaled to 32000 to store as int16s. 
	float sigmar100; ///<Second moments about the mean of/for all particles falling within r100. 
	int16_t sigmar10, sigmar25, sigmar50, sigmar67, sigmar75, sigmar90; ///<Second moments about the mean of/for all particles falling within r_n. 
	
    float   vcirc_max; ///< max velocity 
	int16_t rvcirc_max; ///< radius of max velocity, stored as int16 ratio of r100 scaled by 32000.

	FOFparticle SO_central_particle; ///< Coordinates of the SO central particle (densest particle). 
	FOFloat     SO_central_density;  ///< Density of the SO central particle. 
	FOFloat     SO_radius;           ///< Radius of SO halo (distance to particle furthest from central particle) 

    // The largest (most massive) subhalo center of mass
    float subhalo_x[3];   ///< Center of mass pos of the largest L2 subhalo
    float subhalo_v[3];   ///< Center of mass vel of the largest L2 subhalo
    // The profile properties computed from that center point
    float subhalo_sigmav[3]; ///< sqrt(Eigenvalues) of the velocity tensor of the L1 halo from L2 center

    float subhalo_r100; /// Radius of 100% of mass, relative to L2 center. 
    int16_t subhalo_r10, subhalo_r25, subhalo_r50, subhalo_r67, subhalo_r75, subhalo_r90;
    	///< Radii of this percentage of mass, relative to L2 center. Expressed as ratios of r100 and compressed to int16. 
	float subhalo_sigmar100; ///<Second moments about the mean of/for all particles falling within r100, relative to L2 center. 
	int16_t subhalo_sigmar10, subhalo_sigmar25, subhalo_sigmar50, subhalo_sigmar67, subhalo_sigmar75, subhalo_sigmar90; ///<Second moments about the mean of/for all particles falling within r_n, relative to L2 center. 
	
    float subhalo_vcirc_max, subhalo_rvcirc_max;   ///< max circular velocity and radius thereof, relative to L2 center

	FOFparticle SO_subhalo_central_particle; ///< Coordinates of the SO central particle (densest particle) for the largest L2 subhalo. 
	FOFloat     SO_subhalo_central_density;  ///< Density of the SO central particle of the largest L2 subhalo. 
	FOFloat     SO_subhalo_radius;           ///< Radius of SO halo (distance to particle furthest from central particle) for the largest L2 subhalo
	
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
