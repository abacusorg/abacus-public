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




/// This class stores the position and velocity in a 20+12 bit format.
/// The positions are supplied in the range [-0.5,0.5)
/// For a 2 Gpc/h box, this implies resolution of about 2 kpc/h.
/// The velocities are supplied in km/s and will saturate above +-6000 km/s.
/// This means that velocities have a resolution of 3 km/s, injecting 1 km/s of rms rounding error.
/// This saturation level is motivated by being mildly above the expected
/// escape velocity of the largest clusters.
class RVint {
  public:
    int32_t pv[3];
    inline int32_t pack_pos(float x) {
        int32_t ix = round(x*1048576);
        if (ix<-524288) ix+=1048576;
        if (ix>=524288) ix-=1048576;
        return (ix*4096)&0xfffff000;
    }
    inline int32_t pack_vel(float v) {
        const float velscale = 6000.0;   // km/s
        int iv = round((v/velscale+1.0)*2048.0);
        if (iv<0) iv=0;
        if (iv>4095) iv=4095;    // We just saturate on super-velocity particles.
        return iv;
    }

    RVint(float px, float py, float pz, float vx, float vy, float vz) {
        pv[0] = pack_pos(px)|pack_vel(vx);
        pv[1] = pack_pos(py)|pack_vel(vy);
        pv[2] = pack_pos(pz)|pack_vel(vz);
    }

    /// This is the code to undo the packing above (for one coordinate).
    // This is amenable to vectorization
    unpack(int32_t input, float &pos, float &vel) {
        int iv = input&0xfff;
        const float velscale = 6000.0;   // km/s
        vel = velscale*(iv-2048);   // km/s
        int ix = input-iv;   // Slightly more expensive, but safer if input ended up cast to int64
        // int ix = input&0xfffff000;         // Alternate
        const float posscale= pow(2.0,-32.0); // Return to [-0.5,0.5)
        // Alternatively, one might prefer include the boxsize
        // const float posscale= boxsize*pow(2.0,-32.0); // Return to [-L/2,L/2)
        pos = ix*posscale;   
    }
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
