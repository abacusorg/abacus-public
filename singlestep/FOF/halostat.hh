/** We compute properties of each L1 group and store them in a HaloStat
object.
*/

#ifndef __HALOSTAT_HH
#define __HALOSTAT_HH

// The following file contains the actual struct definition
#include "halostat_cstruct.h"

#define RVfloat RVint

// We could provide C bindings to the following code if we wanted
// (the RVint unpacking in particular)

/// This class stores the position and velocity in a 20+12 bit format.
/// The positions are supplied relative to the first cell of the group, and converted to box units in the range [-0.5,0.5)
/// The range [-0.5,0.5) is mapped to [-500,000,+500,000) so that there is a few % overflow capacity.
/// Periodic wrapping is not supplied; instead, we saturate at +-(2**19-2)
/// For a 2 Gpc/h box, this implies resolution of about 2 kpc/h.
/// The velocities are supplied in km/s and will saturate above +-6000 km/s.
/// This means that velocities have a resolution of 3 km/s, injecting 1 km/s of rms rounding error.
/// This saturation level is motivated by being mildly above the expected
/// escape velocity of the largest clusters.
class RVint {
  public:
    int32_t pv[3];
    inline int32_t pack_pos(float x) {
        int32_t ix = round(x*1000000);  // Round off to human-readable 1 million
        if (ix<-524286) ix = -524286;
        if (ix> 524286) ix = +524286;   // Special saturated values
        return (ix*4096);
    }
    inline int32_t pack_vel(float v) {
        const float velscale = 6000.0;   // km/s
        int iv = round((v/velscale+1.0)*2048.0);
        if (iv<0) iv=0;
        if (iv>4095) iv=4095;    // We just saturate on super-velocity particles.
        return iv;
    }

    RVint(float px, float py, float pz, float vx, float vy, float vz) {

        // pos from global group coming in in units relative to group's first cell.
        // calculate how much first cell is offset from box center, where box goes from -0.5 to 0.5. 
        // then wrap. 
        pv[0] = pack_pos(px)|pack_vel(vx);
        pv[1] = pack_pos(py)|pack_vel(vy);
        pv[2] = pack_pos(pz)|pack_vel(vz);
    }

    /// This is the code to undo the packing above (for one coordinate).
    // This is amenable to vectorization
    void unpack(int32_t input, float &pos, float &vel) {
        int iv = input&0xfff;
        const float velscale = 6000.0/2048.0;   // km/s
        vel = velscale*(iv-2048);   // km/s
        int ix = input-iv;   // Slightly more expensive, but safer if input ended up cast to int64
        // int ix = input&0xfffff000;         // Alternate
        const float posscale= pow(2.0,-12.0)/1.0e6; // Return to [-0.5,0.5)
        // Alternatively, one might prefer include the boxsize
        // const float posscale= boxsize*pow(2.0,-12.0)/1.0e6; // Return to [-L/2,L/2)
        pos = ix*posscale;   
    }

};

class RVFloat {
  public:
    float pos[3];
    float vel[3];
    RVFloat(float px, float py, float pz, float vx, float vy, float vz) {
    	pos[0] = px; pos[1] = py; pos[2] = pz;
    	vel[0] = vx; vel[1] = vy; vel[2] = vz;
    }
};

class TaggedPID {
  public:
    uint64_t _pid;
    uint64_t pid() { return _pid; }
    TaggedPID(auxstruct a) { _pid = a.aux & AUX_PID_TAG_DENS_MASK; }
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

#endif // __HALOSTAT_HH
