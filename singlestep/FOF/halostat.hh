/** We compute properties of each L1 group and store them in a HaloStat
object.
*/

#ifndef __HALOSTAT_HH
#define __HALOSTAT_HH

#define RVfloat RVint


#define N_LARGEST_SUBHALOS 5

class HaloStat {
  public:
    uint64_t id;	///< A unique halo number.
    uint64_t npstartA;	///< Where to start counting in the particle output for subsample A
    uint64_t npstartB;  ///< Where to start counting in the particle output for subsample B
    uint32_t npoutA;	///< Number of taggable particles pos/vel/aux written out in subsample A
    uint32_t npoutB;    ///< Number of taggable particles pos/vel/aux written out in subsample B
    uint32_t ntaggedA;	    ///< Number of tagged particle PIDs written out in subsample A. A particle is tagged if it is taggable and is in the largest L2 halo for a given L1 halo. 
    uint32_t ntaggedB; 
    uint32_t N;	///< The number of particles in this halo
    uint32_t L2cntr_N[N_LARGEST_SUBHALOS];   ///< The number of particles in the largest L2 subhalos
    uint32_t L0_N;    ///< The number of particles in the L0 parent group

    float x[3];      ///< Center of mass position
    float v[3];      ///< Center of mass velocity
    float sigmav3d;  ///< Sum of eigenvalues/3
    float r100; ///<Radius of 100% of mass 
    float vcirc_max; ///< max velocity 
    float SO_central_particle[3]; ///< Coordinates of the SO central particle (densest particle). 
    float SO_central_density;  ///< Density of the SO central particle. 
    float SO_radius;           ///< Radius of SO halo (distance to particle furthest from central particle) 
    float L2cntr_x[3];   ///< Center of mass pos of the largest L2 subhalo
    float L2cntr_v[3];   ///< Center of mass vel of the largest L2 subhalo
    float L2cntr_sigmav3d;  ///< Sum of eigenvalues/3 
    float L2cntr_r100; /// Radius of 100% of mass, relative to L2 center. 
    float L2cntr_vcirc_max; 
    float SO_L2max_central_particle[3]; ///< Coordinates of the SO central particle (densest particle) for the largest L2 subhalo. 
    float SO_L2max_central_density;  ///< Density of the SO central particle of the largest L2 subhalo. 
    float SO_L2max_radius;           ///< Radius of SO halo (distance to particle furthest from central particle) for the largest L2 subhalo

    int16_t sigmavMin_to_sigmav3d; ///< sigmav_z / sigmavSum, compressed

    int16_t sigmavMax_to_sigmav3d; ///< sigmav_z / sigmavSum, compressed
    // int16_t sigmav3d_to_sigmavMaj; ///< sigmav_x / sigmavSum, compressed
    //!!!int16_t sigmav_eigenvecs;  ///<Eigenvectors of the velocity dispersion tensor, compressed into 16 bits. Compression format TBD.
    // B.H.
    int16_t sigmavrad_to_sigmav3d; ///< sigmav_rad / sigmavSum, compressed
    int16_t sigmavtan_to_sigmav3d; ///< sigmav_tan / sigmavSum, compressed
    //float r25, r50, r75, r90;   ///< Radii of this percentage of mass

	int16_t r10, r25, r33, r50, r67, r75, r90; ///<Expressed as ratios of r100, and scaled to 32000 to store as int16s. 
    int16_t sigmar[3]; ///<sqrt( Eigenvalues of the moment of inertia tensor ) 
    uint16_t sigmar_eigenvecs;  ///<Eigenvectors of the moment of inertia tensor, compressed into 16 bits. Compression format TBD. 
	int16_t rvcirc_max; ///< radius of max velocity, stored as int16 ratio of r100 scaled by 32000.

	//!!!float     SO_radius;           ///< Radius of SO halo (distance to particle furthest from central particle) 

    // The largest (most massive) subhalo center of mass    
    int16_t L2cntr_sigmavMin_to_sigmav3d; ///< sigmav_z / sigmavSum, compressed
    int16_t L2cntr_sigmavMax_to_sigmav3d; ///< sigmav_z / sigmavSum, compressed
    // int16_t L2cntr_sigmav3d_to_sigmavMaj; ///< sigmav_x / sigmavSum, compressed
    int16_t L2cntr_sigmavrad_to_sigmav3d; ///< sigmav_rad / sigmavSum, compressed
    int16_t L2cntr_sigmavtan_to_sigmav3d; ///< sigmav_tan / sigmavSum, compressed
    //!!!int16_t L2cntr_sigmav_eigenvecs;  ///<Eigenvectors of the velocity dispersion tensor, compressed into 16 bits. Compression format TBD. 
    int16_t L2cntr_r10, L2cntr_r25, L2cntr_r33, L2cntr_r50, L2cntr_r67, L2cntr_r75, L2cntr_r90;
    	///< Radii of this percentage of mass, relative to L2 center. Expressed as ratios of r100 and compressed to int16. 

	int16_t L2cntr_sigmar[3]; 
    uint16_t L2cntr_sigmar_eigenvecs;
    int16_t L2cntr_rvcirc_max;   ///< max circular velocity and radius thereof, relative to L2 center

};




/// This class stores the position and velocity in a 20+12 bit format.
/// The positions are supplied relative to the first cell of the group, and converted to box units in the range [-0.5,0.5)
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
        const float posscale= pow(2.0,-32.0); // Return to [-0.5,0.5)
        // Alternatively, one might prefer include the boxsize
        // const float posscale= boxsize*pow(2.0,-32.0); // Return to [-L/2,L/2)
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

#endif
