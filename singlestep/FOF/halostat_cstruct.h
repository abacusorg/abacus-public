/* This file contains the pure C parts of our HaloStat definition.
   This facilitates interoperability with C and Python.
   C++ applications should not include this directly, but instead
   include "halostat.hh", which includes this file.

   Note: can't use #include guards here for CFFI compatibility!
*/

#define N_LARGEST_SUBHALOS 5

struct HaloStat {
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
    uint16_t sigmav_eigenvecs;  ///<Eigenvectors of the velocity dispersion tensor, compressed into 16 bits. 
    int16_t sigmavrad_to_sigmav3d; ///< sigmav_rad / sigmavSum, compressed
    int16_t sigmavtan_to_sigmav3d; ///< sigmav_tan / sigmavSum, compressed

	int16_t r10, r25, r33, r50, r67, r75, r90; ///<Expressed as ratios of r100, and scaled to 32000 to store as int16s. 
    int16_t sigmar[3]; ///<sqrt( Eigenvalues of the moment of inertia tensor ) 
    int16_t sigman[3]; ///<sqrt( Eigenvalues of the weighted moment of inertia tensor ) 
    uint16_t sigmar_eigenvecs;  ///<Eigenvectors of the moment of inertia tensor, compressed into 16 bits. Compression format TBD. 
    uint16_t sigman_eigenvecs;  ///<Eigenvectors of the weighted moment of inertia tensor, compressed into 16 bits. Compression format TBD. 
	int16_t rvcirc_max; ///< radius of max velocity, stored as int16 ratio of r100 scaled by 32000.

    // The largest (most massive) subhalo center of mass    
    int16_t L2cntr_sigmavMin_to_sigmav3d; ///< sigmav_z / sigmavSum, compressed
    int16_t L2cntr_sigmavMax_to_sigmav3d; ///< sigmav_z / sigmavSum, compressed
    uint16_t L2cntr_sigmav_eigenvecs;  ///<Eigenvectors of the velocity dispersion tensor, compressed into 16 bits. 
    int16_t L2cntr_sigmavrad_to_sigmav3d; ///< sigmav_rad / sigmavSum, compressed
    int16_t L2cntr_sigmavtan_to_sigmav3d; ///< sigmav_tan / sigmavSum, compressed
    int16_t L2cntr_r10, L2cntr_r25, L2cntr_r33, L2cntr_r50, L2cntr_r67, L2cntr_r75, L2cntr_r90;
    	///< Radii of this percentage of mass, relative to L2 center. Expressed as ratios of r100 and compressed to int16. 

	int16_t L2cntr_sigmar[3]; 
	int16_t L2cntr_sigman[3]; 
    uint16_t L2cntr_sigmar_eigenvecs;
    uint16_t L2cntr_sigman_eigenvecs;
    int16_t L2cntr_rvcirc_max;   ///< max circular velocity and radius thereof, relative to L2 center
};
