/* This file contains the pure C parts of our HaloStat definition.
   This facilitates interoperability with C and Python.
   C++ applications should not include this directly, but instead
   include "halostat.hh", which includes this file.

   Note: can't use #include guards here for CFFI compatibility!

   Quantities are typically computed either around the L1 center of mass (*_com)
       or around the center of mass of the largest L2 subhalo (*_L2com).
   Center of masses are based on 100% of the particles in the halo.
   All velocity stats are computed relative to the COM velocity.
   Second moments are computed using only the inner 90% of particles (by radius),
       except for a couple based on inner 50%.
   The weighted moment of inertia weighs the particles by 1/r, 
       so it is the second moment of the radial unit vectors.
*/

#define N_LARGEST_SUBHALOS 5

struct HaloStat {
    uint64_t id;    ///< A unique halo number.
    uint64_t npstartA;  ///< Where to start counting in the particle output for subsample A
    uint64_t npstartB;  ///< Where to start counting in the particle output for subsample B
    uint32_t npoutA;    ///< Number of taggable particles pos/vel/aux written out in subsample A
    uint32_t npoutB;    ///< Number of taggable particles pos/vel/aux written out in subsample B
    uint32_t ntaggedA;      ///< Number of tagged particle PIDs written out in subsample A. A particle is tagged if it is taggable and is in the largest L2 halo for a given L1 halo. 
    uint32_t ntaggedB; 
    uint32_t N; ///< The number of particles in this halo
    uint32_t L2_N[N_LARGEST_SUBHALOS];   ///< The number of particles in the largest L2 subhalos
    uint32_t L0_N;    ///< The number of particles in the L0 parent group

    float x_com[3];      ///< Center of mass position
    float v_com[3];      ///< Center of mass velocity
    float sigmav3d_com;  ///< Sum of eigenvalues
    float meanSpeed_com;  ///< Mean speed (the norm of the velocity vector)
    float sigmav3d_r50_com;  ///< Velocity dispersion of the inner 50% of particles
    float meanSpeed_r50_com;  ///< Mean speed of the inner 50% of particles
    float r100_com; ///<Radius of 100% of mass 
    float vcirc_max_com; ///< max circular velocity, based on the particles in this L1 halo
    float SO_central_particle[3]; ///< Coordinates of the SO central particle
    float SO_central_density;  ///< Density of the SO central particle. 
    float SO_radius;           ///< Radius of SO halo (distance to particle furthest from central particle) 

    float x_L2com[3];   ///< Center of mass pos of the largest L2 subhalo
    float v_L2com[3];   ///< Center of mass vel of the largest L2 subhalo
    float sigmav3d_L2com;  ///< Sum of eigenvalues
    float meanSpeed_L2com;  ///< Mean speed
    float sigmav3d_r50_L2com;  ///< Velocity dispersion of the inner 50% of particles
    float meanSpeed_r50_L2com;  ///< Mean speed of the inner 50% of particles
    float r100_L2com; /// Radius of 100% of mass, relative to L2 center. 
    float vcirc_max_L2com;   ///< max circular velocity, based on the particles in this L1 halo 
    float SO_L2max_central_particle[3]; ///< Coordinates of the SO central particle for the largest L2 subhalo. 
    float SO_L2max_central_density;  ///< Density of the SO central particle of the largest L2 subhalo. 
    float SO_L2max_radius;           ///< Radius of SO halo (distance to particle furthest from central particle) for the largest L2 subhalo

    int16_t sigmavMin_to_sigmav3d_com; ///< Min(sigmav_eigenvalue) / sigmav3d, compressed
    int16_t sigmavMax_to_sigmav3d_com; ///< Max(sigmav_eigenvalue) / sigmav3d, compressed
    uint16_t sigmav_eigenvecs_com;  ///<Eigenvectors of the velocity dispersion tensor, compressed into 16 bits. 
    int16_t sigmavrad_to_sigmav3d_com; ///< sigmav_rad / sigmav3d, compressed
    int16_t sigmavtan_to_sigmav3d_com; ///< sigmav_tan / sigmav3d, compressed

    int16_t r10_com, r25_com, r33_com, r50_com, r67_com, r75_com, r90_com, r95_com, r98_com; ///<Expressed as ratios of r100, and scaled to 32000 to store as int16s. 
    int16_t sigmar_com[3]; ///<sqrt( Eigenvalues of the moment of inertia tensor ), sorted largest to smallest
    int16_t sigman_com[3]; ///<sqrt( Eigenvalues of the weighted moment of inertia tensor ), sorted largest to smallest
    uint16_t sigmar_eigenvecs_com;  ///<Eigenvectors of the moment of inertia tensor, compressed into 16 bits. Compression format Euler16. 
    uint16_t sigman_eigenvecs_com;  ///<Eigenvectors of the weighted moment of inertia tensor, compressed into 16 bits. Compression format Euler16. 
    int16_t rvcirc_max_com; ///< radius of max velocity, stored as int16 ratio of r100 scaled by 32000.

    // The largest (most massive) subhalo center of mass    
    int16_t sigmavMin_to_sigmav3d_L2com; ///< Min(sigmav_eigenvalue) / sigmav3d, compressed
    int16_t sigmavMax_to_sigmav3d_L2com; ///< Max(sigmav_eigenvalue) / sigmav3d, compressed
    uint16_t sigmav_eigenvecs_L2com;  ///<Eigenvectors of the velocity dispersion tensor, compressed into 16 bits. 
    int16_t sigmavrad_to_sigmav3d_L2com; ///< sigmav_rad / sigmav3d, compressed
    int16_t sigmavtan_to_sigmav3d_L2com; ///< sigmav_tan / sigmav3d, compressed
    int16_t r10_L2com, r25_L2com, r33_L2com, r50_L2com, r67_L2com, r75_L2com, r90_L2com, r95_L2com, r98_L2com;
        ///< Radii of this percentage of mass, relative to L2 center. Expressed as ratios of r100 and compressed to int16. 

    int16_t sigmar_L2com[3];  
    int16_t sigman_L2com[3]; 
    uint16_t sigmar_eigenvecs_L2com;   ///< euler16 format
    uint16_t sigman_eigenvecs_L2com;   ///< euler16 format
    int16_t rvcirc_max_L2com;   ///< radius of max circular velocity, stored as ratio to r100, relative to L2 center

};
