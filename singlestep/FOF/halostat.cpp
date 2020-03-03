#include "halostat.hh"
#include "sym3eigenval.cpp"
#include "euler16.c"

inline float3 WrapPosition(float3 a) {
    while (a.x> 0.5) a.x-=1.0;
    while (a.x<-0.5) a.x+=1.0;
    while (a.y> 0.5) a.y-=1.0;
    while (a.y<-0.5) a.y+=1.0;
    while (a.z> 0.5) a.z-=1.0;
    while (a.z<-0.5) a.z+=1.0;
    return a;
}

/// This sorts the eigenvalues and vectors so that the largest is first.
inline void sort_eig(double sigma[3], double sigma_vecs[3][3]) {
    // Remember that the eigenvectors are in columns, so [0..2][j] is the jth.
    #define SWAP_EIG(a,b) std::swap(sigma[a],sigma[b]); \
        std::swap(sigma_vecs[0][a],sigma_vecs[0][b]); \
        std::swap(sigma_vecs[1][a],sigma_vecs[1][b]); \
        std::swap(sigma_vecs[2][a],sigma_vecs[2][b]); 

    if (sigma[0]<sigma[1]) { SWAP_EIG(0,1) }   // Now 0>1
    if (sigma[1]<sigma[2]) { SWAP_EIG(1,2) }   // Now 1>2
    if (sigma[0]<sigma[1]) { SWAP_EIG(0,1) }   // Now 0>1
    return;
    #undef SWAP_EIG
}

/// In addition to returning the eigenvector euler16 packing,
/// this code also sorts the inputs into descending order of eigenvalues.
uint16_t pack_euler16_eig(double sigma[3], double sigma_vecs[3][3]) {
    sort_eig(sigma, sigma_vecs);
    int imaj = 0, imin = 2;
    float major[3], minor[3];
    major[0] = sigma_vecs[0][imaj];
    major[1] = sigma_vecs[1][imaj];
    major[2] = sigma_vecs[2][imaj];
    minor[0] = sigma_vecs[0][imin];
    minor[1] = sigma_vecs[1][imin];
    minor[2] = sigma_vecs[2][imin];
    return pack_euler16(major, minor);
}

#define assign_to_vector(a,b) { a[0] = b.x; a[1] = b.y; a[2] = b.z; }
    /** Fill a HaloStat object and return it.

    We are given the L1 particles as pos/vel/aux from [0,size).
    These are in the original order; we don't care.
    We are also given the L2 FOF results class.
    The particle positions have already been wrapped to the first cell;
    we will wrap to global coords, using offset
    Velocities will be converted to our usual velocity output
    convention of unit-box redshift-space displacements
    */

HaloStat ComputeStats(int size, 
	posstruct *L1pos, velstruct *L1vel, auxstruct *L1aux, 
	#ifdef SPHERICAL_OVERDENSITY
	    SOcell &L2, 
	#else
	    FOFcell &L2, 
	#endif
	posstruct offset) {
    HaloStat h;
    const int16_t INT16SCALE = 32000;
    h.N = size;
    float isize = 1.0/size;

    
    // Compute the center of mass
    double3 com = double3(0.0, 0.0, 0.0);
    for (int p=0; p<size; p++) com += L1pos[p];
    float3 x_com = com*isize;
    com = double3(0.0, 0.0, 0.0);
    for (int p=0; p<size; p++) com += L1vel[p];
    float3 v_com = com*isize/ReadState.VelZSpace_to_Canonical;
    assign_to_vector(h.v_com, v_com);

    // Find the largest L2 subhalos and the largest COM
    // Groups are already in descending order of multiplicity
    for (int j=0; j<N_LARGEST_SUBHALOS; j++) 
	if (j<L2.ngroups) h.L2_N[j] = L2.groups[j].n; 
	    else h.L2_N[j] = 1;
	    // Strictly speaking, there might not be any singlet particles,
	    // but the far more likely case is that there are.
	    // So this output could be in error, i.e., 0 and 1 are not distinguished

	// We are always guaranteed to have some L2 singlet particles,
	// but singlet particles do not have groups!
	int L2_largest_np = 1;
    FOFparticle *start = L2.p;
	if (L2.ngroups > 0){
		start = L2.p + L2.groups[0].start;
		L2_largest_np = L2.groups[0].n;
	}

    com = double3(0.0, 0.0, 0.0);
    for (int p=0; p<L2_largest_np; p++) 
    	com += L1pos[start[p].index()];
    float3 x_L2com = com/h.L2_N[0];
    com = double3(0.0, 0.0, 0.0);
    for (int p=0; p<L2_largest_np; p++) 
    	com += L1vel[start[p].index()];
    float3 v_L2com = com/h.L2_N[0]/ReadState.VelZSpace_to_Canonical;
    assign_to_vector(h.v_L2com, v_L2com);

    // Now we can go through the particles to compute radii and moments. We can use
    // L2.d2_active and L2.d2buffer for scratch space; they are guaranteed to be big enough
    for (int p=0; p<size; p++) {
                posstruct dr = L1pos[p]-x_com;
		L2.d2buffer[p] = dr.norm2();
		L2.d2_active[p] = L2.d2buffer[p];
    }	

    std::sort(L2.d2_active, L2.d2_active+size);

    float r90sq = L2.d2_active[size*9/10]; // radius with respect to which we compute the second moments
    float r50sq = L2.d2_active[size/2]; // radius with respect to which we compute the second moments
    float isize_r90 = 1./(size*9/10); // number of elements within that radius
    float isize_r50 = 1./(size/2); // number of elements within that radius
    
    h.r100_com = sqrt(L2.d2_active[size-1]); 
    // r10, r25, r50, r67, r75, r90, r95, r98: Expressed as ratios of r100, and scaled to 32000 to store as int16s.   
    h.r10_com  = lround(sqrt(L2.d2_active[size/10  ]) / h.r100_com * INT16SCALE); 
    h.r25_com  = lround(sqrt(L2.d2_active[size/4   ]) / h.r100_com * INT16SCALE); 
    h.r33_com  = lround(sqrt(L2.d2_active[size/3   ]) / h.r100_com * INT16SCALE); 
    h.r50_com  = lround(sqrt(L2.d2_active[size/2   ]) / h.r100_com * INT16SCALE); 
    h.r67_com  = lround(sqrt(L2.d2_active[size*2/3 ]) / h.r100_com * INT16SCALE); 
    h.r75_com  = lround(sqrt(L2.d2_active[size*3/4 ]) / h.r100_com * INT16SCALE); 
    h.r90_com  = lround(sqrt(r90sq) / h.r100_com * INT16SCALE);
    h.r95_com  = lround(sqrt(L2.d2_active[size*19/20]) / h.r100_com * INT16SCALE); 	
    h.r98_com  = lround(sqrt(L2.d2_active[size*49/50]) / h.r100_com * INT16SCALE);
    
    double vxx, vxy, vxz, vyy, vyz, vzz, vrr, vtt;
    double rxx, rxy, rxz, ryy, ryz, rzz;
    double nxx, nxy, nxz, nyy, nyz, nzz;
    float vmax, rvmax; 
    double vmean, vmean_r50, vsq_r50;

    vxx = vxy = vxz = vyy = vyz = vzz = 0.0;
    rxx = rxy = rxz = ryy = ryz = rzz = 0.0;
    nxx = nxy = nxz = nyy = nyz = nzz = 0.0;
    vrr = vtt = vmean = vmean_r50 = vsq_r50 = 0.0;
    for (int p=0; p<size; p++) {
        if (L2.d2buffer[p] < r90sq) {
            posstruct dr = L1pos[p]-x_com;
            velstruct dv = L1vel[p]/ReadState.VelZSpace_to_Canonical-v_com;

            vxx += dv.x*dv.x; vxy += dv.x*dv.y; vxz += dv.x*dv.z;
            vyy += dv.y*dv.y; vyz += dv.y*dv.z; vzz += dv.z*dv.z;
            rxx += dr.x*dr.x; rxy += dr.x*dr.y; rxz += dr.x*dr.z;
            ryy += dr.y*dr.y; ryz += dr.y*dr.z; rzz += dr.z*dr.z;
            posstruct n = dr*(1.0/sqrt(dr.norm2()+1e-20));
            float vr = dv.x*n.x+dv.y*n.y+dv.z*n.z;
            vrr += vr*vr; // Accumulate
            float vnorm2 = dv.norm2();
            vtt += vnorm2-vr*vr;  // Accumulate
            nxx += n.x*n.x; nxy += n.x*n.y; nxz += n.x*n.z;
            nyy += n.y*n.y; nyz += n.y*n.z; nzz += n.z*n.z;
            float vnorm = sqrt(vnorm2);
            vmean += vnorm;
            if (L2.d2buffer[p]<r50sq) {
                vsq_r50 += vnorm2;
                vmean_r50 += vnorm;
            }
        }
    }
    // Normalize by the number of particles
    rxx *= isize_r90; rxy *= isize_r90; rxz *= isize_r90;
    ryy *= isize_r90; ryz *= isize_r90; rzz *= isize_r90;
    vxx *= isize_r90; vxy *= isize_r90; vxz *= isize_r90;
    vyy *= isize_r90; vyz *= isize_r90; vzz *= isize_r90;
    nxx *= isize_r90; nxy *= isize_r90; nxz *= isize_r90;
    nyy *= isize_r90; nyz *= isize_r90; nzz *= isize_r90;
    vrr *= isize_r90; vtt *= isize_r90*0.5; vmean *= isize_r90; // Tangential is scaled per dimension
    vmean_r50 *= isize_r50; vsq_r50 *= isize_r50;
     
    double sigmav[3], sigmar[3], sigman[3]; 
    double sigmav_vecs[3][3]; 
    double sigmar_vecs[3][3]; 
    double sigman_vecs[3][3]; 

    FindEigensystem(vxx, vxy, vxz, vyy, vyz, vzz, sigmav, (double * )sigmav_vecs);
    FindEigensystem(rxx, rxy, rxz, ryy, ryz, rzz, sigmar, (double * )sigmar_vecs);
    FindEigensystem(nxx, nxy, nxz, nyy, nyz, nzz, sigman, (double * )sigman_vecs);
    
    h.sigmar_eigenvecs_com = pack_euler16_eig(sigmar, sigmar_vecs);
    h.sigmav_eigenvecs_com = pack_euler16_eig(sigmav, sigmav_vecs);
    h.sigman_eigenvecs_com = pack_euler16_eig(sigman, sigman_vecs);
    // The eigenvalues are now sorted in descending order
	  
    h.sigmav3d_com = sqrt(sigmav[0] + sigmav[1] + sigmav[2]); 
    h.sigmavMin_to_sigmav3d_com = lround( sqrt(sigmav[2])  / h.sigmav3d_com * INT16SCALE ); 
    h.sigmavMax_to_sigmav3d_com = lround( sqrt(sigmav[0])  / h.sigmav3d_com * INT16SCALE ); 
    h.sigmavtan_to_sigmav3d_com = lround(sqrt(vtt)/h.sigmav3d_com * INT16SCALE ); 
    h.sigmavrad_to_sigmav3d_com = lround(sqrt(vrr)/h.sigmav3d_com * INT16SCALE );
    h.meanSpeed_com = vmean;
    h.meanSpeed_r50_com = vmean_r50;
    h.sigmav3d_r50_com  = vsq_r50;

    for(int i = 0; i < 3; i++) h.sigmar_com[i] = lround(sqrt(sigmar[i]) / h.r100_com * INT16SCALE );
    for(int i = 0; i < 3; i++) h.sigman_com[i] = lround(sqrt(sigman[i]) * INT16SCALE );
    
    
#ifdef SPHERICAL_OVERDENSITY
    h.SO_L2max_central_particle[0] = L2.p[0].x;
    h.SO_L2max_central_particle[1] = L2.p[0].y;
    h.SO_L2max_central_particle[2] = L2.p[0].z;
    h.SO_L2max_central_particle[3] = L2.p[0].n;	
    h.SO_central_density  = L2.density[0]; 
	//!!!h.SO_radius           = sqrt(L2.halo_thresh2); 
#endif 	
    // We search for the max of vcirc, which is proportional to sqrt(G*M/R).
    // The 4th power of that is proportional to N^2/R^2.
    vmax = 0.0;
    for (int p=(size<1000?size/10:100); p<size; p++) {
		float v4 = (float)p*p/L2.d2_active[p];
		if (v4>vmax) { vmax = v4; rvmax = L2.d2_active[p]; }
    }
    h.rvcirc_max_com = lround(sqrt(rvmax) / h.r100_com * INT16SCALE );    // Get to radial units and compress into int16. 
    float GMpart = 3*(P.Omega_M-P.Omega_Smooth)*pow(100*ReadState.BoxSizeHMpc,2)/(8*M_PI*P.np*ReadState.ScaleFactor);
    h.vcirc_max_com = sqrt(GMpart*sqrt(vmax))/ReadState.VelZSpace_to_kms;  // This is sqrt(G*M_particle*N/R).
    
    // Repeat this, finding moments and radii around the largest subhalo COM
    for (int p=0; p<size; p++) {
                posstruct dr = L1pos[p]-x_L2com;
		L2.d2buffer[p] = dr.norm2();
		L2.d2_active[p] = L2.d2buffer[p];
    }
    
    std::sort(L2.d2_active, L2.d2_active+size);

    r90sq = L2.d2_active[size*9/10]; // radius within which we compute the second moments
    r50sq = L2.d2_active[size/2]; // radius within which we compute the second moments
    isize_r90 = 1./(size*9/10); // number of elements within that radius
    isize_r50 = 1./(size/2); // number of elements within that radius
    
    h.r100_L2com = sqrt(L2.d2_active[size-1]); 
    // r10, r25, r50, r67, r75, r90, r95, r98 wrt COM of the largest L2. Expressed as ratios of r100, and scaled to 32000 to store as int16s.   
    h.r10_L2com  = lround(sqrt(L2.d2_active[size/10  ]) / h.r100_L2com * INT16SCALE); 
    h.r25_L2com  = lround(sqrt(L2.d2_active[size/4   ]) / h.r100_L2com * INT16SCALE); 
    h.r33_L2com  = lround(sqrt(L2.d2_active[size/3   ]) / h.r100_L2com * INT16SCALE); 
    h.r50_L2com  = lround(sqrt(L2.d2_active[size/2   ]) / h.r100_L2com * INT16SCALE); 
    h.r67_L2com  = lround(sqrt(L2.d2_active[size*2/3 ]) / h.r100_L2com * INT16SCALE); 
    h.r75_L2com  = lround(sqrt(L2.d2_active[size*3/4 ]) / h.r100_L2com * INT16SCALE); 
    h.r90_L2com  = lround(sqrt(r90sq) / h.r100_L2com * INT16SCALE);
    h.r95_L2com  = lround(sqrt(L2.d2_active[size*19/20]) / h.r100_L2com * INT16SCALE); 	
    h.r98_L2com  = lround(sqrt(L2.d2_active[size*49/50]) / h.r100_L2com * INT16SCALE);

    vxx = vxy = vxz = vyy = vyz = vzz = 0.0;
    rxx = rxy = rxz = ryy = ryz = rzz = 0.0;
    nxx = nxy = nxz = nyy = nyz = nzz = 0.0;
    vrr = vtt = vmean = vsq_r50 = vmean_r50 = 0.0;
    for (int p=0; p<size; p++) {
        if (L2.d2buffer[p] < r90sq) {
            posstruct dr = L1pos[p]-x_L2com;
            velstruct dv = L1vel[p]/ReadState.VelZSpace_to_Canonical-v_L2com;

            vxx += dv.x*dv.x; vxy += dv.x*dv.y; vxz += dv.x*dv.z;
            vyy += dv.y*dv.y; vyz += dv.y*dv.z; vzz += dv.z*dv.z;
            rxx += dr.x*dr.x; rxy += dr.x*dr.y; rxz += dr.x*dr.z;
            ryy += dr.y*dr.y; ryz += dr.y*dr.z; rzz += dr.z*dr.z;
            posstruct n = dr*(1.0/sqrt(dr.norm2()+1e-20));
            float vr = dv.x*n.x+dv.y*n.y+dv.z*n.z;
            vrr += vr*vr; // Accumulate
            float vnorm2 = dv.norm2();
            float vnorm = sqrt(vnorm2);
            vtt += vnorm2-vr*vr;  // Accumulate
            nxx += n.x*n.x; nxy += n.x*n.y; nxz += n.x*n.z;
            nyy += n.y*n.y; nyz += n.y*n.z; nzz += n.z*n.z;
            vmean += vnorm;
            if (L2.d2buffer[p]<r50sq) {
                vsq_r50 += vnorm2;
                vmean_r50 += vnorm;
            }
        }
    }
    // Normalize by the number of particles
    rxx *= isize_r90; rxy *= isize_r90; rxz *= isize_r90;
    ryy *= isize_r90; ryz *= isize_r90; rzz *= isize_r90;
    vxx *= isize_r90; vxy *= isize_r90; vxz *= isize_r90;
    vyy *= isize_r90; vyz *= isize_r90; vzz *= isize_r90;
    nxx *= isize_r90; nxy *= isize_r90; nxz *= isize_r90;
    nyy *= isize_r90; nyz *= isize_r90; nzz *= isize_r90;
    vrr *= isize_r90; vtt *= isize_r90*0.5; vmean *= isize_r90; // Tangential is scaled per dimension
    vmean_r50 *= isize_r50; vsq_r50 *= isize_r50;

    FindEigensystem(vxx, vxy, vxz, vyy, vyz, vzz, sigmav, (double * )sigmav_vecs);
    FindEigensystem(rxx, rxy, rxz, ryy, ryz, rzz, sigmar, (double * )sigmar_vecs);
    FindEigensystem(nxx, nxy, nxz, nyy, nyz, nzz, sigman, (double * )sigman_vecs);
    
    h.sigmar_eigenvecs_L2com = pack_euler16_eig(sigmar, sigmar_vecs);
    h.sigmav_eigenvecs_L2com = pack_euler16_eig(sigmav, sigmav_vecs);
    h.sigman_eigenvecs_L2com = pack_euler16_eig(sigman, sigman_vecs);

    // The eigenvalues are now sorted in descending order

    h.sigmav3d_L2com = sqrt(sigmav[0] + sigmav[1] + sigmav[2]); 
    h.sigmavMin_to_sigmav3d_L2com = lround( sqrt(sigmav[2])  / h.sigmav3d_L2com * INT16SCALE ); 
    h.sigmavMax_to_sigmav3d_L2com = lround( sqrt(sigmav[0])  / h.sigmav3d_L2com * INT16SCALE ); 
    h.sigmavtan_to_sigmav3d_L2com = lround(sqrt(vtt)/h.sigmav3d_L2com * INT16SCALE ); 
    h.sigmavrad_to_sigmav3d_L2com = lround(sqrt(vrr)/h.sigmav3d_L2com * INT16SCALE ); 
    h.meanSpeed_L2com = vmean;
    h.meanSpeed_r50_L2com = vmean_r50;
    h.sigmav3d_r50_L2com  = vsq_r50;
    
    for(int i = 0; i < 3; i++) h.sigmar_L2com[i] = lround(sqrt(sigmar[i]) / h.r100_L2com * INT16SCALE );
    for(int i = 0; i < 3; i++) h.sigman_L2com[i] = lround(sqrt(sigman[i]) * INT16SCALE );
    
    // We search for the max of vcirc, which is proportional to sqrt(G*M/R).
    // The 4th power of that is proportional to N^2/R^2.
    vmax = 0.0;
    for (int p=(size<1000?size/10:100); p<size; p++) {
	float v4 = (float)p*p/L2.d2_active[p];
	if (v4>vmax) { vmax = v4; rvmax = L2.d2_active[p]; }
    }
    h.rvcirc_max_L2com = lround(sqrt(rvmax) / h.r100_L2com * INT16SCALE );    // Get to radial units and compress into int16. 
    h.vcirc_max_L2com = sqrt(GMpart*sqrt(vmax))/ReadState.VelZSpace_to_kms;  // This is sqrt(N/R).
	
    x_com += offset; x_com = WrapPosition(x_com);
    x_L2com += offset; x_L2com = WrapPosition(x_L2com);
    assign_to_vector(h.x_com, x_com);
    assign_to_vector(h.x_L2com, x_L2com);
 
    return h;
};



