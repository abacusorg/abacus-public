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

uint16_t pack_euler16_eig(double sigma[3], double sigma_vecs[3][3]) {
    int imaj = 0, imin = 0;
    for (int i=0; i<3; i++) {
        if (sigma[i]<sigma[imaj]) imaj = i;
        if (sigma[i]>sigma[imin]) imin = i;
    }
    float major[3], minor[3];
    major[0] = sigma_vecs[imaj][0];
    major[1] = sigma_vecs[imaj][1];
    major[2] = sigma_vecs[imaj][2];
    minor[0] = sigma_vecs[imin][0];
    minor[1] = sigma_vecs[imin][1];
    minor[2] = sigma_vecs[imin][2];
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
    float3 x = com*isize;
    h.x[0] = com.x; h.x[1] = com.y; h.x[2] = com.z; 
    com = double3(0.0, 0.0, 0.0);
    for (int p=0; p<size; p++) com += L1vel[p];
    float3 v = com*isize/ReadState.VelZSpace_to_Canonical;
    assign_to_vector(h.v, v);

    // Find the largest L2 subhalos and the largest COM
    // Groups are already in descending order of multiplicity
    for (int j=0; j<N_LARGEST_SUBHALOS; j++) 
	if (j<L2.ngroups) h.L2cntr_N[j] = L2.groups[j].n; 
	    else h.L2cntr_N[j] = 1;
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
    float3 L2cntr_x = com/h.L2cntr_N[0];
    com = double3(0.0, 0.0, 0.0);
    for (int p=0; p<L2_largest_np; p++) 
    	com += L1vel[start[p].index()];
    float3 L2cntr_v = com/h.L2cntr_N[0]/ReadState.VelZSpace_to_Canonical;
    assign_to_vector(h.L2cntr_v, L2cntr_v);

    // Now we can go through the particles to compute radii and moments
    // We can use L2.d2buffer for scratch space; it is guaranteed to be big enough
    double vxx, vxy, vxz, vyy, vyz, vzz, vrr, vtt;
    double rxx, rxy, rxz, ryy, ryz, rzz;
    double nxx, nxy, nxz, nyy, nyz, nzz;
    float vmax, rvmax;	

    vxx = vxy = vxz = vyy = vyz = vzz = 0.0;
    rxx = rxy = rxz = ryy = ryz = rzz = 0.0;
    nxx = nxy = nxz = nyy = nyz = nzz = 0.0;
    vrr = vtt = 0.0;
    for (int p=0; p<size; p++) {
		posstruct dr = L1pos[p]-x;
		L2.d2buffer[p] = dr.norm2();

		velstruct dv = L1vel[p]/ReadState.VelZSpace_to_Canonical - v;

		vxx += dv.x*dv.x; vxy += dv.x*dv.y; vxz += dv.x*dv.z;
		vyy += dv.y*dv.y; vyz += dv.y*dv.z; vzz += dv.z*dv.z;
		rxx += dr.x*dr.x; rxy += dr.x*dr.y; rxz += dr.x*dr.z;
		ryy += dr.y*dr.y; ryz += dr.y*dr.z; rzz += dr.z*dr.z;
		posstruct n = dr*(1.0/sqrt(dr.norm2()+1e-20));
		float vr = dv.x*n.x+dv.y*n.y+dv.z*n.z;
		vrr += vr*vr; // Accumulate
		vtt += dv.norm2()-vr*vr;  // Accumulate
		nxx += n.x*n.x; nxy += n.x*n.y; nxz += n.x*n.z;
		nyy += n.y*n.y; nyz += n.y*n.z; nzz += n.z*n.z;
    }
    // Normalize by the number of particles
    rxx *= isize; rxy *= isize; rxz *= isize;
    ryy *= isize; ryz *= isize; rzz *= isize;
    vxx *= isize; vxy *= isize; vxz *= isize;
    vyy *= isize; vyz *= isize; vzz *= isize;
    nxx *= isize; nxy *= isize; nxz *= isize;
    nyy *= isize; nyz *= isize; nzz *= isize;
    vrr *= isize; vtt *= isize*0.5;  // Tangential is scaled per dimension

    //beta = 1-sigma_r^2/sigma_t^2;
    
    
    std::sort(L2.d2buffer, L2.d2buffer+size);
    
    h.r100 = sqrt(L2.d2buffer[size-1]); 
    // r10, r25, r50, r67, r75, r90: Expressed as ratios of r100, and scaled to 32000 to store as int16s.   
    h.r10  = lround(sqrt(L2.d2buffer[size/10  ]) / h.r100 * INT16SCALE); 
    h.r25  = lround(sqrt(L2.d2buffer[size/4   ]) / h.r100 * INT16SCALE); 
    h.r33  = lround(sqrt(L2.d2buffer[size/3   ]) / h.r100 * INT16SCALE); 
    h.r50  = lround(sqrt(L2.d2buffer[size/2   ]) / h.r100 * INT16SCALE); 
    h.r67  = lround(sqrt(L2.d2buffer[size*2/3 ]) / h.r100 * INT16SCALE); 
    h.r75  = lround(sqrt(L2.d2buffer[size*3/4 ]) / h.r100 * INT16SCALE); 
    h.r90  = lround(sqrt(L2.d2buffer[size*9/10]) / h.r100 * INT16SCALE); 
	
	double sigmav[3], sigmar[3], sigman[3]; 
	double sigmav_vecs[3][3]; 
	double sigmar_vecs[3][3]; 
	double sigman_vecs[3][3]; 

	FindEigensystem(vxx, vxy, vxz, vyy, vyz, vzz, sigmav, (double * )sigmav_vecs);
    FindEigensystem(rxx, rxy, rxz, ryy, ryz, rzz, sigmar, (double * )sigmar_vecs);
    FindEigensystem(nxx, nxy, nxz, nyy, nyz, nzz, sigman, (double * )sigman_vecs);

    h.sigmar_eigenvecs = pack_euler16_eig(sigmar, sigmar_vecs);
    h.sigmav_eigenvecs = pack_euler16_eig(sigmav, sigmav_vecs);
    h.sigman_eigenvecs = pack_euler16_eig(sigman, sigman_vecs);

    h.sigmav3d = sqrt(sigmav[0] + sigmav[1] + sigmav[2]); 
    h.sigmavMin_to_sigmav3d = lround( sqrt(sigmav[2])  / h.sigmav3d * INT16SCALE ); 
    h.sigmavMax_to_sigmav3d = lround( sqrt(sigmav[0])  / h.sigmav3d * INT16SCALE ); 
    // h.sigmav3d_to_sigmavMaj = lround( h.sigmav3d/sqrt(sigmav[0]) * INT16SCALE );
    h.sigmavtan_to_sigmav3d = lround(sqrt(vtt)/h.sigmav3d * INT16SCALE ); 
    h.sigmavrad_to_sigmav3d = lround(sqrt(vrr)/h.sigmav3d * INT16SCALE ); 

    for(int i = 0; i < 3; i++) h.sigmar[i] = lround(sqrt(sigmar[i]) / h.r100 * INT16SCALE );
    for(int i = 0; i < 3; i++) h.sigman[i] = lround(sqrt(sigman[i]) * INT16SCALE );

    
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
    for (int p=size/10; p<size; p++) {
		float v4 = p*p/L2.d2buffer[p];
		if (v4>vmax) { vmax = v4; rvmax = L2.d2buffer[p]; }
    }
    h.rvcirc_max = lround(sqrt(rvmax) / h.r100 * INT16SCALE );    // Get to radial units and compress into int16. 
    float GMpart = 3*P.Omega_M*pow(100*ReadState.BoxSizeHMpc,2)/(8*M_PI*P.np*ReadState.ScaleFactor);
    h.vcirc_max = sqrt(GMpart*sqrt(vmax))/ReadState.VelZSpace_to_kms;  // This is sqrt(G*M_particle*N/R).

    // Repeat this, finding moments and radii around the largest subhalo COM
    vxx = vxy = vxz = vyy = vyz = vzz = 0.0;
    rxx = rxy = rxz = ryy = ryz = rzz = 0.0;
    nxx = nxy = nxz = nyy = nyz = nzz = 0.0;
    vrr = vtt = 0.0;
    for (int p=0; p<size; p++) {
		posstruct dr = L1pos[p]-L2cntr_x;
		L2.d2buffer[p] = dr.norm2();
		velstruct dv = L1vel[p]/ReadState.VelZSpace_to_Canonical - L2cntr_v;
		vxx += dv.x*dv.x; vxy += dv.x*dv.y; vxz += dv.x*dv.z;
		vyy += dv.y*dv.y; vyz += dv.y*dv.z; vzz += dv.z*dv.z;
		rxx += dr.x*dr.x; rxy += dr.x*dr.y; rxz += dr.x*dr.z;
		ryy += dr.y*dr.y; ryz += dr.y*dr.z; rzz += dr.z*dr.z;
		posstruct n = dr*(1.0/sqrt(dr.norm2()+1e-20));
		float vr = dv.x*n.x+dv.y*n.y+dv.z*n.z;
		vrr += vr*vr; // Accumulate
		vtt += dv.norm2()-vr*vr;  // Accumulate
		nxx += n.x*n.x; nxy += n.x*n.y; nxz += n.x*n.z;
		nyy += n.y*n.y; nyz += n.y*n.z; nzz += n.z*n.z;
    }
    // Normalize by the number of particles
    rxx *= isize; rxy *= isize; rxz *= isize;
    ryy *= isize; ryz *= isize; rzz *= isize;
    vxx *= isize; vxy *= isize; vxz *= isize;
    vyy *= isize; vyz *= isize; vzz *= isize;
    nxx *= isize; nxy *= isize; nxz *= isize;
    nyy *= isize; nyz *= isize; nzz *= isize;
    vrr *= isize; vtt *= isize*0.5;  // Tangential is scaled per dimension
    std::sort(L2.d2buffer, L2.d2buffer+size);

    // B.H.
    vtt *= 0.5;
    
    h.L2cntr_r100 = sqrt(L2.d2buffer[size-1]);   
    // r10, r25, r50, r67, r75, r90 relative to largest L2 center: Expressed as ratios of r100, and scaled to 32000 to store as int16s. 
    h.L2cntr_r10  = lround(sqrt(L2.d2buffer[size/10  ]) / h.L2cntr_r100 * INT16SCALE); 
    h.L2cntr_r25  = lround(sqrt(L2.d2buffer[size/4   ]) / h.L2cntr_r100 * INT16SCALE); 
    h.L2cntr_r33  = lround(sqrt(L2.d2buffer[size/3   ]) / h.L2cntr_r100 * INT16SCALE); 
    h.L2cntr_r50  = lround(sqrt(L2.d2buffer[size/2   ]) / h.L2cntr_r100 * INT16SCALE); 
    h.L2cntr_r67  = lround(sqrt(L2.d2buffer[size*2/3 ]) / h.L2cntr_r100 * INT16SCALE); 
    h.L2cntr_r75  = lround(sqrt(L2.d2buffer[size*3/4 ]) / h.L2cntr_r100 * INT16SCALE); 
    h.L2cntr_r90  = lround(sqrt(L2.d2buffer[size*9/10]) / h.L2cntr_r100 * INT16SCALE); 


    FindEigensystem(vxx, vxy, vxz, vyy, vyz, vzz, sigmav, (double * )sigmav_vecs);
    FindEigensystem(rxx, rxy, rxz, ryy, ryz, rzz, sigmar, (double * )sigmar_vecs);
    FindEigensystem(nxx, nxy, nxz, nyy, nyz, nzz, sigman, (double * )sigman_vecs);

    h.L2cntr_sigmar_eigenvecs = pack_euler16_eig(sigmar, sigmar_vecs);
    h.L2cntr_sigmav_eigenvecs = pack_euler16_eig(sigmav, sigmav_vecs);
    h.L2cntr_sigman_eigenvecs = pack_euler16_eig(sigman, sigman_vecs);

    h.L2cntr_sigmav3d = sqrt(sigmav[0] + sigmav[1] + sigmav[2]);
    h.L2cntr_sigmavMin_to_sigmav3d = lround( sqrt(sigmav[2])  / h.L2cntr_sigmav3d * INT16SCALE ); 
    h.L2cntr_sigmavMax_to_sigmav3d = lround( sqrt(sigmav[0])  / h.L2cntr_sigmav3d * INT16SCALE ); 
    // h.L2cntr_sigmav3d_to_sigmavMaj = lround( h.L2cntr_sigmav3d/sqrt(sigmav[0]) * INT16SCALE );
    h.L2cntr_sigmavtan_to_sigmav3d = lround(sqrt(vtt)/h.L2cntr_sigmav3d * INT16SCALE ); 
    h.L2cntr_sigmavrad_to_sigmav3d = lround(sqrt(vrr)/h.L2cntr_sigmav3d * INT16SCALE ); 
	
    for(int i = 0; i < 3; i++) h.L2cntr_sigmar[i] = lround(sqrt(sigmar[i]) / h.r100 * INT16SCALE );
    for(int i = 0; i < 3; i++) h.L2cntr_sigman[i] = lround(sqrt(sigman[i]) * INT16SCALE );

    // We search for the max of vcirc, which is proportional to sqrt(G*M/R).
    // The 4th power of that is proportional to N^2/R^2.
    vmax = 0.0;
    for (int p=size/10; p<size; p++) {
	float v4 = p*p/L2.d2buffer[p];
	if (v4>vmax) { vmax = v4; rvmax = L2.d2buffer[p]; }
    }
    h.L2cntr_rvcirc_max = lround(sqrt(rvmax) / h.L2cntr_r100 * INT16SCALE );    // Get to radial units and compress into int16. 
    h.L2cntr_vcirc_max = sqrt(GMpart*sqrt(vmax))/ReadState.VelZSpace_to_kms;  // This is sqrt(N/R).
	
    x += offset; 
    x = WrapPosition(x);
    L2cntr_x += offset; L2cntr_x = WrapPosition(L2cntr_x);
    assign_to_vector(h.x, x);
    assign_to_vector(h.L2cntr_x, L2cntr_x);
 
    return h;
};



