#include "halostat.hh"
#include "sym3eigenval.cpp"

inline float3 WrapPosition(float3 a) {
    while (a.x> 0.5) a.x-=1.0;
    while (a.x<-0.5) a.x+=1.0;
    while (a.y> 0.5) a.y-=1.0;
    while (a.y<-0.5) a.y+=1.0;
    while (a.z> 0.5) a.z-=1.0;
    while (a.z<-0.5) a.z+=1.0;
    return a;
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

    h.N = size;
    // Compute the center of mass
    double3 com = double3(0.0, 0.0, 0.0);
    for (int p=0; p<size; p++) com += L1pos[p];
    float3 x = com/size;
    h.x[0] = com.x; h.x[1] = com.y; h.x[2] = com.z; 
    com = double3(0.0, 0.0, 0.0);
    for (int p=0; p<size; p++) com += L1vel[p];
    float3 v = com/size/ReadState.VelZSpace_to_Canonical;
    assign_to_vector(h.v, v);

    // Find the largest L2 subhalos and the largest COM
    // Groups are already in descending order of multiplicity
    for (int j=0; j<N_LARGEST_SUBHALOS; j++) 
	if (j<L2.ngroups) h.subhalo_N[j] = L2.groups[j].n; 
	    else h.subhalo_N[j] = 1;
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
    float3 subhalo_x = com/h.subhalo_N[0];
    com = double3(0.0, 0.0, 0.0);
    for (int p=0; p<L2_largest_np; p++) 
    	com += L1vel[start[p].index()];
    float3 subhalo_v = com/h.subhalo_N[0]/ReadState.VelZSpace_to_Canonical;
    assign_to_vector(h.subhalo_v, subhalo_v);

    // Now we can go through the particles to compute radii and moments
    // We can use L2.d2buffer for scratch space; it is guaranteed to be big enough
    double vxx, vxy, vxz, vyy, vyz, vzz;
    float vmax, rvmax;
	float dmean = 0.0; 
	
    vxx = vxy = vxz = vyy = vyz = vzz = 0.0;
    for (int p=0; p<size; p++) {
		posstruct dx = L1pos[p]-x;
		L2.d2buffer[p] = dx.norm2();
		d2mean += sqrt(L2.d2buffer[p]);    //NAM doing this partially defeats the point of having the d2buffer-- here we do a sqrt for every particle. Is there a way to avoid this? 
		velstruct dv = L1vel[p]/ReadState.VelZSpace_to_Canonical - v;
		vxx += dv.x*dv.x; vxy += dv.x*dv.y; vxz += dv.x*dv.z;
		vyy += dv.y*dv.y; vyz += dv.y*dv.z; vzz += dv.z*dv.z;
    }
	dmean /= (float)size; 
	
    FindEigenvalues(vxx, vxy, vxz, vyy, vyz, vzz, h.sigmav);
    for(int i = 0; i < 3; i++) h.sigmav[i] = sqrt(h.sigmav[i]);
    std::sort(L2.d2buffer, L2.d2buffer+size);
	
	int16_t INT16SCALE = 32000; 
	h.r100 = sqrt(L2.d2buffer[-1]; 
	// r10, r25, r50, r67, r75, r90: Expressed as ratios of r100, and scaled to 32000 to store as int16s. 	
	h.r10  = trunc(sqrt(L2.d2buffer[size/10  ]) / h.r100 * INT16SCALE); 
	h.r25  = trunc(sqrt(L2.d2buffer[size/4   ]) / h.r100 * INT16SCALE); 
	h.r50  = trunc(sqrt(L2.d2buffer[size/2   ]) / h.r100 * INT16SCALE); 
	h.r67  = trunc(sqrt(L2.d2buffer[size*2/3 ]) / h.r100 * INT16SCALE); 
	h.r75  = trunc(sqrt(L2.d2buffer[size*3/4 ]) / h.r100 * INT16SCALE); 
	h.r90  = trunc(sqrt(L2.d2buffer[size*9/10]) / h.r100 * INT16SCALE); 
	
	float sigmar_accum = 0.0; 
	for (int p=0; p<size; p++){
		sigmar_accum += pow(sqrt(L2.d2buffer[p]) - dmean, 2.0);
		
		//this is redundant but avoids calculating this value every iteration of the loop. think if there's a better way to do this. 
		if      (p == size/10  ) h.sigmar10  = trunc(sigmar_accum / h.r100 * INT16SCALE); 
		else if (p == size/4   ) h.sigmar25  = trunc(sigmar_accum / h.r100 * INT16SCALE);
		else if (p == size/2   ) h.sigmar50  = trunc(sigmar_accum / h.r100 * INT16SCALE);
		else if (p == size*2/3 ) h.sigmar67  = trunc(sigmar_accum / h.r100 * INT16SCALE);
		else if (p == size*3/4 ) h.sigmar75  = trunc(sigmar_accum / h.r100 * INT16SCALE);
		else if (p == size*9/10) h.sigmar90  = trunc(sigmar_accum / h.r100 * INT16SCALE);
		else if (p == size - 1 ) h.sigmar100 = trunc(sigmar_accum / h.r100 * INT16SCALE);
		
	}

	
#ifdef SPHERICAL_OVERDENSITY
	h.SO_central_particle = L2.particles[0];
	h.SO_central_density  = L2.density[0]; 
	h.SO_radius           = sqrt(L2.d2buffer[-1]); 
#endif 
	
    // We search for the max of vcirc, which is proportional to sqrt(G*M/R).
    // The 4th power of that is proportional to N^2/R^2.
    vmax = 0.0;
    for (int p=size/10; p<size; p++) {
		float v4 = p*p/L2.d2buffer[p];
		if (v4>vmax) { vmax = v4; rvmax = L2.d2buffer[p]; }
    }
    h.rvcirc_max = trunc(sqrt(rvmax) / r100 * INT16SCALE );    // Get to radial units and compress into int16. 
    float GMpart = 3*P.Omega_M*pow(100*ReadState.BoxSizeHMpc,2)/(8*M_PI*P.np*ReadState.ScaleFactor);
    h.vcirc_max = sqrt(GMpart*sqrt(vmax))/ReadState.VelZSpace_to_kms;  // This is sqrt(G*M_particle*N/R).

    // Repeat this, finding moments and radii around the largest subhalo COM
	dmean = 0.0; 
    vxx = vxy = vxz = vyy = vyz = vzz = 0.0;
    for (int p=0; p<size; p++) {
		posstruct dx = L1pos[p]-subhalo_x;
		L2.d2buffer[p] = dx.norm2();
		dmean += sqrt(L2.d2buffer[p]); 
		velstruct dv = L1vel[p]/ReadState.VelZSpace_to_Canonical - subhalo_v;
		vxx += dv.x*dv.x; vxy += dv.x*dv.y; vxz += dv.x*dv.z;
		vyy += dv.y*dv.y; vyz += dv.y*dv.z; vzz += dv.z*dv.z;
    }
	dmean /= (float)size; 
	
    FindEigenvalues(vxx, vxy, vxz, vyy, vyz, vzz, h.subhalo_sigmav);
    for(int i = 0; i < 3; i++) h.subhalo_sigmav[i] = sqrt(h.subhalo_sigmav[i]);
    std::sort(L2.d2buffer, L2.d2buffer+size);

	h.subhalo_r100 = sqrt(L2.d2buffer[-1]);   
	// r10, r25, r50, r67, r75, r90 relative to largest L2 center: Expressed as ratios of r100, and scaled to 32000 to store as int16s. 
	h.subhalo_r10  = trunc(sqrt(L2.d2buffer[size/10  ]) / h.subhalo_r100 * INT16SCALE); 
	h.subhalo_r25  = trunc(sqrt(L2.d2buffer[size/4   ]) / h.subhalo_r100 * INT16SCALE); 
	h.subhalo_r50  = trunc(sqrt(L2.d2buffer[size/2   ]) / h.subhalo_r100 * INT16SCALE); 
	h.subhalo_r67  = trunc(sqrt(L2.d2buffer[size*2/3 ]) / h.subhalo_r100 * INT16SCALE); 
	h.subhalo_r75  = trunc(sqrt(L2.d2buffer[size*3/4 ]) / h.subhalo_r100 * INT16SCALE); 
	h.subhalo_r90  = trunc(sqrt(L2.d2buffer[size*9/10]) / h.subhalo_r100 * INT16SCALE); 
	
	sigmar_accum = 0.0; 
	for (int p=0; p<size; p++){
		sigmar_accum += pow(sqrt(L2.d2buffer[p]) - dmean, 2.0);
		
		if      (p == size/10  ) h.subhalo_sigmar10  = trunc(sigmar_accum / h.r100 * INT16SCALE);
		else if (p == size/4   ) h.subhalo_sigmar25  = trunc(sigmar_accum / h.r100 * INT16SCALE);
		else if (p == size/2   ) h.subhalo_sigmar50  = trunc(sigmar_accum / h.r100 * INT16SCALE);
		else if (p == size*2/3 ) h.subhalo_sigmar67  = trunc(sigmar_accum / h.r100 * INT16SCALE);
		else if (p == size*3/4 ) h.subhalo_sigmar75  = trunc(sigmar_accum / h.r100 * INT16SCALE);
		else if (p == size*9/10) h.subhalo_sigmar90  = trunc(sigmar_accum / h.r100 * INT16SCALE);
		else if (p == size - 1 ) h.subhalo_sigmar100 = trunc(sigmar_accum / h.r100 * INT16SCALE);
	}
	
#ifdef SPHERICAL_OVERDENSITY
	h.SO_subhalo_central_particle = L2.particles[0];
	h.SO_subhalo_central_density  = L2.density[0]; 
	h.SO_subhalo_radius           = sqrt(L2.d2buffer[-1]); 
#endif 

    // We search for the max of vcirc, which is proportional to sqrt(G*M/R).
    // The 4th power of that is proportional to N^2/R^2.
    vmax = 0.0;
    for (int p=size/10; p<size; p++) {
	float v4 = p*p/L2.d2buffer[p];
	if (v4>vmax) { vmax = v4; rvmax = L2.d2buffer[p]; }
    }
    h.subhalo_rvcirc_max = sqrt(rvmax);    // Get to radial units
    h.subhalo_vcirc_max = sqrt(GMpart*sqrt(vmax))/ReadState.VelZSpace_to_kms;  // This is sqrt(N/R).
	
    x += offset; 
    x = WrapPosition(x);
    subhalo_x += offset; subhalo_x = WrapPosition(subhalo_x);
    assign_to_vector(h.x, x);
    assign_to_vector(h.subhalo_x, subhalo_x);
    return h;
};
