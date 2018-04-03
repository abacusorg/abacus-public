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

HaloStat ComputeStats(int size, 
	posstruct *L1pos, velstruct *L1vel, auxstruct *L1aux, 
	FOFcell &L2, posstruct offset) {
    // We are given the L1 particles as pos/vel/aux from [0,size).
    // These are in the original order; we don't care.
    // We are also given the L2 FOF results class.
    // Task is to fill a HaloStat object and return it.
    // The particle positions have already been wrapped to the first cell;
    // we will wrap to global coords, using offset
    // Velocities will be converted to our usual velocity output
    // convention of unit-box redshift-space displacements
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
    vxx = vxy = vxz = vyy = vyz = vzz = 0.0;
    for (int p=0; p<size; p++) {
	posstruct dx = L1pos[p]-x;
	L2.d2buffer[p] = dx.norm2();
	velstruct dv = L1vel[p]/ReadState.VelZSpace_to_Canonical - v;
	vxx += dv.x*dv.x; vxy += dv.x*dv.y; vxz += dv.x*dv.z;
	vyy += dv.y*dv.y; vyz += dv.y*dv.z; vzz += dv.z*dv.z;
    }
    FindEigenvalues(vxx, vxy, vxz, vyy, vyz, vzz, h.sigmav);
    for(int i = 0; i < 3; i++) h.sigmav[i] = sqrt(h.sigmav[i]);
    std::sort(L2.d2buffer, L2.d2buffer+size);
    h.r25 = sqrt(L2.d2buffer[size/4]);   h.r50 = sqrt(L2.d2buffer[size/2]);
    h.r75 = sqrt(L2.d2buffer[size*3/4]); h.r90 = sqrt(L2.d2buffer[size*9/10]);
    // We search for the max of vcirc, which is proportional to sqrt(G*M/R).
    // The 4th power of that is proportional to N^2/R^2.
    vmax = 0.0;
    for (int p=size/10; p<size; p++) {
	float v4 = p*p/L2.d2buffer[p];
	if (v4>vmax) { vmax = v4; rvmax = L2.d2buffer[p]; }
    }
    h.rvcirc_max = sqrt(rvmax);    // Get to radial units
    float GMpart = 3*P.Omega_M*pow(100*ReadState.BoxSizeHMpc,2)/(8*M_PI*P.np*ReadState.ScaleFactor);
    h.vcirc_max = sqrt(GMpart*sqrt(vmax))/ReadState.VelZSpace_to_kms;  // This is sqrt(G*M_particle*N/R).

    // Repeat this, finding moments and radii around the largest subhalo COM
    vxx = vxy = vxz = vyy = vyz = vzz = 0.0;
    for (int p=0; p<size; p++) {
	posstruct dx = L1pos[p]-subhalo_x;
	L2.d2buffer[p] = dx.norm2();
	velstruct dv = L1vel[p]/ReadState.VelZSpace_to_Canonical - subhalo_v;
	vxx += dv.x*dv.x; vxy += dv.x*dv.y; vxz += dv.x*dv.z;
	vyy += dv.y*dv.y; vyz += dv.y*dv.z; vzz += dv.z*dv.z;
    }
    FindEigenvalues(vxx, vxy, vxz, vyy, vyz, vzz, h.subhalo_sigmav);
    for(int i = 0; i < 3; i++) h.subhalo_sigmav[i] = sqrt(h.subhalo_sigmav[i]);
    std::sort(L2.d2buffer, L2.d2buffer+size);
    h.subhalo_r25 = sqrt(L2.d2buffer[size/4]);   
    h.subhalo_r50 = sqrt(L2.d2buffer[size/2]);
    h.subhalo_r75 = sqrt(L2.d2buffer[size*3/4]); 
    h.subhalo_r90 = sqrt(L2.d2buffer[size*9/10]);
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
