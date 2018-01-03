#include "halostat.hh"
#include "sym3eigenval.cpp"

inline void WrapPosition(float3 x) {
    while (x.x> 0.5) x.x-=1.0;
    while (x.x<-0.5) x.x+=1.0;
    while (x.y> 0.5) x.y-=1.0;
    while (x.y<-0.5) x.y+=1.0;
    while (x.z> 0.5) x.z-=1.0;
    while (x.z<-0.5) x.z+=1.0;
    return;
}

#define assign_to_vector(a,b) { a[0] = b.x; a[1] = b.y; a[2] = b.z; }

HaloStat ComputeStats(int size, 
	posstruct *L1pos, velstruct *L1vel, auxstruct *L1aux, 
	FOFcell &L2, int ci, int cj, int ck) {
    // We are given the L1 particles as pos/vel/aux from [0,size).
    // These are in the original order; we don't care.
    // We are also given the L2 FOF results class.
    // Task is to fill a HaloStat object and return it.
    // The particle positions have already been wrapped to the first cell,
    // given as ci,cj,ck; we will wrap to global coords
    HaloStat h;

    h.N = size;
    // Compute the center of mass
    double3 com = double3(0.0, 0.0, 0.0);
    for (int p=0; p<size; p++) com += L1pos[p];
    double3 x = com/size;
    h.x[0] = com.x; h.x[1] = com.y; h.x[2] = com.z; 
    com = double3(0.0, 0.0, 0.0);
    for (int p=0; p<size; p++) com += L1vel[p];
    double3 v = com/size;
    assign_to_vector(h.v, v);

    // Find the largest L2 subhalos and the largest COM
    std::sort(L2.groups, L2.groups+L2.ngroups);
    	// Groups now in descending order of multiplicity
    for (int j=0; j<N_LARGEST_SUBHALOS; j++) 
	if (j<L2.ngroups) h.subhalo_N[j] = L2.groups[j].n; 
	    else h.subhalo_N[j] = 0;

    FOFparticle *start = L2.p + L2.groups[0].start;
    com = double3(0.0, 0.0, 0.0);
    for (int p=0; p<L2.groups[0].n; p++) 
    	com += L1pos[start[p].index()];
    double3 subhalo_x = com/h.subhalo_N[0];
    com = double3(0.0, 0.0, 0.0);
    for (int p=0; p<L2.groups[0].n; p++) 
    	com += L1vel[start[p].index()];
    double3 subhalo_v = com/h.subhalo_N[0];
    assign_to_vector(h.subhalo_v, subhalo_v);

    // Now we can go through the particles to compute radii and moments
    // We can use L2.d2buffer for scratch space; it is guaranteed to be big enough
    double vxx, vxy, vxz, vyy, vyz, vzz;
    float vmax, rvmax;
    vxx = vxy = vxz = vyy = vyz = vzz = 0.0;
    for (int p=0; p<size; p++) {
	posstruct dx = L1pos[p]-x;
	L2.d2buffer[p] = dx.norm2();
	velstruct dv = L1vel[p]-v;
	vxx += dv.x*dv.x; vxy += dv.x*dv.y; vxz += dv.x*dv.z;
	vyy += dv.y*dv.y; vyz += dv.y*dv.z; vzz += dv.z*dv.z;
    }
    FindEigenvalues(vxx, vxy, vxz, vyy, vyz, vzz, h.sigmav);
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
    h.vcirc_max = sqrt(sqrt(vmax));  // This is sqrt(N/R).
    	// TODO: Have to figure out what G*Mpart is in code units.

    // Repeat this, finding moments and radii around the largest subhalo COM
    vxx = vxy = vxz = vyy = vyz = vzz = 0.0;
    for (int p=0; p<size; p++) {
	posstruct dx = L1pos[p]-subhalo_x;
	L2.d2buffer[p] = dx.norm2();
	velstruct dv = L1vel[p]-subhalo_v;
	vxx += dv.x*dv.x; vxy += dv.x*dv.y; vxz += dv.x*dv.z;
	vyy += dv.y*dv.y; vyz += dv.y*dv.z; vzz += dv.z*dv.z;
    }
    FindEigenvalues(vxx, vxy, vxz, vyy, vyz, vzz, h.subhalo_sigmav);
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
    h.subhalo_vcirc_max = sqrt(sqrt(vmax));  // This is sqrt(N/R).
    	// TODO: Have to figure out what G*Mpart is in code units.

    // TODO: Need to supply this method.  Might convert test driver to use grid.
    posstruct offset = PP->CellCenter(ci,cj,ck);
    x += offset; WrapPosition(x);
    subhalo_x += offset; WrapPosition(subhalo_x);
    assign_to_vector(h.x, x);
    assign_to_vector(h.subhalo_x, subhalo_x);
    return h;
};
