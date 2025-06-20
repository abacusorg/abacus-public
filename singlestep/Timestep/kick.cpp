// Copyright 2012-2025 The Abacus Developers
// SPDX-License-Identifier: GPL-3.0-or-later

/* kick.cpp
 * Kick all of the particles in a slab by a specified kick factor.
 *
 * Also contains routine to rescale the accelerations.  Of course, 
 * in principle, we could just attach that to the kick factor!
 *
 * Because we have the velocities and accelerations here, we also
 * compute some statistics for each cell suitable for monitoring
 * timestep adequacy.  Note that these are appropriate to the read-state
 * time, not the write-state time!
 *
 * We want to compute the rms of the velocity.  In principle, we'd like to 
 * subtract off the bulk velocity, but perhaps that's not so important?
 * Usually the mean velocity is a fair bit smaller than the rms.
 *
 * We'd like to compute the maximum component of the velocity.  This is for
 * scaling the velocity output.
 * 
 * We'd like to compute the maximum norm of the acceleration.  But perhaps 
 * it is enough to track the max component of the acceleration?
 *
 * If we can phrase these tests in a way that sums over components, then 
 * it would eventually be AVX-able.  That's considerably faster.
*/


inline void KickCell(Cell &c, FLOAT kick1, FLOAT kick2, int set_aux_dens) {
    FLOAT maxvel = 0.0;
    FLOAT maxacc = 0.0;
    FLOAT sumvel2 = 0.0;

    uint32_t N = c.count();
    #pragma omp simd reduction(max:maxvel,maxacc) reduction(+:sumvel2)
    for (uint32_t i=0;i<N;i++) {
        // First half kick, to get synchronous
        c.vel[i] += static_cast<FLOAT3>(c.acc[i]) * kick1;
        // Some simple stats
	sumvel2 += c.vel[i].norm2();
	maxvel = std::max(maxvel, std::abs(c.vel[i].x));
	maxacc = std::max(maxacc, std::abs(c.acc[i].x));
	maxvel = std::max(maxvel, std::abs(c.vel[i].y));
	maxacc = std::max(maxacc, std::abs(c.acc[i].y));
	maxvel = std::max(maxvel, std::abs(c.vel[i].z));
	maxacc = std::max(maxacc, std::abs(c.acc[i].z));
        // Second half kick, to advance to time i+1/2
	c.vel[i] += static_cast<FLOAT3>(c.acc[i]) * kick2;

        if(set_aux_dens) c.aux[i].set_density(c.acc[i].w);
    }
    if (c.count()>0) sumvel2/=c.count();  // Now this has the mean square velocity 
    c.ci->mean_square_velocity = sumvel2;
    c.ci->max_component_acceleration = maxacc;
    c.ci->max_component_velocity = maxvel;
}


void KickSlab(int slab, FLOAT kick1, FLOAT kick2, int set_aux_dens,
void (*KickCell)(Cell &c, FLOAT kick1, FLOAT kick2, int set_aux_dens)) {
    int cpd = CP->cpd;
    //#pragma omp parallel for schedule(static)
    //for (int y=0;y<cpd;y++) {
    NUMA_FOR(y,0,cpd, NO_CLAUSE, FALLBACK_DYNAMIC){
        for (int z = node_z_start; z < node_z_start + node_z_size; z++) {
            Cell c = CP->GetCell(slab, y, z);
            (*KickCell)(c,kick1,kick2,set_aux_dens);
        }
    }
    NUMA_FOR_END;
}

void RescaleAndCoAddAcceleration(int slab) {
    // The accelerations are computed with unit particle mass.
    // We need to rescale them to the correct cosmology.
    FLOAT rescale = -3.0*(P.Omega_M-P.Omega_Smooth)/(8.0*M_PI*P.np);
    accstruct * __restrict__ nacc = (accstruct *) SB->GetSlabPtr(AccSlab,slab);
    acc3struct * __restrict__ facc = (acc3struct *) SB->GetSlabPtr(FarAccSlab,slab);
    
    // Reverse the sign of the acceleration if we are making glass
    if(P.MakeGlass)
        rescale *= -1;

    // no ghost kick
    uint64 N = SS->size(slab);

    #ifdef DIRECTSINGLESPLINE
    FLOAT inv_eps3 = 1./(NFD->SofteningLengthInternal*NFD->SofteningLengthInternal*NFD->SofteningLengthInternal);
    #endif

    #pragma omp parallel
    {
        // the nontemporal clause doesn't seem to help, although icc -qstreaming-stores does
        //#pragma omp for simd schedule(static) nontemporal(nacc,facc)
        #pragma omp for simd schedule(static)
        for (uint64 j=0; j<N;j++) {
            #ifdef DIRECTSINGLESPLINE
            nacc[j] = (nacc[j]*inv_eps3+facc[j])*rescale;
            #else
            nacc[j] = (nacc[j]+facc[j] )*rescale;
            #endif

            #ifdef COMPUTE_FOF_DENSITY
            nacc[j].w -= WriteState.DensityKernelRad2;
            #endif
        }
    }
}


void ZeroAcceleration(int slab,int Slabtype) {
    // Null out the acceleration
    // Note that this is specific to accstruct, not acc3struct!
    accstruct *acc = (accstruct *) SB->GetSlabPtr(Slabtype,slab);
    memset(acc,0,SS->size(slab)*sizeof(accstruct));
}
