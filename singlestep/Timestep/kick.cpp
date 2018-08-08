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


    // TODO: Now that accstruct is not simply float3, we have do this
    // explicitly, so we lose SIMD.  Any tricks needed?
void KickCell(Cell &c, accstruct *cellacc, FLOAT kick1, FLOAT kick2) {
    FLOAT maxvel = 0.0;
    FLOAT maxacc = 0.0;
    FLOAT sumvel2 = 0.0;

    uint32_t N = c.count();
    for (uint32_t i=0;i<N;i++) {
        // First half kick, to get synchronous
        c.vel[i] += TOFLOAT3(c.acc[i]) * kick1;
        // Some simple stats
	sumvel2 += c.vel[i].norm2();
	maxvel = std::max(maxvel, fabs(c.vel[i].x));
	maxacc = std::max(maxacc, fabs(c.acc[i].x));
	maxvel = std::max(maxvel, fabs(c.vel[i].y));
	maxacc = std::max(maxacc, fabs(c.acc[i].y));
	maxvel = std::max(maxvel, fabs(c.vel[i].z));
	maxacc = std::max(maxacc, fabs(c.acc[i].z));
        // Second half kick, to advance to time i+1/2
	c.vel[i] += TOFLOAT3(c.acc[i]) * kick2;
    }
    if (c.count()>0) sumvel2/=c.count();  // Now this has the mean square velocity 
    c.ci->mean_square_velocity = sumvel2;
    c.ci->max_component_acceleration = maxacc;
    c.ci->max_component_velocity = maxvel;
}


void KickSlab(int slab, FLOAT kick1, FLOAT kick2,
void (*KickCell)(Cell &c, accstruct *cellacc, FLOAT kick1, FLOAT kick2)) {
    accstruct *acc = (accstruct *) LBW->ReturnIDPtr(AccSlab,slab);
    int cpd = PP->cpd;
    #pragma omp parallel for schedule(static)
    for (int y=0;y<cpd;y++) {
        for (int z=0;z<cpd;z++) {
            Cell c = PP->GetCell(slab, y, z);
            accstruct *cellacc = acc+c.ci->startindex;
            (*KickCell)(c,cellacc,kick1,kick2);
        }
    }
}

void RescaleAndCoAddAcceleration(int slab) {
    // The accelerations are computed with unit particle mass.
    // We need to rescale them to the correct cosmology.
    FLOAT rescale = -3.0*P.Omega_M/(8.0*M_PI*P.np);
    accstruct *nacc = (accstruct *) LBW->ReturnIDPtr(AccSlab,slab);
    acc3struct *facc = (acc3struct *) LBW->ReturnIDPtr(FarAccSlab,slab);
    
    // Reverse the sign of the acceleration if we are making glass
    if(strcmp(P.ICFormat, "Glass") == 0)
        rescale *= -1;
    
    uint64 N = Slab->size(slab);

    #ifdef DIRECTSINGLESPLINE
    FLOAT inv_eps3 = 1./(JJ->SofteningLengthInternal*JJ->SofteningLengthInternal*JJ->SofteningLengthInternal);
    #endif
    
    #pragma omp parallel for schedule(static)
    // TODO: Because nacc and facc can differ in type, we can't use SIMD.  
    //       Ok?  Perhaps bandwidth limited anyways?
    for (uint64 j=0; j<N;j++) {
        #ifdef DIRECTSINGLESPLINE
        nacc[j] = (nacc[j]*inv_eps3+facc[j])*rescale;
        #else
        nacc[j] = (nacc[j]+facc[j] )*rescale;
        #endif
    }
}



void ZeroAcceleration(int slab,int Slabtype) {
    // Null out the acceleration
    // Note that this is specific to accstruct, not acc3struct!
    accstruct *acc = (accstruct *) LBW->ReturnIDPtr(Slabtype,slab);
    memset(acc,0,Slab->size(slab)*sizeof(accstruct));
}
