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


void KickCell(Cell &c, accstruct *cellacc, FLOAT kick1, FLOAT kick2) {
    FLOAT maxacc = 0.0, maxvel = 0.0, sumvel2 = 0.0;
    FLOAT *vel = (FLOAT *)c.vel; 
    FLOAT *acc = (FLOAT *)cellacc;	// These are flattened arrays, not triples.

    uint32_t N = c.count()*3;
    #pragma simd reduction(max:maxvel) reduction(max:maxacc) reduction(+:sumvel2) assert
    for (uint32_t i=0;i<N;i++) {
        // First half kick, to get synchronous
        vel[i] += kick1 * acc[i];
        // Some simple stats
        FLOAT _vel = fabs(vel[i]);
        FLOAT _acc = fabs(acc[i]);
        if (_vel>maxvel) maxvel = _vel;
        if (_acc>maxacc) maxacc = _acc;
        sumvel2 += _vel*_vel;
        // Second half kick, to advance to time i+1/2
        vel[i] += kick2 * acc[i];
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
    accstruct *nacc = (accstruct *) LBW->ReturnIDPtr(NearAccSlab,slab);
    accstruct *facc = (accstruct *) LBW->ReturnIDPtr(AccSlab,slab);
    FLOAT *naccxyz = (FLOAT *) nacc;
    FLOAT *faccxyz = (FLOAT *) facc;
    
    // Reverse the sign of the acceleration if we are making glass
    if(strcmp(P.ICFormat, "Glass") == 0)
        rescale *= -1;
    
    uint64 N = Slab->size(slab)*3;
    
    #pragma omp parallel for schedule(static)
    #pragma simd assert
    for (uint64 j=0; j<N;j++) {
        faccxyz[j] = (faccxyz[j] + naccxyz[j])*rescale;
    }
}



void ZeroAcceleration(int slab,int Slabtype) {
    // Null out the acceleration
    accstruct *acc = (accstruct *) LBW->ReturnIDPtr(Slabtype,slab);
    memset(acc,0,Slab->size(slab)*sizeof(accstruct));
}
