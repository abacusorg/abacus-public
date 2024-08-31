/* Abacus forces are of sufficient accuracy that we can support an
alternative derivation of higher-order LPT, in which we derive the 
high-order displacements from force calculations. 

The parameter element LagrangianPTOrder should be =1 for Zel'dovich,
=2 for 2LPT, =3 for 3LPT.  This parameter doesn't change.

When we injest initial conditions, we receive Zel'dovich ICs and
build multipoles to form a legal state. We probably should call
that FullStepNumber==0.

An easy convention would be that if 
    LagrangianPTOrder>1 && FullStepNumber<=LagrangianPTOrder,
then we have more LPT work to do.  E.g., if we're going 2LPT, then
FullStepNumber 1 and 2 are not normal step.  This has implication for
group finding.  But it's reasonable because each LPT step really does
involve a full force calculation and rebinning.

Another implication of the LPT implementation is that we must have
access to the Initial Lagrangian Position.  This can be encoded
into the Aux variable as a particle ID number.  We might opt to
overwrite that PID with a group ID later in the code; that's ok.
But it means that the position of the Lagrangian integer index must
be consistent with loadIC.  Best to force that with a function!

*/

// Provide a simple routine to check that this is an LPT IC step.

int LPTStepNumber() {
    // This will be >0 if we're doing LPT IC work.
    // The step number for the LPT work will be returned.
    int step = WriteState.FullStepNumber;
    if (P.LagrangianPTOrder>1 && step<=P.LagrangianPTOrder)
        return step;
    else return 0;
}

inline double3 ZelPos(int64 ix, int64 iy, int64 iz) {
    // This will be the position, in code units, of the initial grid.
    double3 p;
    p.x = (double)ix/WriteState.ppd -.5;
    p.y = (double)iy/WriteState.ppd -.5;
    p.z = (double)iz/WriteState.ppd -.5;
    return p;
}

inline double3 ZelPos(integer3 ijk) {
    return ZelPos(ijk.x, ijk.y, ijk.z);
}

// If we're doing Zel'dovich only, then no problems: the original
// displacement times f is the velocity.  The current IC code does that.

// If we're doing 2LPT, then we will overwrite those initial velocities.

void KickCell_2LPT_1(Cell &c, FLOAT kick1, FLOAT kick2, int _set_aux_dens) {
    int N = c.count();
    FLOAT Canonical_to_VelZSpace = 1.0/WriteState.VelZSpace_to_Canonical;
    FLOAT invvscale = 1./ReadState.LPTVelScale;

#ifdef GLOBALPOS
    // Set cellcenter to zero to return to box-centered positions
    double3 cellcenter = double3(0.0);
#else
    double3 cellcenter = CP->WrapCellCenter(c.ijk);
#endif
    
    for (int i = 0; i < N; i++) {
        double3 displ = double3(c.pos[i]) - (ZelPos(c.aux[i].xyz()) - cellcenter);
        displ -= displ.round();
        velstruct linearvel = displ;
        velstruct v = c.vel[i]*Canonical_to_VelZSpace;

        // save the vel to the aux
        c.aux[i].set_velocity(v - linearvel, invvscale);

        // overwrite the vel with the acc
        c.vel[i] = TOFLOAT3(c.acc[i]);

        // why move vel to aux instead of writing acc in aux?
        // mostly because we were able to compute a strict upper
        // bound on (vel - linearvel) during the IC step, which
        // enables good use of the limited aux bits
    }
}

void DriftCell_2LPT_1(Cell &c, FLOAT driftfactor) {
    assertf(P.is_np_perfect_cube(), "LPT reconstruction requires np ({:d}) to be a perfect cube.\n",P.np);
    int e = c.count();
#ifdef GLOBALPOS
    // Set cellcenter to zero to return to box-centered positions
    double3 cellcenter = double3(0.0);
#else
    double3 cellcenter = CP->WrapCellCenter(c.ijk);
#endif

    for (int b = 0; b<e; b++) {
        // Flip the displacement: q = pos-grid; new = grid-q = grid*2-pos
        // Need to be careful in cell-centered case.  Now pos is cell-centered,
        // and grid is global.  We should subtract the cell-center from ZelPos.
        double3 pos_g = 2.0*(ZelPos(c.aux[b].xyz())-cellcenter) - c.pos[b];
        c.pos[b] = pos_g - pos_g.round();
    }
}

void KickCell_2LPT_2(Cell &c, FLOAT kick1, FLOAT kick2, int _set_aux_dens) {
    // Now we can co-add the first two kicks to isolate the second-order
    // part.  This will be stored in the velocity.
    for (int i=0;i<c.count();i++) {
        c.vel[i] += TOFLOAT3(c.acc[i]);
    }
}

void DriftCell_2LPT_2(Cell &c, FLOAT driftfactor) {
    // Now we have to adjust the positions and velocities
    assertf(P.is_np_perfect_cube(), "LPT reconstruction requires np ({:d}) to be a perfect cube.\n",P.np);
    int e = c.count();
    posstruct displ1, displ2;
    // This is the factor to convert from redshift-space displacements
    // to canonical velocities.
    // HubbleNow is H(z)/H_0.
    double convert_velocity = WriteState.VelZSpace_to_Canonical;
    // WriteState.ScaleFactor*WriteState.ScaleFactor *WriteState.HubbleNow;

#ifdef GLOBALPOS
    // Set cellcenter to zero to return to box-centered positions
    double3 cellcenter = double3(0.0);
#else
    double3 cellcenter = CP->WrapCellCenter(c.ijk);
#endif
    
    double H = 1;
    double fsmooth = P.Omega_Smooth/P.Omega_M;
    double vel_to_displacement = 2.0/(19.0-18.0*fsmooth-sqrt(25.0-24.0*fsmooth))/H/H/P.Omega_M;
        // The second order displacement is vel/7H^2Omega_m
        // With a smooth component, the Omega_m remains the same, but the factor of 7
        // becomes 1+18 f_cluster - sqrt(1+24*f_cluster)
    double vscale = ReadState.LPTVelScale;

    for (int b = 0; b<e; b++) {
        // The first order displacement
        // double3 until we ensure positions are wrapped
        double3 displ1_g = ZelPos(c.aux[b].xyz())-cellcenter-c.pos[b];
        double3 displ2_g = c.vel[b]*vel_to_displacement;
            // displ2 = c.vel[b]/7/H/H/P.Omega_M;
        
        // Internally, everything is in a unit box, so we actually don't need to normalize by BoxSize before rounding
        displ1 = displ1_g - displ1_g.round();
        displ2 = displ2_g - displ2_g.round();

        double3 pos_g = ZelPos(c.aux[b].xyz())-cellcenter + displ1+displ2;
        c.pos[b] = pos_g - pos_g.round();
    
        // Reconstruct the first-order velocity from the linear velocity
        // plus the correction stored in the aux
        velstruct vel1 = (c.aux[b].get_velocity(vscale) + displ1)*convert_velocity;

        // restore aux
        c.aux[b].zero_velocity();
        
        uint64 _sumA = 0, _sumB = 0;
        ICFile::set_taggable_bits(c.aux[b], _sumA, _sumB);
        /*
        // If we were only supplied with Zel'dovich displacements, then construct the linear theory velocity
        // Or, if we're doing 3LPT, then we also want the ZA velocity for the next step
        else if(P.ICFormat == "Zeldovich" || P.LagrangianPTOrder > 2){
            // TODO: where are we applying f_growth in the RVZel path?
            vel1 = WriteState.f_growth*displ1*convert_velocity;
        }*/
        
        // Here's where we add on the second-order velocity.
        // This factor of 2*f_growth is correct even for Smooth components
        // (where f_growth < 1) because the time dependence is still the square of first order.
        c.vel[b] = vel1 + WriteState.f_growth*2*displ2*convert_velocity;
    }
}

// 3LPT kick and drift

void KickCell_2LPT_3(Cell &c, FLOAT kick1, FLOAT kick2, int _set_aux_dens) {
    // We could update the 3LPT velocity from the acceleration here,
    // but then we'd be repeating the same calculations for the position update in the Drift
    return;
}

void DriftCell_2LPT_3(Cell &c, FLOAT driftfactor) {
    // Now we have to adjust the positions and velocities    
    // Usually, Drift doesn't need accelerations, but 3LPT does
    // TODO: Haven't considered Smooth components at 3rd order
    assertf(P.Omega_Smooth==0, "3rd order not implemented with Omega_Smooth!=0\n");
    int slab = c.ijk.x;
    
    int e = c.count();
    posstruct displ12, displ3;
    // This is the factor to convert from redshift-space displacements
    // to canonical velocities.
    // HubbleNow is H(z)/H_0.
    double convert_velocity = WriteState.VelZSpace_to_Canonical;
    // WriteState.ScaleFactor*WriteState.ScaleFactor *WriteState.HubbleNow;

#ifdef GLOBALPOS
    // Set cellcenter to zero to return to box-centered positions
    double3 cellcenter = double3(0.0);
#else
    double3 cellcenter = CP->WrapCellCenter(c.ijk);
#endif
    
    double H = 1;
    for (int b = 0; b<e; b++) {
        // The first+second order displacement
        displ12 = c.pos[b] - (ZelPos(c.aux[b].xyz())-cellcenter);
        displ12 -= displ12.round();
            
        // Third order displacement
        displ3 = (2./(3*H*H*P.Omega_M)*TOFLOAT3(c.acc[b]) - (7./(3*H*WriteState.f_growth*convert_velocity)*c.vel[b] - 4./3*displ12))/6;
        displ3 -= displ3.round();
        assertf(displ3.norm() < displ12.norm(), "Error: 3rd-order LPT displacement ({:f}, {:f}, {:f}) is larger than 1st+2nd order ({:f}, {:f}, {:f})!\n",
               displ3.x, displ3.y, displ3.z, displ12.x, displ12.y, displ12.z);

        c.pos[b] = ZelPos(c.aux[b].xyz())-cellcenter + displ12 + displ3;
        c.pos[b] -= c.pos[b].round();
        
        velstruct vel3 = 3*WriteState.f_growth*displ3*convert_velocity;
        
        // If we were supplied with Zel'dovich velocities and displacements,
        // we want to re-read the 1st order velocity
        if(1){
            velstruct vel1, vel1_ic, vel2;
            //vel1_ic = get_ic_vel(c.aux[b].xyz(), NULL, NULL);
            QUIT("3LPT needs update for new vel unpacking\n");
            vel2 = c.vel[b] - displ12*WriteState.f_growth*convert_velocity;  // Isolate the 2nd order vel from the linear 1st order
            c.vel[b] = vel1_ic + vel2;
            
        }
        // If we were only supplied with Zel'dovich displacements,
        // then nothing special to be done here
        c.vel[b] += vel3;
    }
}
