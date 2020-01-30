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

double3 ZelPos(integer3 ijk) {
    // This will be the position, in code units, of the initial grid.
    double3 p;
    //ijk = CP->WrapCell(ijk);
    p.x = (double)(ijk.x)/WriteState.ppd -.5;
    p.y = (double)(ijk.y)/WriteState.ppd -.5;
    p.z = (double)(ijk.z)/WriteState.ppd -.5;
    return p;
}

// If we're doing Zel'dovich only, then no problems: the original
// displacement times f is the velocity.  The current IC code does that.

// If we're doing 2LPT, then we will overwrite those initial velocities.

void KickCell_2LPT_1(Cell &c, FLOAT kick1, FLOAT kick2) {
    // Just store the acceleration
    int N = c.count();
    for (int i = 0; i < N; i++) {
        c.vel[i] = TOFLOAT3(c.acc[i]);
    }
}

void DriftCell_2LPT_1(Cell &c, FLOAT driftfactor) {
    assertf(P.is_np_perfect_cube(), "LPT reconstruction requires np (%d) to be a perfect cube.\n",P.np);
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
        c.pos[b] = 2.0*(ZelPos(c.aux[b].xyz())-cellcenter) - c.pos[b];
        c.pos[b] -= c.pos[b].round();
    }
}

void KickCell_2LPT_2(Cell &c, FLOAT kick1, FLOAT kick2) {
    // Now we can co-add the first two kicks to isolate the second-order
    // part.  This will be stored in the velocity.
    for (int i=0;i<c.count();i++) {
        c.vel[i] += TOFLOAT3(c.acc[i]);
    }
}

// We have an arena with the full IC 
// Extract the velocity field and fix the units
void unpack_ic_vel_slab(int slab){
    double convert_pos, convert_vel;
    get_IC_unit_conversions(convert_pos, convert_vel);
    unique_ptr<ICFile> ic_file = ICFile::FromFormat(P.ICFormat, slab);

    // Initialize the arena
    velstruct* velslab = (velstruct*) SB->AllocateArena(VelLPTSlab, slab);
    
    uint64 count = ic_file->unpack_to_velslab(velslab, convert_pos, convert_vel);
    
    assertf(count * sizeof(velstruct) == SB->SlabSizeBytes(VelLPTSlab, slab),
        "The size of particle vel data (%d) read from slab %d did not match the arena size (%d)\n",
        count*sizeof(velstruct), slab, SB->SlabSizeBytes(VelLPTSlab, slab));
}

// Returns the IC velocity of the particle with the given PID
// This "random access" lookup can be slow, probably not much we can do though
inline velstruct get_ic_vel(integer3 ijk, velstruct **velslabs, uint64 *velslab_sizes){
    // TODO: can probably simplify math
    int slab = ijk.x*P.cpd / WriteState.ippd;  // slab number
    uint64 slab_offset = ijk.x - ceil(((double)slab)*WriteState.ippd/P.cpd);  // number of planes into slab
    uint64 offset = ijk.z + WriteState.ippd*(ijk.y + WriteState.ippd*slab_offset);  // number of particles into slab
    
    if(velslabs[slab] == NULL){
        assertf(SB->IsSlabPresent(VelLPTSlab, slab),
            "IC vel slab %d not loaded.  Possibly an IC particle crossed more than FINISH_WAIT_RADIUS=%d slab boundaries?\n",
            slab, FINISH_WAIT_RADIUS);
        velslabs[slab] = (velstruct*) SB->GetSlabPtr(VelLPTSlab, slab);
        velslab_sizes[slab] = SB->SlabSizeBytes(VelLPTSlab, slab)/sizeof(velstruct);
    }

    velstruct* velslab = velslabs[slab];
    
    // Ensure that we're not reading past the end of the slab
    assertf(offset < velslab_sizes[slab],
        "Tried to read particle %d from IC slab %d, which only has %d particles\n",
        offset, slab, velslab_sizes[slab]);
    velstruct vel = velslab[offset];
    
    assertf(vel.is_finite(), "vel bad value: (%f,%f,%f)\n", vel.x, vel.y, vel.z);

    return vel;
}

void DriftCell_2LPT_2(Cell &c, FLOAT driftfactor) {
    // Now we have to adjust the positions and velocities
    assertf(P.is_np_perfect_cube(), "LPT reconstruction requires np (%d) to be a perfect cube.\n",P.np);
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

    // We don't want to be fetching the arena pointers for every particle, so cache them here
    // TODO: do we really want to do this for every cell?
    velstruct **velslabs = new velstruct*[P.cpd];
    uint64 *velslab_sizes = new uint64[P.cpd];
    for(int i = 0; i < P.cpd; i++){
        velslabs[i] = NULL;
        velslab_sizes[i] = 0;
    }
    
    double H = 1;
    double fsmooth = P.Omega_Smooth/P.Omega_M;
    double vel_to_displacement = 2.0/(19.0-18.0*fsmooth+sqrt(25.0-24.0*fsmooth)/H/H/P.Omega.M;
        // The second order displacement is vel/7H^2Omega_m
        // With a smooth component, the Omega_m remains the same, but the factor of 7
        // becomes 1+18 f_cluster - sqrt(1+24*f_cluster)

    for (int b = 0; b<e; b++) {
        // The first order displacement
        displ1 = ZelPos(c.aux[b].xyz())-cellcenter-c.pos[b];
        // displ2 = c.vel[b]/7/H/H/P.Omega_M;
        displ2 = c.vel[b]*vel_to_displacement;
        
        // Internally, everything is in a unit box, so we actually don't need to normalize by BoxSize before rounding
        displ1 -= displ1.round();
        displ2 -= displ2.round();

        c.pos[b] = ZelPos(c.aux[b].xyz())-cellcenter + displ1+displ2;
        c.pos[b] -= c.pos[b].round();
    
        // If we were supplied with Zel'dovich velocities and displacements,
        // we want to re-read the IC files to restore the velocities, which were overwritten in the 1st 2LPT step
        velstruct vel1;
        if(WriteState.Do2LPTVelocityRereading){
            vel1 = get_ic_vel(c.aux[b].xyz(), velslabs, velslab_sizes);
        }
        // If we were only supplied with Zel'dovich displacements, then construct the linear theory velocity
        // Or, if we're doing 3LPT, then we also want the ZA velocity for the next step
        else if(strcmp(P.ICFormat, "Zeldovich") == 0 || P.LagrangianPTOrder > 2){
            vel1 = WriteState.f_growth*displ1*convert_velocity;
        }
        // Unexpected IC format; fail.
        else {
            QUIT("Unexpected ICformat \"%s\" in 2LPT code.  Must be one of: Zeldovich, RVdoubleZel, RVZel\n", P.ICFormat);
        }
        
        // Here's where we add on the second-order velocity.
        // This factor of 2*f_growth is correct even for Smooth components
        // (where f_growth < 1) because the time dependence is still the square of first order.
        c.vel[b] = vel1 + WriteState.f_growth*2*displ2*convert_velocity;
    }

    delete[] velslabs;
    delete[] velslab_sizes;
}

// 3LPT kick and drift

void KickCell_2LPT_3(Cell &c, FLOAT kick1, FLOAT kick2) {
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
        assertf(displ3.norm() < displ12.norm(), "Error: 3rd-order LPT displacement (%f, %f, %f) is larger than 1st+2nd order (%f, %f, %f)!\n",
               displ3.x, displ3.y, displ3.z, displ12.x, displ12.y, displ12.z);

        c.pos[b] = ZelPos(c.aux[b].xyz())-cellcenter + displ12 + displ3;
        c.pos[b] -= c.pos[b].round();
        
        velstruct vel3 = 3*WriteState.f_growth*displ3*convert_velocity;
        
        // If we were supplied with Zel'dovich velocities and displacements,
        // we want to re-read the 1st order velocity
        if(WriteState.Do2LPTVelocityRereading){
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
