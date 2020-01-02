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
        c.pos[b] = 2.0*(ZelPos(c.aux[b].xyz())-cellcenter) - c.pos[b]; // Does abacus automatically box wrap this?
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

// Load an IC velocity slab into an arena through the LoadIC module
void load_ic_vel_slab(int slabnum){
    slabnum = CP->WrapSlab(slabnum);
    assertf(!SB->IsSlabPresent(VelLPTSlab, slabnum), "Trying to re-load velocity IC slab %d, which is already loaded.\n", slabnum);
    
    STDLOG(1, "Re-reading velocity IC slab %d\n", slabnum);
    
    // Initialize the arena
    velstruct* slab = (velstruct*) SB->AllocateArena(VelLPTSlab, slabnum);  // Automatically determines the arena size from the IC file size 
    
    // Use the LoadIC module to do the reading
    ICFile* ic_file;
    if(strcmp(P.ICFormat, "RVdoubleZel") == 0){
        ;//ic_file = new ICFile_RVdoubleZel((char*)SB->ReadSlabPath(VelLPTSlab, slabnum).c_str());
    } else if(strcmp(P.ICFormat, "RVZel") == 0){
        ;//ic_file = new ICFile_RVZel((char*)SB->ReadSlabPath(VelLPTSlab, slabnum).c_str());
    }
    
    uint64 count = 0;
    double3 pos;
    velstruct vel;
    auxstruct aux;
    /*while (ic_file->getparticle(&pos, &vel, &aux)) {
        vel *= WriteState.VelZSpace_to_Canonical;
        slab[count++] = vel;
    }*/
    
    assertf(count * sizeof(velstruct) == SB->SlabSizeBytes(VelLPTSlab, slabnum),
        "The size of particle vel data (%d) read from slab %d did not match the arena size (%d)\n", count, slab, SB->SlabSizeBytes(VelLPTSlab, slabnum));
    delete ic_file;
}

// Returns the IC velocity of the particle with the given PID
// TODO: this "random access" lookup is really slow.  Consider whether there's another way to phrase this.
inline velstruct* get_ic_vel(integer3 ijk){
    int slabnum = ijk.x*P.cpd / WriteState.ppd;  // slab number

    uint64 slab_offset = ijk.x - ceil(((double)slabnum)*WriteState.ppd/P.cpd);  // number of planes into slab
    // We know the exact slab number and position of the velocity we want.
    uint64 offset = ijk.z + WriteState.ppd*(ijk.y + WriteState.ppd*slab_offset);
    
    assertf(SB->IsSlabPresent(VelLPTSlab, slabnum), "IC vel slab %d not loaded.  Possibly an IC particle crossed more than FINISH_WAIT_RADIUS=%d slab boundaries?\n", slabnum, FINISH_WAIT_RADIUS);
    
    velstruct* slab = (velstruct*) SB->GetSlabPtr(VelLPTSlab, slabnum);
    
    // Ensure that we're not reading past the end of the slab
    assertf(offset < SB->SlabSizeBytes(VelLPTSlab, slabnum)/sizeof(velstruct),
        "Tried to read particle %d from IC slab %d, which only has %d particles\n", offset, slabnum, SB->SlabSizeBytes(VelLPTSlab, slabnum)/sizeof(velstruct));
    velstruct* vel = slab + offset;
    
    assertf(vel->is_finite(), "vel bad value: (%f,%f,%f)\n", vel->x, vel->y, vel->z);

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
    
    double H = 1;
    for (int b = 0; b<e; b++) {
        // The first order displacement
        displ1 = ZelPos(c.aux[b].xyz())-cellcenter-c.pos[b];
        // The second order displacement is vel/7H^2Omega_m
        displ2 = c.vel[b]/7/H/H/P.Omega_M;
        
        // Internally, everything is in a unit box, so we actually don't need to normalize by BoxSize before rounding
        displ1 -= displ1.round();
        displ2 -= displ2.round();

        c.pos[b] = ZelPos(c.aux[b].xyz())-cellcenter + displ1+displ2;
        c.pos[b] -= c.pos[b].round();
    
        // If we were only supplied with Zel'dovich displacements, then construct the linear theory velocity
        // Or, if we're doing 3LPT, then we also want the ZA velocity for the next step
        velstruct vel1;
        if(strcmp(P.ICFormat, "Zeldovich") == 0 || P.LagrangianPTOrder > 2){
            vel1 = WriteState.f_growth*displ1*convert_velocity;
        }
        // If we were supplied with Zel'dovich velocities and displacements,
        // we want to re-read the IC files to restore the velocities, which were overwritten in the 1st 2LPT step
        else if(WriteState.Do2LPTVelocityRereading){
            vel1 = *get_ic_vel(c.aux[b].xyz());
        }
        // Unexpected IC format; fail.
        else {
            QUIT("Unexpected ICformat \"%s\" in 2LPT code.  Must be one of: Zeldovich, RVdoubleZel, RVZel\n", P.ICFormat);
        }
        
        c.vel[b] = vel1 + WriteState.f_growth*2*displ2*convert_velocity;
    }
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
            vel1_ic = *get_ic_vel(c.aux[b].xyz());
            vel2 = c.vel[b] - displ12*WriteState.f_growth*convert_velocity;  // Isolate the 2nd order vel from the linear 1st order
            c.vel[b] = vel1_ic + vel2;
            
        }
        // If we were only supplied with Zel'dovich displacements,
        // then nothing special to be done here
        c.vel[b] += vel3;
    }
}
