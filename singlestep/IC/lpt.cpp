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


// Our Zel'dovich IC implementation must be coordinated with choices here.
// Therefore, we will codify some items in some simple functions:

uint64 ZelPID(integer3 ijk) {
    // This will set the PID to be stored in aux.
    uint64 ppd = WriteState.ppd;                // Must ensure promotion to 64-bit
    return (ppd*ijk.x+ijk.y)*ppd+ijk.z;
}

double3 ZelPos(integer3 ijk) {
    // This will be the position, in code units, of the initial grid.
    double3 p;
    //ijk = PP->WrapCell(ijk);
    p.x = (float)(ijk.x)/WriteState.ppd -.5;
    p.y = (float)(ijk.y)/WriteState.ppd -.5;
    p.z = (float)(ijk.z)/WriteState.ppd -.5;
    return p;
}

integer3 ZelIJK(uint64 PID) {
    // Given a PID, unwrap it base-PPD and return the grid indices
    integer3 ijk;
    uint64 residual = PID/WriteState.ppd;
    ijk.z = PID-residual*WriteState.ppd;
    ijk.x = residual/WriteState.ppd;
    ijk.y = residual-ijk.x*WriteState.ppd;
    return ijk;
}

double3 ZelPos(uint64 PID) {
    // Given a PID, unwrap it base-PPD and return the position.
    return ZelPos(ZelIJK(PID));
}


// If we're doing Zel'dovich only, then no problems: the original
// displacement times f is the velocity.  The current IC code does that.

// If we're doing 2LPT, then we will overwrite those initial velocities.

void KickCell_2LPT_1(Cell c, accstruct *cellacc, FLOAT kick1, FLOAT kick2) {
    // Just store the acceleration
    int N = c.count();
    for (int i = 0; i < N; i++) {
        c.vel[i] = cellacc[i];
    }
}

void DriftCell_2LPT_1(Cell c, FLOAT driftfactor) {
    assertf(P.is_np_perfect_cube(), "LPT reconstruction requires np (%d) to be a perfect cube.\n",P.np);
    int e = c.count();
#ifdef GLOBALPOS
    // Set cellcenter to zero to return to box-centered positions
    double3 cellcenter = double3(0.0);
#else
    double3 cellcenter = PP->WrapCellCenter(c.ijk);
#endif

    for (int b = 0; b<e; b++) {
        // Flip the displacement: q = pos-grid; new = grid-q = grid*2-pos
        // Need to be careful in cell-centered case.  Now pos is cell-centered,
        // and grid is global.  We should subtract the cell-center from ZelPos.
        c.pos[b] = 2.0*(ZelPos(c.aux[b].pid())-cellcenter) - c.pos[b]; // Does abacus automatically box wrap this?
        for (int i = 0; i < 3; i++){
            c.pos[b][i] -= round(c.pos[b][i]);
        }
    }
}

void KickCell_2LPT_2(Cell c, accstruct *cellacc, FLOAT kick1, FLOAT kick2) {
    // Now we can co-add the first two kicks to isolate the second-order
    // part.  This will be stored in the velocity.
    for (int i=0;i<c.count();i++) {
        c.vel[i] += cellacc[i];
    }
}

// Load an IC velocity slab into an arena through the LoadIC module
void load_ic_vel_slab(int slabnum){
    slabnum = PP->WrapSlab(slabnum);
    assertf(!LBW->IDPresent(VelLPTSlab, slabnum), "Trying to re-load velocity IC slab %d, which is already loaded.\n", slabnum);
    
    assertf(vel_ics[slabnum].n_part == 0, "Trying to load velocity IC slab %d, which should have already been completely re-read for 2LPT.\n", slabnum);
    STDLOG(1, "Re-reading velocity IC slab %d\n", slabnum);
    
    // Initialize the arena
    velstruct* slab = (velstruct*) LBW->AllocateArena(VelLPTSlab, slabnum);  // Automatically determines the arena size from the IC file size 
    
    // Use the LoadIC module to do the reading
    ICfile* ic_file;
    if(strcmp(P.ICFormat, "RVdoubleZel") == 0){
        ic_file = new ICfile_RVdoubleZel((char*)LBW->ReadSlabDescriptorName(VelLPTSlab, slabnum).c_str());
    } else if(strcmp(P.ICFormat, "RVZel") == 0){
        ic_file = new ICfile_RVZel((char*)LBW->ReadSlabDescriptorName(VelLPTSlab, slabnum).c_str());
    }
    
    uint64 count = 0;
    double3 pos;
    velstruct vel;
    auxstruct aux;
    while (ic_file->getparticle(&pos, &vel, &aux)) {
        vel *= WriteState.VelZSpace_to_Canonical;
        slab[count] = vel;
        count++;
    }
    
    // Initialize the VelIC struct, as a signal that the slab is truly ready
    vel_ics[slabnum].n_part = LBW->IDSizeBytes(VelLPTSlab, slabnum) / sizeof(velstruct);
    
    assertf(count == vel_ics[slabnum].n_part, "The number of particles (%d) read from slab %d did not match the number computed from its file size (%d)\n", slabnum, count, vel_ics[slabnum].n_part);
    delete ic_file;
}

// Loads the neighboring IC vel slabs if they are unloaded
void load_ic_vel_neighbors(int slabnum){
    for(int i = slabnum-1; i <= slabnum+1; i++){
        int wi = PP->WrapSlab(i);
        if (vel_ics[wi].n_part == 0){
            load_ic_vel_slab(wi);
        }
    }
}

// Unloaded neighboring IC vel slabs, if they are finished being used
// Note that Drift->done(slabnum) will not be true, even though it is for our purposes
void unload_finished_ic_vel_neighbors(int slabnum, Dependency* Drift){
    // Previous slab
    if(Drift->done(slabnum-1) && Drift->done(slabnum-2)){ // the other neighbor is the current slab
        STDLOG(1, "Finished re-reading all velocities from slab %d; unloading slab.\n", slabnum-1);
        LBW->DeAllocate(VelLPTSlab, slabnum-1);
    }
    // Current slab
    if(Drift->done(slabnum-1) && Drift->done(slabnum+1)){ // the middle is the current slab
        STDLOG(1, "Finished re-reading all velocities from slab %d; unloading slab.\n", slabnum);
        LBW->DeAllocate(VelLPTSlab, slabnum);
    }
    // Next slab
    if(Drift->done(slabnum+1) && Drift->done(slabnum+2)){ // the other neighbor is the current slab
        STDLOG(1, "Finished re-reading all velocities from slab %d; unloading slab.\n", slabnum+1);
        LBW->DeAllocate(VelLPTSlab, slabnum+1);
    }
}

// Returns the IC velocity of particle number "offset" relative to the beginning of the slab
inline velstruct* get_ic_vel(int slabnum, uint64 offset){
    assertf(LBW->IDPresent(VelLPTSlab, slabnum), "IC vel slab %d not loaded.  Possibly an IC particle crossed two slab boundaries?\n", slabnum);
    
    velstruct* slab = (velstruct*) LBW->ReturnIDPtr(VelLPTSlab, slabnum);
    
    assertf(offset < vel_ics[slabnum].n_part, "Tried to read particle %lld from IC slab %d, which only has %d particles\n", offset, slabnum, vel_ics[slabnum].n_part);  // Ensure that we're not reading past the end of the slab
    velstruct* vel = slab + offset;
    
    return vel;
}

void DriftCell_2LPT_2(Cell c, FLOAT driftfactor) {
    // Now we have to adjust the positions and velocities
    // The following probably assumes Omega_m = 1
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
    double3 cellcenter = PP->WrapCellCenter(c.ijk);
#endif
    
    double H = 1;
    for (int b = 0; b<e; b++) {
        // The first order displacement
        displ1 = ZelPos(c.aux[b].pid())-cellcenter-c.pos[b];
        // The second order displacement is vel/7H^2Omega_m
        displ2 = c.vel[b]/7/H/H/P.Omega_M;
        
        for (int i = 0; i < 3; i++){
            // Internally, everything is in a unit box, so we actually don't need to normalize by BoxSize before rounding
            displ1[i] -= round(displ1[i]);
            displ2[i] -= round(displ2[i]);
        }

        c.pos[b] = ZelPos(c.aux[b].pid())-cellcenter + displ1+displ2;
        for(int i = 0; i < 3; i++)
            c.pos[b][i] -= round(c.pos[b][i]);
    
        // If we were only supplied with Zel'dovich displacements, then construct the linear theory velocity:
        velstruct vel1;
        if(strcmp(P.ICFormat, "Zeldovich") == 0){
            vel1 = WriteState.f_growth*displ1*convert_velocity;
        }
        // If we were supplied with Zel'dovich velocities and displacements,
        // we want to re-read the IC files to restore the velocities, which were overwritten above
        else
        if(strcmp(P.ICFormat, "RVdoubleZel") == 0 || strcmp(P.ICFormat, "RVZel") == 0){
            integer3 ijk = ZelIJK(c.aux[b].pid());
            int slab = ijk.x*P.cpd / WriteState.ppd;  // slab number
            
            int slab_offset = ijk.x - ceil(((double)slab)*WriteState.ppd/P.cpd);  // number of planes into slab
            // We know the exact slab number and position of the velocity we want.
            uint64 offset = ijk.z + WriteState.ppd*(ijk.y + WriteState.ppd*slab_offset);
            velstruct* vel = get_ic_vel(slab, offset);
            
            vel1.x = vel->x;
            vel1.y = vel->y;
            vel1.z = vel->z;
            assertf(!isnan(vel1.x), "vel1.x is nan: %f\n", vel1.x);
            assertf(!isnan(vel1.y), "vel1.y is nan: %f\n", vel1.y);
            assertf(!isnan(vel1.z), "vel1.z is nan: %f\n", vel1.z);
        }
        // Unexpected IC format; fail.
        else {
            assertf(false, "Unexpected ICformat in 2LPT code.  Must be one of: Zeldovich, RVdoubleZel, RVZel\n");
        }
        
        c.vel[b] = vel1 + WriteState.f_growth*2*displ2*convert_velocity;
    }
}

// Going to 3LPT requires one more force calculation, which then gets
// used for both velocity and positions.  This requires that the acceleration
// be preserved until the Drift step and then used to adjust both.



