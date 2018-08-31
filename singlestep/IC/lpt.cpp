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

// Bookkeeping for 2LPT velocity re-reading
uint64 *vel_ics_npart;

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
    p.x = (double)(ijk.x)/WriteState.ppd -.5;
    p.y = (double)(ijk.y)/WriteState.ppd -.5;
    p.z = (double)(ijk.z)/WriteState.ppd -.5;
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
    double3 cellcenter = PP->WrapCellCenter(c.ijk);
#endif

    for (int b = 0; b<e; b++) {
        // Flip the displacement: q = pos-grid; new = grid-q = grid*2-pos
        // Need to be careful in cell-centered case.  Now pos is cell-centered,
        // and grid is global.  We should subtract the cell-center from ZelPos.
        c.pos[b] = 2.0*(ZelPos(c.aux[b].pid())-cellcenter) - c.pos[b]; // Does abacus automatically box wrap this?
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

// Allocate 2LPT velocity bookkeeping
void init_2lpt_rereading(){
    vel_ics_npart = new uint64[P.cpd]();
}

// Free 2LPT velocity bookkeeping
void finish_2lpt_rereading(){
    for(int i = 0; i < P.cpd; i++)
        assertf(!LBW->IDPresent(VelLPTSlab, i), "A 2LPT velocity slab was left loaded\n");
    delete[] vel_ics_npart;
}

// Load an IC velocity slab into an arena through the LoadIC module
void load_ic_vel_slab(int slabnum){
    slabnum = PP->WrapSlab(slabnum);
    assertf(!LBW->IDPresent(VelLPTSlab, slabnum), "Trying to re-load velocity IC slab %d, which is already loaded.\n", slabnum);
    
    assertf(vel_ics_npart[slabnum] == 0, "Trying to load velocity IC slab %d, which should have already been completely re-read for 2LPT.\n", slabnum);
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
    vel_ics_npart[slabnum] = LBW->IDSizeBytes(VelLPTSlab, slabnum) / sizeof(velstruct);
    
    assertf(count == vel_ics_npart[slabnum], "The number of particles (%d) read from slab %d did not match the number computed from its file size (%d)\n", slabnum, count, vel_ics_npart[slabnum]);
    delete ic_file;
}

// Returns the IC velocity of the particle with the given PID
inline velstruct* get_ic_vel(uint64 pid){
    integer3 ijk = ZelIJK(pid);
    int slabnum = ijk.x*P.cpd / WriteState.ppd;  // slab number

    uint64 slab_offset = ijk.x - ceil(((double)slabnum)*WriteState.ppd/P.cpd);  // number of planes into slab
    // We know the exact slab number and position of the velocity we want.
    uint64 offset = ijk.z + WriteState.ppd*(ijk.y + WriteState.ppd*slab_offset);
    
    assertf(LBW->IDPresent(VelLPTSlab, slabnum), "IC vel slab %d not loaded.  Possibly an IC particle crossed two slab boundaries?\n", slabnum);
    
    velstruct* slab = (velstruct*) LBW->ReturnIDPtr(VelLPTSlab, slabnum);
    
    assertf(offset < vel_ics_npart[slabnum], "Tried to read particle %lld from IC slab %d, which only has %d particles\n", offset, slabnum, vel_ics_npart[slabnum]);  // Ensure that we're not reading past the end of the slab
    velstruct* vel = slab + offset;
    
    assertf(std::isfinite(vel->x), "vel.x bad value: %f\n", vel->x);
    assertf(std::isfinite(vel->y), "vel.y bad value: %f\n", vel->y);
    assertf(std::isfinite(vel->z), "vel.z bad value: %f\n", vel->z);

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
    double3 cellcenter = PP->WrapCellCenter(c.ijk);
#endif
    
    double H = 1;
    for (int b = 0; b<e; b++) {
        // The first order displacement
        displ1 = ZelPos(c.aux[b].pid())-cellcenter-c.pos[b];
        // The second order displacement is vel/7H^2Omega_m
        displ2 = c.vel[b]/7/H/H/P.Omega_M;
        
        // Internally, everything is in a unit box, so we actually don't need to normalize by BoxSize before rounding
        displ1 -= displ1.round();
        displ2 -= displ2.round();

        c.pos[b] = ZelPos(c.aux[b].pid())-cellcenter + displ1+displ2;
        c.pos[b] -= c.pos[b].round();
    
        // If we were only supplied with Zel'dovich displacements, then construct the linear theory velocity
        // Or, if we're doing 3LPT, then we also want the ZA velocity for the next step
        velstruct vel1;
        if(strcmp(P.ICFormat, "Zeldovich") == 0 || P.LagrangianPTOrder > 2){
            vel1 = WriteState.f_growth*displ1*convert_velocity;
        }
        // If we were supplied with Zel'dovich velocities and displacements,
        // we want to re-read the IC files to restore the velocities, which were overwritten above
        else if(WriteState.Do2LPTVelocityRereading){
            vel1 = *get_ic_vel(c.aux[b].pid());
        }
        // Unexpected IC format; fail.
        else {
            assertf(false, "Unexpected ICformat \"%s\" in 2LPT code.  Must be one of: Zeldovich, RVdoubleZel, RVZel\n", P.ICFormat);
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
    double3 cellcenter = PP->WrapCellCenter(c.ijk);
#endif
    
    double H = 1;
    for (int b = 0; b<e; b++) {
        // The first+second order displacement
        displ12 = c.pos[b] - (ZelPos(c.aux[b].pid())-cellcenter);
        displ12 -= displ12.round();
            
        // Third order displacement
        displ3 = (2./(3*H*H*P.Omega_M)*TOFLOAT3(c.acc[b]) - (7./(3*H*WriteState.f_growth*convert_velocity)*c.vel[b] - 4./3*displ12))/6;
        displ3 -= displ3.round();
        assertf(displ3.norm() < displ12.norm(), "Error: 3rd-order LPT displacement (%f, %f, %f) is larger than 1st+2nd order (%f, %f, %f)!\n",
               displ3.x, displ3.y, displ3.z, displ12.x, displ12.y, displ12.z);

        c.pos[b] = ZelPos(c.aux[b].pid())-cellcenter + displ12 + displ3;
        c.pos[b] -= c.pos[b].round();
        
        velstruct vel3 = 3*WriteState.f_growth*displ3*convert_velocity;
        
        // If we were supplied with Zel'dovich velocities and displacements,
        // we want to re-read the 1st order velocity
        if(WriteState.Do2LPTVelocityRereading){
            velstruct vel1, vel1_ic, vel2;
            vel1_ic = *get_ic_vel(c.aux[b].pid());
            vel2 = c.vel[b] - displ12*WriteState.f_growth*convert_velocity;  // Isolate the 2nd order vel from the linear 1st order
            c.vel[b] = vel1_ic + vel2;
            
        }
        // If we were only supplied with Zel'dovich displacements,
        // then nothing special to be done here
        c.vel[b] += vel3;
    }
}
