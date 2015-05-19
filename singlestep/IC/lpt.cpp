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
    uint64 ppd = WriteState.ppd;		// Must ensure promotion to 64-bit
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

double3 ZelPos(uint64 PID) {
    // Given a PID, unwrap it base-CPD and return the position.
    integer3 ijk;
    uint64 residual = PID/WriteState.ppd;
    ijk.z = PID-residual*WriteState.ppd;
    ijk.x = residual/WriteState.ppd;
    ijk.y = residual-ijk.x*WriteState.ppd;
    return ZelPos(ijk);
}


// If we're doing Zel'dovich only, then no problems: the original
// displacement times f is the velocity.  The current IC code does that.

// If we're doing 2LPT, then we will overwrite those initial velocities.

void KickCell_2LPT_1(Cell c, accstruct *cellacc, FLOAT kick1, FLOAT kick2) {
    // Just store the acceleration
    for (int i=0;i<c.count();i++) {
        c.vel[i] = cellacc[i];
    }
}

void DriftCell_2LPT_1(Cell c, FLOAT driftfactor) {
    assertf(P.is_np_perfect_cube(), "LPT reconstruction requires np to be a perfect cube.\n");
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
        c.pos[b] = 2.0*(ZelPos(c.aux[b].pid())-cellcenter) - c.pos[b];
    }
}

void KickCell_2LPT_2(Cell c, accstruct *cellacc, FLOAT kick1, FLOAT kick2) {
    // Now we can co-add the first two kicks to isolate the second-order
    // part.  This will be stored in the velocity.
    for (int i=0;i<c.count();i++) {
        c.vel[i] += cellacc[i];
    }
}

void DriftCell_2LPT_2(Cell c, FLOAT driftfactor) {
    // Now we have to adjust the positions and velocities
    // The following probably assumes Omega_m = 1
    assertf(P.is_np_perfect_cube(), "LPT reconstruction requires np to be a perfect cube.\n");
    int e = c.count();
    double3 displ1, displ2;
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
    
    for (int b = 0; b<e; b++) {
	// The first order displacement
	displ1 = ZelPos(c.aux[b].pid())-cellcenter-c.pos[b];
	// The second order displacement is vel/7H^2
	// TODO: Is this line ok if we don't have cosm.H0=1?
	displ2 = c.vel[b]/7/WriteState.HubbleNow/WriteState.HubbleNow;
	//
        c.pos[b] = ZelPos(c.aux[b].pid())-cellcenter+displ1+displ2;
	c.vel[b] = WriteState.f_growth*(displ1 + 2*displ2)*convert_velocity;
    }
}

// Going to 3LPT requires one more force calculation, which then gets
// used for both velocity and positions.  This requires that the acceleration
// be preserved until the Drift step and then used to adjust both.


