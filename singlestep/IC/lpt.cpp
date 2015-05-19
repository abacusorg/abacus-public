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
    for (int i=0;i<c.count();i++) {
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

void DriftCell_2LPT_2(Cell c, FLOAT driftfactor) {
    // Now we have to adjust the positions and velocities
    // The following probably assumes Omega_m = 1
    assertf(P.is_np_perfect_cube(), "LPT reconstruction requires np (%d) to be a perfect cube.\n",P.np);
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
    
    //double H = WriteState.HubbleNow*P.H0;
    //double H = WriteState.HubbleNow;
    double H = 1;
    //printf("HubbleNow: %f\ncosm->C.H0: %f\nP.H0: %f\nH: %f\n", WriteState.HubbleNow, cosm->C.H0, P.H0, H);
    for (int b = 0; b<e; b++) {
        // The first order displacement
        displ1 = ZelPos(c.aux[b].pid())-cellcenter-c.pos[b];
        // The second order displacement is vel/7H^2Omega_m
        displ2 = c.vel[b]/7/H/H/P.Omega_M;
        
        for (int i = 0; i < 3; i++){
            // Internally, everything is in a unit box, so we actually don't need to normalize
            //displ1[i] -= P.BoxSize*round(displ1[i]/P.BoxSize);
            //displ2[i] -= P.BoxSize*round(displ2[i]/P.BoxSize);
            displ1[i] -= round(displ1[i]);
            displ2[i] -= round(displ2[i]);
        }

        double3 total_disp = displ1 + displ2;
        for (int i = 0; i < 3; i++)
            total_disp[i] -= round(total_disp[i]);
        double max_comp = total_disp.maxcomponent() > (-total_disp).maxcomponent() ? total_disp.maxcomponent() : -((-total_disp).maxcomponent());
        /*if(max_comp >= 1./33){
            printf("max_comp of total_disp: %g\n", max_comp);
        }*/

        c.pos[b] = ZelPos(c.aux[b].pid())-cellcenter + displ1+displ2;
        for(int i = 0; i < 3; i++)
        c.pos[b][i] -= round(c.pos[b][i]);
    
        // If we were only supplied with Zel'dovich displacements, then construct the linear theory velocity:
        double3 vel1;
        if(strcmp(P.ICFormat, "Zeldovich") == 0){
            vel1 = WriteState.f_growth*displ1;
        }
        // If we were supplied with Zel'dovich velocities and displacements,
        // we want to re-read the IC files to restore the velocities, which were overwritten above
        else
        #pragma omp critical
        if(strcmp(P.ICFormat, "RVdoubleZel") == 0){
            STDLOG(1,"Re-reading initial conditions files to restore 1st order velocities for 2LPT\n");
            integer3 ijk = ZelIJK(c.aux[b].pid());
            int slab = ijk.x*P.cpd / WriteState.ppd;  // slab number
            int slab_offset = ijk.x - ceil(((double)slab)*WriteState.ppd/P.cpd);  // number of planes into slab
            assertf(ceil((double)slab*WriteState.ppd/P.cpd) + slab_offset == ijk.x, "Wrong slab offset!\n");
            // We know the exact slab number and position of the velocity we want.
            long int offset = ijk.z + WriteState.ppd*(ijk.y + WriteState.ppd*slab_offset);
            offset *= sizeof(ICfile_RVdoubleZel::ICparticle);
            int seek_status = fseek(ic_fp[slab], offset, SEEK_SET);
            assertf(seek_status == 0, "Seek to position %ld failed in IC file %d\n", offset, slab);
            ICfile_RVdoubleZel::ICparticle p;
            int count = fread(&p, sizeof(ICfile_RVdoubleZel::ICparticle), 1, ic_fp[slab]);
            if(count != 1){
                perror("Error message");
            }
            assertf(count == 1, "Error reading IC file %d. Count = %d, expected 1\n", slab, count);
            vel1.x = p.vel[0];
            vel1.y = p.vel[1];
            vel1.z = p.vel[2];
            vel1 /= P.BoxSize;
        }
        // Unexpected IC format; fail.
        else {
            assertf(false, "Unexpected ICformat in 2LPT code.  Must be one of: Zeldovich, RVdoubleZel\n");
        }
        c.vel[b] = (vel1 + WriteState.f_growth*2*displ2)*convert_velocity;
    }
}

// Going to 3LPT requires one more force calculation, which then gets
// used for both velocity and positions.  This requires that the acceleration
// be preserved until the Drift step and then used to adjust both.



