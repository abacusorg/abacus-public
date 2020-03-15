Cosmology *InitializeCosmology(double ScaleFactor) {
    // Be warned that all of the Cosmology routines quote time units
    // by what is entered as H0.  The code wants to use H0=1 units.
    // But P.H0 is defined to be in km/s/Mpc!
    // If you want to compare any times or H(z), remember that H0=1.
    // This routine should only be called once.
    MyCosmology cosmo;
    cosmo.Omega_m = P.Omega_M;
    cosmo.Omega_smooth = P.Omega_Smooth;
    cosmo.Omega_K = P.Omega_K;
    cosmo.Omega_DE = P.Omega_DE;
    cosmo.H0 = 1.0;	
    cosmo.w0 = P.w0;
    cosmo.wa = P.wa;
    STDLOG(0,"Initialized Cosmology at a= %6.4f\n",ScaleFactor);
    return new Cosmology(ScaleFactor,cosmo);
    	// This will set cosm->current and cosm->next to epoch a=ScaleFactor.
}

void FillStateWithCosmology(State &S) {
    // Fill up the given state with information about the Cosmology cosm,
    // pulled from epoch cosm->next.
    S.ScaleFactor = cosm->next.a;
    S.Redshift = cosm->next.z;
    S.Time = cosm->next.t;                // In same units as H_0 was given
    S.etaK = cosm->next.etaK;
    S.etaD = cosm->next.etaD;
    S.Growth = cosm->next.growth;
    S.Growth_on_a_n = cosm->next.growth/pow(cosm->next.a, cosm->n);  // D/a^n
    S.f_growth = cosm->next.f_growth;
    S.w = cosm->next.w;
    S.HubbleNow = cosm->next.H;           // In same units as H_0 was given
    S.Htime = cosm->next.H*cosm->next.t;               // Time*H(z), in code units

    double total = cosm->next.OmegaHat_m+cosm->next.OmegaHat_X+cosm->next.OmegaHat_K;
    S.OmegaNow_m = cosm->next.OmegaHat_m/total;
    S.OmegaNow_K = cosm->next.OmegaHat_K/total;
    S.OmegaNow_DE = cosm->next.OmegaHat_X/total;
    S.fsmooth = cosm->C.Omega_smooth/cosm->C.Omega_m;

    S.HubbleTimeHGyr = 9.7782;
        // The Hubble Time is 9.7782 h^-1 Gyr
    S.HubbleTimeGyr = S.HubbleTimeHGyr*(100.0/P.H0);
    	// This is in Gyr
    S.BoxSizeMpc = P.BoxSize*(P.hMpc?(100.0/P.H0):1.0);
    S.BoxSizeHMpc = S.BoxSizeMpc*(P.H0/100.0);
    	// Redundant, but we might as well be explicit.
    S.ParticleMassMsun = 2.7746e11*(P.Omega_M-P.Omega_Smooth)*pow(P.H0/100,2.0)*pow(S.BoxSizeMpc,3.0)/P.np;
    	// This is in Msun.  
        // The critical density is 2.7746e11 h^2 Msun/Mpc^3
    S.ParticleMassHMsun = S.ParticleMassMsun*(P.H0/100);
    	// This is in h^-1 Msun.  

    S.CoordinateDistanceHMpc = (cosm->today.etaK-cosm->next.etaK)*2997.92;
        // This is the coordinate distance from z=0 to this epoch, in Mpc/h

    // The code uses canonical velocities, for the unit box and 1/H_0 time unit.
    // However, our output standard is redshift-space comoving displacements, again
    // in units where the box is one.
    // The conversion is v_canon = v_zspace * a^2 H(z)/H_0
    S.VelZSpace_to_Canonical = S.ScaleFactor * S.ScaleFactor * (S.HubbleNow/cosm->C.H0);

    // After output, one might want to convert to km/s.  
    // Here, we quote the number of km/s across the full box.
    // The proper size of the box is BoxSize/(1+z).
    // The Hubble parameter is (HubbleNow/cosm->C.H0)*H_0.
    // If hMpc is set, then we should use 100 km/s/Mpc instead of H_0.
    S.VelZSpace_to_kms = P.BoxSize*S.ScaleFactor*(S.HubbleNow/cosm->C.H0)*
    	(P.hMpc?100:P.H0);
}



double ChooseTimeStep(){
	// Choose the maximum allowable timestep
	// We start with the absolute maximum timestep allowed by the parameter file,
	// then see if it needs to be shorter.

	// cosm has already been loaded with the ReadState.ScaleFactor.
        // Don't advance time if we are still doing LPT
        STDLOG(0,"LPTStepNumber = %d, FullStepNumber = %d, PTO = %d\n",LPTStepNumber(), WriteState.FullStepNumber, P.LagrangianPTOrder);
        if(LPTStepNumber()>0) return 0.;

	double da = ReadState.ScaleFactor*P.TimeStepDlna;
	STDLOG(0,"da from Hubble Dlna limit is %f\n", da);
	if (da==0.0) return da;

	// TODO: I think below might be simplified if we tried to construct
	// cosm->BuildEpoch(cosm->current, cosm->next, cosm->current.a+da_max);
	// and then did interpolations with that.  Or after each attempt, 
	// call BuildEpoch and *test* whether cosm->next is acceptable,
	// then interpolate down.

	// Perhaps the next output is sooner than this?
	// TimeSlizez array might not be in order!  Look at all of them.
	for (int i = 0; i < P.nTimeSlice; i ++){
		double tsa = 1.0/(1+P.TimeSliceRedshifts[i]);
		if (ReadState.Redshift > P.TimeSliceRedshifts[i]+1e-12 && ReadState.ScaleFactor + da > tsa) {
			// Need to guard against round-off in this comparison
			// Doing that in Redshift, to match in PlanOutput()
			da = tsa - ReadState.ScaleFactor;
			STDLOG(0,"da to reach next output is %f\n", da);
		}
	}
    double finala = 1.0/(1+P.FinishingRedshift());
    if (ReadState.Redshift > P.FinishingRedshift()+1e-12 && ReadState.ScaleFactor + da > finala) {
        da = finala - ReadState.ScaleFactor;
        STDLOG(0,"da to reach finishing redshift is %f\n", da);
    }
	if (da<1e-12) return da;

	// We might have already reached the FinishingRedshift.
	if (ReadState.Redshift < P.FinishingRedshift()+1e-12) {
	    STDLOG(0,"We have reached the Finishing Redshift of %f\n", P.FinishingRedshift());
	    da = 0.0; return da;
	}

    // Do we need to output a subsample sooner?
    for(int i = 0; i < P.nTimeSliceSubsample; i++){
        double L1z = P.TimeSliceRedshifts_Subsample[i];
        double L1a = 1.0/(1+L1z);

        if(ReadState.Redshift > L1z + 1e-12 && ReadState.ScaleFactor + da > L1a){
            da = L1a - ReadState.ScaleFactor;
            STDLOG(0,"da to reach next timeslice subsample output is %f\n", da);

        }
    }

    

	// Particles should not be able to move more than one cell per timestep
	double maxdrift = cosm->DriftFactor(cosm->current.a, da)*ReadState.MaxVelocity;
	maxdrift *= P.cpd;
	STDLOG(1,"Maximum velocity would drift %f cells in this time step\n", maxdrift);
	if (maxdrift>0.8) {
	    da *= 0.8/maxdrift;   // Just linearly interpolate
	    STDLOG(0,"da based on not letting particles drift more than a cell is %f.\n", da);
	}

	// Perhaps the acceleration limits us more than this?
	// dt = eta*sqrt(epsilon/amax) is a time.  So epsilon ~ amax*(dt/eta)^2
	// Since accel beget velocities, which beget positions, we use one Kick and one Drift.
	// But this is not appropriate at early times, when particles are separated by
	// much more than a softening length!

	
	maxdrift = cosm->KickFactor(cosm->current.a, da);
	maxdrift *= cosm->DriftFactor(cosm->current.a, da);
	maxdrift *= ReadState.MaxAcceleration;
	maxdrift /= P.TimeStepAccel*P.TimeStepAccel;
	double da_eona = da;
	if (maxdrift>NFD->SofteningLength) {  // Plummer-equivalent softening length
	    if(maxdrift >1e-12) da_eona *= sqrt(NFD->SofteningLength/maxdrift);
	    STDLOG(0,"da based on sqrt(epsilon/amax) is %f.\n", da_eona);
        // We deliberately do not limit the timestep based on this criterion
	}
	

	// Perhaps the acceleration compared to the velocity is too big?
	// We want amax*dt = eta*vrms, or 
	double maxkick = cosm->KickFactor(cosm->current.a, da);
	double goal = ReadState.MinVrmsOnAmax;
	double goal2;
	if (ReadState.MaxAcceleration!=0.0) 
	    goal2 = ReadState.RMS_Velocity/ReadState.MaxAcceleration;
	else goal2 = 1e10;    // This doesn't exist in the first step.
	STDLOG(1,"Cell-based Vrms/Amax = %f\n", goal);
	STDLOG(1,"Global     Vrms/Amax = %f\n", goal2);
	// We have both a global value and a cell value.  Take the maximum of these,
	// to guard against abnormally cold cells.
	goal = max(goal,goal2) * P.TimeStepAccel;

	if (maxkick>goal) {
	    da *= goal/maxkick;
	    STDLOG(0,"da based on vrms/amax is %f. dlna = %f.\n", da, da/ReadState.ScaleFactor);
	}

    if(P.MicrostepTimeStep > 0)
        MicrostepEpochs = new MicrostepEpochTable(cosm, cosm->current.a, cosm->current.a + da, P.np);

    // Do we need to output a merger tree redshift during this step?
    for(int i = 0; i < P.nTimeSliceL1; i++){
        double L1z = P.L1OutputRedshifts[i];
        double L1a = 1.0/(1+L1z);

        // Need a bit of wiggle room to avoid duplicate outputs, but not so much that we miss outputs
        if(ReadState.Redshift > L1z + 1e-12 && ReadState.ScaleFactor + da > L1a + 1e-10){
            // We don't shorten our timestep to land exactly on a merger tree redshift.
            // Sometimes, this can lead to two consecutive group finding steps, if the L1OutputRedshift
            // we're about to say is "close enough" to the current redshift also appears in either
            // TimeSliceRedshifts or TimeSliceRedshifts_Subsample. To prevent this, check
            // if that's the case; if it is, wait to do GF and output until next step. If not, 
            // sally forth.
            STDLOG(0,"Group finding at this redshift requested by L1OutputRedshifts[%d]\n", i);
            ReadState.DoGroupFindingOutput = 1; 

            for(int i = 0; i < P.nTimeSliceSubsample; i++){
                if (L1z == P.TimeSliceRedshifts_Subsample[i]){
                    STDLOG(0,"...but will hold off, as this redshift %f appears in TimeSliceRedshifts_Subsample[%d].\n", L1z, i);
                    ReadState.DoGroupFindingOutput = 0; 
                    break;
                }
            }

            for(int i = 0; i < P.nTimeSlice; i++){
                if (L1z == P.TimeSliceRedshifts[i]){
                    STDLOG(0,"...but will hold off, as this redshift %f appears in TimeSliceRedshifts[%d].\n", L1z, i);
                    ReadState.DoGroupFindingOutput = 0;
                    break; 
                }
            }
        }
    }
	return da;
}

