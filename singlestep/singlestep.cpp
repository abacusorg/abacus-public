#include "proepi.cpp"

#define BIGNUM 1000000.0
#define DATOLERANCE 1.e-12


void AssertStateLegal() {
    //make sure read state and parameters are compatible
    assertf(ReadState.order_state == P.order, 
	    "ReadState and Parameter order do not match, %d != %d\n", 
	    ReadState.order_state, P.order);
    assertf(ReadState.cpd_state == P.cpd, 
	    "ReadState and Parameter cpd do not match, %d != %d\n", 
	    ReadState.cpd_state, P.cpd);
    assertf(ReadState.np_state == P.np, 
	    "ReadState and Parameter np do not match, %d != %d\n", 
	    ReadState.np_state, P.np);
    assertf(ReadState.MaxCellSize < (2.048e9)/3,
	    "The largest cell has %d particles, exceeding the allowed amount.\n",
	    ReadState.MaxCellSize);
}

Cosmology *InitializeCosmology(double ScaleFactor) {
    // Be warned that all of the Cosmology routines quote time units
    // by what is entered as H0.  The code wants to use H0=1 units.
    // But P.H0 is defined to be in km/s/Mpc!
    // If you want to compare any times or H(z), remember that H0=1.
    // This routine should only be called once.
    MyCosmology cosmo;
    cosmo.Omega_m = P.Omega_M;
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
    S.Growth_on_a = cosm->next.growth/cosm->next.a;
    S.f_growth = cosm->next.f_growth;
    S.w = cosm->next.w;
    S.HubbleNow = cosm->next.H;           // In same units as H_0 was given
    S.Htime = cosm->next.H*cosm->next.t;               // Time*H(z), in code units

    double total = cosm->next.OmegaHat_m+cosm->next.OmegaHat_X+cosm->next.OmegaHat_K;
    S.OmegaNow_m = cosm->next.OmegaHat_m/total;
    S.OmegaNow_K = cosm->next.OmegaHat_K/total;
    S.OmegaNow_DE = cosm->next.OmegaHat_X/total;

    S.HubbleTimeHGyr = 9.7782;
	// The Hubble Time is 9.7782 h^-1 Gyr
    S.HubbleTimeGyr = S.HubbleTimeHGyr*(100.0/P.H0);
    	// This is in Gyr
    S.BoxSizeMpc = P.BoxSize*(P.hMpc?(100.0/P.H0):1.0);
    S.BoxSizeHMpc = S.BoxSizeMpc*(P.H0/100.0);
    	// Redundant, but we might as well be explicit.
    S.ParticleMassMsun = 2.7746e11*P.Omega_M*pow(P.H0/100,2.0)*pow(S.BoxSizeMpc,3.0)/P.np;
    	// This is in Msun.  
	// The critical density is 2.7746e11 h^2 Msun/Mpc^3
    S.ParticleMassHMsun = S.ParticleMassMsun*(P.H0/100);
    	// This is in h^-1 Msun.  

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
    S.RedshiftSpaceConversion = P.BoxSize*S.ScaleFactor*(S.HubbleNow/cosm->C.H0)*
    	(P.hMpc?100:P.H0);
}



double ChooseTimeStep(){
	// Choose the maximum allowable timestep
	// We start with the absolute maximum timestep allowed by the parameter file,
	// then see if it needs to be shorter.

	// cosm has already been loaded with the ReadState.ScaleFactor.


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
	if (da<1e-12) return da;

	// We might have already reached the FinishingRedshift.
	if (ReadState.Redshift < P.FinishingRedshift()+1e-12) {
	    STDLOG(0,"We have reached the Finishing Redshift of %f\n", P.FinishingRedshift());
	    da = 0.0; return da;
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

	/*
	maxdrift = cosm->KickFactor(cosm->current.a, da);
	maxdrift *= cosm->DriftFactor(cosm->current.a, da);
	maxdrift *= ReadState.MaxAcceleration;
	maxdrift /= P.TimeStepAccel*P.TimeStepAccel;
	if (maxdrift>P.SofteningLength) {
	    da *= sqrt(P.SofteningLength/maxdrift);
	    STDLOG(0,"da based on sqrt(epsilon/amax) is %f.\n", da);
	}
	*/

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

	return da;
}



void BuildWriteState(double da){
	STDLOG(0,"Building WriteState for a step from a=%f by da=%f\n", cosm->current.a, da);

	// fill in WriteState from the Parameter file
	WriteState.np_state = P.np;
	WriteState.cpd_state = P.cpd;
	WriteState.order_state = P.order;
	WriteState.ppd = P.ppd();

	// Fill in the logistical reporting fields
#ifdef GITVERSION	
	STDLOG(0,"Git Hash = %s\n", GITVERSION);
	strncpy(WriteState.CodeVersion, GITVERSION, 1024);
#endif
	time_t timet = time(0);
	string now = string(asctime(localtime(&timet)));
	sprintf(WriteState.RunTime,"%s",now.substr(0,now.length()-1).c_str());
	gethostname(WriteState.MachineName,1024);
	STDLOG(0,"Host machine name is %s\n", WriteState.MachineName);

	WriteState.DoublePrecision = (sizeof(FLOAT)==8)?1:0;
	STDLOG(0,"Bytes per float is %d\n", (WriteState.DoublePrecision+1)*4);
	STDLOG(0,"Bytes per auxstruct is %d\n", sizeof(auxstruct));
	STDLOG(0,"Bytes per cellinfo is %d\n", sizeof(cellinfo));
	WriteState.FullStepNumber = ReadState.FullStepNumber+1;
	STDLOG(0,"This is step number %d\n", WriteState.FullStepNumber);

	//get the next timestep and build the cosmology for it
	double nexta = cosm->current.a + da;
	STDLOG(0,"Next scale factor is %f\n", nexta);
	cosm->BuildEpoch(cosm->current, cosm->next, nexta);
	FillStateWithCosmology(WriteState);

	// Get differences between the two cosmological epochs
	WriteState.DeltaTime = cosm->next.t - cosm->current.t;
	WriteState.DeltaRedshift = -(cosm->next.z - cosm->current.z);
	WriteState.DeltaScaleFactor = cosm->next.a - cosm->current.a;

	cosm->t2a(0.5*(cosm->next.t+cosm->current.t));
	STDLOG(0,"Scale factor halfway in between is %f\n", cosm->search.a);
	// cosm->search now has the midpoint epoch.
	WriteState.ScaleFactorHalf = cosm->search.a;
	WriteState.LastHalfEtaKick = 
		cosm->KickFactor(cosm->search.a,WriteState.ScaleFactor-cosm->search.a);
	WriteState.FirstHalfEtaKick = 
		cosm->KickFactor(cosm->current.a,cosm->search.a-cosm->current.a);
	WriteState.DeltaEtaDrift = 
		cosm->DriftFactor(cosm->current.a, cosm->next.a-cosm->current.a);

	// Just truncate some underflow cases
	if (fabs(WriteState.DeltaEtaDrift)   <1e-14*fabs(cosm->current.etaD)) 
		WriteState.DeltaEtaDrift = 0.;
	if (fabs(WriteState.FirstHalfEtaKick)<1e-14*fabs(cosm->current.etaK)) 
		WriteState.FirstHalfEtaKick = 0.;
	if (fabs(WriteState.LastHalfEtaKick) <1e-14*fabs(cosm->current.etaK)) 
		WriteState.LastHalfEtaKick = 0.;

	// Initialize some statistics to accumulate
	WriteState.MaxCellSize = 0;
	WriteState.MinCellSize = 1e9;
	WriteState.StdDevCellSize = 0.0;
	WriteState.MaxVelocity = 0.0;
	WriteState.MaxAcceleration = 0.0;
	WriteState.RMS_Velocity = 0.0;
	WriteState.MinVrmsOnAmax = 1e10;

	// Build the output header.
	// Note we actually will output from ReadState,
	// but we build this to write the write/state file from the same code.
	WriteState.make_output_header();
}

void PlanOutput(bool MakeIC) {
    // Check the time slice and decide whether to do output.
    ReadState.DoTimeSliceOutput = 0;
    ReadState.OutputIsAllowed = 0;
    if (MakeIC) return;   // Do no output on this slice.
    if (LPTStepNumber()>0) return;  // We're doing IC work; no output

    // Build the output header.  The cosmology is from ReadState,
    // but we'd like to use some elements from WriteState.  So we 
    // overwrite some ReadState elements.
    // We will call this output by the WriteState FullStepNumber, as that
    // is what is required for the velocities (and is the run writing the output). 
    strncpy(ReadState.ParameterFileName, WriteState.ParameterFileName, 1024);
    strncpy(ReadState.CodeVersion, WriteState.CodeVersion, 1024);
    strncpy(ReadState.MachineName, WriteState.MachineName, 1024);
    strncpy(ReadState.RunTime,     WriteState.RunTime, 1024);
    ReadState.FullStepNumber = WriteState.FullStepNumber;
    ReadState.make_output_header();

    // Just let later routines know that this is a valid epoch
    // for output, e.g., not a LPT IC epoch.
    ReadState.OutputIsAllowed = 1;

    // Now check whether we're asked to do a TimeSlice.
    for (int nn = 0; nn < P.nTimeSlice; nn++) {
	if (fabs(ReadState.Redshift-P.TimeSliceRedshifts[nn])<1e-12) {
	    STDLOG(0,"Planning to output a TimeSlice, element %d\n", nn);
	    ReadState.DoTimeSliceOutput = 1;
	    char slicedir[128];
	    sprintf(slicedir,"slice%5.3f", ReadState.Redshift);
	    CreateSubDirectory(P.OutputDirectory,slicedir);
	    break;
	}
    }
}


int main(int argc, char **argv) {

    WallClockDirect.Start();
    SingleStepSetup.Start();

    if (argc!=3) {
       // Can't use assertf() or QUIT here: stdlog not yet defined!
       fprintf(stderr, "singlestep(): command line must have 3 parameters given, not %d.\nLegal usage: singlestep <parameter_file> <allow creation of initial conditions 1/0>\n", argc);
       assert(0==99);
    }
    
    int AllowIC = atoi(argv[2]);
    P.ReadParameters(argv[1],1);
    strcpy(WriteState.ParameterFileName, argv[1]);

    // Setup the log
    stdlog_threshold_global = P.LogVerbosity;
    char logfn[1050];
    sprintf(logfn,"%s/lastrun.log", P.LogFileDirectory);
    stdlog.open(logfn);
    STDLOG_TIMESTAMP;
    STDLOG(0,"Read Parameter file %s\n", argv[1]);
    STDLOG(0,"AllowIC = %d\n", AllowIC);

    double da = -1.0;   // If we set this to zero, it will skip the timestep choice
    bool MakeIC; //True if we should make the initial state instead of doing a real timestep

    // Check if ReadStateDirectory is accessible, or if we should 
    // build a new state from the IC file
    char rstatefn[1050];
    sprintf(rstatefn,"%s/state",P.ReadStateDirectory);

    if(access(rstatefn,0) ==-1){
	STDLOG(0,"Can't find ReadStateDirectory %s\n", P.ReadStateDirectory);
    	if(AllowIC != 1){
	    QUIT("Read State Directory ( %s ) is inaccessible and initial state creation is prohibited. Terminating.\n",P.ReadStateDirectory);

    	} else{
	    STDLOG(0,"Generating initial State from initial conditions\n");
	    // We have to fill in a few items, just to bootstrap the rest of the code.
	    ReadState.ScaleFactor = 1.0/(1+P.InitialRedshift);
	    ReadState.FullStepNumber = -1;  
		// So that this number is the number of times forces have been computed.
		// The IC construction will yield a WriteState that is number 0,
		// so our first time computing forces will read from 0 and write to 1.
	    da = 0;
	    MakeIC = true;
	}
    } else {
	// We're doing a normal step
    	CheckDirectoryExists(P.ReadStateDirectory);
    	STDLOG(0,"Reading ReadState from %s\n",P.ReadStateDirectory);
    	ReadState.read_from_file(P.ReadStateDirectory);
	AssertStateLegal();
    	MakeIC = false;
	// Handle some special cases
	if (P.ForceOutputDebug==1) {
	    STDLOG(0,"ForceOutputDebug option invoked; setting time step to 0.\n");
	    da = 0;
	}
    }

    //Check if WriteStateDirectory/state exists, and fail if it does
    char wstatefn[1050];
    sprintf(wstatefn,"%s/state",P.WriteStateDirectory);
    if(access(wstatefn,0) !=-1)
    	QUIT("WriteState exists and would be overwritten. Please move or delete it to continue.\n");

    // Initialize the Cosmology and set up the State epochs and the time step
    cosm = InitializeCosmology(ReadState.ScaleFactor);
    if (MakeIC) FillStateWithCosmology(ReadState);
    if (da!=0) da = ChooseTimeStep();
    STDLOG(0,"Chose Time Step da = %6.4f, dlna = %6.4f\n",da, da/ReadState.ScaleFactor);
    feenableexcept(FE_INVALID | FE_DIVBYZERO);
    BuildWriteState(da);

    // Make a plan for output
    PlanOutput(MakeIC);

    SingleStepSetup.Stop();

    // Now execute the timestep
    Prologue(P,MakeIC);
    if (MakeIC)  timestepIC();
	    else timestep();

    // Let the IO finish, so that it is included in the time log.
    IO_Terminate();
    fedisableexcept(FE_INVALID | FE_DIVBYZERO);

    // Write out the timings.  This must precede the epilogue, because 
    // we need to look inside some instances of classes for runtimes.
    WallClockDirect.Stop();
    if (!MakeIC){
    	char timingfn[1050];
    	sprintf(timingfn,"%s/lastrun.steptiming", P.LogFileDirectory);
    	FILE * timingfile = fopen(timingfn,"w");
    	ReportTimings(timingfile);
    	STDLOG(0,"Wrote Timing File to %s\n",timingfn);
    }

    // The epilogue contains some tests of success.
    Epilogue(P,MakeIC);
    delete cosm;

    // The state should be written last, since that officially signals success.
    WriteState.StdDevCellSize = sqrt(WriteState.StdDevCellSize);
    WriteState.write_to_file(P.WriteStateDirectory);
    STDLOG(0,"Wrote WriteState to %s\n",P.WriteStateDirectory);

    stdlog.close();  
    exit(0);
}
