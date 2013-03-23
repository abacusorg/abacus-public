#include "proepi.cpp"

#define BIGNUM 1000000.0


Cosmology *InitializeCosmology(double ScaleFactor) {
    // Be warned that all of the Cosmology routines quote time units
    // by what is entered as H0.  The code wants to use H0=1 units.
    // But P.H0 is defined to be in km/s/Mpc!
    // If you want to compare any times or H(z), remember that H0=1.
    MyCosmology cosmo;
    cosmo.Omega_m = P.Omega_M;
    cosmo.Omega_K = P.Omega_K;
    cosmo.Omega_DE = P.Omega_DE;
    cosmo.H0 = 1.0;	
    cosmo.w0 = P.w0;
    cosmo.wa = P.wa;
    STDLOG(0,"Initialized Cosmology at a= %6.4f\n",ScaleFactor);
    return new Cosmology(ScaleFactor,cosmo);
}

double ChooseTimeStep(){
	// Choose the maximum allowable timestep
	// We start with the absolute maximum timestep allowed by the parameter file,
	// then see if it needs to be shorter.

	// cosm has already been loaded with the ReadState.ScaleFactor.

	double da = ReadState.ScaleFactor*P.Dlna;
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
		double tsa = 1.0/(1+P.TimeSlicez[i]);
		if (ReadState.ScaleFactor < tsa-1e-12 && ReadState.ScaleFactor + da > tsa) {
			// Need to guard against round-off in this comparison
			da = tsa - ReadState.ScaleFactor;
			STDLOG(0,"da to reach next output is %f\n", da);
		}
	}

	if (da<1e-12) return da;

	// Particles should not be able to move more than one cell per timestep
	double maxdrift = cosm->DriftFactor(cosm->current.a, da)*ReadState.MaxVelocity;
	if (maxdrift>0.8) {
	    da *= 0.8/maxdrift;   // Just linearly interpolate
	    STDLOG(0,"da based on not letting particles drift more than a cell is %f.\n", da);
	}

	// Perhaps the acceleration limits us more than this?
	// dt = eta*sqrt(epsilon/amax) is a time.  So epsilon ~ amax*(dt/eta)^2
	// Since accel beget velocities, which beget positions, we use one Kick and one Drift.

	maxdrift = cosm->KickFactor(cosm->current.a, da);
	maxdrift *= cosm->DriftFactor(cosm->current.a, da);
	maxdrift *= ReadState.MaxAcceleration;
	maxdrift /= P.Eta*P.Eta;
	if (maxdrift>P.SofteningLength) {
	    da *= sqrt(P.SofteningLength/maxdrift);
	    STDLOG(0,"da based on sqrt(epsilon/amax) is %f.\n", da);
	}

	// Perhaps the acceleration compared to the velocity is too big?
	// We want amax*dt = eta*vrms, or 
	double maxkick = cosm->KickFactor(cosm->current.a, da);
	double goal = ReadState.MinVrmsOnAmax;
	// double goal = P.TimeStepVonA*ReadState.MinVrmsOnAmax;
	STDLOG(0,"da based on vrms/amax would differ by %f\n", goal/maxkick);

	/* 
	if (maxkick>goal) {
	    da *= goal/maxkick;
	    STDLOG(0,"da based on vrms/amax is %f.\n", da_acc2);
	}
	*/

	return da;
}


void BuildWriteState(double da){
	STDLOG(0,"Building WriteState for a step from a=%f by da=%f\n", cosm->current.a, da);

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

	//fill in WriteState from the Parameter file
	WriteState.np_state =P.np;
	WriteState.cpd_state = P.cpd;
	WriteState.order_state = P.order;
	WriteState.ppd = P.ppd();
	WriteState.DoublePrecision = (sizeof(FLOAT)==8)?1:0;

	//get the next timestep and build the cosmology for it
	double nexta = cosm->current.a + da;
	STDLOG(0,"Next scale factor is %f\n", nexta);
	cosm->BuildEpoch(cosm->current, cosm->next, nexta);

	WriteState.ScaleFactor = cosm->next.a;
	WriteState.Redshift = cosm->next.z;
	WriteState.Time = cosm->next.t;                // In Gyr or Gyr/h, depending on hMpc flag
	WriteState.etaK = cosm->next.etaK;
	WriteState.etaD = cosm->next.etaD;
	WriteState.Growth = cosm->next.growth;
	WriteState.Growth_on_a = cosm->next.growth/cosm->next.a;
	WriteState.f_growth = cosm->next.f_growth;
	WriteState.w = cosm->next.w;
	WriteState.HubbleNow = cosm->next.H;           // In km/s/Mpc
	WriteState.Htime = cosm->next.H*cosm->next.t;               // Time*H(z)

	double total = cosm->next.OmegaHat_m+cosm->next.OmegaHat_X+cosm->next.OmegaHat_K;
	WriteState.OmegaNow_m = cosm->next.OmegaHat_m/total;
	WriteState.OmegaNow_K = cosm->next.OmegaHat_K/total;
	WriteState.OmegaNow_DE = cosm->next.OmegaHat_DE/total;

	WriteState.ParticleMass = ReadState.ParticleMass; //FIXME: This is just a place holder // In Msun or Msun/h, depending on hMpc flag
	WriteState.RedshiftSpaceConversion = ReadState.RedshiftSpaceConversion ;//FIXME: Another placeholder until the actual math is worked out
	WriteState.LPTstatus = ReadState.LPTstatus; //TODO: Depricated?
	WriteState.FullStepNumber = ReadState.FullStepNumber+1;


	WriteState.DeltaTime = cosm->next.t - cosm->current.t;
	WriteState.DeltaRedshift = -(cosm->next.z - cosm->current.z);
	WriteState.DeltaScaleFactor = cosm->next.a - cosm->current.a;

	// Might also compute the Scale Factor for the halfway time, since
	// we'll need that.
	cosm->t2a(0.5*(cosm->next.t+cosm->current.t));
	STDLOG(0,"Scale factor halfway in between is %f\n", cosm->search.a);
	// cosm->search now has the midpoint epoch.
	WriteState.ScaleFactorHalf = cosm->search.a;
	WriteState.LastHalfEtaKick = cosm->KickFactor(cosm->search.a,WriteState.ScaleFactor-cosm->search.a);
	WriteState.FirstHalfEtaKick = cosm->KickFactor(cosm->current.a,cosm->search.a-cosm->current.a);
	WriteState.DeltaEtaDrift = cosm->DriftFactor(cosm->current.a, cosm->next.a-cosm->current.a);
		// cosm->next.etaD - cosm->current.etaD;
		// WriteState.etaD - ReadState.etaD;

	// Just truncate some underflow cases
	if (fabs(WriteState.DeltaEtaDrift)   <1e-14*fabs(cosm->current.etaD)) WriteState.DeltaEtaDrift = 0.;
	if (fabs(WriteState.FirstHalfEtaKick)<1e-14*fabs(cosm->current.etaK)) WriteState.FirstHalfEtaKick = 0.;
	if (fabs(WriteState.LastHalfEtaKick) <1e-14*fabs(cosm->current.etaK)) WriteState.LastHalfEtaKick = 0.;

	// Some statistics to accumulate
	WriteState.MaxCellSize = ReadState.MaxCellSize;
	WriteState.StdDevCellSize = 0;

	// Build the output header
	WriteState.make_output_header();
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
    stdlog_threshold_global = P.LogVerbosity;
    char logfn[1050];
    sprintf(logfn,"%s/lastrun.log", P.LogFileDirectory);
    stdlog.open(logfn);
    STDLOG_TIMESTAMP;
    STDLOG(0,"Read Parameter file %s\n", argv[1]);

    strcpy(WriteState.ParameterFileName, argv[1]);
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

    	}
    	else{
    		STDLOG(0,"Generating initial State from initial conditions\n");
    		ReadState.ParticleMass = 1.0/P.np; //FIXME: This is just a place holder // In Msun or Msun/h, depending on hMpc flag
    		ReadState.RedshiftSpaceConversion = 1.0 ;//FIXME: Another placeholder until the actual math is worked out
    		ReadState.LPTstatus = P.LagrangianPTOrder;
    		ReadState.FullStepNumber = -1;  
			// So that this number is the number of times
			// forces have been computed.
    		sprintf(ReadState.ParameterFileName,"%s",argv[1]);
        	ReadState.ScaleFactor = 1.0/(1+P.InitialRedshift);
        	da = 0;
        	MakeIC = true;
    	}
    }
    else{
    	CheckDirectoryExists(P.ReadStateDirectory);
	STDLOG(0,"Reading ReadState\n");
    	ReadState.read_from_file(P.ReadStateDirectory);

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
    	STDLOG(0,"Read ReadState from %s\n",P.ReadStateDirectory);

	if (P.ForceOutputDebug==1) {
	    STDLOG(0,"ForceOutputDebug option invoked; setting time step to 0.\n");
	    da = 0;
	}
    	MakeIC = false;
    }

    //Check if WriteStateDirectory/state exists, and fail if it does
    char wstatefn[1050];
    sprintf(wstatefn,"%s/state",P.WriteStateDirectory);
    if(access(wstatefn,0) !=-1)
    	QUIT("WriteState exists and would be overwritten. Please move or delete it to continue.\n");

    cosm = InitializeCosmology(ReadState.ScaleFactor);
    if (da!=0) da = ChooseTimeStep();
    STDLOG(0,"Chose Time Step da = %6.4f\n",da);

    feenableexcept(FE_INVALID | FE_DIVBYZERO);
    BuildWriteState(da);

    SingleStepSetup.Stop();

    Prologue(P,MakeIC);

    if (MakeIC)  timestepIC();
	    else timestep();

    // Let the IO finish, so that it is included in the time log.
    IO_Terminate();

    // The timings need to proceed the epilogue, because they need to look inside 
    // some instances of classes.
    WallClockDirect.Stop();

    fedisableexcept(FE_INVALID | FE_DIVBYZERO);

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
