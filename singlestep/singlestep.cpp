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

    S.HubbleTime = 1.0; 	
    	// The code always uses H_0=1, so here's the true 1/H_0 in Gyr or Gyr/h, depending on hMpc flag
	// TODO: Not yet in place.
    S.ParticleMass = 1.0/P.np;
    	//FIXME: This is just a place holder // In Msun or Msun/h, depending on hMpc flag

    // The code uses canonical velocities, for the unit box and 1/H_0 time unit.
    // However, our output standard is redshift-space comoving displacements, again
    // in units where the box is one.
    // The conversion is v_canon = v_zspace * a^2 H(z)/H_0
    S.VelZSpace_to_Canonical = S.ScaleFactor * S.ScaleFactor * (S.HubbleNow/cosm->C.H0);

    // After output, one might want to convert to km/s.  
    // Here, we quote the number of km/s across the full box.
    S.RedshiftSpaceConversion = 1.0;
    	//FIXME: Another placeholder until the actual math is worked out
	// What are we doing about the hMpc flag?
}



double ChooseTimeStep(){
	//choose the maximum allowable timestep
	//we start with the absolute maximum timestep allowed by the parameter file
	// cosm has already been loaded with the ReadState.ScaleFactor.

	double da_max = ReadState.ScaleFactor*P.Dlna;
	// TODO: I think below might be simplified if we used the kickfactor and
	// driftfactor methods from the cosmology function.  Also, we need to be
	// careful that our velocities and accelerations are in code (canonical) units,
	// so there are redshift factors slinging around.  Using kicks and drifts
	// helps to avoid problems.

	//first we calculate the maximum timesteps in time units, then choose the minimum and convert to da

	//limit imposed by acceleration
	double dt_acc = P.Eta * sqrt(P.SofteningLength/ReadState.MaxAcceleration);
	if (isnan(dt_acc)) dt_acc = BIGNUM;
	//particles should move only one cell per timestep
	double dt_v = 1.0/P.cpd *1.0/(ReadState.MaxVelocity);
	if (isnan(dt_v)) dt_v = BIGNUM;
	if (ReadState.MaxAcceleration <= 0) dt_acc = 10000000.0;
	//limit imposed by acceleration on velocity
	/*double dt_vona = .5 * ReadState.MinVrmsOnAmax;
	if (isnan(dt_vona)) dt_vona = BIGNUM;
	*/
	double dt_vona  = BIGNUM;
//	printf("dt_acc: %f \t dt_v: %f \t dt_vona: %f \n ", dt_acc,dt_v,dt_vona);
	double dtmin = min(min(dt_acc,dt_v),dt_vona);
	double da;
	if (dtmin > 100000.0) da = BIGNUM;
	else da = cosm->t2a(ReadState.Time+dtmin) - ReadState.ScaleFactor; //TODO: This is essentially from Marc's code. Is it kosher?

	da = min(da,da_max);

	for (int i = 0; i < P.nTimeSlice; i ++){
		double tsa = 1.0/(1+P.TimeSlicez[i]);
		if (ReadState.ScaleFactor < tsa && ReadState.ScaleFactor + da > tsa){
			da = tsa - ReadState.ScaleFactor;
			break;
		}
	}

	assertf(da >= 0,"Maximum da:%f was not > 0. dt_acc: %f \t dt_v: %f \n",da,dt_acc,dt_v);
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
	WriteState.FullStepNumber = ReadState.FullStepNumber+1;

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
	WriteState.StdDevCellSize = 0.0;
	WriteState.MaxVelocity = 0.0;
	WriteState.MaxAcceleration = 0.0;
	WriteState.rms_velocity = 0.0;
	WriteState.MinVrmsOnAmax = 1e10;

	// Build the output header
	WriteState.make_output_header();
}

void PlanOutput() {
    // Check the time slice and decide whether to do output.
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
    STDLOG(0,"Chose Time Step da = %6.4f\n",da);
    feenableexcept(FE_INVALID | FE_DIVBYZERO);
    BuildWriteState(da);

    // Make a plan for output
    PlanOutput();

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
