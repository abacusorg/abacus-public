#include "proepi.cpp"

#define BIGNUM 1000000.0


Cosmology * InitializeCosmology(double ScaleFactor) {
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
    return new Cosmology(ScaleFactor,cosmo);
}

double ChooseTimeStep(){
	//choose the maximum allowable timestep
	//we start with the absolute maximum timestep allowed by the parameter file

	MyCosmology cosmo;
	cosmo.Omega_m = P.Omega_M;
	cosmo.Omega_K = P.Omega_K;
	cosmo.Omega_DE = P.Omega_DE;
	cosmo.H0 = 1.0;
	cosmo.w0 = P.w0;
	cosmo.wa = P.wa;
	cosm =  new Cosmology(ReadState.ScaleFactor,cosmo);

	double da_max = ReadState.ScaleFactor*P.Dlna;

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
	delete cosm;
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
	WriteState.np =P.np;
	WriteState.cpd = P.cpd;
	WriteState.order = P.order;
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

    double a;
    double da;
    bool MakeIC; //True if we should make the initial state instead of doing a real timestep
    //Check if ReadStateDirectory is accessible, or if we should build a new state from the IC file

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
        	a = 1.0/(1+P.InitialRedshift);
        	da = 0;
        	MakeIC = true;
    	}
    }
    else{
    	CheckDirectoryExists(P.ReadStateDirectory);
	STDLOG(0,"Reading ReadState\n");
    	readstate(ReadState,P.ReadStateDirectory);
	// Strange :: to stdout during this step.

    	//make sure read state and parameters are compatible
    	assertf(ReadState.order == P.order, 
		"ReadState and Parameter order do not match, %d != %d\n", 
		ReadState.order, P.order);
    	assertf(ReadState.cpd == P.cpd, 
		"ReadState and Parameter cpd do not match, %d != %d\n", 
		ReadState.cpd, P.cpd);
    	assertf(ReadState.np == P.np, 
		"ReadState and Parameter np do not match, %d != %d\n", 
		ReadState.np, P.np);

    	STDLOG(0,"Read ReadState from %s\n",P.ReadStateDirectory);
    	a = ReadState.ScaleFactor;
    	da = ChooseTimeStep();
	if (P.ForceOutputDebug==1) {
	    STDLOG(0,"ForceOutputDebug option invoked; setting time step to 0.\n");
	    da = 0;
	}
    	MakeIC = false;
    }
    feenableexcept(FE_INVALID | FE_DIVBYZERO);

    cosm = InitializeCosmology(a);
    STDLOG(0,"Initialized Cosmology at a= %6.4f\n",a);

    //Check if WriteStateDirectory/state exists, and fail if it does
    char wstatefn[1050];
    sprintf(wstatefn,"%s/state",P.WriteStateDirectory);
    if(access(wstatefn,0) !=-1)
    	QUIT("WriteState exists and would be overwritten. Please move or delete it to continue.\n");


    STDLOG(0,"Chose Time Step da = %6.4f\n",da);
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
    writestate(&WriteState,P.WriteStateDirectory);
    STDLOG(0,"Wrote WriteState to %s\n",P.WriteStateDirectory);

    stdlog.close();  
    exit(0);
}
