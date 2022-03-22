#include "proepi.cpp"
#include "timestep_ic.cpp"
#include "cosmo_setup.cpp"


void BuildWriteState(double da){
	STDLOG(0,"Building WriteState for a step from a=%f by da=%f\n", cosm->current.a, da);

	// fill in WriteState from the Parameter file
	WriteState.np_state = P.np;
	WriteState.cpd_state = P.cpd;
	WriteState.order_state = P.order;
	WriteState.ppd = P.ppd();
    WriteState.ippd = (int64) round(WriteState.ppd);
        
	// Fill in the logistical reporting fields
#ifdef GITVERSION
	STDLOG(0,"Git Hash = %s\n", GITVERSION);
	strncpy(WriteState.CodeVersion, GITVERSION, 1024);
#endif
	strncpy(WriteState.OutputFormatVersion, "v2.2", 1024);
        // v2.0 is the first versioned format (so named because this was the Abacus release version at the time)
        // v2.0 uses RVint with the full 2**20 range; halo_info at 296 bytes.
        // v2.2 uses RVint with a 1,000,000 range; halo_info at 296 bytes.
	time_t timet = time(0);
	string now = string(asctime(localtime(&timet)));
	sprintf(WriteState.RunTime,"%s",now.substr(0,now.length()-1).c_str());
	gethostname(WriteState.MachineName,1024);
    WriteState.NodeRank = MPI_rank;
    WriteState.NodeSize = MPI_size;

	WriteState.DoublePrecision = (sizeof(FLOAT)==8)?1:0;
	STDLOG(0,"Bytes per float is %d\n", sizeof(FLOAT));
	STDLOG(0,"Bytes per auxstruct is %d\n", sizeof(auxstruct));
	STDLOG(0,"Bytes per cellinfo is %d\n", sizeof(cellinfo));
	assert(WriteState.FullStepNumber == ReadState.FullStepNumber+1);  // already set this in load_read_state()
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

    WriteState.GhostRadius = MERGE_GHOST_RADIUS;
}

void BuildWriteStateOutput() {
	// Build the output header.
	// Note we actually will output from ReadState,
	// but we build this to write the write/state file from the same code.
	WriteState.make_output_header();
}

double EvolvingDelta(float z){
    float omegaMz = P.Omega_M * pow(1.0 + z, 3.0) / (P.Omega_DE + P.Omega_M *  pow(1.0 + z, 3.0) );
	// The Bryan & Norman (1998) threshold is in units of (redshift-dependent) critical density;
	// the 1/omegaMz factor makes it relative to the mean density at that redshift
    float Deltaz = (18.0*M_PI*M_PI + 82.0 * (omegaMz - 1.0) - 39.0 * pow(omegaMz - 1.0, 2.0) ) / omegaMz;
    return Deltaz / (18.0*M_PI*M_PI); //Params are given at high-z, so divide by high-z asymptote to find rescaling.
}

void PlanOutput(bool MakeIC) {
    // Check the time slice and decide whether to do output.
    ReadState.DoTimeSliceOutput = 0;
    ReadState.OutputIsAllowed = 0;
    ReadState.DoBinning = 0;
    density = 0;
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


    #ifdef SPHERICAL_OVERDENSITY
    if (P.SO_EvolvingThreshold) {

        float rescale = EvolvingDelta(ReadState.Redshift);

        STDLOG(2, "Rescaling SO Delta as a function of redshift.\n\t\tL0: %f --> %f\n\t\tL1: %f --> %f\n\t\tL2: %f --> %f\n",
                      P.L0DensityThreshold, rescale * P.L0DensityThreshold,
                      P.SODensity[0], rescale * P.SODensity[0],
                      P.SODensity[1], rescale * P.SODensity[1]);


        P.SODensity[0] *= rescale;
        P.SODensity[1] *= rescale;
        P.L0DensityThreshold *= rescale;

        ReadState.SODensityL1 = P.SODensity[0];
        ReadState.SODensityL2 = P.SODensity[1];
        ReadState.L0DensityThreshold = P.L0DensityThreshold;
    }
    else STDLOG(2, "Using constant SO Delta (no redshift-dependent rescaling).\n");
    #endif


    ReadState.make_output_header();

    // Just let later routines know that this is a valid epoch
    // for output, e.g., not a LPT IC epoch.
    ReadState.OutputIsAllowed = 1;

    // Now check whether we're asked to do a TimeSlice.
    for (int nn = 0; nn < P.nTimeSlice || P.OutputEveryStep == 1; nn++) {
	if ((abs(ReadState.Redshift-P.TimeSliceRedshifts[nn])<1e-12) || P.OutputEveryStep == 1) {
	    STDLOG(0,"Planning to output a TimeSlice, element %d\n", nn);
	    ReadState.DoTimeSliceOutput = 1;
	    char slicedir[128];
	    sprintf(slicedir,"slice%5.3f", ReadState.Redshift);
	    CreateSubDirectory(P.OutputDirectory,slicedir);
	    break;
	}
    }

    //check if we should bin
    if (P.PowerSpectrumStepInterval >0 && ReadState.FullStepNumber %P.PowerSpectrumStepInterval == 0){
    	density = new FLOAT[P.PowerSpectrumN1d*P.PowerSpectrumN1d*P.PowerSpectrumN1d];
    	ReadState.DoBinning = 1;
    }

}


int main(int argc, char **argv) {
    //signal(SIGUSR1, graceful_exit_signal_handler);

    SET_LOG_REFERENCE_TIME;  // Establish the base timestamp
        // but logging still not available until setup_log() below!

    //Enable floating point exceptions
    feenableexcept(FE_INVALID | FE_DIVBYZERO);

    WallClockDirect.Start();
    SingleStepSetup.Start();

    if (argc!=3) {
       // Can't use assertf() or QUIT here: stdlog not yet defined!
       fprintf(stderr, "singlestep(): command line must have 3 parameters given, not %d.\nLegal usage: singlestep PARAMETER_FILE MAKE_IC\n\tPARAMETER_FILE: path to parameter file (usually called abacus.par)\n\tMAKE_IC: 0 or 1, whether this is an IC step.\n", argc);
       return 1;
    }

    // Set up MPI
    StartMPI();  // call MPI_Init as early as possible
    P.ReadParameters(argv[1],0);

    int MakeIC = atoi(argv[2]);
    load_read_state(MakeIC);  // Load the bare minimum (step num, etc) to get the log up and running

    setup_log(); // STDLOG and assertf now available
    STDLOG(0,"Read Parameter file %s\n", argv[1]);
    STDLOG(0,"MakeIC = %d\n", MakeIC);

    InitializePipelineWidths();
    InitializeParallel();  // MPI_rank now available

    SetupLocalDirectories(MakeIC);

    // Set up OpenMP
    init_openmp();

    // Decide what kind of step to do
    double da = -1.0;   // If we set this to zero, it will skip the timestep choice

    check_read_state(MakeIC, da);

    // Initialize the Cosmology and set up the State epochs and the time step
    cosm = InitializeCosmology(ReadState.ScaleFactor);
    if (MakeIC) FillStateWithCosmology(ReadState);

    // Set some WriteState values before ChooseTimeStep()
    // This also sets the SofteningLength, needed by the NFD constructor
    InitWriteState(MakeIC, "singlestep", argv[1]);

    // Set up the major classes (including NFD)
    Prologue(P,MakeIC);

    // Check if WriteStateDirectory/state exists, and fail if it does
    char wstatefn[1050];
    int ret = snprintf(wstatefn, 1050, "%s/state", P.WriteStateDirectory);
    assert(ret >= 0 && ret < 1050);
    if(access(wstatefn,0) !=-1 && !WriteState.OverwriteState)
        QUIT("WriteState \"%s\" exists and would be overwritten. Please move or delete it to continue.\n", wstatefn);

    if (da!=0) da = ChooseTimeStep();

    // da *= -1;  // reverse the time step TODO: make parameter
    double dlna = da/ReadState.ScaleFactor;
    STDLOG(0,"Chose Time Step da = %6.4f, dlna = %6.4f\n", da, dlna);
    if(dlna > 0){
        STDLOG(0, "\t\tAt the current rate, this implies %d more steps to z_final=%f\n", (int64)ceil(log(1./ReadState.ScaleFactor/(1. + P.FinishingRedshift()))/dlna), P.FinishingRedshift());
    }

    BuildWriteState(da);
    assertf(P.NLightCones<=NUMLIGHTCONES, "Parameter file requests %d light cones, but AUX data model supports only %d\n", P.NLightCones, NUMLIGHTCONES);
    LCOrigin = (double3 *) malloc(NUMLIGHTCONES*sizeof(double3));  // max NUMLIGHTCONES light cones
    // TODO: Why is this NUMLIGHTCONES instead of P.NLightCones?
    for(int i = 0; i < NUMLIGHTCONES; i++)
        LCOrigin[i] = ((double3*) P.LightConeOrigins)[i]/P.BoxSize;  // convert to unit-box units

    // Make a plan for output
    PlanOutput(MakeIC);

    // Set up the Group Finding concepts and decide if Group Finding output is requested.
    InitGroupFinding(MakeIC);
    BuildWriteStateOutput();    // Have to delay this until after GFC is made

    SingleStepSetup.Stop();

    // Now execute the timestep
    if (MakeIC)  timestepIC();
	    else timestep();

    fedisableexcept(FE_INVALID | FE_DIVBYZERO);

    // Let the IO finish, so that it is included in the time log.
    IOFinish.Start();
    IO_Terminate();
    IOFinish.Stop();
    WallClockDirect.Stop();

    SingleStepTearDown.Start();

    // While we still have global objects, log their timings
    GatherTimings();

    // The epilogue contains some tests of success.
    Epilogue(P,MakeIC);

    delete cosm;
    free(LCOrigin);

    // Print out some final stats
    FinalizeWriteState();

    // The state should be written last, since that officially signals success.
    if (WriteState.NodeRank==0) {
        WriteState.write_to_file(P.WriteStateDirectory);
        STDLOG(0,"Wrote WriteState to %s\n",P.WriteStateDirectory);
    }

    if (!MakeIC && P.ProfilingMode){
        STDLOG(0,"ProfilingMode is active. Removing the write state in %s\n",P.LocalWriteStateDirectory);
        char command[1024];
        int printret = snprintf(command, 1024, "rm -rf %s/*", P.LocalWriteStateDirectory);
        assert(printret >= 0 && printret < 1024);
        int ret = system(command);  // hacky!
    }

    // Delete the read state and move write to read
    if(!P.ProfilingMode)
        MoveLocalDirectories();

    FinalizeParallel();  // This may be the last synchronization point?
    SingleStepTearDown.Stop();

    // Finished cleanup.  Now we can log that time too.
    ReportTimings();

    stdlog.close();
    
    return 0;
}
