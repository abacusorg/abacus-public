#include "proepi.cpp"
#include "timestep_ic.cpp"
#include "cosmo_setup.cpp"


void BuildWriteState(double da){
	STDLOG(0,"Building WriteState for a step from a={:f} by da={:f}\n", cosm->current.a, da);

	// fill in WriteState from the Parameter file
	WriteState.np_state = P.np;
	WriteState.cpd_state = P.cpd;
	WriteState.order_state = P.order;
	WriteState.ppd = P.ppd();
    WriteState.ippd = (int64) round(WriteState.ppd);
    // DoTimeSliceOutput, DoSubsampleOutput set in ChooseTimeStep()
    WriteState.VelIsSynchronous = da == 0;
        
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
    WriteState.NodeRankX = MPI_rank_x;
    WriteState.NodeRankZ = MPI_rank_z;
    WriteState.NodeSizeX = MPI_size_x;
    WriteState.NodeSizeZ = MPI_size_z;

	WriteState.DoublePrecision = (sizeof(FLOAT)==8)?1:0;
	STDLOG(0,"Bytes per float is {:d}\n", sizeof(FLOAT));
	STDLOG(0,"Bytes per auxstruct is {:d}\n", sizeof(auxstruct));
	STDLOG(0,"Bytes per cellinfo is {:d}\n", sizeof(cellinfo));
	assert(WriteState.FullStepNumber == ReadState.FullStepNumber+1);  // already set this in load_read_state()
	STDLOG(0,"This is step number {:d}\n", WriteState.FullStepNumber);

	//get the next timestep and build the cosmology for it
	double nexta = cosm->current.a + da;
	STDLOG(0,"Next scale factor is {:f}\n", nexta);
	cosm->BuildEpoch(cosm->current, cosm->next, nexta);
	FillStateWithCosmology(WriteState);

	// Get differences between the two cosmological epochs
	WriteState.DeltaTime = cosm->next.t - cosm->current.t;
	WriteState.DeltaRedshift = -(cosm->next.z - cosm->current.z);
	WriteState.DeltaScaleFactor = cosm->next.a - cosm->current.a;

	cosm->t2a(0.5*(cosm->next.t+cosm->current.t));
	STDLOG(0,"Scale factor halfway in between is {:f}\n", cosm->search.a);
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

    // carry forward
    WriteState.LPTVelScale = ReadState.LPTVelScale;

    // Set the WriteState softening, will be used by the next time step

    // Decrease the softening length if we are doing a 2LPT step
    // This helps ensure that we are using the true 1/r^2 force
    /*if(LPTStepNumber()>0){
        WriteState.SofteningLengthNow = P.SofteningLength / 1e4;  // This might not be in the growing mode for this choice of softening, though
        STDLOG(0,"Reducing softening length from {:f} to {:f} because this is a 2LPT step.\n", P.SofteningLength, WriteState.SofteningLengthNow);
    }
    else{
        WriteState.SofteningLengthNow = P.SofteningLength;
    }*/

    // Is the softening fixed in proper coordinates?
    if(P.ProperSoftening){
        WriteState.SofteningLengthNow = min(P.SofteningLength/WriteState.ScaleFactor, P.SofteningMax);
        STDLOG(1, "Adopting a comoving softening of {:g}, fixed in proper coordinates\n", WriteState.SofteningLengthNow);
    }
    else{
        WriteState.SofteningLengthNow = P.SofteningLength;
        STDLOG(1, "Adopting a comoving softening of {:g}, fixed in comoving coordinates\n", WriteState.SofteningLengthNow);
    }

    // Now scale the softening to match the minimum Plummer orbital period
#if defined DIRECTCUBICSPLINE
    strcpy(WriteState.SofteningType, "cubic_spline");
    WriteState.SofteningLengthNowInternal = WriteState.SofteningLengthNow * 1.10064;
#elif defined DIRECTSINGLESPLINE
    strcpy(WriteState.SofteningType, "single_spline");
    WriteState.SofteningLengthNowInternal = WriteState.SofteningLengthNow * 2.15517;
#elif defined DIRECTCUBICPLUMMER
    strcpy(WriteState.SofteningType, "cubic_plummer");
    WriteState.SofteningLengthNowInternal = WriteState.SofteningLengthNow * 1.;
#else
    strcpy(WriteState.SofteningType, "plummer");
    WriteState.SofteningLengthNowInternal = WriteState.SofteningLengthNow;
#endif
}

void BuildOutputHeaders() {
	// Build the output headers.
	// Note that outputs at this epoch use the ReadState header
    ReadState.make_output_header();
	WriteState.make_output_header();
}

void PlanOutput(bool MakeIC) {
    // Set up some metadata for time slice output. The decision is made in ChooseTimeStep().
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
    ReadState.ParameterFileName = WriteState.ParameterFileName;
    ReadState.CodeVersion = WriteState.CodeVersion;
    ReadState.MachineName = WriteState.MachineName;
    ReadState.RunTime = WriteState.RunTime;
    ReadState.FullStepNumber = WriteState.FullStepNumber;
    
    ReadState.NodeRankX = WriteState.NodeRankX;
    ReadState.NodeRankZ = WriteState.NodeRankZ;
    ReadState.NodeSizeX = WriteState.NodeSizeX;
    ReadState.NodeSizeZ = WriteState.NodeSizeZ;

    // Just let later routines know that this is a valid epoch
    // for output, e.g., not a LPT IC epoch.
    ReadState.OutputIsAllowed = 1;
    if(P.OutputEveryStep) ReadState.DoTimeSliceOutput = 1;

    // Now check whether we're asked to do a TimeSlice.
	if (ReadState.DoTimeSliceOutput) {
	    STDLOG(0,"Planning to output a TimeSlice\n");
        fs::create_directory(P.OutputDirectory / fmt::format("slice{:5.3f}", ReadState.Redshift));
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

    WallClockDirect.Start();
    SingleStepSetup.Start();

    if (argc!=3) {
       // Can't use assertf() or QUIT here: stdlog not yet defined!
       fmt::print(stderr, "singlestep(): command line must have 3 parameters given, not {:d}.\nLegal usage: singlestep PARAMETER_FILE MAKE_IC\n\tPARAMETER_FILE: path to parameter file (usually called abacus.par)\n\tMAKE_IC: 0 or 1, whether this is an IC step.\n", argc);
       return 1;
    }

    // Set up MPI
    StartMPI();  // call MPI_Init as early as possible

    //Enable floating point exceptions
    feenableexcept(FE_INVALID | FE_DIVBYZERO);
    
    P.ReadParameters(argv[1],0);

    InitializeParallelTopology();  // MPI_rank now available
    P.ProcessStateDirectories();  // needs MPI_rank

    int MakeIC = atoi(argv[2]);
    load_read_state(MakeIC);  // Load the bare minimum (step num, etc) to get the log up and running
    int Do2DGroupFinding = MPI_size_z > 1 && ReadState.DoGroupFindingOutput && ReadState.HaveAuxDensity && ReadState.VelIsSynchronous;
    int NoForces = MakeIC || Do2DGroupFinding;

    setup_log(); // STDLOG and assertf now available
    STDLOG(0,"Read Parameter file {:s}\n", argv[1]);
    STDLOG(0,"MakeIC = {:d}\n", MakeIC);
    STDLOG(0,"NoForces = {:d}\n", NoForces);

    InitializeForceRadius(NoForces);
    InitializeParallelDomain();  // needs ReadState
    
    SetupLocalDirectories(MakeIC);

    // Set up OpenMP
    init_openmp();

    check_read_state(MakeIC);

    // Initialize the Cosmology and set up the State epochs and the time step
    cosm = InitializeCosmology(ReadState.ScaleFactor);
    if (MakeIC) FillStateWithCosmology(ReadState);

    // Set some WriteState values before ChooseTimeStep()
    // This also sets the SofteningLength, needed by the NFD constructor
    InitWriteState(MakeIC, "singlestep", argv[1]);

    // Set up the major classes (including NFD)
    Prologue(P,MakeIC,NoForces);

    double da = ChooseTimeStep(NoForces);

    // da *= -1;  // reverse the time step TODO: make parameter
    double dlna = da/ReadState.ScaleFactor;
    STDLOG(0,"Chose Time Step da = {:6.4f}, dlna = {:6.4f}\n", da, dlna);
    if(dlna > 0){
        STDLOG(0, "\t\tAt the current rate, this implies {:d} more steps to z_final={:f}\n", (int64)ceil(log(1./ReadState.ScaleFactor/(1. + P.FinishingRedshift()))/dlna), P.FinishingRedshift());
    }

    BuildWriteState(da);
    InitializeLightCones();

    // Set up the Group Finding concepts and decide if Group Finding output is requested.
    InitGroupFinding(MakeIC);

    // Set up output metadata
    PlanOutput(MakeIC);
    
    InitializePipelineWidths(MakeIC);  // needs to know if this step will do group finding
    InitializeParallelMergeDomain();  // needs to know if the *next* step will do group finding
    LogParallelTopology();

    BuildOutputHeaders();    // needs group finding and merge domain

    SingleStepSetup.Stop();

    // Now execute the timestep
    if (MakeIC)  timestepIC();
	    else timestep(NoForces);

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
    FinalizeLightCones();
    free_dependencies();

    // Print out some final stats
    FinalizeWriteState();

    // The state should be written last, since that officially signals success.
    if (MPI_rank==0) {
        WriteState.write_to_file(P.WriteStateDirectory);
        STDLOG(0,"Wrote WriteState to {}\n",P.WriteStateDirectory);
    }

    if (!MakeIC && P.ProfilingMode){
        STDLOG(0,"ProfilingMode is active. Removing the write state in {}\n",P.LocalWriteStateDirectory);
        fs::remove_all(P.LocalWriteStateDirectory);
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
