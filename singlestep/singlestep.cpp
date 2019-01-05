#include "proepi.cpp"

#include "cosmo_setup.cpp"

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
    #ifdef PARALLEL
        // MPI_Comm_rank(MPI_COMM_WORLD, &WriteState.NodeRank);
    #else 
        WriteState.NodeRank = 0;
    #endif

	WriteState.DoublePrecision = (sizeof(FLOAT)==8)?1:0;
	STDLOG(0,"Bytes per float is %d\n", sizeof(FLOAT));
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
}

void BuildWriteStateOutput() {
	// Build the output header.
	// Note we actually will output from ReadState,
	// but we build this to write the write/state file from the same code.
	WriteState.make_output_header();
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


void InitGroupFinding(int MakeIC){
    int do_output;
    // Request output of L1 groups and halo/field subsamples if:
        // - By going from ReadState to WriteState we are crossing a L1Output_dlna checkpoint
        // - We are doing a TimeSlice output
    // We may not end up outputting group if GFC is not initialized below
    if(P.L1Output_dlna >= 0)
        do_output = log(WriteState.ScaleFactor) - log(ReadState.ScaleFactor) >= P.L1Output_dlna ||
                    fmod(log(WriteState.ScaleFactor), P.L1Output_dlna) < fmod(log(ReadState.ScaleFactor), P.L1Output_dlna);
    else
        do_output = 0;
    do_output |= ReadState.DoTimeSliceOutput;

    WriteState.DensityKernelRad2 = 0.0;   // Don't compute densities

    // We need to enable group finding if:
        // - We are doing microstepping
        // - We are outputting groups
    // But we can't enable it if:
        // - AllowGroupFinding is disabled
        // - ForceOutputDebug is enabled
        // - This is an IC step
    // ForceOutputDebug outputs accelerations as soon as we compute them
    // i.e. before GroupFinding has a chance to rearrange them
    if((P.MicrostepTimeStep > 0 || do_output) &&
        !(!P.AllowGroupFinding || P.ForceOutputDebug || MakeIC)){
        STDLOG(1, "Setting up group finding\n");
        
        ReadState.DoGroupFindingOutput = do_output;

        GFC = new GroupFindingControl(P.FoFLinkingLength[0]/pow(P.np,1./3),
                    #ifdef SPHERICAL_OVERDENSITY
                    P.SODensity[0], P.SODensity[1],
                    #else
                    P.FoFLinkingLength[1]/pow(P.np,1./3),
                    P.FoFLinkingLength[2]/pow(P.np,1./3),
                    #endif
                    P.cpd, P.GroupRadius, P.MinL1HaloNP, P.np);

        #ifdef COMPUTE_FOF_DENSITY
        #ifdef CUDADIRECT   // For now, the CPU doesn't compute FOF densities, so signal this by leaving Rad2=0.
        WriteState.DensityKernelRad2 = GFC->linking_length;
        WriteState.DensityKernelRad2 *= WriteState.DensityKernelRad2*(1.0+1.0e-5); 
        // We use square radii.  The radius is padded just a little
        // bit so we don't risk underflow with 1 particle at r=b
        // in comparison to the self-count.
        #endif
        #endif
        STDLOG(1,"Using DensityKernelRad2 = %f (%f of interparticle)\n", WriteState.DensityKernelRad2, sqrt(WriteState.DensityKernelRad2)*pow(P.np,1./3.));

    } else{
        STDLOG(1, "Group finding not enabled for this step.\n");
    }
}

void InitializeParallel() {
    #ifdef PARALLEL
         // TODO: MPI_Init() and other items
         int size = 0, rank = 0;
         // MPI_Init(NULL, NULL);
         STDLOG(0,"Initializing MPI.");   
         // MPI_Comm_size(MPI_COMM_WORLD, &size);
         // MPI_Comm_rank(MPI_COMM_WORLD, &rank);
         STDLOG(0,"Node %d of %d total\n", rank, size);
         sprintf(NodeString,".%04d",rank);
    #else
    #endif
    return;
}

void FinalizeParallel() {
    #ifdef PARALLEL
         // TODO: MPI_Finalize() and other items
         // MPI_Finalize();
         STDLOG(0,"Calling MPI_Finalize()");
    #else
    #endif
}


int main(int argc, char **argv) {
    WallClockDirect.Start();
    SingleStepSetup.Start();

    if (argc!=3) {
       // Can't use assertf() or QUIT here: stdlog not yet defined!
       fprintf(stderr, "singlestep(): command line must have 3 parameters given, not %d.\nLegal usage: singlestep PARAMETER_FILE MAKE_IC\n\tPARAMETER_FILE: path to parameter file (usually called abacus.par)\n\tMAKE_IC: 0 or 1, whether this is an IC step.\n", argc);
       return 1;
    }
    InitializeParallel();
    
    int MakeIC = atoi(argv[2]);
    P.ReadParameters(argv[1],0);
    strcpy(WriteState.ParameterFileName, argv[1]);

    setup_log(); // STDLOG and assertf now available
    STDLOG(0,"Read Parameter file %s\n", argv[1]);
    STDLOG(0,"MakeIC = %d\n", MakeIC);

    SetupStateDirectories(P, MakeIC);
    
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
    InitWriteState(MakeIC);

    // Set up the major classes (including NFD)
    Prologue(P,MakeIC);

    // Check if WriteStateDirectory/state exists, and fail if it does
    char wstatefn[1050];
    sprintf(wstatefn,"%s/state", P.WriteStateDirectory);
    if(access(wstatefn,0) !=-1 && !WriteState.OverwriteState)
        QUIT("WriteState \"%s\" exists and would be overwritten. Please move or delete it to continue.\n", wstatefn);
    
    if (da!=0) da = ChooseTimeStep();
    double dlna = da/ReadState.ScaleFactor;
    STDLOG(0,"Chose Time Step da = %6.4f, dlna = %6.4f\n", da, dlna);
    if(dlna > 0){
        STDLOG(0, "\t\tAt the current rate, this implies %d more steps to z_final=%f\n", (int)ceil(log(1./ReadState.ScaleFactor/(1. + P.FinishingRedshift()))/dlna), P.FinishingRedshift());
    }
    
    //Enable floating point exceptions unless we are doing profiling (where they tend to break the profilier)
    feenableexcept(FE_INVALID | FE_DIVBYZERO);
    //if(!P.ProfilingMode) feenableexcept(FE_INVALID | FE_DIVBYZERO);
    //else fedisableexcept(FE_INVALID | FE_DIVBYZERO);
    
    BuildWriteState(da);
    LCOrigin = (double3 *) malloc(8*sizeof(double3));  // max 8 light cones
    for(int i = 0; i < 8; i++)
        LCOrigin[i] = ((double3*) P.LightConeOrigins)[i]/P.BoxSize;  // convert to unit-box units
    
    // Make a plan for output
    PlanOutput(MakeIC);

    // Set up the Group Finding concepts
    InitGroupFinding(MakeIC);
    BuildWriteStateOutput();    // Have to delay this until after GFC is made

    SingleStepSetup.Stop();

    // Now execute the timestep
    if (MakeIC)  timestepIC();
	    else timestep();

    // Let the IO finish, so that it is included in the time log.
    SingleStepTearDown.Start();
    IO_Terminate();
    fedisableexcept(FE_INVALID | FE_DIVBYZERO);

    SingleStepTearDown.Stop();
    WallClockDirect.Stop();

    // The epilogue contains some tests of success.
    Epilogue(P,MakeIC);
    delete cosm;
    free(LCOrigin);

    // Print out some final stats
    FinalizeWriteState();

    // The state should be written last, since that officially signals success.
    // TODO: For MPI, this should be the Global Write State.
    // TODO: For MPI, only rank 0 node does this
    if (WriteState.NodeRank==0) {
        WriteState.write_to_file(P.WriteStateDirectory);
        STDLOG(0,"Wrote WriteState to %s\n",P.WriteStateDirectory);
    }
    
    if (!MakeIC && P.ProfilingMode){
        STDLOG(0,"ProfilingMode is active. Removing the write state in %s\n",P.LocalWriteStateDirectory);
        char command[1024];
        sprintf(command, "rm -rf %s/*", P.LocalWriteStateDirectory);
        int ret = system(command);  // hacky!
    }
    stdlog.close();

    FinalizeParallel();  // This may be the last synchronization point?
    exit(0);
}
