/*  Run the on-the-fly group finder on a state or time slice outputs.
 *  This executable is primarily invoked via the `Abacus/standalone_fof.py` script.
 *  TODO: add functionality to run on states.
 */
 
#include "proepi.cpp"
#include "timestep_standalone_fof.cpp"
 
 int main(int argc, char **argv) {
     //Enable floating point exceptions
    feenableexcept(FE_INVALID | FE_DIVBYZERO);
     
    WallClockDirect.Start();
    SingleStepSetup.Start();

    if(argc < 3) {
       fprintf(stderr, "%s: command line must have at least 2 parameters given.\nUsage: %s <time slice directory> <parameter file>\n", argv[0], argv[0]);
       exit(1);
    }

    // We will cheat and load this as both the parameters and state file
    HeaderStream hs(argv[2]);
    ReadState.ReadHeader(hs);
    hs.Close();
    P.ReadParameters(argv[2], 0);
    strcpy(ReadState.Pipeline, "standalone_fof");
    strcpy(ReadState.ParameterFileName, argv[1]);

    // Override any parameters that don't make sense in this context
    P.AllowGroupFinding = 1;
    P.L1Output_dlna = 0.;  // force group finding output
    WriteState.ScaleFactor = ReadState.ScaleFactor;
    // TODO: is it okay to leave the WriteState uninitialized?

    setup_log(); // STDLOG and assertf now available
    STDLOG(0,"Beginning standalone_fof\n");
    STDLOG(0,"Read Parameter file {:s}\n", argv[2]);
    char hostname[1024];
    gethostname(hostname,1024);
    STDLOG(0,"Host machine name is {:s}\n", hostname);

    STDLOG(1,"Slab header indicates CPD {:d}, ppd {:d}\n", P.cpd, (int) ReadState.ppd);
    STDLOG(1,"Using group radius {:d} and MinL1HaloNP {:d}\n", P.GroupRadius, P.MinL1HaloNP);
    STDLOG(1,"Running with {:d} threads\n", omp_get_max_threads());
     
    bool MakeIC = false;
     
    SetupLocalDirectories(MakeIC);

    init_openmp();
    
    double da = 0.;  // If we set this to zero, it will skip the timestep choice

    check_read_state(MakeIC, da);

    // Initialize the Cosmology and set up the State epochs and the time step
    cosm = InitializeCosmology(ReadState.ScaleFactor);
    // Set some WriteState values before ChooseTimeStep()
    // This also sets the SofteningLength, needed by the NFD constructor
    InitWriteState(MakeIC, "singlestep", argv[1]);
    
    SingleStepSetup.Stop();

    // Now execute the timestep
    Prologue(P,MakeIC);

    // There are some things we set up in the Prologue that would be used in an IC step, but not here
    delete MF;
    MF = NULL;
     
    //BuildWriteState(da);
     
     // Make a plan for output
    PlanOutput(MakeIC);

    // Set up the Group Finding concepts and decide if Group Finding output is requested.
    InitGroupFinding(MakeIC);
    BuildWriteStateOutput();    // Have to delay this until after GFC is made

    // Do we need a SlabSizePack14?  Right now we populate SS on-the-fly
    //SS = new SlabSize(P);
    timestepStandaloneFOF(argv[1]);

    // Let the IO finish, so that it is included in the time log.
    SingleStepTearDown.Start();
    IO_Terminate();
    SingleStepTearDown.Stop();
    WallClockDirect.Stop();
     
    // While we still have global objects, log their timings
    GatherTimings();

    // Stop Epilogue from trying to write slabsize
    delete SS;
    SS = NULL;

    // The epilogue contains some tests of success.
    Epilogue(P,MakeIC);
     
     // Print out some final stats
    FinalizeWriteState();
     
    // Finished cleanup.  Now we can log that time too.
    ReportTimings();
     
    stdlog.close();
    
    fprintf(stderr, "%s completed successfully!\n", argv[0]);
    exit(0);
}
