/*  Run the on-the-fly group finder on a state or time slice outputs.
 *  This executable is primarily invoked via the `Abacus/standalone_fof.py` script.
 *  TODO: add functionality to run on states.
 */
 
 #include "proepi.cpp"
 
 int main(int argc, char **argv) {
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

    // Override any parameters that it don't make sense in this context
    P.AllowGroupFinding = 1;
    P.L1Output_dlna = 0.;  // force group finding output
    WriteState.ScaleFactor = ReadState.ScaleFactor;
    // TODO: is it okay to leave the WriteState uninitialized?

    setup_log(); // STDLOG and assertf now available
    STDLOG(0,"Beginning standalone_fof\n");
    STDLOG(0,"Read Parameter file %s\n", argv[2]);

    STDLOG(1,"Slab header indicates CPD %d, ppd %d\n", P.cpd, (int) ReadState.ppd);
    STDLOG(1,"Using group radius %d and MinL1HaloNP %d\n", P.GroupRadius, P.MinL1HaloNP);
    STDLOG(1,"Running with %d threads\n", omp_get_max_threads());

    init_openmp();

    bool MakeIC = true;  // For most purposes, pretend we're making ICs (e.g. no forces)
    
    SingleStepSetup.Stop();

    // Now execute the timestep
    Prologue(P,MakeIC);

    // There are some things we set up in the Prologue that would be used in an IC step, but not here
    delete MF;
    MF = NULL;

    // Prologue didn't do this because it thought it was an IC step
    GFC = new GroupFindingControl(P.FoFLinkingLength[0]/pow(P.np,1./3),
                                  P.FoFLinkingLength[1]/pow(P.np,1./3),
                                  P.FoFLinkingLength[2]/pow(P.np,1./3),
                                  P.cpd, PP->invcpd, P.GroupRadius, P.MinL1HaloNP, P.np);
    // Do we need a SlabSizePack14?
    //Slab = new SlabSize(P.cpd);
    //load_slabsize(P);
    timestepStandaloneFOF(argv[1]);

    // Let the IO finish, so that it is included in the time log.
    SingleStepTearDown.Start();
    IO_Terminate();
    SingleStepTearDown.Stop();
    WallClockDirect.Stop();

    // The epilogue contains some tests of success.
    Epilogue(P,MakeIC);
    stdlog.close();
    
    fprintf(stderr, "%s completed successfully!\n", argv[0]);
    exit(0);
}
