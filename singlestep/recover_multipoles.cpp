/*  This is a multipole recovery utility that simply loads positions
 *  from the read state and outputs multipoles to the multipole directory.
 *  
 *  This executable is primarily invoked by `Abacus/abacus.py` when it
 *  detects missing multipoles.
 */
 
 #include "proepi.cpp"
 
 int main(int argc, char **argv) {
    WallClockDirect.Start();
    SingleStepSetup.Start();

    if (argc!=2) {
       // Can't use assertf() or QUIT here: stdlog not yet defined!
       fprintf(stderr, "recover_multipoles: command line must have 2 parameters given, not %d.\nLegal usage: recover_multipoles <parameter_file>\n", argc);
       assert(0==99);
    }
    
    P.ReadParameters(argv[1],0);
    strcpy(WriteState.ParameterFileName, argv[1]);
    strcpy(WriteState.Pipeline, "recover_multipoles");

    setup_log(); // STDLOG and assertf now available
    STDLOG(0,"Read Parameter file %s\n", argv[1]);

    init_openmp();

    bool MakeIC;
    double da;
    check_read_state(MakeIC, da);
    MakeIC = true;  // For most purposes, pretend we're making ICs
    da = 0.; // never used, but let's be clear
    
    SingleStepSetup.Stop();

    // Informs some of our directory structure
    // Also sets up SlabSize
    InitWriteState(MakeIC);

    // Now execute the timestep
    Prologue(P,MakeIC);
    SS->load_from_params(P);  // normally done during Prologue, but we pretended this was an IC step
    timestepMultipoles();

    // Let the IO finish, so that it is included in the time log.
    SingleStepTearDown.Start();
    IO_Terminate();
    SingleStepTearDown.Stop();
    WallClockDirect.Stop();

    // The epilogue contains some tests of success.
    Epilogue(P,MakeIC);
    stdlog.close();
        
    exit(0);
}
