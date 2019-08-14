/*  This is a multipole recovery utility that simply loads positions
 *  from the read state and outputs multipoles to the multipole directory.
 *  
 *  This executable is primarily invoked by `Abacus/abacus.py` when it
 *  detects missing multipoles.
 */
 
 #include "proepi.cpp"

#ifdef PARALLEL
void RecoverReadStateFiles(Parameters &P){
	STDLOG(0, "Beginning read state file recovery:\n");
	
	STDLOG(1, "\t nodeslabs\n");
	
	STDLOG(1, "\t state\n");
	
	STDLOG(1, "\t slabsize\n");
	
	STDLOG(0, "Completed read state file recovery.\n");
	
}
#endif

 
 int main(int argc, char **argv) {
    WallClockDirect.Start();
    SingleStepSetup.Start();

#ifdef PARALLEL 
	int num_args = 4; 
#else
	int num_args = 2;
#endif 
	
    if (argc!=num_args) {
       // Can't use assertf() or QUIT here: stdlog not yet defined!
       fprintf(stderr, "recover_multipoles: command line must have 2 (serial) or 4 (parallel) parameters given, not %d.\nLegal usage: recover_multipoles <parameter_file> (and, if parallel: <reconstruct_read_files> <reconstruct_multipoles>)\n", argc);
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
	
#ifdef PARALLEL
	int reconstruct_read_files = atoi(argv[3]); 
	int reconstruct_multipoles = atoi(argv[4]);

	Prologue(P, MakeIC, reconstruct_read_files);
	STDLOG(1, "Calling timestep for multipole and/or read state file recovery. reconstruct_multipoles = %d.\n", reconstruct_multipoles);
	
	timestepMultipoles(reconstruct_multipoles);
	if (reconstruct_read_files) RecoverReadStateFiles(P); 		
	
    // The epilogue contains some tests of success.
    Epilogue(P, MakeIC, reconstruct_read_files);
	
	stdlog.close();
	
    FinalizeParallel();  // This may be the last synchronization point?
	
	exit(0); 
	}
#else

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
	
#endif
}
