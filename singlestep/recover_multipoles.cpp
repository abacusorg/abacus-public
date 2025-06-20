// Copyright 2012-2025 The Abacus Developers
// SPDX-License-Identifier: GPL-3.0-or-later

/*  This is a multipole recovery utility that simply loads positions
 *  from the read state and outputs multipoles to the multipole directory.
 *  
 *  This executable is primarily invoked by `Abacus/abacus.py` when it
 *  detects missing multipoles.
 */
 
#include "proepi.cpp"
#include "timestep_recover_multipoles.cpp"

#ifdef PARALLEL
void RecoverReadStateFiles(Parameters &P){
	STDLOG(0, "RecoverReadStateFiles not yet implemented.");
}
#endif

 
 int main(int argc, char **argv) {
    WallClockDirect.Start();
    SingleStepSetup.Start();

	int num_args = 2;
	
    if (argc!=num_args) {
       // Can't use assertf() or QUIT here: stdlog not yet defined!
       fmt::print(stderr, "recover_multipoles: command line must have 2 parameters given, not {:d}.\nLegal usage: recover_multipoles <parameter_file> (and, if parallel: <reconstruct_read_files> <reconstruct_multipoles>)\n", argc);
       assert(0==99);
    }  

    InitializeParallel(MPI_size, MPI_rank);
	
    P.ReadParameters(argv[1],0);

    bool MakeIC;
    MakeIC = true;  // For most purposes, pretend we're making ICs

    load_read_state(0);  // Load the bare minimum (step num, etc) to get the log up and running

    setup_log(); // STDLOG and assertf now available
    STDLOG(0,"Read Parameter file {:s}\n", argv[1]);

    SetupLocalDirectories(0);

    init_openmp();

    double da;
    check_read_state(MakeIC, da);
    da = 0.; // never used, but let's be clear
    
    SingleStepSetup.Stop();

    // Informs some of our directory structure
    // Also sets up SlabSize
    InitWriteState("recover_multipoles", argv[1]);
	Prologue(P, MakeIC);
    SS->load_from_params(P);  // normally done during Prologue, but we pretended this was an IC step
	
	timestepMultipoles();
	
#ifdef PARALLEL	
	RecoverReadStateFiles(P); 
#endif
	
    // TODO: is the following parallel-unsafe for any reason?		
    // Let the IO finish, so that it is included in the time log.
    SingleStepTearDown.Start();
    IO_Terminate();
    SingleStepTearDown.Stop();
    WallClockDirect.Stop();

    // Stop Epilogue from trying to write slabsize
    delete SS;
    SS = NULL;

    // The epilogue contains some tests of success.
    Epilogue(P);
	
	stdlog.close();
    FinalizeParallel(); 
	
	exit(0); 
}
