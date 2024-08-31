/*  This is an IO benchmark utility that simply reads the state files
 *  and writes them back out.
 */
 
#include "proepi.cpp"
#include "timestep_benchmark_io.cpp"
 
 int main(int argc, char **argv) {
    WallClockDirect.Start();
    SingleStepSetup.Start();

    if(argc < 2) {
       fmt::print(stderr, "{:s}: command line must have at least 1 parameter given.\nUsage: {:s} <parameter_file> [nslabs]\n", argv[0], argv[0]);
       exit(1);
    }
    
    P.ReadParameters(argv[1],0);
    strcpy(WriteState.ParameterFileName, argv[1]);

    setup_log(); // STDLOG and assertf now available
    STDLOG(0,"Read Parameter file {:s}\n", argv[1]);

    init_openmp();

    bool MakeIC;
    double da;
    check_read_state(0, MakeIC, da);
    MakeIC = true;  // For most purposes, pretend we're making ICs
    da = 0.; // never used, but let's be clear
    
    SingleStepSetup.Stop();

    // Now execute the timestep
    Prologue(P,MakeIC);
    SS = new SlabSize(P);
    int nslab = argc >= 3 ? atoi(argv[2]) : -1;
    timestepBenchmarkIO(nslab);

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