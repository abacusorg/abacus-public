#include "proepi.cpp"

    
int main(int argc, char **argv) {
    WallClockDirect.Clear();
    WallClockDirect.Start();

    if(argc!=4) {
        printf("usage: singlestep <sampleparameterfile> kickfactor driftfactor\n");
        assert(1==0);
    }
    
    P.ReadParameters(argv[1],1);

    GlobalKickFactor    = atof(argv[2]);
    GlobalDriftFactor   = atof(argv[3]);

    Prologue(P,0);

    timestep(); 

    WallClockDirect.Stop();
    ReportTimings();

    Epilogue(P,0);
    exit(0);
}
