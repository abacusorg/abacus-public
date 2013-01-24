#include "proepi.cpp"

    
int main(int argc, char **argv) {
    WallClockDirect.Clear();
    WallClockDirect.Start();

    if(argc!=4) {
        printf("usage: singlestep <sampleparameterfile> kickfactor driftfactor\n");
        assert(1==0);
    }
    
    P.ReadParameters(argv[1],1);
    stdlog.open("mylog");   // Need a real name for this.
    STDLOG("Read Parameter file %s\n", argv[1]);

    GlobalKickFactor    = atof(argv[2]);
    GlobalDriftFactor   = atof(argv[3]);
    STDLOG("GlobalKickFactor = %f\n", GlobalKickFactor);
    STDLOG("GlobalDriftFactor = %f\n", GlobalDriftFactor);

    Prologue(P,0);

    timestep(); 

    WallClockDirect.Stop();
    ReportTimings();

    Epilogue(P,0);
    stdlog.close();  
    exit(0);
}
