#include "proepi.cpp"

state ReadState, WriteState;


Cosmology * InitializeCosmology() {
    MyCosmology cosmo;
    cosmo.Omega_m = P.Omega_M;
    cosmo.Omega_K = P.Omega_K;
    cosmo.Omega_DE = P.Omega_DE;
    cosmo.H0 = P.H0;
    cosmo.w0 = P.w0;
    cosmo.wa = P.wa;
    return new Cosmology(ReadState.ScaleFactor,cosmo);
}

double ChooseTimeStep(){
	return .01;
}

void BuildWriteState(double da){
	//get the next timestep and build the cosmology for it
	double nexta = cosm->current.a + da;
	cosm->BuildEpoch(cosm->current, cosm->next, nexta);

	//fill in WriteState

	WriteState.ScaleFactor = cosm->next.a;
	WriteState.Redshift = cosm->next.z;
	WriteState.Time = cosm->next.t;                // In Gyr or Gyr/h, depending on hMpc flag
	WriteState.etaK = cosm->next.etaK;
	WriteState.etaD = cosm->next.etaD;
	WriteState.Growth = cosm->next.growth;
	WriteState.Growth_on_a = cosm->next.growth/cosm->next.a;
	WriteState.f_growth = cosm->next.f_growth;
	WriteState.w = cosm->next.w;
	WriteState.HubbleNow = cosm->next.H;           // In km/s/Mpc
	WriteState.Htime = cosm->next.H*cosm->next.t;               // Time*H(z)

	double total = cosm->next.OmegaHat_m+cosm->next.OmegaHat_X+cosm->next.OmegaHat_K;
	WriteState.OmegaNow_m = cosm->next.OmegaHat_m/total;
	WriteState.OmegaNow_K = cosm->next.OmegaHat_K/total;
	WriteState.OmegaNow_DE = cosm->next.OmegaHat_DE/total;

	WriteState.ParticleMass = ReadState.ParticleMass; //FIXME: This is just a place holder // In Msun or Msun/h, depending on hMpc flag
	WriteState.RedshiftSpaceConversion =ReadState.RedshiftSpaceConversion;//FIXME: Another placeholder until the actual math is worked out
	WriteState.DeltaTime = cosm->next.t - cosm->current.t;
	WriteState.DeltaRedshift = -(cosm->next.z - cosm->current.z);
	WriteState.DeltaScaleFactor = cosm->next.a - cosm->current.a;

	// Might also compute the Scale Factor for the halfway time, since
	// we'll need that.
	cosm->t2a(0.5*(cosm->next.t+cosm->current.t));
	// cosm->search now has the midpoint epoch.
	WriteState.ScaleFactorHalf = cosm->search.a;
	WriteState.LastHalfEtaKick = cosm->KickFactor(cosm->search.a,WriteState.ScaleFactor-cosm->search.a);
	WriteState.FirstHalfEtaKick = cosm->KickFactor(cosm->current.a,cosm->search.a-cosm->current.a);

	WriteState.LPTstatus = ReadState.LPTstatus;
	WriteState.FullStepNumber = ReadState.FullStepNumber +1;
}


int main(int argc, char **argv) {
    WallClockDirect.Clear();
    WallClockDirect.Start();

    if(argc!=3) {
        printf("usage: singlestep <parameter file name> <allow creation of initial conditions 1/0>\n");
        assert(1==0);
    }
    
    int AllowIC = atoi(argv[2]);

    P.ReadParameters(argv[1],1);
    stdlog.open("mylog");   // Need a real name for this.
    STDLOG("Read Parameter file %s\n", argv[1]);


    //Check if ReadStateDirectory is accessible, or if we should build a new state from the IC file
    if(access(P.ReadStateDirectory,0) ==-1){
    	if(AllowIC != 1){
    		STDLOG("ERROR: Read State Directory ( %s ) is inaccessible and initial state creation is prohibited. Terminating.\n",P.ReadStateDirectory);
    		fprintf(stderr,"ERROR: Read State Directory ( %s ) is inaccessible and initial state creation is prohibited. Terminating.\n",P.ReadStateDirectory);
    	}
    	else{
    		//Make initial conditions
    	}
    }
    else{
    	CheckDirectoryExists(P.ReadStateDirectory);
    	readstate(ReadState,P.ReadStateDirectory);
    	STDLOG("Read ReadState from %s\n",P.ReadStateDirectory);
    }


    //Check if WriteStateDirectory exists, and fail if it does
    assert(access(P.WriteStateDirectory,0) ==-1);

    cosm = InitializeCosmology();
    STDLOG("Initialized Cosmology at a= %4.2f\n",ReadState.ScaleFactor);

    double da = ChooseTimeStep();
    STDLOG("Chose Time Step da = $5.4f\n",da);
    BuildWriteState(da);

    GlobalKickFactor    = cosm->KickFactor(ReadState.ScaleFactor,da);
    GlobalDriftFactor   = cosm->DriftFactor(ReadState.ScaleFactor,da);
    STDLOG("GlobalKickFactor = %f\n", GlobalKickFactor);
    STDLOG("GlobalDriftFactor = %f\n", GlobalDriftFactor);

    Prologue(P,0);

    timestep(); 

    WallClockDirect.Stop();
    ReportTimings();

    Epilogue(P,0);

    WriteState(WriteState,P.WriteStateDirectory);
    STDLOG("Wrote WriteState to %s\n",P.WriteStateDirectory);

    stdlog.close();  
    exit(0);
}
