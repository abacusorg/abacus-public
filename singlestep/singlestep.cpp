#include "proepi.cpp"




Cosmology * InitializeCosmology(double ScaleFactor) {
    MyCosmology cosmo;
    cosmo.Omega_m = P.Omega_M;
    cosmo.Omega_K = P.Omega_K;
    cosmo.Omega_DE = P.Omega_DE;
    cosmo.H0 = P.H0;
    cosmo.w0 = P.w0;
    cosmo.wa = P.wa;
    return new Cosmology(ScaleFactor,cosmo);
}

double ChooseTimeStep(){
	return .01;
}

void BuildWriteState(double da){

	//make the WriteState directory
	char cmd[2048];
	sprintf(cmd,"mkdir %s",P.WriteStateDirectory);
	system(cmd);
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

	WriteState.ParticleMass = 1.0/P.np; //FIXME: This is just a place holder // In Msun or Msun/h, depending on hMpc flag
	WriteState.RedshiftSpaceConversion = 1.0 ;//FIXME: Another placeholder until the actual math is worked out
	WriteState.LPTstatus = 0; //TODO: Should this come from parameters?
	WriteState.FullStepNumber = 0;


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
    stdlog.open("mylog");   // TODO:Need a real name for this.
    STDLOG("Read Parameter file %s\n", argv[1]);
    stdlog.flush();

    double a;
    double da;
    bool MakeIC; //True if we should make the initial state instead of doing a real timestep
    //Check if ReadStateDirectory is accessible, or if we should build a new state from the IC file
    if(access(P.ReadStateDirectory,0) ==-1){
    	if(AllowIC != 1){
    		fprintf(stderr,"ERROR: Read State Directory ( %s ) is inaccessible and initial state creation is prohibited. Terminating.\n",P.ReadStateDirectory);
    		QUIT("Read State Directory ( %s ) is inaccessible and initial state creation is prohibited. Terminating.\n",P.ReadStateDirectory);

    	}
    	else{
    		STDLOG("ReadState not found... Generating initial State from initial conditions\n");
    		ReadState.ParticleMass = 1.0/P.np; //FIXME: This is just a place holder // In Msun or Msun/h, depending on hMpc flag
    		ReadState.RedshiftSpaceConversion = 1.0 ;//FIXME: Another placeholder until the actual math is worked out
    		ReadState.LPTstatus = P.LagrangianPTOrder;
    		ReadState.FullStepNumber = -1;
        	a = 1.0/(1+P.InitialRedshift);
        	da = 0;
        	MakeIC = true;
    	}
    }
    else{
    	CheckDirectoryExists(P.ReadStateDirectory);
    	readstate(ReadState,P.ReadStateDirectory);
    	STDLOG("Read ReadState from %s\n",P.ReadStateDirectory);
    	a = ReadState.ScaleFactor;
    	da = ChooseTimeStep();
    	MakeIC = false;
    }


    cosm = InitializeCosmology(a);
    STDLOG("Initialized Cosmology at a= %4.2f\n",a);

    //Check if WriteStateDirectory exists, and fail if it does
    if(access(P.WriteStateDirectory,0) !=-1) QUIT("WriteState exists and would be overwritten. Please move or delete it to continue\n");


    STDLOG("Chose Time Step da = %5.4f\n",da);
    BuildWriteState(da);

    Prologue(P,0);

    if (MakeIC) timestepIC();
    else timestep();

    WallClockDirect.Stop();
    //ReportTimings();

    Epilogue(P,0);

    writestate(&WriteState,P.WriteStateDirectory);
    STDLOG("Wrote WriteState to %s\n",P.WriteStateDirectory);

    stdlog.close();  
    exit(0);
}
