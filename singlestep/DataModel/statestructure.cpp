/* statestructure.cpp

This defines the State class.

States contain the time-dependent global information about the
simulation.  The code is evolving from the ReadState at time j to
the WriteState at time j+1.  The state class is written (in ASCII)
to the state file in the state subdirectory, so as to track the 
time-varying aspects of the simulation.

Time-independent instructions for the simulation should be in the
Parameter class.

State uses ParseHeader to handle the reading of the ASCII file.
However, it must provide its own routine to write the new file!

*/

#ifndef STATESTRUCTURE_CPP
#define STATESTRUCTURE_CPP
#include "ParseHeader.hh"
#include "Parameters.cpp"
#include <cassert>
#include <ctime>
#include <unistd.h>
#include <iostream>
#include <fstream>

class State:public ParseHeader{
public:

    long long int np_state;
    int cpd_state;
    int order_state;
    
    char ParameterFileName[1024];   // State must contain a pointer to the Parameter file
    char CodeVersion[1024];
    char RunTime[1024];
    char MachineName[1024];
    double ppd;		// Particles per dimension
    int DoublePrecision;  // =1 if code is using double precision positions
    char SofteningType[128];
    double SofteningLength;  // Effective Plummer length, used for timestepping.  Same units as BoxSize.
    double SofteningLengthInternal;  // The equivalent length for the current softening technique.  Same units as BoxSize.

    double ScaleFactor;
    int FullStepNumber; // Counting the full steps
    int LPTStepNumber; // Counting the full steps

    double BoxSizeMpc; 		// In Mpc
    double BoxSizeHMpc;		// In h^-1 Mpc
    double HubbleTimeGyr;      // In Gyr 
    double HubbleTimeHGyr;      // In h^-1 Gyr 
    double ParticleMassMsun;    // In Msun 
    double ParticleMassHMsun;    // In Msun/h

    double VelZSpace_to_kms; // The number of km/s across the full box.
    double VelZSpace_to_Canonical;  // Converting to code canonical units 
    	// from redshift-space displacements in unit box units.

    // =1..3 if we're deriving LPT initial conditions
    double LastHalfEtaKick; // What one should kick by to get the velocities to synchronous
    double ScaleFactorHalf;
    double FirstHalfEtaKick;
    double DeltaEtaDrift;

     // Derived quantities for the cosmology; these are all in H_0=1 units
    double Redshift;
    double Time;        // In H_0=1 units
    double etaK;
    double etaD;
    double Growth;
    double Growth_on_a;
    double f_growth;
    double w;
    double HubbleNow;   // In H_0=1 units
    double Htime;       // Time*H(z)
    double OmegaNow_m;
    double OmegaNow_K;
    double OmegaNow_DE;

    // More description of the last time step
    double DeltaTime;
    double DeltaScaleFactor;
    double DeltaRedshift;

    // Some statistics about the particle distribution that might be useful.
    double MaxVelocity;
    double MaxAcceleration;
    double MinVrmsOnAmax;
    int MaxCellSize;
    int MinCellSize;
    double StdDevCellSize;
    double RMS_Velocity;

    // The variables below are not intended for output.  Just used in the code.

    // We will write a lot of state information into the header of output
    // files.  We will collect that information here.
    std::string output_header;
    void make_output_header();
    const char *header() { 
    	return output_header.c_str(); 	// Standard C-style char[] string
    }
    int DoTimeSliceOutput;
    int OutputIsAllowed;
    int DoBinning;
    
    int Do2LPTVelocityRereading;

    void read_from_file(const char *fn);
    void write_to_file(const char *dir);
    
    State();
    
    void AssertStateLegal(Parameters &P);

};

    State::State(){
    	installscalar("np_state",np_state,MUST_DEFINE);
    	installscalar("cpd_state",cpd_state,MUST_DEFINE);
    	installscalar("order_state",order_state,MUST_DEFINE);
    	installscalar("ParameterFileName",ParameterFileName,DONT_CARE);
    	installscalar("ppd",ppd,DONT_CARE);
        installscalar("SofteningType", SofteningType,DONT_CARE);
        installscalar("SofteningLength", SofteningLength,DONT_CARE);
        installscalar("SofteningLengthInternal", SofteningLengthInternal,DONT_CARE);

    	sprintf(CodeVersion,"version_not_defined");
    	installscalar("CodeVersion",CodeVersion,DONT_CARE);
	// These will now be set in BuildWriteState();
	// Don't bother loading these in ReadState
    	// time_t timet = time(0);
    	// string now = string(asctime(localtime(&timet)));
    	// sprintf(RunTime,"%s",now.substr(0,now.length()-1).c_str()); //valid
    	// gethostname(MachineName,1024); //valid
    	installscalar("DoublePrecision",DoublePrecision, DONT_CARE);

	installscalar("ScaleFactor",ScaleFactor, MUST_DEFINE);
	installscalar("FullStepNumber",FullStepNumber,MUST_DEFINE);
    installscalar("LPTStepNumber",LPTStepNumber,DONT_CARE);
	installscalar("BoxSizeMpc",BoxSizeMpc,DONT_CARE);
	installscalar("BoxSizeHMpc",BoxSizeHMpc,DONT_CARE);
	installscalar("HubbleTimeGyr",HubbleTimeGyr,DONT_CARE);
	installscalar("HubbleTimeHGyr",HubbleTimeHGyr,DONT_CARE);
	installscalar("ParticleMassHMsun",ParticleMassHMsun,DONT_CARE);
	installscalar("ParticleMassMsun",ParticleMassMsun,DONT_CARE);
	installscalar("VelZSpace_to_kms",VelZSpace_to_kms,DONT_CARE);
	installscalar("VelZSpace_to_Canonical",VelZSpace_to_Canonical,DONT_CARE);
	installscalar("LastHalfEtaKick",LastHalfEtaKick,MUST_DEFINE);
    	installscalar("ScaleFactorHalf",ScaleFactorHalf,DONT_CARE);
    	installscalar("FirstHalfEtaKick",FirstHalfEtaKick,DONT_CARE);
    	installscalar("DeltaEtaDrift",DeltaEtaDrift,DONT_CARE);
    	installscalar("Redshift",Redshift,DONT_CARE);
	installscalar("Time",Time,MUST_DEFINE);        // In Gyr or Gyr/h, depending on hMpc flag
	installscalar("etaK",etaK,MUST_DEFINE);
    	installscalar("etaD",etaD,DONT_CARE);
    	installscalar("Growth",Growth,DONT_CARE);
    	installscalar("Growth_on_a",Growth_on_a,DONT_CARE);
    	installscalar("f_growth",f_growth,DONT_CARE);
    	installscalar("w",w,DONT_CARE);
    	installscalar("HubbleNow",HubbleNow,DONT_CARE);   // In km/s/Mpc
    	installscalar("Htime",Htime,DONT_CARE);       // Time*H(z)
    	installscalar("OmegaNow_m",OmegaNow_m,DONT_CARE);
    	installscalar("OmegaNow_K", OmegaNow_K,DONT_CARE);
    	installscalar("OmegaNow_DE",OmegaNow_DE,DONT_CARE);
    	installscalar("DeltaTime",DeltaTime,DONT_CARE);
    	installscalar("DeltaScaleFactor",DeltaScaleFactor,DONT_CARE);
    	installscalar("DeltaRedshift",DeltaRedshift,DONT_CARE);
	
	MaxVelocity = 0;
    	installscalar("MaxVelocity",MaxVelocity,DONT_CARE);
	MaxAcceleration = 0;
    	installscalar("MaxAcceleration",MaxAcceleration,DONT_CARE);
        MinVrmsOnAmax = 1e10;
    	installscalar("MinVrmsOnAmax",MinVrmsOnAmax,DONT_CARE);
	MaxCellSize = 0;
    	installscalar("MaxCellSize",MaxCellSize,DONT_CARE);
	MinCellSize = 1e9;
    	installscalar("MinCellSize",MinCellSize,DONT_CARE);
	StdDevCellSize = 0.0;
    	installscalar("StdDevCellSize",StdDevCellSize,DONT_CARE);
	RMS_Velocity = 0.0;
    	installscalar("RMS_Velocity",RMS_Velocity,DONT_CARE);

	// Initialize helper variables
	DoTimeSliceOutput = 0;
     	OutputIsAllowed = 0;
        
    Do2LPTVelocityRereading = 0;
    }


void State::read_from_file(const char *fn) {
    char statefn[1050];
    sprintf(statefn,"%s/state",fn);
    HeaderStream hs(statefn);
    ReadHeader(hs);
    hs.Close();
}



#define FSYM 20.15f
#define ESYM 25.15e
#define ISYM d
#define SSYM s


#define PRQUOTEME(X) #X
#define WPR(X,XSYM) sprintf(tmp, PRQUOTEME(%22s = %XSYM\n), PRQUOTEME(X), X); ss << tmp
#define WPRS(X,XSYM) sprintf(tmp, "%22s = \"%s\" \n", PRQUOTEME(X), X); ss << tmp

void State::make_output_header() {
    // We're going to output most, but not all of the fields, into a 
    // nice header.  This will be prepended to many output files.
    // It also will get used to write the state.
    // Only those items that are setup in BuildWriteState should be in here.
    char tmp[1024];
    std::stringstream ss;

    WPRS(ParameterFileName        , s);
    WPRS(CodeVersion              , s);
    WPRS(RunTime                  , s);
    WPRS(MachineName              , s);
    WPR(DoublePrecision          , ISYM);
    WPR(ppd                      , FSYM);

    WPR(FullStepNumber           , ISYM);
    WPR(LPTStepNumber            , ISYM);
    WPR(ScaleFactor              , ESYM);
    WPR(BoxSizeMpc               , FSYM);
    WPR(BoxSizeHMpc              , FSYM);
    WPR(HubbleTimeGyr            , FSYM);
    WPR(HubbleTimeHGyr           , FSYM);
    WPR(ParticleMassMsun         , ESYM);
    WPR(ParticleMassHMsun        , ESYM);
    WPR(VelZSpace_to_kms         , ESYM);
    WPR(VelZSpace_to_Canonical   , FSYM);

    WPR(Redshift                 , FSYM);
    WPR(Time                     , FSYM);
    WPR(etaK                     , FSYM);
    WPR(etaD                     , FSYM);
    WPR(Growth                   , FSYM);
    WPR(Growth_on_a              , FSYM);
    WPR(f_growth                 , FSYM);
    WPR(w                        , FSYM);
    WPR(HubbleNow                , FSYM);
    WPR(Htime                    , FSYM);
    WPR(OmegaNow_m               , FSYM);
    WPR(OmegaNow_K               , FSYM);
    WPR(OmegaNow_DE              , FSYM);
    
    WPRS(SofteningType            , s);
    WPR(SofteningLength          , ESYM);
    WPR(SofteningLengthInternal  , ESYM);

    WPR(DeltaTime                , FSYM);
    WPR(DeltaScaleFactor         , FSYM);
    WPR(DeltaRedshift            , FSYM);

    WPR(DeltaEtaDrift            , ESYM);
    WPR(FirstHalfEtaKick         , ESYM);
    WPR(LastHalfEtaKick          , ESYM);
    WPR(ScaleFactorHalf          , ESYM);
    
    WPR(Do2LPTVelocityRereading  , ISYM);

    output_header = ss.str();
}

#undef WPR
#undef WPRS
#define WPR(X,XSYM) fprintf(statefp, PRQUOTEME(%22s = %XSYM\n), PRQUOTEME(X), X); 
#define WPRS(X,XSYM) fprintf(statefp, "%22s = \"%s\" \n", PRQUOTEME(X), X); 
// #define WPR(X,XSYM) fprintf(statefp, PRQUOTEME(X = %XSYM\n), X)
// #define WPRS(X,XSYM) fprintf(statefp, PRQUOTEME(X) " = \"%s\" \n", X)

void State::write_to_file(const char *dir) {
    char statefn[1050];
    sprintf(statefn,"%s/state",dir);
    FILE *statefp;
    statefp = fopen(statefn,"wb");
    assertf(statefp!=NULL, "Couldn't open file %s to write state\n", statefn);

    WPR(np_state                       , llu);
    WPR(cpd_state                      , ISYM);
    WPR(order_state                    , ISYM);

    fprintf(statefp,"%s", header());

    WPR(MaxVelocity		 , FSYM);
    WPR(MaxAcceleration		 , FSYM);
    WPR(RMS_Velocity             , FSYM);
    WPR(MinVrmsOnAmax		 , FSYM);
    WPR(MaxCellSize              , ISYM);
    WPR(MinCellSize              , ISYM);
    WPR(StdDevCellSize           , FSYM);

    time_t now  = time(0);
    fprintf(statefp,"#State written:%s\n",asctime(localtime(&now)) );
    FinalizeHeader(statefp);

    fclose(statefp);
}

void State::AssertStateLegal(Parameters &P) {
    //make sure read state and parameters are compatible
    assertf(order_state == P.order, 
	    "State and Parameter order do not match, %d != %d\n", 
	    order_state, P.order);
    assertf(cpd_state == P.cpd, 
	    "State and Parameter cpd do not match, %d != %d\n", 
	    cpd_state, P.cpd);
    assertf(np_state == P.np, 
	    "State and Parameter np do not match, %d != %d\n", 
	    np_state, P.np);
    assertf(MaxCellSize < (2.048e9)/3,
	    "The largest cell has %d particles, exceeding the allowed amount.\n",
	    MaxCellSize);
}


#endif //STATESTRUCTURECPP
