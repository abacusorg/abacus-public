// Copyright 2012-2025 The Abacus Developers
// SPDX-License-Identifier: GPL-3.0-or-later

/* statestructure.cpp
 *
 * This defines the State class.
 *
 * States contain the time-dependent global information about the
 * simulation.  The code is evolving from the ReadState at time j to
 * the WriteState at time j+1.  The state class is written (in ASCII)
 * to the state file in the state subdirectory, so as to track the 
 * time-varying aspects of the simulation.
 * 
 * Time-independent instructions for the simulation should be in the
 * Parameter class.
 * 
 * State uses ParseHeader to handle the reading of the ASCII file.
 * However, it must provide its own routine to write the new file!
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

    uint64 np_state;
    uint64 np_with_ghost_state;
    uint64 np_subA_state;
    uint64 np_subB_state;
    int cpd_state;
    int order_state;
    
    fs::path ParameterFileName;   // State must contain a pointer to the Parameter file
    std::string CodeVersion;
    std::string OutputFormatVersion;
    std::string RunTime;
    std::string MachineName;
    int NodeRankX;   // The MPI X rank, 0 if serial
    int NodeRankZ;   // The MPI Z rank, 0 if serial or 1D
    int NodeSizeX;   // The MPI X size, 1 if serial
    int NodeSizeZ;   // The MPI Z size, 1 if serial or 1D
    double ppd;		// Particles per dimension
    int64 ippd;     // The closest integer value to NP^(1/3)
    int DoublePrecision;  // =1 if code is using double precision positions
    std::string SofteningType;  // The force law.  This is here because it's a compile-time parameter.
    double SofteningLengthNow;  // Effective Plummer length, used for timestepping.  Same units as BoxSize.
    double SofteningLengthNowInternal;  // The equivalent length for the current softening technique.  Same units as BoxSize.

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
    double Growth_on_a_n;  // D/a^n
    double f_growth;
    double w;
    double HubbleNow;   // In H_0=1 units
    double Htime;       // Time*H(z)
    double OmegaNow_m;
    double OmegaNow_K;
    double OmegaNow_DE;
    double fsmooth;     // This is time-independent: Omega_Smooth/Omega_M
    double CoordinateDistanceHMpc;   // In Mpc/h

    // More description of the last time step
    double DeltaTime;
    double DeltaScaleFactor;
    double DeltaRedshift;

    // The FOF density scale being used (in code units)
    // This matters because 0 indicates that it was not computed.
    double DensityKernelRad2;
    // Unit density values for our kernel. Used in group finding and density aux packing
    FLOAT FOFunitdensity;
    FLOAT invFOFunitdensity;
    // The density threshold for L0 particle eligibility (units of cosmic mean)
    double L0DensityThreshold;
    double SODensityL1; //density threshold for SO L1 groups.
    double SODensityL2; //density threshold for SO L2 groups. 

    // Some statistics about the particle distribution that might be useful.
    double MaxVelocity;
    double MaxAcceleration;
    double MinVrmsOnAmax;
    int MaxCellSize;
    int MinCellSize;
    double StdDevCellSize;
    double RMS_Velocity;
    int MaxGroupDiameter; 
    int MaxL0GroupSize;
    double DirectsPerParticle;
    // The variables below are not intended for output.  Just used in the code.

    // We will write a lot of state information into the header of output
    // files.  We will collect that information here.
    std::string output_header;
    void make_output_header();
    std::string header() const { 
    	return output_header;
    }
    int DoTimeSliceOutput;
    int OutputIsAllowed;
    int DoBinning;
    int DoGroupFindingOutput;
    int DoSubsampleOutput; 
    int VelIsSynchronous;
    int HaveAuxDensity;
    int SetAuxDensity;          // Should the Kick store the density in the aux?
    int DidGroupFindingOutput;  // did we already do group output on the positions in this state?
    int LastTimeSliceOutput;  // index of the most recent time slice
    int LastSubsampleOutput;  // index of the most recent time slice subsample
    std::string GroupFindingDensitySource;
    
    double LPTVelScale;  // normalization for the aux compression of the LPT vel

    int OverwriteState;
    int OverwriteConvState;
    int StripeState;
    int StripeConvState;

    std::string Pipeline;

    fs::path LogDirectory;  // step-numbered log directory

    int64 np_lightcone;

    int GhostRadius;

    void read_from_file(const fs::path &fn);
    void write_to_file(const fs::path &dir, const fs::path &fname);
    void write_to_file(const fs::path &dir) { write_to_file(dir,""); }
    std::string get_state_string() const;
    
    State();
    
    void AssertStateLegal(Parameters &P);

};

    State::State(){
        // These state values are not written to file headers, only to the statefile proper
        // But we might want to load a file header as a state for analysis!
        // So make these optional, but they will fail validation against the Parameters if missing in a simulation context
    	installscalar("np_state",np_state,DONT_CARE);
        np_with_ghost_state = 0;
        installscalar("np_with_ghost_state",np_with_ghost_state,DONT_CARE);
        np_subA_state = 0;
        installscalar("np_subA_state",np_subA_state,DONT_CARE);
        np_subB_state = 0;
        installscalar("np_subB_state",np_subB_state,DONT_CARE);
    	installscalar("cpd_state",cpd_state,DONT_CARE);
    	installscalar("order_state",order_state,DONT_CARE);

        installscalar("Pipeline",Pipeline,DONT_CARE);
    	installscalar("ParameterFileName",ParameterFileName,DONT_CARE);
    	installscalar("ppd",ppd,DONT_CARE);
        installscalar("SofteningType", SofteningType,DONT_CARE);
        installscalar("SofteningLengthNow", SofteningLengthNow, MUST_DEFINE);
        installscalar("SofteningLengthNowInternal", SofteningLengthNowInternal,MUST_DEFINE);

    	CodeVersion = "version_not_defined";
    	installscalar("CodeVersion",CodeVersion,DONT_CARE);
    	OutputFormatVersion = "version_not_defined";
    	installscalar("OutputFormatVersion",OutputFormatVersion,DONT_CARE);
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
    	installscalar("Growth_on_a_n",Growth_on_a_n,DONT_CARE);  // D/a^n
    	installscalar("f_growth",f_growth,DONT_CARE);
    	installscalar("w",w,DONT_CARE);
    	installscalar("HubbleNow",HubbleNow,DONT_CARE);   // In km/s/Mpc
    	installscalar("Htime",Htime,DONT_CARE);       // Time*H(z)
    	installscalar("OmegaNow_m",OmegaNow_m,DONT_CARE);
    	installscalar("OmegaNow_K", OmegaNow_K,DONT_CARE);
    	installscalar("OmegaNow_DE",OmegaNow_DE,DONT_CARE);
    	installscalar("fsmooth",fsmooth,DONT_CARE);
    	installscalar("CoordinateDistanceHMpc",CoordinateDistanceHMpc,DONT_CARE);
    	installscalar("DeltaTime",DeltaTime,DONT_CARE);
    	installscalar("DeltaScaleFactor",DeltaScaleFactor,DONT_CARE);
    	installscalar("DeltaRedshift",DeltaRedshift,DONT_CARE);
    	installscalar("DensityKernelRad2",DensityKernelRad2,DONT_CARE);
        installscalar("SODensityL1",SODensityL1,DONT_CARE);
        installscalar("SODensityL2",SODensityL2,DONT_CARE);
        
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
        MaxGroupDiameter = 0; 
        installscalar("MaxGroupDiameter",MaxGroupDiameter,DONT_CARE);

        MaxL0GroupSize = 0;
        installscalar("MaxL0GroupSize",MaxL0GroupSize,DONT_CARE);
        DirectsPerParticle = 0.0;
        installscalar("DirectsPerParticle",DirectsPerParticle,DONT_CARE);

        GhostRadius = 0;
        installscalar("GhostRadius",GhostRadius,DONT_CARE);

        // Initialize helper variables
        DoTimeSliceOutput = 0;
        installscalar("DoTimeSliceOutput",DoTimeSliceOutput,DONT_CARE);

        OutputIsAllowed = 0;
        
        DoGroupFindingOutput = 0;
        installscalar("DoGroupFindingOutput",DoGroupFindingOutput,DONT_CARE);
        
        DoSubsampleOutput = 0;
        installscalar("DoSubsampleOutput",DoSubsampleOutput,DONT_CARE);

        VelIsSynchronous = 0;
        installscalar("VelIsSynchronous",VelIsSynchronous,DONT_CARE);

        HaveAuxDensity = 0;
        installscalar("HaveAuxDensity",HaveAuxDensity,DONT_CARE);

        SetAuxDensity = 0;
        installscalar("SetAuxDensity",SetAuxDensity,DONT_CARE);

        DidGroupFindingOutput = 0;
        installscalar("DidGroupFindingOutput",DidGroupFindingOutput,DONT_CARE);

        LastTimeSliceOutput = -1;
        installscalar("LastTimeSliceOutput",LastTimeSliceOutput,DONT_CARE);

        LastSubsampleOutput = -1;
        installscalar("LastSubsampleOutput",LastSubsampleOutput,DONT_CARE);

        GroupFindingDensitySource[0] = '\0';
        installscalar("GroupFindingDensitySource",GroupFindingDensitySource,DONT_CARE);

        LPTVelScale = 0.;
        installscalar("LPTVelScale",LPTVelScale,DONT_CARE);

        // These will be set in InitWriteState() based on StateIOMode and Conv_IOMode in the parameters file
        OverwriteState = 0;
        OverwriteConvState = 0;
        StripeState = 0;
        StripeConvState = 0;

        np_lightcone = 0;
    }


void State::read_from_file(const fs::path &fn) {
    fs::path statefn = fn / "state";
    HeaderStream hs(statefn);
    ReadHeader(hs);
    hs.Close();
}



#define FSYM 20.15f
#define ESYM 25.15e
#define ISYM d
#define SSYM s


#define PRQUOTEME(X) #X
#define WPR(X,XSYM) {ss << fmt::format(PRQUOTEME({:>26} = {:XSYM}\n), PRQUOTEME(X), X);}
#define WPRS(X,XSYM) {ss << fmt::format("{:>26} = \"{}\" \n", PRQUOTEME(X), X);}

void State::make_output_header() {
    // We're going to output most, but not all of the fields, into a 
    // nice header.  This will be prepended to many output files.
    // It also will get used to write the state.
    
    // Note that this is written *during* the timestep, so only properties
    // that are known at the beginning of the timestep (i.e. those setup in
    // BuildWriteState) should be in here.
    
    std::stringstream ss;

    WPRS(Pipeline                 , s);
    WPRS(ParameterFileName        , s);
    WPRS(CodeVersion              , s);
    WPRS(OutputFormatVersion      , s);
    WPRS(RunTime                  , s);
    WPRS(MachineName              , s);
    WPR(NodeRankX                , ISYM);
    WPR(NodeRankZ                , ISYM);
    WPR(NodeSizeX                , ISYM);
    WPR(NodeSizeZ                , ISYM);
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
    WPR(Growth_on_a_n            , FSYM);
    WPR(f_growth                 , FSYM);
    WPR(w                        , FSYM);
    WPR(HubbleNow                , FSYM);
    WPR(Htime                    , FSYM);
    WPR(OmegaNow_m               , FSYM);
    WPR(OmegaNow_K               , FSYM);
    WPR(OmegaNow_DE              , FSYM);
    WPR(fsmooth                  , FSYM);
    WPR(CoordinateDistanceHMpc   , FSYM);
    
    WPRS(SofteningType           , s);
    WPR(SofteningLengthNow       , ESYM);
    WPR(SofteningLengthNowInternal, ESYM);

    WPR(DeltaTime                , FSYM);
    WPR(DeltaScaleFactor         , FSYM);
    WPR(DeltaRedshift            , FSYM);

    WPR(DeltaEtaDrift            , ESYM);
    WPR(FirstHalfEtaKick         , ESYM);
    WPR(LastHalfEtaKick          , ESYM);
    WPR(ScaleFactorHalf          , ESYM);
    
    WPR(DensityKernelRad2        , FSYM);
    WPR(L0DensityThreshold       , FSYM);
    WPR(SODensityL1              , FSYM);
    WPR(SODensityL2              , FSYM);

    WPR(GhostRadius              , ISYM);
    
    WPR(DoTimeSliceOutput        , ISYM);
    WPR(DoSubsampleOutput        , ISYM);
    WPR(VelIsSynchronous         , ISYM);
    WPR(DoGroupFindingOutput     , ISYM);
    WPR(HaveAuxDensity           , ISYM);
    WPR(SetAuxDensity            , ISYM);
    WPR(DidGroupFindingOutput    , ISYM);
    WPR(LastTimeSliceOutput      , ISYM);
    WPR(LastSubsampleOutput      , ISYM);

    WPRS(GroupFindingDensitySource, s);

    output_header = ss.str();
}

#undef WPR
#undef WPRS
#define WPR(X,XSYM) {ss << fmt::format(PRQUOTEME({:>26s} = {:XSYM}\n), PRQUOTEME(X), X);}
#define WPRS(X,XSYM) {ss << fmt::format("{:>26s} = \"{}\" \n", PRQUOTEME(X), X);}

std::string State::get_state_string() const {
    std::stringstream ss;
    
    // The quantities known at the beginning of the timestep
    ss << header();

    // The quantities known at the end of the timestep
    WPR(np_state, ISYM);
    WPR(np_with_ghost_state, ISYM);
    WPR(np_subA_state, ISYM);
    WPR(np_subB_state, ISYM);
    WPR(cpd_state, ISYM);
    WPR(order_state, ISYM);
    WPR(np_lightcone, ISYM);

    WPR(MaxVelocity, FSYM);
    WPR(MaxAcceleration, FSYM);
    WPR(RMS_Velocity, FSYM);
    WPR(MinVrmsOnAmax, FSYM);
    WPR(MaxCellSize, ISYM);
    WPR(MinCellSize, ISYM);
    WPR(StdDevCellSize, FSYM);
    WPR(MaxGroupDiameter, ISYM); 
    WPR(MaxL0GroupSize, ISYM); 
    WPR(DirectsPerParticle, FSYM);
    WPR(LPTVelScale, ESYM);

    time_t now = time(0);
    ss << fmt::format("#State written:{:s}\n", asctime(localtime(&now)));
    
    return ss.str();
}

void State::write_to_file(const fs::path &dir, const fs::path &suffix) {
    fs::path statefn = (dir / "state").string() + suffix.string();
    FILE *statefp;
    statefp = fopen(statefn.c_str(), "wb");
    assertf(statefp != NULL, "Couldn't open file {} to write state\n", statefn);

    std::string state_content = get_state_string();
    
    fwrite(state_content.c_str(), sizeof(char), state_content.length(), statefp);
    FinalizeHeader(statefp);
    
    fclose(statefp);
}

void State::AssertStateLegal(Parameters &P) {
    //make sure read state and parameters are compatible
    assertf(order_state == P.order, 
            "State and Parameter order do not match, {:d} != {:d}\n", 
            order_state, P.order);
    assertf(cpd_state == P.cpd, 
            "State and Parameter cpd do not match, {:d} != {:d}\n", 
            cpd_state, P.cpd);
    assertf(np_state == P.np, 
            "State and Parameter np do not match, {:d} != {:d}\n", 
            np_state, P.np);
    assertf(MaxCellSize < (2.048e9)/3,
            "The largest cell has {:d} particles, exceeding the allowed amount.\n",
            MaxCellSize);
}


#endif //STATESTRUCTURECPP
