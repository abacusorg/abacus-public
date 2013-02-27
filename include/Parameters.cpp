/* Parameters.cpp

The Parameters class contains the time-independent global variables
for the simulation.  These get read from an ASCII parameter file
via ParseHeader.

NB: When adding parameters, you should add:

1) a variable to the class definition,
2) an installscalar/vector line to the constructor, 
3) (optional) a validation check in ValidateParameters.

*/



#ifdef OLD_STUFF
#define QUOTEME(X) #X
#define SCANLINE(X, XSYM)   ret += sscanf(line, QUOTEME(X  =  %XSYM),   &X);
#define STRSCANLINE(X)      ret += sscanf(line, QUOTEME(X  =  %s), X);
#define DUMPVAR(X, XSYM)    fprintf(stderr, QUOTEME(X  =  %XSYM\n), X);

#define TESTSTRINGUNDEFINED(X) if(strcmp(X,"NotDefined")==0) { fprintf(stderr,QUOTEME(You didnt define X\n)); assert(1==0); }

// This had been testing against 0, but that would exclude legal negative entries!
// This code will crash any entries less than -1 million.
#define NO_ENTRY -1237654
#define TestVariablePresent(variable) if(variable<=NO_ENTRY+0.1) { \
     printf("Can't find paramter <"#variable"> in the parameter file\n");  \
     assert(1==0); }
#endif // OLD_STUFF




// TODO: Is it ok that this uses floats not doubles?
#define MAX_LINE_LENGTH 1024
#define STRUNDEF "NONE"

#include "ParseHeader.hh"

class Parameters: public ParseHeader {
public:
    
	char RunName[1024]; //What to call this run

    int np;
    int cpd;
    int order;
    int ppd; //Exact only if np is a perfect cube

    int NearFieldRadius;    // Radius of cells in the near-field
    float SofteningLength; // Softening length in units of interparticle spacing

    int  DerivativeExpansionRadius;
    int  MAXConvolutionRAMMB;
    int  ConvolutionCacheSizeMB;

    int  DirectNewtonRaphson;  // 0 or 1 
    int  DirectDoublePrecision; // 0 or 1 

    char DerivativesDirectory[1024];

    char InitialConditionsDirectory[1024];   // The initial condition file name
    char ICFormat[1024];		// The format of the IC files
    float ICPositionRange;		// The box size of the IC positions, 
    	// in file units.  If ==0, then will default to BoxSize;
    float ICVelocity2Displacement;	// The conversion factor from file velocities
    	// to redshift-space comoving displacements (at the IC redshift!).
    float NumSlabsInsertList;		
         // The amount of space to allocate for the insert list, in units 
	 // of np/cpd particles.  Set =0 to allocate the full np particles.
	 // Default is 2.
    float NumSlabsInsertListIC;		
         // The amount of space to allocate for the insert list, in units 
	 // of np/cpd particles.  Set =0 to allocate the full np particles.
	 // This parameter applies only to the IC step.
	 // Recommend either 4 or 0.

    char ReadStateDirectory[1024];  // Where the input State lives
    char WriteStateDirectory[1024]; // Where the output State lives
    char PastStateDirectory[1024];  // Where the old input State lives
    char WorkingDirectory[1024];	// If Read/Write/Past not specified, where to put the states
    char LogFileDirectory[1024];
    char OutputDirectory[1024];     // Where the outputs go
    char BaseDistributionDirectory[1024];

    char DumpFilePrefix[1024];      // What the outputs are called
    char GroupFilePrefix[1024];     // What the group outputs are called
    char LightFilePrefix[1024];

    float Dumpz[1024];
    int nDumpz;

    float H0;          // The Hubble constant in km/s/Mpc
    float Omega_M;
    float Omega_DE;
    float Omega_K;
    float w0;          // w(z) = w_0 + (1-a)*w_a
    float wa;

    float BoxSize;
    int hMpc;           // =1 if we're using Mpc/h units.  =0 if Mpc units
    float InitialRedshift;
    int LagrangianPTOrder;  // =1 for Zel'dovich, =2 for 2LPT, =3 for 3LPT

    int GroupRadius;        // Maximum size of a group, in units of cell sizes
    float Eta;         // Time-step parameter based on accelerations
    float Dlna;        // Maximum time step in d(ln a)
    // Could have microstepping instructions
    // Could have group finding or coevolution set instructions

    int  NLightCones;
    double LightConeOrigins[24];
    char LightConeDirectory[1024];

    int StoreForces; // If 1, store the accelerations
    int ForcesOnly; //If 1, do not drift or kick
    int ForceOutputDebug; // If 1, output near and far forces seperately. Should only be set if ForcesOnly is also set


    Parameters() {

    	installscalar("NP",np, MUST_DEFINE);
    	installscalar("CPD",cpd,MUST_DEFINE);
    	installscalar("Order",order,MUST_DEFINE);

    	installscalar("NearFieldRadius",NearFieldRadius,MUST_DEFINE);    // Radius of cells in the near-field
    	installscalar("SofteningLength", SofteningLength,MUST_DEFINE); // Softening length in units of interparticle spacing

    	installscalar("DerivativeExpansionRadius", DerivativeExpansionRadius,MUST_DEFINE);
    	installscalar("MAXConvolutionRAMMB", MAXConvolutionRAMMB,MUST_DEFINE);
    	installscalar("ConvolutionCacheSizeMB", ConvolutionCacheSizeMB,MUST_DEFINE);

    	installscalar("DirectNewtonRaphson",DirectNewtonRaphson,MUST_DEFINE);  // 0 or 1
    	installscalar("DirectDoublePrecision",DirectDoublePrecision,MUST_DEFINE); // 0 or 1

    	installscalar("DerivativesDirectory",DerivativesDirectory,MUST_DEFINE);

    	installscalar("InitialConditionsDirectory",InitialConditionsDirectory,MUST_DEFINE);   // The initial condition file name
    	installscalar("ICFormat",ICFormat,MUST_DEFINE);   // The initial condition file format
    	installscalar("ICPositionRange",ICPositionRange,MUST_DEFINE);   // The initial condition file position convention
    	installscalar("ICVelocity2Displacement",ICVelocity2Displacement,MUST_DEFINE);   // The initial condition file velocity convention

	NumSlabsInsertList = 2.0;
    	installscalar("NumSlabsInsertList",NumSlabsInsertList,DONT_CARE);   
	NumSlabsInsertListIC = 0.0;
    	installscalar("NumSlabsInsertListIC",NumSlabsInsertListIC,DONT_CARE);   

    	sprintf(ReadStateDirectory,STRUNDEF);
    	sprintf(WriteStateDirectory,STRUNDEF);
    	sprintf(PastStateDirectory,STRUNDEF);
    	sprintf(WorkingDirectory,STRUNDEF);
    	installscalar("ReadStateDirectory",ReadStateDirectory,DONT_CARE);  // Where the input State lives
    	installscalar("WriteStateDirectory",WriteStateDirectory,DONT_CARE); // Where the output State lives
    	installscalar("PastStateDirectory",PastStateDirectory,DONT_CARE);  // Where the old input State lives
    	installscalar("WorkingDirectory",WorkingDirectory,DONT_CARE);
    	installscalar("LogFileDirectory",LogFileDirectory,MUST_DEFINE);
    	installscalar("OutputDirectory",OutputDirectory,MUST_DEFINE);     // Where the outputs go

    	installscalar("BaseDistributionDirectory",BaseDistributionDirectory,MUST_DEFINE);

    	installscalar("DumpFilePrefix",DumpFilePrefix,MUST_DEFINE);      // What the outputs are called
    	installscalar("GroupFilePrefix",GroupFilePrefix,MUST_DEFINE);     // What the group outputs are called
    	installscalar("LightConeDirectory",LightConeDirectory,MUST_DEFINE); //Where the lightcones go. Generally will be the same as the Output directory
    	installscalar("NLightCones",NLightCones,DONT_CARE); //if not set, we assume 0

    	installvector("Dumpz",Dumpz,1024,1,MUST_DEFINE);
    	installscalar("nDumpz",nDumpz,MUST_DEFINE);

    	installscalar("H0", H0, MUST_DEFINE);
    	installscalar("Omega_M", Omega_M, MUST_DEFINE);
    	installscalar("Omega_DE", Omega_DE, MUST_DEFINE);
    	installscalar("Omega_K", Omega_K, MUST_DEFINE);
    	installscalar("w0", w0, MUST_DEFINE);
    	installscalar("wa", wa, MUST_DEFINE);

    	installscalar("BoxSize",BoxSize,MUST_DEFINE);
    	installscalar("hMpc",hMpc,MUST_DEFINE);           // =1 if we're using Mpc/h units.  =0 if Mpc units
    	installscalar("InitialRedshift",InitialRedshift,MUST_DEFINE);
    	installscalar("LagrangianPTOrder",LagrangianPTOrder,MUST_DEFINE);  // =1 for Zel'dovich, =2 for 2LPT, =3 for 3LPT

    	installscalar("GroupRadius",GroupRadius,MUST_DEFINE);        // Maximum size of a group, in units of cell sizes
    	installscalar("Eta",Eta,MUST_DEFINE);         // Time-step parameter based on accelerations
    	installscalar("Dlna",Dlna,MUST_DEFINE);        // Maximum time step in d(ln a)

    	StoreForces = 1;
    	installscalar("StoreForces",StoreForces, DONT_CARE);
    	ForcesOnly = 0;
    	installscalar("ForcesOnly",ForcesOnly, DONT_CARE);
    	ForceOutputDebug = 0;
    	installscalar("ForceOutputDebug",ForceOutputDebug,DONT_CARE);
    	installscalar("RunName",RunName,MUST_DEFINE);

    	double npcr = pow(np,1.0/3.0);
    	ppd = (long long int) floor(npcr+0.5);
    	assert(ppd >0);


    }

    void ReadParameters(char *paramaterfile, int icflag);
    // void DumpParameters(void);
    // void CheckVariablesPresent(void);
    void ValidateParameters(void);
    void register_vars();
private:
    void ProcessStateDirectories();
};

void Parameters::ProcessStateDirectories(){
	if (strcmp(WorkingDirectory,STRUNDEF) !=0){
		if ( strcmp(ReadStateDirectory,STRUNDEF)!=0 || strcmp(WriteStateDirectory,STRUNDEF)!=0 || strcmp(PastStateDirectory,STRUNDEF)!=0   ){
			QUIT("If WorkingDirectory is defined, Read/Write/PastStateDirectory should be undefined. Terminating\n")
		}
		else{
			sprintf(ReadStateDirectory,"%s/read",WorkingDirectory);
			sprintf(WriteStateDirectory,"%s/write",WorkingDirectory);
			sprintf(PastStateDirectory,"%s/past",WorkingDirectory);
		}
	}
}


void Parameters::ReadParameters(char *parameterfile, int icflag) {
    HeaderStream hs(parameterfile);
    ReadHeader(hs);
    ProcessStateDirectories();
    if(!icflag) ValidateParameters();
}

void Parameters::ValidateParameters(void) {
    // Warning: Can't use STDLOG(), QUIT(), or assertf() calls in here:
    // We haven't opened the stdlog file yet!

    if(np<0) {
        fprintf(stderr,
            "[ERROR] np = %d must be greater than zero!\n", np);
        assert(1==0);
    }


    if(cpd<0) {
        fprintf(stderr,
            "[ERROR] cpd = %d must be greater than zero!\n", cpd);
        assert(1==0);
    }

    if( !( (DirectNewtonRaphson == 1) ||  (DirectNewtonRaphson == 0) ) ) {
        fprintf(stderr,"DirectNewtonRapson must be 0 or 1\n");
        assert(1==0);
    }

    if( !( (DirectDoublePrecision == 1) || (DirectDoublePrecision==0) ) ) {
        fprintf(stderr,"DirectDoublePrecision must be 0 or 1 \n");
        assert(1==0);
    }

    if(cpd%2==0) {
        fprintf(stderr,
            "[ERROR] cpd = %d  must be odd!\n", cpd);
        assert(1==0);
    }

    if(NearFieldRadius<=0) {
        fprintf(stderr, "[ERROR] NearFieldRadius = %d  must be greater than 0\n",
            NearFieldRadius );
        assert(1==0);
    }

    if(NearFieldRadius>(cpd-1)/2) {
        fprintf(stderr,
            "[ERROR] NearFieldRadius = %d must be less than (cpd-1)/2 = %d\n",
                    NearFieldRadius, (cpd-1)/2  );
        assert(1==0);
    }

    if(GroupRadius<=0) {
        fprintf(stderr, "[ERROR] GroupRadius = %d  must be greater than 0\n",
            GroupRadius );
        assert(1==0);
    }

    if(order>16 || order < 2 ) {
        fprintf(stderr,
            "[ERROR] order = %d must be less than or equal to 16 and greater than 1\n",
            order);
        assert(1==0);
    }

    if(MAXConvolutionRAMMB<0) {
        fprintf(stderr,
            "[ERROR] MAXConvolutionRAMMB = %d must be greater than 0 \n",
                MAXConvolutionRAMMB);
        assert(1==0);
    }

    if(ConvolutionCacheSizeMB<=0) {
        fprintf(stderr,
            "[ERROR] ConvolutionCacheSizeMB = %d must be greater than 0 \n",
                ConvolutionCacheSizeMB);
        assert(1==0);
    }


    if( (DerivativeExpansionRadius!=8) &&
        (DerivativeExpansionRadius!=16) &&
        (DerivativeExpansionRadius!=32) ) {

        fprintf(stderr,
            "[ERROR] DerivativeExpansionRadius = %d has to be 8 or 16 or 32\n",
                DerivativeExpansionRadius);
        assert(1==0);
    }

    if( (SofteningLength <= 0) || (SofteningLength>1) ) {
        fprintf(stderr,
            "[ERROR] SofteningLength = %e has to be in (0,1)\n",
                SofteningLength);
        assert(1==0);
    }

    if (NumSlabsInsertList<0.0 || NumSlabsInsertList>cpd) {
        fprintf(stderr,
            "[ERROR] NumslabsInsertList = %e must be in range [0..CPD]\n",
                NumSlabsInsertList);
        assert(1==0);
    }

    if (NumSlabsInsertListIC<0.0 || NumSlabsInsertListIC>cpd) {
        fprintf(stderr,
            "[ERROR] NumslabsInsertListIC = %e must be in range [0..CPD]\n",
                NumSlabsInsertListIC);
        assert(1==0);
    }

    // Note: Not putting a requirement that ICVelocity2Displacement>0,
    // because one might use 0 or negative values to play with the ICs.

    // Illegal ICFormat's will crash in loadIC.cpp; no need to crash here.

    ExpandPathName(DerivativesDirectory);
    ExpandPathName(ReadStateDirectory);
    ExpandPathName(WriteStateDirectory);
    ExpandPathName(PastStateDirectory);
    ExpandPathName(OutputDirectory);
    ExpandPathName(LogFileDirectory);
    ExpandPathName(InitialConditionsDirectory);

    CheckDirectoryExists(DerivativesDirectory);
    CheckDirectoryExists(ReadStateDirectory);
    CheckDirectoryExists(WriteStateDirectory);
    CheckDirectoryExists(PastStateDirectory);
    CheckDirectoryExists(OutputDirectory);
    CheckDirectoryExists(LogFileDirectory);
    CheckDirectoryExists(BaseDistributionDirectory);

    CheckDirectoryExists(InitialConditionsDirectory);


    char dfn[1024];

    for(int i=0;i<(cpd+1)/2;i++) {
        sprintf(dfn,"%s/fourierspace_%d_%d_%d_%d_%d",
            DerivativesDirectory,
            cpd,
            order,
            NearFieldRadius,
            DerivativeExpansionRadius,i );

        CheckFileExists(dfn);
    }

    if (ForceOutputDebug && !StoreForces){
    	QUIT("ForcesOutputDebug set to on, but StoreForces was not set. This is not supported.\n")
    }
    if (ForceOutputDebug && !ForcesOnly){
    	QUIT("ForcesOutputDebug set to on, but ForcesOnly was not set. This is not supported.\n")
    }



}
Parameters P;
