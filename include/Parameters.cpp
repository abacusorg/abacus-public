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




#define MAX_LINE_LENGTH 1024
#define STRUNDEF "NONE"

#include "ParseHeader.hh"

class Parameters: public ParseHeader {
public:
    
    char SimName[1024]; //What to call this run
    // TODO: Rename this to SimName

    long long int np;
    int cpd;
    int order;

    int NearFieldRadius;    // Radius of cells in the near-field
    double SofteningLength; // Softening length in units of interparticle spacing

    int  DerivativeExpansionRadius;
    int  MAXConvolutionRAMMB;
    int  ConvolutionCacheSizeMB;
    int RamDisk;	// ==0 for a normal disk, ==1 for a ramdisk (which don't have DIO support)
    int ForceBlockingIO;   // ==1 if you want to force all IO to be blocking.

    int  DirectNewtonRaphson;  // 0 or 1 

    char DerivativesDirectory[1024];

    char InitialConditionsDirectory[1024];   // The initial condition file name
    char ICFormat[1024];		// The format of the IC files
    double ICPositionRange;		// The box size of the IC positions, 
    	// in file units.  If ==0, then will default to BoxSize;
    double ICVelocity2Displacement;	// The conversion factor from file velocities
    	// to redshift-space comoving displacements (at the IC redshift!).
	// =-1 to instead supply file velocities in km/s.
	// =0 to force the initial velocities to zero!
    double NumSlabsInsertList;		
         // The amount of space to allocate for the insert list, in units 
	 // of np/cpd particles.  Set =0 to allocate the full np particles.
	 // Default is 2.
    double NumSlabsInsertListIC;		
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

    char TimeSliceFilePrefix[1024];      // What the outputs are called
    char GroupFilePrefix[1024];     // What the group outputs are called
    char LightFilePrefix[1024];

    char OutputFormat[1024];		// The format of the Output files
    int  OmitOutputHeader;		// =1 if you want to skip the ascii header

    double TimeSlicez[1024];
    int nTimeSlice;

    double H0;          // The Hubble constant in km/s/Mpc
    double Omega_M;
    double Omega_DE;
    double Omega_K;
    double w0;          // w(z) = w_0 + (1-a)*w_a
    double wa;

    double BoxSize;
    int hMpc;           // =1 if we're using Mpc/h units.  =0 if Mpc units
    double InitialRedshift;
    int LagrangianPTOrder;  // =1 for Zel'dovich, =2 for 2LPT, =3 for 3LPT

    int GroupRadius;        // Maximum size of a group, in units of cell sizes
    double Eta;         // Time-step parameter based on accelerations
    	// TODO: Rename to TimeStepAccel
    double Dlna;        // Maximum time step in d(ln a)
    	// Rename to TimeStepDlna

    // Could have microstepping instructions
    // Could have group finding or coevolution set instructions

    int  NLightCones;
    double LightConeOrigins[24];
    char LightConeDirectory[1024];

    int LogVerbosity;   // If 0, production-level log; higher numbers are more verbose
    int StoreForces; // If 1, store the accelerations
    int ForceOutputDebug; // If 1, output near and far forces seperately. 


    Parameters() {

    	installscalar("NP",np, MUST_DEFINE);
    	installscalar("CPD",cpd,MUST_DEFINE);
    	installscalar("Order",order,MUST_DEFINE);

    	installscalar("NearFieldRadius",NearFieldRadius,MUST_DEFINE);    // Radius of cells in the near-field
    	installscalar("SofteningLength", SofteningLength,MUST_DEFINE); // Softening length in units of interparticle spacing

    	installscalar("DerivativeExpansionRadius", DerivativeExpansionRadius,MUST_DEFINE);
    	installscalar("MAXConvolutionRAMMB", MAXConvolutionRAMMB,MUST_DEFINE);
    	installscalar("ConvolutionCacheSizeMB", ConvolutionCacheSizeMB,MUST_DEFINE);
	RamDisk = 0;
    	installscalar("RamDisk",RamDisk,DONT_CARE);
	ForceBlockingIO = 0;
    	installscalar("ForceBlockingIO",ForceBlockingIO,DONT_CARE);

    	installscalar("DirectNewtonRaphson",DirectNewtonRaphson,MUST_DEFINE);  // 0 or 1

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


    	installscalar("TimeSliceFilePrefix",TimeSliceFilePrefix,MUST_DEFINE);      // What the outputs are called
    	installscalar("GroupFilePrefix",GroupFilePrefix,MUST_DEFINE);     // What the group outputs are called
    	installscalar("LightConeDirectory",LightConeDirectory,MUST_DEFINE); //Where the lightcones go. Generally will be the same as the Output directory
    	installscalar("NLightCones",NLightCones,DONT_CARE); //if not set, we assume 0

    	installvector("TimeSlicez",TimeSlicez,1024,1,MUST_DEFINE);
    	installscalar("nTimeSlice",nTimeSlice,MUST_DEFINE);

	// strcpy(OutputFormat,"RVdouble");
	strcpy(OutputFormat,"Packed");
    	installscalar("OutputFormat",OutputFormat,DONT_CARE);
	OmitOutputHeader = 0;
    	installscalar("OmitOutputHeader",OmitOutputHeader,DONT_CARE);

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

    	LogVerbosity = 1;
    	installscalar("LogVerbosity",LogVerbosity, DONT_CARE);
    	StoreForces = 0;
    	installscalar("StoreForces",StoreForces, DONT_CARE);
    	ForceOutputDebug = 0;
    	installscalar("ForceOutputDebug",ForceOutputDebug,DONT_CARE);
    	installscalar("SimName",SimName,MUST_DEFINE);
	hs = NULL;
    }

    // We're going to keep the HeaderStream, so that we can output it later.
    HeaderStream *hs;
    ~Parameters(void) { delete hs; }
    char *header() { 
	assert(hs!=NULL); assert(hs->buffer!=NULL);
        return hs->buffer;	// This is just a standard C-style string.
    }


    void ReadParameters(char *paramaterfile, int icflag);
    void ValidateParameters(void);
    void register_vars();

    double ppd() {
        // return the cube root of np, but be careful to avoid round-off 
	// of perfect cubes.
	double _ppd = pow((double)np, 1.0/3.0);
	if ( fabs(_ppd-floor(_ppd+0.1))<1e-10 ) _ppd = floor(_ppd+0.1);
	return _ppd;
    }
    int is_np_perfect_cube() {
	// Return 1 if np is a perfect cube.
        int n = floor(ppd());
	if (n*n*n!=np) return 1; else return 0;
    }

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
    hs = new HeaderStream(parameterfile);
    ReadHeader(*hs);
    hs->Close();
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

    if (ForceOutputDebug && StoreForces) {
    	fprintf(stderr,"ForcesOutputDebug and StoreForces both set. This is not supported.\n");
	assert(1==0);
    }

    if (LogVerbosity<0) {
        fprintf(stderr,"LogVerbosity must be >=0\n");
	assert(1==0);
    }



}
Parameters P;
