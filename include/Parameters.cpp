/* Parameters.cpp

The Parameters class contains the time-independent global variables
for the simulation.  These get read from an ASCII parameter file
via ParseHeader.

NB: When adding parameters, you should add:

1) a variable to the class definition,
2) an installscalar/vector line to the constructor, 
3) (optional) a validation check in ValidateParameters.

*/

#ifndef __PARAMETERS_CPP
#define __PARAMETERS_CPP

#define MAX_LINE_LENGTH 1024
#define STRUNDEF "NONE"

#include "ParseHeader.hh"

#include <sys/sysinfo.h>

class Parameters: public ParseHeader {
public:
    
    char SimName[1024]; //What to call this run
    long long int np;
    int cpd;
    int order;

    int NearFieldRadius;    // Radius of cells in the near-field
    double SofteningLength; // Softening length in the same units as BoxSize

    int  DerivativeExpansionRadius;
    int  MAXRAMMB;
    int  ConvolutionCacheSizeMB; // Set to manually override the detected cache size
    int RamDisk;	// ==0 for a normal disk, ==1 for a ramdisk (which don't have DIO support)
    int OverwriteState; // 0 for normal, separate read and write states; 1 to overwrite the read state to save space
    int ForceBlockingIO;   // ==1 if you want to force all IO to be blocking.
    
    int OMP_NUM_THREADS;  // Number of OpenMP threads.  0 does not modify the system value (usually OMP_NUM_THREADS, or all threads).
                        // Negative values use that many fewer than the max.

    int  DirectNewtonRaphson;  // 0 or 1 

    char DerivativesDirectory[1024];

    char InitialConditionsDirectory[1024];   // The initial condition file name
    char ICFormat[1024];		// The format of the IC files
    int FlipZelDisp;                    // If non-zero and using ICFormat = Zeldovich, flip the Zeldovich displacements
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
    char MultipoleDirectory[1024];
    char TaylorDirectory[1024];
    char WorkingDirectory[1024];	// If Read/Write/Past not specified, where to put the states
    char LogDirectory[1024];
    char OutputDirectory[1024];     // Where the outputs go
    int OutputEveryStep; //Force timeslices to be output every step if 1
    char OutputFormat[1024];		// The format of the Output files
    int  OmitOutputHeader;		// =1 if you want to skip the ascii header

    double FinalRedshift;	// When to stop.  This will override TimeSliceRedshifts.
    double TimeSliceRedshifts[1024];
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
    double TimeStepAccel;         // Time-step parameter based on accelerations
    double TimeStepDlna;        // Maximum time step in d(ln a)

    // Could have microstepping instructions
    // Could have group finding or coevolution set instructions

    int  NLightCones;

    double LightConeOrigins[24];  // Same units as BoxSize

    char LightConeDirectory[1024];

    int PowerSpectrumStepInterval;
    int PowerSpectrumN1d; //1D number of bins to use in the powerspectrum

    int LogVerbosity;   // If 0, production-level log; higher numbers are more verbose
    int StoreForces; // If 1, store the accelerations
    int ForceOutputDebug; // If 1, output near and far forces seperately. 

    int ForceCPU; //IF 1, force directs to be executed exclusively on cpu even if cuda is available
    int GPUMinCellSinks;// If AVX directs are compiled, cells with less than this many particles go to cpu
    int ProfilingMode;//If 1, enable profiling mode, i.e. delete the write-state after creating it to run repeatedly on same dat

    // in MB
    unsigned int getCacheSize(){
        unsigned int cache_size = 0;
        FILE *fp = 0;
        fp = fopen("/sys/devices/system/cpu/cpu0/cache/index3/size", "r");  // L3 cache size in KB
        if(fp){
            int nscan = fscanf(fp, "%dK", &cache_size);
            assert(nscan == 1);
            fclose(fp);
        }
        cache_size /= 1024; // to MB
        return cache_size;
    }
    
    // in MB
    int getRAMSize(){
        struct sysinfo info;
        int status = sysinfo(&info);
        assert(status == 0);
        return (int) ((double) info.totalram*info.mem_unit / 1024./1024.);
    }

    Parameters() {

    	installscalar("NP",np, MUST_DEFINE);
    	installscalar("CPD",cpd,MUST_DEFINE);
    	installscalar("Order",order,MUST_DEFINE);

    	installscalar("NearFieldRadius",NearFieldRadius,MUST_DEFINE);    // Radius of cells in the near-field
    	installscalar("SofteningLength", SofteningLength,MUST_DEFINE); // Softening length in the same units as BoxSize
    	installscalar("DerivativeExpansionRadius", DerivativeExpansionRadius,MUST_DEFINE);
        MAXRAMMB = getRAMSize();
    	installscalar("MAXRAMMB", MAXRAMMB, DONT_CARE);
        ConvolutionCacheSizeMB = getCacheSize();
    	installscalar("ConvolutionCacheSizeMB", ConvolutionCacheSizeMB, DONT_CARE);
        RamDisk = 0;
    	installscalar("RamDisk",RamDisk,DONT_CARE);
        OverwriteState = 0;
    	installscalar("OverwriteState",OverwriteState,DONT_CARE);
        ForceBlockingIO = 0;
    	installscalar("ForceBlockingIO",ForceBlockingIO,DONT_CARE);
        
        OMP_NUM_THREADS = 0;
        installscalar("OMP_NUM_THREADS",OMP_NUM_THREADS,DONT_CARE);

        DirectNewtonRaphson = 1;
    	installscalar("DirectNewtonRaphson",DirectNewtonRaphson,DONT_CARE);  // 0 or 1

    	installscalar("DerivativesDirectory",DerivativesDirectory,MUST_DEFINE);

    	installscalar("InitialConditionsDirectory",InitialConditionsDirectory,MUST_DEFINE);   // The initial condition file name
    	installscalar("ICFormat",ICFormat,MUST_DEFINE);   // The initial condition file format
    	installscalar("ICPositionRange",ICPositionRange,MUST_DEFINE);   // The initial condition file position convention
    	installscalar("ICVelocity2Displacement",ICVelocity2Displacement,MUST_DEFINE);   // The initial condition file velocity convention
        FlipZelDisp = 0;
        installscalar("FlipZelDisp",FlipZelDisp,DONT_CARE);   // Flip Zeldovich ICs

    	NumSlabsInsertList = 2.0;
    	installscalar("NumSlabsInsertList",NumSlabsInsertList,DONT_CARE);   
    	NumSlabsInsertListIC = 4.0;
    	installscalar("NumSlabsInsertListIC",NumSlabsInsertListIC,DONT_CARE);   

    	sprintf(ReadStateDirectory,STRUNDEF);
    	sprintf(WriteStateDirectory,STRUNDEF);
    	sprintf(PastStateDirectory,STRUNDEF);
        sprintf(WorkingDirectory,STRUNDEF);
        
    	installscalar("ReadStateDirectory",ReadStateDirectory,DONT_CARE);  // Where the input State lives
    	installscalar("WriteStateDirectory",WriteStateDirectory,DONT_CARE); // Where the output State lives
    	installscalar("PastStateDirectory",PastStateDirectory,DONT_CARE);  // Where the old input State lives
        installscalar("WorkingDirectory",WorkingDirectory,DONT_CARE);
        installscalar("MultipoleDirectory",MultipoleDirectory,MUST_DEFINE);
        installscalar("TaylorDirectory",TaylorDirectory,MUST_DEFINE);
    	installscalar("LogDirectory",LogDirectory,MUST_DEFINE);
    	installscalar("OutputDirectory",OutputDirectory,MUST_DEFINE);     // Where the outputs go
    	OutputEveryStep = 0;
    	installscalar("OutputEveryStep",OutputEveryStep,DONT_CARE);

    	installscalar("LightConeDirectory",LightConeDirectory,MUST_DEFINE); //Where the lightcones go. Generally will be the same as the Output directory
    	installscalar("NLightCones",NLightCones,DONT_CARE); //if not set, we assume 0
        installvector("LightConeOrigins",LightConeOrigins,24,1,DONT_CARE);

	FinalRedshift = -2.0;	// If <-1, then we will cascade back to the minimum of the TimeSliceRedshifts list
    	installscalar("FinalRedshift",FinalRedshift,DONT_CARE);
    	installvector("TimeSliceRedshifts",TimeSliceRedshifts,1024,1,MUST_DEFINE);
    	installscalar("nTimeSlice",nTimeSlice,MUST_DEFINE);

	strcpy(OutputFormat,"RVdouble");
	// strcpy(OutputFormat,"Packed");
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
    	installscalar("TimeStepAccel",TimeStepAccel,MUST_DEFINE);         // Time-step parameter based on accelerations
    	installscalar("TimeStepDlna",TimeStepDlna,MUST_DEFINE);        // Maximum time step in d(ln a)

    	LogVerbosity = 10000;
    	installscalar("LogVerbosity",LogVerbosity, DONT_CARE);
    	StoreForces = 0;
    	installscalar("StoreForces",StoreForces, DONT_CARE);
    	ForceOutputDebug = 0;

    	installscalar("ForceOutputDebug",ForceOutputDebug,DONT_CARE);
        ForceCPU = 0;
        installscalar("ForceCPU",ForceCPU,DONT_CARE);
        GPUMinCellSinks = 0;
        installscalar("GPUMinCellSinks",GPUMinCellSinks,DONT_CARE);
        ProfilingMode=0;
        installscalar("ProfilingMode",ProfilingMode,DONT_CARE);

    	installscalar("SimName",SimName,MUST_DEFINE);
    	installscalar("PowerSpectrumStepInterval",PowerSpectrumStepInterval,DONT_CARE);
    	installscalar("PowerSpectrumN1d",PowerSpectrumN1d,DONT_CARE);
    	PowerSpectrumStepInterval = -1; //Do not calculate OTF powerspectra
    	PowerSpectrumN1d = 1;
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
        long long int n = floor(ppd());
	if (n*n*n==np) return 1; else return 0;
    }

    double FinishingRedshift() {
	// Return the redshift where the code should halt.
	// FinalRedshift is the controlling item, but if that's
	// not set (value<=-1), then we will run to the minimum of the TimeSliceRedshifts list.
	// If no TimeSlices are requested, then z=0.
        if (FinalRedshift>-1) return FinalRedshift;
	if (nTimeSlice) {
	    double minz = 1e10;
	    for (int i=0; i<nTimeSlice; i++)
		if (TimeSliceRedshifts[i]<minz) minz = TimeSliceRedshifts[i];
	    return minz;
	}
	return 0.0;
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
            if(OverwriteState){
                strcpy(WriteStateDirectory, ReadStateDirectory);
            }
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
            "[ERROR] np = %lld must be greater than zero!\n", np);
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

    if(MAXRAMMB<0) {
        fprintf(stderr,
            "[ERROR] MAXRAMMB = %d must be greater than 0 \n",
                MAXRAMMB);
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

    if( (SofteningLength < 0) || (SofteningLength>BoxSize) ) {
        fprintf(stderr,
            "[ERROR] SofteningLength = %e has to be in [0,BoxSize)\n",
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

    if (nTimeSlice<0) {
        fprintf(stderr,"nTimeSlice must be >=0\n");
	assert(1==0);
    }
    
    if(abs(Omega_M + Omega_DE + Omega_K - 1.) > 1e-6){
        fprintf(stderr,"Omega_M + Omega_DE + Omega_K must equal 1, but is %g\n", Omega_M + Omega_DE + Omega_K);
        assert(1==0);
    }

    // Illegal ICFormat's will crash in loadIC.cpp; no need to crash here.
    /*
    ExpandPathName(DerivativesDirectory);
    ExpandPathName(ReadStateDirectory);
    ExpandPathName(WriteStateDirectory);
    ExpandPathName(PastStateDirectory);
    ExpandPathName(OutputDirectory);
    ExpandPathName(MultipoleDirectory);
    ExpandPathName(LogDirectory);
    ExpandPathName(InitialConditionsDirectory);

    CheckDirectoryExists(DerivativesDirectory);
    CheckDirectoryExists(ReadStateDirectory);
    CheckDirectoryExists(WriteStateDirectory);
    CheckDirectoryExists(PastStateDirectory);
    CheckDirectoryExists(OutputDirectory);
    CheckDirectoryExists(LogDirectory);

    CheckDirectoryExists(InitialConditionsDirectory);
    */


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

#endif
