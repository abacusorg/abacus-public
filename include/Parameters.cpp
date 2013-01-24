#define FSYM "f"
#define ESYM "e"
#define ISYM "d"
#define MAX_LINE_LENGTH 1024

#define QUOTEME(X) #X
#define SCANLINE(X, XSYM)   ret += sscanf(line, QUOTEME(X  =  %XSYM),   &X);
#define STRSCANLINE(X)      ret += sscanf(line, QUOTEME(X  =  %s), X);
#define DUMPVAR(X, XSYM)    fprintf(stderr, QUOTEME(X  =  %XSYM\n), X);

#define TESTSTRINGUNDEFINED(X) if(strcmp(X,"NotDefined")==0) { fprintf(stderr,QUOTEME(You didnt define X\n)); assert(1==0); }

#define TestVariablePresent(variable) if(variable<0) { \
     printf("Can't find paramter <"#variable"> in the parameter file\n");  \
     assert(1==0); }

class Parameters {
public:
    
    int np;
    int cpd;
    int order;

    int NearFieldRadius;    // Radius of cells in the near-field
    float SofteningLength; // Softening length in units of interparticle spacing

    int  DerivativeExpansionRadius;
    int  MAXConvolutionRAMMB;
    int  ConvolutionCacheSizeMB;

    int  DirectNewtonRaphson;  // 0 or 1 
    int  DirectDoublePrecision; // 0 or 1 

    char DerivativesDirectory[1024];

    char InitialConditionsFile[1024];   // The initial condition file name

    char ReadStateDirectory[1024];  // Where the input State lives
    char WriteStateDirectory[1024]; // Where the output State lives
    char PastStateDirectory[1024];  // Where the old input State lives
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


    Parameters() {

        np                                  = -1;
        cpd                                 = -1;
        order                               = -1;

        NearFieldRadius                     = -1;
        SofteningLength                     = -1.0;

        DirectDoublePrecision               = -1;
        DirectNewtonRaphson                 = -1;

        DerivativeExpansionRadius           = -1;
        MAXConvolutionRAMMB                 = -1;
        ConvolutionCacheSizeMB              = -1;
        sprintf(DerivativesDirectory,"NotDefined");

        sprintf(InitialConditionsFile,"NotDefined");

        sprintf(ReadStateDirectory,"NotDefined");
        sprintf(WriteStateDirectory,"NotDefined");
        sprintf(PastStateDirectory,"NotDefined");
        sprintf(LogFileDirectory,"NotDefined");
        sprintf(OutputDirectory,"NotDefined");
        sprintf(BaseDistributionDirectory,"NotDefined");

        sprintf(DumpFilePrefix,"NotDefined");
        sprintf(GroupFilePrefix,"NotDefined");
        sprintf(LightFilePrefix,"NotDefined");
    
        nDumpz                              = -1;
    
        H0                                  = -1.0;
        Omega_M                             = -1.0;
        Omega_DE                            = -1.0;
        Omega_K                             = -1.0;
        w0                                  = -1.0;
        wa                                  = -1.0;

        BoxSize                             = -1.0;
        hMpc                                = -1;
        InitialRedshift                     = -1.0;
        LagrangianPTOrder                   = -1.0;

        GroupRadius                         = 1;
        Eta                                 = -1.0;
        Dlna                                = -1.0;

    }

    void ReadParameters(char *paramaterfile, int icflag);
    void DumpParameters(void);
    void CheckVariablesPresent(void);
    void ValidateParameters(void);

};

void Parameters::ReadParameters(char *parameterfile, int icflag) {
    char *rawline;
    char *line;

    rawline = new char[ MAX_LINE_LENGTH ];

    char *rawlinebase = rawline;

    FILE *fp;
    fp = fopen(parameterfile,"r");
    if(fp==NULL) { 
        printf("Cannot find parameterfile %s\n", parameterfile);
        assert(1==0);
    }
    assert(fp!=NULL);

    int ret=0;

    int lineno = 1;

    while (fgets(rawline, MAX_LINE_LENGTH, fp) != NULL) {
        ret = 0;

        line = rawline;
        int l = strlen(rawline);
        int i =0;
        while( (*rawline == ' ' || *rawline == '\t' || *rawline == '\n') && (i < l) ) {
            rawline++;
            i++;
        }
        line = rawline;
        
        SCANLINE(np,d);
        SCANLINE(cpd,d);
        SCANLINE(order,d);

        SCANLINE(NearFieldRadius,d);
        SCANLINE(SofteningLength,f);

        SCANLINE(DirectDoublePrecision,d);
        SCANLINE(DirectNewtonRaphson,d);

        SCANLINE(DerivativeExpansionRadius,d);
        SCANLINE(MAXConvolutionRAMMB,d);
        SCANLINE(ConvolutionCacheSizeMB,d);
        STRSCANLINE(DerivativesDirectory);

        STRSCANLINE(InitialConditionsFile);

        STRSCANLINE(ReadStateDirectory);
        STRSCANLINE(WriteStateDirectory);
        STRSCANLINE(PastStateDirectory);
        STRSCANLINE(LogFileDirectory);
        STRSCANLINE(OutputDirectory);
        STRSCANLINE(BaseDistributionDirectory);

        STRSCANLINE(DumpFilePrefix);
        STRSCANLINE(GroupFilePrefix);
        STRSCANLINE(LightFilePrefix);

        SCANLINE(nDumpz,d);
        SCANLINE(H0,f);
        SCANLINE(Omega_M,f);
        SCANLINE(Omega_DE,f);
        SCANLINE(Omega_K,f);
        SCANLINE(w0,f);
        SCANLINE(wa,f);

        SCANLINE(BoxSize,f);
        SCANLINE(hMpc,d);
        SCANLINE(InitialRedshift,f);
        SCANLINE(LagrangianPTOrder,d);

        SCANLINE(GroupRadius,d);
        SCANLINE(Eta,f);
        SCANLINE(Dlna,f);

        int allwhitespace=0;
        int whitespace = 0;
        int sl = strlen(line);
        for(int i=0;i<sl;i++) whitespace += (line[i] == ' ');
        if(whitespace==sl-1) allwhitespace = 1;

        char firstnonspacecharacter;
        for(i=0;i<sl;i++) if(line[i]!=' ') break;
        firstnonspacecharacter = line[i];

        if( (ret==0) && !allwhitespace ) {
            if( firstnonspacecharacter == '#' ) { }
            else {
                fprintf(stderr, 
                   "[ERROR] line<#%d>:%s   was NOT interpreted\n", lineno, line);
                assert(1==0);
            }
        }

        lineno++;
    }

    delete[] rawlinebase;
    fclose(fp);

    CheckVariablesPresent();
    if(!icflag) ValidateParameters();
}

void Parameters::CheckVariablesPresent(void) {

    TestVariablePresent( np                           );
    TestVariablePresent( cpd                          );
    TestVariablePresent( order                        );

    TestVariablePresent( NearFieldRadius              );
    TestVariablePresent( SofteningLength              );

    TestVariablePresent( DirectDoublePrecision        );
    TestVariablePresent( DirectNewtonRaphson          );

    TestVariablePresent( DerivativeExpansionRadius    );
    TestVariablePresent( MAXConvolutionRAMMB          );
    TestVariablePresent( ConvolutionCacheSizeMB       );

    TestVariablePresent( H0                           );
    TestVariablePresent( Omega_M                      );
    TestVariablePresent( Omega_DE                     );
    TestVariablePresent( Omega_K                      );
    TestVariablePresent( w0                           );
    TestVariablePresent( wa                           );

    TestVariablePresent( BoxSize                      );
    TestVariablePresent( hMpc                         );
    TestVariablePresent( InitialRedshift              );
    TestVariablePresent( LagrangianPTOrder            );
    TestVariablePresent( GroupRadius                  );
    TestVariablePresent( Eta                          );
    TestVariablePresent( Dlna                         );

    TESTSTRINGUNDEFINED(DerivativesDirectory);
    TESTSTRINGUNDEFINED(InitialConditionsFile);

    TESTSTRINGUNDEFINED(ReadStateDirectory);
    TESTSTRINGUNDEFINED(WriteStateDirectory);
    TESTSTRINGUNDEFINED(PastStateDirectory);
    TESTSTRINGUNDEFINED(LogFileDirectory);
    TESTSTRINGUNDEFINED(OutputDirectory);
    TESTSTRINGUNDEFINED(BaseDistributionDirectory);

    TESTSTRINGUNDEFINED(DumpFilePrefix);
    TESTSTRINGUNDEFINED(GroupFilePrefix);
    TESTSTRINGUNDEFINED(LightFilePrefix);

}

void Parameters::ValidateParameters(void) {

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
        printf("DirectNewtonRapson must be 0 or 1\n");
        assert(1==0);
    }

    if( !( (DirectDoublePrecision == 1) || (DirectDoublePrecision==0) ) ) {
        printf("DirectDoublePrecision must be 0 or 1 \n");
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

    if(ConvolutionCacheSizeMB<0) {
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

    ExpandPathName(DerivativesDirectory);
    ExpandPathName(ReadStateDirectory);
    ExpandPathName(WriteStateDirectory);
    ExpandPathName(PastStateDirectory);
    ExpandPathName(OutputDirectory);
    ExpandPathName(LogFileDirectory);
    ExpandPathName(InitialConditionsFile);

    CheckDirectoryExists(DerivativesDirectory);
    CheckDirectoryExists(ReadStateDirectory);
    CheckDirectoryExists(WriteStateDirectory);
    CheckDirectoryExists(PastStateDirectory);
    CheckDirectoryExists(OutputDirectory);
    CheckDirectoryExists(LogFileDirectory);
    CheckDirectoryExists(BaseDistributionDirectory);

    CheckFileExists(InitialConditionsFile);


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


}


void Parameters::DumpParameters(void) {

    DUMPVAR( np                           , d);
    DUMPVAR( cpd                          , d);
    DUMPVAR( order                        , d);

    DUMPVAR( NearFieldRadius              , d);
    DUMPVAR( SofteningLength              , e);

    DUMPVAR( DirectNewtonRaphson          , d);
    DUMPVAR( DirectDoublePrecision        , d);

    DUMPVAR( DerivativeExpansionRadius    , d);
    DUMPVAR( MAXConvolutionRAMMB          , d);
    DUMPVAR( ConvolutionCacheSizeMB       , d);
    DUMPVAR( DerivativesDirectory         , s);

    DUMPVAR( InitialConditionsFile        , s);

    DUMPVAR( ReadStateDirectory           , s);
    DUMPVAR( WriteStateDirectory          , s);
    DUMPVAR( PastStateDirectory           , s);
    DUMPVAR( LogFileDirectory             , s);
    DUMPVAR( OutputDirectory              , s);
    DUMPVAR( BaseDistributionDirectory    , s);

    DUMPVAR( DumpFilePrefix               , s);
    DUMPVAR( GroupFilePrefix              , s);
    DUMPVAR( LightFilePrefix              , s);

    DUMPVAR( nDumpz                       , d);
    
    DUMPVAR( H0                           , e);
    DUMPVAR( Omega_M                      , e);
    DUMPVAR( Omega_DE                     , e);
    DUMPVAR( Omega_K                      , e);
    DUMPVAR( w0                           , e);
    DUMPVAR( wa                           , e);

    DUMPVAR( BoxSize                      , e);
    DUMPVAR( hMpc                         , d);
    DUMPVAR( InitialRedshift              , e);
    DUMPVAR( LagrangianPTOrder            , d);
    DUMPVAR( GroupRadius                  , d);
    DUMPVAR( Eta                          , e);
    DUMPVAR( Dlna                         , e);

    
}

Parameters P;
