#include "header.cpp"
#include "threevector.hh"


#include "STimer.h"
#include "PTimer.h"

STimer TotalWallClock;
STimer Setup;
STimer ConvolutionWallClock;

#include "file.h"
#include "factorial.h"

#include "stdlog.cc"

#include "Parameters.cpp"

#include "OutofCoreConvolution.h"

#include <fenv.h>

void dumpstats(OutofCoreConvolution *OCC, char *fn) {

    FILE *fp;
    fp = fopen(fn,"w");
    assert(fp!=NULL);


    double accountedtime  = OCC->CS.ConvolutionArithmetic;
           accountedtime += OCC->CS.ForwardZFFTMultipoles + OCC->CS.InverseZFFTTaylor;
           accountedtime += OCC->CS.ReadDerivatives + OCC->CS.ReadMultipoles + OCC->CS.WriteTaylor;
           accountedtime += OCC->CS.ArraySwizzle;
    double discrepency = OCC->CS.ConvolveWallClock - accountedtime;

    int computecores = omp_get_num_procs()/2;
    fprintf(fp,"Convolution parameters:  RamAllocated = %dMB CacheSizeMB = %dMB nreal_cores=%d blocksize=%d zwidth=%d cpd=%d order=%d\n\n",
        (int) (OCC->CS.totalMemoryAllocated/(1<<20)), OCC->CS.runtime_ConvolutionCacheSizeMB, computecores, OCC->CS.blocksize, OCC->CS.zwidth, OCC->CP.runtime_cpd, OCC->CP.runtime_order);

    fprintf(fp,"\t ConvolutionWallClock:  %2.2e seconds \n", OCC->CS.ConvolveWallClock );
    fprintf(fp,"\t \t \t %50s : %1.1e seconds\n", "Array Swizzling", OCC->CS.ArraySwizzle );
    double e = OCC->CS.ReadDerivativesBytes/OCC->CS.ReadDerivatives/(1.0e+6);
    fprintf(fp,"\t \t \t %50s : %1.1e seconds --> rate was %4d MB/s\n", "ReadDiskDerivatives", OCC->CS.ReadDerivatives, (int) (e) );
    e = OCC->CS.ReadMultipolesBytes/OCC->CS.ReadMultipoles/(1.0e+6);
    fprintf(fp,"\t \t \t %50s : %1.1e seconds --> rate was %4d MB/s\n", "ReadDiskMultipoles", OCC->CS.ReadMultipoles, (int) (e) );
    e = OCC->CS.WriteTaylorBytes/OCC->CS.WriteTaylor/(1.0e+6);
    fprintf(fp,"\t \t \t %50s : %1.1e seconds --> rate was %4d MB/s\n", "WriteDiskTaylor", OCC->CS.WriteTaylor, (int) (e) );
    double Gops = ((double) OCC->CS.ops)/(1.0e+9);
    fprintf(fp,"\t \t \t %50s : %1.1e seconds for %5.3f billion double precision operations\n", "Convolution Arithmetic", OCC->CS.ConvolutionArithmetic, Gops );
    fprintf(fp,"\t \t \t %50s : %1.1e seconds\n", "Forward FFT Z Multipoles", OCC->CS.ForwardZFFTMultipoles );
    fprintf(fp,"\t \t \t %50s : %1.1e seconds\n", "Inverse FFT Z Taylor",         OCC->CS.InverseZFFTTaylor );
    fprintf(fp,"\t \t %50s : %1.1e seconds which is %d%% \n", "Unaccounted remaining wallclock time", discrepency, (int) (discrepency/OCC->CS.ConvolveWallClock*100) );

    double cae = OCC->CS.ConvolutionArithmetic;
    int farithp   = cae/OCC->CS.ConvolveWallClock*100;
    int ffftp     = (OCC->CS.ForwardZFFTMultipoles + OCC->CS.InverseZFFTTaylor)/OCC->CS.ConvolveWallClock*100;
    int fiop      = (OCC->CS.ReadDerivatives + OCC->CS.ReadMultipoles + OCC->CS.WriteTaylor )/OCC->CS.ConvolveWallClock*100;
    int swzp      = OCC->CS.ArraySwizzle/OCC->CS.ConvolveWallClock*100;

    fprintf(fp,"\t Summary: Fourier Transforms = %d%%     Convolution Arithmetic = %d%%    DiskIO = %d%%   Array Swizzle = %d%% \n", ffftp, farithp, fiop, swzp );
    fprintf(fp,"\t          Arithmetic rate = %d DGOPS --> rate per core = %1.1f DGOPS\n", (int) (Gops/cae), Gops/cae/computecores );
    fprintf(fp,"\t          [DGOPS == Double Precision Billion operations per second]\n");
    fprintf(fp,"\n");

    fclose(fp);
}


int main(int argc, char ** argv){
	TotalWallClock.Start();
	Setup.Start();

	if (argc!=2) {
	       // Can't use assertf() or QUIT here: stdlog not yet defined!
	       fprintf(stderr, "singlestep(): command line must have 2 parameters given, not %d.\nLegal usage: convolve <parameter_file>\n", argc);
	       assert(0==99);
	    }

	    P.ReadParameters(argv[1],1);

	    // Setup the log
	    stdlog_threshold_global = P.LogVerbosity;
	    char logfn[1050];
	    sprintf(logfn,"%s/lastconvolution.log", P.LogDirectory);
	    stdlog.open(logfn);
	    OutofCoreConvolution OCC;
	    STDLOG(1,"Read parameter file\n");
	    Setup.Stop();

	    ConvolutionParameters p;
	    p.runtime_ConvolutionCacheSizeMB = P.ConvolutionCacheSizeMB;
	    p.runtime_DerivativeExpansionRadius = P.DerivativeExpansionRadius;
	    strcpy(p.runtime_DerivativesDirectory,P.DerivativesDirectory);
	    p.runtime_DiskBufferSizeKB = 1LL<<23;
	    p.runtime_IsRamDisk = P.RamDisk;
	    p.runtime_MaxConvolutionRAMMB = P.MAXRAMMB;
	    strcpy(p.runtime_MultipoleDirectory,P.ReadStateDirectory);
	    sprintf(p.runtime_MultipolePrefix, "Multipoles");
	    p.runtime_NearFieldRadius = P.NearFieldRadius;
	    strcpy(p.runtime_TaylorDirectory,P.WriteStateDirectory);
	    p.runtime_cpd = P.cpd;
	    p.runtime_order = P.order;
	    sprintf(p.runtime_TaylorPrefix, "Taylor");

	    ConvolutionWallClock.Start();
	    STDLOG(1,"Starting Convolution\n");
	    OCC.Convolve(p);
	    STDLOG(1,"Convolution Complete\n");
	    ConvolutionWallClock.Stop();
	    TotalWallClock.Stop();
	    char timingfn[1050];
	    sprintf(timingfn,"%s/lastconvolution.timing",P.LogDirectory);
	    dumpstats(&OCC,timingfn);
	    stdlog.close();
	        exit(0);



 
}
