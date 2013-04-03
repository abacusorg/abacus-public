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

    double accountedtime  = OCC->ConvolutionArithmetic.Elapsed();
           accountedtime += OCC->ForwardZFFTMultipoles.Elapsed() + OCC->InverseZFFTTaylor.Elapsed();
           accountedtime += OCC->ReadDerivatives.Elapsed() + OCC->ReadMultipoles.Elapsed() + OCC->WriteTaylor.Elapsed();
           accountedtime += OCC->ArraySwizzle.Elapsed();
    double discrepency = OCC->ConvolveWallClock.Elapsed() - accountedtime;

    int computecores = omp_get_num_procs()/2;
    fprintf(fp,"Convolution parameters:  RamAllocated = %dMB CacheSizeMB = %dMB nreal_cores=%d blocksize=%d zwidth=%d cpd=%d order=%d\n\n",
        (int) (OCC->totalMemoryAllocated/(1<<20)), OCC->runtime_ConvolutionCacheSizeMB, computecores, OCC->blocksize, OCC->zwidth, OCC->runtime_cpd, OCC->runtime_order);

    fprintf(fp,"\t ConvolutionWallClock:  %2.2e seconds \n", OCC->ConvolveWallClock.Elapsed() );
    fprintf(fp,"\t \t \t %50s : %1.1e seconds\n", "Array Swizzling", OCC->ArraySwizzle.Elapsed() );
    double e = OCC->ReadDerivativesBytes/OCC->ReadDerivatives.Elapsed()/(1.0e+6);
    fprintf(fp,"\t \t \t %50s : %1.1e seconds --> rate was %4d MB/s\n", "ReadDiskDerivatives", OCC->ReadDerivatives.Elapsed(), (int) (e) );
    e = OCC->ReadMultipolesBytes/OCC->ReadMultipoles.Elapsed()/(1.0e+6);
    fprintf(fp,"\t \t \t %50s : %1.1e seconds --> rate was %4d MB/s\n", "ReadDiskMultipoles", OCC->ReadMultipoles.Elapsed(), (int) (e) );
    e = OCC->WriteTaylorBytes/OCC->WriteTaylor.Elapsed()/(1.0e+6);
    fprintf(fp,"\t \t \t %50s : %1.1e seconds --> rate was %4d MB/s\n", "WriteDiskTaylor", OCC->WriteTaylor.Elapsed(), (int) (e) );
    double Gops = ((double) OCC->ops)/(1.0e+9);
    fprintf(fp,"\t \t \t %50s : %1.1e seconds for %5.3f billion double precision operations\n", "Convolution Arithmetic", OCC->ConvolutionArithmetic.Elapsed(), Gops );
    fprintf(fp,"\t \t \t %50s : %1.1e seconds\n", "Forward FFT Z Multipoles", OCC->ForwardZFFTMultipoles.Elapsed() );
    fprintf(fp,"\t \t \t %50s : %1.1e seconds\n", "Inverse FFT Z Taylor",         OCC->InverseZFFTTaylor.Elapsed() );
    fprintf(fp,"\t \t %50s : %1.1e seconds which is %d%% \n", "Unaccounted remaining wallclock time", discrepency, (int) (discrepency/OCC->ConvolveWallClock.Elapsed()*100) );

    double cae = OCC->ConvolutionArithmetic.Elapsed();
    int farithp   = cae/OCC->ConvolveWallClock.Elapsed()*100;
    int ffftp     = (OCC->ForwardZFFTMultipoles.Elapsed() + OCC->InverseZFFTTaylor.Elapsed())/OCC->ConvolveWallClock.Elapsed()*100;
    int fiop      = (OCC->ReadDerivatives.Elapsed() + OCC->ReadMultipoles.Elapsed() + OCC->WriteTaylor.Elapsed() )/OCC->ConvolveWallClock.Elapsed()*100;
    int swzp      = OCC->ArraySwizzle.Elapsed()/OCC->ConvolveWallClock.Elapsed()*100;

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
	    sprintf(logfn,"%s/lastconvolution.log", P.LogFileDirectory);
	    stdlog.open(logfn);
	    OutofCoreConvolution OCC;
	    STDLOG(1,"Read parameter file\n")
	    Setup.Stop();
	    ConvolutionWallClock.Start();
	    STDLOG(1,"Starting Convolution\n");
	    OCC.Convolve( P.cpd, P.order, P.NearFieldRadius, P.DerivativeExpansionRadius,
	                      P.RamDisk, 1024,  P.ConvolutionCacheSizeMB, P.MAXConvolutionRAMMB,
	                      P.DerivativesDirectory, P.ReadStateDirectory, P.WriteStateDirectory,
	                      "Multipoles", "Taylor" );
	    STDLOG(1,"Convolution Complete\n")
	    ConvolutionWallClock.Stop();
	    TotalWallClock.Stop();
	    char timingfn[1050];
	    sprintf(timingfn,"%s/lastconvolution.timing",P.LogFileDirectory);
	    dumpstats(&OCC,timingfn);
	    stdlog.close();
	        exit(0);




}
