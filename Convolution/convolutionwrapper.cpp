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


    double accountedtime  = OCC->CS.ConvolutionArithmetic.Elapsed();
           accountedtime += OCC->CS.ForwardZFFTMultipoles.Elapsed() + OCC->CS.InverseZFFTTaylor.Elapsed();
           accountedtime += OCC->CS.ReadDerivatives.Elapsed() + OCC->CS.ReadMultipoles.Elapsed() + OCC->CS.WriteTaylor.Elapsed();
           accountedtime += OCC->CS.ArraySwizzle.Elapsed();
    double discrepency = OCC->CS.ConvolveWallClock.Elapsed() - accountedtime;

    int computecores = omp_get_num_procs()/2;
    fprintf(fp,"Convolution parameters:  RamAllocated = %dMB CacheSizeMB = %dMB nreal_cores=%d blocksize=%d zwidth=%d cpd=%d order=%d\n\n",
        (int) (OCC->CS.totalMemoryAllocated/(1<<20)), OCC->CS.runtime_ConvolutionCacheSizeMB, computecores, OCC->CS.blocksize, OCC->CS.zwidth, OCC->CP.runtime_cpd, OCC->CP.runtime_order);

    fprintf(fp,"\t ConvolutionWallClock:  %2.2e seconds \n", OCC->CS.ConvolveWallClock.Elapsed() );
    fprintf(fp,"\t \t \t %50s : %1.1e seconds\n", "Array Swizzling", OCC->CS.ArraySwizzle.Elapsed() );
    double e = OCC->CS.ReadDerivativesBytes/OCC->CS.ReadDerivatives.Elapsed()/(1.0e+6);
    fprintf(fp,"\t \t \t %50s : %1.1e seconds --> rate was %4d MB/s\n", "ReadDiskDerivatives", OCC->CS.ReadDerivatives.Elapsed(), (int) (e) );
    e = OCC->CS.ReadMultipolesBytes/OCC->CS.ReadMultipoles.Elapsed()/(1.0e+6);
    fprintf(fp,"\t \t \t %50s : %1.1e seconds --> rate was %4d MB/s\n", "ReadDiskMultipoles", OCC->CS.ReadMultipoles.Elapsed(), (int) (e) );
    e = OCC->CS.WriteTaylorBytes/OCC->CS.WriteTaylor.Elapsed()/(1.0e+6);
    fprintf(fp,"\t \t \t %50s : %1.1e seconds --> rate was %4d MB/s\n", "WriteDiskTaylor", OCC->CS.WriteTaylor.Elapsed(), (int) (e) );
    double Gops = ((double) OCC->CS.ops)/(1.0e+9);
    fprintf(fp,"\t \t \t %50s : %1.1e seconds for %5.3f billion double precision operations\n", "Convolution Arithmetic", OCC->CS.ConvolutionArithmetic.Elapsed(), Gops );
    fprintf(fp,"\t \t \t %50s : %1.1e seconds\n", "Forward FFT Z Multipoles", OCC->CS.ForwardZFFTMultipoles.Elapsed() );
    fprintf(fp,"\t \t \t %50s : %1.1e seconds\n", "Inverse FFT Z Taylor",         OCC->CS.InverseZFFTTaylor.Elapsed() );
    fprintf(fp,"\t \t %50s : %1.1e seconds which is %d%% \n", "Unaccounted remaining wallclock time", discrepency, (int) (discrepency/OCC->CS.ConvolveWallClock.Elapsed()*100) );

    double cae = OCC->CS.ConvolutionArithmetic.Elapsed();
    int farithp   = cae/OCC->CS.ConvolveWallClock.Elapsed()*100;
    int ffftp     = (OCC->CS.ForwardZFFTMultipoles.Elapsed() + OCC->CS.InverseZFFTTaylor.Elapsed())/OCC->CS.ConvolveWallClock.Elapsed()*100;
    int fiop      = (OCC->CS.ReadDerivatives.Elapsed() + OCC->CS.ReadMultipoles.Elapsed() + OCC->CS.WriteTaylor.Elapsed() )/OCC->CS.ConvolveWallClock.Elapsed()*100;
    int swzp      = OCC->CS.ArraySwizzle.Elapsed()/OCC->CS.ConvolveWallClock.Elapsed()*100;

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
	    STDLOG(1,"Read parameter file\n");
	    Setup.Stop();

	    ConvolutionParameters p;
	    p.runtime_ConvolutionCacheSizeMB = P.ConvolutionCacheSizeMB;
	    p.runtime_DerivativeExpansionRadius = P.DerivativeExpansionRadius;
	    strcpy(P.DerivativesDirectory,p.runtime_DerivativesDirectory);
	    p.runtime_DiskBufferSizeKB = 4;
	    p.runtime_IsRamDisk = P.RamDisk;
	    p.runtime_MaxConvolutionRAMMB = P.MAXConvolutionRAMMB;
	    strcpy(P.ReadStateDirectory,p.runtime_MultipoleDirectory);
	    p.runtime_MultipolePrefix = "Multipoles";
	    p.runtime_NearFieldRadius = P.NearFieldRadius;
	    strcpy(P.WriteStateDirectory,p.runtime_TaylorDirectory);
	    p.runtime_cpd = P.cpd;
	    p.runtime_order = P.order;


	    ConvolutionWallClock.Start();
	    STDLOG(1,"Starting Convolution\n");
	    OCC.Convolve(p);
	    STDLOG(1,"Convolution Complete\n");
	    ConvolutionWallClock.Stop();
	    TotalWallClock.Stop();
	    char timingfn[1050];
	    sprintf(timingfn,"%s/lastconvolution.timing",P.LogFileDirectory);
	    dumpstats(&OCC,timingfn);
	    stdlog.close();
	        exit(0);




}
