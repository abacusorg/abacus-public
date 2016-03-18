/* reporting.cpp

This writes timing outputs.  It should provide a nicely formatted
summary of the timings from the various portions of the code.


TODO: Need to uncomment the Multipole and Taylor portions.

*/

int DOUBLECLOSE(double a, double b, double eps) {
    return ( fabs(a-b) < eps*fabs(a+b) );
}

void DumpMultipoleStats(double totalwallclock,int np) {
#ifdef NOT_READY_YET
    double a, g, q, f1, f2;

    g = 100.0/ComputeMultipoles.Elapsed();
    a = MF->FFTMultipole.Elapsed();         printf("\t \t \t \t \t        MultipoleFFT: %2.2e seconds (%2.2f %%)\n", a, a*g);
    q = ComputeMultipoles.Elapsed() - MF->FFTMultipole.Elapsed();
    f1 = MF->MultipoleASM.Elapsed()/(MF->MultipoleASM.Elapsed() + MF->MultipoleC2R.Elapsed() +1e-15);
    f2 = 1.0-f1;
    printf("\t \t \t \t \t        MultipoleASM: %2.2e seconds (%2.2f %%) ---> %e pps \n", f1*q, f1*q*g, np/(f1*q) );
    printf("\t \t \t \t \t        MultipoleC2R: %2.2e seconds (%2.2f %%)\n", f2*q, f2*q*g);

#endif
}

void DumpTaylorStats(double totalwallclock,int np) {

#ifdef NOT_READY_YET
    double a, g, t, f1, f2, f3, z;

    g = 100.0/TaylorForce.Elapsed();
    a = TY->FFTTaylor.Elapsed();        printf("\t \t \t \t      FFTTaylor: %e seconds (%2.2f %%)\n", a, a*g );

    t = TY->TaylorRedlack.Elapsed() + TY->TaylorR2C.Elapsed() + TY->TaylorASM.Elapsed()+1e-15;
    f1 = TY->TaylorRedlack.Elapsed()/t;
    f2 = TY->TaylorR2C.Elapsed()/t;
    f3 = TY->TaylorASM.Elapsed()/t;

    z = TaylorForce.Elapsed() - TY->FFTTaylor.Elapsed()+1e-15;

     printf("\t \t \t \t      TaylorASM: %2.2e seconds (%2.2f %%) ---> %e pps \n", f3*z, f3*z*g, np/(f3*z) );
     printf("\t \t \t \t  TaylorR2C: %2.2e seconds (%2.2f %%)\n", f2*z, f2*z*g );
     printf("\t \t \t \t  TaylorRedlack: %2.2e seconds (%2.2f %%)\n", f1*z, f1*z*g);
 #endif
}

#ifdef DONOTCOMPILE
void DumpFinishStats(void) { 

    double a,g;
    g = 100.0/(Finish.Elapsed()+1e-15);
    a = MemCpy.Elapsed();               printf("\t \t \t \t                 Memcpy: %2.2e seconds (%2.2f %%)\n", a, a*g );
    a = WriteMergeSlab.Elapsed();       printf("\t \t \t \t         WriteMergeSlab: %2.2e seconds (%2.2f %%)\n", a, a*g );
    a = WriteMultipoleSlab.Elapsed();   printf("\t \t \t \t     WriteMultipoleSlab: %2.2e seconds (%2.2f %%)\n", a, a*g );
    a = WriteCellInfo.Elapsed();        printf("\t \t \t \t          WriteCellInfo: %2.2e seconds (%2.2f %%)\n", a, a*g );
    a = FillCellInfo.Elapsed();         printf("\t \t \t \t           FillCellInfo: %2.2e seconds (%2.2f %%)\n", a, a*g );
}



void ReportTimings(FILE * timingfile) {
    double a,f,g,q;

    double w = FetchSlabs.Elapsed() + FetchCellInfo.Elapsed() + NearForce.Elapsed() + TaylorFetch.Elapsed() +
                    TaylorForce.Elapsed() + Kick.Elapsed() + Group.Elapsed() + Output.Elapsed() + Drift.Elapsed() + Finish.Elapsed();

    assert( DOUBLECLOSE(w,TimeStepWallClock.Elapsed(),0.05) );

    f = 100.0/(WallClockDirect.Elapsed()+1e-15);

    fprintf(timingfile,"DirectWallClock:  %2.2e seconds\n", WallClockDirect.Elapsed() );
    fprintf(timingfile,"\t \t       Prologue : %2.2e seconds (%2.2f %%)\n", prologue.Elapsed(), prologue.Elapsed() * f);
    fprintf(timingfile,"\t \t       Epilogue : %2.2e seconds (%2.2f %%)\n", epilogue.Elapsed(), epilogue.Elapsed() * f);
    fprintf(timingfile,"\t \t       TimeStep : %2.2e seconds (%2.2f %%)\n", w, w*f);
    q = WallClockDirect.Elapsed()-epilogue.Elapsed()-prologue.Elapsed()-w ;
    fprintf(timingfile,"\t \t    Discrepancy : %2.2e seconds (%2.2f %%)\n", q, q*f );
    fprintf(timingfile,"\n");


    a = FetchSlabs.Elapsed();       fprintf(timingfile,"\t \t \t    FetchSlabs: %2.2e seconds (%2.2f %%)\n", a, a*f );
    a = FetchCellInfo.Elapsed();    fprintf(timingfile,"\t \t \t FetchCellInfo: %2.2e seconds (%2.2f %%)\n", a, a*f );
    a = NearForce.Elapsed();        fprintf(timingfile,"\t \t \t     NearForce: %2.2e seconds (%2.2f %%) ---> %e DIPS \n", a, a*f,  P.np/(NearForce.Elapsed()+1e-15) );
    a = TaylorFetch.Elapsed();      fprintf(timingfile,"\t \t \t   TaylorFetch: %2.2e seconds (%2.2f %%)\n", a, a*f );
    a = TaylorForce.Elapsed();     fprintf(timingfile,"\t \t \t   TaylorForce: %2.2e seconds (%2.2f %%) ---> %e pps \n", a, a*f, P.np/(a+1e-15) );


    DumpTaylorStats( WallClockDirect.Elapsed() , P.np); 

    a = Kick.Elapsed();             fprintf(timingfile,"\t \t \t          Kick: %2.2e seconds (%2.2f %%)\n", a, a*f );
    a = Group.Elapsed();            fprintf(timingfile,"\t \t \t         Group: %2.2e seconds (%2.2f %%)\n", a, a*f );
    a = Output.Elapsed();           fprintf(timingfile,"\t \t \t        Output: %2.2e seconds (%2.2f %%)\n", a, a*f );
    a = Drift.Elapsed();            fprintf(timingfile,"\t \t \t         Drift: %2.2e seconds (%2.2f %%)\n", a, a*f );

    fprintf(timingfile,"\n");
    a = Finish.Elapsed();               fprintf(timingfile,"\t \t \t        Finish: %2.2e seconds (%2.2f %%)\n", a, a*f );
    g = 100.0/(Finish.Elapsed()+1e-15);

    a = MemCpy.Elapsed();               fprintf(timingfile,"\t \t \t \t                 Memcpy: %e seconds (%2.2f %%)\n", a, a*g );
    a = WriteMergeSlab.Elapsed();       fprintf(timingfile,"\t \t \t \t         WriteMergeSlab: %e seconds (%2.2f %%)\n", a, a*g );
    a = WriteMultipoleSlab.Elapsed();   fprintf(timingfile,"\t \t \t \t     WriteMultipoleSlab: %e seconds (%2.2f %%)\n", a, a*g );
    a = WriteCellInfo.Elapsed();        fprintf(timingfile,"\t \t \t \t          WriteCellInfo: %e seconds (%2.2f %%)\n", a, a*g );
    a = FillCellInfo.Elapsed();         fprintf(timingfile,"\t \t \t \t           FillCellInfo: %e seconds (%2.2f %%)\n", a, a*g );

    a = ComputeMultipoles.Elapsed();    fprintf(timingfile,"\t \t \t \t      ComputeMultipoles: %e seconds (%2.2f %%) ---> %e pps \n", a, a*g, P.np/(a+1e-15) );
    DumpMultipoleStats(WallClockDirect.Elapsed(), P.np);

    a = FinishSort.Elapsed();           fprintf(timingfile,"\t \t \t \t             FinishSort: %e seconds (%2.2f %%)\n", a, a*g );
    a = FinishPartition.Elapsed();      fprintf(timingfile,"\t \t \t \t        FinishPartition: %e seconds (%2.2f %%)\n", a, a*g );

}
#endif


#define REPORT(tabs, str, a) \
     do { \
	 fprintf(timingfile, "\n"); \
	 for (int i=0; i<tabs; i++) fprintf(timingfile, "    "); \
	 thistime = a; fprintf(timingfile, "%-30s: %10.2e sec (%5.2f%%) ", str, thistime, 100*thistime/(denom+1e-15)); \
     } while(0)

void ReportTimings(FILE * timingfile) {
    double thistime, denom, total;
    denom = WallClockDirect.Elapsed();
    // fprintf(timingfile,"DirectWallClock:  %2.2e seconds\n", WallClockDirect.Elapsed() );
    REPORT(0, "Total Wall Clock Time", WallClockDirect.Elapsed()); 
    fprintf(timingfile,"---> %6.3f Mpart/sec",P.np/(thistime+1e-15)/1e6);
    fprintf(timingfile,"\n");

    total = 0.0;
    REPORT(0, "Prologue", prologue.Elapsed()); total += thistime;
    // REPORT(0, "Epilogue", epilogue.Elapsed()); total += thistime;
    REPORT(0, "TimeStep", TimeStepWallClock.Elapsed()); total += thistime;
    REPORT(0, "SingleStep Setup", SingleStepSetup.Elapsed()); total += thistime;
    REPORT(0, "SingleStep TearDown", SingleStepTearDown.Elapsed()); total += thistime;
    REPORT(0, "Unaccounted", WallClockDirect.Elapsed()-total);

    fprintf(timingfile, "\n");
    REPORT(0, "Total Blocking Disk Reads", BlockingIOReadTime.Elapsed());
	fprintf(timingfile, "---> %6.2f MB/sec on %6.3f GB", blocking_read_bytes/(thistime+1e-15)/1e6, blocking_read_bytes/1e9);
    REPORT(0, "Total Blocking Disk Writes", BlockingIOWriteTime.Elapsed());
	fprintf(timingfile, "---> %6.2f MB/sec on %6.3f GB", blocking_write_bytes/(thistime+1e-15)/1e6, blocking_write_bytes/1e9);

	fprintf(timingfile, "\n");
	REPORT(0, "Total Non-blocking Disk Reads", NonBlockingIOReadTime.Elapsed());
	    fprintf(timingfile, "---> %6.2f MB/sec on %6.3f GB", non_blocking_read_bytes/(thistime+1e-15)/1e6, non_blocking_read_bytes/1e9);
	REPORT(0, "Total Non-blocking Disk Writes", NonBlockingIOWriteTime.Elapsed());
	    fprintf(timingfile, "---> %6.2f MB/sec on %6.3f GB", non_blocking_write_bytes/(thistime+1e-15)/1e6, non_blocking_write_bytes/1e9);

    fprintf(timingfile, "\n\nBreakdown of TimeStep: ");
    total = 0.0;
    REPORT(1, "FetchSlabs", FetchSlabs.Elapsed()); total += thistime;
#ifdef CUDADIRECT
    int NGPU = GetNGPU();
    REPORT(1, "NearForce [blocking]", NearForce.Elapsed()); total += thistime;
    fprintf(timingfile,"---> %6.3f Mpart/sec", P.np/(thistime+1e-15)/1e6 );
    REPORT(1, "NearForce [non-blocking]", GPUTimers[NGPU].Elapsed()); //total += thistime;
    double gdi_gpu = JJ->DirectInteractions_GPU/1e9;
    fprintf(timingfile,"---> %6.3f GDIPS, %6.3f Gdirects, %6.3f Mpart/sec", gdi_gpu/(thistime+1e-15), gdi_gpu, P.np/(thistime+1e-15)/1e6);
    double total_di = (JJ->DirectInteractions_CPU +JJ->DirectInteractions_GPU)/1e9;
	//fprintf(timingfile,"---> %6.3f GDIPS, %6.3f Gdirects, %6.3f Mpart/sec", total_di/(thistime+1e-15), total_di, P.np/(thistime+1e-15)/1e6);
    REPORT(1, "Spinning while waiting for GPU", WaitingForGPU.Elapsed()); total += thistime;
#else
    REPORT(1, "NearForce", NearForce.Elapsed()); total += thistime;
    fprintf(timingfile,"---> %6.3f Mpart/sec", P.np/(thistime+1e-15)/1e6 );
#endif
    REPORT(1, "TaylorForce", TaylorForce.Elapsed()); total += thistime;
	fprintf(timingfile,"---> %6.3f Mpart/sec", P.np/(thistime+1e-15)/1e6 );
    REPORT(1, "Kick", Kick.Elapsed()); total += thistime;
	fprintf(timingfile,"---> %6.3f Mpart/sec", P.np/(thistime+1e-15)/1e6 );
    REPORT(1, "Output", Output.Elapsed()); total += thistime;
    fprintf(timingfile,"---> %6.3f Mpart/sec",P.np/(thistime+1e-15)/1e6);
    REPORT(1, "Group", Group.Elapsed()); total += thistime;
    REPORT(1, "Drift", Drift.Elapsed()); total += thistime;
	fprintf(timingfile,"---> %6.3f Mpart/sec", P.np/(thistime+1e-15)/1e6 );
    REPORT(1, "Finish", Finish.Elapsed()); total += thistime;
	fprintf(timingfile,"---> %6.3f Mpart/sec", P.np/(thistime+1e-15)/1e6 );
    REPORT(1, "Unaccounted (spinning?)", TimeStepWallClock.Elapsed()-total);

    fprintf(timingfile, "\n\n Breakdown per slab (Wall Clock)");
    double slabforcetimemean = 0;
    double slabforcetimesigma = 0;
    double slabforcemaxtime = 0;
    double slabforcemintime =  1e9;
    double slabforcelatencymean = 0;
    double slabforcelatencysigma = 0;
    double slabforcemaxlatency = 0;
    double slabforceminlatency = 1e9;

    char fn[1024];
    sprintf(fn,"%s/lastrun.slabtimes",P.LogDirectory);
    FILE* slabtimefile = fopen(fn,"wb");
    for(int i =0; i < P.cpd;i++){
	double slabtime = SlabForceTime[i].Elapsed();
    	slabforcetimemean += SlabForceTime[i].Elapsed()/P.cpd;
    	slabforcetimesigma+= SlabForceTime[i].Elapsed()*SlabForceTime[i].Elapsed()/P.cpd;
    	if(SlabForceTime[i].Elapsed() > slabforcemaxtime) slabforcemaxtime = SlabForceTime[i].Elapsed();
    	if(SlabForceTime[i].Elapsed() < slabforcemintime) slabforcemintime = SlabForceTime[i].Elapsed();

    	double slablatency = SlabForceLatency[i].Elapsed() -SlabFarForceTime[i].Elapsed();
    	slabforcelatencymean += slablatency/P.cpd;
    	slabforcelatencysigma+= slablatency*slablatency/P.cpd;
    	if(slablatency > slabforcemaxlatency) slabforcemaxlatency = slablatency;
    	if(slablatency < slabforceminlatency) slabforceminlatency = slablatency;

	fwrite(&slabtime,sizeof(double),1,slabtimefile);
	fwrite(&slablatency,sizeof(double),1,slabtimefile);
    }
    fclose(slabtimefile);
    slabforcetimesigma =  sqrt(slabforcetimesigma    - slabforcetimemean*slabforcetimemean);
    slabforcelatencysigma=sqrt(slabforcelatencysigma - slabforcelatencymean*slabforcelatencymean);

    denom = WallClockDirect.Elapsed()/P.cpd;
    REPORT(1,"Mean Force Computation",slabforcetimemean);
    fprintf(timingfile,"\n Sigma(STD): %e s\t Min: %e s\t Max: %e s ",slabforcetimesigma,slabforcemintime,slabforcemaxtime);
    REPORT(1,"Mean Force Latency",slabforcelatencymean);
    fprintf(timingfile,"\n Sigma(STD): %e s\t Min: %e s\t Max: %e s ",slabforcelatencysigma,slabforceminlatency,slabforcemaxlatency);

    fprintf(timingfile, "\n");
    fprintf(timingfile, "\n\n Breakdown of Output Step:");
    denom = Output.Elapsed();
    REPORT(2, "TimeSlice", OutputTimeSlice.Elapsed());
    fprintf(timingfile,"---> %6.3f Mpart/sec",P.np/(thistime+1e-15)/1e6);
    REPORT(2, "LightCone", OutputLightCone.Elapsed());
    fprintf(timingfile,"---> %6.3f Mpart/sec",P.np/(thistime+1e-15)/1e6);
    REPORT(2, "Binning", OutputBin.Elapsed());
    fprintf(timingfile,"---> %6.3f Mpart/sec",P.np/(thistime+1e-15)/1e6);

    fprintf(timingfile, "\n\n Breakdown of Finish Step:");
    denom = Finish.Elapsed();
    REPORT(2, "Partition Insert List", FinishPartition.Elapsed());
    REPORT(2, "Sort Insert List", FinishSort.Elapsed());
    fprintf(timingfile,"---> %6.3f Mitem/sec (%lld items)",IL->TotalLength/(thistime+1e-15)/1e6, IL->TotalLength);
    REPORT(2, "Index Cells", FinishCellIndex.Elapsed());
    REPORT(2, "Merge", FinishMerge.Elapsed());
    REPORT(2, "Compute Multipoles", ComputeMultipoles.Elapsed());
	fprintf(timingfile,"---> %6.3f Mpart/sec", P.np/(thistime+1e-15)/1e6 );
    REPORT(2, "Write Particles", WriteMergeSlab.Elapsed());
    REPORT(2, "Write Multipoles", WriteMultipoleSlab.Elapsed());

    fprintf(timingfile, "\n\n Subdivisions of Drift:");
    denom = Drift.Elapsed();
    REPORT(2, "Move & Rebin (S)", DriftMoveRebin.Elapsed());
    REPORT(2, "Inserting (S)",    DriftInsert.Elapsed());
    REPORT(2, "LPT IC Re-reading (S)",   LPTDriftICReRead.Elapsed());
    REPORT(2, "Move (P)",         DriftMove.Elapsed());
    REPORT(2, "Rebin (P)",        DriftRebin.Elapsed());
    
    fprintf(timingfile, "\n\n Subdivisions of Near Force:");
#ifdef CUDADIRECT
    denom = NearForce.Elapsed();
    
    REPORT(2, "Blocking (Setup/Memcpy/Output/CPU/Launch GPU)", NearForce.Elapsed());
    REPORT(3, "Cell bookkeeping", CellBookkeeping.Elapsed());
    REPORT(3, "Sink counting", CountSinks.Elapsed());
    REPORT(3, "Find unpinnable slabs", FindUnpin.Elapsed());
    REPORT(3, "Setup GPU", SetupGPUTimers[0].Elapsed());
    //REPORT(4, "Waiting for GPU", SetupGPUTimers[3].Elapsed());
    REPORT(4, "Init and Launch Setup Kernels", SetupGPUTimers[1].Elapsed());
    REPORT(5, "Launch memcpy kernels", SetupGPUTimers[4].Elapsed());
    REPORT(5, "cudaFree (implicit sync)", SetupGPUTimers[5].Elapsed());
    REPORT(4, "Waiting for GPU", SetupGPUTimers[2].Elapsed());
    REPORT(3, "Unpin slabs", GPUUnpinTimer.Elapsed());
    REPORT(3, "GPU Directs Launch", GPUDirectsLaunchTimer.Elapsed());
    REPORT(3, "CPU Fallback", CPUFallbackTimer.Elapsed());
    double gdi_cpu = JJ->DirectInteractions_CPU/1e9;
    fprintf(timingfile,"---> %6.3f GDIPS, %6.3f Gdirects, %6.3f Mpart/sec", gdi_cpu/(thistime+1e-15), gdi_cpu, JJ->NSink_CPU/(thistime+1e-15)/1e6);
    
    REPORT(2, "Non-blocking (GPU Directs)", GPUTimers[NGPU].Elapsed());
    //double gdi_gpu = JJ->DirectInteractions_GPU()/1e9;
    fprintf(timingfile,"---> %6.3f GDIPS, %6.3f Gdirects, %6.3f Mpart/sec", gdi_gpu/(thistime+1e-15), gdi_gpu, P.np/(thistime+1e-15)/1e6);
    
    for(int g = 0; g < NGPU; g++){
        char buffer[100];
        const char* fmt = "GPU %d Directs:";
        sprintf(buffer, fmt, g);
        REPORT(3, buffer, GPUTimers[g].Elapsed());
        double gdi_gpu = JJ->DirectInteractions_GPU/1e9;
        fprintf(timingfile,"---> %6.3f GDIPS, %6.3f Gdirects, %6.3f Mpart/sec", gdi_gpu/(thistime+1e-15), gdi_gpu, P.np/(thistime+1e-15)/1e6);
    }

    denom = NearForce.Elapsed();
    REPORT(2, "NearForce Setup", JJ->Construction);
    denom = JJ->Construction;
    REPORT(3, "NearForce Fill Source", JJ->FillSourceLists);
    denom = JJ->FillSourceLists;
    REPORT(4, "Count Sources",JJ->CountSources);
    REPORT(4, "Calc Source Offset",JJ->CalcSourceBlocks);
    REPORT(4, "Fill Sources",JJ->FillSources);
    denom = JJ->Construction;
    REPORT(3, "NearForce Fill Sink", JJ->FillSinkLists);
    denom = JJ->FillSinkLists;
    REPORT(4, "Count Sinks",JJ->CountSinks);
    REPORT(4, "Calc Sink Offset",JJ->CalcSinkBlocks);
    REPORT(4, "Fill Sinks",JJ->FillSinks);
    denom = JJ->Construction;
    REPORT(3, "NearForce Fill Interaction", JJ->FillInteractionList);
    denom = NearForce.Elapsed();
    REPORT(2, "SIC Execution", JJ->SICExecute.Elapsed());
    //REPORT(2, "NearForce Core", JJ->DeviceAcceleration.Elapsed());
    //REPORT(2, "NearForce Thread Sort", JJ->GPUThreadSort);
    //REPORT(2, "NearForce Thread Spin", JJ->GPUThreadSpin);
    //fprintf(timingfile,"---> Mean: %e s\t +- %e s",JJ->GPUThreadSpinAverage,JJ->GPUThreadSpinSigma);
    //REPORT(2, "NearForce Thread Device kernel + copy", JJ->GPUDeviceDirectsPlusCopy);
    //REPORT(2, "NearForce Thread blockify/deblockify", JJ->GPUThreadBlockifyInto64Blocks);
    //denom = JJ->GPUThreadBlockifyInto64Blocks;
    //REPORT(3, "Count Source Blocks",JJ->SourceBlockCount);
    //REPORT(3, "Fill Source Blocks",JJ->SourceBlockFill);
    //REPORT(3, "Count Sink Blocks",JJ->SinkBlockCount);
    //REPORT(3, "Fill Sink Blocks",JJ->SinkBlockFill);
    //REPORT(3, "Memcpy and Free",JJ->MemcpyAndFree);
    //REPORT(3, "Deblockify",JJ->DeblockifyAndCoAdd);
    //denom = NearForce.Elapsed();
    //REPORT(2, "NearForce Total", JJ->Total.Elapsed());

    //fprintf(timingfile, "\n\t\t Near Force Device Breakdown:");
    //fprintf(timingfile, "\n\t\t\t Total Data Transfered to GPU: %6.3f GB",JJ->TotalGPUDataSize/1e9);
    //fprintf(timingfile, "\n\t\t\t -->Inferred Minimum GPU Bandwith: %6.3f GB/s",JJ->TotalGPUDataSize/(TimeStepWallClock.Elapsed()*1e9));
    //fprintf(timingfile, "\n\t\t\t Average Data Transfered Per GPU kernel: %6.4f GB",JJ->MeanGPUDataSizePerKernel/1e9);
    //fprintf(timingfile, "\n\t\t\t -->Inferred PCI-E Bandwith (~8GB/s Min) Limited Kernel Time: %6.3f s",JJ->MeanGPUDataSizePerKernel/1e9 *(1/8.0));
    //fprintf(timingfile, "\n\t\t\t GPU Data Transfer Size Range: %6.3f GB",JJ->GPUDataSizeRange/1e9);
    //fprintf(timingfile, "\n\t\t\t Kernels Per Device Range: %6.3f",JJ->GPUKernelsRange);
#endif

    fprintf(timingfile, "\n\n Subdivisions of Compute Multipole:");
    denom = ComputeMultipoles.Elapsed();
    REPORT(2, "Multipole ASM (P)", MF->MultipoleASM.Elapsed());
    REPORT(2, "Multipole C to R (S)", MF->MultipoleC2R.Elapsed());
    REPORT(2, "Multipole FFT (S)", MF->FFTMultipole.Elapsed());

    if (TY!=NULL) {    // Not in IC steps
	denom = TaylorCompute.Elapsed();
	fprintf(timingfile, "\n\n Subdivisions of Taylor Evaluate:");
	REPORT(2, "Taylor Computation (S)", TaylorCompute.Elapsed());
	fprintf(timingfile,"---> %6.3f Mpart/sec", P.np/(thistime+1e-15)/1e6 );
	REPORT(2, "Taylor FFT (S)", TY->FFTTaylor.Elapsed());
	REPORT(2, "Taylor R to C (P)", TY->TaylorR2C.Elapsed());
	REPORT(2, "Taylor ASM (P)", TY->TaylorASM.Elapsed());
	REPORT(2, "Taylor Redlack (P)", RL->TaylorRedlack.Elapsed());
    }

    return;
}

