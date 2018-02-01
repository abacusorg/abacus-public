/* reporting.cpp
 *
 * This writes timing outputs.  It should provide a nicely formatted
 * summary of the timings from the various portions of the code.
 *
 *
 * TODO: Need to uncomment the Multipole and Taylor portions.
 *
 * FIXME: DWF recommends a major refactor of how we handle reports. 
 *  Ideally, we should have a top level class TimingReport that 
 *  holds a registry of ModuleTimings objects (an interface
 *  that is implemented by each "section" of the timings.
 *  These ModuleTimings are responsible for supplying
 *  formatted strings for each timing they are responsible.
 *  These strings can inlcude any other information (like counts)
 *  held by the timer class.
 *
 *  The TimingReport object itself is responsible for producing
 *  the final timing report. Each ModuleTimings is registered along
 *  with its format options---these should at least include desired
 *  indentation level, but could also support some sort of nested timing
 *  system.
 *
 *  This will allow timing reports to be defined closer to the code they are 
 *  reporting on. This makes adding new sections to the timing file much easier, 
 *  and allows seperation of timer setup from main code. This is much better
 *  than the current solution, where the timings report requires knowledge
 *  of every other class's deep implementation details. This makes the timing
 *  report tend to lag the rest of the code, since changing it is a major process.
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

    a = IL->FinishSort.Elapsed();           fprintf(timingfile,"\t \t \t \t             FinishSort: %e seconds (%2.2f %%)\n", a, a*g );
    a = IL->FinishPartition.Elapsed();      fprintf(timingfile,"\t \t \t \t        FinishPartition: %e seconds (%2.2f %%)\n", a, a*g );

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
    //REPORT(0, "Epilogue", epilogue.Elapsed()); total += thistime;
    REPORT(0, "TimeStep", TimeStepWallClock.Elapsed()); total += thistime;
    REPORT(0, "SingleStep Setup", SingleStepSetup.Elapsed()); total += thistime;
    REPORT(0, "SingleStep TearDown", SingleStepTearDown.Elapsed()); total += thistime;
    REPORT(0, "Unaccounted", WallClockDirect.Elapsed()-total);
    fprintf(timingfile, "\n");
    
    // Collect stats on IO
    
#ifdef IOTHREADED
    #define PRINT_IO_THREAD(dir) fprintf(timingfile, " [thread %d]", GetIOThread(dir));
#else
    #define PRINT_IO_THREAD(dir)
#endif
    
#define IOSTATS(NAME, DESC)\
    do { \
        double total_##NAME##_time = 0, total_##NAME##_bytes = 0; \
        for(auto &iter : NAME##Time){ \
            total_##NAME##_time += iter.second.Elapsed(); \
            total_##NAME##_bytes += NAME##Bytes[iter.first]; \
        } \
        denom = WallClockDirect.Elapsed(); \
        REPORT(0, DESC, total_##NAME##_time); \
        denom = thistime+1e-15; \
        fprintf(timingfile, "---> %6.1f MB/sec on %6.2f GB", total_##NAME##_bytes/(thistime+1e-15)/1e6, total_##NAME##_bytes/1e9); \
        for(auto &iter : NAME##Time){ \
            double time = iter.second.Elapsed(); \
            const char* dir = iter.first.c_str(); \
            auto bytes = NAME##Bytes[iter.first]; \
            REPORT(1, dir, time); \
            fprintf(timingfile, "---> %6.1f MB/sec on %6.2f GB", bytes/(time+1e-15)/1e6, bytes/1e9); \
            PRINT_IO_THREAD(dir); \
        } \
    } while (0)
    
    IOSTATS(BlockingIORead, "Blocking Disk Reads");
    IOSTATS(BlockingIOWrite, "Blocking Disk Writes");
#ifdef IOTHREADED
    IOSTATS(NonBlockingIORead, "Non-Blocking Disk Reads");
    IOSTATS(NonBlockingIOWrite, "Non-Blocking Disk Writes");
#endif

    fprintf(timingfile, "\n\nBreakdown of TimeStep: ");
    total = 0.0;
    denom = TimeStepWallClock.Elapsed();
    REPORT(1, "FetchSlabs", FetchSlabs.Elapsed()); total += thistime;
    
    // Add up the time spent spinning
    double spinning = Dependency::global_spin_timer.Elapsed();
    
#ifdef CUDADIRECT
    REPORT(1, "Transpose Positions", TransposePos.Elapsed()); total += thistime;
    
    int NGPU = GetNGPU();
    
    // Compute some GPU timing totals
    double total_copy_time = 0, total_execution_time = 0, total_copyback_time = 0, total_gpu_time = 0;
    double total_GB_to = 0, total_GB_from = 0, total_sinks = 0, total_sources = 0;
    for(int g = 0; g < NGPU*DirectBPD; g++){
        total_GB_to += JJ->GB_to_device[g];
        total_GB_from += JJ->GB_from_device[g];
        total_sinks += JJ->DeviceSinks[g];
        total_sources += JJ->DeviceSources[g];
    }
    
    // Measures the total amount of time we have at least one GPU thread running
    assertf(!SetInteractionCollection::GPUThroughputTimer.timeron, "GPU throughput timer still on! Atomic thread counting failure?\n");
    double GPUThroughputTime = SetInteractionCollection::GPUThroughputTimer.Elapsed();
    
    REPORT(1, "NearForce [blocking]", NearForce.Elapsed()); total += thistime;
    fprintf(timingfile,"---> %6.3f Mpart/sec", P.np/(thistime+1e-15)/1e6 );
    REPORT(1, "NearForce [non-blocking]", GPUThroughputTime); //total += thistime;
    double gdi_gpu = JJ->TotalDirectInteractions_GPU/1e9;
    fprintf(timingfile,"---> %6.3f effective GDIPS, %6.3f Gdirects, %6.3f Mpart/sec", gdi_gpu/(thistime+1e-15), gdi_gpu, P.np/(thistime+1e-15)/1e6);
    double total_di = (JJ->DirectInteractions_CPU +JJ->TotalDirectInteractions_GPU)/1e9;
#else
    REPORT(1, "NearForce", NearForce.Elapsed()); total += thistime;
    fprintf(timingfile,"---> %6.3f Mpart/sec", P.np/(thistime+1e-15)/1e6 );
#endif
    REPORT(1, "TaylorForce", TaylorForce.Elapsed()); total += thistime;
        fprintf(timingfile,"---> %6.3f Mpart/sec", P.np/(thistime+1e-15)/1e6 );
    REPORT(1, "Kick", Kick.Elapsed()); total += thistime;
        fprintf(timingfile,"---> %6.3f Mpart/sec", P.np/(thistime+1e-15)/1e6 );
	double gf_total = MakeCellGroups.Elapsed() + FindCellGroupLinks.Elapsed() + DoGlobalGroups.Elapsed();
    if (GFC != NULL){
        REPORT(1, "Group Finding", gf_total); total += thistime;
        fprintf(timingfile,"---> %6.3f Mpart/sec",P.np/(thistime+1e-15)/1e6);
    }
    REPORT(1, "Output", Output.Elapsed()); total += thistime;
    fprintf(timingfile,"---> %6.3f Mpart/sec",P.np/(thistime+1e-15)/1e6);
    if (GFC != NULL){
        REPORT(1, "Microstep", Microstep.Elapsed()); total += thistime;
        fprintf(timingfile,"---> %6.3f Mpart/sec",P.np/(thistime+1e-15)/1e6);
    }
   
    if(WriteState.Do2LPTVelocityRereading){
        REPORT(1, "Velocity Re-reading for LPT", LPTVelocityReRead.Elapsed()); total += thistime;
        fprintf(timingfile,"---> %6.3f Mpart/sec", P.np/(thistime+1e-15)/1e6 );
    }
    REPORT(1, "Drift", Drift.Elapsed()); total += thistime;
        fprintf(timingfile,"---> %6.3f Mpart/sec", P.np/(thistime+1e-15)/1e6 );
    REPORT(1, "Finish", Finish.Elapsed()); total += thistime;
        fprintf(timingfile,"---> %6.3f Mpart/sec", P.np/(thistime+1e-15)/1e6 );
    REPORT(1, "Spinning", spinning); total += thistime;
    REPORT(1, "Unaccounted", TimeStepWallClock.Elapsed()-total);

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
    fprintf(timingfile,"\n\t\tSigma(STD): %.2g s\t Min: %.2e s\t Max: %.2e s ",slabforcetimesigma,slabforcemintime,slabforcemaxtime);
    REPORT(1,"Mean Force Latency",slabforcelatencymean);
    fprintf(timingfile,"\n\t\tSigma(STD): %.2e s\t Min: %.2e s\t Max: %.2e s ",slabforcelatencysigma,slabforceminlatency,slabforcemaxlatency);
    
    fprintf(timingfile, "\n\n Subdivisions of Near Force:");
    double gdi_cpu = JJ->DirectInteractions_CPU/1e9;  // Measure per-core load balancing?
#ifdef CUDADIRECT
    fprintf(timingfile, "\n\t Notes about non-blocking timing:\n");
    fprintf(timingfile, "\t -\"Directs Throughput\" is the wall clock time while at least one GPU thread is running (copy or compute).\n");
    fprintf(timingfile, "\t -\"Effective\" GDIPS is based on this throughput.\n");
    denom = NearForce.Elapsed();
    REPORT(1, "Blocking", NearForce.Elapsed());
    REPORT(2, "Calculate Direct Splits", JJ->CalcSplitDirects.Elapsed());
    
    REPORT(2, "Construct Pencils", JJ->Construction);
        denom = JJ->Construction;
        REPORT(3, "Plan Sinks", JJ->FillSinkLists);
            denom = JJ->FillSinkLists;
            REPORT(4, "Count Sinks",JJ->CountSinks);
            REPORT(4, "Calc Sink Blocks",JJ->CalcSinkBlocks);
            REPORT(4, "Allocate Acceleration Memory", JJ->AllocAccels);
        denom = JJ->Construction;
        REPORT(3, "Plan Sources", JJ->FillSourceLists);
            denom = JJ->FillSourceLists;
            REPORT(4, "Count Sources",JJ->CountSources);
            REPORT(4, "Calc Source Blocks",JJ->CalcSourceBlocks);
        denom = JJ->Construction;
        REPORT(3, "Fill Interaction", JJ->FillInteractionList);
    denom = NearForce.Elapsed();
    REPORT(2, "Dispatch Interaction", JJ->SICExecute.Elapsed());
    REPORT(2, "CPU Fallback", JJ->CPUFallbackTimer.Elapsed());
            fprintf(timingfile,"---> %6.3f GDIPS, %6.3f Gdirects, %6.3f Mpart/sec", gdi_cpu/(thistime+1e-15), gdi_cpu, JJ->NSink_CPU/(thistime+1e-15)/1e6);
    
    denom = JJ->DeviceThreadTimer;
    char str[1024];  sprintf(str, "Non-Blocking (thread-seconds, %d threads)", NGPU*DirectBPD);
    REPORT(1, str, JJ->DeviceThreadTimer);
        REPORT(2, "Fill Sinks", JJ->FillSinks);
                fprintf(timingfile,"---> %6.1f MB/s, %6.3f MSink/sec", total_sinks*sizeof(FLOAT3)/1e6/thistime, total_sinks/1e6/thistime);
        REPORT(2, "Fill Sources", JJ->FillSources);
                fprintf(timingfile,"---> %6.1f MB/s, %6.3f MSource/sec", total_sources*sizeof(FLOAT3)/1e6/thistime, total_sources/1e6/thistime);
        REPORT(2, "Launch Kernels", JJ->LaunchDeviceKernels);
        REPORT(2, "Wait for GPU Result", JJ->WaitForResult);
        REPORT(2, "Copy Accel from Pinned", JJ->CopyAccelFromPinned);
                fprintf(timingfile,"---> %6.1f MB/s, %6.3f MSink/sec", total_sinks*sizeof(FLOAT3)/1e6/thistime, total_sinks/1e6/thistime);  // same number of accels as sinks
        
    denom = GPUThroughputTime;
    REPORT(1, "Non-Blocking Throughput (Wall Clock)", GPUThroughputTime);
            fprintf(timingfile,"\n\t\t\t\t---> %6.3f effective GDIPS, %6.3f Gdirects, %6.3f Mpart/sec, %6.3f Msink/sec", gdi_gpu/(thistime+1e-15), gdi_gpu, P.np/(thistime+1e-15)/1e6, total_sinks/(thistime+1e-15)/1e6);
            fprintf(timingfile,"\n\t\t\t\t---> with %d device threads, estimate %.1f%% thread concurrency", NGPU*DirectBPD, (JJ->DeviceThreadTimer - GPUThroughputTime)/(JJ->DeviceThreadTimer - JJ->DeviceThreadTimer/(NGPU*DirectBPD))*100);
            
        fprintf(timingfile, "\n    Device stats:\n");
        for(int g = 0; g < NGPU*DirectBPD; g++){
            fprintf(timingfile, "        Device thread %d (GPU %d):", g, g % NGPU);
            fprintf(timingfile, " %.2f GB to device, %.2f GB from device, %.2f Msink, %.2f Gdirects\n", JJ->GB_to_device[g], JJ->GB_from_device[g], JJ->DeviceSinks[g]/1e6, JJ->DirectInteractions_GPU[g]/1e9);
        }
#else
    REPORT(1, "CPU directs", NearForce.Elapsed());
        fprintf(timingfile,"---> %6.3f GDIPS, %6.3f Gdirects, %6.3f Mpart/sec", gdi_cpu/(thistime+1e-15), gdi_cpu, JJ->NSink_CPU/(thistime+1e-15)/1e6);

#endif
    
    if (TY!=NULL) {    // Not in IC steps
        denom = TaylorCompute.Elapsed();
        fprintf(timingfile, "\n Subdivisions of Taylor Evaluate:");
        REPORT(1, "Taylor Computation", TaylorCompute.Elapsed());
        fprintf(timingfile,"---> %6.3f Mpart/sec", P.np/(thistime+1e-15)/1e6 );
        REPORT(1, "Compute Cell Offsets", TY->ConstructOffsets.Elapsed());
        REPORT(1, "Taylor FFT", TY->FFTTaylor.Elapsed());
        REPORT(1, "Taylor R to C", TY->TaylorR2C.Elapsed());
        REPORT(1, "Taylor ASM", TY->TaylorASM.Elapsed());
        REPORT(1, "Taylor Redlack", RL->TaylorRedlack.Elapsed());
    }

    fprintf(timingfile, "\n\n Subdivisions of Kick:");
    denom = Kick.Elapsed();
    REPORT(1, "Finalize accelerations", JJ->FinalizeTimer.Elapsed());
    denom = thistime;
    REPORT(2, "Bookkeeping/Fetch Timings", JJ->FinalizeBookkeeping.Elapsed());
    REPORT(2, "Copy Pencil to Slab", JJ->CopyPencilToSlab.Elapsed());
    denom = Kick.Elapsed();
    REPORT(1, "Add Near + Far Accel", AddAccel.Elapsed());
    REPORT(1, "Kick Cell", KickCellTimer.Elapsed());
    
	fprintf(timingfile, "\n\n Breakdown of Group Finding:");
	denom = gf_total;
	REPORT(1, "MakeCellGroups", MakeCellGroups.Elapsed());
    fprintf(timingfile,"---> %6.3f Mpart/sec",P.np/(thistime+1e-15)/1e6);
    REPORT(1, "FindCellGroupLinks", FindCellGroupLinks.Elapsed());
    fprintf(timingfile,"---> %6.3f Mpart/sec",P.np/(thistime+1e-15)/1e6);
    REPORT(1, "DoGlobalGroups", DoGlobalGroups.Elapsed());
    fprintf(timingfile,"---> %6.3f Mpart/sec",P.np/(thistime+1e-15)/1e6);
	
	// Now write some detailed multiplicity and timing stats to lastrun.grouplog
	if(GFC != NULL)
		GFC->report();
	
    fprintf(timingfile, "\n\n Breakdown of Output Step:");
    denom = Output.Elapsed();
    REPORT(1, "TimeSlice", OutputTimeSlice.Elapsed());
    fprintf(timingfile,"---> %6.3f Mpart/sec",P.np/(thistime+1e-15)/1e6);
    REPORT(1, "LightCone", OutputLightCone.Elapsed());
    fprintf(timingfile,"---> %6.3f Mpart/sec",P.np/(thistime+1e-15)/1e6);
    REPORT(1, "Binning", OutputBin.Elapsed());
    fprintf(timingfile,"---> %6.3f Mpart/sec",P.np/(thistime+1e-15)/1e6);

    fprintf(timingfile, "\n\n Breakdown of Microstep Step:");
    denom = Microstep.Elapsed();
    REPORT(1, "Slab Kick", MicrostepSlabKick.Elapsed());
    fprintf(timingfile,"---> %6.3f Mpart/sec",P.np/(thistime+1e-15)/1e6);
    if(GFC != NULL){
        REPORT(1, "Group Kick", MicrostepGroupKick.Elapsed());
        fprintf(timingfile,"---> %6.3f M_group_part/sec", GFC->L0stats.tot/(thistime+1e-15)/1e6);
        REPORT(1, "CPU Microsteps", MicrostepCPU.Elapsed());
        fprintf(timingfile,"---> %6.3f M_group_part/sec",GFC->L0stats.tot/(thistime+1e-15)/1e6);
        REPORT(1, "Scatter Groups", GFC->ScatterGroups.Elapsed());
        fprintf(timingfile,"---> %6.3f M_group_part/sec",GFC->L0stats.tot/(thistime+1e-15)/1e6);
    }

    fprintf(timingfile, "\n\n Subdivisions of Drift:");
    denom = Drift.Elapsed();
    REPORT(1, "Move",         DriftMove.Elapsed());
    REPORT(1, "Rebin",        DriftRebin.Elapsed());
    REPORT(1, "Inserting",    DriftInsert.Elapsed());

    fprintf(timingfile, "\n\n Subdivisions of Compute Multipole:");
    denom = ComputeMultipoles.Elapsed();
    REPORT(1, "Compute Cell Offsets", MF->ConstructOffsets.Elapsed());
    REPORT(1, "Multipole ASM", MF->MultipoleASM.Elapsed());
    REPORT(1, "Multipole C to R", MF->MultipoleC2R.Elapsed());
    REPORT(1, "Multipole FFT", MF->FFTMultipole.Elapsed());
    
    fprintf(timingfile, "\n\n Breakdown of Finish Step:");
    denom = Finish.Elapsed();
    REPORT(1, "Partition Insert List", IL->FinishPartition.Elapsed());
    REPORT(1, "Sort Insert List", IL->FinishSort.Elapsed());
    fprintf(timingfile,"---> %6.3f Mitems/sec (%.2g items)",IL->n_sorted/(thistime+1e-15)/1e6, (double) IL->n_sorted);
    REPORT(1, "Index Cells", FinishCellIndex.Elapsed());
    REPORT(1, "Merge", FinishMerge.Elapsed());
    REPORT(1, "Compute Multipoles", ComputeMultipoles.Elapsed());
    fprintf(timingfile,"---> %6.3f Mpart/sec", P.np/(thistime+1e-15)/1e6 );
    REPORT(1, "Write Particles", WriteMergeSlab.Elapsed());
    REPORT(1, "Write Multipoles", WriteMultipoleSlab.Elapsed());
    
    denom = TimeStepWallClock.Elapsed();
    REPORT(0, "\nAllocate Arena Memory", LBW->ArenaMalloc.Elapsed());
    
    fprintf(timingfile, "\n\n Reasons for Spinning:");
    fprintf(timingfile, "\n\t Note: may add up to >100%% if there are multiple simultaneous reasons for spinning");
    denom = spinning;
    REPORT(1, "Not enough RAM to load slabs", Dependency::spin_timers[NOT_ENOUGH_RAM].Elapsed());
    REPORT(1, "Waiting for slab IO", Dependency::spin_timers[WAITING_FOR_IO].Elapsed());
    REPORT(1, "Waiting for GPU", Dependency::spin_timers[WAITING_FOR_GPU].Elapsed());
    
    return;
}
