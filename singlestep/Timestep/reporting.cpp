/* reporting.cpp
 *
 * This writes timing outputs.  It should provide a nicely formatted
 * summary of the timings from the various portions of the code.
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

#define REPORT(tabs, str, a) \
     do { \
         fprintf(timingfile, "\n"); \
         for (int i=0; i<tabs; i++) fprintf(timingfile, "    "); \
         thistime = a; fprintf(timingfile, "%-30s: %10.2e sec (%5.2f%%) ", str, thistime, denom ? 100*thistime/denom : 0.); \
     } while(0)

#define REPORT_RATE() \
     do { \
        fprintf(timingfile,"---> %6.2f Mpart/sec", thistime ? P.np/thistime/1e6 : 0.); \
     } while(0)

void ReportTimings(FILE * timingfile) {
    double thistime, denom, total;
    denom = WallClockDirect.Elapsed();
    // fprintf(timingfile,"DirectWallClock:  %2.2e seconds\n", WallClockDirect.Elapsed() );
    REPORT(0, "Total Wall Clock Time", WallClockDirect.Elapsed()); 
    fprintf(timingfile,"---> %6.3f Mpart/sec", thistime ? P.np/thistime/1e6 : 0.);
    fprintf(timingfile,"\n");

    total = 0.0;
    REPORT(0, "SingleStep Setup", SingleStepSetup.Elapsed()); total += thistime;
    REPORT(0, "Prologue", prologue.Elapsed()); total += thistime;
    REPORT(0, "TimeStep", TimeStepWallClock.Elapsed()); total += thistime;
    REPORT(0, "SingleStep TearDown", SingleStepTearDown.Elapsed()); total += thistime;
    //REPORT(0, "Epilogue", epilogue.Elapsed()); total += thistime;  // The timings get written at the start of the Epilogue. Can't time what hasn't happened!
    REPORT(0, "Unaccounted", WallClockDirect.Elapsed()-total);
    fprintf(timingfile, "\n");
    
    // Collect stats on IO
#ifndef IOTHREADED
    #define niothreads 1
    #define GetIOThread(a) (1)
#endif
// iter.first is the dir, iter.second is the timer
#define ACCUMULATE_THREAD_TOTALS(NAME, RW)\
    do{\
        for (auto &iter : NAME##Time){\
            if(GetIOThread(iter.first.c_str())-1 == i){\
                total_##RW##_time += iter.second.Elapsed();\
                total_##RW##_bytes += NAME##Bytes[iter.first];\
            }\
        }\
    } while(0)

#define REPORT_DIR_IOSTATS(NAME, RW, BLOCKING)\
    do{\
        for (auto &iter : NAME##Time){\
            if(GetIOThread(iter.first.c_str())-1 == i){\
                double this_io_time = iter.second.Elapsed();\
                double this_io_bytes = NAME##Bytes[iter.first];\
                REPORT(1, iter.first.c_str(), this_io_time);\
                fprintf(timingfile, "---> %6.1f MB/s " #RW " on %6.2f GB ["  #BLOCKING "]", this_io_time ? this_io_bytes/this_io_time/1e6 : 0., this_io_bytes/1e9);\
            }\
        }\
    } while(0)

    for(int i = 0; i < niothreads; i++){
        double total_read_time = 0., total_read_bytes = 0.;
        double total_write_time = 0., total_write_bytes = 0.;

        ACCUMULATE_THREAD_TOTALS(BlockingIORead, read);
        ACCUMULATE_THREAD_TOTALS(BlockingIOWrite, write);
        ACCUMULATE_THREAD_TOTALS(NonBlockingIORead, read);
        ACCUMULATE_THREAD_TOTALS(NonBlockingIOWrite, write);

        double total_time = total_read_time + total_write_time;
        double total_bytes = total_read_bytes + total_write_bytes;

        denom = WallClockDirect.Elapsed();
#ifdef IOTHREADED
        std::string threadname = "IO Thread " + std::to_string(i+1ll);
        REPORT(0, threadname.c_str(), total_time);
#else
        REPORT(0, "Blocking IO", total_time);
#endif
        fprintf(timingfile, "---> %6.1f MB/s read, %6.1f MB/s write", total_read_time ? total_read_bytes/total_read_time/1e6 : 0., total_write_time ? total_write_bytes/total_write_time/1e6 : 0.);
        denom = thistime;
        REPORT_DIR_IOSTATS(BlockingIORead, read, blocking);
        REPORT_DIR_IOSTATS(BlockingIOWrite, write, blocking);
        REPORT_DIR_IOSTATS(NonBlockingIORead, read, non-blocking);
        REPORT_DIR_IOSTATS(NonBlockingIOWrite, write, non-blocking);
    }
#undef niothreads

    fprintf(timingfile, "\n\nBreakdown of TimeStep: ");
    total = 0.0;
    denom = TimeStepWallClock.Elapsed();
    REPORT(1, "FetchSlabs", FetchSlabs.Elapsed()); total += thistime;
    
    // Add up the time spent spinning
    double spinning = Dependency::global_spin_timer.Elapsed();
    
#ifdef CUDADIRECT
    REPORT(1, "Transpose Positions", TransposePos.Elapsed()); total += thistime;
    
    int NGPU = GetNGPU();

    if(JJ){
        REPORT(1, "NearForce [blocking]", NearForce.Elapsed()); total += thistime;
        REPORT(1, "NearForce [non-blocking]", JJ->GPUThroughputTime); //total += thistime;
        fprintf(timingfile,"---> %6.2f effective GDIPS, %6.2f Gdirects, %.2f Mpart/sec", thistime ? JJ->gdi_gpu/thistime : 0, JJ->gdi_gpu, thistime ? P.np/thistime/1e6 : 0.);
        double total_di = (JJ->DirectInteractions_CPU +JJ->TotalDirectInteractions_GPU)/1e9;
    }
#else
    REPORT(1, "NearForce", NearForce.Elapsed()); total += thistime;
    REPORT_RATE();
#endif
    REPORT(1, "TaylorForce", TaylorForce.Elapsed()); total += thistime;
    REPORT_RATE();
    REPORT(1, "Kick", Kick.Elapsed()); total += thistime;
    REPORT_RATE();
    double gf_total = MakeCellGroups.Elapsed() + FindCellGroupLinks.Elapsed() + DoGlobalGroups.Elapsed() + FinishGroups.Elapsed();
    if (GFC != NULL){
        REPORT(1, "Group Finding", gf_total); total += thistime;
        REPORT_RATE();
    }
    REPORT(1, "Output", Output.Elapsed()); total += thistime;
    REPORT_RATE();
    if (GFC != NULL){
        REPORT(1, "Microstep", Microstep.Elapsed()); total += thistime;
        REPORT_RATE();
    }
   
    if(WriteState.Do2LPTVelocityRereading){
        REPORT(1, "Velocity Re-reading for LPT", LPTVelocityReRead.Elapsed()); total += thistime;
        REPORT_RATE();
    }
    REPORT(1, "Drift", Drift.Elapsed()); total += thistime;
        REPORT_RATE();
    REPORT(1, "Finish", Finish.Elapsed()); total += thistime;
        REPORT_RATE();
    REPORT(1, "Spinning", spinning); total += thistime;
    REPORT(1, "Unaccounted", TimeStepWallClock.Elapsed()-total);

    fprintf(timingfile, "\n\nBreakdown per slab (Wall Clock)");
    double slabforcetimemean = 0;
    double slabforcetimesigma = 0;
    double slabforcemaxtime = 0;
    double slabforcemintime =  1e9;
    double slabforcelatencymean = 0;
    double slabforcelatencysigma = 0;
    double slabforcemaxlatency = 0;
    double slabforceminlatency = 1e9;

    if(JJ){
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
    }
    slabforcetimesigma =  sqrt(slabforcetimesigma    - slabforcetimemean*slabforcetimemean);
    slabforcelatencysigma=sqrt(slabforcelatencysigma - slabforcelatencymean*slabforcelatencymean);

    denom = WallClockDirect.Elapsed()/P.cpd;
    REPORT(1,"Mean Force Computation",slabforcetimemean);
    fprintf(timingfile,"\n\t\tSigma(STD): %.2g s\t Min: %.2e s\t Max: %.2e s ",slabforcetimesigma,slabforcemintime,slabforcemaxtime);
    REPORT(1,"Mean Force Latency",slabforcelatencymean);
    fprintf(timingfile,"\n\t\tSigma(STD): %.2e s\t Min: %.2e s\t Max: %.2e s ",slabforcelatencysigma,slabforceminlatency,slabforcemaxlatency);
    
    if(JJ){
        fprintf(timingfile, "\n\nBreakdown of Near Force:");
        double gdi_cpu = JJ->DirectInteractions_CPU/1e9;  // Measure per-core load balancing?
#ifdef CUDADIRECT
        fprintf(timingfile, "\n\tNotes about non-blocking timing:\n");
        fprintf(timingfile, "\t- \"Directs Throughput\" is the wall clock time while at least one GPU thread is running (copy or compute).\n");
        fprintf(timingfile, "\t- \"Effective\" GDIPS is based on this throughput.\n");
        denom = NearForce.Elapsed();
        REPORT(1, "Blocking", NearForce.Elapsed());
        REPORT(2, "Calculate Splits", JJ->CalcSplitDirects.Elapsed());
        REPORT(2, "Construct Pencils", JJ->SICConstruction.Elapsed());
        REPORT(2, "Dispatch Interaction", JJ->SICExecute.Elapsed());
        REPORT(2, "CPU Fallback", JJ->CPUFallbackTimer.Elapsed());
                fprintf(timingfile,"---> %6.2f GDIPS, %6.2f Gdirects, %6.2f Mpart/sec\n", thistime ? gdi_cpu/thistime : 0., gdi_cpu, thistime ? JJ->NSink_CPU/thistime/1e6 : 0.);
        
        denom = JJ->DeviceThreadTimer;
        char str[1024];  sprintf(str, "Non-Blocking (thread-seconds, %d threads)", NGPU*DirectBPD);
        REPORT(1, str, JJ->DeviceThreadTimer);
            REPORT(2, "Fill Sinks", JJ->FillSinks);  // filling is unpadded
                    fprintf(timingfile,"---> %6.1f MB/s, %6.2f MSink/sec", thistime ? JJ->total_sinks*sizeof(posstruct)/1e6/thistime : 0., thistime ? JJ->total_sinks/1e6/thistime : 0.);
            REPORT(2, "Fill Sources", JJ->FillSources);
                    fprintf(timingfile,"---> %6.1f MB/s, %6.2f MSource/sec", thistime ? JJ->total_sources*sizeof(posstruct)/1e6/thistime : 0., thistime ? JJ->total_sources/1e6/thistime : 0.);
            REPORT(2, "Launch Kernels", JJ->LaunchDeviceKernels);
            REPORT(2, "Wait for GPU Result", JJ->WaitForResult);
            REPORT(2, "Copy Accel from Pinned", JJ->CopyAccelFromPinned);
                    fprintf(timingfile,"---> %6.1f MB/s, %6.2f MSink/sec\n", thistime ? JJ->total_sinks*sizeof(accstruct)/1e6/thistime : 0., thistime ? JJ->total_sinks/1e6/thistime : 0.);  // same number of accels as sinks
            
        denom = JJ->GPUThroughputTime;
        REPORT(1, "Non-Blocking Throughput (Wall Clock)", JJ->GPUThroughputTime);
                fprintf(timingfile,"\n\t\t\t\t---> %6.2f effective GDIPS, %6.2f Mpart/sec, %6.2f Msink/sec", thistime ? JJ->gdi_gpu/thistime : 0., thistime ? P.np/thistime/1e6 : 0., thistime ? JJ->total_sinks/thistime/1e6 : 0.);
                fprintf(timingfile,"\n\t\t\t\t---> %6.2f Gdirects, %6.2f padded Gdirects", JJ->gdi_gpu, JJ->gdi_padded_gpu);
                fprintf(timingfile,"\n\t\t\t\t---> with %d device threads, estimate %.1f%% thread concurrency", NGPU*DirectBPD, (JJ->DeviceThreadTimer - JJ->GPUThroughputTime)/(JJ->DeviceThreadTimer - JJ->DeviceThreadTimer/(NGPU*DirectBPD))*100);
                
            fprintf(timingfile, "\n    Device stats:\n");
            for(int g = 0; g < NGPU*DirectBPD; g++){
                fprintf(timingfile, "        Device thread %d (GPU %d):", g, g % NGPU);
                fprintf(timingfile, " %.2f GB to device, %.2f GB from device, %.2f Msink, %.2f Gdirects\n", JJ->GB_to_device[g], JJ->GB_from_device[g], JJ->DeviceSinks[g]/1e6, JJ->DirectInteractions_GPU[g]/1e9);
            }
#else
        REPORT(1, "CPU directs", NearForce.Elapsed());
            fprintf(timingfile,"---> %6.2f GDIPS, %6.2f Gdirects, %6.2f Mpart/sec", thistime ? gdi_cpu/thistime : 0., gdi_cpu, thistime ? JJ->NSink_CPU/thistime/1e6 : 0.);

#endif
        } // if(JJ)
    
    if (TY!=NULL) {    // Not in IC steps
        fprintf(timingfile, "\nBreakdown of Taylor Evaluate:");
        //REPORT(1, "Taylor Computation", TaylorCompute.Elapsed());
        denom = TimeStepWallClock.Elapsed();
        REPORT(1, "Taylor Computation", TaylorForce.Elapsed());
            REPORT_RATE();
        denom = thistime;
            REPORT(2, "Compute Cell Offsets", TY->ConstructOffsets.Elapsed());
            REPORT(2, "Taylor FFT", TY->FFTTaylor.Elapsed());
            REPORT(2, "Taylor R to C", TY->TaylorR2C.Elapsed());
            REPORT(2, "Taylor ASM", TY->TaylorASM.Elapsed());
            REPORT(2, "Taylor Redlack", RL->TaylorRedlack.Elapsed());
    }

    fprintf(timingfile, "\n\nBreakdown of Kick:");
    denom = TimeStepWallClock.Elapsed();
    REPORT(1, "Kick", Kick.Elapsed());
        REPORT_RATE();
    denom = thistime;
        if(JJ){
            REPORT(2, "Accumulate Pencil Stats", JJ->FinalizeTimer.Elapsed());
        }
        REPORT(2, "Add Near + Far Accel", AddAccel.Elapsed());
            fprintf(timingfile,"---> %6.2f GB/sec", thistime ? P.np/thistime*3*sizeof(accstruct)/1e9 : 0.);
        REPORT(2, "Kick Cell", KickCellTimer.Elapsed());
            REPORT_RATE();
    
    if(GFC != NULL){
        fprintf(timingfile, "\n\nBreakdown of Group Finding:");
        denom = TimeStepWallClock.Elapsed();
        REPORT(1, "Group Finding", gf_total);
            REPORT_RATE();
        denom = thistime;
            REPORT(2, "MakeCellGroups", MakeCellGroups.Elapsed());
                REPORT_RATE();
            REPORT(2, "FindCellGroupLinks", FindCellGroupLinks.Elapsed());
                REPORT_RATE();
            REPORT(2, "DoGlobalGroups", DoGlobalGroups.Elapsed());
                REPORT_RATE();
            denom = thistime;
            REPORT(3, "Create (Sort + Index + Find)", GFC->SortLinks.Elapsed() + GFC->IndexLinks.Elapsed() + GFC->FindGlobalGroupTime.Elapsed());
            REPORT(3, "Gather particles/Scatter aux", GFC->IndexGroups.Elapsed() + GFC->GatherGroups.Elapsed() + GFC->ScatterAux.Elapsed());
            REPORT(3, "Find L1 Halos", GFC->ProcessLevel1.Elapsed());
            REPORT(3, "Output L1 Halos", GFC->OutputLevel1.Elapsed());
            denom = gf_total;
            REPORT(2, "FinishGroups", FinishGroups.Elapsed());
            denom = thistime;
            REPORT(3, "Scatter Groups", GFC->ScatterGroups.Elapsed());
                fprintf(timingfile,"---> %6.2f M_group_part/sec",thistime ? GFC->L0stats.tot/thistime/1e6 : 0.);

        // Now write some detailed multiplicity and timing stats to lastrun.grouplog
        GFC->report();
    }
    
    fprintf(timingfile, "\n\nBreakdown of Output:");
    denom = TimeStepWallClock.Elapsed();
    REPORT(1, "Output", Output.Elapsed());
        REPORT_RATE();
    denom = thistime;
        REPORT(2, "TimeSlice", OutputTimeSlice.Elapsed());
            REPORT_RATE();
        REPORT(2, "LightCone", OutputLightCone.Elapsed());
            REPORT_RATE();
        REPORT(2, "Binning", OutputBin.Elapsed());
            REPORT_RATE();

    if(GFC != NULL){
        fprintf(timingfile, "\n\nBreakdown of Microstep:");
        denom = TimeStepWallClock.Elapsed();
        REPORT(1, "Microstep", Microstep.Elapsed());
            REPORT_RATE();
        denom = thistime;
        REPORT(2, "CPU Microsteps", MicrostepCPU.Elapsed());
            fprintf(timingfile,"---> %6.2f M_group_part/sec", thistime ? GFC->L0stats.tot/thistime/1e6 : 0.);
    }

    fprintf(timingfile, "\n\nBreakdown of Drift:");
    denom = TimeStepWallClock.Elapsed();
    REPORT(1, "Drift", Drift.Elapsed());
        REPORT_RATE();
    denom = thistime;
        REPORT(2, "Move",         DriftMove.Elapsed());
            REPORT_RATE();
        REPORT(2, "Rebin",        DriftRebin.Elapsed());
            REPORT_RATE();
        REPORT(2, "Collect Insert List Gaps",    DriftInsert.Elapsed());

    if(MF != NULL){
        fprintf(timingfile, "\n\nBreakdown of Compute Multipole:");
        denom = Finish.Elapsed();
        REPORT(1, "Compute Multipoles", ComputeMultipoles.Elapsed());
            REPORT_RATE();
        denom = thistime;
            REPORT(2, "Compute Cell Offsets", MF->ConstructOffsets.Elapsed());
            REPORT(2, "Multipole ASM", MF->MultipoleASM.Elapsed());
            REPORT(2, "Multipole C to R", MF->MultipoleC2R.Elapsed());
            REPORT(2, "Multipole FFT", MF->FFTMultipole.Elapsed());
    }
    
    fprintf(timingfile, "\n\nBreakdown of Finish:");
    denom = TimeStepWallClock.Elapsed();
    REPORT(1, "Finish", Finish.Elapsed());
        REPORT_RATE();
    denom = Finish.Elapsed();
        REPORT(2, "Partition Insert List", IL->FinishPartition.Elapsed());
        REPORT(2, "Sort Insert List", IL->FinishSort.Elapsed());
            fprintf(timingfile,"---> %6.2f Mitems/sec (%.2g items)", thistime ? IL->n_sorted/thistime/1e6 : 0., (double) IL->n_sorted);
        REPORT(2, "Index Cells", FinishCellIndex.Elapsed());
        REPORT(2, "Merge", FinishMerge.Elapsed());
            REPORT_RATE();
        REPORT(2, "Compute Multipoles", ComputeMultipoles.Elapsed());
            REPORT_RATE();
        REPORT(2, "Write Particles", WriteMergeSlab.Elapsed());
        REPORT(2, "Write Multipoles", WriteMultipoleSlab.Elapsed());
    
    // Misc global timings
    denom = TimeStepWallClock.Elapsed();
    REPORT(0, "\nAllocate Arena Memory", LBW->ArenaMalloc.Elapsed());
    REPORT(0, "Free Arena Memory", LBW->ArenaFreeTime());
    REPORT(0, "Free SlabAccum Variables", SlabAccumFree.Elapsed());
    
    fprintf(timingfile, "\n\nReasons for Spinning:");
    fprintf(timingfile, "\n\t Note: may add up to >100%% if there are multiple simultaneous reasons for spinning");
    denom = spinning;
    REPORT(1, "Not enough RAM to load slabs", Dependency::spin_timers[NOT_ENOUGH_RAM].Elapsed());
    REPORT(1, "Waiting for slab IO", Dependency::spin_timers[WAITING_FOR_IO].Elapsed());
    REPORT(1, "Waiting for GPU", Dependency::spin_timers[WAITING_FOR_GPU].Elapsed());
    
    return;
}
