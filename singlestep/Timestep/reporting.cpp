/* reporting.cpp
 *
 * This writes timing outputs.  It should provide a nicely formatted
 * summary of the timings from the various portions of the code.
 *
 * Note: we use the decimal convention for bytes; i.e. 1 KB = 1000 B
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

#define REPORT_BUFFER_SIZE (sizeof(char) * 256*1024)

char *reportbuffer;
FILE *reportfp;

void InitializeReport() {
    reportbuffer = (char *) malloc(REPORT_BUFFER_SIZE);  //  allocate a 256 KB string buffer for the timings file
    reportfp = fmemopen(reportbuffer, REPORT_BUFFER_SIZE, "w");
    return;
}

#define REPORT(tabs, str, a) \
     do { \
         fprintf(reportfp, "\n"); \
         for (int _i=0; _i<tabs; _i++) fprintf(reportfp, "    "); \
         thistime = a; fprintf(reportfp, "%-30s: %10.2e sec (%5.2f%%) ", str, thistime, denom ? 100*thistime/denom : 0.); \
     } while(0)
		 
#define REPORT_RATE(dependency) \
     do { \
        fprintf(reportfp,"---> %6.2f Mpart/sec", thistime ? dependency->num_particles/thistime/1e6 : 0.); \
     } while(0)		 


std::string cpu_model_str(){
    // https://stackoverflow.com/a/24957090
    int cpuInfo[4] = {-1};
    char CPUBrandString[0x40];

    memset(CPUBrandString, 0, sizeof(CPUBrandString));

    __cpuid(cpuInfo, 0x80000002);
    memcpy(CPUBrandString, cpuInfo, sizeof(cpuInfo));

    __cpuid(cpuInfo, 0x80000003);
    memcpy(CPUBrandString + 16, cpuInfo, sizeof(cpuInfo));

    __cpuid(cpuInfo, 0x80000004);
    memcpy(CPUBrandString + 32, cpuInfo, sizeof(cpuInfo));

    return std::string(CPUBrandString);
}


/* This function gathers timings from major global classes before they're destroyed in the Epilogue.
 * They're printed to a string buffer (via a file-like interface); this string buffer saves a few spots
 * for teardown timings that happen after this function runs.  Those spots get filled in in ReportTimings(),
 * which runs after the Epilogue.
 */
void GatherTimings() {
    InitializeReport();
    
    if(NFD)
        NFD->AggregateStats();

    double thistime, denom, total;
    denom = WallClockDirect.Elapsed();
    REPORT(0, "Total Wall Clock Time", WallClockDirect.Elapsed()); 

    // What particle count do we use to report the overall rate?  FinishParticles works for IC and normal steps.
    int64 np_node = FinishParticles->num_particles;
    int64 np_node_with_ghost = FinishParticles->num_particles_with_ghost;
    fprintf(reportfp, "---> %6.3f Mpart/sec\n", thistime ? np_node/thistime/1e6 : 0.);
    fprintf(reportfp, "                                       % " PRId64 " particles ", np_node);
    double ghostfrac = (double)(np_node_with_ghost - np_node) / np_node;
    if(MPI_size_z > 1) fprintf(reportfp,"(+%.3g%% ghost) ", ghostfrac*100);
	
    //TODO : consider reporting number of particles microstepped here as well. 
    fprintf(reportfp,"finished by this node.\n\n");

    int NGPU = NFD ? NFD->NGPU : 0;

    fprintf(reportfp, "| %s (%d^3, CPD=%d); ", P.SimName, (int) roundf(P.ppd()), P.cpd);
#ifdef PARALLEL
    fprintf(reportfp, "%dx%d MPI ranks\n", MPI_size_x, MPI_size_z);
#else
    fprintf(reportfp, "no MPI\n");
#endif
    fprintf(reportfp, "| %d x %s\n", omp_get_num_procs(), cpu_model_str().c_str());
    if(NFD) fprintf(reportfp, "| %d x %s\n", NGPU, NFD->GPUName.c_str());

    total = 0.0;
    REPORT(0, "SingleStep Setup", SingleStepSetup.Elapsed()); total += thistime;
#ifdef PARALLEL
    REPORT(0, "Convolution", ConvolutionWallClock.Elapsed()); total += thistime;
#endif
    REPORT(0, "TimeStep", TimeStepWallClock.Elapsed()); total += thistime;
    REPORT(0, "Finish IO", IOFinish.Elapsed()); total += thistime;
    REPORT(0, "Unaccounted", WallClockDirect.Elapsed()-total);
    fprintf(reportfp, "\n");

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
                fprintf(reportfp, "---> %6.1f MB/s " #RW " on %6.2f GB ["  #BLOCKING "]", this_io_time ? this_io_bytes/this_io_time/1e6 : 0., this_io_bytes/1e9);\
            }\
        }\
    } while(0)

    for(int i = 0; i < niothreads; i++){
        // TODO: could probably move some of this inside the IO thread destructor
        double total_read_time = 0., total_read_bytes = 0.;
        double total_write_time = 0., total_write_bytes = 0.;

        ACCUMULATE_THREAD_TOTALS(BlockingIORead, read);
        ACCUMULATE_THREAD_TOTALS(BlockingIOWrite, write);
        ACCUMULATE_THREAD_TOTALS(NonBlockingIORead, read);
        ACCUMULATE_THREAD_TOTALS(NonBlockingIOWrite, write);

        // Though not a measure of disk IO performance, the checksumming does affect how fast we can write data
        total_write_time += ChecksumTime[i];

        double total_time = total_read_time + total_write_time;
        double total_bytes = total_read_bytes + total_write_bytes;

        denom = WallClockDirect.Elapsed();
#ifdef IOTHREADED
        std::string threadname = "IO Thread " + std::to_string(i+1ll);
        REPORT(0, threadname.c_str(), total_time);
#else
        REPORT(0, "Blocking IO", total_time);
#endif
        fprintf(reportfp, "---> %6.1f MB/s read, %6.1f MB/s write", total_read_time ? total_read_bytes/total_read_time/1e6 : 0., total_write_time ? total_write_bytes/total_write_time/1e6 : 0.);
        denom = thistime;
        REPORT_DIR_IOSTATS(BlockingIORead, read, blocking);
        REPORT_DIR_IOSTATS(BlockingIOWrite, write, blocking);
        REPORT_DIR_IOSTATS(NonBlockingIORead, read, non-blocking);
        REPORT_DIR_IOSTATS(NonBlockingIOWrite, write, non-blocking);
        if(ChecksumBytes[i]){
            REPORT(1, "Calculate Checksum", ChecksumTime[i]);
            fprintf(reportfp, "---> %6.1f MB/s on %6.2f GB", ChecksumTime[i] ? ChecksumBytes[i]/ChecksumTime[i]/1e6 : 0., ChecksumBytes[i]/1e9);
        }
    }
#undef niothreads

    // Total time checking for MPI completion and freeing buffers
    double manifest_check_time = 0.;
    double RManifestTime = 0.;
    double SManifestTime = 0.;
    double total_mpi_check = 0.;  // things not timed in other actions

#ifdef PARALLEL
    for(int i = 0; i < nManifest; i++){
        manifest_check_time += _SendManifest[i].CheckCompletion.Elapsed() + 
                                _ReceiveManifest[i].CheckCompletion.Elapsed();
        RManifestTime += _ReceiveManifest[i].Load.Elapsed() +
                         _ReceiveManifest[i].Transmit.Elapsed();
        SManifestTime += _SendManifest[i].Load.Elapsed() +
                         _SendManifest[i].Transmit.Elapsed();  // in Finish
    }
#endif

    fprintf(reportfp, "\n\nBreakdown of TimeStep: ");
    total = 0.0;
    denom = TimeStepWallClock.Elapsed();
    REPORT(1, "FetchSlabs", FetchSlabs->Elapsed()); total += thistime;

    if(NFD){
        REPORT(1, "Transpose Positions", TransposePos->Elapsed()); total += thistime;

    #ifdef CUDADIRECT
        REPORT(1, "NearForce [blocking]", NearForce->Elapsed()); total += thistime;
        REPORT(1, "NearForce [non-blocking]", NFD->GPUThroughputTime);
        fprintf(reportfp,"---> %6.2f effective GDIPS, %.2f Mpart/sec", thistime ? NFD->gdi_gpu/thistime : 0, thistime ? NearForce->num_particles/thistime/1e6 : 0.);
        double total_di = (NFD->DirectInteractions_CPU +NFD->TotalDirectInteractions_GPU)/1e9;
    #else
        REPORT(1, "NearForce", NearForce->Elapsed()); total += thistime;
        REPORT_RATE(NearForce);
    #endif
    }

    if(TY){
        double ty_total = TaylorTranspose->Elapsed() + TaylorForce->Elapsed();
        REPORT(1, "Taylors", ty_total); total += thistime;
        REPORT_RATE(TaylorForce);
        REPORT(1, "Kick", Kick->Elapsed()); total += thistime;
        REPORT_RATE(Kick);
    }
    
    double gf_total = 0, gf_precon = 0;
    if (GFC != NULL){
        gf_total = MakeCellGroups->Elapsed() + FindCellGroupLinks->Elapsed() + DoGlobalGroups->Elapsed() + FinishGroups->Elapsed();
        gf_precon = MakeCellGroups->ElapsedPrecon() + FindCellGroupLinks->ElapsedPrecon() + DoGlobalGroups->ElapsedPrecon() + FinishGroups->ElapsedPrecon();
    
        REPORT(1, "Group Finding", gf_total); total += thistime;
        REPORT_RATE(DoGlobalGroups);
    }

    if(Output){
        REPORT(1, "Output", Output->Elapsed()); total += thistime;
        REPORT_RATE(Output);
        if (GFC != NULL){
            REPORT(1, "Microstep", Microstep->Elapsed()); total += thistime;
            REPORT_RATE(Microstep);
        }
    }
   
    if(WriteState.Do2LPTVelocityRereading){
        REPORT(1, "Unpack IC Velocities for LPT", UnpackLPTVelocity->Elapsed()); total += thistime;
        REPORT_RATE(UnpackLPTVelocity);
    }
    REPORT(1, "Drift", Drift->Elapsed()); total += thistime;
        REPORT_RATE(Drift);

    double neighbor_tot;
    if(NeighborSend && MPI_size_z > 1){
        neighbor_tot = NeighborSend->Elapsed() + DoNeighborRecv->Elapsed();
        REPORT(1, "Neighbor Exchange", neighbor_tot);
        total += neighbor_tot;
    }

    double finish_total = FinishParticles->Elapsed() + FinishMultipoles->Elapsed();
    REPORT(1, "Finish", finish_total); total += thistime;
        REPORT_RATE(FinishParticles);

#ifdef PARALLEL
    double multipoles_mpi_check = CheckForMultipoles->Elapsed();
    if(Check2DMultipoleMPI) multipoles_mpi_check += Check2DMultipoleMPI->Elapsed();
    if(Check2DTaylorMPI) multipoles_mpi_check += Check2DTaylorMPI->Elapsed();
    REPORT(1, "Check Multipoles MPI", multipoles_mpi_check); total += thistime;
#endif


#ifdef PARALLEL
    double manifest_total = DoReceiveManifest->Elapsed() + DoSendManifest->Elapsed();
    REPORT(1, "Manifest", manifest_total); total += thistime;
    REPORT(1, "MPI Barrier", BarrierWallClock.Elapsed()); total += thistime;
#endif 
    
    REPORT(1, "Spinning", Dependency::SpinTimer.Elapsed()); total += thistime;
    REPORT(1, "Unaccounted", TimeStepWallClock.Elapsed()-total);

    // breakdown per-slab
    double slabforcetimemean = 0;
    double slabforcetimesigma = 0;
    double slabforcemaxtime = 0;
    double slabforcemintime =  1e9;
    double slabforcelatencymean = 0;
    double slabforcelatencysigma = 0;
    double slabforcemaxlatency = 0;
    double slabforceminlatency = 1e9;
    char *slabtimesbuffer = NULL;

    if(NFD){
        //char fn[1024];
        //int ret = snprintf(fn, 1024, "%s/lastrun%s.slabtimes",P.LogDirectory, NodeString);
        //assert(ret >= 0 && ret < 1024);
        //FILE* slabtimefile = fopen(fn,"wb");
        slabtimesbuffer = (char *) malloc(REPORT_BUFFER_SIZE);  //  allocate a 128 KB string buffer for the timings file
        FILE *slabtimesfp = fmemopen(slabtimesbuffer, REPORT_BUFFER_SIZE, "w");
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
            if (slabtime>0.0) 
                fprintf(slabtimesfp, "%4d %#10.4f %#10.4f\n", i, slabtime, slablatency);

            //fwrite(&slabtime,sizeof(double),1,slabtimefile);
            //fwrite(&slablatency,sizeof(double),1,slabtimefile);
        }
        //fclose(slabtimefile);
        fclose(slabtimesfp);
    }
    slabforcetimesigma =  sqrt(slabforcetimesigma    - slabforcetimemean*slabforcetimemean);
    slabforcelatencysigma=sqrt(slabforcelatencysigma - slabforcelatencymean*slabforcelatencymean);

    /*
    fprintf(reportfp, "\n\nBreakdown per slab (Wall Clock)");
    denom = TimeStepWallClock.Elapsed()/P.cpd;
    REPORT(1,"Mean Force Computation",slabforcetimemean);
    fprintf(reportfp,"\n\t\tSigma(STD): %.2g s\t Min: %.2e s\t Max: %.2e s ",slabforcetimesigma,slabforcemintime,slabforcemaxtime);
    REPORT(1,"Mean Force Latency",slabforcelatencymean);
    fprintf(reportfp,"\n\t\tSigma(STD): %.2e s\t Min: %.2e s\t Max: %.2e s ",slabforcelatencysigma,slabforceminlatency,slabforcemaxlatency);*/
    
    if(NFD){
        denom = TimeStepWallClock.Elapsed();
        fprintf(reportfp, "\n\nBreakdown of Near Force:");
        double gdi_cpu = NFD->DirectInteractions_CPU/1e9;  // Measure per-core load balancing?
#ifdef CUDADIRECT
        /*fprintf(reportfp, "\n\tNotes about non-blocking timing:\n");
        fprintf(reportfp, "\t- \"Directs Throughput\" is the wall clock time while at least one GPU thread is running (copy or compute).\n");
        fprintf(reportfp, "\t- \"Effective\" GDIPS is based on this throughput.\n");*/
        REPORT(1, "Blocking", NearForce->Elapsed());
            denom = NearForce->Elapsed();
            REPORT(2, "Precondition", NearForce->ElapsedPrecon());
            REPORT(2, "Calculate Splits", NFD->CalcSplitDirects.Elapsed());
            REPORT(2, "Construct Pencils", NFD->SICConstruction.Elapsed());
            REPORT(2, "Dispatch Interaction", NFD->SICExecute.Elapsed());
            REPORT(2, "CPU Fallback", NFD->CPUFallbackTimer.Elapsed());
                fprintf(reportfp,"---> %6.2f GDIPS, %6.2f Gdirects, %6.2f Mpart/sec\n", thistime ? gdi_cpu/thistime : 0., gdi_cpu, thistime ? NFD->NSink_CPU/thistime/1e6 : 0.);
        
        denom = TimeStepWallClock.Elapsed();
        char str[1024];  sprintf(str, "Non-Blocking (thread-seconds, %d threads)", NFD->NBuffers);
        REPORT(1, str, NFD->DeviceThreadTimer);
            denom = NFD->DeviceThreadTimer;
            REPORT(2, "Fill Sinks", NFD->FillSinks);  // filling is unpadded
                    fprintf(reportfp,"---> %6.1f MB/s, %6.2f MSink/sec", thistime ? NFD->total_sinks*sizeof(posstruct)/1e6/thistime : 0., thistime ? NFD->total_sinks/1e6/thistime : 0.);
            REPORT(2, "Fill Sources", NFD->FillSources);
                    fprintf(reportfp,"---> %6.1f MB/s, %6.2f MSource/sec", thistime ? NFD->total_sources*sizeof(posstruct)/1e6/thistime : 0., thistime ? NFD->total_sources/1e6/thistime : 0.);
            REPORT(2, "Launch Kernels", NFD->LaunchDeviceKernels);
            REPORT(2, "Wait for GPU Result", NFD->WaitForResult);
            REPORT(2, "Copy Accel from Pinned", NFD->CopyAccelFromPinned);
                    fprintf(reportfp,"---> %6.1f MB/s, %6.2f MSink/sec\n", thistime ? NFD->total_sinks*sizeof(accstruct)/1e6/thistime : 0., thistime ? NFD->total_sinks/1e6/thistime : 0.);  // same number of accels as sinks
            
        denom = TimeStepWallClock.Elapsed();
        REPORT(1, "Non-Blocking Throughput (Wall Clock)", NFD->GPUThroughputTime);
                denom = NFD->GPUThroughputTime;
                fprintf(reportfp,"\n\t\t\t\t---> %6.2f effective GDIPS, %6.2f Mpart/sec, %6.2f Msink/sec", thistime ? NFD->gdi_gpu/thistime : 0., thistime ? NearForce->num_particles/thistime/1e6 : 0., thistime ? NFD->total_sinks/thistime/1e6 : 0.);
                fprintf(reportfp,"\n\t\t\t\t---> %6.2f Gdirects, %6.2f padded Gdirects", NFD->gdi_gpu, NFD->gdi_padded_gpu);
                fprintf(reportfp,"\n\t\t\t\t---> with %d device threads, estimate %.1f%% thread concurrency", NFD->NBuffers, (NFD->DeviceThreadTimer - NFD->GPUThroughputTime)/(NFD->DeviceThreadTimer - NFD->DeviceThreadTimer/NFD->NBuffers)*100);
                
            
            /*// Such detailed stats are not very useful
            fprintf(reportfp, "\n    Device stats:\n");
            for(int g = 0; g < NFD->NBuffers; g++){
                fprintf(reportfp, "        Device thread %d (GPU %d):", g, g % NGPU);
                fprintf(reportfp, " %.2f GB to device, %.2f GB from device, %.2f Msink, %.2f Gdirects\n", NFD->GB_to_device[g], NFD->GB_from_device[g], NFD->DeviceSinks[g]/1e6, NFD->DirectInteractions_GPU[g]/1e9);
            }*/

            // max(directs) - mean(directs) is an estimate of the load balancing
            double maxdirects = *std::max_element(NFD->DirectInteractions_GPU, NFD->DirectInteractions_GPU + NFD->NBuffers);
            fprintf(reportfp,"\n\t\t\t\t---> directs load imbalance is %.3g%%\n", (maxdirects*NFD->NBuffers - 1e9*NFD->gdi_gpu)/(1e9*NFD->gdi_gpu)*100);
#else
        REPORT(1, "CPU directs", NearForce->Elapsed());
            fprintf(reportfp,"---> %6.2f GDIPS, %6.2f Gdirects, %6.2f Mpart/sec", thistime ? gdi_cpu/thistime : 0., gdi_cpu, thistime ? NFD->NSink_CPU/thistime/1e6 : 0.);

#endif
        } // if(NFD)
    
    if (TY!=NULL) {    // Not in IC steps
        fprintf(reportfp, "\nBreakdown of Taylors:");
        //REPORT(1, "Taylor Computation", TaylorCompute.Elapsed());
        denom = TimeStepWallClock.Elapsed();
        REPORT(1, "Taylor Computation", TaylorForce->Elapsed());
            REPORT_RATE(TaylorForce);
        denom = thistime;
            REPORT(2, "Precondition", TaylorForce->ElapsedPrecon());
            REPORT(2, "Compute Cell Offsets", TY->ConstructOffsets.Elapsed());
            if(MPI_size_z > 1) REPORT(2, "Unpack MPI Buffer", TY->UnpackRecvBuf.Elapsed());
            REPORT(2, "Taylor FFT", TY->FFTTaylor.Elapsed());
            REPORT(2, "Taylor R to C", TY->TaylorR2C.Elapsed());
            REPORT(2, "Taylor Kernel", TY->TaylorKernel.Elapsed());
            REPORT(2, "Taylor Redlack", RL->TaylorRedlack.Elapsed());

        if(MPI_size_z > 1){
            denom = TimeStepWallClock.Elapsed();
            REPORT(1, "Taylor Transpose", TaylorTranspose->Elapsed());
            denom = thistime;
                REPORT(2, "Precondition", TaylorTranspose->ElapsedPrecon());
                REPORT(2, "Taylor Z-FFT", TY->FFTZTaylor.Elapsed());
                REPORT(2, "Fill MPI Buffers", TY->FillMPIBufs.Elapsed());
                REPORT(2, "Launch MPI All-to-all", TY->AllToAll.Elapsed());
        }
    }

    if(Kick){
        fprintf(reportfp, "\n\nBreakdown of Kick:");
        denom = TimeStepWallClock.Elapsed();
        REPORT(1, "Kick", Kick->Elapsed());
            REPORT_RATE(Kick);
        denom = thistime;
            REPORT(2, "Precondition", Kick->ElapsedPrecon());
            if(NFD){
                REPORT(2, "Accumulate Pencil Stats", NFD->FinalizeTimer.Elapsed());
            }
            REPORT(2, "Add Near + Far Accel", AddAccel.Elapsed());
                fprintf(reportfp,"---> %6.2f GB/sec", thistime ? Kick->num_particles/thistime*3*sizeof(accstruct)/1e9 : 0.);
            REPORT(2, "Kick Cell", KickCellTimer.Elapsed());
                REPORT_RATE(Kick);
        
        if(GFC != NULL){
            fprintf(reportfp, "\n\nBreakdown of Group Finding:");
            denom = TimeStepWallClock.Elapsed();
            REPORT(1, "Group Finding", gf_total);
                REPORT_RATE(DoGlobalGroups);
            denom = thistime;
                REPORT(2, "Preconditions", gf_precon);
                REPORT(2, "MakeCellGroups", MakeCellGroups->Elapsed());
                    REPORT_RATE(MakeCellGroups);
                REPORT(2, "FindCellGroupLinks", FindCellGroupLinks->Elapsed());
                    REPORT_RATE(FindCellGroupLinks);
                REPORT(2, "DoGlobalGroups", DoGlobalGroups->Elapsed());
                    REPORT_RATE(DoGlobalGroups);
                denom = thistime;
                REPORT(3, "Create (Sort + Index + Find)", GFC->SortLinks.Elapsed() + GFC->IndexLinks.Elapsed() + GFC->FindGlobalGroupTime.Elapsed());
                REPORT(3, "Gather particles/Scatter aux", GFC->IndexGroups.Elapsed() + GFC->GatherGroups.Elapsed() + GFC->ScatterAux.Elapsed());
                REPORT(3, "Find L1 Halos", GFC->ProcessLevel1.Elapsed());
                REPORT(3, "Output L1 Halos", GFC->OutputLevel1.Elapsed());
                denom = gf_total;
                REPORT(2, "FinishGroups", FinishGroups->Elapsed());
                denom = thistime;
                REPORT(3, "Scatter Groups", GFC->ScatterGroups.Elapsed());
                    fprintf(reportfp,"---> %6.2f M_group_part/sec",thistime ? GFC->L0stats.tot/thistime/1e6 : 0.);
        }
    }
    
    if(Output){
        fprintf(reportfp, "\n\nBreakdown of Output:");
        denom = TimeStepWallClock.Elapsed();
        REPORT(1, "Output", Output->Elapsed());
            REPORT_RATE(Output);
        denom = thistime;
            REPORT(2, "Precondition", Output->ElapsedPrecon());
            REPORT(2, "Output Time Slice", OutputTimeSlice.Elapsed());
                REPORT_RATE(Output);
            REPORT(2, "Output Light Cone", OutputLightCone.Elapsed());
                REPORT_RATE(Output);
            REPORT(3, "Output Light Cone Search", OutputLightConeSearch.Elapsed());
                REPORT_RATE(Output);
            REPORT(2, "Output Bin", OutputBin.Elapsed());
                REPORT_RATE(Output);
    }

    if(GFC != NULL){
        fprintf(reportfp, "\n\nBreakdown of Microstep:");
        denom = TimeStepWallClock.Elapsed();
        REPORT(1, "Microstep", Microstep->Elapsed());
            REPORT_RATE(Microstep);
        denom = thistime;
        REPORT(2, "Precondition", Microstep->ElapsedPrecon());
        REPORT(2, "CPU Microsteps", MicrostepCPU.Elapsed());
            fprintf(reportfp,"---> %6.2f M_group_part/sec", thistime ? GFC->L0stats.tot/thistime/1e6 : 0.);
    }

    fprintf(reportfp, "\n\nBreakdown of Drift:");
    denom = TimeStepWallClock.Elapsed();
    REPORT(1, "Drift", Drift->Elapsed());
        REPORT_RATE(Drift);
    denom = thistime;
        REPORT(2, "Precondition", Drift->ElapsedPrecon());
        REPORT(2, "Move",         DriftMove.Elapsed());
        REPORT(2, "Rebin",        DriftRebin.Elapsed());

#ifdef PARALLEL
    if(NeighborSend && MPI_size_z > 1){
        fprintf(reportfp, "\n\nBreakdown of Neighbor Exchange:");
        denom = TimeStepWallClock.Elapsed();
        REPORT(1, "Neighbor Send", NeighborSend->Elapsed());
        REPORT(1, "Neighbor Receive", DoNeighborRecv->Elapsed());
    }
#endif

    if(MF != NULL){
        fprintf(reportfp, "\n\nBreakdown of Compute Multipole:");
        denom = FinishParticles->Elapsed();
        REPORT(1, "Compute Multipoles", ComputeMultipoles.Elapsed());
            REPORT_RATE(FinishParticles);
        denom = thistime;
            REPORT(2, "Compute Cell Offsets", MF->ConstructOffsets.Elapsed());
            REPORT(2, "Multipole Kernel", MF->MultipoleKernel.Elapsed());
            REPORT(2, "Multipole C to R", MF->MultipoleC2R.Elapsed());
            REPORT(2, "Multipole FFT", MF->FFTMultipole.Elapsed());

            if(MPI_size_z > 1){
                REPORT(2, "Fill MPI Buffers", MF->FillMPIBufs.Elapsed());
                REPORT(2, "Launch MPI All-to-all", MF->AllToAll.Elapsed());
            }
    }
    
    fprintf(reportfp, "\n\nBreakdown of Finish:");
    denom = TimeStepWallClock.Elapsed();
    REPORT(1, "Finish Particles", FinishParticles->Elapsed());
        REPORT_RATE(FinishParticles);
    denom = thistime;
        REPORT(2, "Precondition", FinishParticles->ElapsedPrecon());
		REPORT(2, "Free ReadState Slabs", FinishFreeSlabs.Elapsed());
        REPORT(2, "Partition Insert List", IL->FinishPartition.Elapsed());
        REPORT(2, "Sort Insert List", IL->FinishSort.Elapsed());
            fprintf(reportfp,"---> %6.2f Mitems/sec (%.2g items)", thistime ? IL->n_sorted/thistime/1e6 : 0., (double) IL->n_sorted);
        REPORT(2, "Index Cells", FinishCellIndex.Elapsed());
        REPORT(2, "Merge", FinishMerge.Elapsed());
        	REPORT_RATE(FinishParticles);
        REPORT(2, "Compute Multipoles", ComputeMultipoles.Elapsed());
        	REPORT_RATE(FinishParticles);
        REPORT(2, "Write Particles", WriteMergeSlab.Elapsed());
#ifdef PARALLEL
        REPORT(2, "Queuing Send Manifest", SManifestTime);
#endif
    denom = TimeStepWallClock.Elapsed();
    REPORT(1, "Finish Multipoles", FinishMultipoles->Elapsed());
    denom = thistime;
        REPORT(2, "Precondition", FinishMultipoles->ElapsedPrecon());
    #ifdef PARALLEL
        if (MPI_size_z > 1){
            REPORT(2, "Unpack MPI Buffer", MF->UnpackRecvBuf.Elapsed());
            REPORT(2, "Multipole Z-FFT", MF->FFTZMultipole.Elapsed());
            REPORT(2, "Transpose into Slab", MF->FFTZTranspose.Elapsed());
        }
        REPORT(2, "Queuing Multipole MPI", QueueMultipoleMPI.Elapsed());
    #endif
    REPORT(2, "Write Multipoles", WriteMultipoleSlab.Elapsed());
    REPORT(2, "Release Free Memory To Kernel", ReleaseFreeMemoryTime.Elapsed());
#ifdef PARALLEL
    denom = TimeStepWallClock.Elapsed();
    REPORT(1, "Free MPI Multipoles", CheckForMultipoles->Elapsed());
#endif


#ifdef PARALLEL
    fprintf(reportfp, "\n\nBreakdown of various MPI work:");
    denom = TimeStepWallClock.Elapsed();
    REPORT(1, "Check Manifest MPI", manifest_check_time);
    REPORT(1, "Queuing Receive Manifest", RManifestTime);
    REPORT(1, "Queuing Send Manifest [Finish]", SManifestTime);
    REPORT(1, "Queuing Multipole MPI [Finish]", QueueMultipoleMPI.Elapsed());
    REPORT(1, "Check MPI Multipoles [precon]", CheckForMultipoles->ElapsedPrecon());
    REPORT(1, "Free MPI Multipoles", CheckForMultipoles->Elapsed());
    if(MPI_size_z > 1 && Check2DMultipoleMPI) REPORT(1, "Check 2D Multiple MPI", Check2DMultipoleMPI->Elapsed());
    if(MPI_size_z > 1 && Check2DTaylorMPI) REPORT(1, "Check 2D Taylor MPI", Check2DTaylorMPI->Elapsed());

    /*
    // TOOD: can't use *Manifest->, need to sum up _*Manifest
    fprintf(reportfp, "\n\nBreakdown of Manifest:");
    denom = TimeStepWallClock.Elapsed();
    REPORT(1, "Receive Manifest", RManifestTime);
    denom = RManifestTime;
    //REPORT(2, "SendManifest Check (spinning)", SendManifest->CheckCompletion.Elapsed());
    REPORT(2, "RecvManifest Load", ReceiveManifest->Load.Elapsed());
    REPORT(2, "RecvManifest Transmit", ReceiveManifest->Transmit.Elapsed());
    //REPORT(2, "RecvManifest Check (spinning)", ReceiveManifest->CheckCompletion.Elapsed()); 
    */
#endif
    
    // Misc global timings
    double spinning = Dependency::SpinTimer.Elapsed();  // `global_spin_timer2` may overlap non-action work; avoid
    fprintf(reportfp, "\n\nBreakdown of Spinning:");
    fprintf(reportfp, "\n\t Note: may add up to >100%% if there are multiple simultaneous reasons for spinning");
    denom = TimeStepWallClock.Elapsed();
    REPORT(1, "Spinning", spinning);
    denom = spinning;
    REPORT(2, "Checking preconditions", SlabDependency::global_precon_timer.Elapsed());
    REPORT(2, "Not enough RAM to load slabs", SlabDependency::spin_timers[NOT_ENOUGH_RAM].Elapsed());
    REPORT(2, "Waiting for slab IO", SlabDependency::spin_timers[WAITING_FOR_IO].Elapsed());
    REPORT(2, "Waiting for GPU", SlabDependency::spin_timers[WAITING_FOR_GPU].Elapsed());
#ifdef PARALLEL
    REPORT(2, "Waiting for MPI", SlabDependency::spin_timers[WAITING_FOR_MPI].Elapsed());
#endif
	

    denom = TimeStepWallClock.Elapsed();
    double arena_malloc, arena_free, munmap_thread;
    SB->GetMallocFreeTimes(&arena_malloc, &arena_free, &munmap_thread);

    REPORT(0, "\nAllocate Arena Memory", arena_malloc);
    REPORT(0, "Free Arena Memory", arena_free);
    REPORT(0, "Munmap Thread", munmap_thread);
    REPORT(0, "Free SlabAccum Variables", SlabAccumFree.Elapsed());

    fprintf(reportfp,"\n\nMinCellSize = %d, MaxCellSize = %d, RMS Fractional Overdensity = %10.4e\n",
        WriteState.MinCellSize, WriteState.MaxCellSize, WriteState.StdDevCellSize);
    fprintf(reportfp,"Rms |v| in simulation is %f.\n", WriteState.RMS_Velocity);

    if (GFC!=NULL) {
        // Now write some detailed multiplicity and timing stats to lastrun.grouplog
        fprintf(reportfp, "\n\n========================================================================\n\n");
        GFC->report(reportfp);
    }

    #ifdef PARALLEL
    if (convtimebuffer!=NULL) {
        fprintf(reportfp, "\n\n========================================================================\n\n");
        fputs(convtimebuffer, reportfp);
        free(convtimebuffer);
    }
    #endif

    if (NFD) {
        fprintf(reportfp, "\n\n========================================================================\n\n");
        fprintf(reportfp, "GPU Timings\nSlab   Time   Latency\n");
        fputs(slabtimesbuffer, reportfp);
        free(slabtimesbuffer);
    }
	
}

/* This function writes the timing report to disk.
 * It runs after the Epilogue has destroyed all the global objects; all the global object timings were collected in GatherTimings.
 * So we just have to fill in how long the Epilogue/teardown itself took.
 */
void ReportTimings(){
    // Now add a row with the teardown timings
    double denom, thistime;
    denom = TimeStepWallClock.Elapsed();

    fprintf(reportfp, "\n");
    REPORT(0, "SingleStep TearDown", SingleStepTearDown.Elapsed());
    fprintf(reportfp, " [not included in Total Wall Clock Time]\n");

    // Not actually closing the file; just the file-like interface to the string buffer
    fclose(reportfp);

    // and write the whole buffer to disk
    char timingfn[1050];
    sprintf(timingfn,"%s/step%04d%s.time",
        WriteState.LogDirectory, WriteState.FullStepNumber, NodeString);
    FILE *timingfp = fopen(timingfn,"w");
    assertf(timingfp != NULL, "Couldn't open timing file \"%s\"\n", timingfn);
    fputs(reportbuffer, timingfp);

    free(reportbuffer);
    fclose(timingfp);

    STDLOG(0, "Wrote timings to %s\n", timingfn);
}
