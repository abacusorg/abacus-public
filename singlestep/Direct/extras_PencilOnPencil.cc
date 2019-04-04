// This is extra code used to check PencilOnPencil.cpp by providing
// CPU equivalents.  It is not normally executed.

void SetInteractionCollection::CPUExecute(){
    // Currently Plummer only
    // TODO: Also, this does not compute FOF densities
#ifndef DIRECTPLUMMER
    QUIT("Error: executing CPU pencils with non-Plummer softening is not currently supported.\n");
#endif
    
    int WIDTH = nfwidth;
    
    FLOAT cpu_eps = NFD->SofteningLengthInternal*NFD->SofteningLengthInternal;
    
    // This hasn't been allocated yet, so do it here.
    assert(posix_memalign((void **) &SinkSetAccelerations, 4096, sizeof(accstruct) * PaddedSinkTotal) == 0);

    // Copy the sources and sinks into position
    FillSinks.Start();
    List3<FLOAT> *SinkSetPositions = new List3<FLOAT>(NSinkBlocks*NFBlockSize);
    #pragma omp parallel for schedule(dynamic,1)
    for (int j=0; j<NSinkSets; j++) {
        SinkPlan[j].copy_into_pinned_memory(*SinkSetPositions, SinkSetStart[j], SinkSetCount[j], SinkPosSlab, nfradius, Nslab[nfradius]);
    }
    FillSinks.Stop();

    FillSources.Start();
    List3<FLOAT> *SourceSetPositions = new List3<FLOAT>(NSourceBlocks*NFBlockSize);
    #pragma omp parallel for schedule(dynamic,1)
    for (int j=0; j<NSourceSets; j++) {
        SourcePlan[j].copy_into_pinned_memory(*SourceSetPositions, SourceSetStart[j], SourceSetCount[j], SourcePosSlab, nfradius, Nslab);
    }
    FillSources.Stop();

    for(int blockIdx = 0; blockIdx < NSinkBlocks; blockIdx++){
        #pragma omp parallel for schedule(dynamic,1)
        for(int threadIdx = 0; threadIdx < NFBlockSize; threadIdx++){

            int id = NFBlockSize*blockIdx + threadIdx;
            int sinkIdx = SinkBlockParentPencil[blockIdx];

            FLOAT sinkX, sinkY, sinkZ;
            if(id <  SinkSetIdMax[sinkIdx]){
                sinkX = SinkSetPositions->X[id];
                sinkY = SinkSetPositions->Y[id];
                sinkZ = SinkSetPositions->Z[id];
            } else {
                 sinkX =0;
                 sinkY =0;
                 sinkZ =0;
            }


            FLOAT3 a(0,0,0);

            int InteractionStart = sinkIdx*WIDTH;
            int InteractionMax = InteractionStart + WIDTH;

            for(int c = InteractionStart; c < InteractionMax; c++){
                int sourceIdx = SinkSourceInteractionList[c];
                int sourceStart = SourceSetStart[sourceIdx];
                int sourceCount = SourceSetCount[sourceIdx];
                FLOAT yoffset = SinkSourceYOffset[c];
                assert(sourceCount > 0);
                int nB = sourceCount/NFBlockSize;

                for(int b = 0; b < nB; b+=1){
                    int idx = sourceStart + b*NFBlockSize;

                    for(int i = 0; i < NFBlockSize; i++){
                        FLOAT sourceX = SourceSetPositions->X[idx + i];
                        FLOAT sourceY = SourceSetPositions->Y[idx + i];
                        FLOAT sourceZ = SourceSetPositions->Z[idx + i];

                        FLOAT drx, dry, drz, r;
                        drx = sourceX - sinkX;
                        dry = sourceY - sinkY + yoffset;
                        drz = sourceZ - sinkZ;

                        r = (FLOAT)1./std::sqrt( drx*drx + dry*dry + drz*drz  + cpu_eps);
                        r *=r*r;
                        a.x -= r * drx;
                        a.y -= r * dry;
                        a.z -= r * drz;
                    }
                }

                int remaining = sourceCount%NFBlockSize;

                //if(threadIdx < remaining){
                    int idx = sourceStart + nB*NFBlockSize;

                    for(int i = 0; i < remaining; i++){
                        FLOAT sourceX = SourceSetPositions->X[idx + i];
                        FLOAT sourceY = SourceSetPositions->Y[idx + i];
                        FLOAT sourceZ = SourceSetPositions->Z[idx + i];

                        FLOAT drx, dry, drz, r;
                        drx = sourceX - sinkX;
                        dry = sourceY - sinkY + yoffset;
                        drz = sourceZ - sinkZ;

                        r = (FLOAT)1./std::sqrt( drx*drx + dry*dry + drz*drz  + cpu_eps);
                        r *=r*r;
                        a.x -= r * drx;
                        a.y -= r * dry;
                        a.z -= r * drz;
			// TODO: Need to compute the FOFneighbor count
                    }
                //}

            }

            if(id < SinkSetIdMax[sinkIdx]){
                assert(isfinite(a.x));
                assert(isfinite(a.y));
                assert(isfinite(a.z));
                assert(id < PaddedSinkTotal);
                SinkSetAccelerations[id].x = a.x;
                SinkSetAccelerations[id].y = a.y;
                SinkSetAccelerations[id].z = a.z;
		#ifdef COMPUTE_FOF_DENSITY
		SinkSetAccelerations[id].w = 0.0;
		#endif
            }
        }
    }

    //for(int i = 0; i < PaddedSinkTotal; i++)
    //    assert(TOFLOAT3(SinkSetAccelerations[i]).is_finite());

    // Need to copy back the SinkSetAccel
    CopyAccelFromPinned.Start();
    #pragma omp parallel for schedule(dynamic,1)
    for (int j=0; j<j_width; j++) {
	// We need each Y skewer to be done by one thread
	for (int k=0; k<cpd; k++) {
	    int idx = j*Nk+k;
	    SinkPlan[idx].copy_from_pinned_memory((void *)SinkSetAccelerations, SinkSetStart[idx], SinkSetCount[idx], (void *)SinkAccSlab, idx, nfradius, Nslab[nfradius]);
	}

    }
    CopyAccelFromPinned.Stop();

    SetCompleted();
    free(SinkSetAccelerations);
    delete SinkSetPositions;
    delete SourceSetPositions;
}





/// This routine removed.  Old code in 1bde1c52074dcfdf6355f03ea3c11cdff3f68206
void SetInteractionCollection::PrintInteractions(){
}

/// This routine removed.  Old code in 1bde1c52074dcfdf6355f03ea3c11cdff3f68206
void SetInteractionCollection::AddInteractionList( std::vector<uint64> ** il){
}
