// This is extra code used to check PencilOnPencil.cpp by providing
// CPU equivalents.  It is not normally executed.

void SetInteractionCollection::CPUExecute(){
    // Currently Plummer only
    // TODO: Also, this does not compute FOF densities
#ifndef DIRECTPLUMMER
    QUIT("Error: executing CPU pencils with non-Plummer softening is not currently supported.\n");
#endif
    
    int WIDTH = 2*P.NearFieldRadius + 1;
    
    FLOAT cpu_eps = JJ->SofteningLengthInternal*JJ->SofteningLengthInternal;

    // Copy the sources and sinks into position
    FillSinks.Start();
    List3<FLOAT> *SinkSetPositions = new List3<FLOAT>(NSinkBlocks*NFBlockSize);
    #pragma omp parallel for schedule(dynamic,1)
    for (int j=0; j<NSinkSets; j++) {
        SinkPlan[j].copy_into_pinned_memory(*SinkSetPositions, SinkSetStart[j], SinkSetCount[j], SinkPosSlab);
    }
    FillSinks.Stop();

    FillSources.Start();
    List3<FLOAT> *SourceSetPositions = new List3<FLOAT>(NSourceBlocks*NFBlockSize);
    #pragma omp parallel for schedule(dynamic,1)
    for (int j=0; j<NSourceSets; j++) {
        SourcePlan[j].copy_into_pinned_memory(*SourceSetPositions, SourceSetStart[j], SourceSetCount[j], SourcePosSlab);
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
                SinkSetAccelerations[id].x = a.x;
                SinkSetAccelerations[id].y = a.y;
                SinkSetAccelerations[id].z = a.z;
            }
        }
    }
    SetCompleted();
    delete SinkSetPositions;
    delete SourceSetPositions;
}

void SetInteractionCollection::PrintInteractions(){

    if(Nk * nfwidth + k_mod < P.cpd)
        Nk++;
    int nfr = nfradius;

    printf("SIC for slab: %d w: %d, k: %d - %d\n",SlabId,k_mod,j_low,j_high);

    printf("\t%d Sink pencils and %d Source pencils\n\n",NSinkSets,NSourceSets);

    printf("\tSink Pencils:\n");

    // WARNING: This routine didn't swap the meaning of j & k
    for(int i = 0; i < NSinkSets; i++){
        int j = i%Nk;
        int k = i/Nk;
        int jj = k_mod + j * nfwidth;
        int zmid = PP->WrapSlab(jj + nfr);
        printf("\t\t%d: %d, %d, %d - %d\n", i, SlabId, k + j_low, zmid-nfr,zmid+nfr);
    }

    printf("\tSourcePencils:\n");
    for(int i = 0; i < NSourceSets; i++){
        int j = i%Nk;
        int jj = k_mod + j * nfwidth;
        int zmid = PP->WrapSlab(jj + nfr);
        int y = i/Nk + j_low - nfr;
        printf("\t\t%d: %d - %d, %d, %d\n", i, SlabId-nfr,SlabId+nfr, y, zmid);
    }


    printf("\n\tInteractionList (sinkIdx<-sourceIdx ||sink | source):\n");

    for(int i = 0; i < InteractionCount; i++){
        int sinkIdx = i/nfwidth;
        int sinkj = sinkIdx%Nk;
        int sinkk = sinkIdx/Nk;
        int sinkjj = k_mod + sinkj * nfwidth;
        int sinkzmid = PP->WrapSlab(sinkjj + nfr);

        int sourceIdx = SinkSourceInteractionList[i];
        int sourcej = sourceIdx%Nk;
        int sourcey = sourceIdx/Nk + j_low - nfr;
        int sourcejj =  k_mod + sourcej * nfwidth;
        int sourcezmid =  PP->WrapSlab(sourcejj + nfr);
        FLOAT yoffset = SinkSourceYOffset[i];
                
        printf("\t\t%d: %d<-)%d|| %d: %d, %d, %d - %d | %d: %d - %d, %d, %d, %f\n",
            i, sinkIdx, sourceIdx,
            sinkIdx,SlabId, sinkk + j_low, sinkzmid-nfr,sinkzmid+nfr,
            sourceIdx,SlabId-nfr,SlabId+nfr, sourcey, sourcezmid, yoffset);
    }    
}

void SetInteractionCollection::AddInteractionList( std::vector<uint64> ** il){
    // TODO: Not sure what this routine is doing!  Not sure how to adjust it.

    if(Nk*nfwidth + k_mod < P.cpd)
        Nk++;
    int nfr = nfradius;

    // WARNING: This routine didn't swap the meaning of j & k
    for(int i = 0; i < InteractionCount; i++){
        int sinkIdx = i/nfwidth;
        int sinkj = sinkIdx%Nk;
        int sinkk = sinkIdx/Nk;
        int sinkjj = k_mod + sinkj * nfwidth;
        int sinkzmid = PP->WrapSlab(sinkjj + nfr);

        int sourceIdx = SinkSourceInteractionList[i];
        int sourcej = sourceIdx%Nk;
        int sourcey = sourceIdx/Nk + j_low - nfr;
        int sourcejj =  k_mod + sourcej * nfwidth;
        int sourcezmid =  PP->WrapSlab(sourcejj + nfr);

        for(int sinkz = sinkzmid - nfr; sinkz <= sinkzmid+nfr; sinkz++){
            for(int sourcex = SlabId - nfr; sourcex <= SlabId + nfr; sourcex++){
                int sinky = PP->WrapSlab(sinkk + j_low);
                int sinkzw = PP->WrapSlab(sinkz);

                int sourcexw = PP->WrapSlab(sourcex);
                int sourceyw = PP->WrapSlab(sourcey);
                int sourcez = sourcezmid;

                uint64 sinkSlabCellIdx = cpd*sinky + sinkzw;
                uint64 sourceCellIdx = cpd*cpd*sourcexw + cpd*sourceyw + sourcez;

                il[sinkSlabCellIdx]->push_back(sourceCellIdx);

            }
        }
    }
}
