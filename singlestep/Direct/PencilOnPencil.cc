//code related to direct pencil interaction creation for directdriver

//Collection of Interactions of all 5 cell pencils in a specified region of a slab
//TODO:This is similar enough to the group multistep interaction list that they could share an interface/superclass

// TODO: There is heavy use of signed 32-bit integers in this package.
// For now, DJE has simply added a lot of asserts.

#include "StructureOfLists.cc"
#include "SetInteractionCollection.hh"



inline int SetInteractionCollection::NumPaddedBlocks(int nparticles) {
    // Given the number of particles, compute and return the number of blocks
    // required.  1..N means 1, N+1..2N means 2, etc.  0 means 0.
    return (nparticles+NFBlockSize-1)/NFBlockSize;
}

inline int SetInteractionCollection::PaddedSinkCount(int sinkindex) {
    // The GPU thread needs to figure out how much space has been allocated.
    // No point to allocate a big array when it's so easy to compute.
    return NFBlockSize*NumPaddedBlocks(SinkSetCount[sinkindex]);
}

inline int SetInteractionCollection::PaddedSourceCount(int sourceindex) {
    // The GPU thread needs to figure out how much space has been allocated
    return NFBlockSize*NumPaddedBlocks(SourceSetCount[sourceindex]);
}

SetInteractionCollection::SetInteractionCollection(int slab, int w, int k_low, int k_high){
    Construction.Start();

    eps = JJ->eps;
    
    //set known class variables
    CompletionFlag = 0;
    K_low = k_low;
    K_high = k_high;
    W = w;
    SlabId = slab;
    
    //usefull construction constants
    int width = 2*P.NearFieldRadius+1;
    int k_width = k_high-k_low;
    int Nj = (P.cpd - w)/width;
    if(Nj * width + w < P.cpd)
        Nj++;
    
    assertf( (uint64)k_width*Nj < INT32_MAX, 
            "The number of sink sets will overflow a 32-bit signed int");
    NSinkList = k_width * Nj;
    assertf(NSinkList <= MaxNSink, "NSinkList (%d) larger than allocated space (MaxNSink = %d)\n", NSinkList, MaxNSink);
    
    assert(posix_memalign((void **) &SinkSetStart, 4096, sizeof(int) * NSinkList) == 0);
    assert(posix_memalign((void **) &SinkSetCount, 4096, sizeof(int) * NSinkList) == 0);
    assert(posix_memalign((void **) &SinkPlan, 4096, sizeof(SinkPencilPlan) * NSinkList) == 0);

    // TODO: Need SinkPlan and SourcePlan in the class

    assert(posix_memalign((void **) &SinkSetIdMax, 4096, sizeof(int) * NSinkList) == 0);
    
    int localSinkTotal = 0;

    assertf( (uint64)Nj * (k_width + width) < INT32_MAX,
        "The number of source sets will overflow a 32-bit signed int");
    NSourceSets = Nj * (k_width + width);
    
    assertf(NSourceSets <= MaxNSource, "NSourceSets (%d) larger than allocated space (MaxNSource = %d)\n", NSourceSets, MaxNSource);
    assertf(NSourceSets <= P.cpd*(P.cpd+width), "NSourceSets (%d) exceeds SourceSet array allocation (%d)\n", NSourceSets, P.cpd*(P.cpd+width));
    
    assert(posix_memalign((void **) &SourceSetStart, 4096, sizeof(int) * NSourceSets) == 0);  // can this be of size NSourceSets?
    assert(posix_memalign((void **) &SourceSetCount, 4096, sizeof(int) * NSourceSets) == 0);
    assert(posix_memalign((void **) &SourcePlan, 4096, sizeof(SourcePencilPlan) * NSourceSets) == 0);
    
    int localSourceTotal = 0;

    DirectTotal = 0;

    //Fill in sink data
    FillSinkLists.Clear(); FillSinkLists.Start();
    //Count the Sinks
    CountSinks.Clear(); CountSinks.Start();

    uint64 skewer_blocks[k_width+width];   // Number of blocks in this k-skewer
            // This is oversized because we'll re-use for Sources

    #pragma omp parallel for schedule(static) reduction(+:localSinkTotal)
    for(int k = 0; k < k_width; k++){
        int this_skewer_blocks = 0;   // Just to have a local variable
        for(int j = 0; j < Nj; j ++) {
            int jj = w + j * width;
            int zmid = PP->WrapSlab(jj + P.NearFieldRadius);
            int sinkindex = k * Nj + j;
            // This loads all of the position pointers into the Plan,
            // but also returns the total size.
            int pencilsize = SinkPlan[sinkindex].load(slab, k + k_low, zmid);
            SinkSetCount[sinkindex] = pencilsize;
            localSinkTotal += pencilsize;
            this_skewer_blocks += NumPaddedBlocks(pencilsize);
        }
        skewer_blocks[k] = this_skewer_blocks;
    }
    SinkTotal = localSinkTotal;  // OpenMP can't do reductions directly on member variables
    CountSinks.Stop();
    
    CalcSinkBlocks.Clear(); CalcSinkBlocks.Start();

    // Cumulate the number of blocks in each skewer, so we know how 
    // to start enumerating.
    uint64 skewer_blocks_start[k_width+1+width]; 
            // This is oversized because we'll re-use for Sources
    skewer_blocks_start[0] = 0;
    for (int k=0; k<k_width; k++) 
            skewer_blocks_start[k+1] = skewer_blocks_start[k]+skewer_blocks[k];

    NSinkBlocks = skewer_blocks_start[k_width];   // The total for all skewers
    assertf(skewer_blocks_start[k_width]*NFBlockSize < INT32_MAX, 
            "The number of padded sink particles will overflow a 32-bit signed int");

    assertf(NSinkBlocks <= MaxSinkBlocks, "NSinkBlocks (%d) larger than allocated space (MaxSinkBlocks = %d)\n", NSinkBlocks, MaxSinkBlocks);
    assert(posix_memalign((void **) &SinkBlockParentPencil, 4096, sizeof(int) * NSinkBlocks) == 0);

    // Loop over the pencils to set these:
    //   SinkSetStart[set] points to where the padded particles starts
    //   SinkSetIdMax[set] points to where the padded particles end
    //   SinkBlockParentPencil[block] points back to the (k,j) sinkset enumeration
    //
    
    #pragma omp parallel for schedule(static) 
    for(int k = 0; k < k_width; k++){
        // Figure out where to start the enumeration in this skewer
        int block_start = skewer_blocks_start[k];
        int NPaddedSinks = NFBlockSize*block_start; 

        for(int j = 0; j < Nj; j ++) {
            int sinkset = k * Nj + j;
            int sinklength = SinkSetCount[sinkset];    // num particles
            SinkSetStart[sinkset] = NPaddedSinks;
            SinkSetIdMax[sinkset] = SinkSetStart[sinkset] + sinklength;
            int nblocks = NumPaddedBlocks(sinklength);
            NPaddedSinks += nblocks*NFBlockSize;

            int block_end = block_start+nblocks;
            for (int b=block_start; b<block_end; b++)
                SinkBlockParentPencil[b] = sinkset;
            block_start=block_end;
        }
    }

    int NPaddedSinks = NFBlockSize*NSinkBlocks;   
            // The total padded number of particles
    assertf(NPaddedSinks <= MaxSinkSize, "NPaddedSinks (%d) larger than allocated space (MaxSinkSize = %d)\n", NPaddedSinks, MaxSinkSize);

//    SinkSetPositions = new List3<FLOAT>(NPaddedSinks);

    assert(posix_memalign((void **) &SinkSetAccelerations, 4096, sizeof(FLOAT3) * NPaddedSinks) == 0);
    
    CalcSinkBlocks.Stop();
    
//    FillSinks.Start();
//    // Pencil filling should use static scheduling
//    // so threads aren't trying to write next to each other
//    #pragma omp parallel for schedule(static)
//    for(int sinkset = 0; sinkset < NSinkList; sinkset++) {
//        int j = sinkset%Nj;
//        int k = sinkset/Nj + k_low;
//        int jj = w + j * width;
//        int zmid = PP->WrapSlab(jj + P.NearFieldRadius);
//        CreateSinkPencil(slab, k, zmid, SinkSetStart[sinkset] );
//    }
//    FillSinks.Stop();

    FillSinkLists.Stop();

    // All done with the Sinks.  Now for the Sources.
    FillSourceLists.Start();
    CountSources.Start();

    // Notice that the Source Pencils go from [k_low-NFR, k_high+NFR).

    // We once again need to precompute the enumeration of the padded
    // blocks, totaled by skewer.  
    // DJE chose to reverse this loop ordering, to match the sink case.

    #pragma omp parallel for schedule(static) reduction(+:localSourceTotal)
    for(int k = 0; k<k_width+width; k++) {
        int this_skewer_blocks = 0;   // Just to have a local variable
        for(int j = 0; j < Nj; j ++) {
            int sourceindex = k * Nj + j;
            int jj = w + j * width;
            int zmid = PP->WrapSlab(jj + P.NearFieldRadius);
            // This loads the Plan and returns the length
            int sourcelength = SourcePlan[sourceindex].load(slab,
                            k+k_low-P.NearFieldRadius, zmid);
            SourceSetCount[sourceindex] = sourcelength;
            localSourceTotal += sourcelength;
            this_skewer_blocks += NumPaddedBlocks(sourcelength);
        }
        skewer_blocks[k] = this_skewer_blocks;
    }
    SourceTotal = localSourceTotal;  // OpenMP can't do reductions directly on member variables
    CountSources.Stop();

    CalcSourceBlocks.Clear(); CalcSourceBlocks.Start();

    // Cumulate the number of blocks in each skewer, so we know how 
    // to start enumerating.
    skewer_blocks_start[0] = 0;
    for (int k=0; k<k_width+width; k++) 
            skewer_blocks_start[k+1] = skewer_blocks_start[k]+skewer_blocks[k];

    NSourceBlocks = skewer_blocks_start[k_width+width];   
            // The total for all skewers
    assertf(skewer_blocks_start[k_width+width]*NFBlockSize < INT32_MAX, 
            "The number of padded source particles will overflow a 32-bit signed int");

    // SourceSetStart[set] holds the padded particle index
    #pragma omp parallel for schedule(static) 
    for(int k = 0; k<k_width+width; k++) {
        int block_start = skewer_blocks_start[k];
        int NPaddedSources = NFBlockSize*block_start; 
        for(int j = 0; j < Nj; j ++) {
            int sourceset = k * Nj + j;
            int sourcelength = SourceSetCount[sourceset];
            SourceSetStart[sourceset] = NPaddedSources;
            int nblocks = NumPaddedBlocks(sourcelength);
            NPaddedSources += nblocks*NFBlockSize;
        }
    }
    
    int NPaddedSources = NFBlockSize*NSourceBlocks;
            // The total number of padded sources

//    SourceSetPositions = new List3<FLOAT>(NPaddedSources);
    assertf(NPaddedSources <= MaxSourceSize, "NPaddedSources (%d) larger than allocated space (MaxSourceSize = %d)\n", NPaddedSources, MaxSourceSize);
    CalcSourceBlocks.Stop();

//    FillSources.Start();
//    #pragma omp parallel for schedule(static)
//    for(int sourceIdx = 0; sourceIdx < NSourceSets; sourceIdx ++){
//        int j = sourceIdx%Nj;
//        int jj = w + j * width;
//        int zmid = PP->WrapSlab(jj + P.NearFieldRadius);
//        int y = sourceIdx/Nj + k_low - P.NearFieldRadius;
//        CreateSourcePencil(slab, y, zmid, SourceSetStart[sourceIdx]);
//    }
//    FillSources.Stop();

    FillSourceLists.Stop();


    //fill the interaction lists for the sink sets
    FillInteractionList.Start();

    InteractionCount = width*k_width*Nj;
    assertf(InteractionCount <= MaxNSink*width, "InteractionCount (%d) larger than allocated space (MaxNSink * width = %d)\n", InteractionCount, MaxNSink * width);
    assertf((uint64)width*k_width*Nj < INT32_MAX, 
            "Interaction Count exceeds 32-bit signed int");
    assert(posix_memalign((void **) &SinkSourceInteractionList, 4096, sizeof(int) * InteractionCount) == 0);
    assert(posix_memalign((void **) &SinkSourceYOffset, 4096, sizeof(FLOAT) * InteractionCount) == 0);
    FLOAT cellsize = PP->invcpd;
    
    uint64 localDirectTotal = 0;
    #pragma omp parallel for schedule(static) reduction(+:localDirectTotal)
    for(int k = 0; k < k_width; k++){
        int g = omp_get_thread_num();
        assertf(k*Nj + Nj <= NSinkList, "SinkSetCount array access at %d would exceed allocation %d\n", k*Nj + Nj, NSinkList);
        for(int j=0; j < Nj; j++) {
            int jj = w + j * width;
            int zmid = PP->WrapSlab(jj + P.NearFieldRadius);

            int sinkindex = k*Nj + j;
            int l = width * sinkindex;

            assertf(l + width <= InteractionCount, "SinkSourceInteractionList array access at %d would exceed allocation %d\n", l + width, InteractionCount);
            assertf((k + width)*(Nj) +  j <= P.cpd*(P.cpd+width), "SourceSetCount array access at %d would exceed allocation %d\n", (k + width)*(Nj), P.cpd*(P.cpd+width));
            for(int y=0;y<width;y++) {
                int sourceindex = (k + y)*(Nj) +  j;
                SinkSourceInteractionList[l+y] = sourceindex;
                #ifndef GLOBAL_POS
                    SinkSourceYOffset[l+y] = (y-width/2)*cellsize;
                #else
                    // The y coordinate is (k_low+k+y-width/2)
                    int tmpy = k_low+k+y-width/2;
                    SinkSourceYOffset[l+y] = (tmpy-PP->WrapSlab(tmpy))*cellsize;
                #endif
                localDirectTotal += SinkSetCount[sinkindex] * SourceSetCount[sourceindex];
            }
        }
    }
    DirectTotal = localDirectTotal;

    FillInteractionList.Stop();
    Construction.Stop();
}

SetInteractionCollection::~SetInteractionCollection(){
    free(SinkSetStart);
    free(SinkSetCount);
    free(SinkPlan);
    free(SinkSetIdMax);
    free(SourceSetStart);
    free(SourceSetCount);
    free(SourcePlan);
    free(SinkBlockParentPencil);
    free(SinkSetAccelerations);
    free(SinkSourceInteractionList);
    free(SinkSourceYOffset);
}

void SetInteractionCollection::SetCompleted(){
    STDLOG(1,"Completed SIC for slab %d w: %d k: %d - %d\n",SlabId,W,K_low,K_high); 

    // delete SinkSetPositions;
    // delete SourceSetPositions;

    CompletionFlag = 1;
}

//inline int SetInteractionCollection::SinkPencilCount(int sinkx, int sinky, int sinkz) {
//    int s = 0;
//    int zmax = sinkz + P.NearFieldRadius;
//    for(int z = sinkz - P.NearFieldRadius; z <= zmax; z++) {
//        s += PP->NumberParticle(sinkx,sinky,z);
//    }
//    return s;
//}


//inline int SetInteractionCollection::SourcePencilCount(int slab, int ymid, int zz) { 
//    int s = 0;
//    int xmax = slab + P.NearFieldRadius;
//    for(int x = slab - P.NearFieldRadius; x <= xmax; x++) {
//        s += PP->NumberParticle(x,ymid,zz);
//    }
//    return s;
//}


/*inline void SetInteractionCollection::CreateSinkPencil(int sinkx, int sinky, int sinkz, uint64 start){
    int zmax = sinkz+P.NearFieldRadius;
    for(int z=sinkz-P.NearFieldRadius;z<=zmax;z++) {
        int count = PP->NumberParticle(sinkx,sinky,z);

        if(count>0) {
            List3<FLOAT> pos = PP->PosXYZCell(sinkx,sinky,z);
            
            FLOAT3 newsinkcenter = PP->CellCenter(sinkx, sinky, z);
            FLOAT * X = &(SinkSetPositions->X[start]);
            FLOAT * Y = &(SinkSetPositions->Y[start]);
            FLOAT * Z = &(SinkSetPositions->Z[start]);
            
            #pragma simd assert
            for(int p = 0; p < count; p++)
                X[p] = pos.X[p] + newsinkcenter.x;
            #pragma simd assert
            for(int p = 0; p < count; p++)
                Y[p] = pos.Y[p] + newsinkcenter.y;
            #pragma simd assert
            for(int p = 0; p < count; p++)
                Z[p] = pos.Z[p] + newsinkcenter.z;
            
            start+=count;
        }
    }
    assertf(start <= SinkSetPositions->N, "Wrote to particle %d. This overruns SinkSetPositions (length %d)\n", start, SinkSetPositions->N);
}

inline void SetInteractionCollection::CreateSourcePencil(int sx, int sy, int nz, uint64 start){
    int xmax = sx+P.NearFieldRadius;
    for(int x=sx-P.NearFieldRadius;x<=xmax;x++)  {
        int count = PP->NumberParticle(x,sy,nz);
        if(count>0) {
            List3<FLOAT> pos = PP->PosXYZCell(x,sy,nz);
            
            FLOAT3 newsourcecenter = PP->CellCenter(x,sy,nz);
            FLOAT * X = &(SourceSetPositions->X[start]);
            FLOAT * Y = &(SourceSetPositions->Y[start]);
            FLOAT * Z = &(SourceSetPositions->Z[start]);
            
            #pragma simd assert
            for(int p = 0; p < count; p++)
                X[p] = pos.X[p] + newsourcecenter.x;
            #pragma simd assert
            for(int p = 0; p < count; p++)
                Y[p] = pos.Y[p] + newsourcecenter.y;
            #pragma simd assert
            for(int p = 0; p < count; p++)
                Z[p] = pos.Z[p] + newsourcecenter.z;
            
            start += count;
        }
    }
    assertf(start <= SourceSetPositions->N, "Wrote to particle %d. This overruns SourceSetPositions (length %d)\n", start, SourceSetPositions->N);
}*/

int SetInteractionCollection::CheckCompletion(){
    return CompletionFlag;
}

// Currently Plummer only
void SetInteractionCollection::CPUExecute(){
    int WIDTH = 2*P.NearFieldRadius + 1;
    
    FLOAT cpu_eps = JJ->SofteningLengthInternal*JJ->SofteningLengthInternal;

    // Copy the sources and sinks into position
    FillSinks.Start();
    List3<FLOAT> *SinkSetPositions = new List3<FLOAT>(NSinkBlocks*NFBlockSize);
    #pragma omp parallel for schedule(dynamic,1)
    for (int j=0; j<NSinkList; j++) {
        SinkPlan[j].copy_into_pinned_memory(*SinkSetPositions, SinkSetStart[j], SinkSetCount[j]);
    }
    FillSinks.Stop();

    FillSources.Start();
    List3<FLOAT> *SourceSetPositions = new List3<FLOAT>(NSourceBlocks*NFBlockSize);
    #pragma omp parallel for schedule(dynamic,1)
    for (int j=0; j<NSourceSets; j++) {
        SourcePlan[j].copy_into_pinned_memory(*SourceSetPositions, SourceSetStart[j], SourceSetCount[j]);
    }
    FillSources.Stop();

    for(int blockIdx = 0; blockIdx < NSinkBlocks; blockIdx++){
        #pragma omp parallel for schedule(dynamic,1)
        for(int threadIdx = 0; threadIdx < NFBlockSize; threadIdx++){

            int id = NFBlockSize*blockIdx + threadIdx;
            int sinkIdx = SinkBlockParentPencil[blockIdx];

            FLOAT sinkX, sinkY, sinkZ;
            if(id <  SinkSetStart[sinkIdx] + SinkSetCount[sinkIdx]){
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

                        r = RSQRT( drx*drx + dry*dry + drz*drz  + cpu_eps);
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

                        r = RSQRT( drx*drx + dry*dry + drz*drz  + cpu_eps);
                        r *=r*r;
                        a.x -= r * drx;
                        a.y -= r * dry;
                        a.z -= r * drz;
                    }
                //}

            }

            if(id < SinkSetStart[sinkIdx] + SinkSetCount[sinkIdx]){
                assert(isfinite(a.x));
                assert(isfinite(a.y));
                assert(isfinite(a.z));
                SinkSetAccelerations[id] = a;
            }
        }
    }
    SetCompleted();
    delete SinkSetPositions;
    delete SourceSetPositions;
}

void SetInteractionCollection::PrintInteractions(){

    int width = 2*P.NearFieldRadius+1;
    int k_width = K_high-K_low;
    int Nj = (P.cpd - W)/width;
    if(Nj * width + W < P.cpd)
        Nj++;
    int nfr = P.NearFieldRadius;

    printf("SIC for slab: %d w: %d, k: %d - %d\n",SlabId,W,K_low,K_high);

    printf("\t%d Sink pencils and %d Source pencils\n\n",NSinkList,NSourceSets);

    printf("\tSink Pencils:\n");

    for(int i = 0; i < NSinkList; i++){
        int j = i%Nj;
        int k = i/Nj;
        int jj = W + j * width;
        int zmid = PP->WrapSlab(jj + nfr);
        printf("\t\t%d: %d, %d, %d - %d\n", i, SlabId, k + K_low, zmid-nfr,zmid+nfr);
    }

    printf("\tSourcePencils:\n");
    for(int i = 0; i < NSourceSets; i++){
        int j = i%Nj;
        int jj = W + j * width;
        int zmid = PP->WrapSlab(jj + nfr);
        int y = i/Nj + K_low - nfr;
        printf("\t\t%d: %d - %d, %d, %d\n", i, SlabId-nfr,SlabId+nfr, y, zmid);
    }


    printf("\n\tInteractionList (sinkIdx<-sourceIdx ||sink | source):\n");

    for(int i = 0; i < InteractionCount; i++){
        int sinkIdx = i/width;
        int sinkj = sinkIdx%Nj;
        int sinkk = sinkIdx/Nj;
        int sinkjj = W + sinkj * width;
        int sinkzmid = PP->WrapSlab(sinkjj + nfr);

        int sourceIdx = SinkSourceInteractionList[i];
        int sourcej = sourceIdx%Nj;
        int sourcey = sourceIdx/Nj + K_low - nfr;
        int sourcejj =  W + sourcej * width;
        int sourcezmid =  PP->WrapSlab(sourcejj + nfr);
        FLOAT yoffset = SinkSourceYOffset[i];
                
        printf("\t\t%d: %d<-)%d|| %d: %d, %d, %d - %d | %d: %d - %d, %d, %d, %f\n",
            i, sinkIdx, sourceIdx,
            sinkIdx,SlabId, sinkk + K_low, sinkzmid-nfr,sinkzmid+nfr,
            sourceIdx,SlabId-nfr,SlabId+nfr, sourcey, sourcezmid, yoffset);
    }    
}

void SetInteractionCollection::AddInteractionList( std::vector<uint64> ** il){
    // TODO: Not sure what this routine is doing!  Not sure how to adjust it.

    int width = 2*P.NearFieldRadius+1;
    int k_width = K_high-K_low;
    int Nj = (P.cpd - W)/width;
    if(Nj*width + W < P.cpd)
        Nj++;
    int nfr = P.NearFieldRadius;
    int cpd = P.cpd;

    for(int i = 0; i < InteractionCount; i++){
        int sinkIdx = i/width;
        int sinkj = sinkIdx%Nj;
        int sinkk = sinkIdx/Nj;
        int sinkjj = W + sinkj * width;
        int sinkzmid = PP->WrapSlab(sinkjj + nfr);

        int sourceIdx = SinkSourceInteractionList[i];
        int sourcej = sourceIdx%Nj;
        int sourcey = sourceIdx/Nj + K_low - nfr;
        int sourcejj =  W + sourcej * width;
        int sourcezmid =  PP->WrapSlab(sourcejj + nfr);

        for(int sinkz = sinkzmid - nfr; sinkz <= sinkzmid+nfr; sinkz++){
            for(int sourcex = SlabId - nfr; sourcex <= SlabId + nfr; sourcex++){
                int sinky = PP->WrapSlab(sinkk + K_low);
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

#ifndef CUDADIRECT
// If we're not compiling the GPU code, need a stub for this function
void SetInteractionCollection::GPUExecute(int){
    assertf(0, "Abacus was not compiled with CUDA.  Try running ./configure again?\n");
}
#endif

volatile int SetInteractionCollection::ActiveThreads = 0;
pthread_mutex_t SetInteractionCollection::GPUTimerMutex = PTHREAD_MUTEX_INITIALIZER;
STimer SetInteractionCollection::GPUThroughputTimer;

