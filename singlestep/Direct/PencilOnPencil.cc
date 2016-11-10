//code related to direct pencil interaction creation for directdriver

//Collection of Interactions of all 5 cell pencils in a specified region of a slab
//TODO:This is similar enough to the group multistep interaction list that they could share an interface/superclass
#include "StructureOfLists.cc"
#include "SetInteractionCollection.hh"

SetInteractionCollection::SetInteractionCollection(int slab, int w, int k_low, int k_high){
    Construction.Start();
    CopyTime = 0.;
    ExecutionTime = 0.;
    CopybackTime = 0.;
    TotalTime = 0.;

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
    
    NSinkList = k_width * Nj;
    SinkSetStart = (int *) malloc(sizeof(int) * NSinkList);
    SinkSetCount = (int *) malloc(sizeof(int) * NSinkList);
    SinkTotal = 0;

    NSourceSets = Nj * (k_width + width);
    assertf(NSourceSets <= P.cpd*(P.cpd+width), "NSourceSets (%d) exceeds SourceSet array allocation (%d)\n", NSourceSets, P.cpd*(P.cpd+width));
    SourceSetStart = (int *) malloc(sizeof(int) * P.cpd*(P.cpd+width));
    SourceSetCount = (int *) malloc(sizeof(int) * P.cpd*(P.cpd+width));

    DirectTotal = 0;

    //Fill in sink data
    FillSinkLists.Clear(); FillSinkLists.Start();
    //Count the Sinks
    CountSinks.Clear(); CountSinks.Start();
    #pragma omp parallel for schedule(dynamic,1)
    for(int k = 0; k < k_width; k++){
        for(int j = 0; j < Nj; j ++) {
            int jj = w + j * width;
            int zmid = PP->WrapSlab(jj + P.NearFieldRadius);
            int sinkindex = k * Nj + j;
            int pencilsize = SinkPencilCount(slab, k + k_low, zmid );
            SinkSetCount[sinkindex] = pencilsize;
            SinkTotal += pencilsize;
        }
    }
    CountSinks.Stop();
    
    CalcSinkBlocks.Clear(); CalcSinkBlocks.Start();
    //Count the total number of sink blocks we will require
    //TODO: Another non-trivialy parallelizable sum
    int NPaddedSinks = 0;
    NSinkBlocks = 0;
    for(int sinkset=0; sinkset < NSinkList; sinkset++) {
        int sinklength = SinkSetCount[sinkset];
        SinkSetStart[sinkset] = NPaddedSinks; 

        int PencilBlocks = sinklength/NFBlockSize;
        int remaining = (NFBlockSize - (sinklength%NFBlockSize))%NFBlockSize;
        PencilBlocks+= remaining&&remaining;
        NSinkBlocks += PencilBlocks;
        NPaddedSinks+= NFBlockSize*PencilBlocks;
    }
    assert(NPaddedSinks==NFBlockSize*NSinkBlocks);

    SinkBlockParentPencil = (int *) malloc(sizeof(int) * NSinkBlocks);
    for(int sinkset=0; sinkset < NSinkList; sinkset++) {
        if(SinkSetCount[sinkset] == 0) continue;
        int PencilBlocks = SinkSetCount[sinkset]/NFBlockSize;
        int remaining = (NFBlockSize - (SinkSetCount[sinkset]%NFBlockSize))%NFBlockSize;
        PencilBlocks+= remaining&&remaining;
        int b0 =  SinkSetStart[sinkset]/NFBlockSize;
        assertf(b0 + PencilBlocks <= NSinkBlocks, "SinkBlockParentPencil array access at index %d will exceed allocation %d.\n",b0 + PencilBlocks, NSinkBlocks);
        for(int b = b0; b < b0 + PencilBlocks; b++)
            SinkBlockParentPencil[b] = sinkset;
    }


    SinkSetPositions = new List3<FLOAT>(NPaddedSinks);

    SinkSetAccelerations = (FLOAT3 *) malloc(sizeof(FLOAT3) * NPaddedSinks);
    //to make valgrind more accurate, we only zero the regions we use
    //memset(SinkSetAccelerations, 0 , sizeof(FLOAT3)*NPaddedSinks);
    CalcSinkBlocks.Stop();

    FillSinks.Start();

    #pragma omp parallel for schedule(dynamic,1) 
    for(int sinkset = 0; sinkset < NSinkList; sinkset++) {
        int j = sinkset%Nj;
        int k = sinkset/Nj + k_low;
        int jj = w + j * width;
        int zmid = PP->WrapSlab(jj + P.NearFieldRadius);
        CreateSinkPencil(slab, k, zmid, SinkSetStart[sinkset] );
        memset(&(SinkSetAccelerations[SinkSetStart[sinkset]]), 0 , sizeof(FLOAT3)*SinkSetCount[sinkset]);
    }

    FillSinks.Stop();
    FillSinkLists.Stop();

    //Fill in source data
    FillSourceLists.Start();
    CountSources.Start();
    memset(SourceSetCount, 0, sizeof(int) * NSourceSets );
    #pragma omp parallel for schedule(dynamic,1)
    for(int j = 0; j < Nj; j ++) {
        int jj = w + j * width;
        int zmid = PP->WrapSlab(jj + P.NearFieldRadius);
        for(int y = k_low; y < k_high + width; y++) {
            int sourceindex = (y - k_low) * Nj + j;
            int sourcelength = SourcePencilCount( slab, y - P.NearFieldRadius, zmid);
            SourceSetCount[sourceindex] = sourcelength;
        }
    }
    CountSources.Stop();

    CalcSourceBlocks.Clear(); CalcSourceBlocks.Start();
    //Count the total number of sink blocks we will require
    //TODO: Another non-trivialy parallelizable sum
    int NPaddedSources = 0;
    NSourceBlocks = 0;
    for(int sourceset=0; sourceset < NSourceSets; sourceset++) {
        int sourcelength = SourceSetCount[sourceset];
        SourceSetStart[sourceset] = NPaddedSources;

        NPaddedSources += sourcelength;
        NSourceBlocks  += sourcelength/NFBlockSize;

        int remaining = (NFBlockSize - (sourcelength%NFBlockSize))%NFBlockSize;
        NPaddedSources += remaining;
        NSourceBlocks += remaining&&remaining;
    }
    assert(NPaddedSources==NFBlockSize*NSourceBlocks);

    SourceSetPositions = new List3<FLOAT>(NPaddedSources);

    CalcSourceBlocks.Stop();

    FillSources.Start();
    #pragma omp parallel for schedule(dynamic,1)
    for(int sourceIdx = 0; sourceIdx < NSourceSets; sourceIdx ++){
        int j = sourceIdx%Nj;
        int jj = w + j * width;
        int zmid = PP->WrapSlab(jj + P.NearFieldRadius);
        int y = sourceIdx/Nj + k_low - P.NearFieldRadius;
        CreateSourcePencil(slab, y, zmid, SourceSetStart[sourceIdx]);
    }
    FillSources.Stop();

    FillSourceLists.Stop();


    //fill the interaction lists for the sink sets
    FillInteractionList.Start();

    InteractionCount = width*k_width*Nj;
    SinkSourceInteractionList = (int *) malloc(sizeof(int) * InteractionCount);
    
    int cpd = P.cpd;
    int nprocs = omp_get_max_threads();
    uint64 threadDI[nprocs];
    memset(threadDI,0,sizeof(uint64)*nprocs);
    #pragma omp parallel for schedule(dynamic,1)
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
                threadDI[g] += SinkSetCount[sinkindex] * SourceSetCount[sourceindex];
            }
        }
    }
    DirectTotal = 0;
    for(int g = 0; g< nprocs; g++) DirectTotal+= threadDI[g];

    FillInteractionList.Stop();
    Construction.Stop();
}

SetInteractionCollection::~SetInteractionCollection(){
    free(SinkSetStart);
    free(SinkSetCount);
    free(SourceSetStart);
    free(SourceSetCount);
    free(SinkBlockParentPencil);
    free(SinkSetAccelerations);
    free(SinkSourceInteractionList);
}

void SetInteractionCollection::SetCompleted(){
    STDLOG(1,"Completed SIC for slab %d w: %d k: %d - %d\n",SlabId,W,K_low,K_high); 
    //free(SourceSetStart);
    //free(SourceSetCount);
    //free(SinkSourceInteractionList);

    delete SinkSetPositions;
    delete SourceSetPositions;

    CompletionFlag = 1;
}

inline int SetInteractionCollection::SinkPencilCount(int sinkx, int sinky, int sinkz) {
    int s = 0;
    for(int z=sinkz-P.NearFieldRadius;z<=sinkz+P.NearFieldRadius;z++) {
        s += PP->NumberParticle(sinkx,sinky,z);
    }
    return s;
}


inline int SetInteractionCollection::SourcePencilCount(int slab, int ymid, int zz) { 
    int s = 0;
    int xmax = slab+P.NearFieldRadius;
    for(int x=slab-P.NearFieldRadius;x<=xmax;x++)  {
        s += PP->NumberParticle(x,ymid,zz);
    }
    return s;
}

inline void SetInteractionCollection::CreateSinkPencil(int sinkx, int sinky, int sinkz, uint64 start){
    int zmax = sinkz+P.NearFieldRadius;
    for(int z=sinkz-P.NearFieldRadius;z<=zmax;z++) {
        posstruct * pos = PP->PosCell(sinkx,sinky,z);
        int count = PP->NumberParticle(sinkx,sinky,z);

        if(count>0) {
            //memcpy(&(sinkpositions[offset],&pos,sizeof(FLOAT3)*count));
            FLOAT3 newsinkcenter = PP->CellCenter(sinkx, sinky, z);
            FLOAT * X = &(SinkSetPositions->X[start]);
            FLOAT * Y = &(SinkSetPositions->Y[start]);
            FLOAT * Z = &(SinkSetPositions->Z[start]);

            for(int p = 0; p < count; p++){
                X[p] = pos[p].x; 
                Y[p] = pos[p].y;
                Z[p] = pos[p].z;
            }

            #pragma simd
            for(int p=0;p<count;p++){
                X[p] += newsinkcenter.x;
                Y[p] += newsinkcenter.y;
                Z[p] += newsinkcenter.z;
            }
            start+=count;
        }
    }
    assertf(start <= SinkSetPositions->N, "Wrote %d particles. This overruns SinkSetPositions (length %d)\n", start, SinkSetPositions->N);
}

inline void SetInteractionCollection::CreateSourcePencil(int sx, int sy, int nz, uint64 start){
    int offset = 0;
    int xmax = sx+P.NearFieldRadius;
    for(int x=sx-P.NearFieldRadius;x<=xmax;x++)  {
        posstruct * pos = PP->PosCell(x,sy,nz);
        int count = PP->NumberParticle(x,sy,nz);
        if(count>0) {
            //memcpy(&(sourcepositions[offset],&pos,sizeof(FLOAT3)*count));
            FLOAT3 newsourcecenter = PP->CellCenter(x,sy,nz);
            FLOAT * X = &(SourceSetPositions->X[start]);
            FLOAT * Y = &(SourceSetPositions->Y[start]);
            FLOAT * Z = &(SourceSetPositions->Z[start]);

            for(int p = 0; p < count; p++){
                X[p+offset] = pos[p].x;
                Y[p+offset] = pos[p].y;
                Z[p+offset] = pos[p].z;
            }

            #pragma simd
            for(int p=offset;p<count + offset;p++){
                X[p] += newsourcecenter.x;
                Y[p] += newsourcecenter.y;
                Z[p] += newsourcecenter.z;
            }
            offset+=count;
        }
    }
    assertf(start + offset <= SourceSetPositions->N, "Wrote %d particles starting at %d. This overruns SourceSetPositions (length %d)\n", offset, start, SourceSetPositions->N);
}

int SetInteractionCollection::CheckCompletion(){
    return CompletionFlag;
}

// Currently Plummer only
void SetInteractionCollection::CPUExecute(){
    int WIDTH = 2*P.NearFieldRadius + 1;
    
    FLOAT cpu_eps = JJ->SofteningLengthInternal*JJ->SofteningLengthInternal;

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
                        dry = sourceY - sinkY;
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
                        dry = sourceY - sinkY;
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
                
        printf("\t\t%d: %d<-)%d|| %d: %d, %d, %d - %d | %d: %d - %d, %d, %d\n",
            i, sinkIdx, sourceIdx,
            sinkIdx,SlabId, sinkk + K_low, sinkzmid-nfr,sinkzmid+nfr,
            sourceIdx,SlabId-nfr,SlabId+nfr, sourcey, sourcezmid);
    }    
}

void SetInteractionCollection::AddInteractionList( std::vector<uint64> ** il){

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

