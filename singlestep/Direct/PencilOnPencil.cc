//code related to direct pencil interaction creation for directdriver

//Collection of Interactions of all 5 cell pencils in a specified region of a slab
//TODO:This is simmilar enough to the group multistep interaction list that they could share an interface/superclass
#include "StructureOfLists.cc"
#include "SetInteractionCollection.hh"

SetInteractionCollection::SetInteractionCollection(int slab, int w, int k_low, int k_high){
    eps2 = P.SofteningLength;
    int k_width = k_high-k_low;
    CompletionFlag = 0;
    //Allocate regions of known size
    K_low = k_low;
    K_high = k_high;
    SinkSetStart = (int *) malloc(sizeof(int) * P.cpd*k_width);
    SinkSetCount = (int *) malloc(sizeof(int) * P.cpd*k_width);
    memset(SinkSetCount,0,sizeof(int) * P.cpd*P.cpd);

    SinkSetInteractionListStart = (int *) malloc(sizeof(int) * P.cpd*k_width);

    int width = 2*P.NearFieldRadius+1;
    
    SourceSetStart = (int *) malloc(sizeof(int) * P.cpd*(P.cpd+width));
    SourceSetCount = (int *) malloc(sizeof(int) * P.cpd*(P.cpd+width));

    DirectTotal = 0;

    //Fill in sink data
    FillSinkLists.Clear(); FillSinkLists.Start();
    //Count the Sinks
    CountSinks.Clear(); CountSinks.Start();
    #pragma omp parallel for schedule(dynamic,1)
    for(int k = K_low; k < K_high; k++){
        for(int j = w; j < P.cpd; j += width) {
            int zmid = PP->WrapSlab(j + P.NearFieldRadius);
            int sinkindex = k * P.cpd + zmid;
            int pencilsize = SinkPencilCount(slab, k, zmid );
            SinkSetCount[sinkindex] = pencilsize;
        }
    }
    CountSinks.Stop();
    
    //Calculate the starting index for each sink set's positions 
    //TODO:This is a classic prefix sum and can be parallelized if it is a time sink
    CalcSinkOffset.Clear(); CalcSinkOffset.Start();
    int nsinksets = 0;
    int CountSinkPencilPositions = 0;
    for(int sinkindex = 0; sinkindex < P.cpd*k_width; sinkindex++) {
        int sinklength = SinkSetCount[sinkindex];          
        nsinksets++;
        SinkSetStart[sinkindex] = CountSinkPencilPositions;
        CountSinkPencilPositions += SinkSetCount[sinkindex];
    }
    NSinkList = nsinksets;
    CalcSinkOffset.Stop();

    CalcSinkBlocks.Clear(); CalcSinkBlocks.Start();
    //Count the total number of sink blocks we will require
    //TODO: Another non-trivialy parallelizable sum
    int NPaddedSinks = 0;
    NSinkBlocks = 0;
    for(int sinkset=0; sinkset < NSinkList; sinkset++) {
        int sinkindex = sinkset;
        int sinklength = SinkSetCount[sinkindex];
        SinkSetStart[sinkindex] = NPaddedSinks; 

        NPaddedSinks += sinklength;
        int PencilBlocks = sinklength/NFBlockSize;

        int remaining = (NFBlockSize - (sinklength%NFBlockSize))%NFBlockSize;
        NPaddedSinks += remaining;
        //SinkSetCount[sinkindex] += remaining;
        PencilBlocks+= remaining&&remaining;
        NSinkBlocks += PencilBlocks;
    }
    assert(NPaddedSinks==NFBlockSize*NSinkBlocks);

    SinkBlockParentPencil = (int *) malloc(sizeof(int) * NSinkBlocks);
    for(int sinkset=0; sinkset < NSinkList; sinkset++) {
        int PencilBlocks = SinkSetCount[sinkset]/NFBlockSize;
        int remaining = (NFBlockSize - (SinkSetCount[sinkset]%NFBlockSize))%NFBlockSize;
        PencilBlocks+= remaining&&remaining;
        for(int b = SinkSetStart[sinkset]/NFBlockSize; b < PencilBlocks; b++)
            SinkBlockParentPencil[b] = sinkset;
    }


    SinkSetPositions = new List3<FLOAT>(NPaddedSinks);
    SinkSetPositions->SetZero(0,NPaddedSinks);

    SinkSetAccelerations = (FLOAT3 *) malloc(sizeof(FLOAT3) * CountSinkPencilPositions);
    memset(SinkSetAccelerations, 0 , sizeof(FLOAT3)*CountSinkPencilPositions);
    CalcSinkBlocks.Stop();

    FillSinks.Start();

    #pragma omp parallel for schedule(dynamic,1) 
    for(int sinkset = 0; sinkset < nsinksets; sinkset++) {
        int sinkindex = sinkset;
        int zmid = sinkindex%P.cpd;
        int k = (sinkindex - zmid)/P.cpd;
        CreateSinkPencil(slab, k, zmid, SinkSetStart[sinkindex] );
    }

    FillSinks.Stop();
    FillSinkLists.Stop();

    //Fill in source data
    FillSourceLists.Start();
    CountSources.Start();
    NSourceSets = P.cpd * (k_width + width);
    memset(SourceSetCount, 0, sizeof(int) * NSourceSets );
    #pragma omp parallel for schedule(dynamic,1)
    for(int j = w; j < P.cpd; j += width) {
        int zmid = PP->WrapSlab(j + P.NearFieldRadius);
        for(int y = K_low - P.NearFieldRadius; y <= K_high + P.NearFieldRadius; y++) {
            int sourceindex = (y + P.NearFieldRadius) * (P.cpd) + (zmid);
            int sourcelength = SourcePencilCount( slab, y, zmid);
            SourceSetCount[sourceindex] = sourcelength;
        }
    }
    CountSources.Stop();

    CalcSourceOffset.Start();
    //TODO: As with the sinks, this is a prefix sum which can be done efficiently in parallel
    int CountSourcePencilPositions = 0;
    for(int y = K_low - P.NearFieldRadius; y <= K_high + P.NearFieldRadius; y++) {
        for(int j = w; j < P.cpd; j += width){
            int zmid = PP->WrapSlab(j + P.NearFieldRadius);
            int sourceindex = (y + P.NearFieldRadius)*(P.cpd) + (zmid);
            SourceSetStart[sourceindex] = CountSourcePencilPositions;
            CountSourcePencilPositions += SourceSetCount[sourceindex];
        }
    }
    CalcSourceOffset.Stop();

    CalcSourceBlocks.Clear(); CalcSourceBlocks.Start();
    //Count the total number of sink blocks we will require
    //TODO: Another non-trivialy parallelizable sum
    int NPaddedSources = 0;
    NSourceBlocks = 0;
    for(int sourceset=0; sourceset < NSourceSets; sourceset++) {
        int sourcelength = SourceSetCount[sourceset];
        SinkSetStart[sourceset] = NPaddedSources;

        NPaddedSources += sourcelength;
        NSourceBlocks  += sourcelength/NFBlockSize;

        int remaining = (NFBlockSize - (sourcelength%NFBlockSize))%NFBlockSize;
        NPaddedSources += remaining;
        //SourceSetCount[sourceindex] += remaining;
        NSourceBlocks += remaining&&remaining;
    }
    assert(NPaddedSources==NFBlockSize*NSourceBlocks);

    SourceSetPositions = new List3<FLOAT>(NPaddedSources);
    SourceSetPositions->SetZero(0,NPaddedSources);

    CalcSourceBlocks.Stop();

    FillSources.Start();
    #pragma omp parallel for schedule(dynamic,1)
    for(int j = w; j < P.cpd; j+= 2 * P.NearFieldRadius + 1) {
        for(int y= K_low - P.NearFieldRadius; y <= K_high + P.NearFieldRadius; y++) {
            int zmid = PP->WrapSlab(j + P.NearFieldRadius);
            int sourceindex = (y + P.NearFieldRadius) * (P.cpd) + zmid;
            CreateSourcePencil(slab, y, zmid, SourceSetStart[sourceindex]);
        }
    }
    FillSources.Stop();

    FillSourceLists.Stop();


    //fill the interaction lists for the sink sets
    FillInteractionList.Start();

    #pragma omp parallel for schedule(dynamic,1)
    for(int k = K_low; k < K_high; k++){
        for(int j = w; j < P.cpd; j += width) {
            int zmid = PP->WrapSlab(j + P.NearFieldRadius);
            int sinkindex = k*P.cpd + zmid;
            SinkSetInteractionListStart[sinkindex] = width*((k-k_low)*((P.cpd-w)/width+1) +((j-w)/width))  ;
        }
    }

    int InteractionCount = width*((k_high-k_low)*((P.cpd-w)/width+1));
    SinkSourceInteractionList = (int *) malloc(sizeof(int) * InteractionCount);
    
    int nfr = P.NearFieldRadius;
    int cpd = P.cpd;
    int nprocs = omp_get_max_threads();
    uint64 threadDI[nprocs];
    memset(threadDI,0,sizeof(uint64)*nprocs);
    #pragma omp parallel for schedule(dynamic,1)
    for(int k = K_low; k <= K_high; k++){
        int g = omp_get_thread_num();
        for(int j=w;j<cpd;j+=2*nfr+1) {
            int zmid = PP->WrapSlab(j+nfr);
            int sinkindex = k*P.cpd + zmid;
            int l = SinkSetInteractionListStart[ sinkindex ];
            
            for(int y=0;y<width;y++) {
                int yy = y + k-nfr;
                int sourceindex = (yy+nfr)*(cpd) +  zmid;
                SinkSourceInteractionList[l+y] = sourceindex;
                threadDI[g] += SinkSetCount[sinkindex] * SourceSetCount[sourceindex];
            }
        }
    }
    DirectTotal = 0;
    for(int g = 0; g< nprocs; g++) DirectTotal+= threadDI[g];

    FillInteractionList.Stop();
}

SetInteractionCollection::~SetInteractionCollection(){
    free(SinkSetStart);
    free(SinkSetCount);

    free(SinkSetAccelerations);
}

void SetInteractionCollection::SetCompleted(){
    free(SinkSetInteractionListStart);
    free(SourceSetStart);
    free(SourceSetCount);
    free(SinkSourceInteractionList);

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
    int offset = 0;
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
                X[p+offset] = pos[p].x;
                Y[p+offset] = pos[p].y;
                Z[p+offset] = pos[p].z;
            }

            #pragma simd
            for(int p=offset;p<count + offset;p++){
                X[p] += newsinkcenter.x;
                Y[p] += newsinkcenter.y;
                Z[p] += newsinkcenter.z;
            }
            offset+=count;
        }
    }
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
}

int SetInteractionCollection::CheckCompletion(){
    return CompletionFlag;
}
