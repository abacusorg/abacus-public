// This file contains parts of directdriver.cpp that are relevant
// to CPU computation and comparison of the GPU & CPU code.
// These routine are not normally executed.

void NearFieldDriver::CheckGPUCPU(int slabID){
// Computes the CPU result to compare to the GPU result
// but does not overwrite the GPU forces.
    size_t len = Slab->size(slabID) *sizeof(accstruct);
    accstruct * a_cpu = (accstruct *)malloc(len);
    accstruct * a_tmp = (accstruct *)malloc(len);
    accstruct * a_gpu = (accstruct *) LBW->ReturnIDPtr(AccSlab,slabID);
    auxstruct * aux = (auxstruct *)LBW->ReturnIDPtr(AuxSlab,slabID);
    memcpy(a_tmp,a_gpu,len);  // Save the GPU result in tmp
    memset(a_gpu,0,len);
    ExecuteSlabCPU(slabID);  // Compute the CPU result
    memcpy(a_cpu,a_gpu,len);  // Save the CPU result
    memcpy(a_gpu,a_tmp,len); // Load the GPU result
    #ifdef DOUBLEPRECISION
    FLOAT target = 1e-10;
    #else
    FLOAT target = 1e-1;
    #endif

    for(int i = 0; i < Slab->size(slabID);i++){
        acc3struct ai_g = TOFLOAT3(a_gpu[i]);
        acc3struct ai_c = TOFLOAT3(a_cpu[i]);
        if(ai_g.norm() == 0. && ai_c.norm() == 0.)
            continue;
        FLOAT delta =2* (ai_g-ai_c).norm()/(ai_g.norm() + ai_c.norm());
        if(!(delta < target)){
            printf("Error in slab %d:\n\ta_gpu[%d]: (%5.4f,%5.4f,%5.4f)\n\ta_cpu[%d]: (%5.4f,%5.4f,%5.4f)\n\tdelta:%f\n",
                    slabID,i,ai_g.x,ai_g.y,ai_g.z,i,ai_c.x,ai_c.y,ai_c.z,delta);
            //assert(delta < target);
        }
        assert(isfinite(ai_g.x));
        assert(isfinite(ai_g.y));
        assert(isfinite(ai_g.z));

        assert(isfinite(ai_c.x));
        assert(isfinite(ai_c.y));
        assert(isfinite(ai_c.z));
    }
    free(a_cpu);
    free(a_tmp);
    
    STDLOG(1,"GPU-CPU comparison passed for slab %d\n", slabID);
}

void NearFieldDriver::ExecuteSlabCPU(int slabID){
    ExecuteSlabCPU(slabID,(int *)NULL);
}


void NearFieldDriver::ExecuteSlabCPU(int slabID, int * predicate){
    CPUFallbackTimer.Start();
    if(!LBW->IDPresent(AccSlab, slabID))
        LBW->AllocateArena(AccSlab,slabID);
    ZeroAcceleration(slabID,AccSlab);
    
    #ifdef DIRECTSINGLESPLINE
    FLOAT inv_eps3 = 1./(SofteningLengthInternal*SofteningLengthInternal*SofteningLengthInternal);
    #endif

    uint64 DI_slab = 0;
    uint64 NSink_CPU_slab = 0;
    #pragma omp parallel for schedule(dynamic,1) reduction(+:DI_slab,NSink_CPU_slab)
    for(int y = 0; y < P.cpd; y++){
        int g = omp_get_thread_num();
        //STDLOG(1,"Executing directs on pencil y=%d in slab %d, in OMP thread %d on CPU %d (nprocs: %d)\n", y, slabID, g, sched_getcpu(), omp_get_num_procs());
        for(int z = 0; z < P.cpd; z++){
            if(predicate != NULL && !predicate[y*P.cpd +z]) continue;
            
            // We can use PosCell instead of PosXYZ here because it's for the sink positions
            posstruct * sink_pos = PP->PosCell(slabID,y,z);
            accstruct * sink_acc = PP->NearAccCell(slabID,y,z);
            uint64 np_sink = PP->NumberParticle(slabID,y,z);
            NSink_CPU_slab += np_sink;
            if (np_sink == 0) continue;
            for(int i = slabID - RADIUS; i <= slabID + RADIUS; i++){
                for(int j = y - RADIUS; j <= y + RADIUS; j++){
                    for(int k = z - RADIUS; k <= z + RADIUS; k++){
                        uint64 np_source = PP->NumberParticle(i,j,k);
                        
                        // We assume that PosXYZ is used for all sources
                        // This lets us drift PosSlab much earlier
                        List3<FLOAT> source_pos_xyz = PP->PosXYZCell(i,j,k);
                        posstruct *source_pos = new posstruct[np_source];
                        for(uint64 ii = 0; ii < np_source; ii++){
                            source_pos[ii].x = source_pos_xyz.X[ii];
                            source_pos[ii].y = source_pos_xyz.Y[ii];
                            source_pos[ii].z = source_pos_xyz.Z[ii];
                        }
                        
                        FLOAT3 delta = PP->CellCenter(slabID,y,z)-PP->CellCenter(i,j,k);
			// TODO: At present, the b2 parameter is not passed
			// into the CPU directs, so there is no FOF neighbor
			// computation.
			// TODO: sink_acc may now be a float4, but the CPU routines
			// want a float3.  We'll overload this space and fix it later
                        if(np_source >0) DD[g].AVXExecute(sink_pos,source_pos,np_sink,np_source,
                                delta,eps,(FLOAT3 *)sink_acc);
                        delete[] source_pos;
                        
                        DI_slab += np_sink*np_source;
                    }
                }
            }
	    // All done with this cell.  Fix the float4 to float3 issue
	    acc3struct *acc3 = (acc3struct *)sink_acc;
	    for (int64_t i=np_sink-1; i>=0; i--) 
	    	sink_acc[i] = accstruct(acc3[i]);
        }
    }
    DirectInteractions_CPU += DI_slab;
    NSink_CPU += NSink_CPU_slab;
    CPUFallbackTimer.Stop();
}

#include <vector>
#include <algorithm>
void NearFieldDriver::CheckInteractionList(int slab){

    vector<uint64> ** il = new std::vector<uint64> *[P.cpd*P.cpd];

    for(int i = 0; i < P.cpd*P.cpd; i++){
        il[i] = new vector<uint64>();
        il[i]->reserve(WIDTH*WIDTH*WIDTH);
    }

    SetInteractionCollection ** Slices = SlabInteractionCollections[slab];
    int NSplit = SlabNSplit[slab];

    for(int s = 0; s < WIDTH*NSplit; s++)
        Slices[s]->AddInteractionList(il);


    //check number of interactions
    for(int i = 0; i < P.cpd*P.cpd; i++){
        if (il[i]->size() != WIDTH*WIDTH*WIDTH){
            assertf(il[i]->size() == WIDTH*WIDTH*WIDTH,
            "Error: cell %d has %ld =/= %d interactions\n",
	    i,il[i]->size(),WIDTH*WIDTH*WIDTH);
        }
    }

    //Now we do a full check of the interaction list
    
    //build a correct interaction list for comparison
    
    vector<uint64> ** il_test = new std::vector<uint64> *[P.cpd*P.cpd];
    for(int i = 0; i < P.cpd*P.cpd; i++){
        il_test[i] = new vector<uint64>();
        il_test[i]->reserve(WIDTH*WIDTH*WIDTH);
    }
    for(int y = 0; y < P.cpd; y++){
        for(int z = 0; z < P.cpd; z++){
            for(int i = slab - RADIUS; i <= slab + RADIUS; i++){
                for(int j = y - RADIUS; j <= y + RADIUS; j++){
                    for(int k = z - RADIUS; k <= z + RADIUS; k++){
                        uint64 sinkCellId = y*P.cpd + z;
                        uint64 sourceCellId = P.cpd*P.cpd * PP->WrapSlab(i) +
                            P.cpd*PP->WrapSlab(j) +
                            PP->WrapSlab(k);
                        il_test[sinkCellId]->push_back(sourceCellId);
                    }
                }
            }
        }
    }
    
    
    //sort the lists and check them
    for(int i=0; i <P.cpd*P.cpd; i++){
        std::sort(il[i]->begin(),il[i]->end());
        std::sort(il_test[i]->begin(),il_test[i]->end());

        assert(il[i]->size() == il_test[i]->size());

        for(int j = 0; j < il[i]->size(); j ++)
                assertf(il[i]->at(j) == il_test[i]->at(j), "Interaction list %d at %d failed in slab %d\n", i, j, slab);
        delete il[i];
        delete il_test[i];
    }

    delete[] il;
    delete[] il_test;
    
    STDLOG(1,"Checking the interaction list for slab %d passed.\n", slab);
    //Checking the interaction list is time-consuming
    //it should never be done in an actual run
}
