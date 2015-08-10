//Benchmark the direct module using "realistic" fake slabs
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <pthread.h>
#include <algorithm>
#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>

#include <cmath>


#include "threevector.hh"
#include "STimer.cc"

#ifdef DOUBLEPRECISION
#define FLOAT double
#define FLOAT4 double4
#define FLOAT3 double3
#else
#define FLOAT float
#define FLOAT4 float4
#define FLOAT3 float3
#define posstruct FLOAT3
#define accstruct FLOAT3
#endif

#ifdef CUDADIRECT
#include "DeviceFunctions.h"
#endif

#include "header.cpp"
#include "direct.h"
#include "direct.cpp"


#ifndef NFRADIUS
    #define NFRADIUS 2
#endif
#define NFR NFRADIUS

#define WIDTH (2*NFRADIUS+1)



#define check_tol(a,b,tol)  (!(fabs(a-b)/((a+b) +1e-10) < tol))




inline int _wrap(int i,int cpd) {
    if(i>=cpd) i-=cpd;
    if(i<0) i+= cpd; 
    return i;
}   
    
#define wrap(a) _wrap(a,cpd)

long long int CreateInput(int cpd, double mean_ppc, int slab,
            FLOAT3  **   pos, 
            int     **   c_start,
            int     **   c_np){
    //set up the rng
    const gsl_rng_type * T;
            
    gsl_rng_env_setup();
    T = gsl_rng_default; 
    int seed = 12345 + slab;
    gsl_rng * rng = gsl_rng_alloc (T);
    gsl_rng_set(rng,seed);
    
    gsl_rng * pos_rng = gsl_rng_alloc(T);
    gsl_rng_set(pos_rng,seed);
    
    long long int total = 0;
    int maxppc = 0;
    int minppc = mean_ppc;
    for(long long int i = 0; i < cpd*cpd;i++){
        (c_start[slab])[i] = total;
        (c_np[slab])[i] = gsl_ran_poisson(rng,mean_ppc);
        maxppc = std::max((c_np[slab])[i],maxppc);
        minppc = std::min((c_np[slab])[i],minppc);
        int y = i/cpd;
        int z = i % cpd;
        total+= (c_np[slab])[i];
    }   
    pos[slab] = (FLOAT3 *) malloc(total*sizeof(FLOAT3));
        
    for(long long int i = 0; i <cpd*cpd;i++){
        int y = i/cpd;
        int z = i % cpd;
        for(int j = (c_start[slab])[i]; j < (c_start[slab])[i] + (c_np[slab])[i]; j++){
            (pos[slab])[j].x = gsl_rng_uniform(pos_rng)*1.0f/cpd;
            (pos[slab])[j].y = gsl_rng_uniform(pos_rng)*1.0f/cpd;
            (pos[slab])[j].z = gsl_rng_uniform(pos_rng)*1.0f/cpd;
        }   
    }       
            
    printf("\t\tTotal NP for slab %d: %lld | Max PPC: %d | Min PPC %d\n",slab,total,maxppc,minppc);
    gsl_rng_free(pos_rng);
    gsl_rng_free(rng);
    return total;    
}

long long unsigned int CalcDI(int ** c_n,long long unsigned int cpd, int slab){
    long long unsigned  int di = 0;
    for(int j = 0; j < cpd; j++){
        for(int k = 0; k < cpd; k++){
            long long unsigned int nsinks = (c_n[slab])[j*cpd+k];
            for(int dx = -NFR; dx <=NFR; dx++){
                for(int dy = -NFR; dy <= NFR; dy++){
                    for(int dz = -NFR; dz <= NFR; dz++){
                        int sourceidx = wrap(j+dy)*cpd + wrap(k+dz);
                        long long unsigned int nsources = (c_n[wrap(slab + dx)])[sourceidx];
                        di+= nsources*nsinks; 
                    }}}
        }}              
    return di;          
}

Direct *DD;

void ExecuteSlabCPU(int slabID, int * predicate,
        FLOAT3 ** pos, int** c_start, int** c_np,
        int cpd, FLOAT3 * acc, FLOAT eps){
    long long unsigned int DI_slab = 0;
    #pragma omp parallel for schedule(dynamic,1) reduction(+:DI_slab)
    for(int y = 0; y < cpd; y++){
        int g = omp_get_thread_num();
        for(int z = 0; z < cpd; z++){
            if(predicate != NULL && !predicate[y*cpd +z]) continue;
            int sink_offset = (c_start[slabID])[y*cpd+z];
            FLOAT3 * sink_pos = &((pos[slabID])[sink_offset]);
            FLOAT3 * sink_acc = &(acc[sink_offset]);
            uint64 np_sink = (c_np[slabID])[y*cpd+z];
            if (np_sink == 0) continue;
            for(int i =             slabID  - NFRADIUS; i <= slabID + NFRADIUS; i++){
                for(int j =         y       - NFRADIUS; j <= y      + NFRADIUS; j++){
                    for(int k =     z       - NFRADIUS; k <= z      + NFRADIUS; k++){
                        int wi = wrap(i);
                        int wj = wrap(j);
                        int wk = wrap(k);
                        int source_offset = (c_start[wi])[wj*cpd+wk];
                        FLOAT3 *source_pos = &((pos[wi])[source_offset]);
                        FLOAT3 delta;
                        delta.x = (1.0/cpd) *(slabID-i);
                        delta.y = (1.0/cpd)*(y-j);
                        delta.z = (1.0/cpd)*(z-k);
                        int np_source = (c_np[wi])[wj*cpd+wk];
                        if(np_source >0) DD[g].AVXExecute(sink_pos,source_pos,np_sink,np_source,
                                delta,eps*eps,sink_acc);
                        DI_slab+=np_sink*np_source;
                    }   
                }               
            }           
        }           
    }
    //if(predicate == NULL) printf("\t\tCPU reports %lld interactions for slab %d\n",DI_slab,slabID);

}
int slabcomplete[WIDTH];

int GPUMinCellSinks = 0;

void ExecuteSlabGPU(FLOAT3 ** pos, int ** c_start, int ** c_np,int cpd,int slabID, int *NPslab, FLOAT3 * acc,FLOAT eps){
    #ifdef CUDADIRECT
    int * CellStart[WIDTH];
    int * CellNP[WIDTH];
    #ifdef AVXDIRECT
    int pred[cpd*cpd];
    #else
    int * pred = NULL;
    #endif
    int SlabNP[WIDTH];
    int SlabIDs[WIDTH];
    FLOAT3 *SlabPos[WIDTH];

    for(int i = 0; i < WIDTH; i++){
        int sid = wrap(slabID+i-NFRADIUS);
        SlabIDs[i] = sid;
        SlabNP[i] = NPslab[sid];
        SlabPos[i] = pos[sid];
        PinSlab(SlabPos[i],sid,SlabNP[i],cpd);

        CellStart[i] = GetCellStart(cpd,i);
        CellNP[i] = GetCellNP(cpd,i);
        #pragma omp parallel for schedule(dynamic,1)
        for(int y = 0; y <cpd*cpd;y++){
            (CellStart[i])[y] = (c_start[sid])[y];
            (CellNP[i])[y] = (c_np[sid])[y];
            #ifdef AVXDIRECT
            pred[y] = ((CellNP[i])[y] < GPUMinCellSinks);
            #endif
        }
    }
    SetupGPU(SlabIDs,SlabPos,SlabNP,CellStart,CellNP,cpd);

    DeviceAcceleration(acc, slabID, NPslab[slabID],
            cpd, eps*eps, slabcomplete, 0,pred);
    if(pred != NULL) ExecuteSlabCPU(slabID,pred,pos,c_start,c_np,cpd,acc,eps);
    #endif
}

void BenchmarkGPU(FLOAT3 ** pos, int ** c_start, int ** c_np, FLOAT3 ** acc,
        int cpd, int nslabs,
        int *NPslab,
        FLOAT eps ){
#ifdef CUDADIRECT
    printf("\tRunning GPU benchmark...\n");
    STimer wc;
    wc.Clear();
    wc.Start();

    for(int i = NFR; i < nslabs - NFR; i++){

        ExecuteSlabGPU(pos,c_start,c_np,cpd,i,NPslab, acc[i], eps);
        UnpinSlab(pos[i],i);
    }
    long long unsigned int GPUDI = DeviceGetDI();
    //printf("\t\tGPU reports %lld DI for slabs %d - %d\n",GPUDI,NFR,nslabs-NFR-1);
    wc.Stop();
    printf("\t\tDone.\n");
#else
    printf("\tNot compiled with CUDA directs; skipping GPU benchmark.\n");
#endif
}

void BenchmarkCPU(FLOAT3 ** pos, int ** c_start, int ** c_np, FLOAT3 ** acc,
        int cpd, int nslabs,
        int *NPslab,
        FLOAT eps){

    printf("\tRunning CPU benchmark...\n");
    STimer wc;
    wc.Clear();
    wc.Start();

    for(int i = NFR; i < nslabs - NFR; i++)
        ExecuteSlabCPU(i,NULL,pos,c_start,c_np,cpd,acc[i],eps);
    wc.Stop();
}

void CheckGPUCPU(int nslab, int *NPslab, FLOAT3 ** a_gpu, FLOAT3 ** a_cpu){

    double delta_m = 0.0;
    int NP = 0;
    printf("\tCross-checking accelerations...");

    for(int i = NFR; i < nslab-NFR; i++){
        NP += NPslab[i];
        for(int n = 0; n < NPslab[i]; n++){
            FLOAT3 delta = (a_gpu[i])[n]-(a_cpu[i])[n];
            FLOAT3 sigma = (a_gpu[i])[n]+(a_cpu[i])[n];
            if( abs(delta.norm()/sigma.norm()) > 1e-1){
                printf("Bad value for particle %d in slab %d\n",n,i);
                printf("\tGPU: ( %f, %f, %f )\n", (a_gpu[i])[n].x,(a_gpu[i])[n].y,(a_gpu[i])[n].z);
                printf("\tCPU: ( %f, %f, %f )\n", (a_cpu[i])[n].x,(a_cpu[i])[n].y,(a_cpu[i])[n].z);
                printf("\tdel: ( %f, %f, %f )\n", delta.x,delta.y,delta.z);
                assert(delta.norm()/sigma.norm() <= 1e-1);
            }
            delta_m += abs(delta.norm()/sigma.norm());
        }
    }
    printf("\tPassed.\n\tMean relative error: %f\n", delta_m/NP);

}



int main(int argc, char ** argv){
    if(argc != 5){
        printf("Usage: %s [cpd] [mean ppc] [n slabs] [mode]\n",argv[0]);
        printf("\tModes:\n\t\t0: Run both, no check\n\t\t1: Run both, check\n\t\t2:GPU only\n\t\t3:CPUOnly\n");
        return 1;
    }
    int cpd = atoi(argv[1]);
    double ppc = atof(argv[2]);
    int slabs = atoi(argv[3]);
    int check = atoi(argv[4]);

    printf("Benchmarking directs with %d slabs (%d will be executed) of %d cells with a mean %f particles per cell\n",
            slabs,slabs - 2*NFR, cpd*cpd,ppc);
    FLOAT3  * pos[slabs];
    FLOAT3  * acc_cpu[slabs];
    FLOAT3 * acc_gpu[slabs];
    int     * c_start[slabs];
    int     * c_np[slabs];
    int       NPslab[slabs];
    long long unsigned int DIslab[slabs];
    int NPTotal = 0;
    long long unsigned int DITotal = 0;

    double GDirect;
    int nthread = omp_get_max_threads();
    DD = new Direct[nthread];
    
    printf("\tCreating Test slabs...\n");
    STimer input; input.Clear(); input.Start();
    for(int i = 0; i < slabs; i++)
    {
        c_start[i] = (int *)malloc(cpd*cpd*sizeof(int));
        c_np[i]    = (int *)malloc(cpd*cpd*sizeof(int));
        NPslab[i]  = CreateInput(cpd, ppc, i, pos, c_start, c_np);
    }

    for(int i = NFR; i < slabs - NFR; i ++){
        NPTotal += NPslab[i];
        if(check < 3){
            acc_gpu[i] = (FLOAT3 *)malloc(NPslab[i]*sizeof(FLOAT3));
            memset(acc_gpu[i],0,NPslab[i]*sizeof(FLOAT3));
        }
        if(check!= 2) {
            acc_cpu[i] = (FLOAT3 *)malloc(NPslab[i]*sizeof(FLOAT3));
            memset(acc_cpu[i],0,NPslab[i]*sizeof(FLOAT3));
        }
        DIslab[i] = CalcDI(c_np,cpd,i);
        DITotal += DIslab[i];
    }
    GDirect = DITotal/1e9;
    double Mp = NPTotal/1e6;
    double eps = .1/(ppc *cpd);

    input.Stop();
    printf("\tDone. Created %f million effective particles and %f Gdirect interactions at eps = %e\n\t\t\t--->%f Mp/s\n", Mp, GDirect,eps,Mp/input.Elapsed());


    STimer gpu; gpu.Clear(); gpu.Start();
    if(check < 3) 
        BenchmarkGPU(pos,c_start,c_np,acc_gpu,cpd,slabs,NPslab,eps);
    gpu.Stop();
    double gpu_time = gpu.Elapsed();

    if(check < 3) 
        printf("\tGPU: %f s for %f Gd/s (%f Mp/s)\n\n",gpu_time,GDirect/gpu_time,Mp/gpu_time);

    STimer cpu; cpu.Clear(); cpu.Start();
    if(check!= 2) 
        BenchmarkCPU(pos,c_start,c_np,acc_cpu,cpd,slabs,NPslab,eps);
    cpu.Stop();
    double cpu_time = cpu.Elapsed();

    if(check!= 2) 
        printf("\tCPU: %f s for %f Gd/s (%f Mp/s)\n\n",cpu_time,GDirect/cpu_time,Mp/cpu_time);

    if(check == 1) 
        CheckGPUCPU(slabs,NPslab,acc_gpu, acc_cpu);
}
