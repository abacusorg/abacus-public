/* slabbuffer.cpp

This defines all of the types of slabs that we will use
It provides routines to give the read and write filenames
as well as the sizes of these arenas.

It provides additional wrappings of SlabAllocator that depend
on file names or sizes.

The class SlabBuffer is defined here.  The global SB object
is an instance of this class and is the primary interface
for slab memory allocation and IO.

*/ 

#ifndef INCLUDE_SLABBUFFER
#define INCLUDE_SLABBUFFER

#include "particlestruct.cpp"
#include "arenaallocator.cpp"
#include "io_interface.h"

enum SlabType { CellInfoSlab,           //0
                PosSlab,                //1
                PosXYZSlab,             //2
                VelSlab,                //3
                AuxSlab,                //4
                AccSlab,                //5
                MergeCellInfoSlab,      //6
                MergePosSlab,           //7
                MergeVelSlab,           //8
                MergeAuxSlab,           //9
                TaylorSlab,             //10
                MultipoleSlab,          //11
                InsertCellInfoSlab,     //12
                NearAccSlab,            //13
                FarAccSlab,             //14
                FieldTimeSlice,         //15
                FieldTimeSlicePIDs,     //16
                VelLPTSlab,             //17
                CellGroupArena,         //18
                NearField_SIC_Slab,     //19
                
                L1halosSlab,            //20
                HaloRVSlabA,            //21
                HaloRVSlabB,            //22
                HaloPIDsSlabA,          //23
                HaloPIDsSlabB,          //24
                FieldRVSlabA,           //25
				FieldRVSlabB,           //26
                FieldPIDSlabA,          //27
                FieldPIDSlabB,          //28
                L0TimeSlice,            //29
                L0TimeSlicePIDs,        //30
                
                LightCone0,             //31
                LightCone1,             //32
                LightCone2,             //33

                LightCone0PID,          //34
                LightCone1PID,          //35
                LightCone2PID,          //36

                ICSlab,                 //37

                NUMTYPES
                };

enum SlabIntent { READSLAB,
                  WRITESLAB,
                  INMEMORY
                };

class SlabBuffer {
private:
    ArenaAllocator *AA;
    int order;                // The multipole order, used to determine ArenaSizes
    int cpd;

    // AA manages arenas using IDs, not slab types, so provide the mapping
    inline int TypeSlab2ID(int type, int s) {
        int ws = Grid->WrapSlab(s);
        // +1 for Reuse slab
        return type*(cpd+1) + ws;
    }
    
    // Every slab type has at most one "reuse slab"
    inline int ReuseID(int type) {
        return type*(cpd+1)+cpd;
    }

    // Check if this slab type is supposed to be directly allocated in shared memory
    int IsRamdiskSlab(int type, int hint=RAMDISK_AUTO);

public:

    SlabBuffer(int _cpd, int _order, uint64 max_allocations) {
        order = _order;
        cpd = _cpd;

        AA = new ArenaAllocator((_cpd+1)*NUMTYPES, max_allocations, P.UseMunmapThread, P.MunmapThreadCore);
    }

    ~SlabBuffer(void) {
        // Clean up the Reuse buffers
        for (int t = 0; t < NUMTYPES; t++) {
            int id = ReuseID(t);
            if (AA->IsArenaPresent(id)) {
                AA->DeAllocateArena(id, id);
            }
        }

        delete AA;
    }

    // Determine whether a slab is used for reading or writing by default
    int GetSlabIntent(int type);

    // Determine whether a slab is an output slab
    int IsOutputSlab(int type);
    // Determine whether a slab should have its checksum recorded
    int WantChecksum(int type);
    
    // The write and read names for the files of this SlabType.
    // The slab number will be wrapped.
    // Program will abort if one specifies a type not in this function;
    // this is a feature to block incorrect I/O.
    std::string WriteSlabPath(int type, int slab);
    std::string ReadSlabPath(int type, int slab);

    // The amount of memory to be allocated for the specified arena.
    uint64 ArenaSize(int type, int slab);
    
    char *AllocateArena(int type, int slab, int ramdisk = RAMDISK_AUTO);
    char *AllocateSpecificSize(int type, int slab, uint64 sizebytes, int ramdisk = RAMDISK_AUTO);

    // "Write" commands write the arena but do not delete it
    void WriteArenaBlockingWithoutDeletion(int type, int slab);
    void WriteArenaNonBlockingWithoutDeletion(int type, int slab);
    
    // "Store" commands will write the arena and then delete it.
    void StoreArenaBlocking(int type, int slab);
    void StoreArenaNonBlocking(int type, int slab);

    // "Read" commands read into an already-allocated arena
    void ReadArenaBlocking(int type, int slab);
    void ReadArenaNonBlocking(int type, int slab);
    
    // "Load" commands allocate the arena and then read it.
    void LoadArenaBlocking(int type, int slab);
    void LoadArenaNonBlocking(int type, int slab);

    // Use of these lower-level interfaces is discouraged in favor of the above.
    void ReadArena(int type, int slab, int blocking);
    void WriteArena(int type, int slab, int deleteafter, int blocking);
    void ReadArena(int type, int slab, int blocking, const char *fn);
    void WriteArena(int type, int slab, int deleteafter, int blocking, const char *fn);

    void DeAllocate(int type, int slab, int delete_file=0);

    inline char* GetSlabPtr(int type, int slab) {
        int id = TypeSlab2ID(type,slab);
        assertf(AA->IsArenaPresent(id), "Error: slab %d of type %d was not present, when pointer requested.\n", slab, type);
        return AA->GetArenaPtr(id);
    }

    int IsSlabPresent(int type, int slab) {
        int id = TypeSlab2ID(type,slab);
        return AA->IsArenaPresent(id);
    }

    void MarkSlabUnavailable(int type, int slab) {
        int id = TypeSlab2ID(type,slab);
        AA->MarkArenaUnavailable(id);
        return;
    }

    uint64 SlabSizeBytes(int type, int slab) {
        int id = TypeSlab2ID(type,slab);
        return AA->ArenaSizeBytes(id);
    }

    void ResizeSlab(int type, int slab, uint64 size) {
        STDLOG(2,"Resizing slab %d of type %d to size %l\n", slab, type, size);
        int id = TypeSlab2ID(type,slab);
        AA->ResizeArena(id, size);
    }

    int IsIOCompleted(int type, int slab) {
        int id = TypeSlab2ID(type,slab);
        if (!AA->IsArenaPresent(id)) return 0;
        return AA->IsIOCompleted(id);
    }

    void SetIOCompleted(int type, int slab) {
        int id = TypeSlab2ID(type,slab);
        AA->SetIOCompletedArena(id);
    }

    void GetMallocFreeTimes(double *malloc_time, double *free_time, double *disposal_thread_munmap){
        *malloc_time = AA->ArenaMalloc.Elapsed();
        //*free_time = AA->ArenaFree->Elapsed();
		*free_time = AA->ArenaFree_elapsed; 
        *disposal_thread_munmap = AA->DisposalThreadMunmap.Elapsed();
    }

    void report(){
        AA->report();
    }
};


/** Decide whether the given slab lives on the ramdisk
  * If so, we want to inform ArenaAllocator of this with the RAMDISK_READSLAB or RAMDISK_WRITESLAB flags,
  * which will be returned by this function as appropriate.
  * The hint allows one to override normal assumptions about a slab being a "read type" or "write type".
  * Valid values of `hint` are the RAMDISK_* flags.
  */
int SlabBuffer::IsRamdiskSlab(int type, int hint){

    // It's still possible to write a slab to ramdisk even if it's not listed below
    // The IO module should detect that it's a ramdisk path and use fopen instead of DIO

    // Does the hint indicate that we already decided the fate of this slab?
    if(hint >= RAMDISK_NO && hint < RAMDISK_AUTO)
        return hint;

    // If not, compare the file path for the slab to the system ramdisk path

    const int ReadHint = NUMTYPES;
    const int WriteHint = NUMTYPES+1;

    int intent_or_hint = GetSlabIntent(type);
    assert(intent_or_hint < ReadHint);  // for sanity

    switch(hint){
        case RAMDISK_AUTO_READSLAB:
            intent_or_hint = ReadHint;
            break;
        case RAMDISK_AUTO_WRITESLAB:
            intent_or_hint = WriteHint;
            break;
        case RAMDISK_AUTO:
            break;
        default:
            QUIT("Did not understand ramdisk hint %d\n", hint);
    }

    switch(intent_or_hint){
        case READSLAB:
        case ReadHint:
            if(is_path_on_ramdisk(ReadSlabPath(type, -1)))
                return RAMDISK_READSLAB;
            break;

        case WRITESLAB:
        case WriteHint:
            if(is_path_on_ramdisk(WriteSlabPath(type, -1)))
                return RAMDISK_WRITESLAB;
            break;

        // Absence of a slab type from both lists means that it should not explicitly use shared memory
        // An example is TimeSlice, which usually requires slab resizing, which is not trivial with shared memory
        // If the TimeSlice path lands on /dev/shm anyway, io_thread knows to use fopen instead of direct IO,
        // so it should be safe, just slower.
        case INMEMORY:
        default:
            break;
    }

    return RAMDISK_NO;
}

/* Determine whether a given slab type (e.g. PosSlab) is something that we typically read from disk
 * or write to disk (or neither).
 *
 * There are at least two scenarios when we need this information: the ramdisk code needs to know
 * whether to overwrite an existing file or not, and the manifest code needs to know whether a slab
 * has a corresponding file on disk.
*/
int SlabBuffer::GetSlabIntent(int type){
    switch(type){
        case CellInfoSlab:
        case PosSlab:
        case VelSlab:
        case AuxSlab:
        case TaylorSlab: 
            return READSLAB;

        case MultipoleSlab:
        case MergeCellInfoSlab:
        case MergePosSlab:
        case MergeVelSlab:
        case MergeAuxSlab:
            return WRITESLAB;

        default:
            return INMEMORY;
    }

    return -1;
}

/* Determine whether a given slab type (e.g. TimeSliceSlab) is an output type
 * that we want to record the checksum of when we write it to disk.
 *
 * P.NoChecksum disables all checksumming.
*/
int SlabBuffer::IsOutputSlab(int type){
   switch(type){
        case L0TimeSlice:
        case FieldTimeSlice:
        case L1halosSlab:
        case HaloRVSlabA:
        case HaloRVSlabB:
        case HaloPIDsSlabA:
        case HaloPIDsSlabB:
        case FieldRVSlabA:
        case FieldRVSlabB:
        case FieldPIDSlabA:
        case FieldPIDSlabB:
        case L0TimeSlicePIDs:
        case FieldTimeSlicePIDs:

        // Do we want to checksum light cones?  They're probably being merged/repacked very soon.
        // Maybe same for the halos?  But if we move either of these to another filesystem
        // for post-processing, then we'd like to be able to verify the checksums.
        case LightCone0:
        case LightCone1:
        case LightCone2:
            return 1;
     }
     return 0;
}

int SlabBuffer::WantChecksum(int type){
    if(P.NoChecksum)
        return 0;
    return IsOutputSlab(type);
}


std::string SlabBuffer::WriteSlabPath(int type, int slab) {
    slab = Grid->WrapSlab(slab);

    char slabnum[8]; sprintf(slabnum,"%04d",slab);
    char stepnum[8]; sprintf(stepnum,"%04d",ReadState.FullStepNumber);
    char redshift[16]; sprintf(redshift, "%5.3f", ReadState.Redshift);
    std::stringstream ss;
    std::string s;

    switch(type) {
#ifdef PARALLEL
        case TaylorSlab     : {
            // Read odd taylors from TaylorDirectory2, if it's defined
            if (WriteState.StripeConvState && slab % 2 == 1) {
                ss << P.TaylorDirectory2 << "/Taylor_" << slabnum;
                break;
            }
            ss << P.TaylorDirectory << "/Taylor_"     << slabnum; break;
        }
#endif
        case MultipoleSlab       : {
						
			
            // Send odd multipoles to MultipoleDirectory2, if it's defined
            if (WriteState.StripeConvState && slab % 2 == 1) {
                ss << P.MultipoleDirectory2 << "/Multipoles_" << slabnum;
                break;
            }
            ss << P.MultipoleDirectory << "/Multipoles_" << slabnum; break;
						
        }
        case MergeCellInfoSlab   : { ss << P.LocalWriteStateDirectory << "/cellinfo_"   << slabnum; break; }
        case MergePosSlab        : { ss << P.LocalWriteStateDirectory << "/position_"   << slabnum; break; }
        case MergeVelSlab        : { ss << P.LocalWriteStateDirectory << "/velocity_"   << slabnum; break; }
        case MergeAuxSlab        : { ss << P.LocalWriteStateDirectory << "/auxillary_"  << slabnum; break; }
        case AccSlab             : { ss << P.OutputDirectory << "/acc_"            << slabnum; break; }
        case NearAccSlab         : { ss << P.OutputDirectory << "/nearacc_"        << slabnum; break; }
        case FarAccSlab          : { ss << P.OutputDirectory << "/faracc_"         << slabnum; break; }
       
        case L1halosSlab           : { ss << P.GroupDirectory << "/Step" << stepnum << "_z" << redshift << "/halo_info_"           << slabnum; break;}
		
        case HaloPIDsSlabA    : { ss << P.GroupDirectory << "/Step" << stepnum << "_z" << redshift << "/halo_pids_A_"    << slabnum; break;}
        case HaloPIDsSlabB    : { ss << P.GroupDirectory << "/Step" << stepnum << "_z" << redshift << "/halo_pids_B_"    << slabnum; break;}
        case FieldPIDSlabA    : { ss << P.GroupDirectory << "/Step" << stepnum << "_z" << redshift << "/field_pids_A_"   << slabnum; break;}
        case FieldPIDSlabB    : { ss << P.GroupDirectory << "/Step" << stepnum << "_z" << redshift << "/field_pids_B_"   << slabnum; break;}
		
        case HaloRVSlabA : { ss << P.GroupDirectory << "/Step" << stepnum << "_z" << redshift << "/halo_rv_A_"      << slabnum; break;}
        case HaloRVSlabB : { ss << P.GroupDirectory << "/Step" << stepnum << "_z" << redshift << "/halo_rv_B_"      << slabnum; break;}
        case FieldRVSlabA       : { ss << P.GroupDirectory << "/Step" << stepnum << "_z" << redshift << "/field_rv_A_"     << slabnum; break;}
        case FieldRVSlabB       : { ss << P.GroupDirectory << "/Step" << stepnum << "_z" << redshift << "/field_rv_B_"     << slabnum; break;}

        case L0TimeSlice          : { ss << P.OutputDirectory << "/slice" << redshift << "/" << P.SimName << ".z" << redshift << ".slab" << slabnum << ".L0_pack9.dat"; break; }
        case FieldTimeSlice            : { ss << P.OutputDirectory << "/slice" << redshift << "/" << P.SimName << ".z" << redshift << ".slab" << slabnum << ".field_pack9.dat"; break; }
        case L0TimeSlicePIDs      : { ss << P.OutputDirectory << "/slice" << redshift << "/" << P.SimName << ".z" << redshift << ".slab" << slabnum << ".L0_pack9_pids.dat"; break; }
        case FieldTimeSlicePIDs        : { ss << P.OutputDirectory << "/slice" << redshift << "/" << P.SimName << ".z" << redshift << ".slab" << slabnum << ".field_pack9_pids.dat"; break; }
		
        default:
            QUIT("Illegal type %d given to WriteSlabPath()\n", type);
    }

    s = ss.str();
    return s;
}

std::string SlabBuffer::ReadSlabPath(int type, int slab) {
    slab = Grid->WrapSlab(slab);

    char slabnum[8]; sprintf(slabnum,"%04d",slab);
    std::stringstream ss;
    std::string s;

    switch(type) {
        case TaylorSlab     : {
            // Read odd taylors from TaylorDirectory2, if it's defined
            if (WriteState.StripeConvState && slab % 2 == 1) {
                ss << P.TaylorDirectory2 << "/Taylor_" << slabnum;
                break;
            }
            ss << P.TaylorDirectory << "/Taylor_"     << slabnum; break;
        }

        // The Merge slabs are usually write slabs, but can be used as read slabs in multipole recovery mode
        case MergeCellInfoSlab :
        case CellInfoSlab  : { ss << P.LocalReadStateDirectory << "/cellinfo_"   << slabnum; break; }

        case MergePosSlab :
        case PosSlab       : { ss << P.LocalReadStateDirectory << "/position_"   << slabnum; break; }

        case VelSlab       : { ss << P.LocalReadStateDirectory << "/velocity_"   << slabnum; break; }
        case AuxSlab       : { ss << P.LocalReadStateDirectory << "/auxillary_"  << slabnum; break; }
        case VelLPTSlab    : { ss << P.InitialConditionsDirectory << "/ic_" << slab; break; }
        case FieldTimeSlice     : { ss << WriteSlabPath(type, slab); }  // used for standalone FOF

        case ICSlab    : { ss << P.InitialConditionsDirectory << "/ic_" << slab; break; }

        default:
            QUIT("Illegal type %d given to ReadSlabPath()\n", type);
    }

    s = ss.str();
    return s;
}

uint64 SlabBuffer::ArenaSize(int type, int slab) {
    int rml = (order+1)*(order+1);
    uint64 lcpd = cpd;  // Just to assure that we don't spill 32-bits

    switch(type) {
        case CellInfoSlab        : { return sizeof(cellinfo)*lcpd*lcpd; }
        case MergeCellInfoSlab   : { return sizeof(cellinfo)*lcpd*lcpd; }
        case InsertCellInfoSlab  : { return sizeof(cellinfo)*lcpd*lcpd; }
        case MultipoleSlab  : { return lcpd*(lcpd+1)/2*rml*sizeof(MTCOMPLEX); }
        case TaylorSlab     : { return lcpd*(lcpd+1)/2*rml*sizeof(MTCOMPLEX); }
        case PosXYZSlab     :  // let these fall through to PosSlab
        case PosSlab        : { 
            return SS->size(slab)*sizeof(posstruct);
        }
        case MergePosSlab   : { 
            return SS->size(slab)*sizeof(posstruct);
        }
        case VelSlab        : { 
            return SS->size(slab)*sizeof(velstruct);
        }
        case AuxSlab        : { 
            return SS->size(slab)*sizeof(auxstruct);
        }
        case AccSlab        : {
            return SS->size(slab)*sizeof(accstruct);
        }
        case FarAccSlab     : {
            return SS->size(slab)*sizeof(acc3struct);
        }
        case NearAccSlab    : {
            return SS->size(slab)*sizeof(accstruct);
        }
        case ICSlab     : {
            return fsize(ReadSlabPath(ICSlab,slab).c_str());
        }
        case VelLPTSlab : {
            return ICFile::FromFormat(P.ICFormat, slab)->Npart*sizeof(velstruct);
        }
        case FieldTimeSlice : {
            return fsize(ReadSlabPath(FieldTimeSlice,slab).c_str());
        }
        
        /* // Ideally this is how we would allocate group finding arenas
                // but there's an annoying circular depdendency: we don't know about GFC yet, but GFC has fields that require slabbuffer.
                // Not easily solved with a forward declaration.
                // TODO: better way to do this?
                case L1halosSlab           : { return GFC->globalslabs[slab]->L1halos.get_slab_bytes(); }
        case TaggedPIDsSlab        : { return GFC->globalslabs[slab]->TaggedPIDs.get_slab_bytes(); }
        case L1ParticlesSlab       : { return GFC->globalslabs[slab]->L1Particles.get_slab_bytes(); }
        case HaloPIDsSlab            : { return GFC->globalslabs[slab]->L1PIDs.get_slab_bytes(); }
        case TaggableFieldSlab     : { 
            uint64 maxsize = P.np*P.HaloTaggableFraction*1.05;  // TODO: better heuristic? what will happen in very small sims?  Also technically HaloTaggableFraction is only used in the IC step
            return maxsize*sizeof(RVfloat);
        }
        case TaggableFieldPIDSlab  : {
            uint64 maxsize = P.np*P.HaloTaggableFraction*1.05;
            return maxsize*sizeof(TaggedPID);
        }*/

        default            : QUIT("Illegal type %d given to ArenaSize()\n", type);
    }
    return -1; //should be unreachable
}

char *SlabBuffer::AllocateArena(int type, int slab, int ramdisk) {	
    slab = Grid->WrapSlab(slab);
    uint64 s = ArenaSize(type, slab);
    return AllocateSpecificSize(type, slab, s, ramdisk);
}

char *SlabBuffer::AllocateSpecificSize(int type, int slab, uint64 sizebytes, int ramdisk) {
    // Most slabs are happy with RAMDISK_AUTO
    ramdisk = IsRamdiskSlab(type, ramdisk);
		
    std::string spath;
    switch(ramdisk){
        case RAMDISK_READSLAB:
            spath = ReadSlabPath(type, slab);
            break;
        case RAMDISK_WRITESLAB:
            spath = WriteSlabPath(type, slab);
            break;
        case RAMDISK_NO:
            spath = "";
            break;
        default:
            QUIT("Unexpected value %d of ramdisk\n", ramdisk);
            break;
    }
	const char *ramdisk_fn = spath.c_str();
	STDLOG(3, "Allocating slab %d of type %d to size %l (ramdisk = %d), total %5.3f GB\n",
                slab, type, sizebytes, ramdisk, AA->total_allocation/1024./1024./1024.);

    int id = TypeSlab2ID(type,slab);
    AA->Allocate(id, sizebytes, ReuseID(type), ramdisk, ramdisk_fn);

    return GetSlabPtr(type, slab);
}

void SlabBuffer::WriteArenaBlockingWithoutDeletion(int type, int slab) {
    // This is writing without deletion
    WriteArena(type, slab, IO_KEEP, IO_BLOCKING);
}

void SlabBuffer::WriteArenaNonBlockingWithoutDeletion(int type, int slab) {
    // This is writing without deletion
    WriteArena(type, slab, IO_KEEP, IO_NONBLOCKING);
}

void SlabBuffer::StoreArenaBlocking(int type, int slab) {
    // This is writing with deletion
    WriteArena(type, slab, IO_DELETE, IO_BLOCKING);
}

void SlabBuffer::StoreArenaNonBlocking(int type, int slab) {
    // This is writing with deletion
    WriteArena(type, slab, IO_DELETE, IO_NONBLOCKING);
}

void SlabBuffer::WriteArena(int type, int slab, int deleteafter, int blocking){
    std::string spath = WriteSlabPath(type,slab);
    const char *path = spath.c_str();
    // Determine the actual, allocated Ramdisk type of the current slab
    // If it was allocated on Ramdisk, then the writing is already done by definition!
    if(AA->ArenaRamdiskType(TypeSlab2ID(type,slab)) != RAMDISK_NO){
        STDLOG(2, "Skipping explicit write of type %d Ramdisk slab \"%s\"\n", type, path);

        // still might have to deallocate
        if(deleteafter == IO_DELETE){
#ifdef PARALLEL
			if (type == MultipoleSlab)
                STDLOG(1, "Uh oh. DeAllocating Multipole slab %d prematurely\n", slab);
#endif
            DeAllocate(type, slab);
		}
				
        return;
    }	

    WriteArena(type, slab, deleteafter, blocking, path);
}


void SlabBuffer::ReadArenaBlocking(int type, int slab) {
    ReadArena(type, slab, IO_BLOCKING);
}

void SlabBuffer::ReadArenaNonBlocking(int type, int slab) {
    ReadArena(type, slab, IO_NONBLOCKING);
}

void SlabBuffer::ReadArena(int type, int slab, int blocking){
    std::string spath = ReadSlabPath(type,slab);
    const char *path = spath.c_str();

    // Determine the actual, allocated Ramdisk type of the current slab
    // If it was allocated from existing shared memory (i.e. RAMDISK_READSLAB),
    // we probably don't want to read into it
    if(AA->ArenaRamdiskType(TypeSlab2ID(type,slab)) == RAMDISK_READSLAB){
        STDLOG(2, "Skipping explicit read of Ramdisk slab %d \"%s\"\n", slab, path);
        SetIOCompleted(type, slab);
        return;
    }

    ReadArena(type, slab, blocking, path);
}

void SlabBuffer::LoadArenaBlocking(int type, int slab) {
    AllocateArena(type, slab, RAMDISK_AUTO_READSLAB);
    ReadArenaBlocking(type, slab);

    assertf(SlabSizeBytes(type, slab) == fsize(ReadSlabPath(type,slab).c_str()),
            "Unexpected file size for slab %d of type %d\n", slab, type);
}

void SlabBuffer::LoadArenaNonBlocking(int type, int slab) {
    AllocateArena(type, slab, RAMDISK_AUTO_READSLAB);
    ReadArenaNonBlocking(type, slab);

    assertf(SlabSizeBytes(type, slab) == fsize(ReadSlabPath(type,slab).c_str()),
            "Unexpected file size for slab %d of type %d\n", slab, type);
}


void SlabBuffer::ReadArena(int type, int slab, int blocking, const char *fn) {
    // This will read into an arena.
    // This always triggers a real read, even if the path is on the ramdisk!
    // The read is always ordered to be the full usable size of the arena.

    if (P.ForceBlockingIO!=0)
        blocking = IO_BLOCKING;

    STDLOG(1,"Reading slab %d of type %d from file %s with blocking %d.\n", slab, type, fn, blocking);
    assertf(IsSlabPresent(type, slab), 
        "Type %d and Slab %d doesn't exist\n", type, slab);
    assertf(FileExists(fn),
        "File %s does not exist\n", fn);
    assertf(fsize(fn) >= 0 && (uint64) fsize(fn) >= SlabSizeBytes(type,slab),
        "File %s appears to be too small (%d bytes) compared to the arena (%d)\n",
            fn, fsize(fn), SlabSizeBytes(type,slab));
    
    ReadFile( GetSlabPtr(type,slab), 
        SlabSizeBytes(type,slab),
        type, slab,
        fn, 
        0,      // No file offset
        blocking);
}

void SlabBuffer::WriteArena(int type, int slab, int deleteafter, int blocking, const char *fn) {
    // This will write into an arena.
    // This always triggers a real write, even if the path is on the ramdisk!

    if (P.ForceBlockingIO!=0)
        blocking = IO_BLOCKING;

    STDLOG(1,"Writing slab %d of type %d to file %s with blocking %d and delete status %d.\n",
        slab, type, fn, blocking, deleteafter);
    assertf(IsSlabPresent(type, slab), 
        "Type %d and Slab %d doesn't exist\n", type, slab);

    WriteFile( GetSlabPtr(type,slab), 
        SlabSizeBytes(type,slab),
        type, slab,
        fn, 
        0,      // No file offset
        deleteafter, 
        blocking,
        WantChecksum(type));
}

// Free the memory associated with this slab.  Might stash the arena as a reuse slab!
// One can request deletion of the file that corresponds to this slab, too.
void SlabBuffer::DeAllocate(int type, int slab, int delete_file) {
    STDLOG(2,"Deallocating slab %d of type %d.\n", slab, type);
    AA->DeAllocateArena(TypeSlab2ID(type,slab), ReuseID(type));

    if(delete_file){
        int intent = GetSlabIntent(type);
        std::string path;
        switch(intent){
            case READSLAB:
                path = ReadSlabPath(type, slab);
                break;
            case WRITESLAB:
                path = WriteSlabPath(type, slab);
                break;
            case INMEMORY:
                return;  // No file; nothing to delete from disk!
            default:
                QUIT("Did not understand slab intent %d\n", intent);
        }

        char buffer[1024];
        strcpy(buffer, path.c_str());
        char *tmp = dirname(buffer);
					
        ExpandPathName(tmp);
				

        if(!IsTrueLocalDirectory(tmp)){
            STDLOG(2,"Not deleting slab file \"%s\" because it is in a global directory\n", path);
            return;
        }

        STDLOG(2,"Deleting slab file \"%s\"\n", path);
        int ret = unlink(path.c_str());
        if(ret != 0){
            assertf(errno == ENOENT, "Failed to remove path \"%s\" for reason %d: %s\n", path, errno, strerror(errno));
            // TODO: is it really safe to fail to remove the file?
            STDLOG(2, "Failed to remove path \"%s\"; does not exist. Continuing.\n", path);
        }
    }
}

#endif // INCLUDE_SLABBUFFER
