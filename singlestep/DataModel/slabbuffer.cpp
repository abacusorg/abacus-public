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

enum SlabType { CellInfoSlab,
                PosSlab,
                PosXYZSlab,
                VelSlab,
                AuxSlab,
                AccSlab,
                MergeCellInfoSlab,
                MergePosSlab,
                MergeVelSlab,
                MergeAuxSlab,
                TaylorSlab, 
                MultipoleSlab, 
                InsertCellInfoSlab,
                NearAccSlab,
                FarAccSlab,
                TimeSlice,
                VelLPTSlab,
                NearField_SIC_Slab, 
                
                L1halosSlab,
                TaggedPIDsSlab,
                L1ParticlesSlab,
                L1PIDsSlab,
                TaggableFieldSlab,
                TaggableFieldPIDSlab,
                L0TimeSlice,
                
                LightCone0,
                LightCone1,
                LightCone2,
                LightCone3,
                LightCone4,
                LightCone5,
                LightCone6,
                LightCone7,

                NUMTYPES
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

        AA = new ArenaAllocator((_cpd+1)*NUMTYPES, max_allocations);
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
    
    // The write and read names for the files of this SlabType.
    // The slab number will be wrapped.
    // Program will abort if one specifies a type not in this function;
    // this is a feature to block incorrect I/O.
    std::string WriteSlabPath(int type, int slab);
    std::string ReadSlabPath(int type, int slab);

    // The amount of memory to be allocated for the specified arena.
    uint64 ArenaSize(int type, int slab);
    
    char *AllocateArena(int type, int slab, int ramdisk = RAMDISK_AUTO);
    void AllocateSpecificSize(int type, int slab, uint64 sizebytes, int ramdisk = RAMDISK_AUTO);

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

    inline char* GetSlabPtr(int type, int slab) {
        int id = TypeSlab2ID(type,slab);
        return AA->GetArenaPtr(id);
    }

    int IsSlabPresent(int type, int slab) {
        int id = TypeSlab2ID(type,slab);
        return AA->IsArenaPresent(id);
    }

    uint64 SlabSizeBytes(int type, int slab) {
        int id = TypeSlab2ID(type,slab);
        return AA->ArenaSizeBytes(id);
    }

    void ResizeSlab(int type, int slab, uint64 size) {
        STDLOG(1,"Resizing slab %d of type %d to size %l\n", slab, type, size);
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

    void DeAllocate(int type, int slab) {
        STDLOG(1,"Deallocating slab %d of type %d.\n", slab, type);
        AA->DeAllocateArena(TypeSlab2ID(type,slab), ReuseID(type));
    }

    void GetMallocFreeTimes(double *malloc_time, double *free_time){
        *malloc_time = AA->ArenaMalloc.Elapsed();
        *free_time = AA->ArenaFree->Elapsed();
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
    assert(type < ReadHint);  // for sanity

    int type_or_hint = type;
    switch(hint){
        case RAMDISK_AUTO_READSLAB:
            type_or_hint = ReadHint;
            break;
        case RAMDISK_AUTO_WRITESLAB:
            type_or_hint = WriteHint;
            break;
        case RAMDISK_AUTO:
            break;
        default:
            QUIT("Did not understand ramdisk hint %d\n", hint);
    }

    switch(type_or_hint){
        // Here we express defaults about a slab being a "read type" or "write type"
        case TaylorSlab:
        case CellInfoSlab:
        case PosSlab:
        case VelSlab:
        case AuxSlab:
        case ReadHint:
            if(is_path_on_ramdisk(ReadSlabPath(type, -1)))
                return RAMDISK_READSLAB;
            break;

        case MultipoleSlab:
        case MergeCellInfoSlab:
        case MergePosSlab:
        case MergeVelSlab:
        case MergeAuxSlab:
        case WriteHint:
            if(is_path_on_ramdisk(WriteSlabPath(type, -1)))
                return RAMDISK_WRITESLAB;
            break;

        // Absence of a slab type from both lists means that it should not explicitly use shared memory
        // An example is TimeSlice, which usually requires slab resizing, which is not trivial with shared memory
        // If the TimeSlice path lands on /dev/shm anyway, io_thread knows to use fopen instead of direct IO,
        // so it should be safe, just slower.
        default:
            break;
    }

    return RAMDISK_NO;
}


std::string SlabBuffer::WriteSlabPath(int type, int slab) {
    slab = Grid->WrapSlab(slab);

    char slabnum[8]; sprintf(slabnum,"%04d",slab);
    char stepnum[8]; sprintf(stepnum,"%04d",ReadState.FullStepNumber);
    char redshift[16]; sprintf(redshift, "%5.3f", ReadState.Redshift);
    std::stringstream ss;
    std::string s;

    switch(type) {
        case MultipoleSlab       : {
            // Send odd multipoles to MultipoleDirectory2, if it's defined
            if (WriteState.StripeConvState && slab % 2 == 1) {
                ss << P.MultipoleDirectory2 << "/Multipoles_" << slabnum;
                break;
            }
            ss << P.MultipoleDirectory << "/Multipoles_" << slabnum; break;
        }
        case MergeCellInfoSlab   : { ss << P.WriteStateDirectory << "/cellinfo_"   << slabnum; break; }
        case MergePosSlab        : { ss << P.WriteStateDirectory << "/position_"   << slabnum; break; }
        case MergeVelSlab        : { ss << P.WriteStateDirectory << "/velocity_"   << slabnum; break; }
        case MergeAuxSlab        : { ss << P.WriteStateDirectory << "/auxillary_"  << slabnum; break; }
        case AccSlab             : { ss << P.OutputDirectory << "/acc_"            << slabnum; break; }
        case NearAccSlab         : { ss << P.OutputDirectory << "/nearacc_"        << slabnum; break; }
        case FarAccSlab          : { ss << P.OutputDirectory << "/faracc_"         << slabnum; break; }
        
        case L1halosSlab           : { ss << P.GroupDirectory << "/Step" << stepnum << "_z" << redshift << "/halos_"                   << slabnum; break;}
        case TaggedPIDsSlab        : { ss << P.GroupDirectory << "/Step" << stepnum << "_z" << redshift << "/tagged_pids_"             << slabnum; break;}
        case L1ParticlesSlab       : { ss << P.GroupDirectory << "/Step" << stepnum << "_z" << redshift << "/halo_particles_"          << slabnum; break;}
        case L1PIDsSlab            : { ss << P.GroupDirectory << "/Step" << stepnum << "_z" << redshift << "/halo_particle_ids_"       << slabnum; break;}
        case TaggableFieldSlab     : { ss << P.GroupDirectory << "/Step" << stepnum << "_z" << redshift << "/taggable_particles_"      << slabnum; break;}
        case TaggableFieldPIDSlab  : { ss << P.GroupDirectory << "/Step" << stepnum << "_z" << redshift << "/taggable_particle_ids_"   << slabnum; break;}

        case L0TimeSlice : { ss << P.OutputDirectory << "/slice" << redshift << "/" << P.SimName << ".z" << redshift << ".slab" << slabnum << ".L0.dat"; break; }
        case TimeSlice   : { ss << P.OutputDirectory << "/slice" << redshift << "/" << P.SimName << ".z" << redshift << ".slab" << slabnum << ".dat"; break; }

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
        case CellInfoSlab  : { ss << P.ReadStateDirectory << "/cellinfo_"   << slabnum; break; }

        case MergePosSlab :
        case PosSlab       : { ss << P.ReadStateDirectory << "/position_"   << slabnum; break; }

        case VelSlab       : { ss << P.ReadStateDirectory << "/velocity_"   << slabnum; break; }
        case AuxSlab       : { ss << P.ReadStateDirectory << "/auxillary_"  << slabnum; break; }
        case VelLPTSlab    : { ss << P.InitialConditionsDirectory << "/ic_" << slab; break; }
        case TimeSlice     : { ss << WriteSlabPath(type, slab); }  // used for standalone FOF

        default:
            QUIT("Illegal type %d given to ReadSlabPath()\n", type);
    }

    s = ss.str();
    return s;
}

uint64 SlabBuffer::ArenaSize(int type, int slab) {
    int rml = (order+1)*(order+1);
    uint64 lcpd = cpd;  // Just to assure that we don't spill 32-bits

    // TODO: can we use the SlabSize info instead of fsize for these?
    switch(type) {
        case CellInfoSlab        : { return sizeof(cellinfo)*lcpd*lcpd; }
        case MergeCellInfoSlab   : { return sizeof(cellinfo)*lcpd*lcpd; }
        case InsertCellInfoSlab  : { return sizeof(cellinfo)*lcpd*lcpd; }
        case MultipoleSlab  : { return lcpd*(lcpd+1)/2*rml*sizeof(MTCOMPLEX); }
        case TaylorSlab     : { return lcpd*(lcpd+1)/2*rml*sizeof(MTCOMPLEX); }
        case PosXYZSlab     :  // let these fall through to PosSlab
        case PosSlab        : { 
            return fsize(ReadSlabPath(PosSlab,slab).c_str());
        }
        case MergePosSlab   : { 
            return fsize(ReadSlabPath(MergePosSlab,slab).c_str());
        }
        case VelSlab        : { 
            return fsize(ReadSlabPath(VelSlab,slab).c_str());
        }
        case AuxSlab        : { 
            return fsize(ReadSlabPath(AuxSlab,slab).c_str());
        }
        case AccSlab        : {
            return fsize(ReadSlabPath(PosSlab,slab).c_str())/sizeof(posstruct) * sizeof(accstruct);
        }
        case FarAccSlab        : {
                    return fsize(ReadSlabPath(PosSlab,slab).c_str())/sizeof(posstruct) * sizeof(acc3struct);
                }
        case NearAccSlab        : {
                // TODO: We might be able to remove this
                    return fsize(ReadSlabPath(PosSlab,slab).c_str())/sizeof(posstruct) * sizeof(accstruct);
                }
        case VelLPTSlab        : {
            if(strcmp(P.ICFormat, "RVZel") == 0)
                // TODO: when we re-route LoadIC through SB to be non-blocking for 2LPT velocity IO, get rid of these magic numbers
                //return fsize(ReadSlabPath(VelLPTSlab,slab).c_str())/sizeof(ICfile_RVZel::ICparticle) * sizeof(velstruct);
                return fsize(ReadSlabPath(VelLPTSlab,slab).c_str())/32 * sizeof(velstruct);
            else if(strcmp(P.ICFormat, "RVdoubleZel") == 0)
                //return fsize(ReadSlabPath(VelLPTSlab,slab).c_str())/sizeof(ICfile_RVdoubleZel::ICparticle) * sizeof(velstruct);
                return fsize(ReadSlabPath(VelLPTSlab,slab).c_str())/56 * sizeof(velstruct);
            else
                QUIT("Unknown ICFormat for re-reading LPT IC velocities\n");
        }
        case TimeSlice : {
            return fsize(ReadSlabPath(TimeSlice,slab).c_str());
        }
        
        /* // Ideally this is how we would allocate group finding arenas
                // but there's an annoying circular depdendency: we don't know about GFC yet, but GFC has fields that require slabbuffer.
                // Not easily solved with a forward declaration.
                // TODO: better way to do this?
                case L1halosSlab           : { return GFC->globalslabs[slab]->L1halos.get_slab_bytes(); }
        case TaggedPIDsSlab        : { return GFC->globalslabs[slab]->TaggedPIDs.get_slab_bytes(); }
        case L1ParticlesSlab       : { return GFC->globalslabs[slab]->L1Particles.get_slab_bytes(); }
        case L1PIDsSlab            : { return GFC->globalslabs[slab]->L1PIDs.get_slab_bytes(); }
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
    AllocateSpecificSize(type, slab, s, ramdisk);

    return GetSlabPtr(type, slab);
}

void SlabBuffer::AllocateSpecificSize(int type, int slab, uint64 sizebytes, int ramdisk) {
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

    STDLOG(1, "Allocating slab %d of type %d to size %l (ramdisk = %d), total %5.3f GB\n",
                slab, type, sizebytes, ramdisk, AA->total_allocation/1024./1024./1024.);

    int id = TypeSlab2ID(type,slab);
    AA->Allocate(id, sizebytes, ReuseID(type), ramdisk, ramdisk_fn);
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
        STDLOG(1, "Skipping explicit write of Ramdisk slab \"%s\"", path);

        // still might have to deallocate
        if(deleteafter == IO_DELETE)
            DeAllocate(type, slab);
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
        STDLOG(1, "Skipping explicit read of Ramdisk slab \"%s\"", path);
        SetIOCompleted(type, slab);
        return;
    }

    // TODO: a more flexible way to do this might be to explicitly pass paths down the stack into the ArenaAllocator
    // That way, any arena can get attached to any shared memory path
    // Right now, we're married to allocations only coming from the canonical Descriptor path
    ReadArena(type, slab, blocking, path);
}

void SlabBuffer::LoadArenaBlocking(int type, int slab) {
    AllocateArena(type, slab, RAMDISK_AUTO_READSLAB);
    ReadArenaBlocking(type, slab);
}

void SlabBuffer::LoadArenaNonBlocking(int type, int slab) {
    AllocateArena(type, slab, RAMDISK_AUTO_READSLAB);
    ReadArenaNonBlocking(type, slab);
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

    STDLOG(0,"Writing slab %d of type %d to file %s with blocking %d and delete status %d.\n",
        slab, type, fn, blocking, deleteafter);
    assertf(IsSlabPresent(type, slab), 
        "Type %d and Slab %d doesn't exist\n", type, slab);

    WriteFile( GetSlabPtr(type,slab), 
        SlabSizeBytes(type,slab),
        type, slab,
        fn, 
        0,      // No file offset
        deleteafter, 
        blocking);
}

#endif // INCLUDE_SLABBUFFER
