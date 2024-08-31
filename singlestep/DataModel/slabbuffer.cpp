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
                CellGroupArena,         //17
                NearField_SIC_Slab,     //18

                L1halosSlab,            //19
                HaloRVSlabA,            //20
                HaloRVSlabB,            //21
                HaloPIDsSlabA,          //22
                HaloPIDsSlabB,          //23
                FieldRVSlabA,           //24
				FieldRVSlabB,           //25
                FieldPIDSlabA,          //26
                FieldPIDSlabB,          //27
                L0TimeSlice,            //28
                L0TimeSlicePIDs,        //29

                ICSlab,                 //39

                // The start of the LC slabs
                // SB->NumTypes is the total number of types
                LightCone0RV,
                LightCone0PID,
                LightCone0Heal
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

    uint64 ReadOffset(int type, int slab);

public:

    int NumTypes;

    // The write and read names for the files of this SlabType.
    // The slab number will be wrapped.
    // Program will abort if one specifies a type not in this function;
    // this is a feature to block incorrect I/O.
    fs::path WriteSlabPath(int type, int slab);
    fs::path ReadSlabPath(int type, int slab);

    SlabBuffer(int _cpd, int _order, int num_lc) {
        NumTypes = LightCone0RV + 3*num_lc;

        order = _order;
        cpd = _cpd;

        AA = new ArenaAllocator((_cpd+1)*NumTypes, P.UseMunmapThread, P.MunmapThreadCore);
    }

    ~SlabBuffer(void) {
        // Clean up the Reuse buffers
        for (int t = 0; t < NumTypes; t++) {
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
    void ReadArena(int type, int slab, int blocking, const fs::path &fn);
    void WriteArena(int type, int slab, int deleteafter, int blocking, const fs::path &fn);

    void DeAllocate(int type, int slab, int delete_file=0);

    inline char* GetSlabPtr(int type, int slab) {
        int id = TypeSlab2ID(type,slab);
        assertf(AA->IsArenaPresent(id), "Error: slab {:d} of type {:d} was not present, when pointer requested.\n", slab, type);
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
        STDLOG(2,"Resizing slab {:d} of type {:d} to size {:d}\n", slab, type, size);
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

    void report_peak(int force=0){
        AA->report_peak(force);
    }

    void report_current(){
        AA->report_current();
    }
};


/** Decide whether the given slab lives on the ramdisk
  * If so, we want to inform ArenaAllocator of this with the RAMDISK_READSLAB or RAMDISK_WRITESLAB flags,
  * which will be returned by this function as appropriate.
  * The hint allows one to override normal assumptions about a slab being a "read type" or "write type".
  * Valid values of `hint` are the RAMDISK_* flags.
  */
int SlabBuffer::IsRamdiskSlab(int type, int hint){

    //return 0;

    // It's still possible to write a slab to ramdisk even if it's not listed below
    // The IO module should detect that it's a ramdisk path and use fopen instead of DIO

    // Does the hint indicate that we already decided the fate of this slab?
    if(hint >= RAMDISK_NO && hint < RAMDISK_AUTO)
        return hint;

    // If not, compare the file path for the slab to the system ramdisk path

    int intent = GetSlabIntent(type);

    switch(hint){
        case RAMDISK_AUTO_READSLAB:
            intent = READSLAB;
            break;
        case RAMDISK_AUTO_WRITESLAB:
            intent = WRITESLAB;
            break;
        case RAMDISK_AUTO:
            break;
        default:
            QUIT("Did not understand ramdisk hint {:d}\n", hint);
    }

    switch(intent){
        case READSLAB:
            if(is_path_on_ramdisk(ReadSlabPath(type, -1)))
                return RAMDISK_READSLAB;
            break;

        case WRITESLAB:
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

/* Determine whether a given slab type (e.g. TimeSliceSlab) is an
 * output type that we do not want to send in the manifest.
*/
int SlabBuffer::IsOutputSlab(int type){
    if (type >= LightCone0RV && type < NumTypes){
        return 1;
    }

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
            return 1;
     }
     return 0;
}


/* Determine whether a given slab type (e.g. TimeSliceSlab) is an type
 * that we want to record the checksum of when we write it to disk.
 *
 * P.NoChecksum disables all checksumming.
*/
int SlabBuffer::WantChecksum(int type){
    if(P.NoChecksum)
        return 0;
    return IsOutputSlab(type);
}


fs::path SlabBuffer::WriteSlabPath(int type, int slab) {
    slab = Grid->WrapSlab(slab);

    std::string slabstr;
    if(P.NumZRanks > 1)
        slabstr = fmt::format("{:04d}.{:03d}", slab, MPI_rank_z);
    else
        slabstr = fmt::format("{:04d}", slab);
    std::string stepstr = fmt::format("{:04d}", ReadState.FullStepNumber);
    std::string zstr = fmt::format("{:5.3f}", ReadState.Redshift);

    if (type >= LightCone0RV && type < NumTypes) {
        int lcn = (type - LightCone0RV) / 3;
        int lct = (type - LightCone0RV) % 3;
        std::string lct_name = (lct == 0) ? "rv" : (lct == 1) ? "pid" : "heal";
        return P.LightConeDirectory / ("Step" + stepstr) / fmt::format("LightCone{:d}_{:s}{:s}", lcn, lct_name, NodeString);
    }

    switch(type) {
#ifdef PARALLEL
        case TaylorSlab     : {
            // Read odd taylors from TaylorDirectory2, if it's defined
            if (WriteState.StripeConvState && slab % 2 == 1) {
                return P.TaylorDirectory2 / ("Taylor_" + slabstr);
            }
            return P.TaylorDirectory / ("Taylor_"     + slabstr); break;
        }
#endif
        case MultipoleSlab       : {


            // Send odd multipoles to MultipoleDirectory2, if it's defined
            if (WriteState.StripeConvState && slab % 2 == 1) {
                return P.MultipoleDirectory2 / ("Multipoles_" + slabstr);
                break;
            }
            return P.MultipoleDirectory / ("Multipoles_" + slabstr); break;

        }
        case MergeCellInfoSlab   : { return P.LocalWriteStateDirectory / ("cellinfo_"   + slabstr); break; }
        case MergePosSlab        : { return P.LocalWriteStateDirectory / ("position_"   + slabstr); break; }
        case MergeVelSlab        : { return P.LocalWriteStateDirectory / ("velocity_"   + slabstr); break; }
        case MergeAuxSlab        : { return P.LocalWriteStateDirectory / ("auxillary_"  + slabstr); break; }
        case AccSlab             : { return P.OutputDirectory / "acc" / ("Step" + stepstr) / ("acc_" + slabstr); break; }
        case NearAccSlab         : { return P.OutputDirectory / "acc" / ("Step" + stepstr) / ("nearacc_" + slabstr); break; }
        case FarAccSlab          : { return P.OutputDirectory / "acc" / ("Step" + stepstr) / ("faracc_" + slabstr); break; }

        case L1halosSlab           : { return P.GroupDirectory / ("Step" + stepstr) / ("_z" + zstr) / ("halo_info_" + slabstr); break;}

        case HaloPIDsSlabA    : { return P.GroupDirectory / ("Step" + stepstr) / ("_z" + zstr) / ("halo_pids_A_" + slabstr); break;}
        case HaloPIDsSlabB    : { return P.GroupDirectory / ("Step" + stepstr) / ("_z" + zstr) / ("halo_pids_B_" + slabstr); break;}
        case FieldPIDSlabA    : { return P.GroupDirectory / ("Step" + stepstr) / ("_z" + zstr) / ("field_pids_A_" + slabstr); break;}
        case FieldPIDSlabB    : { return P.GroupDirectory / ("Step" + stepstr) / ("_z" + zstr) / ("field_pids_B_" + slabstr); break;}

        case HaloRVSlabA : { return P.GroupDirectory / ("Step" + stepstr) / ("_z" + zstr) / ("halo_rv_A_" + slabstr); break;}
        case HaloRVSlabB : { return P.GroupDirectory / ("Step" + stepstr) / ("_z" + zstr) / ("halo_rv_B_" + slabstr); break;}
        case FieldRVSlabA       : { return P.GroupDirectory / ("Step" + stepstr) / ("_z" + zstr) / ("field_rv_A_" + slabstr); break;}
        case FieldRVSlabB       : { return P.GroupDirectory / ("Step" + stepstr) / ("_z" + zstr) / ("field_rv_B_"+ slabstr); break;}

        case L0TimeSlice          : { return P.OutputDirectory / ("slice" + zstr) / (P.SimName + ".z" + zstr + ".slab" + slabstr + ".L0_pack9.dat"); break; }
        case FieldTimeSlice            : { return P.OutputDirectory / ("slice" + zstr) / (P.SimName + ".z" + zstr + ".slab" + slabstr + ".field_pack9.dat"); break; }
        case L0TimeSlicePIDs      : { return P.OutputDirectory / ("slice" + zstr) / (P.SimName + ".z" + zstr + ".slab" + slabstr + ".L0_pack9_pids.dat"); break; }
        case FieldTimeSlicePIDs        : { return P.OutputDirectory / ("slice" + zstr) / (P.SimName + ".z" + zstr + ".slab" + slabstr + ".field_pack9_pids.dat"); break; }

        default:
            QUIT("Illegal type {:d} given to WriteSlabPath()\n", type);
    }

    return "";
}

fs::path SlabBuffer::ReadSlabPath(int type, int slab) {
    slab = Grid->WrapSlab(slab);

    std::string slabstr;
    if(P.NumZRanks > 1)
        slabstr = fmt::format("{:04d}.{:03d}", slab, MPI_rank_z);
    else
        slabstr = fmt::format("{:04d}", slab);

    switch(type) {
        case TaylorSlab     : {
            // Read odd taylors from TaylorDirectory2, if it's defined
            if (WriteState.StripeConvState && slab % 2 == 1) {
                return P.TaylorDirectory2 / ("Taylor_" + slabstr);
                break;
            }
            return P.TaylorDirectory / ("Taylor_" + slabstr); break;
        }

        // The Merge slabs are usually write slabs, but can be used as read slabs in multipole recovery mode
        case MergeCellInfoSlab :
        case CellInfoSlab  : { return P.LocalReadStateDirectory / ("cellinfo_"   + slabstr); break; }

        case MergePosSlab :
        case PosSlab       : { return P.LocalReadStateDirectory / ("position_"   + slabstr); break; }

        case VelSlab       : { return P.LocalReadStateDirectory / ("velocity_"   + slabstr); break; }
        case AuxSlab       : { return P.LocalReadStateDirectory / ("auxillary_"  + slabstr); break; }
        
        // used for standalone FOF
        case FieldTimeSlice     : { return WriteSlabPath(type, slab); break; }
        case L0TimeSlice     : { return WriteSlabPath(type, slab); break; }
        case FieldTimeSlicePIDs     : { return WriteSlabPath(type, slab); break; }
        case L0TimeSlicePIDs     : { return WriteSlabPath(type, slab); break; }
        
        case ICSlab    : {
            if (P.NumZRanks > 1){
                return P.InitialConditionsDirectory / "2D" / fmt::format("ic2D_{:d}", slab);
            } else {
                return P.InitialConditionsDirectory / fmt::format("ic_{:d}", slab);
            }
            break;
        }

        default:
            QUIT("Illegal type {:d} given to ReadSlabPath()\n", type);
    }

    return "";
}

uint64 SlabBuffer::ArenaSize(int type, int slab) {
    int rml = (order+1)*(order+1);
    uint64 lcpd = cpd;  // Just to assure that we don't spill 32-bits

    switch(type) {
        case CellInfoSlab        : { return sizeof(cellinfo)*lcpd*node_z_size_with_ghost; }
        case MergeCellInfoSlab   : { return sizeof(cellinfo)*lcpd*(node_z_size + 2*MERGE_GHOST_RADIUS); }
        case InsertCellInfoSlab  : { return sizeof(cellinfo)*lcpd*(node_z_size + 2*MERGE_GHOST_RADIUS); }
        case MultipoleSlab  : { return lcpd*node_cpdp1half*rml*sizeof(MTCOMPLEX); }
        case TaylorSlab     : { return lcpd*node_cpdp1half*rml*sizeof(MTCOMPLEX); }
        case PosXYZSlab     :  // let these fall through to PosSlab
        case PosSlab        : {
            return SS->size_with_ghost(slab)*sizeof(posstruct);
        }
        case MergePosSlab   : {
            return SS->size_with_ghost(slab)*sizeof(posstruct);
        }
        case VelSlab        : {
            return SS->size_with_ghost(slab)*sizeof(velstruct);
        }
        case AuxSlab        : {
            return SS->size_with_ghost(slab)*sizeof(auxstruct);
        }
        case AccSlab        : {
            // no ghost for acc
            return SS->size(slab)*sizeof(accstruct);
        }
        case FarAccSlab     : {
            return SS->size(slab)*sizeof(acc3struct);
        }
        case NearAccSlab    : {
            return SS->size(slab)*sizeof(accstruct);
        }
        case ICSlab     : {
            return ICFile::FromFormat(P.ICFormat, slab)->fbytes;
        }
        case FieldTimeSlice : {
            return fs::file_size(ReadSlabPath(FieldTimeSlice,slab));
        }
        case L0TimeSlice : {
            return fs::file_size(ReadSlabPath(L0TimeSlice,slab));
        }
        case L0TimeSlicePIDs : {
            return fs::file_size(ReadSlabPath(L0TimeSlicePIDs,slab));
        }
        case FieldTimeSlicePIDs : {
            return fs::file_size(ReadSlabPath(FieldTimeSlicePIDs,slab));
        }
        default            : QUIT("Illegal type {:d} given to ArenaSize()\n", type);
    }
    return -1; //should be unreachable
}

uint64 SlabBuffer::ReadOffset(int type, int slab) {
    // Byte offset for file IO.
    // Only used for the ICs in the 2D code.

    switch(type) {
        case ICSlab: {
            if(MPI_size_z > 1){
                return ICFile::FromFormat(P.ICFormat, slab)->foffset;
            }
            break;
        }
        default: break;
    }
    return 0;
}

char *SlabBuffer::AllocateArena(int type, int slab, int ramdisk) {
    slab = Grid->WrapSlab(slab);
    uint64 s = ArenaSize(type, slab);
    return AllocateSpecificSize(type, slab, s, ramdisk);
}

char *SlabBuffer::AllocateSpecificSize(int type, int slab, uint64 sizebytes, int ramdisk) {
    // Most slabs are happy with RAMDISK_AUTO
    ramdisk = IsRamdiskSlab(type, ramdisk);

    fs::path spath;
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
            QUIT("Unexpected value {:d} of ramdisk\n", ramdisk);
            break;
    }
	STDLOG(3, "Allocating slab {:d} of type {:d} to size {:d} (ramdisk = {:d}), total {:5.3f} GB\n",
                slab, type, sizebytes, ramdisk, AA->total_allocation/1024./1024./1024.);

    int id = TypeSlab2ID(type,slab);
    AA->Allocate(id, sizebytes, ReuseID(type), ramdisk, spath);

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
    fs::path spath = WriteSlabPath(type,slab);
    // Determine the actual, allocated Ramdisk type of the current slab
    // If it was allocated on Ramdisk, then the writing is already done by definition!
    if(AA->ArenaRamdiskType(TypeSlab2ID(type,slab)) != RAMDISK_NO){
        STDLOG(2, "Skipping explicit write of type {:d} Ramdisk slab \"{}\"\n", type, spath);

        // still might have to deallocate
        if(deleteafter == IO_DELETE){
#ifdef PARALLEL
			if (type == MultipoleSlab)
                STDLOG(1, "Uh oh. DeAllocating Multipole slab {:d} prematurely\n", slab);
#endif
            DeAllocate(type, slab);
		}

        return;
    }

    WriteArena(type, slab, deleteafter, blocking, spath);
}


void SlabBuffer::ReadArenaBlocking(int type, int slab) {
    ReadArena(type, slab, IO_BLOCKING);
}

void SlabBuffer::ReadArenaNonBlocking(int type, int slab) {
    ReadArena(type, slab, IO_NONBLOCKING);
}

void SlabBuffer::ReadArena(int type, int slab, int blocking){
    fs::path spath = ReadSlabPath(type,slab);

    // Determine the actual, allocated Ramdisk type of the current slab
    // If it was allocated from existing shared memory (i.e. RAMDISK_READSLAB),
    // we probably don't want to read into it
    if(AA->ArenaRamdiskType(TypeSlab2ID(type,slab)) == RAMDISK_READSLAB){
        STDLOG(2, "Skipping explicit read of Ramdisk slab {:d} \"{}\"\n", slab, spath);
        SetIOCompleted(type, slab);
        return;
    }

    ReadArena(type, slab, blocking, spath);
}

void SlabBuffer::LoadArenaBlocking(int type, int slab) {
    AllocateArena(type, slab, RAMDISK_AUTO_READSLAB);
    ReadArenaBlocking(type, slab);

    if (type != ICSlab)
        assertf(SlabSizeBytes(type, slab) == fs::file_size(ReadSlabPath(type,slab)),
            "Unexpected file size for slab {:d} of type {:d}\n", slab, type);
}

void SlabBuffer::LoadArenaNonBlocking(int type, int slab) {
    AllocateArena(type, slab, RAMDISK_AUTO_READSLAB);
    ReadArenaNonBlocking(type, slab);

    if (type != ICSlab)
        assertf(SlabSizeBytes(type, slab) == fs::file_size(ReadSlabPath(type,slab)),
            "Unexpected file size for slab {:d} of type {:d}\n", slab, type);
}


void SlabBuffer::ReadArena(int type, int slab, int blocking, const fs::path &fn) {
    // This will read into an arena.
    // This always triggers a real read, even if the path is on the ramdisk!
    // The read is always ordered to be the full usable size of the arena.

    if (P.ForceBlockingIO!=0)
        blocking = IO_BLOCKING;

    STDLOG(1,"Reading slab {:d} of type {:d} from file {} with blocking {:d}.\n", slab, type, fn, blocking);
    assertf(IsSlabPresent(type, slab),
        "Type {:d} and Slab {:d} doesn't exist\n", type, slab);
    assertf(fs::exists(fn),
        "File {} does not exist\n", fn);
    assertf((uint64) fs::file_size(fn) >= SlabSizeBytes(type,slab),
        "File {} appears to be too small ({:d} bytes) compared to the arena ({:d})\n",
            fn, fs::file_size(fn), SlabSizeBytes(type,slab));

    uint64 offset = ReadOffset(type, slab);

    ReadFile( GetSlabPtr(type,slab),
        SlabSizeBytes(type,slab),
        type, slab,
        fn,
        offset,
        blocking);
}

void SlabBuffer::WriteArena(int type, int slab, int deleteafter, int blocking, const fs::path &fn) {
    // This will write into an arena.
    // This always triggers a real write, even if the path is on the ramdisk!

    if (P.ForceBlockingIO!=0)
        blocking = IO_BLOCKING;

    STDLOG(1,"Writing slab {:d} of type {:d} to file {} with blocking {:d} and delete status {:d}.\n",
        slab, type, fn, blocking, deleteafter);
    assertf(IsSlabPresent(type, slab),
        "Type {:d} and Slab {:d} doesn't exist\n", type, slab);

    int use_fp = type >= LightCone0RV && type < NumTypes;
    
    WriteFile( GetSlabPtr(type,slab),
        SlabSizeBytes(type,slab),
        type, slab,
        fn,
        0,      // No file offset
        deleteafter,
        blocking,
        WantChecksum(type),
        use_fp);
}

// Free the memory associated with this slab.  Might stash the arena as a reuse slab!
// One can request deletion of the file that corresponds to this slab, too.
void SlabBuffer::DeAllocate(int type, int slab, int delete_file) {
    STDLOG(2,"Deallocating slab {:d} of type {:d}.\n", slab, type);
    AA->DeAllocateArena(TypeSlab2ID(type,slab), ReuseID(type));

    if(delete_file){
        int intent = GetSlabIntent(type);
        fs::path path;
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
                QUIT("Did not understand slab intent {:d}\n", intent);
        }

        if(!IsTrueLocalDirectory(fs::absolute(path).parent_path())){
            STDLOG(2,"Not deleting slab file \"{}\" because it is in a global directory\n", path);
            return;
        }

        STDLOG(2,"Deleting slab file \"{}\"\n", path);
        fs::remove(path);
    }
}

#endif // INCLUDE_SLABBUFFER
