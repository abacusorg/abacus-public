/*

This class contains implementations of the ICFile abstract base
class from IC_base.h.  This is where we define our IC formats.

New IC formats must implement at least the unpack() method and
a constructor that stores the number of particles in Npart, and
register themselves in ICFile::FromFormat() at the bottom
of this file.  See IC_base.h for the full interface, and look at the
formats below like ICFile_RVZel for examples.

The UnpacktoIL() function at the bottom of this file is the entry point
for timestepIC.

Particles don't really need to be part of the slab indicated by their filename.
They just need to be within FINISH_WAIT_RADIUS~2 slabs.

Several parameters control IC processing. They are:

- ICFormat: will determine the format of the IC file.
- ICPositionRange: The initial positions will be rescaled by ICPositionRange,
   e.g., they might be between 0 and ICPositionRange.
   If ICPositionRange<=0, then it is assumed to be BoxSize.
- ICVelocity2Displacement: The initial velocities will be multiplied by
   ICVelocity2Displacement to get to redshift space displacements
   at the initial redshift in the box units of ICPositionRange in size.
   If ICVelocity2Displacement==-1, then the velocities are to be given in km/s.

*/


// ===================================================================

#include <gsl/gsl_rng.h>
#include "particle_subsample.cpp"

#include "IC_base.h"

/* Each IC format should define a new class that inherits from the ICFile
 * abstract base class.  Each class should also have an entry in the
 * ICFile::FromFormat() factory function at the bottom of this file.
 *
 * Each ICFile class must define two routines: unpack() and a constructor.
 * unpack() pushes the particles to the insert list or (if given) velslab.
 * The constructor must record the number of particles in Npart.
 */

/****************
 * ICFile
 ****************/

// Unpack the particles in the arena and push them to the Insert List
uint64 ICFile::unpack_to_IL(double convert_pos, double convert_vel){
    return unpack(convert_pos, convert_vel);
}

// Start a read via SB
void ICFile::read_nonblocking(){
    SB->LoadArenaNonBlocking(ICSlab, slab);
}

// If a read is in progress, return 0
int ICFile::check_read_done(){
    if(!SB->IsIOCompleted(ICSlab, slab)){
        if(SB->IsSlabPresent(ICSlab, slab))
            SlabDependency::NotifySpinning(WAITING_FOR_IO);
        return 0;
    }
    return 1;
}

ICFile::ICFile(int _slab, int _zsplit){
    slab = _slab;
    zsplit = _zsplit;

    // Subclasses should init this in the constructor
    Npart = 0;
    foffset = 0;
    fbytes = 0;

    lpt_vel_scale = 0.;
}

inline void ICFile::set_taggable_bits(auxstruct &aux, uint64 &sumA, uint64 &sumB){
    // Set the 'taggable' bit for particles as a function of their PID
    // this bit will be used for group finding and merger trees
    int subsample = is_subsample_particle(aux.pid(), P.ParticleSubsampleA, P.ParticleSubsampleB);
    if (subsample == SUBSAMPLE_A) // set all particles 0-ParticleSubsampleA as taggable in set A. 
        { aux.set_taggable_subA(); sumA +=1; }

    else if (subsample == SUBSAMPLE_B) //set all particles A BBig as taggable in set B. 
        { aux.set_taggable_subB(); sumB +=1; }
}

uint64 ICFile::next_pid;  // static next_pid counter


// ICFile_OnDisk<T> is a convenience intermediary class to reduce boilerplate.
// Note that T is actually the type of the child class. So the child's inheritance
// of this class is recursive! This is the "Curiously Recursive Template Pattern" (CRTP).
template<class T>
class ICFile_OnDisk : public ICFile {
public:
    ICFile_OnDisk(int _slab, int _zsplit) : ICFile(_slab, _zsplit){
        // SlabSizeBytes only works if ICSlab is already in memory
        // This avoids an extra file stat
        /*if(SB->IsSlabPresent(ICSlab,slab))
            fbytes = SB->SlabSizeBytes(ICSlab, slab);
        else*/
        fbytes = fs::file_size(SB->ReadSlabPath(ICSlab,slab));

        Npart = fbytes/sizeof(typename T::ICparticle);
        foffset = 0;
        
        assertf(Npart*sizeof(typename T::ICparticle) == fbytes,
            "Size of IC slab {} not divisible by particle size {:d}!\n",
            SB->ReadSlabPath(ICSlab, slab), sizeof(typename T::ICparticle));
    }
};

// ===================================================================
// Simple RVdouble initial conditions: lists of positions & velocities 
// in double precision.  The positions and velocities will be 
// translated by ICPositionRange and ICVelocity2Displacement.

class ICFile_RVdouble: public ICFile_OnDisk<ICFile_RVdouble> {
public:
    class ICparticle {
    public:
        double pos[3];
        double vel[3];
    };

    using ICFile_OnDisk<ICFile_RVdouble>::ICFile_OnDisk;  // inherit the constructor

private:
    // If velslab != NULL, unpack the velocities into there and do nothing else
    uint64 unpack(double convert_pos, double convert_vel) {
        ICparticle *particles = (ICparticle *) SB->GetSlabPtr(ICSlab, slab);

        uint64 sumA = 0, sumB = 0;
        #pragma omp parallel for schedule(static) reduction(+:sumA,sumB)
        for(uint64 i = 0; i < Npart; i++){
            ICparticle p = particles[i];

            velstruct vel(p.vel[0], p.vel[1], p.vel[2]);
            vel *= convert_vel;

            double3 pos(p.pos[0], p.pos[1], p.pos[2]);
            auxstruct aux;

            pos *= convert_pos;

            #ifdef GLOBALPOS
            integer3 newcell = CP->WrapPosition2Cell(&pos);
            #else
            // make pos cell-centered:
            integer3 newcell = CP->LocalPosition2Cell(&pos);
            #endif

            // make pos float3
            posstruct _pos(pos);

            aux.clear();
            aux.packpid(next_pid+i); // use packpid to pack our "linear" pid into the distributed aux format
            set_taggable_bits(aux, sumA, sumB);

            IL->WrapAndPush(&_pos, &vel, &aux, newcell);
        }

        // Store the local reductions in the class variables
        // Curiously, although OpenMP usually supports reduction into class variables,
        // it doesn't like something about the inheritance or templating here.
        NsubsampleA = sumA;
        NsubsampleB = sumB;
        next_pid += Npart;

        // All done with these particles
        SB->DeAllocate(ICSlab, slab);
        
        return Npart;
    }
};

/*
Almost the same binary format as RVdouble, but with a field for the PID.
Positions are global.

This is a template that is typedef'd to RVPID and RVdoublePID below.
*/
template<typename T>
class ICFile_RVPIDTemplate: public ICFile_OnDisk<ICFile_RVPIDTemplate<T>> {
public:
    class ICparticle {
    public:
        T pos[3];
        T vel[3];
        uint64 pid;
    };

    using ICFile_OnDisk<ICFile_RVPIDTemplate<T>>::ICFile_OnDisk;  // inherit the constructor

private:
    // If velslab != NULL, unpack the velocities into there and do nothing else
    uint64 unpack(double convert_pos, double convert_vel) {
        ICparticle *particles = (ICparticle *) SB->GetSlabPtr(ICSlab, this->slab);

        uint64 sumA = 0, sumB = 0;
        #pragma omp parallel for schedule(static) reduction(+:sumA,sumB)
        for(uint64 i = 0; i < this->Npart; i++){
            ICparticle p = particles[i];

            velstruct vel(p.vel[0], p.vel[1], p.vel[2]);
            vel *= convert_vel;

            double3 pos(p.pos[0], p.pos[1], p.pos[2]);
            auxstruct aux;

            pos *= convert_pos;

            #ifdef GLOBALPOS
            integer3 newcell = CP->WrapPosition2Cell(&pos);
            #else
            // make pos cell-centered:
            integer3 newcell = CP->LocalPosition2Cell(&pos);
            #endif

            // make pos float3
            posstruct _pos(pos);

            aux.clear();
            // In this format, we're given a PID that we must then pack into a set of non-contiguous bits
            // in the aux field.  The original PID can still be reconstructed, but beware that naively it
            // will not look the same.
            aux.packpid(p.pid);
            this->set_taggable_bits(aux, sumA, sumB);

            IL->WrapAndPush(&_pos, &vel, &aux, newcell);
        }

        // Store the local reductions in the class variables
        this->NsubsampleA = sumA;
        this->NsubsampleB = sumB;

        // All done with these particles
        SB->DeAllocate(ICSlab, this->slab);
        
        return this->Npart;
    }
};

using ICFile_RVPID = ICFile_RVPIDTemplate<float>;
using ICFile_RVdoublePID = ICFile_RVPIDTemplate<double>;

/*
A template class for RVZel and RVdoubleZel, which are typedef'd with the "using" keyword below
*/
template<typename T>
class ICFile_RVZelTemplate : public ICFile_OnDisk<ICFile_RVZelTemplate<T>> {
private:
    double Canonical_to_VelZSpace;
    int prep_2lpt;

public:
    class ICparticle {
    public:
        unsigned short i,j,k;
        T pos[3];
        T vel[3];
    };

    ICFile_RVZelTemplate(int _slab, int _zsplit)
        : ICFile_OnDisk<ICFile_RVZelTemplate<T>>(_slab, _zsplit) {
        
        // Npart in the parent constructor is the total in the file
        // but this constructor will set it to the per-zrank value.

        uint64 ppd = WriteState.ippd;
        uint64 xplanes = this->Npart / (ppd*ppd);
        assertf(xplanes*ppd*ppd == this->Npart,
            "IC file does not contain a whole number of ppd^2 particle planes?\n");
        
        if(MPI_size_z > 1){
            uint64 zstart = MPI_rank_z*ppd/MPI_size_z;
            uint64 znext = (MPI_rank_z + 1)*ppd/MPI_size_z;
            this->Npart = xplanes*ppd*(znext - zstart);
            this->fbytes = this->Npart*sizeof(ICparticle);
            this->foffset = zstart*ppd*xplanes*sizeof(ICparticle);
        } else {
            this->foffset = 0;
        }

        Canonical_to_VelZSpace = 1./WriteState.VelZSpace_to_Canonical;
        prep_2lpt = P.LagrangianPTOrder > 1;
    }

private:
    // If velslab != NULL, unpack the velocities into there and do nothing else
    uint64 unpack(double convert_pos, double convert_vel) {

        ICparticle *particles = (ICparticle *) SB->GetSlabPtr(ICSlab, this->slab);

        uint64 sumA = 0, sumB = 0;
        FLOAT max_vel_delta = 0.;
        #pragma omp parallel for schedule(static) reduction(+:sumA,sumB) reduction(max:max_vel_delta)
        for(uint64 i = 0; i < this->Npart; i++){
            ICparticle p = particles[i];

            velstruct vel(p.vel[0], p.vel[1], p.vel[2]);
            vel *= convert_vel;

            double3 pos(p.pos[0], p.pos[1], p.pos[2]);
            auxstruct aux;

            pos *= convert_pos;

            if(prep_2lpt){
                // Track max(vel-linearvel)
                // TODO: have we been neglecting f_growth in the zeldovich code?
                velstruct delta = vel*Canonical_to_VelZSpace - pos;
                max_vel_delta = std::max(delta.maxabscomponent(), max_vel_delta);
            }

            integer3 ijk(p.i, p.j, p.k);
            assert(ijk.x >= 0 && ijk.y >= 0 && ijk.z >= 0);
            assert(ijk.x < WriteState.ppd && ijk.y < WriteState.ppd && ijk.z < WriteState.ppd);
            // make pos global:
            pos += ZelPos(ijk);

            #ifdef GLOBALPOS
            integer3 newcell = CP->WrapPosition2Cell(&pos);
            #else
            // make pos cell-centered:
            integer3 newcell = CP->LocalPosition2Cell(&pos);
            #endif

            // make pos float3
            posstruct _pos(pos);

            aux.clear();
            aux.setpid(ijk); // Set aux too.
            this->set_taggable_bits(aux, sumA, sumB);

            IL->WrapAndPush(&_pos, &vel, &aux, newcell);
        }

        // Store the local reductions in the class variables
        this->NsubsampleA = sumA;
        this->NsubsampleB = sumB;

        this->lpt_vel_scale = max_vel_delta;

        // All done with these particles
        SB->DeAllocate(ICSlab, this->slab);
        
        return this->Npart;
    }
};

using ICFile_RVZel = ICFile_RVZelTemplate<float>;
using ICFile_RVdoubleZel = ICFile_RVZelTemplate<double>;


// ===================================================================
// Zel'dovich displacement initial conditions.
// Here we're given the integer-based Lagrangian position and
// the displacement.  The displacement should be in units of 
// P.ICPositionMax, which could plausibly be 1 or BoxSize.
// The velocities will be derived from the displacements.
// The Lagrangian positions will be ijk = 0..CPD-1.  
// We will store i alone, and then jk = j*CPD+k, just to avoid 64-bit math
// and to keep the file 4-byte aligned.

// Our LPT implementation must be coordinated with choices here.
// Therefore, we will codify some items in some simple functions,
// given in lpt.cpp.
// double3 ZelPos(integer3 ijk): Returns the position in code-units 
//             of the initial grid point.


class ICFile_Zeldovich: public ICFile_OnDisk<ICFile_Zeldovich> {
public:
    class ICparticle {
    public:
        unsigned short i,j,k;
        double pos[3];
    };

    using ICFile_OnDisk<ICFile_Zeldovich>::ICFile_OnDisk;  // inherit the constructor

private:
    uint64 unpack(double convert_pos, double convert_vel) {
        ICparticle *particles = (ICparticle *) SB->GetSlabPtr(ICSlab, slab);

        uint64 sumA = 0, sumB = 0;
        #pragma omp parallel for schedule(static) reduction(+:sumA,sumB)
        for(uint64 i = 0; i < Npart; i++){
            ICparticle p = particles[i];

            double3 pos(p.pos[0], p.pos[1], p.pos[2]);
            double3 vel = pos;
            vel *= convert_vel;

            pos *= convert_pos;

            auxstruct aux;

            integer3 ijk(p.i, p.j, p.k);
            assert(ijk.x >= 0 && ijk.y >= 0 && ijk.z >= 0);
            assert(ijk.x < WriteState.ppd && ijk.y < WriteState.ppd && ijk.z < WriteState.ppd);
            // make pos global:
            pos += ZelPos(ijk);

            #ifdef GLOBALPOS
            integer3 newcell = CP->WrapPosition2Cell(&pos);
            #else
            // make pos cell-centered:
            integer3 newcell = CP->LocalPosition2Cell(&pos);
            #endif

            // make float3
            posstruct _pos(pos);
            velstruct _vel(vel);

            aux.clear();
            aux.setpid(ijk); // Set aux too.
            set_taggable_bits(aux, sumA, sumB);

            IL->WrapAndPush(&_pos, &_vel, &aux, newcell);
        }

        // Store the local reductions in the class variables
        NsubsampleA = sumA;
        NsubsampleB = sumB;

        // All done with these particles
        SB->DeAllocate(ICSlab, slab);
        
        return Npart;
    }
};

class ICFile_Lattice: public ICFile {
private:
    int64_t ppdy;  // particles per Y dim

    int firstx, firstz;  // first plane indices in this slab
    int lastx, lastz;
    
public:
    ICFile_Lattice(int _slab, int _zsplit) : ICFile(_slab, _zsplit) {

        ppdy = WriteState.ippd;
        
        firstx = slab*WriteState.ippd/P.cpd;
        lastx = (slab+1)*WriteState.ippd/P.cpd;
        firstz = zsplit*WriteState.ippd/MPI_size_z;
        lastz = (zsplit+1)*WriteState.ippd/MPI_size_z;

        Npart = ppdy*(lastz - firstz)*(lastx - firstx);
    }

    void read_nonblocking(){
        // Lattice is in-memory; nothing to read!
        return;
    }

    int check_read_done(){
        return 1;  // nothing to read!
    }

    uint64 unpack(double convert_pos [[maybe_unused]], double convert_vel) {
        uint64 sumA = 0, sumB = 0;
        double ppd = WriteState.ppd;
        velstruct vel(convert_vel);  // inherit the vel from P.ICVelocity2Displacement
        
        #pragma omp parallel for schedule(static) reduction(+:sumA,sumB)
        for(int64 iy = 0; iy < ppdy; iy++){
            double3 gpos;
            gpos.y = (double)iy/ppd - 0.5;
            for(int64 ix = firstx; ix < lastx; ix++){
                gpos.x = (double)ix/ppd - 0.5;
                for(int64 iz = firstz; iz < lastz; iz++){
                    gpos.z = (double)iz/ppd - 0.5;

                    double3 pos(gpos);
                    
                    #ifdef GLOBALPOS
                    integer3 newcell = CP->WrapPosition2Cell(&pos);
                    #else
                    // make pos cell-centered:
                    integer3 newcell = CP->LocalPosition2Cell(&pos);
                    #endif

                    // make pos float3
                    posstruct _pos(pos);

                    auxstruct aux;
                    aux.clear();
                    aux._setpid(ix, iy, iz); // Set aux too.
                    set_taggable_bits(aux, sumA, sumB);

                    // One might expect to use Push() instead of WrapAndPush() here.
                    // But LocalPosition2Cell() uses doubles, which may spill a cell
                    // once converted to float.
                    IL->WrapAndPush(&_pos, &vel, &aux, newcell);
                }
            }
        }

        // Store the local reductions in the class variables
        NsubsampleA = sumA;
        NsubsampleB = sumB;
        next_pid += Npart;
        
        return Npart;
    }
};

// Generate the position and velocity conversion factors based on the parameter file
// This is called in two places: this module, and lpt.cpp
void get_IC_unit_conversions(double &convert_pos, double &convert_vel){
    if (P.ICPositionRange>0)
        convert_pos = 1.0/P.ICPositionRange;
    else
        convert_pos = 1.0/P.BoxSize;
    if (P.FlipZelDisp)
        convert_pos *= -1;
    if (P.ICVelocity2Displacement!=-1) // Should always be 1 for IC from zeldovich
        convert_vel = P.ICVelocity2Displacement*convert_pos;
    else
        convert_vel = 1.0/ReadState.VelZSpace_to_kms;
    // This gets to the unit box in redshift-space displacement.

    // The getparticle() routine should return velocities in 
    // redshift-space comoving displacements in length units where
    // the box is unit length.  
    // We need to convert to canonical velocities.
    // Code time unit is 1/H_0, but H(z) gets output with H0=1.
    // The conversion is v_code = v_zspace * a^2 H(z)/H_0
    convert_vel *= WriteState.VelZSpace_to_Canonical;
}


// An alternative to this macro approach would be a C++ type registry
#define REGISTER_ICFORMAT(fmt) if(format == #fmt){\
    ic = std::unique_ptr<ICFile_##fmt>(new ICFile_##fmt(_slab, MPI_rank_z));\
} else

// This is a factory function to instantiate ICFile objects of a given format and slab number
// It returns a C++ unique_ptr; no need to call delete on it later.  The object will be deleted
// when it leaves scope.
std::unique_ptr<ICFile> ICFile::FromFormat(const std::string& format, int _slab){

    std::unique_ptr<ICFile> ic;

    // Do we ever need Poisson instead of Lattice?
    //REGISTER_ICFORMAT(Poisson)

    // This is actually a big if-else chain; no semicolons!
    REGISTER_ICFORMAT(RVdouble)
    REGISTER_ICFORMAT(RVdoubleZel)
    REGISTER_ICFORMAT(RVZel)
    REGISTER_ICFORMAT(RVPID)
    REGISTER_ICFORMAT(RVdoublePID)
    REGISTER_ICFORMAT(Zeldovich)
    REGISTER_ICFORMAT(Lattice)
    {
        // We weren't given a legal format name.
        QUIT("Unrecognized case: ICFormat = {:s}\n", format);
    }

    return ic;
}


uint64 UnpackICtoIL(int slab) {

    // Set up unit conversions
    double convert_pos, convert_vel;
    get_IC_unit_conversions(convert_pos, convert_vel);

    STDLOG(1,"Using IC format {:s}\n", P.ICFormat);
    std::unique_ptr<ICFile> ic = ICFile::FromFormat(P.ICFormat, slab);

    // Unpacks the whole slab directly to the insert list
    uint64 count = ic->unpack_to_IL(convert_pos, convert_vel);

    STDLOG(0,"Read {:d} particles from IC slab {:d}\n", count, slab);
    STDLOG(1,"Slab {:d} has {:d} subsample A particles, {:d} subsample B particles.\n", slab, ic->NsubsampleA, ic->NsubsampleB);

    WriteState.LPTVelScale = std::max((double)std::abs(ic->lpt_vel_scale), WriteState.LPTVelScale);

    return count; 
}
