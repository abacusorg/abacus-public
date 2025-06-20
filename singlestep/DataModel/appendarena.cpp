// Copyright 2012-2025 The Abacus Developers
// SPDX-License-Identifier: GPL-3.0-or-later

/* appendarena.cpp
 *
 * \brief Handles outputs in a variety of formats
 *
 * AppendArena is an abstraction layer on top of any number of desired formats for outputting
 * snapshots of particles. All details of how these outputs are done is completely hidden from
 * external code
 *
 * FIXME: There is a significant amount of commented out or non-compiled code in this file to be
 * cleaned up.
 */

// TODO: Ideally we'd turn this into some kind of run-time choice
// between formats, with the PACKED class supplying the pack()
// and pack_cell() methods.

#include "threevector.hh"
#include <math.h>

#ifdef COMPILE
// The following lines are just to check the compiling
// g++ -c -I../include -DCOMPILE appendarena.cpp
#define posstruct double3
#define velstruct double3
#define auxstruct double3
#define SlabType int
#define assertf(a,b) do {} while(0)
#endif

// We need the cell_header class
#include "cell_header.h"
#include "packN_storage.cpp"

struct ArenaPencil {
    char *start;   // Pointer to the first byte of data for this pencil
    char *next;    // Pointer to where we'd write next
    cell_header current_cell;
    char empty[CACHE_LINE_SIZE-2*sizeof(char *)-sizeof(cell_header)];   // Avoid cache line contention
};

class AppendArena {
  private:
    // This is an Abstract Base Class.
    // The following 4 routines must be provided in the derived class.
    // Of course, the derived class must also provide a constructor and destructor.

    //FIXME: Since appending particles is a quick function we are calling many times,
    //  we may want to make the polymorphism compile-time (using templates) rather than
    //  looking up virtual functions at runtime

    virtual void appendparticle(char *c, cell_header current_cell, posstruct pos, velstruct vel, auxstruct aux) =0;
    // Place one particle at the location c
    virtual void appendcell(char *c, cell_header &current_cell, integer3 ijk, float vscale) =0;
    // Place one cell at the location c.  
    // This must also set the cell header current_cell.

    // Note that it is legal for these routines to write no data for cells,
    // but we still need to set the current_cell, since that might be needed
    // to interpret the cell-centered positions.
    // Also note that it is legal to ignore vscale in writing the velocities.
    // But vscale must be >0 to be allowed to add particles.

  public:
    virtual int sizeof_particle() =0;
    // Return the number of bytes to be written for a particle
    virtual int sizeof_cell() =0;
    // Return the number of bytes to be written for a cell

    virtual ~AppendArena(void) {};		// Virtual destructor

  public:
    // The following routines are provided as general services of the ABC.
    // Note that when we write into the Arena, we do not append \0 -- these are not strings!

    // The AppendArena class does know about the Arenas, but we opt to
    // leave the arena Allocation and Write commands
    // in the domain of the calling program.

    // We want to injest velocities and positions in the code units.
    // The base AppendArena addparticle() method then converts the velocities
    // to redshift-space displacements, in the code position units (i.e., unit box).
    // The derived classes are responsible for conversions to global positions and velocities
    // because they might differ between output formats.

    // The derived class constructor might want to use ReadState to
    // figure out the conversion factor from code units to output units.

    // We need to be able to handle multiple threads writing into this.
    // We organize this by pencil; a pencil must be writen consecutively by a single thread.

    char *startarena;  // The first byte of data, including the header
    char *startdata;  // The first byte of data, after the header
    char *endarena;  // The first byte beyond the allowed arena

    ArenaPencil *pencil;  // Pointers for the start and end of each pencil's data
    // The size needed for the pencil is just the difference of these two pointers

    long long int bytes;  // The total number of bytes used
    int cpd;		// CPD to be used in mapping particle positions
    FLOAT VelCanonical_to_ZSpace; 	// Typically 1/ReadState.VelZSpace_to_Canonical

    void initialize(SlabType type, int slab, int _cpd, double VelZSpace_to_Canonical) {
        STDLOG(4, "Initializing...\n");
        startarena = startdata = SB->GetSlabPtr(type,slab);
        long long int size_available = SB->SlabSizeBytes(type,slab);
        endarena = startarena+size_available;
        bytes = 0;
        cpd = _cpd;
        VelCanonical_to_ZSpace = 1.0/VelZSpace_to_Canonical;
        STDLOG(4, "Done Initializing...\n");

        int ret;
        ret = posix_memalign((void **)&pencil, CACHE_LINE_SIZE, sizeof(ArenaPencil)*cpd); assert(ret==0);
        for (int j=0; j<cpd; j++) {
            pencil[j].start = pencil[j].next = NULL;
            endcell(j);
        }
    }

    // ~AppendArena(void) { }

    long long int bytes_written() {
        return bytes;
    }

    void start_pencil(int j, long long int offset) {
        // We're ready to start a pencil with this thread, at an offset location given in bytes.
        // Note: we cannot check that this doesn't overlap another pencils portion of the buffer!
        // So this passes a substantial responsibility to the calling function.
        ArenaPencil *p = pencil+j;
        p->next = p->start = startdata+offset;
    	p->current_cell.vscale = 0; 
    }

    inline void addcell(int j, integer3 ijk, float vscale) {
        // Append this cell and make this cell active.
        int size = sizeof_cell();
        ArenaPencil *p = pencil+j;
        assertf(p->next>=startdata && p->next<endarena, "AppendArena() has run out of space!\n");
        appendcell(p->next, p->current_cell, ijk, vscale);
        p->next+=size; 
    }
    inline void endcell(int j) { 
        // This makes it illegal to add more particles
        ArenaPencil *p = pencil+j;
    	p->current_cell.vscale = 0; 
    }

    inline void addparticle(int j, posstruct pos, velstruct vel, auxstruct aux) {
        // Append this particle
        ArenaPencil *p = pencil+j;
        assertf(p->current_cell.islegal(), "addparticle() called without an active cell\n");
        int size = sizeof_particle();
        assertf(p->next>=startdata && p->next<endarena, "AppendArena() has run out of space!\n");
        // We were given the particle in code units, need to convert to Zspace.
        vel *= VelCanonical_to_ZSpace;
        // The appendparticle() function is responsible for any further conversions,
        // e.g., from local to global positions or to km/s or from unit box to unit cell..
        appendparticle(p->next, p->current_cell, pos, vel, aux);
        p->next+=size; 
    }

    inline void addheader(const std::string &c) {
        // Write the given string to the arena, without the final \0.
        int size = c.size();
        memcpy(startdata, c.c_str(), size);
        startdata += size; bytes += size;
    }

    void finalize_header() {
        // terminate the ParseHeader header.
        // Pad out the header to a multiple of 4096 bytes just to keep the data DIO aligned.
        long long int size = bytes_written()+2;  // including the header
        int pad = 4096-size%4096; 
        if (pad>0) { startdata[0] = '\n'; }
        for (int n=1; n<pad; n++) { startdata[n] = ' '; }
        startdata += pad; bytes += pad;
        startdata[0] = ''; startdata[1] = '\n'; startdata+=2; bytes+=2;
    }

    long long int finalize_arena() {
        // Copy the pencils inside the buffer so that they are contiguous.
        // This must be called, as it cleans up the allocations internal to the class.
        // Return the number of bytes in the arena (including header), since that is often needed.
        char *next = startdata;
        for (int j=0; j<cpd; j++) {
            if (next==pencil[j].start) {
                // Buffer is already where it should be
                next = pencil[j].next;
                continue;   
            }
            if (pencil[j].next == pencil[j].start) continue;  // Pencil was empty
            // Otherwise, we need to copy down 
            long long int size = pencil[j].next-pencil[j].start;
            assertf(size>0, "ArenaPencils appear to be in non-increasing order!  Pencil {:d}\n", j);
            if (pencil[j].start<next+size) {
                // The two regions overlap, so we need to be careful
                memmove(next, pencil[j].start, size);
                // This doesn't happen much, so we just eat the single-threaded case
            } else {
                // The two regions don't overlap, so let's multithread
                #pragma omp parallel for schedule(static)
                for (int i=0; i<size; i++) next[i] = pencil[j].start[i];
            }
            next += size;  // Ready for the next pencil
        }
        bytes += next-startdata;  // This is now the total size of the arena
        free(pencil);
        return bytes_written();
    }
};

// =================================================================
// Here's an example of a derived class that would use the packed functions.
//
// FIXME: The different output formats should probably be in their own files

template <int N>
class OutputPacked: public AppendArena {
  private:
    inline void appendparticle(char *c, cell_header current_cell, posstruct pos, velstruct vel, auxstruct aux) {
        packN<N> *p = (packN<N> *) c;
        // p->pack_global(pos, vel, aux.pid(), current_cell);
        p->pack(pos, vel, aux.pid(), current_cell);
    }
    inline void appendcell(char *c, cell_header &current_cell, integer3 ijk, float vscale) {
        // We're given vscale in Zspace unit-box units, same as velocities.
        // But we need to hand it to the pack14 method in unit-cell units
        packN<N> *p = (packN<N> *) c;
        int vs = ceil(vscale*cpd);
        if (vs<=0) vs = 10;    // Probably just not initialized correctly
        current_cell = p->pack_cell(ijk, cpd, vs);
    }
    float velocity_conversion;

  public:
    int sizeof_cell()     { return sizeof(packN<N>); }
    int sizeof_particle() { return sizeof(packN<N>); }

    OutputPacked() {
        // Use ReadState to figure out the correct conversion of the
        // velocity!  The pack14 class assumes velocites are in
        // redshift-space units; it *then* applies the given vscale.
        velocity_conversion = 1.0;
    }
    ~OutputPacked(void) { }
};


// And here's an example to put things into simple RVdouble

class OutputRVdouble: public AppendArena {
  private:
    struct ICparticle {
        double pos[3];		// Global position, unit box
        double vel[3];		// Zspace, unit box
    };

    float velocity_conversion;
    void appendparticle(char *c, cell_header current_cell, posstruct pos, velstruct vel, auxstruct aux [[maybe_unused]]) {
        struct ICparticle *p = (struct ICparticle *)c;
#ifdef GLOBALPOS
        // This is for global box-referenced positions
        double3 cc = double3(0.0);
#else
        // If we instead have cell-referenced positions, then:
        double3 cc = CP->WrapCellCenter(current_cell.cellid());
#endif
        p->pos[0] = (double) pos.x+cc.x;
        p->pos[1] = (double) pos.y+cc.y;
        p->pos[2] = (double) pos.z+cc.z;
        p->vel[0] = vel.x*velocity_conversion;
        p->vel[1] = vel.y*velocity_conversion;
        p->vel[2] = vel.z*velocity_conversion;
    }
    void appendcell(char *c [[maybe_unused]], cell_header &current_cell, integer3 ijk, float vscale [[maybe_unused]]) {
        current_cell = cell_header(ijk, cpd, 1);
    }

  public:
    int sizeof_cell()     { return 0; }
    int sizeof_particle() { return sizeof(struct ICparticle); }

    OutputRVdouble() {
        velocity_conversion = 1.0;  // Leave it in Zspace, unit box units
    }
    ~OutputRVdouble(void) { }
};


class OutputPID: public AppendArena {
  private:
    //uint64_t pid;

    void appendparticle(char *c, cell_header current_cell [[maybe_unused]], posstruct pos [[maybe_unused]], velstruct vel [[maybe_unused]], auxstruct a) {
        uint64_t *p = (uint64_t *)c;
        *p  = a.get_aux_pid_dens_tagged(); 
    }

    void appendcell(char *c [[maybe_unused]], cell_header &current_cell, integer3 ijk, float vscale [[maybe_unused]]) {
        current_cell = cell_header(ijk, cpd, 1);
    }

  public:
    int sizeof_cell()     { return 0; }
   // int sizeof_particle() { return sizeof(pid); }
    int sizeof_particle() { return sizeof(uint64_t); }

    OutputPID() { }
    ~OutputPID(void) { }
};


// =================================================================
// And here's an example to put things back into the Heitmann format.

class OutputHeitmann: public AppendArena {
  private:
    class ICparticle {
        public:
        float xv[6];	// Global position, unit box; Zspace velocity, unit box
        float mass;
        unsigned int tag;
    };
    float velocity_conversion;
    
    void appendparticle(char *c, cell_header current_cell, posstruct pos, velstruct vel, auxstruct aux) {
        ICparticle *p = (ICparticle *)c;
#ifdef GLOBALPOS
        // This is for global box-referenced positions
        double3 cc = double3(0.0);
#else
        // If we instead have cell-referenced positions, then:
        double3 cc = CP->WrapCellCenter(current_cell.cellid());
#endif
        p->xv[0] = (double) pos.x+cc.x;
        p->xv[2] = (double) pos.y+cc.y;
        p->xv[4] = (double) pos.z+cc.z;
        p->xv[1] = vel[0]*velocity_conversion;
        p->xv[3] = vel[1]*velocity_conversion;
        p->xv[5] = vel[2]*velocity_conversion;
        p->mass = 1.0;
        p->tag = (unsigned int) aux.pid();
        	// This could overflow, as we allow 40-bit ids
    }
    void appendcell(char *c [[maybe_unused]], cell_header &current_cell, integer3 ijk, float vscale [[maybe_unused]]) {
        current_cell = cell_header(ijk, cpd, 1);
    }

  public:
    int sizeof_cell()     { return 0; }
    int sizeof_particle() { return sizeof(ICparticle); }

    OutputHeitmann() {
        velocity_conversion = 1.0;  // Leave it in Zspace, unit box units
    }
    ~OutputHeitmann(void) { }
};


class OutputRVdoublePID: public AppendArena {
  private:
    struct ICparticle {
        double pos[3];		// Global position, unit box
            double vel[3];		// Zspace, unit box
            uint64 tag;
    };

    void appendparticle(char *c, cell_header current_cell, posstruct pos, velstruct vel, auxstruct aux) {
        struct ICparticle *p = (struct ICparticle *)c;
#ifdef GLOBALPOS
        // This is for global box-referenced positions
        double3 cc = double3(0.0);
#else
        // If we instead have cell-referenced positions, then:
        double3 cc = CP->WrapCellCenter(current_cell.cellid());
#endif
        p->pos[0] = (double) pos.x+cc.x;
        p->pos[1] = (double) pos.y+cc.y;
        p->pos[2] = (double) pos.z+cc.z;
        p->vel[0] = vel.x;  // Leave it in Zspace, unit box units
        p->vel[1] = vel.y;
        p->vel[2] = vel.z;
        p->tag = aux.pid();

    }
    void appendcell(char *c [[maybe_unused]], cell_header &current_cell, integer3 ijk, float vscale [[maybe_unused]]) {
        current_cell = cell_header(ijk, cpd, 1);
    }

  public:
    int sizeof_cell()     { return 0; }
    int sizeof_particle() { return sizeof(struct ICparticle); }

    OutputRVdoublePID() {
            STDLOG(1,"Particle size: {:d}\n",sizeof_particle());
    }
    ~OutputRVdoublePID(void) { }
};


class OutputRVZel: public AppendArena {
  private:
    struct ICparticle {
        unsigned short i,j,k;
        float disp[3];      // Lagrangian displacement, unit box
        float vel[3];      // Zspace velocity, unit box
    };

    void appendparticle(char *c, cell_header current_cell, posstruct pos, velstruct vel, auxstruct aux) {
    struct ICparticle *p = (struct ICparticle *)c;
#ifdef GLOBALPOS
    // This is for global box-referenced positions
    double3 cc = double3(0.0);
#else
    // If we instead have cell-referenced positions, then:
    double3 cc = CP->WrapCellCenter(current_cell.cellid());
#endif
    integer3 ijk = aux.xyz();
    double3 zelpos = ZelPos(ijk);
    float3 disp = double3(pos) + cc - zelpos;
    disp -= disp.round();  // wrap to [-.5,+.5)

    p->disp[0] = disp.x;
    p->disp[1] = disp.y;
    p->disp[2] = disp.z;
    p->vel[0] = vel.x;  // Leave it in Zspace, unit box units
    p->vel[1] = vel.y;
    p->vel[2] = vel.z;
    p->i = ijk.x;
    p->j = ijk.y;
    p->k = ijk.z;

    }
    void appendcell(char *c [[maybe_unused]], cell_header &current_cell, integer3 ijk, float vscale [[maybe_unused]]) {
    current_cell = cell_header(ijk, cpd, 1);
    }

  public:
    int sizeof_cell()     { return 0; }
    int sizeof_particle() { return sizeof(struct ICparticle); }

    OutputRVZel() {}
    ~OutputRVZel(void) { }
};


#ifdef DONOTCOMPILE
/* This class takes a blank buffer and starts filling it with packed cells and particles.

Sample use:
*/

    SB->AllocateSpecificSize(FieldTimeSlice, slab, (SS->size(slab)+CP->cpd*(CP->cpd))*sizeof(PACKED));
    char *start = SB->GetSlabPtr(FieldTimeSlice,slab);
    AppendArena AA(start, SB->SlabSizeBytes(FieldTimeSlice,slab));

    FLOAT vel_to_zspace = 1.0;   // Factor to convert from canonical to zspace
    FLOAT kickfactor = 1.0;	   // Amount to unkick.

    for (integer3 ijk(slab,0,0); ijk.y<CP->cpd; ijk.y++)
        for (ijk.z=0;ijk.z<CP->cpd;ijk.z++) {
            Cell c = CP->GetCell(ijk);
            accstruct *acc = CP->AccCell(ijk);
            // Decide on vscale
            for (int p=0,float3 max = 0, float3 tmp; p<c.count(); p++) {
        	if ((tmp.x=fabs(c.vel[p].x))>max.x) max.x=tmp.x;
        	if ((tmp.y=fabs(c.vel[p].y))>max.y) max.y=tmp.y;
        	if ((tmp.z=fabs(c.vel[p].z))>max.z) max.z=tmp.z;
            }
            if (max.y>max.x) max.x = max.y;
            if (max.z>max.x) max.x = max.z;
            max.x*=vel_to_zspace;
            int vscale = ceiling(max.x);

            // Now pack the particles
            AA.addcell(ijk, vscale);
            for (p=0;p<c.count();p++) {
        	FLOAT3 vel = (c.vel[p]-acc[p]*kickfactor)*vel_to_zspace;
        	AA.addparticle(c.pos[p], vel, c.aux[p]);
            }
            AA.endcell();
        }

    SB->ResizeSlab(FieldTimeSlice, slab, AA.bytes_written());

    fs::path filename = ?; 	// Make the file name!
    SB->WriteArena(FieldTimeSlice, slab, IO_DELETE, IO_NONBLOCKING, filename);

#endif
