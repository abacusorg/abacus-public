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
#include "pack14_storage.cpp"


class AppendArena {
  private:
    // This is an Abstract Base Class.
    // The following 4 routines must be provided in the derived class.
    // Of course, the derived class must also provide a constructor and destructor.
    
    //FIXME: Since appending particles is a quick function we are calling many times, 
    //  we may want to make the polymorphism compile-time (using templates) rather than
    //  looking up virtual functions at runtime

    virtual void appendparticle(char *c, posstruct pos, velstruct vel, auxstruct aux) =0;
    // Place one particle at the location c
    virtual void appendcell(char *c, integer3 ijk, float vscale) =0;
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

    char *arena;	// Where we will write next
    long long int size_available, bytes;
    int cpd;		// CPD to be used in mapping particle positions
    FLOAT VelCanonical_to_ZSpace; 	// Typically 1/ReadState.VelZSpace_to_Canonical

    cell_header current_cell;

    void initialize(SlabType type, int slab, int _cpd, double VelZSpace_to_Canonical) {
        arena = LBW->ReturnIDPtr(type,slab);
	size_available = LBW->IDSizeBytes(type,slab);
	bytes = 0;
	cpd = _cpd;
	VelCanonical_to_ZSpace = 1.0/VelZSpace_to_Canonical;
	endcell();
    }

    // ~AppendArena(void) { }

    char *ptr_to_end() {
	// In case you need to know where you are!
        return arena;
    }
    long long int bytes_written() {
        return bytes;
    }

    void addcell(integer3 ijk, float vscale) {
	// Append this cell and make this cell active.
	int size = sizeof_cell();
	assertf(bytes+size<=size_available, "AppendArena() has run out of space!\n");
	appendcell(arena, ijk, vscale);
	arena+=size; bytes+=size;
    }
    void endcell() { 
	// This makes it illegal to add more particles
    	current_cell.vscale = 0; 
    }

    void addparticle(posstruct pos, velstruct vel, auxstruct aux) {
	// Append this particle
	assertf(current_cell.islegal(), "addparticle() called without an active cell\n");
	int size = sizeof_particle();
	assertf(bytes+size<=size_available, "AppendArena() has run out of space!\n");
	// We were given the particle in code units, need to convert to Zspace.
	vel *= VelCanonical_to_ZSpace;
	// The appendparticle() function is responsible for any further conversions,
	// e.g., from local to global positions or to km/s or from unit box to unit cell..
        appendparticle(arena, pos, vel, aux);
	arena+=size; bytes+=size;
    }

    void addheader(const char *c) {
        // Write the given string to the arena, without the final \0.
	int size = 0;
	while (*c!='\0') {
	    *arena = *c;
	    arena++; c++; bytes++;
	}
    }
    void addheader(char *c) {
        addheader((const char *)c);
    }

    void finalize_header() {
        // terminate the ParseHeader header.  
	// Pad out the header to a multiple of 4096 bytes just to keep the data DIO aligned.
	long long int size = bytes_written()+2;  // including the header
	int pad = 4096-size%4096; 
	if (pad>0) { arena[0] = '\n'; }
	for (int n=1; n<pad; n++) { arena[n] = ' '; }
	arena += pad; bytes += pad;
	arena[0] = ''; arena[1] = '\n'; arena+=2; bytes+=2;
    }
};

// =================================================================
// Here's an example of a derived class that would use the packed functions.
//
// FIXME: The different output formats should probably be in their own files

#define PACKED pack14
class OutputPacked: public AppendArena {
  private:
    void appendparticle(char *c, posstruct pos, velstruct vel, auxstruct aux) {
        PACKED *p = (PACKED *) c;
	// p->pack_global(pos, vel, aux.pid(), current_cell);
	p->pack(pos, vel, aux.pid(), current_cell);
    }
    void appendcell(char *c, integer3 ijk, float vscale) {
	// We're given vscale in Zspace unit-box units, same as velocities.
	// But we need to hand it to the pack14 method in unit-cell units
        PACKED *p = (PACKED *) c;
	int vs = ceil(vscale*cpd); 
	if (vs<=0) vs = 10;    // Probably just not initialized correctly
	current_cell = p->pack_cell(ijk, cpd, vs);
    }
    float velocity_conversion;

  public:
    int sizeof_cell()     { return sizeof(PACKED); }
    int sizeof_particle() { return sizeof(PACKED); }

    OutputPacked() {
        // Use ReadState to figure out the correct conversion of the
	// velocity!  The pack14 class assumes velocites are in 
	// redshift-space units; it *then* applies the given vscale.
	velocity_conversion = 1.0;
    }
    ~OutputPacked(void) { }
};
#undef PACKED

// =================================================================
// And here's an example to put things into simple RVdouble

class OutputRVdouble: public AppendArena {
  private:
    struct ICparticle {
        double pos[3];		// Global position, unit box
	double vel[3];		// Zspace, unit box
    };

    float velocity_conversion;
    
    void appendparticle(char *c, posstruct pos, velstruct vel, auxstruct aux) {
	struct ICparticle *p = (struct ICparticle *)c;
#ifdef GLOBALPOS
	// This is for global box-referenced positions
	// p->pos[0] = pos.x;
	// p->pos[1] = pos.y;
	// p->pos[2] = pos.z;
	double3 cc = double3(0.0);
#else
	// If we instead have cell-referenced positions, then:
	double3 cc = PP->WrapCellCenter(current_cell.cellid());
#endif
	p->pos[0] = (double) pos.x+cc.x;
	p->pos[1] = (double) pos.y+cc.y;
	p->pos[2] = (double) pos.z+cc.z;
	p->vel[0] = vel.x*velocity_conversion;
	p->vel[1] = vel.y*velocity_conversion;
	p->vel[2] = vel.z*velocity_conversion;
    }
    void appendcell(char *c, integer3 ijk, float vscale) {
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


//=====================================================================
//The format we will use for lightcones
class OutputLightcone: public AppendArena {
  private:
    struct ICparticle {
        float pos[3];		// Global position, unit box
        float vel[3];		// Zspace, unit box
        auxstruct aux;
    };

    struct CellHeader{
    	long long int ijk[3];
    	long long int np; //number of particles in this cell
    };
    long long int * np;

    float velocity_conversion;
  public:
    void appendparticle(char *c, posstruct pos, velstruct vel, auxstruct aux) {
	struct ICparticle *p = (struct ICparticle *)c;
	p->pos[0] = (float) pos.x;
	p->pos[1] = (float) pos.y;
	p->pos[2] = (float) pos.z;
	p->vel[0] = vel.x ;
	p->vel[1] = vel.y;
	p->vel[2] = vel.z;
	p->aux = aux;
	assert(np !=0);
	(*np)++;
    }


    void appendcell(char *c, integer3 ijk, float vscale) {

    	current_cell = cell_header(ijk, cpd, 1);
    	CellHeader *cc = (CellHeader *) c;
    	cc->ijk[ 0] = ijk.x;
    	cc->ijk[ 1] = ijk.y;
    	cc->ijk[ 2] = ijk.z;
    	np = &(cc->np);
    	*np = 0;
    }




    int sizeof_cell()     { return sizeof(CellHeader); }
    int sizeof_particle() { return sizeof(struct ICparticle); }

    OutputLightcone() {
	velocity_conversion = 1.0;  // Leave it in Zspace, unit box units
	np =0;
    }
    ~OutputLightcone(void) { }
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
    
    void appendparticle(char *c, posstruct pos, velstruct vel, auxstruct aux) {
	ICparticle *p = (ICparticle *)c;
#ifdef GLOBALPOS
	// This is for global box-referenced positions
	// p->xv[0] = pos[0];
	// p->xv[2] = pos[1];
	// p->xv[4] = pos[2];
	double3 cc = double3(0.0);
#else
	// If we instead have cell-referenced positions, then:
	double3 cc = PP->WrapCellCenter(current_cell.cellid());
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
    void appendcell(char *c, integer3 ijk, float vscale) {
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


class OutputRVdoubleTag: public AppendArena {
  private:
    struct ICparticle {
        double pos[3];		// Global position, unit box
	    double vel[3];		// Zspace, unit box
	    uint64 tag;
    };

    void appendparticle(char *c, posstruct pos, velstruct vel, auxstruct aux) {
	struct ICparticle *p = (struct ICparticle *)c;
#ifdef GLOBALPOS
	// This is for global box-referenced positions
	double3 cc = double3(0.0);
#else
	// If we instead have cell-referenced positions, then:
	double3 cc = PP->WrapCellCenter(current_cell.cellid());
#endif
	p->pos[0] = (double) pos.x+cc.x;
	p->pos[1] = (double) pos.y+cc.y;
	p->pos[2] = (double) pos.z+cc.z;
	p->vel[0] = vel.x;  // Leave it in Zspace, unit box units
	p->vel[1] = vel.y;
	p->vel[2] = vel.z;
	p->tag = aux.pid();
	
    }
    void appendcell(char *c, integer3 ijk, float vscale) {
	current_cell = cell_header(ijk, cpd, 1);
    }

  public:
    int sizeof_cell()     { return 0; }
    int sizeof_particle() { return sizeof(struct ICparticle); }

    OutputRVdoubleTag() {
	    STDLOG(0,"Particle size: %d\n",sizeof_particle());
    }
    ~OutputRVdoubleTag(void) { }
};





#ifdef DONOTCOMPILE
/* This class takes a blank buffer and starts filling it with packed cells and particles.

Sample use:
*/

    LBW->AllocateSpecificSize(TimeSlice, slab, (Slab->size(slab)+PP->cpd*(PP->cpd))*sizeof(PACKED));
    char *start = LBW->ReturnIDPtr(TimeSlice,slab);
    AppendArena AA(start, LBW->IDSizeBytes(TimeSlice,slab));

    FLOAT vel_to_zspace = 1.0;   // Factor to convert from canonical to zspace
    FLOAT kickfactor = 1.0;	   // Amount to unkick.

    for (integer3 ijk(slab,0,0); ijk.y<PP->cpd; ijk.y++) 
	for (ijk.z=0;ijk.z<PP->cpd;ijk.z++) {
	    Cell c = PP->GetCell(ijk);
	    accstruct *acc = PP->AccCell(ijk);
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

    LBW->ResizeSlab(TimeSlice, slab, AA.bytes_written());

    char filename[1024]; 	// Make the file name!
    LBW->WriteArena(TimeSlice, slab, IO_DELETE, IO_NONBLOCKING, filename);

#endif
