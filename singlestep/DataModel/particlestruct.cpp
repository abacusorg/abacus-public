/* particlestruct.cpp

Base structures for the particle and cell data.  This is where
we choose between single and double precision; we provide 
types FLOAT and FLOAT3 for use throughout the code.

The auxillary class requires a fair number of methods to pack and
unpack it.

The number of particles in a cell must fit into a 32-bit signed integer,
and moreover we have routines that refer to the components individually.
So we actually are limited to about 680 million particles in a cell.

*/

#ifndef INCLUDE_PARTICLESTRUCT
#define INCLUDE_PARTICLESTRUCT

// So that we can easily adjust to double precision for kinematics
#ifdef DOUBLEPRECISION
#define FLOAT double
#define FLOAT3 double3		
#else 
#define FLOAT float	
#define FLOAT3 float3
#endif

#define posstruct FLOAT3
#define velstruct FLOAT3
#define accstruct FLOAT3

#define AUXPID (uint64)  0xffffffffff	// The lower 5 bytes, bits 0..39
#define AUXLCZEROBIT 40		// The LC bits are 40..47
#define AUXLC  (uint64)0xff0000000000	// The next byte
#define AUXTAGABLEBIT 48llu //Can the particle be tagged.
#define AUXTAGGEDBIT 49llu //Has the particle been tagged
#define AUXINL0BIT 50llu //Is the particle in a level 0 group
#define AUXINL1BIT 51llu //Is the particle in a levl 1 group

class auxstruct {
public:
    uint64 aux;
    // unsigned char lightcones; //1 bit for each of 8 lightcones

    // Methods to extact items, e.g.,
    void clear() { aux = 0; }   // lightcones =0;}

    // We will provide at least one ID, suitable for the particles.
    // This is required for the LPT implementation.
    // The PID will fill the lower 5 bytes.
    // The upper 3 bytes may be used for other items.
    uint64 pid() { return aux&AUXPID; }
    uint64 key() { return pid(); }
    void setpid(uint64 pid) { 
    	assertf(pid<=AUXPID,
		"PID %d is too big\n", pid); 
        aux = pid + (aux&~AUXPID);
    }

    // We will provide a group ID too; this may overwrite the PID.
    uint64 gid() { return pid(); }
    void setgid(uint64 gid) { setpid(gid); }

    // We expose lightconemask() publicly because one might 
    // not want to construct for every particle.

    // Light cones need 1 byte
    static uint64 lightconemask(int number) {
        assertf(number<8 && number>=0, "Lightcone number lcn = %d must satisfy 0 <= lcn < 8.", number);
        return (uint64)1 << (number+AUXLCZEROBIT);
    }

    bool lightconedone(uint64 mask) {
        assert (mask<=AUXLC && mask >= AUXPID);  // better way to do this...
        return (aux&mask);
    }
    bool lightconedone(int number) {
        return lightconedone(lightconemask(number));
    }

    void setlightconedone(uint64 mask) {
        aux |= mask;
    }
    void setlightconedone(int number) {
        setlightconedone(lightconemask(number));
    }


/* OLD CODE
    bool lightconedone(int number) {
    	assert (number < 8 && number >=0);
    	unsigned char mask = 1<<number;
    	return ((mask & lightcones) & mask);
    }
    void setlightconedone(int number){
    	unsigned char mask = 1<<number;
    	lightcones |= mask;
    }
*/

    
};

class cellinfo { 
public:
    uint64 startindex;     // The particle-count offset in the slab list of particles
    int count;          // The total number of particles in the cell
    int active;         // The number of active particles
    FLOAT mean_square_velocity;     // The mean of |vel|^2 (no adjustment for mean v)
    FLOAT max_component_velocity;      // The maximum x,y,z component of v
    FLOAT max_component_acceleration;  // The maximum x,y,z component of a

    void makenull() {
        memset(this,0,sizeof(cellinfo));
        startindex = 0;
        count = 0;
        active = 0;
        mean_square_velocity = 0;
        max_component_velocity  = 0;
        max_component_acceleration = 0;
    }
    int legalvalue(uint64 slabsize) {
        // Do a sanity check on the cellinfo values; return 1 if ok, 0 if not.
	// slabsize is the number of particles in the slab.
	if (!isfinite(startindex)){
        STDLOG(0, "Bad 'startindex' in cellinfo: %llu\n", startindex);
        return 0;
    }
    if(count<0){
        STDLOG(0, "Bad 'count' in cellinfo: %d\n", count);
        return 0;
    }
    if(active<0){
        STDLOG(0, "Bad 'active' in cellinfo: %llu\n", active);
        return 0;
    }
	if (startindex+count > slabsize){
        STDLOG(0, "'startindex+count' (%llu + %d = %lld) > 'slabsize' (%llu) in cellinfo\n", startindex, count, startindex+count, slabsize);
        return 0;
    }
	if (active>count){
        STDLOG(0, "'active' (%d) > 'count' (%d) in cellinfo\n", active, count);
        return 0;
    }
	return 1;
    }
};


// Define the cell class to have a compact way to handle the particles.
// Unfortunately, this does not include accelerations and group lists,
// because we don't always have those.

class Cell {
public:
    // The cell index number (useful for converting from local to global positions)
    integer3 ijk;
    // Pointer to the cell info (not a copy of it!)
    cellinfo *ci;
    // Pointers to the lists of positions, velocities, and auxillaries
    // These lists run 0 .. ci->count-1
    posstruct *pos;
    velstruct *vel;
    auxstruct *aux;
    accstruct *acc = NULL;

    inline int count() { return ci->count; }
    inline int active() { return ci->active; }

    inline void swap(int a, int b) {
        // Swap two particles in the list
        /*assertf(a>=0&&a<ci->count,
            "Particle a = %d is out of range %d\n", a, ci->count);
        assertf(b>=0&&b<ci->count,
            "Particle b = %d is out of range %d\n", b, ci->count);*/
        posstruct ptmp; ptmp = pos[b]; pos[b] = pos[a]; pos[a] = ptmp; 
        velstruct vtmp; vtmp = vel[b]; vel[b] = vel[a]; vel[a] = vtmp; 
        auxstruct atmp; atmp = aux[b]; aux[b] = aux[a]; aux[a] = atmp;
        if(acc != NULL){
            accstruct actmp; actmp = acc[b]; acc[b] = acc[a]; acc[a] = actmp;
        }
    }

    // Routines to copy information in and out of the lists
    inline void getparticle(int n, posstruct *p, velstruct *v, auxstruct *a) {
        *p = pos[n]; *v = vel[n]; *a = aux[n];
    }
    inline void putparticle(int n, posstruct *p, velstruct *v, auxstruct *a) {
        pos[n] = *p; vel[n] = *v; aux[n] = *a;
    }
};


#endif // INCLUDE_PARTICLESTRUCT
