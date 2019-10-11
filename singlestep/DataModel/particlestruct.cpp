/* particlestruct.cpp

Base structures for the particle and cell data.  This is where
we choose between single and double precision; we provide 
types FLOAT and FLOAT3 for use throughout the code.

The auxillary class requires a fair number of methods to pack and
unpack it.

The number of particles in a cell must fit into a 32-bit signed integer,
and moreover we have routines that refer to the components individually.
So we actually are limited to about 680 million particles in a cell.


z PID: bits 0-14
y PID: 16-30
z PID: 32-46
Tagged: 48
Density: 49-58

L0: 15
L1: 31
L2: 47 (if we even use this?)
SubA: 59
SubB: 60
LightCone (3): 61-63

When we store into the output TaggedPID format, we apply a bit mask to zero all of the second set of bits. That's 0x07ff 7fff 7fff 7fff.

*/

#ifndef INCLUDE_PARTICLESTRUCT
#define INCLUDE_PARTICLESTRUCT

#define posstruct FLOAT3
#define velstruct FLOAT3
#define acc3struct FLOAT3
        // The far-field code operates on float3's.

#ifdef COMPUTE_FOF_DENSITY
#define accstruct FLOAT3p1
// We will pass the accelerations around as a float3p1, in which
// the 4th element carries density information and hence should
// not be normalized in the same way.  Adding a float3p1 to a float3
// produces a float3p1 without change to the w element.  
// Typecasting to (FLOAT3) will strip the w element.
#else
#define accstruct FLOAT3
#endif

#ifndef uint16
#define uint16 uint16_t
#endif

#define AUXXPID  (uint64)  0x7fff	// bits 0 - 14 (right-most)
#define AUXYPID  (uint64)  0x7fff0000  // bits 16 - 30
#define AUXZPID  (uint64)  0x7fff00000000  // bits 32 - 46

#define AUXTAGGEDBIT 48llu //Has the particle been tagged
#define AUXDENSITYZEROBIT 49llu //The density bits are 49-58. 
#define AUXDENSITY (uint64) 0x7fe000000000000 //bits 49-58
#define AUXINL0BIT 15llu //Is the particle in a level 0 group
#define AUXINL1BIT 31llu //Is the particle in a level 1 group
#define AUXINL2BIT 47llu //Is the particle in a level 2 group

#define AUXTAGGABLE_A_BIT 59llu //Can this particle be tagged in subsample A? 
#define AUXTAGGABLE_B_BIT 60llu //Can this particle be tagged in subsample B? 


#define AUXLCZEROBIT 61llu		// The LC bits are 61,62,63.
#define AUXLC  (uint64)0xe000000000000000	// The next three bits.

#define AUX_PID_TAG_DENS_MASK 0x07ff7fff7fff7fff //this masks out everything except for bits 0-14 (x pid), 16-30 (y pid), 32-46 (z pid), 48 (tagged), 49-58 (density).
#define AUXPIDMASK 0x7fff7fff7fff

class auxstruct {
public:
    uint64 aux;
    // unsigned char lightcones; //1 bit for each of 8 lightcones

    // Methods to extact items, e.g.,
   void clear() { aux = 0; }   // lightcones =0;}

    // We will provide at least one ID, suitable for the particles.
    // This is required for the LPT implementation.
    uint64 pid() { return aux&AUXPIDMASK; }

    integer3 xyz() { 
        integer3 xyz; 
        xyz.x = (pid() & AUXXPID); 
        xyz.y = (pid() & AUXYPID) >> 16;
        xyz.z = (pid() & AUXZPID) >> 32;
        return xyz; 
    }

    void setpid(uint16 _pid) { 
        uint16 pid[3] = {_pid, _pid, _pid};
        setpid(pid); 
    }

    void setpid(integer3 _pid) { 
        uint16 max = (uint16) AUXXPID; 
            assert(_pid.x <= max and _pid.y <= max and _pid.z <= max);
        uint16 pid[3] = {(uint16) _pid.x, (uint16) _pid.y, (uint16) _pid.z}; 
        setpid(pid); 
    }

    void setpid(uint16 * _pid) { 
        uint16 max = (uint16) AUXXPID; 
           assert(_pid[0] <= max and _pid[1] <= max and _pid[2] <= max);
        setpid((uint64) _pid[0] | (uint64) _pid[1]<<16| (uint64) _pid[2] <<32);
    }

    void setpid(uint64 _pid) {
        aux = _pid | (aux &~ AUXPIDMASK); 
    }
    // We will provide a group ID too; this may overwrite the PID.
    uint64 gid() { return pid(); }
    void setgid(uint64 gid) { setpid(gid); }

    // We expose lightconemask() publicly because one might 
    // not want to construct for every particle.

    // Light cones need 1 byte
    inline static uint64 lightconemask(int number) {
        assertf(number<8 && number>=0, "Lightcone number lcn = %d must satisfy 0 <= lcn < 8.", number);
        return (uint64)1 << (number+AUXLCZEROBIT);
    }

    inline bool lightconedone(uint64 mask) {
        assert (mask<=AUXLC && mask >= AUXPIDMASK);  // better way to do this...
        return (aux&mask);
    }
    inline bool lightconedone(int number) {
        return lightconedone(lightconemask(number));
    }

    inline void setlightconedone(uint64 mask) {
        aux |= mask;
    }
    inline void setlightconedone(int number) {
        setlightconedone(lightconemask(number));
    }

    inline void set_taggable_subA() {
        // The TAGGABLE SUBA bit should be set at the beginning of the sim and not changed.
        aux |= ((uint64)1 << AUXTAGGABLE_A_BIT);
    }

    // Group and subsample related bits
    inline void set_taggable_subB() {
        // The TAGGABLE SUBB bit should be set at the beginning of the sim and not changed.
        aux |= ((uint64)1 << AUXTAGGABLE_B_BIT);
    }

    #define TAGGABLE_SUB_A 1
    #define TAGGABLE_SUB_B 2
    inline int is_taggable() { // >> and mask
        return ((aux >> AUXTAGGABLE_A_BIT) & 0x3); //returns > 0 if something is taggable, 0 otherwise. 
    }
    inline void set_tagged() {
        // The TAGGED bit is a lasting tag, once set.
        aux |= ((uint64)1 << AUXTAGGEDBIT);
    }
    inline bool is_tagged() {
        return aux & ((uint64)1 << AUXTAGGEDBIT);
    }
    inline void set_density(uint64 _density){
        assert(_density < (AUXDENSITY >> AUXDENSITYZEROBIT)); 
        aux |= ((uint64) _density << AUXDENSITYZEROBIT);

    }
    inline void reset_L01_bits() {
        // We need to be able to unset these bits each time we run groupfinding
        uint64 mask = ((uint64)1 << AUXINL0BIT) + ((uint64)1 << AUXINL1BIT);
        aux &= ~mask;
    }

    inline void reset_L1_bit() {
        // We need to be able to unset the L1 bit to output L1 particles in halos that are too small. 
        uint64 mask = ((uint64)1 << AUXINL1BIT);
        aux &= ~mask;
    }

    inline void set_L0() {
        aux |= ((uint64)1 << AUXINL0BIT);
    }
    inline bool is_L0() {
        return aux & ((uint64)1 << AUXINL0BIT);
    }
    inline void set_L1() {
        aux |= ((uint64)1 << AUXINL1BIT);
    }
    inline bool is_L1() {
        return aux & ((uint64)1 << AUXINL1BIT);
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
        STDLOG(0, "Bad 'startindex' in cellinfo: %u\n", startindex);
        return 0;
    }
    if(count<0){
        STDLOG(0, "Bad 'count' in cellinfo: %d\n", count);
        return 0;
    }
    if(active<0){
        STDLOG(0, "Bad 'active' in cellinfo: %u\n", active);
        return 0;
    }
        if (startindex+count > slabsize){
        STDLOG(0, "'startindex+count' (%u + %d = %d) > 'slabsize' (%u) in cellinfo\n", startindex, count, startindex+count, slabsize);
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
// Unfortunately, this does not include group lists,
// because we don't always have those.

class Cell {
public:
    // The cell index number (useful for converting from local to global positions)
    integer3 ijk;
    // Pointer to the cell info (not a copy of it!)
    cellinfo *ci;
    // Pointers to the lists of positions, velocities, and auxiliaries
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
