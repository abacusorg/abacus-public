/* particlestruct.cpp

Base structures for the particle and cell data.  This is where
we choose between single and double precision; we provide 
types FLOAT and FLOAT3 for use throughout the code.

The auxillary class requires a fair number of methods to pack and
unpack it.

The number of particles in a cell must fit into a 32-bit signed integer,
and moreover we have routines that refer to the components individually.
So we actually are limited to about 680 million particles in a cell.


x PID: bits 0-14
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

density: 64-95 (full-precision dens for group finding)

For the LPT steps, we repurpose the following bits:
vx: 48-63
vy: 64-79
vz: 80-95

When we store into the output TaggedPID format, we only output
from the first 64 bits, and apply a bit mask to zero all of the
second set of bits. That's 0x07ff 7fff 7fff 7fff.

*/
#include <bitset>


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

#define NUMLIGHTCONES 3

#define AUXLCZEROBIT 61llu		// The LC bits are 61,62,63.
#define AUXLC  (uint64)0xe000000000000000	// The next three bits.

#define AUX_PID_TAG_DENS_MASK 0x07ff7fff7fff7fff //this masks out everything except for bits 0-14 (x pid), 16-30 (y pid), 32-46 (z pid), 48 (tagged), 49-58 (density).
#define AUXPIDMASK 0x7fff7fff7fff

#define NAUXPIDBITS 45  // total of 45 bits used for PIDs
#define PIDBITGAP 1  // one bit in between each PID segment


#define AUX_LPTVXZEROBIT 48
#define AUX2_LPTVYZEROBIT 0
#define AUX2_LPTVZZEROBIT 16
#define AUX_LPTVX ((uint64) 0xFFFF << AUX_LPTVXZEROBIT)
#define AUX2_LPTVY ((uint64) 0xFFFF << AUX2_LPTVYZEROBIT)
#define AUX2_LPTVZ ((uint64) 0xFFFF << AUX2_LPTVZZEROBIT)

class __attribute__((packed)) auxstruct {
private:
    uint64 aux;

    union {
        uint32 aux2;
        float dens;
    };

public:

   void clear() { aux = 0; aux2 = 0; }

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

    void setpid(integer3 _pid) { 
        assert(_pid.x <= AUXXPID && _pid.y <= AUXXPID && _pid.z <= AUXXPID);
        _setpid(_pid.x, _pid.y, _pid.z);
    }

    void _setpid(uint64 px, uint64 py, uint64 pz) { 
        _setpidbits(px | py<<16 | pz<<32);
    }

    void _setpidbits(uint64 _pid) {
        aux = _pid | (aux & ~AUXPIDMASK);
    }

    // Take a pid and distribute its values to the three segements used for PIDs
    void packpid(uint64 _pid){
        assert(_pid <= ((uint64) 1 << NAUXPIDBITS));
        uint64 pid = (_pid & AUXXPID) | (_pid << PIDBITGAP & AUXYPID) | (_pid << 2*PIDBITGAP & AUXZPID);
        _setpidbits(pid);
    }


    // We will provide a group ID too; this may overwrite the PID.
    uint64 gid() { return pid(); }
    void setgid(uint64 gid) { _setpidbits(gid); }

    // We expose lightconemask() publicly because one might 
    // not want to construct for every particle.

    // Light cones need 1 byte
    inline static uint64 lightconemask(int number) {
        assertf(number<NUMLIGHTCONES && number>=0, "Lightcone number lcn = %d must satisfy 0 <= lcn < %d.", number, NUMLIGHTCONES);
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

    inline void clearLightCone() {
        uint64 mask = AUXLC;
        aux &= ~mask;
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

    #define TAGGABLE_SUB_A 0x1
    #define TAGGABLE_SUB_B 0x2
    inline int is_taggable() { // >> and mask
        return ((aux >> AUXTAGGABLE_A_BIT) & 0x3); //returns > 0 if something is taggable, 0 otherwise
    }
    inline void set_tagged() {
        // The TAGGED bit is a lasting tag, once set.
        aux |= ((uint64)1 << AUXTAGGEDBIT);
    }
    inline bool is_tagged() {
        return aux & ((uint64)1 << AUXTAGGEDBIT);
    }

    // Store a lossy compression of the density from the raw value (in code units)
    // There is a chance that a mismatch in precision between the GPU and CPU codes 
    // could lead to a negative value of the density when the only particle is the self-particle.
    // Hence, we take the absolute value.
    inline void set_compressed_density(FLOAT rawdensity){
        _pack_density((uint64) round(std::sqrt(std::abs(rawdensity) * WriteState.invFOFunitdensity)));
    }

    inline void _pack_density(uint64 _density){
        assert(_density < (AUXDENSITY >> AUXDENSITYZEROBIT)); 
        aux = ( (uint64) _density << AUXDENSITYZEROBIT )  | (aux &~ AUXDENSITY); 
    }

    inline FLOAT get_compressed_density(){
        uint64 d = _unpack_density();
        return d*d * WriteState.FOFunitdensity;
    }

    inline uint64 _unpack_density(){
        return (aux | AUXDENSITY) >> AUXDENSITYZEROBIT;
    }

    // Routines to set/get the full-precision density for group finding
    inline void set_density(FLOAT rawdensity){
        dens = (float) rawdensity;
    }

    inline float get_density(){
        return dens;
    }
    
    inline void reset_L01_bits() {
        // We need to be able to unset these bits each time we run groupfinding
        uint64 mask = ((uint64)1 << AUXINL0BIT) + ((uint64)1 << AUXINL1BIT);
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

    inline uint64 get_aux_pid_dens_tagged(){
        return aux & AUX_PID_TAG_DENS_MASK;
    }

    // Routines to support 2LPT, where we save the velocity in the aux

    inline void set_velocity(velstruct delta, FLOAT invscale){
        // vel is the velocity to store, usually as an offset
        uint16 vx = _pack_float(delta.x, invscale);
        uint16 vy = _pack_float(delta.y, invscale);
        uint16 vz = _pack_float(delta.z, invscale);
        
        aux = (aux & ~AUX_LPTVX) | ((uint64) vx << AUX_LPTVXZEROBIT);
        aux2 = (aux2 & ~(AUX2_LPTVY | AUX2_LPTVZ)) |
                ((uint64) vy << AUX2_LPTVYZEROBIT) |
                ((uint64) vz << AUX2_LPTVZZEROBIT);
    }

    uint16 _pack_float(FLOAT f, FLOAT invfscale){
        // packs f as a 16 bit int, stored as a ratio relative to scale
        // f must be in [-scale,scale]
        FLOAT ratio = f*invfscale; // [-1,1]
        int iscale = (1<<15)-1;
        int enc = (int) (ratio*iscale) + iscale;  // [0,1<<16-2]
        
        // clamp modest overflow
        if(enc == -1) enc = 0;
        if(enc == 2*iscale + 1) enc = 2*iscale;  // we could represent this overflow, but let's be symmetric
        
        assertf(enc >= 0 && enc < 1<<16, "Cannot pack float %g with invscale %g\n", f, invfscale);
        return (uint16) enc;
    }

    inline void zero_velocity(){
        aux &= ~AUX_LPTVX;
        aux2 &= ~(AUX2_LPTVY | AUX2_LPTVZ);
    }

    inline velstruct get_velocity(FLOAT scale){
        FLOAT vx = _unpack_float((aux & AUX_LPTVX) >> AUX_LPTVXZEROBIT, scale);
        FLOAT vy = _unpack_float((aux2 & AUX2_LPTVY) >> AUX2_LPTVYZEROBIT, scale);
        FLOAT vz = _unpack_float((aux2 & AUX2_LPTVZ) >> AUX2_LPTVZZEROBIT, scale);

        return velstruct(vx, vy, vz);
    }

    FLOAT _unpack_float(uint16 enc, FLOAT fscale){
        int iscale = (1<<15)-1;
        FLOAT inviscale = 1.0/iscale;
        FLOAT f = ((int) enc - iscale) * inviscale;
        return f*fscale;
    }

    std::string tostring(){
        std::stringstream stream;
        stream << "0x" << std::hex << aux << std::hex << aux2;
        return stream.str();
    }    
};

static_assert(sizeof(auxstruct) == 12, "unexpected auxstruct size");

class cellinfo {
    /* With regard to ghost cells, our CellInfo slabs will always be sized to hold
    ghosts.  But some of our slabs, like AccSlab, won't. So each cellinfo struct
    has two counts: one for the starting offset in slabs with ghosts, and one
    for the slabs without.

    The alternatives would be to make a new CellInfoWithGhosts slab type, or
    simply require all slabs to have ghosts, even if they are null.
    */
public:
    uint64 startindex;     // The particle-count offset in the slab list of particles
    uint64 startindex_with_ghost;     // The offset in slabs that also hold ghost cells
    int count;          // The total number of particles in the cell
    int active;         // The number of active particles
    FLOAT mean_square_velocity;     // The mean of |vel|^2 (no adjustment for mean v)
    FLOAT max_component_velocity;      // The maximum x,y,z component of v
    FLOAT max_component_acceleration;  // The maximum x,y,z component of a

    void makenull() {
        memset(this,0,sizeof(cellinfo));
        startindex = 0;
        startindex_with_ghost = 0;
        count = 0;
        active = 0;
        mean_square_velocity = 0;
        max_component_velocity  = 0;
        max_component_acceleration = 0;
    }
    int legalvalue(uint64 slabsize, uint64 slabsize_with_ghost, int isghost) {
        // Do a sanity check on the cellinfo values; return 1 if ok, 0 if not.
        // slabsize is the number of particles in the slab.
        if(count<0){
            STDLOG(0, "Bad 'count' in cellinfo: %d\n", count);
            return 0;
        }
        if(active<0){
            STDLOG(0, "Bad 'active' in cellinfo: %u\n", active);
            return 0;
        }
        if(mean_square_velocity<0){
            STDLOG(0, "Bad 'mean_square_velocity' in cellinfo: %u\n", mean_square_velocity);
            return 0;
        }

        // startindex is meaningless in ghost cells
        if (!isghost && startindex+count > slabsize){
            STDLOG(0, "'startindex+count' (%u + %d = %d) > 'slabsize' (%u) in cellinfo\n", startindex, count, startindex+count, slabsize);
            return 0;
        }
        
        if (startindex_with_ghost + count > slabsize_with_ghost){
            STDLOG(0, "'startindex_with_ghost+count' (%u + %d = %d) > 'slabsize_with_ghost' (%u) in cellinfo\n",
                startindex_with_ghost, count, startindex_with_ghost+count, slabsize_with_ghost);
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
