/** \file This file contains the routines to do FOF on a consecutive list of 
particles, either from a cell or from a L0 or L1 gathering.
Optionally, this can perform a permutation on the particle set.

There are several different algorithms in this file, as well as various
AVX tricks to accelerate the computation of distances.

As written, the input particle set cannot exceed 2**23 particles.
*/

//#define FOFTimer DummyTimer
#define FOFTimer STimer

#ifdef TEST

#include <stdlib.h>
#include <cmath>
#include <algorithm>
#include <vector>
#include <cstdint>
#include <sys/time.h>
#include "STimer.cc"
#ifdef OMP
    #include <omp.h>
#else
    int omp_get_num_threads() { return 1; }
    int omp_get_thread_num() { return 0; }
#endif
#include <unistd.h>
#include <string.h>
#include "stdlog.cc"

#undef FOFTimer
#define FOFTimer STimer

#include "promote_numeric.h"
#include "threevector.hh"

typedef float3 posstruct;
typedef float3 velstruct;
typedef float3 accstruct;
typedef unsigned long long int  auxstruct;

#endif 

#if defined(AVXFOF) or defined(AVX512DIRECT)
#include <immintrin.h>
#endif

class DummyTimer {
  public:
    inline void Start() { }
    inline void Stop() { }
    inline double Elapsed() { return 0.0; }
    inline struct timespec get_timer() { struct timespec retval; return retval; }
    inline void increment(struct timespec dt) { }
};

//  Now we get into the specifics of the FOF code.

#ifdef AVXFOF
    typedef float FOFloat;
    // It is assumed in the AVX code that FOFparticle is 16 bytes!
    // Importantly, this will override the DOUBLE_PRECISION flag.
#else
    typedef FLOAT FOFloat;
    // For now, revert to the decision of the rest of the code.
#endif

#define FOF_RESCALE 1e15

/** We will copy the positions and index number into a float4, then
do all of the rest of the work swapping this compact object.
Permutation of the original array is done only optionally.

  We are storing the cell index as a float.  Scaling the positions
  up by a huge number, so that we can do positional differences
  as a float4 norm with SSE, trusting to underflow the indices. 

  OLD: As written, this limits the number of particles in a cell to not
  overflow a 23-bit int.  However, this could be adjusted by setting
  the floating exponent sign bit to be negative and then using bits
  in the exponent.  Not implemented at present!

  NEW: We're now just shoving the bits of the int32 into the 4th float;
  for numbers 0 to 600M, this is <1e-16 because of the way the float exponents are handled.
  Keeping the scaling up by FOF_RESCALE, even though this could have been skipped.

  We usually do this in single precision!
  Particularly if we want to consider some AVX calls for diff2.
  Doing this as aligned float4's will also allow AVX for min/max
  computation of bounding boxes.
*/

union floatint {
    FOFloat w;
    int n;
};


class alignas(16) FOFparticle {

  public:
    FOFloat x,y,z; 
    floatint fi;
    FOFparticle(posstruct p, int _n) {
	x = p.x*FOF_RESCALE;
	y = p.y*FOF_RESCALE;
	z = p.z*FOF_RESCALE;
    fi.n = _n;
	// n = (FOFloat)(_n);
    }
    FOFparticle() { }
    // inline int index() { return (int)(n); }
    inline int index() { return fi.n; }
    inline FOFloat diff2v(FOFparticle *q) {
        FOFloat dx = q->x-x;
        FOFloat dy = q->y-y;
        FOFloat dz = q->z-z;
        FOFloat dn = q->fi.w-fi.w;
	return dx*dx+dy*dy+dz*dz+dn*dn;
    }
    inline FOFloat diff2(FOFparticle *q) {
        FOFloat dx = q->x-x;
        FOFloat dy = q->y-y;
        FOFloat dz = q->z-z;
	return dx*dx+dy*dy+dz*dz;
    }
    inline FOFloat diff2(FOFparticle *q, FOFparticle *offset) {
	// We add offset to q before computing the distance
        FOFloat dx = q->x+offset->x-x;
        FOFloat dy = q->y+offset->y-y;
        FOFloat dz = q->z+offset->z-z;
	return dx*dx+dy*dy+dz*dz;
    }
    inline posstruct FOF_to_pos(){
        posstruct pos; 
        double inv_FOF_RESCALE = 1.0/ FOF_RESCALE;
        pos.x = x * inv_FOF_RESCALE; 
        pos.y = y * inv_FOF_RESCALE;
        pos.z = z * inv_FOF_RESCALE;
        return pos;
    }
    // Hopefully the compiler is smart enough to convert the following to SSE!
    // If not, we could include the SSE primatives...
    inline void min(FOFparticle &q) {
	x = std::min(x,q.x); y = std::min(y,q.y);
	z = std::min(z,q.z); fi.w = std::min(fi.w,q.fi.w);
    }
    inline void max(FOFparticle &q) {
	x = std::max(x,q.x); y = std::max(y,q.y);
	z = std::max(z,q.z); fi.w = std::max(fi.w,q.fi.w);
    }
    inline void mult(FOFloat f) {
        x *= f; y *= f; z *= f; fi.w *= f;
    }
    
    // Provide a sort operand based on the index
    // bool operator< (const FOFparticle& b) const { return (n<b.n); }
    bool operator< (const FOFparticle& b) const { return (fi.n<b.fi.n); }
};

#ifdef AVX512DIRECT  // change this to AVX512FOF if we ever decide it's faster
#include "avx512_calls.h"
// This implementation is slower than AVX; AVX-512 is missing hadd
inline void diff2avx512_4(float *r, float *p, float *a) {
	// 4 copies of the primary float4
	AVX512_FLOATS primary = _mm512_broadcast_f32x4(*(__m128*)p);
	// Load the 4 secondary float4s
	AVX512_FLOATS abcd = AVX512_LOAD_FLOATS_ALIGNED(a);
	abcd = AVX512_SUBTRACT_FLOATS(abcd,primary);
	abcd = AVX512_SQUARE_FLOAT(abcd);

    __m256 ab = _mm256_hadd_ps(_mm512_castps512_ps256(abcd), _mm512_extractf32x8_ps(abcd, 1));
    ab = _mm256_hadd_ps(ab, ab);
    __m128 dists = _mm_unpacklo_ps(_mm256_castps256_ps128(ab),
                                  _mm256_extractf128_ps(ab,1));
    _mm_store_ps(r,dists);

	// This is obviously awful
	/*r[0] = AVX512_MASK_HORIZONTAL_SUM_FLOATS((AVX512_MASK)0b1111, abcd);
	r[1] = AVX512_MASK_HORIZONTAL_SUM_FLOATS((AVX512_MASK)0b11110000, abcd);
	r[2] = AVX512_MASK_HORIZONTAL_SUM_FLOATS((AVX512_MASK)0b111100000000, abcd);
	r[3] = AVX512_MASK_HORIZONTAL_SUM_FLOATS((AVX512_MASK)0b1111000000000000, abcd);*/
}

// Could make a masked version for < 16 particles
// Still marginally slower than the AVX version
inline void diff2avx512_16(float *r, float *p, float *a) {
	// 16 copies of the primary float4
	AVX512_FLOATS px = AVX512_SET_FLOAT(p[0]);
	AVX512_FLOATS py = AVX512_SET_FLOAT(p[1]);
	AVX512_FLOATS pz = AVX512_SET_FLOAT(p[2]);

	// Load the 16 secondary FLOAT4s as List3s
	AVX512_FLOATS ax, ay, az;
	for(int j = 0; j < AVX512_NVEC; j++){
		ax[j] = a[4*j];
		ay[j] = a[4*j + 1];
		az[j] = a[4*j + 2];
	}
	
	AVX512_FLOATS dx = AVX512_SUBTRACT_FLOATS(ax, px);
	AVX512_FLOATS dy = AVX512_SUBTRACT_FLOATS(ay, py);
	AVX512_FLOATS dz = AVX512_SUBTRACT_FLOATS(az, pz);

	AVX512_FLOATS dr2 = AVX512_SQUARE_FLOAT(dx);
    dr2 = AVX512_FMA_ADD_FLOATS(dy, dy, dr2);
    dr2 = AVX512_FMA_ADD_FLOATS(dz, dz, dr2);

    _mm512_store_ps(r, dr2);
}

// Try some manual loop unrolling
inline void diff2avx512_32(float *r, float *p, float *a) {
	// 16 copies of the primary float4
	AVX512_FLOATS px = AVX512_SET_FLOAT(p[0]);
	AVX512_FLOATS py = AVX512_SET_FLOAT(p[1]);
	AVX512_FLOATS pz = AVX512_SET_FLOAT(p[2]);

	// Load the 32 secondary FLOAT4s as List3s
	AVX512_FLOATS ax, ay, az;
	AVX512_FLOATS ax2, ay2, az2;
	for(int j = 0; j < AVX512_NVEC; j++){
		ax[j] = a[4*j];
		ay[j] = a[4*j + 1];
		az[j] = a[4*j + 2];
	}
	
	AVX512_FLOATS dx = AVX512_SUBTRACT_FLOATS(ax, px);
	AVX512_FLOATS dy = AVX512_SUBTRACT_FLOATS(ay, py);
	AVX512_FLOATS dz = AVX512_SUBTRACT_FLOATS(az, pz);

	for(int j = 0; j < AVX512_NVEC; j++){
		ax2[j] = a[4*(AVX512_NVEC + j)];
		ay2[j] = a[4*(AVX512_NVEC + j) + 1];
		az2[j] = a[4*(AVX512_NVEC + j) + 2];
	}

	AVX512_FLOATS dr2 = AVX512_SQUARE_FLOAT(dx);
    dr2 = AVX512_FMA_ADD_FLOATS(dy, dy, dr2);
    dr2 = AVX512_FMA_ADD_FLOATS(dz, dz, dr2);

	AVX512_FLOATS dx2 = AVX512_SUBTRACT_FLOATS(ax2, px);
	AVX512_FLOATS dy2 = AVX512_SUBTRACT_FLOATS(ay2, py);
	AVX512_FLOATS dz2 = AVX512_SUBTRACT_FLOATS(az2, pz);

	_mm512_store_ps(r, dr2);

    AVX512_FLOATS dr22 = AVX512_SQUARE_FLOAT(dx2);
    dr22 = AVX512_FMA_ADD_FLOATS(dy2, dy2, dr22);
    dr22 = AVX512_FMA_ADD_FLOATS(dz2, dz2, dr22);

    _mm512_store_ps(r + AVX512_NVEC, dr22);
}
#endif // AVX512DIRECT


#ifndef AVXFOF

/// This takes 4 FLOAT4's as a[0-3] and computes |p-a|^2
/// There's no additional assumption of alignment.
inline void diff2by4(FOFloat *r, FOFparticle *p, FOFparticle *a) {
    r[0] = (*p).diff2v(a);
    r[1] = (*p).diff2v(a+1);
    r[2] = (*p).diff2v(a+2);
    r[3] = (*p).diff2v(a+3);
}

#else

/** This takes 4 float4's and computes |p-a|^2 by AVX

p and r are required to be 16-byte aligned.
a is required to be 32-byte aligned (unless we want to use loadu),
so this is alignment to two float4's.
We may not want to cross the 64-byte cache line.
*/
inline void diff2avx4(float *r, float *p, float *a) {
    __m256 primary = _mm256_broadcast_ps((__m128 *)p);
    __m256 ab = _mm256_load_ps(a);
    __m256 cd = _mm256_load_ps(a+8);
    ab = _mm256_sub_ps(ab,primary);
    cd = _mm256_sub_ps(cd,primary);
    ab = _mm256_mul_ps(ab,ab);
    cd = _mm256_mul_ps(cd,cd);
    ab = _mm256_hadd_ps(ab,cd);
    ab = _mm256_hadd_ps(ab,ab);
    __m128 abcd = _mm_unpacklo_ps(_mm256_castps256_ps128(ab),
                                  _mm256_extractf128_ps(ab,1));
    _mm_store_ps(r,abcd);
    return;
}



// Checking whether this version is faster, but current tests say slower.
inline void diff2avx8(float *r, float *p, float *a) {
    // p and r are required to be 16-byte aligned.
    // a is required to be 32-byte aligned (unless we want to use loadu),
    // so this is alignment to two float4's.
    // We may not want to cross the 64-byte cache line.
    // This takes 8 float4's and computes |p-a|^2
    __m256 primary = _mm256_broadcast_ps((__m128 *)p);
    __m256 ab = _mm256_load_ps(a);
    __m256 cd = _mm256_load_ps(a+8);
    __m256 ef = _mm256_load_ps(a+16);
    __m256 gh = _mm256_load_ps(a+24);
    ab = _mm256_sub_ps(ab,primary);
    cd = _mm256_sub_ps(cd,primary);
    ef = _mm256_sub_ps(ef,primary);
    gh = _mm256_sub_ps(gh,primary);
    ab = _mm256_mul_ps(ab,ab);
    cd = _mm256_mul_ps(cd,cd);
    ef = _mm256_mul_ps(ef,ef);
    gh = _mm256_mul_ps(gh,gh);
    ab = _mm256_hadd_ps(ab,cd);
    ef = _mm256_hadd_ps(ef,gh);
    ab = _mm256_hadd_ps(ab,ab);
    ef = _mm256_hadd_ps(ef,ef);
    __m128 abcd = _mm_unpacklo_ps(_mm256_castps256_ps128(ab),
                                  _mm256_extractf128_ps(ab,1));
    __m128 efgh = _mm_unpacklo_ps(_mm256_castps256_ps128(ef),
                                  _mm256_extractf128_ps(ef,1));
    _mm_store_ps(r,abcd);
    _mm_store_ps(r+4,efgh);
    return;
}

#endif

/** Return a pointer to a zero-indexed list of floats containing the 
square distances between primary and list[0,nlist).

Note that this will blindly overshoot the end of the list by a multiple
of 4; that's considered ok, so d2buf must have space
If AVX is used, then d2buf is 16-byte aligned space that we supply for the output,
but the returned pointer may not point to the start of it,
since we have to attend to alignment.
*/

FOFloat *compute_d2(FOFparticle *primary, FOFparticle *list, int nlist, FOFloat *d2buf, long long &numdists) {
    FOFloat *d2use;
    numdists += nlist;

    #ifndef AVXFOF
	d2use = d2buf;
	for (int a=0; a<nlist; a+=4)
	    diff2by4(d2use+a, primary, list+a);
    #else
    // This depends on the registration of unassigned.
    // AVX requires the float4 list to be 32-byte aligned
    // We assume that p is 32-byte aligned.
    // if (((list-p)&1)==0) {
    if (((uintptr_t)list)%32==0) {
	// We have an even registration
	d2use = d2buf;
	//#pragma unroll
	for (int a=0; a<nlist; a+=4)
		    diff2avx4(d2use+a, (float *)(primary), (float *)(list+a));
		    //d2use[a] = primary->diff2(list+a);  // with ICC, this is actually the same speed
    } else {
	// We have an odd registration.  Do the first object special.
	d2use = d2buf+3;
	//? d2use[0] = primary->diff2(list);
	d2use[0] = primary->diff2v(list);
	//#pragma unroll
	for (int a=1; a<nlist; a+=4)
		    diff2avx4(d2use+a, (float *)(primary), (float *)(list+a));
		    //d2use[a] = primary->diff2(list+a);
    }
    #endif
    return d2use;
}


/** This is the class to contain one FOF group, including the 
particle location and the Bounding Boxes.

Assume that FOFgroup is divisible by 16 and that the __m128 
elements will be 16-byte aligned for SSE.
*/

class alignas(16) FOFgroup {
  public:
    FOFparticle BBmin, BBmax;
    FOFloat halo_radius2; // The threshold distance squared for each halo center 
    int center_particle;  // The index of the central particle
    int start, n;   // Starting index and Number of particles
    uint8_t tmp[CACHE_LINE_SIZE-8-sizeof(FOFparticle)];     // For alignment padding
    // During calculation, these are in FOF units, but at the
	// end we restore them to input code units
    FOFgroup(int _start, int _n, FOFloat _halo_radius2, int _center) { start = _start; n = _n; halo_radius2 = _halo_radius2; center_particle = _center;}
    FOFgroup(int _start, int _n, FOFparticle _BBmin, FOFparticle _BBmax) { 
        BBmin = _BBmin; BBmax = _BBmax;
    	start = _start; n = _n; 
        halo_radius2=0.0; center_particle = 0;
    }
    // We're setting up the comparand to sort in *descending* order
    bool operator< (const FOFgroup& c) const { return (c.n<n); }
};



/** This is the main class to perform the FOF.  It establishes
a single-threaded workspace that can be invoked many times.
*/

class FOFcell {
  public:
    // These lists provide some work space, consistent between 
    FOFparticle *p;	///< The copied version of the particles, to permute
    FOFloat *d2buffer;	///< A matching array of distances
    FOFloat *d2_active;	///< Matches an array name in SO; not used
    accstruct *permutebuf;   ///< Space for doing permutations
    int *index;		///< The new indices in the original order
    int np;		///< The current size of the group
  
    int maxsize;	///< The maximum size of these lists

    // Primary outputs, aside from the re-ordering of the input particles
    FOFgroup *groups;	///< The multiplet ranges
    int ngroups;	///< The number of groups
    long long numdists;		///< Number of distances computed
    long long numcenters;	///< The number of partitions used
    long long numsorts; 	///< The number of sorts (not used in FOF)
    int nmultiplets; 	///< Number of particles in multiplets
    int nsinglet_boundary;  ///< Plus number of particles near the edge
    		// So boundary singlets are in [nmultiplets, nsinglet_boundary)

    // Some internal timings (not all are implemented, for speed reasons)
    FOFTimer time_partition;
    FOFTimer time_permute;
    FOFTimer time_fof;
    FOFTimer time_d2;
    FOFTimer time_close;
    FOFTimer time_copy;
    FOFTimer time_total;

    // The linking lengths that we've been setup to use.
    FOFloat linking_length;   ///< The original entry, in code units
    FOFloat boundary;   ///< The location of the boundary, in FOF units
    FOFloat b;		///< The linking length, in FOF units
    FOFloat b2;		///< The linking length squared, in FOF units

    uint8_t pad[CACHE_LINE_SIZE];	// Just to ensure that we don't have cache line fights 
    		// within an array of FOFCell's.

    inline void reset(int _size) {
	int ret;
	np = ngroups = 0;
	if (_size+16<maxsize) return;    // Do nothing if we have enough space
	maxsize = _size*1.9+16;   // Oversize to limit the number of re-allocations
	assertf(maxsize<8e6, "Maxsize {:d} is too large\n", maxsize);	
	    // This is required because FOFparticle is storing the 
	    // indices in a float.  We could in principle work around
	    // this, by using additional bits in a negative exponent.
	    // But current Abacus plans do not produce cells this big.
	// fmt::print("Allocating FOFCell to maxsize = {:d}\n", maxsize);
        if (p!=NULL) free(p);
	ret = posix_memalign((void **)&p, CACHE_LINE_SIZE, sizeof(FOFparticle)*maxsize);  assert(ret == 0);
	memset(p, 0, sizeof(FOFparticle)*maxsize);

        if (permutebuf!=NULL) free(permutebuf);
	ret = posix_memalign((void **)&permutebuf, CACHE_LINE_SIZE, sizeof(accstruct)*maxsize);  assert(ret == 0);

        if (d2buffer!=NULL) free(d2buffer);
	ret = posix_memalign((void **)&d2buffer, CACHE_LINE_SIZE, sizeof(FOFloat)*maxsize);  assert(ret == 0);
	if (d2_active!=NULL) free(d2_active);
	ret = posix_memalign((void **)&d2_active, CACHE_LINE_SIZE, sizeof(FOFloat)*maxsize);  assert(ret == 0);
        if (groups!=NULL) free(groups);
	ret = posix_memalign((void **)&groups, CACHE_LINE_SIZE, sizeof(FOFgroup)*maxsize);  assert(ret == 0);
        if (index!=NULL) free(index);
	ret = posix_memalign((void **)&index, CACHE_LINE_SIZE, sizeof(int)*maxsize);  assert(ret == 0);
    }

    /// We have a null constructor, since we'll define one for each thread.
    /// Be sure to call setup() for each instance.
    FOFcell() { 
	p = NULL;
	d2buffer = NULL;
	d2_active = NULL;
	permutebuf = NULL;
	groups = NULL;
	index = NULL;
	maxsize = -1;
	numdists = 0;
	reset(32768);
	// We will have terrible bugs if these aren't true!
	#ifdef AVXFOF
	    assertf(sizeof(FOFparticle)==16, "FOFparticle is illegal size!  Can't mix AVXFOF and double precision FOFloat");
	#endif
	assertf(sizeof(FOFgroup)%16==0, "FOFgroup is illegal size!");  
        return; 
    }   

    void destroy() {
        if (p!=NULL) free(p); p = NULL;
        if (d2buffer!=NULL) free(d2buffer); d2buffer = NULL;
	if (d2_active!=NULL) free(d2_active); d2_active = NULL;
        if (permutebuf!=NULL) free(permutebuf); permutebuf = NULL;
        if (groups!=NULL) free(groups); groups = NULL;
        if (index!=NULL) free(index); index = NULL;
    }
    ~FOFcell() {
	destroy();
    }

    /// _b and _cellsize should be supplied in posstruct units.
    /// Positions are assumed to be cell-centered.
    /// We will change the length units to be single-precision,
    /// with the lower cell edge being at 0 and the units of linking lengths
    void setup (FOFloat _b, FOFloat _boundary) {
	linking_length = _b;
	boundary = _boundary*FOF_RESCALE;
	b = _b*FOF_RESCALE;
	b2 = b*b;
	// invb = 1.0/_b; 		// The length rescaling factor
        // b = b2 = 1.0;		// We've scaled to unit linking length
	return;
    }

    // ======================  Routine to permute at the end ==============

    inline bool near_boundary(FOFparticle p) {
	// printf("%f %f %f\n", p.x/boundary, p.y/boundary, p.z/boundary);
        if (fabs(p.x)>boundary || fabs(p.y)>boundary || fabs(p.z)>boundary) return true;
		else return false;
    }

        /** This performs the permutation of the original arrays.
	
	After we have the FOFparticles sorted into groups, we need to 
	copy the pos/vel/aux into that order.  But we also want to 
	collect the singlets at the end.  Recall that only multiplets
	have FOFgroups at all!
	
	We have the old indices in the new order, but it's more 
	helpful to have the new indices in the old order, 
	so that we can scatter rather than gather.
	*/

    inline void permute(posstruct *pos, velstruct *vel, auxstruct *aux, accstruct *acc, int n) {
	// We also need to undo the scaling of the positions
	FOFloat inv_rescale = 1.0/(FOF_RESCALE);

	// Count the number of particles in the multiplet
	nmultiplets = 0;
	for (int g=0; g<ngroups; g++) nmultiplets += groups[g].n;
	nsinglet_boundary = nmultiplets;

	int j = 0;    // The current particle number
	int g = 0;    // The current multiplet
	int s = np-1;	// The current interior singlet location
	int m = 0;	// The current multiplet location
	int nextstart;
	while (j<np) {
	    // Process any singlets before the next group
	    // If there are no more groups, we just want to 
	    // process to the end.
	    if (g<ngroups) nextstart = groups[g].start;
		else nextstart = np;
	    for (;j<nextstart;j++) {
		// We imagine that this particle was actually at 
		// location s in our list.
		if (near_boundary(p[j])) {
		     index[p[j].index()] = nsinglet_boundary++;
		}
		else index[p[j].index()] = s--;
	    }
	    if (g==ngroups) break;

	    // Process the next group
	    nextstart += groups[g].n;
	    groups[g].start = m;    // Change the starting location!
	    for (;j<nextstart;j++,m++) {
		// We imagine that this particle was actually at 
		// location m in our list.
		index[p[j].index()] = m;
	    }
	    // AVXAVX
	    groups[g].BBmin.mult(inv_rescale);
	    groups[g].BBmax.mult(inv_rescale);
	    g++;
	}

	// Now we need to apply this
	// This code used to re-use p for permutebuf, but that wasn't enough
	// space when accstruct is double4 and FOFparticle is only float4.
	// Not only could it overflow the buffer, but it was wrecking the 
	// AVX code by writing illegally large values into ghost particles.
	{
	posstruct *tmp = (posstruct *)permutebuf;
	for (j=0;j<n;j++) tmp[index[j]] = pos[j];
	memcpy(pos, tmp, n*sizeof(posstruct));
	}
	{
	velstruct *tmp = (velstruct *)permutebuf;
	for (j=0;j<n;j++) tmp[index[j]] = vel[j];
	memcpy(vel, tmp, n*sizeof(velstruct));
	}
	{
	auxstruct *tmp = (auxstruct *)permutebuf;
	for (j=0;j<n;j++) tmp[index[j]] = aux[j];
	memcpy(aux, tmp, n*sizeof(auxstruct));
	}
	if(acc != NULL){
	accstruct *tmp = (accstruct *)permutebuf;
	for (j=0;j<n;j++) tmp[index[j]] = acc[j];
	memcpy(acc, tmp, n*sizeof(accstruct));
	}

    }

    // ============  Some helper routines for the FOF work ============

    inline void add_to_pending(FOFparticle *&b, FOFparticle *a, 
    		FOFparticle &BBmin, FOFparticle &BBmax) {
	// We're going to add 'a' to our pending queue, 
	// moving it to position b and incrementing b.
	FOFparticle new_member = *a;
	BBmin.min(new_member);
	BBmax.max(new_member);
	*a = *b;
	*b = new_member;
	b++;
	return;
    }

    inline void add_to_pending(FOFparticle *&c, FOFparticle *&b, FOFparticle *a, 
    		FOFparticle &BBmin, FOFparticle &BBmax) {
	// We're going to add 'a' to our pending queue, 
	// moving it to position c and incrementing c.
	// Via position b, which is incremented.
	FOFparticle new_member = *a;
	BBmin.min(new_member);
	BBmax.max(new_member);
	*a = *b;
	*b = *c;
	*c = new_member;
	b++; c++;
	return;
    }

    inline void add_to_pending(FOFparticle *&d, FOFparticle *&c, FOFparticle *&b, FOFparticle *a, 
    		FOFparticle &BBmin, FOFparticle &BBmax) {
	// We're going to add 'a' to our pending queue, 
	// moving it to position d and incrementing d.
	// Via position b and c, which are incremented.
	FOFparticle new_member = *a;
	BBmin.min(new_member);
	BBmax.max(new_member);
	*a = *b;
	*b = *c;
	*c = *d;
	*d = new_member;
	b++; c++; d++;
	return;
    }

    inline void add_to_pending(FOFparticle *&e, FOFparticle *&d, FOFparticle *&c, 
    	FOFparticle *&b, FOFparticle *a, FOFparticle &BBmin, FOFparticle &BBmax) {
	// We're going to add 'a' to our pending queue, 
	// moving it to position e and incrementing e.
	// Via position b, c, and d, which are incremented.
	FOFparticle new_member = *a;
	BBmin.min(new_member);
	BBmax.max(new_member);
	*a = *b;
	*b = *c;
	*c = *d;
	*d = *e;
	*e = new_member;
	b++; c++; d++; e++;
	return;
    }


    // ====================  Alternative modified N^2 algorithm ==========

	/** Do the FOF in the Cell, without using any list partitioning.

	This tries to accelerate the N^2 calculation by partitioning 
	the unassigned particle list as we sweep through on one primary,
	so that all of the particles remaining in the pending list can
	be searched more quickly.  

	In practice, we're going collect particles in a 'core' out to
	distance R, and then a 'skin' out to R+b.  The pending queue
	will separate particles in the core, including newly added 
	particles, so that we can complete the search on all of these
	primaries.  

	Setting the distance R will require a heuristic.  We certainly
	want R>=b, but probably we want to tie to the number of particles
	in the cell.  When there are few particles, we want a big radius,
	so that if we complete a group with yet more unassigned particles
	in the core, then we could continue on.
	*/

    void quickFOF() {
    // TODO: is this safe?  do we rely on other side-effects from this routine?
    if (b == 0)
        return;

	time_fof.Start();

	FOFloat core = 2.5;     // Units of linking lengths to start
	int core_size = 5;     // The number of particles found in the core.

	FOFparticle *start; 
	FOFparticle *primary;
	FOFparticle *pending_notcore;
	FOFparticle *unassigned_core;
	FOFparticle *unassigned_skin;
	FOFparticle *unassigned_far;
	FOFparticle *end;
	// The current group is [start,unassigned_core)
	// Core pending is [primary+1,pending_notcore)
	// The unassigned list is [unassigned_core,end)
	// Completed singlets will get moved to the end.

	primary = p;
	end = p+np;
	start = primary;
	unassigned_core = primary+1;
	unassigned_skin = primary;   // Set this to force an initial partitioning

	FOFparticle BBmin, BBmax;
	// AVXAVX
	BBmin = *primary; BBmax = *primary;

	while (primary<end) {
	    // We have a primary to search.  

	    if (primary==start && primary<unassigned_skin) {
		// This is a new group seed, but it is within
		// the core of the previous partition.  Therefore, we
		// can search it efficiently.

		pending_notcore = unassigned_core;
		numcenters++;
		FOFloat *d2use = compute_d2(primary, unassigned_core, 
	    		unassigned_far-unassigned_core, d2buffer, numdists);
		for (FOFparticle *a = unassigned_core; a<unassigned_skin; a++, d2use++) {
		    // numdists++;
		    if (*d2use<b2) 
			// Add this to the pending core
			add_to_pending(pending_notcore, unassigned_core, a, BBmin, BBmax);
		}
		for (FOFparticle *a = unassigned_skin; a<unassigned_far; a++, d2use++) {
		    // numdists++;
		    if (*d2use<b2) 
			// Add this to the pending notcore
			add_to_pending(unassigned_core, unassigned_skin, a, BBmin, BBmax);
		}
	    } else {
		// We need to repartition.
		// Use a heuristic to try to get a few new particles in the
		// core.  We want some, so that we can reuse the partition.
		// TODO: Need to tune this more!
		if (core_size>8 && core>2.5) core /= 1.2;
		if (core_size<3) core *= 1.44;
		FOFloat core_radius = core*b;    
		FOFloat skin_radius = core_radius+b;
		core_radius *= core_radius;
		skin_radius *= skin_radius;

		// We need to compute the distance to all remaining particles
		numcenters++;
		FOFloat *d2use = compute_d2(primary, primary+1,
				    end-(primary+1), d2buffer, numdists);

		// Sweep the partition list to partition core vs notcore.
		pending_notcore = unassigned_core;
		FOFloat *low = d2use;
		FOFloat *high = d2use+(int)(unassigned_core-(primary+1));
		for (FOFparticle *a = primary+1; a<pending_notcore; a++, low++) {
		    // numdists++;
		    if (*low>core_radius) {
			do {
			    pending_notcore--; high--;
			    if (*high<core_radius) {
				std::swap(*a, *pending_notcore);
				break;
			    }
			} while (a<pending_notcore);
		    }
		}

		// Sweep the unassigned list
		unassigned_far = unassigned_skin = unassigned_core;
		d2use += unassigned_core-(primary+1);
		for (FOFparticle *a = unassigned_core; a<end;a++,d2use++) {
		    // numdists++;
		    if (*d2use<b2) {
			// We found a neighbor!  Add to pending core.
			add_to_pending(pending_notcore, unassigned_core, unassigned_skin, unassigned_far, a, BBmin, BBmax);
		    } else if (*d2use<core_radius) {
			// Found an unassigned core particle
			FOFparticle tmp = *a;
			*a = *unassigned_far;
			*(unassigned_far++) = *unassigned_skin;
			*(unassigned_skin++) = tmp;
		    } else if (*d2use<skin_radius) {
			// Found an unassigned skin particle
			FOFparticle tmp = *a;
			*a = *unassigned_far;
			*(unassigned_far++) = tmp;
		    }
		}

		core_size = unassigned_skin-unassigned_core;
		// We're saving this for the next sweep.
	    }

	    // Advance primary, searching core+skin, to exhaust pending_core
	    primary++;
	    while (primary<pending_notcore) {
		// We don't increment numcenters here, because we're 
		// not touching the bulk of the particles.
		FOFloat *d2use = compute_d2(primary, unassigned_core, 
	    		unassigned_far-unassigned_core, d2buffer, numdists);
		// Now consider this particle
		for (FOFparticle *a = unassigned_core; a<unassigned_skin; a++, d2use++) {
		    // numdists++;
		    if (*d2use<b2) 
			// Add this to the pending core
			add_to_pending(pending_notcore, unassigned_core, a, BBmin, BBmax);
		}
		for (FOFparticle *a = unassigned_skin; a<unassigned_far; a++, d2use++) {
		    // numdists++;
		    if (*d2use<b2) 
			// Add this to the pending notcore
			add_to_pending(unassigned_core, unassigned_skin, a, BBmin, BBmax);
		}
		primary++;
	    }
	    // At this point, primary has walked forward to pending_notcore.
	    // If there are pending notcore particles, we'll proceed to those.

	    // Check if the group is now closed (no more pending particles)
	    if (primary==unassigned_core) {
		if (unassigned_core != start+1) {
		    // We've closed a group and it's a multiplet
		    groups[ngroups++] = FOFgroup(start-p, unassigned_core-start, BBmin, BBmax);
		}
		start = primary; unassigned_core = start+1;
		// AVXAVX
		BBmin = *primary; BBmax = *primary;
	    }
	    // If the group is open, then we'll have 
	    // primary=pending_notcore<unassigned_core.
	    // We want to form a new partitioning and keep building the group.
	    // If the group is closed, then the new primary particle
	    // is unassigned.  It may be still in the core, in which
	    // case we don't need a new partition.  This is indicated by
	    // primary==start&&primary<unassigned_skin.
	}
	time_fof.Stop();
	return;
    }

    // ====================  Simple N^2 algorithms ==========

    /// Here's the basic code, using AVX for distances and do bounding boxes.

    void simpleFOFavx() {
	time_fof.Start();
	FOFparticle *start = p;    // The group start
	FOFparticle *unassigned = p+1;    // Group is [start, unassigned)
	FOFparticle *end = p+np;   // Unassigned particles are [unassigned,end)
	FOFparticle *primary = start;   // The current search point
	FOFparticle BBmin, BBmax;
	// AVXAVX
	BBmin = *primary; BBmax = *primary;

	while (primary<end) {
	    // First, we're going to compute the distances to [unassigned,end).
	    // These will start in d2use[0].

	    numcenters++;
	    FOFloat *d2use = compute_d2(primary, unassigned, end-unassigned, d2buffer, numdists);
	    for (FOFparticle *a = unassigned; a<end; a++,d2use++) {
		// numdists++;
		if (*d2use<b2) add_to_pending(unassigned, a, BBmin, BBmax);
	    }
	    primary++;
	    if (primary==unassigned) {
		// All done with this group
		if (unassigned!=start+1) { 
		    // Record the multiplet
		    groups[ngroups++] = FOFgroup(start-p, unassigned-start, BBmin, BBmax);
		}
		// Set up the next group
		start = primary; unassigned=primary+1;
		// AVXAVX
		BBmin = *primary; BBmax = *primary;
	    }
	}
	time_fof.Stop();
	return;
    }

    // This is a reference code.  Doesn't do bounding boxes.
    void simpleFOF() {
	time_fof.Start();
	FOFparticle *start = p;    // The group start
	FOFparticle *unassigned = p+1;    // Group is [start, unassigned)
	FOFparticle *end = p+np;   // Unassigned particles are [unassigned,end)
	FOFparticle *primary = start;   // The current search point

	while (primary<end) {
	    for (FOFparticle *a = unassigned; a<end; a++) {
		// numdists++;
		if (a->diff2(primary)<b2) {
		    if (unassigned!=a) std::swap(*unassigned,*a);
		    unassigned++;
		}
	    }
	    primary++;
	    if (primary==unassigned) {
		// All done with this group
		if (unassigned!=start+1) { 
		    // Record the multiplet
          //groups[ngroups++] = FOFgroup(start-p, unassigned-start);
          groups[ngroups++] = FOFgroup(start-p, unassigned-start,0., 0.); //B.H. doesn't matter-reference
		}
		// Set up the next group
		start = primary; unassigned=primary+1;
	    }
	}
	time_fof.Stop();
	return;
    }

    // ============  And here's the main driver routine

    /** The main driver routine for doing FOF on a set of particles.
    
    This takes in all of the positions, velocities, and aux for 
    a cell of size n.
    It copies pos[0,n) into the FOFparticle array.
    Then perform the FOF and permutes the input arrays
    Return the number of multiplet groups.

    If vel==NULL, then the permutation is skipped; the calling routine
    is left responsible for dealing with the internal variables of the class.
    pos[] will be unchanged.
    */

    inline int findgroups(posstruct *pos, velstruct *vel, auxstruct *aux, accstruct *acc, int n) {
	time_total.Start();
	reset(n);
	np = n;

	// Load the particles
	time_copy.Start();
	for (int j=0; j<np; j++) { 
	    p[j] = FOFparticle(pos[j],j);
	}
	time_copy.Stop();

	if (np>70)     // We get slightly different performances
	    quickFOF();
	else 
	    simpleFOFavx();

	// simpleFOF();

	// Now we want to re-order the input arrays into the multiplets+singlet order
	if (vel!=NULL) {
	    time_permute.Start();
	    permute(pos,vel,aux,acc,np);
	    time_permute.Stop();
	}
	time_total.Stop();

	return ngroups;
    }

    /// This routine does a trivial pass-through: all particles are assigned to
    /// one single group.  This is useful if two levels of the FOF hierarchy have
    /// the same linking length.
    inline int assign_to_one_group(posstruct *pos, velstruct *vel, auxstruct *aux, accstruct *acc, int n) {
        time_total.Start();
        FOFparticle BBmin, BBmax;
        reset(n);
        np = n;
        for (int j=0; j<np; j++) {
            p[j] = FOFparticle(pos[j],j);
            if (j==0) {
                BBmin = *p; BBmax = *p;
            } else {
                BBmin.min(p[j]);
                BBmax.max(p[j]);
            }
        }
        // Now mark one group
        groups[ngroups++] = FOFgroup(0, n, BBmin, BBmax);
        time_total.Stop();
        return ngroups;
    }

};    // End FOFcell









#ifdef TEST


// ========================  COMPARISON CODE ==================


inline void particle_swap(int i, int j, posstruct *pos, velstruct *vel, auxstruct *aux) {
    std::swap(pos[i],pos[j]);
    std::swap(vel[i],vel[j]);
    std::swap(aux[i],aux[j]);
    return;
}



int FOF(posstruct *pos, velstruct *vel, auxstruct *aux, int N, FOFloat b2) {
    int ngroup = 0;     // The number of groups closed so far, zero indexed

    int idxstart = 0;            // The first particle in the current group
    int idxend = 1;	         // The group is [idxstart,idxend)
    int idxprimary = idxstart;   // The next primary particle to search

    while (idxprimary<N) {
	// We have unmatched primary particles remaining.
	// Test the primary against a range [idxend,N)
	posstruct primary = pos[idxprimary];
	for (int j=idxend; j<N; j++) {
	    posstruct res = pos[j]-primary;
	    if (res.norm2()<b2) {
		// We found a new friend
		if (idxend!=j) particle_swap(idxend,j,pos,vel,aux);   
		idxend++;
	    }
	}
	// All done with that primary, move to the next one.
	// But if the group has consumed all remaining particles, then can shortcut to the end.
	if (idxend<N) idxprimary++;
	    else idxprimary = N;     

	if (idxprimary==idxend) {
	    // We've reached the end of the group, so it's closed.
	    // Group is the range [idxstart,idxend), so multiplicity is idxend-idxstart.
	    int thissize = idxend-idxstart;
	    if (thissize>1) ngroup++;  // Only count N>1 groups
	    // ngroup++;  
	    idxstart = idxend; idxend++;   // Ready for the next one
	}
    }
    return ngroup;
}


// ========================  DRIVER CODE ==================



int main() {
    stdlog_threshold_global = 2;
    stdlog.open("/tmp/foftest.log");
    STDLOG_TIMESTAMP;
    int cellsize = 80000;
    int Ncell = 1000, Nindep = 10;
    // int Ncell = 1, Nindep = 1;
    FOFloat b = 0.03;
    srand48(124);
    posstruct *pos = (posstruct *)malloc(sizeof(posstruct)*cellsize*Ncell);
    velstruct *vel = (velstruct *)malloc(sizeof(velstruct)*cellsize*Ncell);
    accstruct *acc = (accstruct *)malloc(sizeof(accstruct)*cellsize*Ncell);
    auxstruct *aux = (auxstruct *)malloc(sizeof(auxstruct)*cellsize*Ncell);
    for (int j=0; j<cellsize*Nindep; j++) {
	// This causes the particles to be clustered (cubicly) near the origin
	FOFloat rescale = drand48();
        pos[j].x = (drand48()-0.5)*rescale;
        pos[j].y = (drand48()-0.5)*rescale;
        pos[j].z = (drand48()-0.5)*rescale;
	vel[j].x = vel[j].y = vel[j].z = 0.0;
	acc[j].x = acc[j].y = acc[j].z = 0.0;
	aux[j] = 0;
    }
    #pragma omp parallel for schedule(static)
    for (int j=cellsize*Nindep; j<cellsize*Ncell; j++) {
	int k = j%(cellsize*Nindep);
        pos[j] = pos[k];
	vel[j].x = vel[j].y = vel[j].z = 0.0;
	acc[j].x = acc[j].y = acc[j].z = 0.0;
	aux[j] = 0;
    }

    // Now to setup the FOF and then execute it
    
    STDLOG(1,"Starting FoF\n"); fflush(NULL);
    FOFcell doFOF[64];
    #pragma omp parallel for schedule(static)
    for (int g=0; g<omp_get_num_threads(); g++) {
	doFOF[g].setup(b, 0.5-b);
	if (g==0) STDLOG(3,"Using {:d} threads\n", omp_get_num_threads());
    }

    STimer FOFtime;
    FOFtime.Start();
    int ngroup=0;
    FOFloat b2 = b*b;
    #pragma omp parallel for schedule(static) reduction(+:ngroup)
    for (int j=0; j<Ncell; j++) {
        // ngroup += FOF(pos+j*cellsize, vel+j*cellsize, aux+j*cellsize, acc+j*cellsize, cellsize, b2);
        ngroup += doFOF[omp_get_thread_num()].findgroups(pos+j*cellsize, vel+j*cellsize, aux+j*cellsize, acc+j*cellsize, cellsize);
	// printf("\n");
    }
    FOFtime.Stop();


    // Output the timing!

    // For timing, the quantity of interest is probably pairs per second.
    float Npair = (float)Ncell*cellsize*cellsize;

    STDLOG(3,"Found {:f} groups per cell, on average\n", (float)ngroup/Ncell);
    STDLOG(3,"Used {:d} pairwise calculations\n", doFOF[0].numdists);
    STDLOG(3,"Time to do {:d} particles in {:d} cells: {:f} ({:f} Mp/sec, {:f} Gpair/sec)\n", 
    	cellsize, Ncell, FOFtime.Elapsed(), cellsize*Ncell/FOFtime.Elapsed()/1e6,
	doFOF[0].numdists/doFOF[0].time_total.Elapsed()/1e9);
    STDLOG(3,"Copy:      {:f}\n", doFOF[0].time_copy.Elapsed());
    STDLOG(3,"   Distances: {:f}\n", doFOF[0].time_d2.Elapsed());
    STDLOG(3,"   Partition: {:f}\n", doFOF[0].time_partition.Elapsed());
    STDLOG(3,"   Close Groups: {:f}\n", doFOF[0].time_close.Elapsed());
    STDLOG(3,"FOF:       {:f}\n", doFOF[0].time_fof.Elapsed());
    STDLOG(3,"Permute:   {:f}\n", doFOF[0].time_permute.Elapsed());
    STDLOG(3,"Total:     {:f}\n", doFOF[0].time_total.Elapsed());
    return 0;
}

#endif
