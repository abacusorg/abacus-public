/*

This class contains implementations of the ICFile abstract base
class from IC_base.h.  This is where we define our IC formats.

New IC formats must implement at least the unpack() and Nparticles()
methods, and register themselves in ICFile::FromFormat() at the bottom
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
 * Each ICFile class must define two routines: unpack() and Nparticles().
 * unpack() pushes the particles to the insert list or (if given) velslab.
 * Nparticles() returns the particle count.
 */

/****************
 * ICFile
 ****************/

// Unpack the particles in the arena and push them to the Insert List
uint64 ICFile::unpack_to_IL(double convert_pos, double convert_vel){
    return unpack(NULL, convert_pos, convert_vel);
}

// Unpack the particle velocities in the arena and store them in velslab
uint64 ICFile::unpack_to_velslab(velstruct *velslab, double convert_pos, double convert_vel){
    assertf(velslab != NULL, "Passed NULL velslab in unpack_to_velslab?\n");
    return unpack(velslab, convert_pos, convert_vel);
}

void ICFile::read_vel_nonblocking(){
    read_nonblocking(1);
}

// Start a read via SB
void ICFile::read_nonblocking(int vel){
    // Base implementation doesn't distinguish
    // between regular reads and vel-only reads
    SB->LoadArenaNonBlocking(ICSlab, slab);
}

int ICFile::check_vel_read_done(){
    return check_read_done(1);
}

// If a read is in progress, return 0
int ICFile::check_read_done(int vel){
    if(!SB->IsIOCompleted(ICSlab, slab)){
        if(SB->IsSlabPresent(ICSlab, slab))
            Dependency::NotifySpinning(WAITING_FOR_IO);
        return 0;
    }
    return 1;
}

ICFile::ICFile(int _slab){
    slab = _slab;

    // Ideally we'd call Nparticles(), but virtual calls are not allowed inside constructors.
    // So we punt to the factory.
    Npart = 0;
}

ICFile::~ICFile(void) {
    next_pid += Npart;
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

// // ===================================================================
// // Simple initial conditions: lists of positions & velocities 
// // in double precision.  The positions and velocities will be 
// // translated by ICPositionRange and ICVelocity2Displacement.

// class ICFile_RVdouble: public ICFile {
// private:
//     double convert_pos, convert_vel;
//     static uint64 pid;
//     FILE *fp;
//     class ICparticle {
//     public:
//         double pos[3];
//         double vel[3];
//     };
//     int read(ICparticle *p) {
//         // Must return 0 if read is unsuccessful
//         // Return the number of particles read otherwise
//         int n = fread(p, sizeof(ICparticle), 1, fp);
//         return n;
//     }
//     int readheader() {
//         return 0;  // All's well
//     }

// public:

//     ICFile_RVdouble(char *filename) {
//         fp = fopen(filename,"r");
//         assertf(fp!=NULL,"Couldn't open IC file %s\n", filename);
//         int err = readheader(); 
//         assertf(err==0,"Error reading header of IC file %s\n", filename);
//         if (P.ICPositionRange>0) convert_pos = 1.0/P.ICPositionRange;
//         else convert_pos = 1.0/P.BoxSize;
//         if (P.ICVelocity2Displacement>-0.99) 
//         convert_vel = P.ICVelocity2Displacement*convert_pos;
//         else convert_vel = 1.0/ReadState.VelZSpace_to_kms;
//         // This gets to the unit box in redshift-space displacement.
//     }

//     ~ICFile_RVdouble(void) {
//         fclose(fp);
//     }

//     inline int getparticle(double3 *global_pos, velstruct *vel, auxstruct *aux) {
//         // This routine should return 1 if a particle has been read and 
//         // successfully loaded in the given variables.
//         // Return 0 otherwise, which will signal the end of this file.
//         ICparticle p;
//         if (read(&p)==0) return 0;
//         global_pos->x = p.pos[0];
//         global_pos->y = p.pos[1];
//         global_pos->z = p.pos[2];
//         *global_pos = (*global_pos*convert_pos);  // -double3(0.5,0.5,0.5);
//         // This maps to -0.5..+0.5
//         vel->x = p.vel[0];
//         vel->y = p.vel[1];
//         vel->z = p.vel[2];
//         *vel = *vel * convert_vel;
//         aux->clear(); aux->setpid((uint64) 0); // Set aux too.
//         aux->setpid(pid++);  // Somewhat useful for debugging
        
//         return 1;
//     }
// };
// uint64 ICFile_RVdouble::pid = 0;

// /*
// Almost the same binary format as RVdouble,
// but with additional fields for the Zeldovich grid index
// and positions are interpreted as Zel-format displacements.
// Velocity is only divided by the BoxSize, just like the displacements.
// */
// class ICFile_RVdoubleZel: public ICFile {
// public:
//     class ICparticle {
//     public:
//         unsigned short i,j,k;
//         double pos[3];
//         double vel[3];
//     };
// private:
//     double convert_pos, convert_vel;
//     FILE *fp;
//     int read(ICparticle *p) {
//         // Must return 0 if read is unsuccessful
//         // Return the number of particles read otherwise
//         int n = fread(p, sizeof(ICparticle), 1, fp);
//         return n;
//     }
//     int readheader() {
//         return 0;  // All's well
//     }

// public:

//     ICFile_RVdoubleZel(char *filename) {
//         fp = fopen(filename,"r");
//         assertf(fp!=NULL,"Couldn't open IC file %s\n", filename);
//         int err = readheader();
//         assertf(err==0,"Error reading header of IC file %s\n", filename);
//         if (P.ICPositionRange>0) convert_pos = 1.0/P.ICPositionRange;
//         else convert_pos = 1.0/P.BoxSize;
//         if (P.FlipZelDisp)
//             convert_pos *= -1;
//         if (P.ICVelocity2Displacement>-0.99) // Should always be 1 for IC from zel.cpp
//             convert_vel = P.ICVelocity2Displacement*convert_pos;
//         else convert_vel = 1.0/ReadState.VelZSpace_to_kms;
//         // This gets to the unit box in redshift-space displacement.
//     }

//     ~ICFile_RVdoubleZel(void) {
//         fclose(fp);
//     }

//     inline int getparticle(double3 *global_pos, velstruct *vel, auxstruct *aux) {
//         // This routine should return 1 if a particle has been read and 
//         // successfully loaded in the given variables.
//         // Return 0 otherwise, which will signal the end of this file.
//         ICparticle p;
//         if (read(&p)==0) return 0;
//         global_pos->x = p.pos[0];
//         global_pos->y = p.pos[1];
//         global_pos->z = p.pos[2];
//         *global_pos = (*global_pos*convert_pos);  // -double3(0.5,0.5,0.5);
//         // This maps to -0.5..+0.5
//         vel->x = p.vel[0];
//         vel->y = p.vel[1];
//         vel->z = p.vel[2];
//         *vel = *vel * convert_vel;
        
//         integer3 ijk;
//         ijk.x = p.i;
//         ijk.y = p.j;
//         ijk.z = p.k;
//         *global_pos += ZelPos(ijk);
        
//         aux->clear(); aux->setpid(ijk); // Set aux too.
//         return 1;
//     }
// };

// /*
// Almost the same binary format as RVdouble, but with a field for the PID.
// Positions are global.
// */
// class ICFile_RVdoubleTag: public ICFile {
// public:
//     class ICparticle {
//     public:
//         double pos[3];
//         double vel[3];
//         uint64 pid;
//     };
// private:
//     double convert_pos, convert_vel;
//     FILE *fp;
//     int read(ICparticle *p) {
//         // Must return 0 if read is unsuccessful
//         // Return the number of particles read otherwise
//         int n = fread(p, sizeof(ICparticle), 1, fp);
//         return n;
//     }
//     int readheader() {
//         return 0;  // All's well
//     }

// public:

//     ICFile_RVdoubleTag(char *filename) {
//         fp = fopen(filename, "r");
//         assertf(fp != NULL, "Couldn't open IC file %s\n", filename);
//         int err = readheader();
//         assertf(err == 0, "Error reading header of IC file %s\n", filename);
//         if (P.ICPositionRange>0)
//             convert_pos = 1.0/P.ICPositionRange;
//         else
//             convert_pos = 1.0/P.BoxSize;
//         if (P.ICVelocity2Displacement>-0.99) // Should always be 1 for IC from zel.cpp
//             convert_vel = P.ICVelocity2Displacement*convert_pos;
//         else
//             convert_vel = 1.0/ReadState.VelZSpace_to_kms;
//         // This gets to the unit box in redshift-space displacement.
//     }

//     ~ICFile_RVdoubleTag(void) {
//         fclose(fp);
//     }

//     inline int getparticle(double3 *global_pos, velstruct *vel, auxstruct *aux) {
//         // This routine should return 1 if a particle has been read and 
//         // successfully loaded in the given variables.
//         // Return 0 otherwise, which will signal the end of this file.
//         ICparticle p;
//         if (read(&p)==0)
//             return 0;
//         global_pos->x = p.pos[0];
//         global_pos->y = p.pos[1];
//         global_pos->z = p.pos[2];
//         *global_pos = (*global_pos*convert_pos);  // -double3(0.5,0.5,0.5);
//         // This maps to -0.5..+0.5
//         vel->x = p.vel[0];
//         vel->y = p.vel[1];
//         vel->z = p.vel[2];
//         *vel = *vel * convert_vel;

//         aux->clear();
//         aux->setpid(p.pid);

//         return 1;
//     }
// };

/*
The same format as RVdoubleZel, except in single precision
*/
class ICFile_RVZel: public ICFile {
public:
    class ICparticle {
    public:
        unsigned short i,j,k;
        float pos[3];
        float vel[3];
    };

    using ICFile::ICFile;  // inherit the constructor

private:
    // If velslab != NULL, unpack the velocities into there and do nothing else
    uint64 unpack(velstruct *velslab, double convert_pos, double convert_vel) {

        ICparticle *particles = (ICparticle *) SB->GetSlabPtr(ICSlab, slab);

        Npart = Nparticles();  // will already be initialized if created from the factory

        uint64 sumA = 0, sumB = 0;
        #pragma omp parallel for schedule(static) reduction(+:sumA,sumB)
        for(uint64 i = 0; i < Npart; i++){
            ICparticle p = particles[i];

            velstruct vel(p.vel[0], p.vel[1], p.vel[2]);
            vel *= convert_vel;

            if(velslab != NULL){
                // Fill the velslab; don't touch the pos,aux
                velslab[i] = vel;
                continue;
            }

            double3 pos(p.pos[0], p.pos[1], p.pos[2]);
            auxstruct aux;

            pos *= convert_pos;

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
            set_taggable_bits(aux, sumA, sumB);

            IL->Push(&_pos, &vel, &aux, newcell);
        }

        // Store the local reductions in the class variables
        NsubsampleA = sumA;
        NsubsampleB = sumB;
        
        return Npart;
    }

    // Called by the factory to initialize this->Npart
    uint64 Nparticles(){
        // SlabSizeBytes only works if ICSlab is already in memory, which it always is in our usage
        // This avoids an extra file stat
        uint64 b = SB->SlabSizeBytes(ICSlab, slab);
        uint64 npart = b/sizeof(ICparticle);
        
        assertf(npart*sizeof(ICparticle) == b, "Size of IC slab %s not divisible by particle size %d!\n",
            SB->ReadSlabPath(ICSlab, slab), sizeof(ICparticle));

        return npart;
    }
};

// /*
// Similar to RVZel, except the position is interpreted as an absolute position,
// not a displacement.
// */
// class ICFile_RVTag: public ICFile {
// public:
//     class ICparticle {
//     public:
//         float pos[3];
//         float vel[3];
//         uint64 tag;
//     };
// private:
//     double convert_pos, convert_vel;
//     FILE *fp;
//     int read(ICparticle *p) {
//         // Must return 0 if read is unsuccessful
//         // Return the number of particles read otherwise
//         int n = fread(p, sizeof(ICparticle), 1, fp);
//         return n;
//     }
//     int readheader() {
//         return 0;  // All's well
//     }

// public:

//     ICFile_RVTag(char *filename) {
//         fp = fopen(filename,"r");
//         assertf(fp!=NULL,"Couldn't open IC file %s\n", filename);
//         int err = readheader();
//         assertf(err==0,"Error reading header of IC file %s\n", filename);
//         if (P.ICPositionRange>0) convert_pos = 1.0/P.ICPositionRange;
//         else convert_pos = 1.0/P.BoxSize;
//         if (P.ICVelocity2Displacement>-0.99) // Should always be 1 for IC from zel.cpp
//             convert_vel = P.ICVelocity2Displacement*convert_pos;
//         else convert_vel = 1.0/ReadState.VelZSpace_to_kms;
//         // This gets to the unit box in redshift-space displacement.
//     }

//     ~ICFile_RVTag(void) {
//         fclose(fp);
//     }

//     inline int getparticle(double3 *global_pos, velstruct *vel, auxstruct *aux) {
//         // This routine should return 1 if a particle has been read and 
//         // successfully loaded in the given variables.
//         // Return 0 otherwise, which will signal the end of this file.
//         ICparticle p;
//         if (read(&p)==0) return 0;
//         global_pos->x = p.pos[0];
//         global_pos->y = p.pos[1];
//         global_pos->z = p.pos[2];
//         *global_pos = (*global_pos*convert_pos);  // -double3(0.5,0.5,0.5);
//         // This maps to -0.5..+0.5
//         vel->x = p.vel[0];
//         vel->y = p.vel[1];
//         vel->z = p.vel[2];
//         *vel = *vel * convert_vel;
        
//         aux->clear(); aux->setpid(p.tag); // Set aux too.
//         return 1;
//     }
// };


// // ===================================================================
// // Zel'dovich displacement initial conditions.
// // Here we're given the integer-based Lagrangian position and
// // the displacement.  The displacement should be in units of 
// // P.ICPositionMax, which could plausibly be 1 or BoxSize.
// // The velocities will be derived from the displacements.
// // The Lagrangian positions will be ijk = 0..CPD-1.  
// // We will store i alone, and then jk = j*CPD+k, just to avoid 64-bit math
// // and to keep the file 4-byte aligned.


// // Our LPT implementation must be coordinated with choices here.
// // Therefore, we will codify some items in some simple functions,
// // given in lpt.cpp.
// // double3 ZelPos(integer3 ijk): Returns the position in code-units 
// //             of the initial grid point.


// class ICFile_Zel: public ICFile {
// private:
//     double convert_pos, convert_vel;
//     FILE *fp;
//     class ICparticle {
//     public:
//         unsigned short i,j,k;
//         double displ[3];
//     };
//     int read(ICparticle *p) {
//         // Must return 0 if read is unsuccessful
//         // Return the number of particles read otherwise
//         int n = fread(p, sizeof(ICparticle), 1, fp);
//         //        printf("(%d,%d,%d) : (%f, %f, %f)\t \t sizeof(ICparticle) = %d\n",p->i,p->j,p->k, p->displ[0],p->displ[1],p->displ[2],sizeof(ICparticle));
//         return n;
//     }
//     int readheader() {
//         char buf[2];
//         int len = 1;
//         fread(&(buf[1]),1,1,fp);
//         do {
//             buf[0] = buf[1];
//             if(fread(&(buf[1]), 1, 1, fp)<=0) QUIT("Error reading header for zeldovich input file. Exiting\n");
//             len++;
//         } while(!(buf[0]==0x2 && buf[1]==0x2));
//         return 0;  // All's well
//     }

// public:

//     ICFile_Zel(char *filename) {
//         fp = fopen(filename,"r");

//         assertf(fp!=NULL,"Couldn't open IC file %s\n", filename);
//         //int err = readheader();
//         //assertf(err==0,"Error reading header of IC file %s\n", filename);
//         if (P.ICPositionRange>0) convert_pos = 1.0/P.ICPositionRange;
//         else convert_pos = 1.0/P.BoxSize;
//         // This will convert the displacement to the code units
//         convert_vel = WriteState.f_growth;
//         // TODO: Check that this is the correct physics.
//         // Zel'dovich velocities in redshift-space displacements units
//         // are just f*position_dispacement.
//     }

//     ~ICFile_Zel(void) {
//         fclose(fp);
//     }

//     inline int getparticle(double3 *global_pos, velstruct *vel, auxstruct *aux) {
//         // This routine should return 1 if a particle has been read and 
//         // successfully loaded in the given variables.
//         // Return 0 otherwise, which will signal the end of this file.
//         ICparticle p;
//         if (read(&p)==0) return 0;
//         global_pos->x = p.displ[0];
//         global_pos->y = p.displ[1];
//         global_pos->z = p.displ[2];
//         if (P.FlipZelDisp)
//             *global_pos *= -1;
//         *global_pos = *global_pos * convert_pos;
//         // This is the displacement in code units (unit box)
//         *vel = *global_pos * convert_vel;
//         // The velocity is related to the displacement.
//         // Now add on the initial Lagrangian position
//         integer3 ijk;
//         ijk.x = p.i;
//         ijk.y = p.j;
//         ijk.z = p.k;
//         *global_pos += ZelPos(ijk);

//         aux->clear(); aux->setpid(ijk); // Set aux too.
//         return 1;
//     }
// };


// class ICFile_Poisson: public ICFile {
// private:
//     static gsl_rng **rng; //The random number generator, NULL by default
//     int slab;  // slab index
//     int count;  // number of particles generated in this slab

//     // x-offset in global coordiantes
//     double slab_start;
//     double slab_width;
    
// public:
//     ICFile_Poisson(int _slab) : slab(_slab), count(0) {
//         maxthreads = omp_get_max_threads();

//         // Initialize once for all slabs
//         // TODO: we will not produce the same ICs for different numbers of threads.  Maybe we don't care?
//         if(rng == NULL){
//             rng = new gsl_rng*[maxthreads];
//             for(int i = 0; i < maxthreads; i++){
//                 rng[i] = gsl_rng_alloc(gsl_rng_mt19937);
//                 gsl_rng_set(rng[i], 999 + i);  // Seed the RNG
//             }
//         }
        
//         this_NP = P.np/P.cpd; // Assign the same number to every slab
//         int64 extra = P.np%P.cpd;  // Number of leftover particles
//         // There are at most CPD-1 leftover particles
//         // So we can assign one extra particle to as many slabs as necessary
//         if(slab < extra)
//             this_NP++;

//         slab_start = -.5 + ((double)slab)/P.cpd;
//         slab_width = 1./P.cpd;
//     }

//     ~ICFile_Poisson(void) {
//         //gsl_rng_free(rng);  // Can't free, because we need to share the state between objects

        
//     }

//     inline int getparticle(double3 *global_pos, velstruct *vel, auxstruct *aux) {
//         // This routine should return 1 if a particle has been read and
//         // successfully loaded in the given variables.
//         // Return 0 otherwise, which will signal the end of this file.
//         // TODO: need this in the serial version, but count is expensive in the parallel version!
//         //if(count >= this_NP)
//         //    return 0;

//         int t = omp_get_thread_num();
        
//         // Uniform random on [0,1)
//         double rx = gsl_rng_uniform(rng[t]);
//         double ry = gsl_rng_uniform(rng[t]);
//         double rz = gsl_rng_uniform(rng[t]);
        
//         // x is slab-width
//         double x = rx*slab_width + slab_start;
        
//         // y,z are always box-width
//         double y = ry - .5;
//         double z = rz - .5;
        
//         global_pos->x = x;
//         global_pos->y = y;
//         global_pos->z = z;
        
//         // Start the particles from 0 velocity
//         vel->x = 0.;
//         vel->y = 0.;
//         vel->z = 0.;
    
//         // This is too slow!
//         //uint64 pid;
//         //#pragma omp atomic capture
//         //    pid = next_pid++;

//         // Don't increment next_pid as we go; that's a race condition
//         // Instead, use it as the starting point and add an offset of i
//         // We passed i via aux
//         assert(0);  // TOOD: new PID scheme
//         aux->setpid(ICFile::next_pid + aux->aux);
        
//         return 1;
//     }
// };

// class ICFile_Lattice: public ICFile {
// private:
//     int slab;  // slab index
//     //int count;  // number of particles generated in this slab
//     int64_t ppd;  // particles per dimension
//     int64_t ppd2;  // particles per dimension squared

//     int first_plane;  // first plane index in this slab
    
// public:
//     ICFile_Lattice(int _slab) : slab(_slab) {
//         maxthreads = omp_get_max_threads();

//         ppd = WriteState.ppd;
//         ppd2 = ppd*ppd;
//         double ppd_per_slab = (double) ppd / P.cpd;

//         // planes [first,last) are in this slab
//         first_plane = (int) ceil(ppd_per_slab*slab);
//         int last_plane = (int) ceil(ppd_per_slab*(slab+1));

//         this_NP = ppd2*(last_plane - first_plane);
//     }

//     ~ICFile_Lattice(void) {

//     }

//     inline int getparticle(double3 *global_pos, velstruct *vel, auxstruct *aux) {
//         // This routine should return 1 if a particle has been read and
//         // successfully loaded in the given variables.
//         // Return 0 otherwise, which will signal the end of this file.
//         // TODO: need this in the serial version, but count is expensive in the parallel version!
//         //if(count >= this_NP)
//         //    return 0;
        
//         // We passed the index via aux
//         // Map it to a i,j,k tuple and add the plane offset
//         integer3 ijk;
//         ijk.x = aux->aux / ppd2;
//         ijk.y = (aux->aux - ppd2*ijk.x)/ppd;
//         ijk.z = (aux->aux - ppd2*ijk.x) - ppd*ijk.y;
//         ijk.x += first_plane;

//         // then set the PID from that
//         aux->setpid(ijk);
        
//         *global_pos = ZelPos(aux->xyz());
        
//         // Start the particles from 0 velocity
//         vel->x = 0.;
//         vel->y = 0.;
//         vel->z = 0.;
        
//         return 1;
//     }
// };


// gsl_rng **ICFile_Poisson::rng;

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


// An alternative to this macro-heavy approach would be a C++ type registry
#define REGISTER_ICFORMAT(fmt) if(strcasecmp(format, #fmt) == 0){\
    STDLOG(1,"Using format " #fmt "\n");\
    ic = unique_ptr<ICFile_##fmt>(new ICFile_##fmt(_slab));\
} else

// This is a factory function to instantiate ICFile objects of a given format and slab number
// It returns a C++ unique_ptr; no need to call delete on it later.  The object will be deleted
// when it leaves scope.
unique_ptr<ICFile> ICFile::FromFormat(const char *format, int _slab){

    unique_ptr<ICFile> ic;

    // This is actually a big if-else chain; no semicolons!
    /*REGISTER_ICFORMAT(RVdouble)
    REGISTER_ICFORMAT(RVdoubleZel)
    REGISTER_ICFORMAT(RVZel)
    REGISTER_ICFORMAT(RVTag)
    REGISTER_ICFORMAT(RVdoubleTag)
    REGISTER_ICFORMAT(Zeldovich)
    REGISTER_ICFORMAT(Heitmann)
    REGISTER_ICFORMAT(Poisson)
    REGISTER_ICFORMAT(Lattice)*/

    REGISTER_ICFORMAT(RVZel)
    {
        // We weren't given a legal format name.
        QUIT("Unrecognized case: ICFormat = %s\n", format);
    }

    // Nparticles is virtual so that each format can count its particles
    // however it likes.  But virtual functions can't be called from constructors,
    // hence this responsibility is delegated to the factory.
    ic->Npart = ic->Nparticles();

    return ic;
}


uint64 UnpackICtoIL(int slab) {

    // Set up unit conversions
    double convert_pos, convert_vel;
    get_IC_unit_conversions(convert_pos, convert_vel);

    unique_ptr<ICFile> ic = ICFile::FromFormat(P.ICFormat, slab);

    // Unpacks the whole slab directly to the insert list
    uint64 count = ic->unpack_to_IL(convert_pos, convert_vel);

    STDLOG(0,"Read %d particles from IC slab %d\n", count, slab);
    STDLOG(1,"Slab %d has %d subsample A particles, %d subsample B particles.\n", slab, ic->NsubsampleA, ic->NsubsampleB);

    return count; 
}
