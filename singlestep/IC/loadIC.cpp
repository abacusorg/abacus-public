// This code should load in the particle notionally in one slab
// and add them to the insert list.

// Whether these particles are really part of exactly this slab is 
// not a requirement.  It just has to be close enough.

// Parameter keywords:
// ICFormat: will determine the format of the IC file.
// ICPositionRange: The initial positions will be rescaled by ICPositionRange,
//    e.g., they might be between 0 and ICPositionRange.
//    If ICPositionRange<=0, then it is assumed to be BoxSize.
// ICVelocity2Displacement: The initial velocities will be multiplied by
//    ICVelocity2Displacement to get to redshift space displacements
//    at the initial redshift in the box units of ICPositionRange in size.
//    If ICVelocity2Displacement<=0, then the velocities are to be given in km/s.



// ===================================================================

// A given kind of IC file should be encapsulated in a class as follows:
// Must provide a constructor, destructor, and getparticle(...) method.
// The getparticle routine is required to parse the input file format 
// and translate into the code internal units.

// The constructor should open the file and deal with the header.
// The destructor should close the file.

class ICfile {         // This is an abstract base class.
public:
    virtual inline int getparticle(double3 *global_pos, velstruct *vel, auxstruct *aux)
    =0;
    // Return one particle.
    // Positions must be scaled to the unit box.
    // Velocities must be comoving redshift-space displacements in the unit box.

    // In addition, the derived class must supply a constructor and
    // destructor that open/close the file and deal with the header.
    virtual ~ICfile(void) {};
    // We need a virtual destructor here too.
};

// ===================================================================
// Simple initial conditions: lists of positions & velocities 
// in double precision.  The positions and velocities will be 
// translated by ICPositionRange and ICVelocity2Displacement.

class ICfile_RVdouble: public ICfile {
private:
    double convert_pos, convert_vel;
    FILE *fp;
    class ICparticle {
    public:
        double pos[3];
        double vel[3];
    };
    int read(ICparticle *p) {
        // Must return 0 if read is unsuccessful
        // Return the number of particles read otherwise
        int n = fread(p, sizeof(ICparticle), 1, fp);
        return n;
    }
    int readheader() {
        return 0;  // All's well
    }

public:

    ICfile_RVdouble(char *filename) {
        fp = fopen(filename,"r");
        assertf(fp!=NULL,"Couldn't open IC file %s\n", filename);
        int err = readheader(); 
        assertf(err==0,"Error reading header of IC file %s\n", filename);
        if (P.ICPositionRange>0) convert_pos = 1.0/P.ICPositionRange;
        else convert_pos = 1.0/P.BoxSize;
        if (P.ICVelocity2Displacement>-0.99) 
        convert_vel = P.ICVelocity2Displacement*convert_pos;
        else convert_vel = 1.0/ReadState.VelZSpace_to_kms;
        // This gets to the unit box in redshift-space displacement.
    }

    ~ICfile_RVdouble(void) {
        fclose(fp);
    }

    inline int getparticle(double3 *global_pos, velstruct *vel, auxstruct *aux) {
        // This routine should return 1 if a particle has been read and 
        // successfully loaded in the given variables.
        // Return 0 otherwise, which will signal the end of this file.
        ICparticle p;
        if (read(&p)==0) return 0;
        global_pos->x = p.pos[0];
        global_pos->y = p.pos[1];
        global_pos->z = p.pos[2];
        *global_pos = (*global_pos*convert_pos);  // -double3(0.5,0.5,0.5);
        // This maps to -0.5..+0.5
        if(global_pos->norm() > .5)
        STDLOG(1, "Large IC input\n");
        vel->x = p.vel[0];
        vel->y = p.vel[1];
        vel->z = p.vel[2];
        *vel = *vel * convert_vel;
        aux->clear(); aux->setpid(0); // Set aux too.
        return 1;
    }
};

/*
Same binary format as RVdouble, but positions are interpreted
as Zel-format displacements.  Velocity is only divided by the BoxSize,
just like the displacements.
*/
class ICfile_RVdoubleZel: public ICfile {
public:
    class ICparticle {
    public:
        unsigned short i,j,k;
        double pos[3];
        double vel[3];
    };
private:
    double convert_pos, convert_vel;
    FILE *fp;
    int read(ICparticle *p) {
        // Must return 0 if read is unsuccessful
        // Return the number of particles read otherwise
        int n = fread(p, sizeof(ICparticle), 1, fp);
        return n;
    }
    int readheader() {
        return 0;  // All's well
    }

public:

    ICfile_RVdoubleZel(char *filename) {
        fp = fopen(filename,"r");
        assertf(fp!=NULL,"Couldn't open IC file %s\n", filename);
        int err = readheader();
        assertf(err==0,"Error reading header of IC file %s\n", filename);
        if (P.ICPositionRange>0) convert_pos = 1.0/P.ICPositionRange;
        else convert_pos = 1.0/P.BoxSize;
        if (P.FlipZelDisp)
        convert_pos *= -1;
        if (P.ICVelocity2Displacement>-0.99) // Should always be 1 for IC from zel.cpp
        convert_vel = P.ICVelocity2Displacement*convert_pos;
        else convert_vel = 1.0/ReadState.VelZSpace_to_kms;
        // This gets to the unit box in redshift-space displacement.
    }

    ~ICfile_RVdoubleZel(void) {
        fclose(fp);
    }

    inline int getparticle(double3 *global_pos, velstruct *vel, auxstruct *aux) {
        // This routine should return 1 if a particle has been read and 
        // successfully loaded in the given variables.
        // Return 0 otherwise, which will signal the end of this file.
        ICparticle p;
        if (read(&p)==0) return 0;
        global_pos->x = p.pos[0];
        global_pos->y = p.pos[1];
        global_pos->z = p.pos[2];
        *global_pos = (*global_pos*convert_pos);  // -double3(0.5,0.5,0.5);
        // This maps to -0.5..+0.5
        vel->x = p.vel[0];
        vel->y = p.vel[1];
        vel->z = p.vel[2];
        *vel = *vel * convert_vel;
        
        integer3 ijk;
        ijk.x = p.i;
        ijk.y = p.j;
        ijk.z = p.k;
        *global_pos += ZelPos(ijk);
        
        aux->clear(); aux->setpid(ZelPID(ijk)); // Set aux too.
        return 1;
    }
};

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
// uint64 ZelPID(integer3 ijk): Converts to a PID.
// double3 ZelPos(integer3 ijk): Returns the position in code-units 
//             of the initial grid point.
// double3 ZelPos(uint64 PID): The same, but starting from a PID.


class ICfile_Zel: public ICfile {
private:
    double convert_pos, convert_vel;
    FILE *fp;
    class ICparticle {
    public:
        unsigned short i,j,k;
        double displ[3];
    };
    int read(ICparticle *p) {
        // Must return 0 if read is unsuccessful
        // Return the number of particles read otherwise
        int n = fread(p, sizeof(ICparticle), 1, fp);
        //        printf("(%d,%d,%d) : (%f, %f, %f)\t \t sizeof(ICparticle) = %d\n",p->i,p->j,p->k, p->displ[0],p->displ[1],p->displ[2],sizeof(ICparticle));
        return n;
    }
    int readheader() {
        char buf[2];
        int len = 1;
        fread(&(buf[1]),1,1,fp);
        do {
            buf[0] = buf[1];
            if(fread(&(buf[1]), 1, 1, fp)<=0) QUIT("Error reading header for zeldovich input file. Exiting\n");
            len++;
        } while(!(buf[0]==0x2 && buf[1]==0x2));
        return 0;  // All's well
    }

public:

    ICfile_Zel(char *filename) {
        fp = fopen(filename,"r");

        assertf(fp!=NULL,"Couldn't open IC file %s\n", filename);
        //int err = readheader();
        //assertf(err==0,"Error reading header of IC file %s\n", filename);
        if (P.ICPositionRange>0) convert_pos = 1.0/P.ICPositionRange;
        else convert_pos = 1.0/P.BoxSize;
        // This will convert the displacement to the code units
        convert_vel = WriteState.f_growth;
        // TODO: Check that this is the correct physics.
        // Zel'dovich velocities in redshift-space displacements units
        // are just f*position_dispacement.
    }

    ~ICfile_Zel(void) {
        fclose(fp);
    }

    inline int getparticle(double3 *global_pos, velstruct *vel, auxstruct *aux) {
        // This routine should return 1 if a particle has been read and 
        // successfully loaded in the given variables.
        // Return 0 otherwise, which will signal the end of this file.
        ICparticle p;
        if (read(&p)==0) return 0;
        global_pos->x = p.displ[0];
        global_pos->y = p.displ[1];
        global_pos->z = p.displ[2];
        if (P.FlipZelDisp)
        *global_pos *= -1;
        *global_pos = *global_pos * convert_pos;
        // This is the displacement in code units (unit box)
        *vel = *global_pos * convert_vel;
        // The velocity is related to the displacement.
        // Now add on the initial Lagrangian position
        integer3 ijk;
        ijk.x = p.i;
        ijk.y = p.j;
        ijk.z = p.k;
        *global_pos += ZelPos(ijk);

        aux->clear(); aux->setpid(ZelPID(ijk)); // Set aux too.
        return 1;
    }
};


class ICfile_Heitmann: public ICfile {
private:
    double convert_pos, convert_vel;
    FILE *fp;
    class ICparticle {
    public:
        float xv[6];
        float mass;
        unsigned int tag;
    };
    int read(ICparticle *p) {
        // Must return 0 if read is unsuccessful
        // Return the number of particles read otherwise
        int n = fread(p, sizeof(ICparticle), 1, fp);
        return n;
    }
    int readheader() {
        return 0;  // All's well
    }

public:

    ICfile_Heitmann(char *filename) {
        fp = fopen(filename,"r");
        assertf(fp!=NULL,"Couldn't open IC file %s\n", filename);
        int err = readheader();
        assertf(err==0,"Error reading header of IC file %s\n", filename);
        if (P.ICPositionRange>0) convert_pos = 1.0/P.ICPositionRange;
        else convert_pos = 1.0/P.BoxSize *(P.H0/100.0);
        if (P.ICVelocity2Displacement>-0.99) 
        convert_vel = P.ICVelocity2Displacement*convert_pos;
        else convert_vel = 1.0/ReadState.VelZSpace_to_kms * (1/(1+ P.InitialRedshift));
        // This gets to the unit box in redshift-space displacement.
    }

    ~ICfile_Heitmann(void) {
        fclose(fp);
    }

    inline int getparticle(double3 *global_pos, velstruct *vel, auxstruct *aux) {
        // This routine should return 1 if a particle has been read and
        // successfully loaded in the given variables.
        // Return 0 otherwise, which will signal the end of this file.
        ICparticle p;
        if (read(&p)==0) return 0;
        global_pos->x = p.xv[0];
        global_pos->y = p.xv[2];
        global_pos->z = p.xv[4];
        *global_pos = (*global_pos*convert_pos);  // -double3(0.5,0.5,0.5);
        // This maps to -0.5..+0.5
        vel->x = p.xv[1];
        vel->y = p.xv[3];
        vel->z = p.xv[5];
        *vel = *vel * convert_vel;
        aux->clear(); aux->setpid(p.tag); // Set aux too.
        return 1;
    }
};

void LoadSlab2IL(int slab) {
    double3 global_pos;
    posstruct pos;
    velstruct vel;
    auxstruct aux;
    
    double convert_velocity = 1.0;
    // The getparticle() routine should return velocities in 
    // redshift-space comoving displacements in length units where
    // the box is unit length.  
    // We need to convert to canonical velocities.
    // Code time unit is 1/H_0, but H(z) gets output with H0=1.
    // The conversion is v_code = v_zspace * a^2 H(z)/H_0
    // TODO: Have to check that the units of H(z) and H_0 are the same!
    convert_velocity = WriteState.VelZSpace_to_Canonical;
    // WriteState.ScaleFactor*WriteState.ScaleFactor*WriteState.HubbleNow;

    char filename[1024];
    sprintf(filename,"%s/ic_%d",P.InitialConditionsDirectory,slab);

    ICfile *ic;    // Abstract class to manipulate the derived class.

    STDLOG(1,"Attempting to read IC file %s\n", filename);
    // TODO: Should make this case-insensitive
    if (strcmp(P.ICFormat,"RVdouble")==0) {
        STDLOG(1,"Using format RVdouble\n");
        ic = new ICfile_RVdouble(filename);
    } else if (strcmp(P.ICFormat,"RVdoubleZel")==0) {
        STDLOG(1,"Using format RVdoubleZel\n");
        ic = new ICfile_RVdoubleZel(filename);
    } else if (strcmp(P.ICFormat,"Zeldovich")==0) {
        STDLOG(1,"Using format Zeldovich\n");
        ic = new ICfile_Zel(filename);
    } else if (strcmp(P.ICFormat,"Heitmann") == 0){
        STDLOG(1,"Using format Heitmann\n");
        ic = new ICfile_Heitmann(filename);
    }
    else {
        // We weren't given a legal format name.
        QUIT("Unrecognized case: ICFormat = %s\n", P.ICFormat);
    }

    uint64 count = 0;
    while (ic->getparticle(&global_pos, &vel, &aux)) {
        // The getparticle routine is responsible for translation
        // to the unit box and comoving z-space velocities.

        // What cell does this go in?
#ifdef GLOBALPOS
        // For box-centered positions:
        integer3 newcell = PP->WrapPosition2Cell(&global_pos);
#else
        // For cell-centered positions:
        integer3 newcell = PP->LocalPosition2Cell(&global_pos);
#endif
        pos = global_pos;    

        vel *= convert_velocity;
        IL->Push(&pos, &vel, &aux, newcell);
        count++;
    }
    delete ic;    // Call the destructor.
    STDLOG(0,"Read %d particles from IC file %s\n", count, filename);
    return;
}
