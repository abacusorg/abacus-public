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
//    If ICVelocity2Displacement<=-1, then the velocities are to be given in km/s.



// ===================================================================

// To resolve a circular dependency issue, all of the IC classes are defined in this header
#include "IC_classes.h"
#include "particle_subsample.cpp"

uint64 LoadSlab2IL(int slab) {
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
    int ret = snprintf(filename, 1024, "%s/ic_%d",P.InitialConditionsDirectory,slab);
    assert(ret >= 0 && ret < 1024);

    ICfile *ic;    // Abstract class to manipulate the derived class.

    STDLOG(1,"Attempting to read IC file %s\n", filename);
    // TODO: Should make this case-insensitive
    if (strcmp(P.ICFormat,"RVdouble")==0) {
        STDLOG(1,"Using format RVdouble\n");
        ic = new ICfile_RVdouble(filename);
    } else if (strcmp(P.ICFormat,"RVdoubleZel")==0) {
        STDLOG(1,"Using format RVdoubleZel\n");
        ic = new ICfile_RVdoubleZel(filename);
    } else if (strcmp(P.ICFormat,"RVZel")==0) {
        STDLOG(1,"Using format RVZel\n");
        ic = new ICfile_RVZel(filename);
    } else if (strcmp(P.ICFormat,"RVTag")==0) {
        STDLOG(1,"Using format RVTag\n");
        ic = new ICfile_RVTag(filename);
    } else if (strcmp(P.ICFormat,"RVdoubleTag")==0) {
        STDLOG(1,"Using format RVdoubleTag\n");
        ic = new ICfile_RVdoubleTag(filename);
    } else if (strcmp(P.ICFormat,"Zeldovich")==0) {
        STDLOG(1,"Using format Zeldovich\n");
        ic = new ICfile_Zel(filename);
    } else if (strcmp(P.ICFormat,"Heitmann") == 0){
        STDLOG(1,"Using format Heitmann\n");
        ic = new ICfile_Heitmann(filename);
    } else if (strcmp(P.ICFormat, "Poisson") == 0){
        STDLOG(1,"Using format Poisson\n");
        STDLOG(1,"Note: ICFormat \"Poisson\" means that we ignore any IC files and generate the random particles in memory.\n");
        ic = new ICfile_Poisson(slab);
        assert(ic->maxthreads > 1);  // LHG TODO
    } else if (strcmp(P.ICFormat, "Lattice") == 0){
        STDLOG(1,"Using format Lattice\n");
        STDLOG(1,"Note: ICFormat \"Lattice\" means that we ignore any IC files and generate the particles in memory.\n");
        ic = new ICfile_Lattice(slab);
        assert(ic->maxthreads > 1);  // LHG TODO
    }
    else {
        // We weren't given a legal format name.
        QUIT("Unrecognized case: ICFormat = %s\n", P.ICFormat);
    }

    uint64 count = 0;

    STDLOG(3, "IC format permits %d thread(s)\n", ic->maxthreads);
    int sumA = 0; 
    int sumB = 0;

    // TODO: it's hard to OMP a while loop, so for now we have two versions to support multi-threaded in-memory IC generation
    if(ic->maxthreads > 1){

        #pragma omp parallel for schedule(static) num_threads(ic->maxthreads) reduction(+: sumA, sumB)
        for(int64 i = 0; i < ic->this_NP; i++){
            double3 global_pos;
            posstruct pos;
            velstruct vel;
            auxstruct aux;

            // Let's be sneaky and pass in i via aux
            aux.aux = i;

            assertf(ic->getparticle(&global_pos, &vel, &aux) == 1, "Expected count wrong for IC file?\n");

    #ifdef GLOBALPOS
            integer3 newcell = CP->WrapPosition2Cell(&global_pos);
    #else
            integer3 newcell = CP->LocalPosition2Cell(&global_pos);
    #endif
            pos = global_pos;

            vel *= convert_velocity;
            
            int subsample = is_subsample_particle(aux.pid(), P.ParticleSubsampleA, P.ParticleSubsampleB);
            if (subsample == SUBSAMPLE_A) // set all particles 0-ParticleSubsampleA as taggable in set A. 
                { aux.set_taggable_subA(); sumA +=1; }

            else if (subsample == SUBSAMPLE_B) //set all particles A BBig as taggable in set B. 
                { aux.set_taggable_subB(); sumB +=1; }
            
            IL->Push(&pos, &vel, &aux, newcell);
        }

        count = ic->this_NP;

    } else {
        double3 global_pos;
        posstruct pos;
        velstruct vel;
        auxstruct aux;

        // The normal, serial IC code
        while (ic->getparticle(&global_pos, &vel, &aux)) {
            // The getparticle routine is responsible for translation
            // to the unit box and comoving z-space velocities.

            // What cell does this go in?
    #ifdef GLOBALPOS
            // For box-centered positions:
            integer3 newcell = CP->WrapPosition2Cell(&global_pos);
    #else
            // For cell-centered positions:
            integer3 newcell = CP->LocalPosition2Cell(&global_pos);
    #endif
            posstruct pos = global_pos;

            vel *= convert_velocity;
            
            // Set the 'taggable' bit for particles as a function of their PID
            // this bit will be used for group finding and merger trees
            int subsample = is_subsample_particle(aux.pid(), P.ParticleSubsampleA, P.ParticleSubsampleB);
            if (subsample == SUBSAMPLE_A) // set all particles 0-ParticleSubsampleA as taggable in set A. 
                { aux.set_taggable_subA(); sumA +=1; }

            else if (subsample == SUBSAMPLE_B) //set all particles A BBig as taggable in set B. 
                { aux.set_taggable_subB(); sumB +=1; }
            
            IL->Push(&pos, &vel, &aux, newcell);
            count++;
        }
    }

    delete ic;    // Call the destructor.
    STDLOG(0,"Read %d particles from IC file %s\n", count, filename);
    STDLOG(1,"Slab %d has %d subsample A particles, %d subsample B particles.\n", slab, sumA, sumB); 

    return count; 
}
