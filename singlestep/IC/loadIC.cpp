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
//    If ICVelocity2Displacement==-1, then the velocities are to be given in km/s.



// ===================================================================

// To resolve a circular dependency issue, all of the IC classes are defined in this header
#include "IC_classes.h"
#include "particle_subsample.cpp"




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


uint64 UnpackICtoIL(int slab) {

    // Set up unit conversions
    double convert_pos, convert_vel;
    get_IC_unit_conversions(convert_pos, convert_vel);

    ICFile *ic = ICFile::FromFormat(P.ICFormat, slab, convert_pos, convert_vel);

    // Unpacks the whole slab directly to the insert list
    uint64 count = ic->unpack_to_IL();

    STDLOG(0,"Read %d particles from IC slab %d\n", count, slab);
    STDLOG(1,"Slab %d has %d subsample A particles, %d subsample B particles.\n", slab, ic->NsubsampleA, ic->NsubsampleB);

    delete ic;    // Call the destructor.

    return count; 
}
