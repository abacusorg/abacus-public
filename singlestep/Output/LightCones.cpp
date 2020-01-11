/*
 * LightCones.cpp
 *
 *  Created on: Oct 1, 2012
 *      Author: dferrer
 *
 *
 *      Identify particles in the lightcone for a given set of origins and output them.
 *
 *      The particles are currently output in the following scheme:
 *      filename --"LCN/StepI/slabJ.lc"
 *      a
 *      da
 *      NCells
 *      cellheader (center, NP)
 *      *particles (pos, vel, z)*
 *      cellheader (center, NP)
 *      *particles (pos, vel, z)*
 *
 */

#define LCTOLERANCE 0.0

double3 *LCOrigin;

#define c_kms 299792.0
#define etaktoHMpc (c_kms/100.)

#include "healpix_shortened.c"

class LightCone {
  private:
    double rmin_tol2;       // Square of rmin-tol
    double rmax_tol2;       // Square of rmax+tol

  public:
    int lcn;        // Light cone number
    double3 origin;  // The observer location, in unit-cell units
    double rmin;     // The minimum distance to the light cone region (i.e., lower redshift)
    double rmax;     // The maximum distance to the light cone region (i.e., higher redshift)
    double tol;      // The tolerance for our searching
    FLOAT DeltaEtaKick;   // The total Delta(eta_Kick) for this step
    FLOAT driftfactor;

    LightCone(int _lcn) {
        lcn = _lcn;
        origin = LCOrigin[lcn];

        // Bounds for light cone, in unit-box units
        rmax = (cosm->today.etaK - cosm->current.etaK)*etaktoHMpc/ReadState.BoxSizeHMpc; // Light cone start
        rmin = (cosm->today.etaK - cosm->next.etaK)*etaktoHMpc/ReadState.BoxSizeHMpc; // Light cone end
    
        // And we need to set some tolerances.
        // A particle have been in front of the light cone in the previous step, but then
        // move so that it is behind the light cone in this step.  Particles can move up to FINISH_RADIUS.
        // And we want to catch this particle even if it is in the corner of the cell, 
        // which means that we have to catch its cell center too.
        // So this is FINISH_WAIT_RADIUS + sqrt(3)/2 cells
        tol = (FINISH_WAIT_RADIUS+sqrt(3.0)/2.0)/P.cpd;   
        // TODO: Could consider stricter bounds, e.g., using MaxVelocity from the previous state.
        rmin_tol2 = rmin-tol; rmin_tol2 *= rmin_tol2;
        rmax_tol2 = rmax+tol; rmax_tol2 *= rmax_tol2;
        driftfactor = WriteState.DeltaEtaDrift;
        DeltaEtaKick = WriteState.FirstHalfEtaKick+WriteState.LastHalfEtaKick;
    }
    ~LightCone() { }

    // This is a simple wrapper for our healpix calls
    inline unsigned int healpixel(double3 pos) {
        int64_t pixel;
        pos -= origin;
        vec2pix_nest64((int64_t)16384, (double *)&pos, &pixel);
        return (unsigned int) pixel;
    }

    inline int isCellInLightCone(double3 pos);
    inline int isParticleInLightCone(double3 cellcenter, posstruct &pos, velstruct &vel, const accstruct acc);
};


// Return whether a CellCenter is in the light cone, including some tolerance
inline int LightCone::isCellInLightCone(double3 pos) {
    double r2 = (pos-origin).norm2();
    return (r2<rmax_tol2) && (r2>rmin_tol2);
}

// pos must be a global position in unit-box coords; vel is in code units
// rmax is the maximum radius, which is the earlier time
// rmin is the minimum radius, which is the later time
// lineofsight should be a unit vector from the LCOrigin to the cell center
//
// Returns 0 if not in LightCone, 1 if it is.
// Further, the position and velocity inputs will be adjusted.
inline int LightCone::isParticleInLightCone(double3 cellcenter, posstruct &pos, velstruct &vel, const accstruct acc) {
    double r0 = (cellcenter-origin+pos).norm();
    posstruct pos1 = pos+vel*driftfactor;   // Take care to match the precision of Drift()
    // Now rebin pos1, matching the precision of Insert()
    double3 cc1 = cellcenter;
    if (pos1.x>CP->halfinvcpd) { pos1.x-=CP->halfinvcpd; cc1.x+=CP->halfinvcpd; }
    if (pos1.y>CP->halfinvcpd) { pos1.y-=CP->halfinvcpd; cc1.y+=CP->halfinvcpd; }
    if (pos1.z>CP->halfinvcpd) { pos1.z-=CP->halfinvcpd; cc1.z+=CP->halfinvcpd; }
    if (pos1.x<-CP->halfinvcpd) { pos1.x+=CP->halfinvcpd; cc1.x-=CP->halfinvcpd; }
    if (pos1.y<-CP->halfinvcpd) { pos1.y+=CP->halfinvcpd; cc1.y-=CP->halfinvcpd; }
    if (pos1.z<-CP->halfinvcpd) { pos1.z+=CP->halfinvcpd; cc1.z-=CP->halfinvcpd; }
    double r1 = (cc1-origin+pos1).norm();

    double frac_step = (rmax-r0)/(rmax-rmin-r0+r1);
        // This is the fraction of the upcoming step when the particle meets the light cone
        // frac_step = 0 means r=rmax, =1 means r-rmin

    if (frac_step<-1.0e-13||frac_step>=1) return 0;
        // We accept the particle into the lightcone only if the two lines cross in
        // the domain of the step.
        // We are accepting a tiny fraction of cases outside the cone, 
        // just in case of floating point math errors

    // The particle is in the light cone!
    // Update the pos and vel to the fractional step (for output).
    pos += vel*driftfactor*frac_step;
    vel += TOFLOAT3(acc)*DeltaEtaKick*frac_step;
    return 1;
}


void makeLightCone(int slab, int lcn){ //lcn = Light Cone Number
    // Use the same format for the lightcones as for the particle subsamples
    SlabAccum<RVfloat>   LightConeRV;     ///< The taggable subset in each lightcone.
    SlabAccum<TaggedPID> LightConePIDs;   ///< The PIDS of the taggable subset in each lightcone.
    SlabAccum<unsigned int> LightConeHealPix;   ///< The Healpix of the particles in each lightcone.

    LightConeRV.setup(  CP->cpd, P.np/P.cpd/30);
    LightConePIDs.setup(CP->cpd, P.np/P.cpd/30);
    LightConeHealPix.setup(CP->cpd, P.np/P.cpd);

    // Here are the Slab numbers
    SlabType lightcone    = (SlabType)((int)(LightCone0 + lcn));
    SlabType lightconePID = (SlabType)((int)LightCone0 + lcn + NUMLIGHTCONES);
    SlabType lightconeHealPix = (SlabType)((int)LightCone0 + lcn + NUMLIGHTCONES*2);

    STDLOG(4, "Making light cone %d, slab num %d, w/ pid slab num %d\n", lcn, lightcone, lightconePID);

    LightCone LC(lcn);
    uint64 mask = auxstruct::lightconemask(lcn);
    double vunits = ReadState.VelZSpace_to_kms/ReadState.VelZSpace_to_Canonical;  // Code to km/s

    integer3 ij(slab,0,0);
    uint64_t slabtotal = 0;
    uint64_t slabtotalsub = 0;
    uint64_t slabtotalcell = 0;
    uint64_t doubletagged = 0;
    #pragma omp parallel for schedule(dynamic,1) reduction(+:slabtotal) reduction (+:slabtotalsub) reduction(+:slabtotalcell) reduction(+:doubletagged)
    for (int y = 0; y < CP->cpd; y ++) {
        integer3 ijk = ij; ijk.y = y;

        PencilAccum<RVfloat>   *pLightConeRV   =   LightConeRV.StartPencil(ijk.y);
        PencilAccum<TaggedPID> *pLightConePIDs = LightConePIDs.StartPencil(ijk.y);
        PencilAccum<unsigned int> *pLightConeHealPix = LightConeHealPix.StartPencil(ijk.y);

        for (ijk.z=0;ijk.z<CP->cpd;ijk.z++) {
            // Check if the cell center is in the lightcone, with some wiggle room
            double3 cc = CP->CellCenter(ijk);
            if(!LC.isCellInLightCone(cc)) continue;  // Skip the rest if too far from the region

            // So you say there's a chance?
            slabtotalcell++;

            Cell c = CP->GetCell(ijk);
            accstruct *acc = CP->AccCell(ijk);
            #ifdef GLOBALPOS
                cc= 0*cc;
            #endif

            // STDLOG(4, "LC: Particles in current cell: %d\n", c.count());
            for (int p=0;p<c.count();p++) {
                // if(!c.aux[p].lightconedone(mask)){   // This particle isn't already in the light cone
                    // Need to unkick by half
                    velstruct vel = c.vel[p] - TOFLOAT3(acc[p])*WriteState.FirstHalfEtaKick;
                    posstruct poscopy = c.pos[p];  // Need a copy, since it will be changed
                    if (LC.isParticleInLightCone(cc, poscopy, vel, acc[p])) { 
                        // Yes, it's in the light cone.  pos and vel were updated.
                        double3 pos = cc+poscopy;    // This is now a global position in double precision

                        pLightConeHealPix->append(LC.healpixel(pos));  // We're outputting all particles for this
                        slabtotal++;

                        if(c.aux[p].is_taggable() or P.OutputFullLightCones){
                            // These output routines take global positions and velocities in km/s
                            pLightConePIDs->append(TaggedPID(c.aux[p]));
                            pLightConeRV->append(RVfloat(pos.x, pos.y, pos.z, vel.x * vunits, vel.y * vunits, vel.z * vunits));
                            slabtotalsub++;

                        }

                        // TODO: For now, we're going to look for particles that get tagged twice.
                        // But maybe we'll find that they are very few, in which case we might stop tagging.
                        if (c.aux[p].lightconedone(mask)) doubletagged++;
                        c.aux[p].setlightconedone(mask);
                    }
                // }
            }  // Done with this particle

            pLightConePIDs->FinishCell();
            pLightConeRV->FinishCell();
            pLightConeHealPix->FinishCell();
        }   // Done with this cell
        pLightConePIDs->FinishPencil();
        pLightConeRV->FinishPencil();
        pLightConeHealPix->FinishPencil();
    }  // Done with this pencil

    STDLOG(1,"Lightcone %d opened %d cells and found %d particles (%d subsampled) in slab %d.  %d double tagged\n",
            lcn,slabtotalcell,slabtotal,slabtotalsub,slab, doubletagged);
    if(slabtotal) {
        #ifdef OLDCODE
        // Find filename for consistency, but writing to pointer anyway
        char filename[1024];
        char headername[1024];
        getLightConeFN(lcn,slab,filename, headername);
        //WriteLightConeHeaderFile(headername.c_str());
        #endif

        // TODO: Someone should write a header for the light cone.

        SB->AllocateSpecificSize(lightcone, slab, LightConeRV.get_slab_bytes());
        LightConeRV.copy_to_ptr((RVfloat *)SB->GetSlabPtr(lightcone, slab));
        SB->StoreArenaNonBlocking(lightcone, slab);

        SB->AllocateSpecificSize(lightconePID, slab, LightConePIDs.get_slab_bytes());
        LightConePIDs.copy_to_ptr((TaggedPID *)SB->GetSlabPtr(lightconePID, slab));
        SB->StoreArenaNonBlocking(lightconePID, slab);

        // TODO: Need to add the HealPix files
        // And in a perfect world, we would *sort* the pixel numbers in the healpix slab before outputing
    }

    LightConeRV.destroy();
    LightConePIDs.destroy();
    LightConeHealPix.destroy();
}
