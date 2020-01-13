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

#ifdef OLD_CODE
void getLightConeFN(int i,int slab, char * fn, char *headerfn){ //ensure the necessary directories exist and return a filename for this light cone.
        char dir1[1080];
        char dir2[1080];
        CheckDirectoryExists(P.LightConeDirectory);
        sprintf(dir1,"%s/LC_raw%.2d",P.LightConeDirectory,i);
        if(!FileExists(dir1)){ //TODO:This assumes the location is not a file. We may want to improve this
            mkdir(dir1,0775);
        }
        sprintf(fn,"%s/LC_raw%.2d/Step%.4d.lc",P.LightConeDirectory,i,WriteState.FullStepNumber);
        sprintf(headerfn,"%s/LC_raw%.2d/header", P.LightConeDirectory,i);
}
#endif

#define c_kms 299792.0
#define etaktoHMpc (c_kms/100.)

#ifdef OLD_CODE
// The positive (rightward) distance between two points (a to b) in a periodic domain
// Essentially implements Python's modulo behavior
inline double positive_dist(double a, double b, double period){
   return fmod(fmod(b - a, period) + period, period);
}

inline int inLightCone(double3 pos, velstruct vel, int lcn, double r1, double r2, double r1tol, double r2tol){ //check if this position is in the current lightcone.
    // Ensure light cone has size
    if (r1 == r2) return false;

    // Get particle position relative to light cone
    double r = (pos-LCOrigin[lcn]).norm();
    double vronc = 0;
    if(r >1e-12) vronc = vel.dot(pos-LCOrigin[lcn])/r * ReadState.VelZSpace_to_kms /c_kms ;
    r = (r-r1) *1/(1+vronc) +r1; //apply velocity correction

    // The lightcone has an absolute output range that is never periodically wrapped
    // We want the lightcone to include particles
    // double drbound = positive_dist(rmin, r, 1.);

    // Check that particle is in light cone
    return (r < r1tol) && (r > r2tol);
}

inline void interpolateParticle(double3 &pos, velstruct &vel, const accstruct acc, const int lcn, double r1, double r2){

    double r = (pos-LCOrigin[lcn]).norm();
    double vronc = 0;
    if(r > 1e-12) vronc = vel.dot(pos-LCOrigin[lcn])/r * ReadState.VelZSpace_to_kms /c_kms;
    double deltar= (r-r1)/(r2-r1) *1/(1+vronc); //apply velocity correction
    if(deltar >=0 or ReadState.DeltaScaleFactor ==0){ //check that we are doing valid interpolation
        double delta_etaK = (WriteState.FirstHalfEtaKick+WriteState.LastHalfEtaKick)*deltar; //how much we need to kick the particle
        double delta_etaD = WriteState.DeltaEtaDrift*deltar; //how much to drift the particle

        pos = pos + delta_etaD * vel;
        vel = vel + delta_etaK * TOFLOAT3(acc);
    }
    else{ //interpolate using the previous timestep information
        double r0 = (cosm->today.etaK -cosm->current.etaK +cosm->KickFactor(ReadState.ScaleFactor -ReadState.DeltaScaleFactor,ReadState.DeltaScaleFactor))*etaktoHMpc/ReadState.BoxSizeMpc;
        double deltar2= (r-r1)/(r0-r1) *1/(1-vronc);

        //assert(deltar2 <=2); //we should not be going back more than one timestep}
        double delta_etaK = (ReadState.FirstHalfEtaKick+ReadState.LastHalfEtaKick)*deltar2; //how much we need to kick the particle
        double delta_etaD = ReadState.DeltaEtaDrift*deltar2; //how much to drift the particle

        pos = pos - delta_etaD * vel;
        vel = vel - delta_etaK * TOFLOAT3(acc);
    }
}
#endif



#include "healpix_shortened.c"

class LightCone {
  private:
    double rmin2;       // Square of rmin
    double rmin_tol2;       // Square of rmin-tol
    double rmax_tol2;       // Square of rmax+tol

  public:
    int lcn;        // Light cone number
    double3 origin;  // The observer location, in unit-cell units
    double rmin;     // The minimum distance to the light cone region (i.e., lower redshift)
    double rmax;     // The maximum distance to the light cone region (i.e., higher redshift)
    double tol;      // The tolerance for our searching
    double DeltaEtaKick;   // The total Delta(eta_Kick) for this step

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
        rmin_tol2 = rmin-tol; rmin_tol2 *= rmin_tol2;
        rmax_tol2 = rmax+tol; rmax_tol2 *= rmax_tol2;
        rmin2 = rmin*rmin;
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
    inline int isParticleInLightCone(double3 &pos, velstruct &vel, const accstruct acc, double3 lineofsight);
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
inline int LightCone::isParticleInLightCone(double3 &pos, velstruct &vel, const accstruct acc, double3 lineofsight) {
    double r = (pos-origin).norm();
    double vr_on_c = vel.dot(lineofsight) * ReadState.VelZSpace_to_kms/c_kms;
    double frac_step = (rmax-r)/(rmax-rmin)*(1.0+vr_on_c);    
        // This is the fraction of the upcoming step when the particle meets the light cone
        // frac_step = 0 means r=rmax, =1 means r-rmin
        // If vr>0, then the particle meets the lightcone earlier

    // Don't let crazy extrapolation happen
    if (frac_step<-0.1) frac_step = -0.1;
    if (frac_step> 1.1) frac_step =  1.1;

    pos += vel*frac_step*WriteState.DeltaEtaDrift;
    // TODO: There should be a way to update r without the full recomputation....
    double r2 = (pos-origin).norm2();
    if (r2>rmax_tol2 || r2<rmin2) return 0;
        // The particle is allowed to be behind the light cone; we want to get it now.
        // But it should be rigorous about the forward edge of the cone.

    // The particle is in the light cone!
    vel += TOFLOAT3(acc)*frac_step*DeltaEtaKick;
    return 1;
}


void makeLightCone(int slab, int lcn){ //lcn = Light Cone Number
    // Use the same format for the lightcones as for the particle subsamples
    if (fabs(cosm->next.etaK-cosm->current.etaK)<1e-12) return;  
            // Nothing to be done, so don't risk divide by zero.

    SlabAccum<RVfloat>   LightConeRV;     ///< The taggable subset in each lightcone.
    SlabAccum<TaggedPID> LightConePIDs;   ///< The PIDS of the taggable subset in each lightcone.
    SlabAccum<unsigned int> LightConeHealPix;   ///< The Healpix of the particles in each lightcone.

    LightConeRV.setup(  CP->cpd, P.np/P.cpd/30);
    LightConePIDs.setup(CP->cpd, P.np/P.cpd/30);
    LightConeHealPix.setup(CP->cpd, P.np/P.cpd);

    // Here are the Slab numbers
    // TODO: These probably should be defined at the bottom, where used.
    SlabType lightcone    = (SlabType)((int)(LightCone0RV + lcn));
    SlabType lightconePID = (SlabType)((int)LightCone0PID + lcn );
    SlabType lightconeHeal = (SlabType)((int)LightCone0Heal + lcn);

    STDLOG(4, "Making light cone %d, slab num %d, w/ pid slab num %d\n", lcn, lightcone, lightconePID);

    LightCone LC(lcn);
    uint64 mask = auxstruct::lightconemask(lcn);
    double vunits = ReadState.VelZSpace_to_kms/ReadState.VelZSpace_to_Canonical;  // Code to km/s

    integer3 ij(slab,0,0);
    uint64_t slabtotal = 0;
    uint64_t slabtotalsub = 0;
    uint64_t slabtotalcell = 0;
    #pragma omp parallel for schedule(dynamic,1) reduction(+:slabtotal) reduction (+:slabtotalsub) reduction(+:slabtotalcell)
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
            // We compute the line of sight direction on a per-cell rather than a per-particle basis.
            // Remember, this is only for a first-order velocity correction to the epoch of the output,
            // and v/c ~ 1e-3.  Plus cells typically subtend 1e-3 radian on the sky, so the radial
            // velocity direction is only changing a tiny bit.
            double3 lineofsight = cc-LCOrigin[lcn];
            lineofsight /= (lineofsight.norm()+1e-15);

            Cell c = CP->GetCell(ijk);
            accstruct *acc = CP->AccCell(ijk);
            #ifdef GLOBALPOS
                cc= 0*cc;
            #endif

            // STDLOG(4, "LC: Particles in current cell: %d\n", c.count());
            for (int p=0;p<c.count();p++) {
                if(!c.aux[p].lightconedone(mask)){   // This particle isn't already in the light cone
                    // Need to unkick by half
                    velstruct vel = c.vel[p] - TOFLOAT3(acc[p])*WriteState.FirstHalfEtaKick;
                    double3 pos = c.pos[p]+cc;  // interpolateParticle takes global positions
                    if (LC.isParticleInLightCone(pos, vel, acc[p], lineofsight)) { 
                        // Yes, it's in the light cone.  pos and vel were updated.

                        pLightConeHealPix->append(LC.healpixel(pos));  // We're outputting all particles for this
                        slabtotal++;

                        if(c.aux[p].is_taggable() or P.OutputFullLightCones){
                            // These output routines take global positions and velocities in km/s
                            pLightConePIDs->append(TaggedPID(c.aux[p]));
                            pLightConeRV->append(RVfloat(pos.x, pos.y, pos.z, vel.x * vunits, vel.y * vunits, vel.z * vunits));
                            slabtotalsub++;

                        }

                        c.aux[p].setlightconedone(mask);
                    }
                }
            }  // Done with this particle

            pLightConePIDs->FinishCell();
            pLightConeRV->FinishCell();
            pLightConeHealPix->FinishCell();
        }   // Done with this cell
        pLightConePIDs->FinishPencil();
        pLightConeRV->FinishPencil();
        pLightConeHealPix->FinishPencil();
    }  // Done with this pencil

    STDLOG(1,"Lightcone %d opened %d cells and found %d particles (%d subsampled) in slab %d\n",
            lcn,slabtotalcell,slabtotal,slabtotalsub,slab);
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

        SB->AllocateSpecificSize(lightconeHeal, slab, LightConeHealPix.get_slab_bytes());
        LightConeHealPix.copy_to_ptr((unsigned int *)SB->GetSlabPtr(lightconeHeal, slab));
        SB->StoreArenaNonBlocking(lightconeHeal, slab);

        // TODO: in a perfect world, we would *sort* the pixel numbers in the healpix slab before outputing
    }

    LightConeRV.destroy();
    LightConePIDs.destroy();
    LightConeHealPix.destroy();
}
