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

void getLightConeFN(int i,int slab, char * fn){ //ensure the necessary directories exist and return a filename for this light cone.
        char dir1[1080];
        char dir2[1080];
        CheckDirectoryExists(P.LightConeDirectory);
        sprintf(dir1,"%s/LC_raw%.2d",P.LightConeDirectory,i);
        if(!FileExists(dir1)){ //TODO:This assumes the location is not a file. We may want to improve this
            mkdir(dir1,0775);
        }
        sprintf(dir2,"%s/LC_raw%.2d/Step%.4d",P.LightConeDirectory,i,WriteState.FullStepNumber);
        if (!FileExists(dir2)){
            mkdir(dir2,0775);
        }
        CheckDirectoryExists(dir2);
        sprintf(fn,"%s/LC_raw%.2d/Step%.4d/slab%.4d.lc",P.LightConeDirectory,i,WriteState.FullStepNumber,slab);

}

// The positive (rightward) distance between two points (a to b) in a periodic domain
// Essentially implements Python's modulo behavior
inline double positive_dist(double a, double b, double period){
   return fmod(fmod(b - a, period) + period, period);
}

#define c_kms 299792.0
#define etaktoMpc (c_kms/100.)
inline int inLightCone(double3 pos, velstruct vel, int lcn, double tol1, double tol2){ //check if this position is in the current lightcone.
    double r1 = (cosm->today.etaK - cosm->current.etaK)*etaktoMpc/ReadState.BoxSizeHMpc;  // lightcone start
    double r2 = (cosm->today.etaK - cosm->next.etaK)*etaktoMpc/ReadState.BoxSizeHMpc;  // lightcone end
    if (r1 ==r2) return false;
    double r = (pos-LCOrigin[lcn]).norm();
    double vronc = 0;
    if(r >1e-12) vronc = vel.dot(pos-LCOrigin[lcn])/r * ReadState.VelZSpace_to_kms /c_kms ;
    r = (r-r1) *1/(1+vronc) +r1; //apply velocity correction
    
    // Expand the valid range by the upper and lower tolerances
    double rmax = r1 + tol1;
    double rmin = r2 - tol2;
    double lcwidth = rmax - rmin;
    // The lightcone has an absolute output range that is never periodically wrapped
    // We want the lightcone to include particles 
    double drbound = positive_dist(rmin, r, 1.);
    
    int result = (r < rmax) && (r > rmin);
    return result;
}

inline void interpolateParticle(double3 &pos, velstruct &vel, const accstruct acc, const int lcn){
    double r1 = (cosm->today.etaK - cosm->current.etaK)*etaktoMpc/ReadState.BoxSizeMpc;
    double r2 = (cosm->today.etaK - cosm->next.etaK)*etaktoMpc/ReadState.BoxSizeMpc;

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
        double r0 = (cosm->today.etaK -cosm->current.etaK +cosm->KickFactor(ReadState.ScaleFactor -ReadState.DeltaScaleFactor,ReadState.DeltaScaleFactor))*etaktoMpc/ReadState.BoxSizeMpc;
        double deltar2= (r-r1)/(r0-r1) *1/(1-vronc);

        //assert(deltar2 <=2); //we should not be going back more than one timestep}
        double delta_etaK = (ReadState.FirstHalfEtaKick+ReadState.LastHalfEtaKick)*deltar2; //how much we need to kick the particle
        double delta_etaD = ReadState.DeltaEtaDrift*deltar2; //how much to drift the particle

        pos = pos - delta_etaD * vel;
        vel = vel - delta_etaK * TOFLOAT3(acc);
    }
}

void makeLightCone(int slab,int lcn){ //lcn = Light Cone Number
    //double r1 = cosm->today.H*(cosm->today.etaK - cosm->current.etaK)*4000/P.BoxSize;
    //double r2 = cosm->today.H*(cosm->today.etaK - cosm->next.etaK)*4000/P.BoxSize;
    //printf("r1: %f r2: %f \n",r1,r2);
    // Use the same format for the lightcones as for the timeslices
    AppendArena *AA = get_AA_by_format(P.OutputFormat);

    SlabType lightcone = (SlabType)((int)(LightCone0 + lcn));
    int headersize = 1024*1024;
    uint64_t slabtotal = 0;
    uint64_t strayslabtotal = 0;

    SB->AllocateSpecificSize(lightcone, slab,
            // allocate enough space for the case where every particle is "stray" and gets its own cell
            SS->size(slab)*(AA->sizeof_particle() + AA->sizeof_cell()) + headersize);
    AA->initialize(lightcone, slab, CP->cpd, ReadState.VelZSpace_to_Canonical);

    if (!P.OmitOutputHeader) {
            AA->addheader((const char *) P.header());
            AA->addheader((const char *) ReadState.header());
            AA->finalize_header();
    }
    uint64 mask = auxstruct::lightconemask(lcn);

    double kickfactor = WriteState.FirstHalfEtaKick;
    velstruct vel;
    double3 pos;  // some positions will be global; ensure double
    integer3 ijk(slab,0,0);
    // if a particle interpolates outside its cell, store its index here
    vector<int> stray_particles;
    stray_particles.reserve(P.np * CP->invcpd3);
    
    for (ijk.y=0; ijk.y<CP->cpd; ijk.y++){
        for (ijk.z=0;ijk.z<CP->cpd;ijk.z++) {
            // Check if the cell center is in the lightcone, with some wiggle room
            if(!inLightCone(CP->CellCenter(ijk),0*pos,lcn,3.0/P.cpd, 3.0/P.cpd))
                continue;
            Cell c = CP->GetCell(ijk);
            // vscale should be the max possible velocity
            // which could be greater than the cell stats would suggest, because we interpolate
            // we could compute the maximum possible interpolated velocity, or just fudge a bit
            float vscale = c.ci->max_component_velocity/ReadState.VelZSpace_to_Canonical;
            accstruct *acc = CP->AccCell(ijk);
            uint64_t cell_np = 0;

            double3 cc = CP->CellCenter(ijk);
#ifdef GLOBALPOS
            cc= 0*cc;
#endif
            for (int p=0;p<c.count();p++) {
                if(!c.aux[p].lightconedone(mask)){
                    vel = (c.vel[p] - TOFLOAT3(acc[p])*kickfactor);
                    if(inLightCone((c.pos[p]+cc),vel,lcn,1.0,0)){
                        pos = c.pos[p]+cc;  // interpolateParticle takes global positions
                        interpolateParticle(pos,vel,acc[p],lcn);
                        // Check which cell the interpolated position falls within
                        // and make the position local to that cell
                        // If the interpolated position has left the cell, deal with it after we've closed the current cell
                        // Same with exceeding the max velocity
                        if(CP->LocalPosition2Cell(&pos) != ijk || vel.maxabscomponent() > c.ci->max_component_velocity){
                            stray_particles.push_back(p);
                            continue;
                        }
                        // Only open the cell here, where we know for sure a particle is going to be written
                        // might have already opened it, though
                        if(!AA->current_cell.islegal())
                            AA->addcell(ijk,vscale);
                        // AA takes cell-centered positions, which LocalPosition2Cell gives
                        AA->addparticle(pos, vel, c.aux[p]);
                        cell_np++;
                        c.aux[p].setlightconedone(mask);
                    }
                }
            }
            AA->endcell();  // can always end cell, even if we never opened one
            slabtotal += cell_np;

            cell_np = 0;
            // Now make a singlet cell for each stray particle
            while(!stray_particles.empty()){
                int p = stray_particles.back();
                stray_particles.pop_back();
                pos = c.pos[p] + cc;  // use the current cell center to make global
                vel = c.vel[p] - TOFLOAT3(acc[p])*kickfactor;
                interpolateParticle(pos,vel,acc[p],lcn);
                
                // find the new cell
                integer3 stray_ijk = CP->LocalPosition2Cell(&pos);  // pos now local
                int slab_distance = abs(stray_ijk.x - slab);
                slab_distance -= P.cpd*round(slab_distance/P.cpd);
                // one of our guarantees should be that particles in a slab file are never actually more than 1 slab away
                // but how can we ensure that without a hard assert?
                //assertf(slab_distance <= 1, "Lightcone particle in slab %d interpolated too many slabs away (to slab %d).\n", slab, stray_ijk.x);
                Cell stray_cell = CP->GetCell(stray_ijk);

                // add one particle and close
                // for a single particle, the vscale can be its own velocity
                float stray_vscale = vel.maxabscomponent()/ReadState.VelZSpace_to_Canonical;
                AA->addcell(stray_ijk, stray_vscale);
                AA->addparticle(pos, vel, c.aux[p]);
                c.aux[p].setlightconedone(mask);
                cell_np++;
                AA->endcell();
            }
            strayslabtotal += cell_np;
            slabtotal += cell_np;
        }
    }
    STDLOG(1,"Slab %d had %d particles (%d stray) in lightcone %d\n",slab,slabtotal,strayslabtotal,lcn);
    if(slabtotal) {
        // Write out this filename
        char filename[1024];
        getLightConeFN(lcn,slab,filename);
        SB->ResizeSlab(lightcone, slab, AA->bytes_written());
        // TODO: make sure this interacts as we want with the ramdisk
        SB->WriteArena(lightcone, slab, IO_DELETE, IO_NONBLOCKING, filename);
    } else {
        // No particles in this lc; don't write anything
        SB->DeAllocate(lightcone, slab);
    }
    delete AA;
}

