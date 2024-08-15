/*
 * LightCones.cpp
 *      Identify particles in the lightcone for a given set of origins and output them.
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

    static void WriteHeaderFile(const char* fn){
        std::ofstream headerfile;
        headerfile.open(fn);
        headerfile << P.header();
        headerfile << ReadState.header();
        headerfile << "\nOutputType = \"LightCone\"\n";
        headerfile.close();
    }
};


// Return whether a CellCenter is in the light cone, including some tolerance
inline int LightCone::isCellInLightCone(double3 pos) {
    double r2 = (pos-origin).norm2();
    return (r2<=rmax_tol2) && (r2>=rmin_tol2);
}

// pos is in cell-centered coords; cellcenter is the position in the box; vel is in code units
// rmax is the maximum radius, which is the earlier time
// rmin is the minimum radius, which is the later time
//
// Returns 0 if not in LightCone, 1 if it is.
// Further, the position and velocity inputs will be adjusted to the epoch where the light cone is crossed.
inline int LightCone::isParticleInLightCone(double3 cellcenter, posstruct &pos, velstruct &vel, const accstruct acc) {
    double r0 = (cellcenter-origin+pos).norm();
    posstruct pos1 = pos+vel*driftfactor;   // Take care to match the precision of Drift()
    /*
    // Now rebin pos1, matching the precision of Insert()
    double3 cc1 = cellcenter;
    if (pos1.x>CP->halfinvcpd) { pos1.x-=CP->halfinvcpd; cc1.x+=CP->halfinvcpd; }
    if (pos1.y>CP->halfinvcpd) { pos1.y-=CP->halfinvcpd; cc1.y+=CP->halfinvcpd; }
    if (pos1.z>CP->halfinvcpd) { pos1.z-=CP->halfinvcpd; cc1.z+=CP->halfinvcpd; }
    if (pos1.x<-CP->halfinvcpd) { pos1.x+=CP->halfinvcpd; cc1.x-=CP->halfinvcpd; }
    if (pos1.y<-CP->halfinvcpd) { pos1.y+=CP->halfinvcpd; cc1.y-=CP->halfinvcpd; }
    if (pos1.z<-CP->halfinvcpd) { pos1.z+=CP->halfinvcpd; cc1.z-=CP->halfinvcpd; }
    // This attempt to improve the precision handling of particles changing cells 
    // didn't actually result in perfectly reproducible answers, so we will rely on the LC bits.
    */
    double r1 = (cellcenter-origin+pos1).norm();

    double frac_step = (rmax-r0)/(rmax-rmin-r0+r1);
        // This is the fraction of the upcoming step when the particle meets the light cone
        // frac_step = 0 means r=rmax, =1 means r-rmin

#ifdef USE_LC_AUX_BITS
    if (frac_step<-1.0e-6||frac_step>=1) return 0;

        // We accept the particle into the lightcone only if the two lines cross in
        // the domain of the step.
        // We are accepting a tiny fraction of cases outside the cone, 
        // just in case of floating point math errors. The aux bits ensure
        // we won't double count.
#else
    if (frac_step<0||frac_step>=1) return 0;

    // But if we're not using aux bits, try to be more precise.
#endif

    // The particle is in the light cone!
    // Update the pos and vel to the fractional step (for output).
    pos += vel*driftfactor*frac_step;
    vel += TOFLOAT3(acc)*DeltaEtaKick*frac_step;
    return 1;
}


/* Geometry recommendation:
If the light cone sphere becomes tangent to the y-z slab, then one will get
more particles selected.  That is probably disfavorable to the pipeline workflow.
This effect is only severe if the observer is far away, so that the radius
of curvature of the sphere is large.  But this disfavors have abutting light
cones where the observers are periodically extended in the x direction.

We are parallelizing over y pencils, so it might be better if more y pencils
have some work.  That would favor observers that are extended in the z direction,
so that the light cone sphere approaches x-y planes.  Such planes will have 
more equal work in x slabs and y pencils, but only for a narrow range of z.

*/


size_t makeLightCone(int slab, int lcn){ //lcn = Light Cone Number
    // Use the same format for the lightcones as for the particle subsamples
    if (fabs(cosm->next.etaK-cosm->current.etaK)<1e-12) return 0;
          // Nothing to be done, so don't risk divide by zero.
    
    OutputLightConeSetup.Start();
    STDLOG(4, "Making light cone %d for slab %d\n", lcn, slab);

    LightCone LC(lcn);

    // Before we start allocating memory, do a quick pass to see if any cells are in the LC
    uint64_t slabtotalcell = 0;

    #pragma omp parallel for schedule(static) reduction(+:slabtotalcell)
    for (int y = 0; y < CP->cpd; y ++) {
        for (int z = node_z_start; z < node_z_start + node_z_size; z++) {
            // Check if the cell center is in the lightcone, with some wiggle room
            double3 cc = CP->CellCenter(slab, y, z);
            if(LC.isCellInLightCone(cc)) slabtotalcell++;
        }
    }

    if(slabtotalcell == 0){
        OutputLightConeSetup.Stop();
        return 0;
    }

    SlabAccum<RVfloat>   LightConeRV;     ///< The taggable subset in each lightcone.
    SlabAccum<TaggedPID> LightConePIDs;   ///< The PIDS of the taggable subset in each lightcone.
    SlabAccum<unsigned int> LightConeHealPix;   ///< The Healpix of the particles in each lightcone.

    // The size estimates are objects per slab.
    // The light cones output only taggables particles for RV and PID, 
    // but all particles for HealPix.
    // What fraction of a slab is the light cone?  It depends on geometry and the time step.
    // But for 1-2 Gpc boxes, it takes a few hundred time steps to cross the box.  We'll guess 1%.
    // That said, if the slab is tangent to the annulus, one can include >1% of the slab.
    LightConeRV.setup(  CP->cpd, node_z_size, P.np/P.cpd*(P.ParticleSubsampleA+P.ParticleSubsampleB)/100);
    LightConePIDs.setup(CP->cpd, node_z_size, P.np/P.cpd*(P.ParticleSubsampleA+P.ParticleSubsampleB)/100);
    LightConeHealPix.setup(CP->cpd, node_z_size, P.np/P.cpd/100);

    uint64 mask = auxstruct::lightconemask(lcn);
    double vunits = ReadState.VelZSpace_to_kms/ReadState.VelZSpace_to_Canonical;  // Code to km/s

    integer3 ij(slab,0,0);
    uint64_t slabtotal = 0;
    uint64_t slabtotalsub = 0;
    uint64_t doubletagged = 0;
    OutputLightConeSetup.Stop();
    OutputLightConeSearch.Start();

    #pragma omp parallel for schedule(dynamic,1) reduction(+:slabtotal) reduction (+:slabtotalsub) reduction(+:doubletagged)
    for (int y = 0; y < CP->cpd; y ++) {
        integer3 ijk = ij; ijk.y = y;

        PencilAccum<RVfloat>   *pLightConeRV   =   LightConeRV.StartPencil(ijk.y);
        PencilAccum<TaggedPID> *pLightConePIDs = LightConePIDs.StartPencil(ijk.y);
        PencilAccum<unsigned int> *pLightConeHealPix = LightConeHealPix.StartPencil(ijk.y);

        for (ijk.z = node_z_start; ijk.z < node_z_start + node_z_size; ijk.z++) {
            // Check if the cell center is in the lightcone, with some wiggle room
            double3 cc = CP->CellCenter(ijk);
            if(!LC.isCellInLightCone(cc)) continue;  // Skip the rest if too far from the region

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
                    posstruct poscopy = c.pos[p];  // Need a copy, since it will be changed
                    if (LC.isParticleInLightCone(cc, poscopy, vel, acc[p])) { 
                        // Yes, it's in the light cone.  pos and vel were updated.
                        double3 pos = cc+poscopy;    // This is now a global position in double precision

                        pLightConeHealPix->append(LC.healpixel(pos));  // We're outputting all particles for this
                        slabtotal++;

                        if(c.aux[p].is_taggable() or P.OutputFullLightCones){
                            // Going to output; pack the density in the aux
                            c.aux[p].set_compressed_density(acc[p].w);
                            
                            // These output routines take global positions and velocities in km/s
                            pLightConePIDs->append(TaggedPID(c.aux[p]));
                            pLightConeRV->append(RVfloat(pos.x, pos.y, pos.z, vel.x * vunits, vel.y * vunits, vel.z * vunits));
                            slabtotalsub++;

                        }

                        // TODO: For now, we're going to look for particles that get tagged twice.
                        // But maybe we'll find that they are very few, in which case we might stop tagging.
                        /*
                        if (c.aux[p].lightconedone(mask)) {
                            doubletagged++;
                            STDLOG(1,"Double tag: (%6.4f %6.4f %6.4f) = %10.7f vs %10.7f %10.7f\n",
                                pos.x, pos.y, pos.z, (pos-LC.origin).norm(), LC.rmin, LC.rmax);
                        }
                        */
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

    OutputLightConeSearch.Stop();
    OutputLightConeTeardown.Start();

    STDLOG(2,"Lightcone %d opened %d cells and found %d particles (%d subsampled) in slab %d.  %d double tagged\n",
            lcn,slabtotalcell,slabtotal,slabtotalsub,slab, doubletagged);
    WriteState.np_lightcone += slabtotal;
    if(slabtotal) {
        // This will create the directory if it doesn't exist (and is parallel safe)
        static int made_dir = 0;
        if (!made_dir) {
            made_dir = 1;
            char dir[32];
            sprintf(dir, "Step%04d", ReadState.FullStepNumber);
            CreateSubDirectory(P.LightConeDirectory, dir);
        }

        SlabType lcrvtype = (SlabType)((int)(LightCone0RV + 3*lcn));
        SlabType lcpidtype = (SlabType)((int)(LightCone0PID + 3*lcn));
        SlabType lchealtype = (SlabType)((int)(LightCone0Heal + 3*lcn));

        SB->AllocateSpecificSize(lcrvtype, slab, LightConeRV.get_slab_bytes());
        SB->AllocateSpecificSize(lcpidtype, slab, LightConePIDs.get_slab_bytes());
        SB->AllocateSpecificSize(lchealtype, slab, LightConeHealPix.get_slab_bytes());

        // TODO: we might want to parallelize the copies themselves, but that may
        // require nested parallelism
        #pragma omp parallel sections
        {
            #pragma omp section
            LightConeRV.copy_to_ptr((RVfloat *)SB->GetSlabPtr(lcrvtype, slab));
            #pragma omp section
            LightConePIDs.copy_to_ptr((TaggedPID *)SB->GetSlabPtr(lcpidtype, slab));
            #pragma omp section
            LightConeHealPix.copy_to_ptr((unsigned int *)SB->GetSlabPtr(lchealtype, slab));
        }

        unsigned int *arenaptr = (unsigned int *) SB->GetSlabPtr(lchealtype, slab);
        OutputLightConeSortHealpix.Start();
        // mm.mmsort(arenaptr, out, slabtotal, UINT32_MAX, 1<<20, 32);
        // tbb::parallel_sort(arenaptr, arenaptr+slabtotal);
        // std::sort(arenaptr, arenaptr+slabtotal);
        // __gnu_parallel::sort(arenaptr, arenaptr+slabtotal);
        // ips4o::parallel::sort(L2.d2_active, L2.d2_active+size, omp_get_max_threads());
        OutputLightConeSortHealpix.Stop();

        SB->StoreArenaNonBlocking(lcrvtype, slab);
        SB->StoreArenaNonBlocking(lcpidtype, slab);
        SB->StoreArenaNonBlocking(lchealtype, slab);
    }

    OutputLightConeFreeSlabAccum.Start();
    LightConeRV.destroy();
    LightConePIDs.destroy();
    LightConeHealPix.destroy();
    OutputLightConeFreeSlabAccum.Stop();

    OutputLightConeTeardown.Stop();

    return slabtotal;
}
