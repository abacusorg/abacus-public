/*
 * LightCones.cpp
 *      Identify particles in the lightcone for a given set of origins and output them.
 */

// Number of healpix weights/maps, besides the density map
#define N_HEAL_WEIGHTS 3
#define LCTOLERANCE 0.0

#define c_kms 299792.0
#define etaktoHMpc (c_kms/100.)

#include "healpix_shortened.c"

template <size_t W>
class HealStruct_ {
public:
    static constexpr uint32_t HALO_MASK = 0x80000000;

    uint32_t id;
    uint32_t N;
    float weights[W];  // TODO: zero-length array may not be allowed

    bool operator<(const HealStruct_<W> &rhs) const {
        return id < rhs.id;
    }

    static constexpr size_t nweights() {
        return W;
    }

    inline HealStruct_<W> &operator+=(const HealStruct_<W> &rhs) {
        // Should have the same IDs, but don't check due to the expense
        N += rhs.N;
        for (size_t i = 0; i < W; i++) {
            weights[i] += rhs.weights[i];
        }
        return *this;
    }

    inline void set_halo() {
        // We'll keep halo particle pixels separate from field.
        // In post, we'll probably make (halo + field) and halo maps.
        id |= HALO_MASK;
    }
};

using HealStruct = HealStruct_<N_HEAL_WEIGHTS>;

// thread-level buffers for healstruct
using HealVec = AbacusVector<HealStruct>;
std::vector<HealVec> healbufs;

class LightCone {
  private:
    double rmin_tol2;       // Square of rmin-tol
    double rmax_tol2;       // Square of rmax+tol

  public:
    double3 origin;  // The observer location, in unit-cell units
    double rmin;     // The minimum distance to the light cone region (i.e., lower redshift)
    double rmax;     // The maximum distance to the light cone region (i.e., higher redshift)
    double tol;      // The tolerance for our searching
    FLOAT DeltaEtaKick;   // The total Delta(eta_Kick) for this step
    FLOAT driftfactor;
    int do_boundary;  // Whether to check a boundary layer of cells from the periodic images
    int nrep;  // The number of box repeats in each direction (0 = normal light cone)
    int lcid;  // Which lightcone this is


    LightCone(double3 origin, int do_boundary, int nrep, int lcid) :
        origin{origin},
        do_boundary{do_boundary},
        nrep{nrep},
        lcid{lcid}
        {

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

    inline void emplace_healstruct(HealVec &out, double3 pos, const velstruct &vel, const auxstruct &aux) const {
        int64_t pixel;
        pos -= origin;
        vec2pix_nest64(P.LCHealpixNside, reinterpret_cast<double *>(&pos), &pixel);

        // pos is already global
        FLOAT3 los = pos / pos.norm();
        FLOAT vel_los = vel.dot(los);

        if constexpr (HealStruct::nweights() == 1) {
            HealStruct &newval = out.emplace_back(static_cast<uint32_t>(pixel), static_cast<uint32_t>(1), vel_los);
            
            if (aux.is_L1_sticky()) {
                newval.set_halo();
            }
        } else if constexpr (HealStruct::nweights() == 3) {
            // Add the two transverse velocities as weights

            // Get phi_hat in cartesian coordinates
            FLOAT3 phi_hat = FLOAT3(0, 0, 1).cross(los);
            phi_hat /= phi_hat.norm();

            // Get theta_hat in cartesian coordinates
            FLOAT3 theta_hat = phi_hat.cross(los);

            FLOAT vel_phi = vel.dot(phi_hat);
            FLOAT vel_theta = vel.dot(theta_hat);

            HealStruct &newval = out.emplace_back(static_cast<uint32_t>(pixel), static_cast<uint32_t>(1), vel_los, vel_theta, vel_phi);

            if (aux.is_L1_sticky()) {
                newval.set_halo();
            }
        } else {
            throw std::runtime_error("Unsupported number of weights for HealStruct");
        }
    }

    inline int isCellInLightCone(double3 pos) const;
    inline int isParticleInLightCone(double3 cellcenter, double3 &dpos, velstruct &vel, const accstruct acc, double3 box_repeat_offset) const;

    static void WriteHeaderFiles(const fs::path &dir) {
        // Write "header_read" and "header_write" files.
        // This is because for light cones, we output continuously during the timestep,
        // so one may want to know cosmological quantities at both the start and end of the step.
        // This is intended to be called at the end of the timestep, which allows us to write
        // not just the output header but the full state(s).
        
        auto write_state_file = [&dir](const State& state, const std::string& suffix) {
            fs::path fn = dir / ("header_" + suffix);
            STDLOG(2, "Writing light cone {} file: {}\n", suffix, fn.string());
            
            FILE *headerfile = fopen(fn.c_str(), "w");
            if (headerfile == nullptr) {
                throw std::runtime_error("Failed to open header file for writing: " + fn.string());
            }
            fmt::print(headerfile, "{}", P.header().c_str());
            fmt::print(headerfile, "{}\n", state.get_state_string().c_str());
            fmt::print(headerfile, "OutputType = \"LightCone\"\n");
            
            std::string HealpixWeightScheme = "";
            if constexpr (HealStruct::nweights() == 1) {
                // LOS vel
                HealpixWeightScheme = "vel";
            } else if constexpr (HealStruct::nweights() == 3) {
                // LOS, theta, phi vel
                HealpixWeightScheme = "vel3";
            } else {
                static_assert(HealStruct::nweights() == 1 || HealStruct::nweights() == 3, "Unsupported number of weights for HealStruct");
            }
            fmt::print(headerfile, "HealpixWeightScheme = {}\n", HealpixWeightScheme);
            fclose(headerfile);
        };
        
        // Write both state files
        write_state_file(ReadState, "read");
        write_state_file(WriteState, "write");
    }
};

// Return whether a CellCenter is in the light cone, including some tolerance
inline int LightCone::isCellInLightCone(double3 pos) const {
    double r2 = (pos-origin).norm2();
    return (r2<=rmax_tol2) && (r2>=rmin_tol2);
}

// pos is in cell-centered coords; cellcenter is the position in the box; vel is in code units
// rmax is the maximum radius, which is the earlier time
// rmin is the minimum radius, which is the later time
//
// Returns 0 if not in LightCone, 1 if it is.
// Further, the position and velocity inputs will be adjusted to the epoch where the light cone is crossed.
inline int LightCone::isParticleInLightCone(double3 cellcenter, double3 &dpos, velstruct &vel, const accstruct acc, double3 box_repeat_offset) const {
    double r0 = (cellcenter-origin+dpos).norm();
    posstruct pos1 = static_cast<posstruct>(dpos) + vel*driftfactor;   // Take care to match the precision of Drift()
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
        // frac_step = 0 means r=rmax, =1 means r=rmin

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
    // and make the pos global
    dpos += vel*driftfactor*frac_step + cellcenter;

    if(this->do_boundary){
        // Check if the interpolated pos is in primary region
        if(!dpos.inrange<IntervalType::HalfOpen>(
            double3(-0.5, -0.5, -0.5) + box_repeat_offset,
            double3(0.5, 0.5, 0.5) + box_repeat_offset
            )) return 0;
    }

    vel += static_cast<FLOAT3>(acc)*DeltaEtaKick*frac_step;
    return 1;
}

std::vector<LightCone> LightCones;

int CountLightCones(){
    assertf(P.LCorigins.size() % 3 == 0, "LCorigins must be specified as a list of 3-tuples\n");

    int NLC = P.LCorigins.size() / 3;
    
    return NLC;
}

void InitializeLightCones(){
    // LightCone objects are persistent through the time step

    int NLC = CountLightCones();
    LightCones.reserve(NLC);
    for(int i = 0; i < NLC; i++){
        LightCones.emplace_back(
            double3(P.LCorigins[3*i], P.LCorigins[3*i+1], P.LCorigins[3*i+2])/P.BoxSize,
            P.LCCheckAcrossWrap,
            P.LCBoxRepeats,
            i
        );
    }

    STDLOG(2, "Initialized Light Cones with {:d} observers\n", LightCones.size());

#ifdef USE_LC_AUX_BITS
    assertf(LightCones.size() <= NUM_LC_AUX_BITS, "Parameter file requests {:d} light cones, but AUX data model supports only {:d}\n", LightCones.size(), NUM_LC_AUX_BITS);
#endif

    if(P.LCHealpixOutputSparseMap && !LightCones.empty()){
        size_t maxcell = ReadState.MaxCellSize;
        // TODO: if healbuf remains useful, track maxpencil in the merge
        // rather than assuming all cells in the pencil have size maxcell
        size_t maxpencil = maxcell * P.cpd * (2 * P.LCBoxRepeats + 1);

        STDLOG(0, "Allocating {:.3g} GB for light cone HealStruct temporary buffers\n",
            sizeof(HealStruct) * maxpencil * omp_get_max_threads() / 1e9
        );
        
        healbufs.resize(omp_get_max_threads());
        #pragma omp parallel
        healbufs[omp_get_thread_num()].reserve(maxpencil);
    }
}

void FinalizeLightCones(){
    if(P.LCHealpixOutputSparseMap && !LightCones.empty()){
        healbufs.clear();
    }

    LightCones.clear();
}

struct LCCellBounds {
    int xstart, xend, ystart, yend, zstart, zend;

    LCCellBounds(int slab, const LightCone &LC){
        if (LC.do_boundary) {
            // Search a boundary layer 1 cell deep.
            // We accept particles from these boundary cells
            // if their interpolated LC crossing is in the primary zone.
            // For y and z, we effectively extend the loop over [0,cpd) to [-1,cpd+1).
            // For x, when we encounter the x=0 and x=cpd-1 slabs,
            // we also consider their periodic images.
            // TODO: as an optimization, if the box face is adjacent to a box repeat,
            // we could skip the boundary layer for that face.

            if (slab == 0) {
                xstart = -1;
                xend = 1;
            } else if (slab == CP->cpd - 1) {
                xstart = CP->cpd - 1;
                xend = CP->cpd + 1;
            } else {
                xstart = slab;
                xend = slab + 1;
            }

            ystart = -1;
            yend = CP->cpd + 1;

            // The z boundaries may not be the box boundaries in 2D
            zstart = node_z_start == 0 ? -1 : node_z_start;
            zend = node_z_size == CP->cpd ? CP->cpd : node_z_start + node_z_size + 1;
        } else {
            xstart = slab;
            xend = slab + 1;
            ystart = 0;
            yend = CP->cpd;
            zstart = node_z_start;
            zend = node_z_start + node_z_size;
        }
    }
};

uint64_t count_cells_in_lightcone(const LightCone &LC, const LCCellBounds &bounds);
void fill_lightcone_slabaccums(
    const int slab,
    const LightCone &LC,
    const LCCellBounds &cell_bounds,
    SlabAccum<RVfloat> &LightConeRV,
    SlabAccum<TaggedPID> &LightConePIDs,
    SlabAccum<HealStruct> &LightConeHealPix,
    size_t &nheal,
    size_t &nrvpid
);
void output_lightcone_slabaccums(
    int slab,
    int lcid,
    SlabAccum<RVfloat> &LightConeRV,
    SlabAccum<TaggedPID> &LightConePIDs,
    SlabAccum<HealStruct> &LightConeHealPix);
void create_sparse_map(const HealStruct *healbuf, size_t Nheal, int slab, int lchealsparsetype);


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

size_t makeLightCone(int slab, size_t lcid){  //lcid = Light Cone ID
    // Use the same format for the lightcones as for the particle subsamples
    if (fabs(cosm->next.etaK-cosm->current.etaK)<1e-12) return 0;
          // Nothing to be done, so don't risk divide by zero.
    
    LCSetup.Start();
    STDLOG(4, "Making light cone {:d} for slab {:d}\n", lcid, slab);

    LightCone &LC = LightCones[lcid];
    LCCellBounds lc_cell_bounds(slab, LC);

    // Before we start allocating memory, do a quick pass to see if any cells are in the LC
    uint64_t slabtotalcell = count_cells_in_lightcone(LC, lc_cell_bounds);

    if(slabtotalcell == 0){
        LCSetup.Stop();
        return 0;
    }

    STDLOG(2,"Lightcone {:d} will open {:d} cells\n", LC.lcid, slabtotalcell);

    SlabAccum<RVfloat>   LightConeRV;     ///< The taggable subset in each lightcone.
    SlabAccum<TaggedPID> LightConePIDs;   ///< The PIDS of the taggable subset in each lightcone.
    SlabAccum<HealStruct> LightConeHealPix;   ///< The Healpix info of the particles in each lightcone.

    size_t nheal = 0, nrvpid = 0;
    fill_lightcone_slabaccums(slab, LC, lc_cell_bounds, LightConeRV, LightConePIDs, LightConeHealPix, nheal, nrvpid);
    WriteState.np_lightcone += nheal;
    
    if(nheal || nrvpid) {
        LCGatherOutput.Start();
        output_lightcone_slabaccums(slab, lcid, LightConeRV, LightConePIDs, LightConeHealPix);
        LCGatherOutput.Stop();
    }

    LCFreeSlabAccum.Start();
    LightConeRV.destroy();
    LightConePIDs.destroy();
    LightConeHealPix.destroy();
    LCFreeSlabAccum.Stop();

    return nheal;
}


uint64_t count_cells_in_lightcone(const LightCone &LC, const LCCellBounds &bounds){
    // This sextuply nested loop illustrates the structure of the light cone search:
    // for each periodic image, loop over each cell in the slab.
    // That's only 5 loops, but we pick up a 6th if we're checking boundary cells.

    uint64_t slabtotalcell = 0;
    #pragma omp parallel for schedule(static) reduction(+:slabtotalcell)
    for (int y = bounds.ystart; y < bounds.yend; y ++) {
        for (int z = bounds.zstart; z < bounds.zend; z++) {
            for(int x = bounds.xstart; x < bounds.xend; x++) {
                int xoff = x == -1 ? CP->cpd : (x == CP->cpd ? -CP->cpd : 0);
                double3 cc_primary = CP->CellCenter(x + xoff, y, z);
                for(int by = -LC.nrep; by <= LC.nrep; by++){
                    for(int bz = -LC.nrep; bz <= LC.nrep; bz++){
                        for(int bx = -LC.nrep; bx <= LC.nrep; bx++){
                            // unit box units
                            double3 box_repeat_offset = double3(bx, by, bz);

                            // Check if the cell center is in the lightcone, with some wiggle room
                            double3 cc_rep = cc_primary + box_repeat_offset;
                            if(LC.isCellInLightCone(cc_rep)) slabtotalcell++;
                        }
                    }
                }
            }
        }
    }
    return slabtotalcell;
}


void compress_healbuf(HealVec &healbuf){
    // Sort, then co-add towards the front, in place.
    // This should be single-threaded.

    ips4o::sort(healbuf.begin(), healbuf.end());

    size_t j = 0;  // where we're writing
    for(size_t i = 1; i < healbuf.size(); i++){
        if(healbuf[i].id == healbuf[j].id){
            healbuf[j] += healbuf[i];
        } else {
            j++;
            healbuf[j] = healbuf[i];
        }
    }

    healbuf.resize(j+1);  // this shouldn't reallocate
}


void fill_lightcone_slabaccums(
    const int slab,
    const LightCone &LC,
    const LCCellBounds &cell_bounds,
    SlabAccum<RVfloat> &LightConeRV,
    SlabAccum<TaggedPID> &LightConePIDs,
    SlabAccum<HealStruct> &LightConeHealPix,
    size_t &nheal,
    size_t &nrvpid
    ){

    // The size estimates are objects per slab.
    // The light cones output only taggables particles for RV and PID, 
    // but all particles for HealPix.
    // What fraction of a slab is the light cone?  It depends on geometry and the time step.
    // But for 1-2 Gpc boxes, it takes a few hundred time steps to cross the box.  We'll guess 1%.
    // That said, if the slab is tangent to the annulus, one can include >1% of the slab.
    // If we're checking boundary cells, we effectively have a slab that is bigger by 2 cells in y and z.
    int pad = LC.do_boundary ? 2 : 0;
    // Could consider zwidth > 1 to keep CPUs on different cache lines, but we may not touch the CellAccums very often
    // TODO: may want to consider smarter defaults here, including LC.nrep scaling
    if (P.LCOutputRVPID){
        LightConeRV.setup(  CP->cpd + pad, 1, P.np/P.cpd/node_z_size*(P.ParticleSubsampleA+P.ParticleSubsampleB)/10);
        LightConePIDs.setup(CP->cpd + pad, 1, P.np/P.cpd/node_z_size*(P.ParticleSubsampleA+P.ParticleSubsampleB)/10);
    }
    if (P.LCHealpixOutputSparseMap){
        LightConeHealPix.setup(CP->cpd + pad, 1, P.np/P.cpd/node_z_size/10);
    }

    uint64 mask = auxstruct::lightconemask(LC.lcid);
    double vunits = ReadState.VelZSpace_to_kms/ReadState.VelZSpace_to_Canonical;  // Code to km/s

    // With P.LCBoxRepeats, part of the goal is to put all the particles from all periodic images
    // (up to the specified number) in the same light cone files and thus the same coordinate system.
    // However, RVfloat (i.e. RVint) can only handle positions in the primary image. With 20 bits, RVint
    // has enough precision that even if we rescale the repeated volume to the primary one, it should
    // still be good enough for many purposes. For high-resolution simulations with many repeats, one
    // might want a different approach.
    // Downstream code will need to know the box repeats to convert the RVfloats to global coordinates.
    double rvfloat_scale = 1./(2*LC.nrep + 1);

    uint64_t doubletagged = 0;

    LCSetup.Stop();
    LCSearch.Start();

    // TODO: we probably don't get any reuse by putting the cell as the outer loop.
    // Can we make the repeats be the outer loop, so we can more easily skip images
    // that are too far away, or too close?
    // How do we parallelize in a way that checking each box repeat is only done once,
    // but we can parallelize over pencils? Ideally over all repeats at once? Tasks?
    // Let's parallelize over y-pencils x all repeats.
    // The valid repeats can be pre-computed for this timestep in InitLightCones()
    // and stored in an array.
    #pragma omp parallel for schedule(dynamic,1) reduction(+:nheal,nrvpid,doubletagged)
    for (int y = cell_bounds.ystart; y < cell_bounds.yend; y++) {
        HealVec *healbuf = nullptr;
        if(P.LCHealpixOutputSparseMap){
            healbuf = &healbufs[omp_get_thread_num()];
        }

        integer3 ijk(slab,0,0);
        integer3 boundary_offset(0,0,0);  // used for shifting the boundary layer of cells
        if (y == -1){
            ijk.y = CP->cpd - 1;
            boundary_offset.y = -CP->cpd;
        }
        else if (y == CP->cpd) {
            ijk.y = 0;
            boundary_offset.y = CP->cpd;
        }
        else {
            ijk.y = y;
            boundary_offset.y = 0;
        }

        PencilAccum<RVfloat>   *pLightConeRV = nullptr;
        PencilAccum<TaggedPID> *pLightConePIDs = nullptr;
        PencilAccum<HealStruct> *pLightConeHealPix = nullptr;

        if (P.LCOutputRVPID){
            pLightConeRV =   LightConeRV.StartPencil(y - cell_bounds.ystart);
            pLightConePIDs = LightConePIDs.StartPencil(y - cell_bounds.ystart);
        }
        if (P.LCHealpixOutputSparseMap) {
            pLightConeHealPix = LightConeHealPix.StartPencil(y - cell_bounds.ystart);
        }

        for (int z = cell_bounds.zstart; z < cell_bounds.zend; z++) {
            if (z == -1){
                ijk.z = CP->cpd - 1;
                boundary_offset.z = -CP->cpd;
            }
            else if (z == CP->cpd) {
                ijk.z = 0;
                boundary_offset.z = CP->cpd;
            }
            else {
                ijk.z = z;
                boundary_offset.z = 0;
            }

            for(int x = cell_bounds.xstart; x < cell_bounds.xend; x++){
                // Same cell, different offset
                boundary_offset.x = x == -1 ? CP->cpd : (x == CP->cpd ? -CP->cpd : 0);
                double3 cc_primary = CP->CellCenter(ijk + boundary_offset);

                for(int by = -LC.nrep; by <= LC.nrep; by++){
                    for(int bz = -LC.nrep; bz <= LC.nrep; bz++){
                        for(int bx = -LC.nrep; bx <= LC.nrep; bx++){
                            // unit box units
                            double3 box_repeat_offset = double3(bx, by, bz);

                            // Check if the cell center is in the lightcone, with some wiggle room
                            double3 cc = cc_primary + box_repeat_offset;
                            if(!LC.isCellInLightCone(cc)) continue;  // Skip the rest if too far from the region

                            Cell c = CP->GetCell(ijk);
                            accstruct *acc = CP->AccCell(ijk);
                            #ifdef GLOBALPOS
                                cc= 0*cc;
                            #endif

                            // STDLOG(4, "LC: Particles in current cell: {:d}\n", c.count());
                            for (int p=0;p<c.count();p++) {
                                if(!c.aux[p].lightconedone(mask) || boundary_offset != integer3(0.)){
                                    // This particle isn't already in the light cone,
                                    // or it interpolated across the periodic wrap to end up in the LC.
                                    // If it did that, then it's almost certainly about to get its LC aux bits reset,
                                    // but that hasn't happened yet.

                                    // Need to unkick by half
                                    // We're going to modify pos and vel, so make sure they're copies
                                    velstruct vel = c.vel[p] - static_cast<FLOAT3>(acc[p])*WriteState.FirstHalfEtaKick;
                                    double3 pos = c.pos[p];
                                    if (LC.isParticleInLightCone(cc, pos, vel, acc[p], box_repeat_offset)) {
                                        // Yes, it's in the light cone.  pos and vel were updated, and the pos made global.
                                        vel *= vunits;
                                        
                                        if (P.LCHealpixOutputSparseMap){
                                            // All particles are part of the healpix map
                                            // We're pushing HealStructs to a temporary cell-wise buffer,
                                            // then co-adding as many as we can before pushing them to the slab accum.
                                            // We expect to have a high hit rate, as many cells will project onto just a few pixels.
                                            LC.emplace_healstruct(*healbuf, pos, vel, c.aux[p]);
                                        }

                                        if(P.LCOutputRVPID && (c.aux[p].is_taggable() || P.OutputFullLightCones)){
                                            // Going to output; pack the density in the aux
                                            c.aux[p].set_compressed_density(acc[p].w);
                                            
                                            // These output routines take global positions and velocities in km/s
                                            pLightConePIDs->append(TaggedPID(c.aux[p]));
                                            pLightConeRV->append(RVfloat(
                                                pos.x * rvfloat_scale, pos.y * rvfloat_scale, pos.z * rvfloat_scale,
                                                vel.x, vel.y, vel.z
                                                ));
                                            nrvpid++;
                                        }

                                        // TODO: For now, we're going to look for particles that get tagged twice.
                                        // But maybe we'll find that they are very few, in which case we might stop tagging.
                                        /*
                                        if (c.aux[p].lightconedone(mask)) {
                                            doubletagged++;
                                            STDLOG(1,"Double tag: ({:6.4f} {:6.4f} {:6.4f}) = {:10.7f} vs {:10.7f} {:10.7f}\n",
                                                pos.x, pos.y, pos.z, (pos-LC.origin).norm(), LC.rmin, LC.rmax);
                                        }
                                        */
                                        c.aux[p].setlightconedone(mask);
                                    }
                                }
                            }  // loop over particles within cell
                        }  // loop over box x repeats
                    }  // loop over box z repeats
                }  // loop over box y repeats
            }  // loop over xoffset
        }  // loop over cell within pencil
        
        // We don't care about cell indexing, so we'll just use one "cell" per pencil

        if (P.LCOutputRVPID){
            pLightConeRV->FinishCell();
            pLightConeRV->FinishPencil();
            
            pLightConePIDs->FinishCell();
            pLightConePIDs->FinishPencil();
        }

        if (P.LCHealpixOutputSparseMap) {
            // compress the healbuf, dump to SlabAccum
            // TODO: we probably ought to adapt SlabAccum for this case rather than use
            // two levels of data structures
            if(!healbuf->empty()){
                compress_healbuf(*healbuf);
                for (const auto &heal : *healbuf) {
                    // This is pretty fast, it seems like <10% of the search time
                    pLightConeHealPix->append(heal);
                }
                nheal += healbuf->size();
                healbuf->clear();
            }

            pLightConeHealPix->FinishCell();            
            pLightConeHealPix->FinishPencil();
        }
    }  // loop over pencils within slab

    LCSearch.Stop();

    STDLOG(2,"Lightcone {:d} opened {:d} cells and found {:d} particles ({:d} subsampled) in slab {:d}.  {:d} double tagged\n",
            LC.lcid,nheal,nrvpid,slab, doubletagged);
}


void output_lightcone_slabaccums(
    int slab,
    int lcid,
    SlabAccum<RVfloat> &LightConeRV,  // TODO: would like this to be const
    SlabAccum<TaggedPID> &LightConePIDs,
    SlabAccum<HealStruct> &LightConeHealPix
    ){

    // This will create the directory if it doesn't exist (and is parallel safe)
    static int made_dir = 0;
    if (!made_dir) {
        made_dir = 1;
        fs::create_directory(P.LCDirectory / fmt::format("Step{:04d}", ReadState.FullStepNumber));
    }

    // TODO: make slabtypes an enum class, and add a function to get the slabtype for a given LC ID
    // Also add IsLCType to enum class
    SlabType lcrvtype = static_cast<SlabType>(static_cast<int>(LightCone0RV) + NumLCTypes*lcid);
    SlabType lcpidtype = static_cast<SlabType>(static_cast<int>(LightCone0PID) + NumLCTypes*lcid);
    SlabType lchealsparsetype = static_cast<SlabType>(static_cast<int>(LightCone0HealSparse) + NumLCTypes*lcid);

    if (P.LCOutputRVPID){
        SB->AllocateSpecificSize(lcrvtype, slab, LightConeRV.get_slab_bytes());
        SB->AllocateSpecificSize(lcpidtype, slab, LightConePIDs.get_slab_bytes());
    }

    // Copy HealStruct slab accum to contiguous temporary buffer, sort, build sparse map, append to disk.
    // Copy particle healpix IDs to arena if needed.
    // If only IDs are needed, then just copy straight into arena and sort.

    if(P.LCHealpixOutputSparseMap){
        LCCreateSparseMap.Start();
        size_t Nheal = LightConeHealPix.get_slab_size();
        HealStruct *healbuf = new HealStruct[Nheal];

        #pragma omp parallel
        #pragma omp single
        LightConeHealPix.copy_to_ptr(healbuf);

        LCSortHealpix.Start();
        ips4o::parallel::sort(healbuf, healbuf + Nheal);
        LCSortHealpix.Stop();

        create_sparse_map(healbuf, Nheal, slab, lchealsparsetype);

        delete[] healbuf;
        LCCreateSparseMap.Stop();
    }

    // TODO: is parallel sections playing nice with taskloop in copy_to_ptr?
    #pragma omp parallel sections
    {
        #pragma omp section
        if (P.LCOutputRVPID) LightConeRV.copy_to_ptr((RVfloat *) SB->GetSlabPtr(lcrvtype, slab));
        #pragma omp section
        if (P.LCOutputRVPID) LightConePIDs.copy_to_ptr((TaggedPID *) SB->GetSlabPtr(lcpidtype, slab));
    }

    if (P.LCOutputRVPID) {
        SB->StoreArenaNonBlocking(lcrvtype, slab);
        SB->StoreArenaNonBlocking(lcpidtype, slab);
    }
    if(P.LCHealpixOutputSparseMap) {
        SB->StoreArenaNonBlocking(lchealsparsetype, slab);
    }
}


void create_sparse_map(const HealStruct *in, size_t Nheal, int slab, int lchealsparsetype){
    // We have an array of {pixel_ID, count, {weights, ...}} structs, sorted by pixel_ID.
    // We want to create a new array where all pixels with the same ID are combined into one entry,
    // with the weights added together.
    // We'll do it out-of-place, into an arena.

    // Our parallelization strategy will be:
    // - Create list of edges, i.e. indices in the input array.
    //   - Each thread can first count its edges to compute a write index in a shared array
    // - Parallelize the writes (statically?) over these indices
    //   - The index of the index is the output index
    //   - The difference with the next index is the input length

    size_t Nthread = omp_get_max_threads();

    // Count how many edges each thread sees
    padded<size_t> thread_offsets[Nthread+1];
    thread_offsets[0] = 0;  // thread_offsets[t+1] holds the counts for thread t
    #pragma omp parallel
    {
        int t = omp_get_thread_num() + 1;
        thread_offsets[t] = 0;

        #pragma omp for schedule(static)
        for(size_t i = 1; i < Nheal; i++){
            if(in[i].id != in[i-1].id) {
                thread_offsets[t]++;
            }
        }
    }

    // prefix sum
    for(size_t i = 1; i <= Nthread; i++){
        thread_offsets[i] += thread_offsets[i-1];
    }
    size_t Nunique = thread_offsets[Nthread];

    // record the input locations of each new pixel
    AbacusVector<size_t> edge_indices(Nunique+1);
    #pragma omp parallel
    {
        size_t j = thread_offsets[omp_get_thread_num()];

        #pragma omp for schedule(static)
        for(size_t i = 1; i < Nheal; i++){
            if(in[i].id != in[i-1].id) {
                edge_indices[j++] = i;
            }
        }
    }
    edge_indices[0] = 0;
    edge_indices[Nunique] = Nheal;

    // Allocate the sparse map
    HealStruct *out = (HealStruct *) SB->AllocateSpecificSize(lchealsparsetype, slab, sizeof(HealStruct) * Nunique);
    // (void) out;

    // TODO: to reduce peak memory usage during post-processing,
    // we could write one array per weight field (plus 1 array for the pixel indexing)

    // Write/compress the pixels, parallelizing over writes
    // There's a potential for load imbalance here, as we haven't tracked how
    // reads each write is doing.

    // Use fewer threads when writing to fewer pixels
    size_t nthread_write = std::max(
        std::min(Nunique * sizeof(HealStruct) / CACHE_LINE_SIZE, Nthread),
        static_cast<size_t>(1)
    );

    #pragma omp parallel for schedule(static) num_threads(nthread_write)
    for(size_t i = 0; i < Nunique; i++){
        size_t j = edge_indices[i];
        size_t jnext = edge_indices[i+1];

        out[i] = in[j];
        j++;

        for(; j < jnext; j++){
            out[i] += in[j];
        }
    }
}
