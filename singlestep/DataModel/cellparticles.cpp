/*
 * This class provides access to particle data in a cell-oriented manner.
 * This is done in a way such that in-order access to the particles should
 * still be easily inlinable in many cases.
 *
 * This class is a global singleton. The name of its instantiation is CP.
 */

#ifndef INCLUDE_PARTICLES
#define INCLUDE_PARTICLES

#include "particlestruct.cpp"


class CellParticles : public grid {
public:
    SlabBuffer *SB;

    CellParticles(int cpd, SlabBuffer *_SB):
            grid(cpd) { SB = _SB; }
    ~CellParticles(void) { }

// ==================================================================
// We have to wrap the 2-d cell index into a linear index.

    inline int _CellID(int y, int z) { return y*cpd+z; }
            // This function does not wrap
    inline int CellID(integer3 xyz) {
        xyz = WrapCell(xyz);
        return _CellID(xyz.y,xyz.z);
    }
    inline int CellID(int y, int z) {
        integer3 xyz = WrapCell(0,y,z);
        return _CellID(xyz.y,xyz.z);
    }
    // These return a linearized ID number for the cell order in these slabs
    // These will wrap; however the _CellID function does not.

// ==================================================================
// Here are routines to be given a pointer to the cellinfo for a given cell,
// or the number of particles in a cell.

    inline cellinfo *_CellInfo(integer3 xyz) {
        // Assumes a wrapped cell
        return ((cellinfo *)SB->GetSlabPtr(CellInfoSlab,xyz.x)) + CellID(xyz);
    }
    inline cellinfo *CellInfo(integer3 xyz) {
        return _CellInfo(WrapCell(xyz)); 
    }
    inline cellinfo *CellInfo(int x, int y, int z) {
        return _CellInfo(WrapCell(x,y,z));
    }
    // Return a pointer to the CellInfo for this cell.
    // Wrap the ID, look up the cellinfo, return pointer to it

    inline int NumberParticle(int x, int y, int z) {
        return CellInfo(x,y,z)->count;
    }
    int NumberParticle(integer3 xyz) {
        return NumberParticle(xyz.x, xyz.y, xyz.z);
    }
    // Return the number of particles in the cell
    // Wrap the ID, look up the cellinfo, return the number of particles



    inline cellinfo *_MergeCellInfo(integer3 xyz) {
        // Assumes a wrapped cell
        return ((cellinfo *)SB->GetSlabPtr(MergeCellInfoSlab,xyz.x)) + CellID(xyz);
    }
    inline cellinfo *MergeCellInfo(integer3 xyz) {
        return _MergeCellInfo(WrapCell(xyz)); 
    }
    inline cellinfo *MergeCellInfo(int x, int y, int z) {
        return _MergeCellInfo(WrapCell(x,y,z));
    }
    // Return a pointer to the MergeCellInfo for this cell.
    // Wrap the ID, look up the cellinfo, return pointer to it

    inline int MergeNumberParticle(int x, int y, int z) {
        return MergeCellInfo(x,y,z)->count;
    }
    inline int MergeNumberParticle(integer3 xyz) {
        return MergeNumberParticle(xyz.x, xyz.y, xyz.z);
    }
    // Return the number of particles in the merged cell
    // Wrap the ID, look up the cellinfo, return the number of particles


    inline cellinfo *_InsertCellInfo(integer3 xyz) {
        // Assumes a wrapped cell
        return ((cellinfo *)SB->GetSlabPtr(InsertCellInfoSlab,xyz.x)) + CellID(xyz);
    }
    inline cellinfo *InsertCellInfo(integer3 xyz) {
        return _InsertCellInfo(WrapCell(xyz)); 
    }
    inline cellinfo *InsertCellInfo(int x, int y, int z) {
        return _InsertCellInfo(WrapCell(x,y,z));
    }
    // Return a pointer to the InsertCellInfo for this cell.
    // Wrap the ID, look up the cellinfo, return pointer to it


// ==================================================================
// We have a Cell class to pass the cellinfo, pos, vel, and aux pointers
// for a given cell in one item.  Here's the routine to get that.

    inline Cell _GetCell(integer3 xyz) {
        // Assumes a wrapped cell
        Cell c;
        c.ijk = xyz;
        c.ci = ((cellinfo *)SB->GetSlabPtr(CellInfoSlab,xyz.x))+CellID(xyz);
        c.pos = (posstruct *)SB->GetSlabPtr(PosSlab,xyz.x)+c.ci->startindex;
        c.vel = (velstruct *)SB->GetSlabPtr(VelSlab,xyz.x)+c.ci->startindex;
        c.aux = (auxstruct *)SB->GetSlabPtr(AuxSlab,xyz.x)+c.ci->startindex;
        if(SB->IsSlabPresent(AccSlab,xyz.x))
            c.acc = (accstruct *) SB->GetSlabPtr(AccSlab,xyz.x)+c.ci->startindex;
        return c;
    }
    inline Cell GetCell(integer3 xyz) {
        return _GetCell(WrapCell(xyz));
    }
    inline Cell GetCell(int x, int y, int z) {
        return _GetCell(WrapCell(x,y,z));
    }
    // These fill an instance of class Cell with the info for the requested cell
    // from the original slabs.  We have to pass a pointer to avoid risks
    // of memory leaks from user error.


    inline Cell _GetMergeCell(integer3 xyz) {
        // Assumes a wrapped cell
        Cell c;
        c.ijk = xyz;
        c.ci = ((cellinfo *)SB->GetSlabPtr(MergeCellInfoSlab,xyz.x))+CellID(xyz);
        c.pos = (posstruct *)SB->GetSlabPtr(MergePosSlab,xyz.x)+c.ci->startindex;
        c.vel = (velstruct *)SB->GetSlabPtr(MergeVelSlab,xyz.x)+c.ci->startindex;
        c.aux = (auxstruct *)SB->GetSlabPtr(MergeAuxSlab,xyz.x)+c.ci->startindex;
        c.acc = NULL;
        return c;
    }
    inline Cell GetMergeCell(integer3 xyz) {
        return _GetMergeCell(WrapCell(xyz));
    }
    inline Cell GetMergeCell(int x, int y, int z) {
        return _GetMergeCell(WrapCell(x,y,z));
    }
    // These fill an instance of class Cell with the info for the requested cell
    // from the merged slabs

// ==================================================================
// Here are routines to get the positions, velocities, auxillary, 
// and acceleration lists separately.


    inline char *_CellPtr(int type, size_t size, integer3 xyz) {
        // Assumes a wrapped cell
        return SB->GetSlabPtr(type,xyz.x) + size*_CellInfo(xyz)->startindex;
    }
    inline char *_MergeCellPtr(int type, size_t size, integer3 xyz) {
        // Assumes a wrapped cell
        return SB->GetSlabPtr(type,xyz.x) + size*_MergeCellInfo(xyz)->startindex;
    }


    inline posstruct *PosCell(integer3 xyz) {
        return (posstruct *)_CellPtr(PosSlab,sizeof(posstruct),WrapCell(xyz));
    }
    inline posstruct *PosCell(int x, int y, int z) {
        return (posstruct *)_CellPtr(PosSlab,sizeof(posstruct),WrapCell(x,y,z));
    }
    
    // These return three separate pointers because x,y,z for a PosXYZ cell is not contiguous
    inline List3<FLOAT> PosXYZCell(integer3 xyz) {
        return PosXYZCell(xyz.x, xyz.y, xyz.z);
    }

    inline List3<FLOAT> PosXYZCell(int x, int y, int z) {
        // This version arranges all slab Xs before slab Ys, etc
        List3<FLOAT> posxyz;
        uint64 Nslab = SS->size(x);
        posxyz.X = (FLOAT *) _CellPtr(PosXYZSlab, sizeof(FLOAT), WrapCell(x,y,z));
        posxyz.Y = posxyz.X + Nslab;  // all Xs are followed by all Ys in a given slab
        posxyz.Z = posxyz.Y + Nslab;
        posxyz.N = (uint64) NumberParticle(x,y,z);
        return posxyz;
    }
    
    // Return the pointer to the start of the position list for this cell.
    // I.e., wrap the id, look up the cellinfo, and offset from PosSlab

    inline velstruct *VelCell(integer3 xyz) {
        return (velstruct *)_CellPtr(VelSlab,sizeof(velstruct),WrapCell(xyz));
    }
    inline velstruct *VelCell(int x, int y, int z) {
        return (velstruct *)_CellPtr(VelSlab,sizeof(velstruct),WrapCell(x,y,z));
    }
    // Return the pointer to the start of the velocity list for this cell.
    // I.e., wrap the id, look up the cellinfo, and offset from VelSlab

    inline auxstruct *AuxCell(integer3 xyz) {
        return (auxstruct *)_CellPtr(AuxSlab,sizeof(auxstruct),WrapCell(xyz));
    }
    inline auxstruct *AuxCell(int x, int y, int z) {
        return (auxstruct *)_CellPtr(AuxSlab,sizeof(auxstruct),WrapCell(x,y,z));
    }
    // Return the pointer to the start of the auxillary list for this cell.
    // I.e., wrap the id, look up the cellinfo, and offset from AuxSlab



    // TODO: Perhaps we no longer use AccCell?
    inline accstruct *AccCell(integer3 xyz) {
        return (accstruct *)_CellPtr(AccSlab,sizeof(accstruct),WrapCell(xyz));
    }
    inline accstruct *NearAccCell(integer3 xyz) {
        return (accstruct *)_CellPtr(AccSlab,sizeof(accstruct),WrapCell(xyz));
    }
    inline accstruct *AccCell(int x, int y, int z) {
        return (accstruct *)_CellPtr(AccSlab,sizeof(accstruct),WrapCell(x,y,z));
    }
    inline accstruct *NearAccCell(int x, int y, int z) {
        return (accstruct *)_CellPtr(AccSlab,sizeof(accstruct),WrapCell(x,y,z));
    }

    // Return the pointer to the start of the acceleration list for this cell.
    // I.e., wrap the id, look up the cellinfo, and offset from AccSlab



    inline posstruct *MergePosCell(integer3 xyz) {
        return (posstruct *)_MergeCellPtr(MergePosSlab,sizeof(posstruct),WrapCell(xyz));
    }
    inline posstruct *MergePosCell(int x, int y, int z) {
        return (posstruct *)_MergeCellPtr(MergePosSlab,sizeof(posstruct),WrapCell(x,y,z));
    }
    // Return the pointer to the start of the position list for this cell.
    // I.e., wrap the id, look up the cellinfo, and offset from PosSlab

    inline velstruct *MergeVelCell(integer3 xyz) {
        return (velstruct *)_MergeCellPtr(MergeVelSlab,sizeof(velstruct),WrapCell(xyz));
    }
    inline velstruct *MergeVelCell(int x, int y, int z) {
        return (velstruct *)_MergeCellPtr(MergeVelSlab,sizeof(velstruct),WrapCell(x,y,z));
    }
    // Return the pointer to the start of the velocity list for this cell.
    // I.e., wrap the id, look up the cellinfo, and offset from VelSlab

    inline auxstruct *MergeAuxCell(integer3 xyz) {
        return (auxstruct *)_MergeCellPtr(MergeAuxSlab,sizeof(auxstruct),WrapCell(xyz));
    }
    inline auxstruct *MergeAuxCell(int x, int y, int z) {
        return (auxstruct *)_MergeCellPtr(MergeAuxSlab,sizeof(auxstruct),WrapCell(x,y,z));
    }
    // Return the pointer to the start of the auxillary list for this cell.
    // I.e., wrap the id, look up the cellinfo, and offset from AuxSlab

};

#endif    // INCLUDE_PARTICLES
