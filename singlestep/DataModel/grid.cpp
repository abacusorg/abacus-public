/* The grid class provides the periodic geometry of the cells and 
the positions (in simulation units).  It is inherited by a number
of other classes.

Cells are indexed as an integer3.  Cell indices are not required to be
periodically wrapped; the grid class will provide the periodic wrapping.

Similarly, positions need not be periodically wrapped in general, but
this grid code will provide the wrapping.  It also defines how cells
and positions are related.

The mapping of Y-Z cell index to a linear cellID for the slabs is the 
responsibility of inherited classes.  It is a requirement that we never
use a linear mapping of the 3-d cells; this will spill a 32-bit int and 
breaks the slab geometry.

All functions that take integer3 as an input should also take 3 integers
x,y,z.  At present, we do not guarantee that wrappings more than 1 box
will function properly.

We will eventually be using cell-centered positions.  These are known as
'local' positions.  They may be single precision.  The box-centered positions
are known as 'global' positions.  These probably should be double precision.

*/

#ifndef INCLUDE_GRID
#define INCLUDE_GRID

#define BOXSIZE  1.0	// Size of box in simulation units
#define HALFBOXSIZE  (BOXSIZE/2.0)   // Half the size of box in simulation units
// Code units have the box at -0.5 to +0.5

class grid {
public:
    // grid(void);
    // grid(int cpd);

    int cpd;		// Must be odd
    int cpdhalf;	// CPD/2 round down
    double invcpd;	// Size of a cell in sim units
    double invcpd3;	// Volume of a cell in sim units
    double halfinvcpd;	// Half the size of a cell in sim units

    grid(int _cpd) {
	assert(_cpd%2==1);  	// check it is odd 
	cpd     = _cpd;
	cpdhalf = (cpd-1)/2;
	invcpd  = BOXSIZE/((double) cpd);   
	invcpd3 = invcpd*invcpd*invcpd;
	halfinvcpd = 0.5*invcpd;
    }



    int IsValidCell(integer3 ijk) {
        if (   ijk.x>0 && ijk.x<cpd
	    && ijk.y>0 && ijk.y<cpd
	    && ijk.z>0 && ijk.z<cpd) return 1; else return 0;
    }
    // Return 0 if a cell index is not in the primary zone; 1 if so.

    int IsValidSlab(integer3 ijk) {
        if (ijk.x>0 && ijk.x<cpd) return 1; else return 0;
    }
    int IsValidSlab(int i) {
        if (i>0 && i<cpd) return 1; else return 0;
    }
    // Return 0 if a slab index is not in the primary zone; 1 if so.


    inline integer3 WrapCell(int i, int j, int k) {
	integer3 c(i,j,k);
        while(c.x>=cpd) c.x -= cpd; while(c.x<0) c.x += cpd;
        while(c.y>=cpd) c.y -= cpd; while(c.y<0) c.y += cpd;
        while(c.z>=cpd) c.z -= cpd; while(c.z<0) c.z += cpd;
    //assertf(c.x < cpd && c.y < cpd && c.z < cpd, "Wrapped cell index (%d, %d, %d) >= cpd (%d)\n", c.x, c.y, c.z, cpd);
    //assertf(c.x >= 0 && c.y >= 0 && c.z >= 0, "Wrapped cell index (%d, %d, %d) < 0\n", c.x, c.y, c.z);
	return c;
    }
    inline integer3 WrapCell(integer3 ijk) {
	return WrapCell(ijk.x, ijk.y, ijk.z);
    }
    // These wrap the given cell index to the primary zone.

    inline int WrapSlab(int i) {
        while (i>=cpd) i-=cpd; while (i<0) i+=cpd;
	return i;
    }
    int WrapSlab(integer3 ijk) {
        return WrapSlab(ijk.x);
    }
    // These return the slab, wrapped to the primary zone.


    inline double3 CellCenter(int i, int j, int k) {
	double3 cc;
	cc.x = (i-cpdhalf)*invcpd;
	cc.y = (j-cpdhalf)*invcpd;
	cc.z = (k-cpdhalf)*invcpd;
	return cc;
    }
    inline double3 CellCenter(integer3 ijk) {
	double3 cc;
	cc.x = (ijk.x-cpdhalf)*invcpd;
	cc.y = (ijk.y-cpdhalf)*invcpd;
	cc.z = (ijk.z-cpdhalf)*invcpd;
	return cc;
    }
    // Returns the center of the cell in global coords.  
    // Cells are centered such that the first cell is -0.5 .. -0.5+invcpd
    // This does *not* wrap to the primary zone.

    inline double3 WrapCellCenter(int i, int j, int k) {
        return CellCenter(WrapCell(i, j, k));
    }
    inline double3 WrapCellCenter(integer3 ijk) {
        return CellCenter(WrapCell(ijk));
    }
    // Returns the center of the cell in global coords.  
    // This does wrap to the primary zone.
    // This is what we probably want for output.

#ifdef GLOBALPOS
    double3 LocalCellCenter(integer3 ijk) { return CellCenter(ijk); }
    double3 LocalCellCenter(int i, int j, int k) { return CellCenter(i,j,k); }
    // For box-centered, this is just the CellCenter.
#else
    inline double3 LocalCellCenter(integer3 ijk) { return double3(0.0,0.0,0.0); }
    inline double3 LocalCellCenter(int i, int j, int k) { return double3(0.0,0.0,0.0); }
    // For cell-centered coords, this is trivial.
#endif
    // The cell center in possibly local coords.


    inline integer3 Position2Cell(float3 p) {
    	integer3 c;
    	c.x = floor((p.x+HALFBOXSIZE)*cpd);
    	c.y = floor((p.y+HALFBOXSIZE)*cpd);
    	c.z = floor((p.z+HALFBOXSIZE)*cpd);
    	return c;
    }
    inline integer3 Position2Cell(double3 p) {
    	integer3 c;
    	c.x = floor((p.x+HALFBOXSIZE)*cpd);
    	c.y = floor((p.y+HALFBOXSIZE)*cpd);
    	c.z = floor((p.z+HALFBOXSIZE)*cpd);
    	return c;
    }
    // This bins a position to a cell with no wrapping
    // Probably only the double precision version should be used.
    // Replaces position2xyz()

    inline integer3 WrapPosition2Cell(float3 *p) {
    	integer3 c = Position2Cell(*p);
    	while (c.x<0) { c.x+=cpd; p->x+=BOXSIZE; }
    	while (c.y<0) { c.y+=cpd; p->y+=BOXSIZE; }
    	while (c.z<0) { c.z+=cpd; p->z+=BOXSIZE; }
    	while (c.x>=cpd) { c.x-=cpd; p->x-=BOXSIZE; }
    	while (c.y>=cpd) { c.y-=cpd; p->y-=BOXSIZE; }
    	while (c.z>=cpd) { c.z-=cpd; p->z-=BOXSIZE; }
    	return c;
    }
    inline integer3 WrapPosition2Cell(double3 *p) {
    	integer3 c = Position2Cell(*p);
    	while (c.x<0) { c.x+=cpd; p->x+=BOXSIZE; }
    	while (c.y<0) { c.y+=cpd; p->y+=BOXSIZE; }
    	while (c.z<0) { c.z+=cpd; p->z+=BOXSIZE; }
    	while (c.x>=cpd) { c.x-=cpd; p->x-=BOXSIZE; }
    	while (c.y>=cpd) { c.y-=cpd; p->y-=BOXSIZE; }
    	while (c.z>=cpd) { c.z-=cpd; p->z-=BOXSIZE; }
    	return c;
    }
    // This bins a position to a cell, then wraps both to the primary zone.
    // The position wrap is done in place.
    // This is the correct way to initialize particles from a global coord.

    inline integer3 LocalPosition2Cell(double3 *p) {
    	integer3 c = Position2Cell(*p);
    	*p = *p-CellCenter(c);
    	return WrapCell(c);
    }
    // This converts a global position, in-place, to a cell-centered position
    // and returns the cell.  Global positions must be double precision; 
    // after this is called, one could reduce to single precision.
    // This is the correct way to initialize particles into cell-centered positions.


};

#endif //  INCLUDE_GRID
