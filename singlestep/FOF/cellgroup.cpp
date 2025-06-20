// Copyright 2012-2025 The Abacus Developers
// SPDX-License-Identifier: GPL-3.0-or-later

/** \file We create a CellGroup for every multiplets and every boundary singlet.
We hide the bounding box test and status bits in the upper 8 bits of the multiplicity
This limits cell groups to 16M particles, which is safe
*/

#define XM_BIT (1<<24)
#define YM_BIT (1<<25)
#define ZM_BIT (1<<26)
#define XP_BIT (1<<27)
#define YP_BIT (1<<28)
#define ZP_BIT (1<<29)
#define DEFER  (1<<30)
#define CLOSED_BIT (1<<31)
    // User: edgebit should be 1<<24 for x_min (i.e., x<min)
    //                         1<<25 for y_min
    //                         1<<26 for z_min
    //                         1<<27 for x_max
    //                         1<<28 for y_max
    //                         1<<29 for z_max

/** The CellGroup class describes a FOF group within a single cell.

We hide information about whether the group is close to a boundary
in the highest 7 bits of the size.
*/

class CellGroup {
  public:
    int start;  ///< The starting index of the particles, zero-indexed in the cell
    int n;      ///< The group multiplicity.  
         ///< Particles are at locations [start, start+n)

    CellGroup(FOFgroup &g, FOFloat boundary) {
        start = g.start;
        n = g.n;
        // The BoundingBoxes are in code units
        if (g.BBmin.x<-boundary) n|= XM_BIT;
        if (g.BBmin.y<-boundary) n|= YM_BIT;
        if (g.BBmin.z<-boundary) n|= ZM_BIT;
        if (g.BBmax.x> boundary) n|= XP_BIT;
        if (g.BBmax.y> boundary) n|= YP_BIT;
        if (g.BBmax.z> boundary) n|= ZP_BIT;
        return;
    }
    CellGroup(int particlenum, posstruct &pos, FOFloat boundary) {
	// This is the constructor to load up a singlet
        start = particlenum;
        n = 1;
        if (pos.x> boundary) n|= XP_BIT;
        if (pos.y> boundary) n|= YP_BIT;
        if (pos.z> boundary) n|= ZP_BIT;
        if (pos.x<-boundary) n|= XM_BIT;
        if (pos.y<-boundary) n|= YM_BIT;
        if (pos.z<-boundary) n|= ZM_BIT;
        return;
    }

    void close_group() { n |= 0x80000000; return; }
    void defer_group() { n |= 0x40000000; return; }
    void clear_deferral() { n &= 0xbfffffff; return; }
    
    bool is_open() { return (n & 0xc0000000)==0; }
        // Open means that are neither closed nor deferred

    int size() { return n&(0x00ffffff); }
    int test(int edgebit) {
        // Returns !=0 if the bit is set, 0 if false
        return n & edgebit;
    }
};

bool test_particle(posstruct pos, int edgebit, FOFloat boundary) {
    // Given a position, decide if it's close to a boundary.
    if (edgebit==XP_BIT) return (pos.x> boundary);
    if (edgebit==YP_BIT) return (pos.y> boundary);
    if (edgebit==ZP_BIT) return (pos.z> boundary);
    if (edgebit==XM_BIT) return (pos.x<-boundary);
    if (edgebit==YM_BIT) return (pos.y<-boundary);
    if (edgebit==ZM_BIT) return (pos.z<-boundary);
    return false;
};

