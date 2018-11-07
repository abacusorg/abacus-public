/* 
The cell-oriented structure of abacus allows for an efficient way
to compress the particle data.  A typical cell size is a few Mpc.
12 bits of precision will yield a part in 4000 precision, which is
about 1 kpc.  This is adequate for all reasonable output tasks.
These files are not intended to be adequate for dynamical restarts!

Abacus always outputs its velocities in redshift-space units, i.e.,
the comoving displacement that the particle would appear to be in
redshift space.  This avoids a need to choose a second set of units
for velocities.  We have v_proper = H*x_proper = H*a*x_comoving,
so v_zspace = v_proper*(1+z)/H(z).

Velocities rarely get above 6000 km/s, so a part in 4000 means 1-2
km/s precision, which is 10 kpc in zspace displacement.

So our standard pack14 class is 12 bits for each of the 6 phase space
components.  We add 5 bytes of ID numbers or group IDs, making a total
of 14 bytes.

We treat these 12 bits as signed and limit the range to +-2000.
This allows us to use the first byte as 0xff to indicate a special
interpretation of the contents as a cell.  Technically, we allow 
the 12-bit values to reach +-2030, i.e., mild overflow.

The global position is then (cellx + particlex/2000.0)*B, where B
is the box size divided by CPD.  

The global velocity is (particlevx/2000.0)*A*B, where A is the
scaling factor chosen when packing the velocities.  A must be
larger than the maximum of the velocities in ZSpace unit-cell units
(not unit box!).  We expect that good A's will be of order 30.  There is
not much value in making A smaller than 1, since this would mean more
resolution in velocity than position!

For the cells, we use the six 12-bit numbers to store the cell 
indicator, the xyz of the cell, the cpd of the problem, and the 
velocity scaling A.  This allows us to reconstruct into positions
and redshift-space displacements according to a unit box length.
We can accomodate cpd up to 4000, and ijk from -30 to cpd+30.

The 12-bit numbers are stored as unsigned.  For the particles, the
positions and velocities are signed quantities, so we add 2048 to
them before packing.

The intent is that with 5-byte IDs and CPD up to 4000, we should
be able to handle trillion particle runs.  If the IDs are group
ids, then we ought to be able to approach 1e13 particle runs with
this format.

Note that because we allow mild (1.5%) excursions beyond the cell
boundaries and the cell indices, it is not a guarantee that the
particles will be returned in the primary box volume.  The user may
wish to apply another periodic wrapping.

TODO:
In addition, we could use the 5 bytes of the id to store the box length
and hMpc as a bonus, so that each cell/particle set is self-contained
to get to full positions.  Or we could store the number of particles in
the cell.  At present, we're just storing 0.

*/

#ifndef PACKED_CPP
#define PACKED_CPP

#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include "threevector.hh"
#include "cell_header.h"




class pack14 {
  public:
    unsigned char c[14];

    pack14() { }
    ~pack14(void) { }

    // Routines to decide if this is a cell or a particle.
  protected:
    void make_cell() { 
        c[0] = 0xff;
    }

  public:
    int iscell() { 
            if (c[0]==0xff) return 1; else return 0;
    }
    int isparticle() { 
            if (c[0]==0xff) return 0; else return 1;
    }

    void print_char() {
         printf("Char binary:  %02x %02x %02x  %02x %02x %02x  %02x %02x %02x  %02x %02x %02x %02x %02x\n", 
         c[0], c[1], c[2], c[3], c[4], c[5], c[6], c[7], c[8], c[9], c[10], c[11], c[12], c[13]);
    }
    void print_short(SHORT s[6]) {
         printf("Short decimal: %d %d %d %d %d %d\n",
             s[0], s[1], s[2], s[3], s[4], s[5]);
         printf("Short binary: %03x %03x  %03x %03x  %03x %03x\n",
             s[0], s[1], s[2], s[3], s[4], s[5]);
        }
                 

  private:
    // Expand from 12-bits to 16-bits, so that the 
    // rest of the math is transparent.  These shorts
    // are signed, so we subtract 2048 from the 12-bit
    // unsigned quantity.
    void expand_to_short(SHORT s[6]) {
        // print_char();
        SHORT s1 = c[1], s4 = c[4], s7 = c[7];
        s[0] = (s1&0xf) | (c[0]<<4);
        s[1] = ((s1&0xf0)<<4)  | (c[2]);
        s[2] = (s4&0xf) | (c[3]<<4);
        s[3] = ((s4&0xf0)<<4)  | (c[5]);
        s[4] = (s7&0xf) | (c[6]<<4);
        s[5] = ((s7&0xf0)<<4)  | (c[8]);
        // print_short(s);
        s[0] -= 2048; s[1] -= 2048; s[2] -= 2048;
        s[3] -= 2048; s[4] -= 2048; s[5] -= 2048;
        return;
    }

    int check_legal(SHORT s[6]) {
        if (s[0]>=-2030 && s[0]<=2030
         && s[1]>=-2030 && s[1]<=2030
         && s[2]>=-2030 && s[2]<=2030
         && s[3]>=-2030 && s[3]<=2030
         && s[4]>=-2030 && s[4]<=2030
         && s[5]>=-2030 && s[5]<=2030) return 1; else return 0;
    }

    void pack_from_short(SHORT s[6]) {
        // These must be between -2030 and +2030.
        assert(check_legal(s));
        s[0] += 2048; s[1] += 2048; s[2] += 2048;
        s[3] += 2048; s[4] += 2048; s[5] += 2048;
        // print_short(s);
        c[0] =  (s[0]&0xff0)>>4;
        c[1] =  (s[0]&0x00f)|((s[1]&0xf00)>>4);
        c[2] =  (s[1]&0x0ff);
        c[3] =  (s[2]&0xff0)>>4;
        c[4] =  (s[2]&0x00f)|((s[3]&0xf00)>>4);
        c[5] =  (s[3]&0x0ff);
        c[6] =  (s[4]&0xff0)>>4;
        c[7] =  (s[4]&0x00f)|((s[5]&0xf00)>>4);
        c[8] =  (s[5]&0x0ff);
        // print_char();
        return;
    }

    void pack_id(uint64_t id) {
        // We also want to pack these carefully, so we don't have endian issues
        // On intel CPUs, the least significant byte is first.  We'll store so
        // that the low 4-bytes of id are a intel uint32_t in bytes 10-13.
        // We put the top byte of the 5-byte id in byte 9.
        c[9]  = (id>>32)&0xff;
        // c[10] =      id &0xff;
        // c[11] = (id>>8) &0xff;
        // c[12] = (id>>16)&0xff;
        // c[13] = (id>>24)&0xff;
        // The above is explicit, but for intel, we can do:
        *((uint32_t *)(c+10)) = id&0xffffffff;
        return;
    }
    uint64_t unpack_id() {
        uint64_t id;
        // Explicit 
        // id = (((((uint64_t)c[9]*256+(uint64_t)c[13])*256+
        //         (uint64_t)c[12])*256+(uint64_t)c[11])*256+(uint64_t)c[10]);
        // But for intel, we can do
        id = ((uint64_t)c[9]<<32) | *((uint32_t *)(c+10));
        return id;
    }

  public:
    // ----------------------------------------------------
    // Routines to unpack and pack cells
    // We subtract 2000 from the cell ijk and cpd, so that 0 is 48 in the file.
    // The check_legal limits us to the range 18..4078
    // This means that we can support slightly negative ijk (to -30).
    // cpd will be limited to 4000, so that we can go slightly above in ijk too (cpd+30).

    cell_header unpack_cell() {
        SHORT s[6]; expand_to_short(s);
        cell_header cell;
        cell.cpd = s[1]+2000;
        cell.vscale = s[2]+2000;
        cell.i = s[3]+2000;
        cell.j = s[4]+2000;
        cell.k = s[5]+2000;
        assert(cell.islegal());   // Negative values are illegal
        return cell;
    }

    cell_header pack_cell(cell_header cell) {
        SHORT s[6]; 
        assert(cell.cpd<=4000);
        c[9] = c[10] = c[11] = c[12] = c[13] = 0;
        s[0] = 0;                // Can't put the make_cell token in, or we'll fail check_legal()
        s[1] = cell.cpd-2000;
        s[2] = cell.vscale-2000;
        s[3] = cell.i-2000;
        s[4] = cell.j-2000;
        s[5] = cell.k-2000;
        pack_from_short(s);
        make_cell();
        return cell;
    }

    cell_header pack_cell(SHORT ijk[3], SHORT cpd, SHORT vscale) {
        cell_header c = cell_header(ijk,cpd,vscale);
        pack_cell(c);
        return c;
    }
    cell_header pack_cell(int ijk[3], int cpd, int vscale) {
        cell_header c = cell_header(ijk,cpd,vscale);
        pack_cell(c);
        return c;
    }
    cell_header pack_cell(int i, int j, int k, int cpd, int vscale) {
        int ijk[3]; ijk[0] = i; ijk[1] = j; ijk[2] = k;
        return pack_cell(ijk,cpd,vscale);
    }

#ifdef __THREEVECTOR_CC__
     cell_header pack_cell(integer3 &ijk, int cpd, int vscale) {
        cell_header c = cell_header(ijk,cpd,vscale);
        pack_cell(c);
        return c;
     }
#endif

    // ----------------------------------------------------
    // Routines to unpack and pack particles
    // Note that single precision is barely enough.  With CPD>2048, this will
    // lose 1-2 bits of precision in the positions.

    void unpack(double pos[3], double vel[3], uint64_t *id, cell_header cell) {
        // These return position and velocities in units of a unit box.
        // Cell centers are half-integer, so position 0 is the left-edge of cell 0
        // User should multiply by Boxsize if desired.
        // We assume that these are in the range -0.5..+0.5
        SHORT s[6]; expand_to_short(s);
        double invcpd = 1.0/cell.cpd;
        pos[0] = (double)(s[0]/2000.0+cell.i+0.5)*invcpd-0.5;
        pos[1] = (double)(s[1]/2000.0+cell.j+0.5)*invcpd-0.5;
        pos[2] = (double)(s[2]/2000.0+cell.k+0.5)*invcpd-0.5;
        vel[0] = (double)s[3]/2000.0*cell.vscale*invcpd;
        vel[1] = (double)s[4]/2000.0*cell.vscale*invcpd;
        vel[2] = (double)s[5]/2000.0*cell.vscale*invcpd;
        *id = unpack_id();
        return;
    }

    void unpack(float pos[3], float vel[3], uint64_t *id, cell_header cell) {
        // These return position and velocities in units of a unit box.
        // Cell centers are half-integer, so position 0 is the left-edge of cell 0
        // User should multiply by Boxsize if desired.
        // We assume that these are in the range -0.5..+0.5
        SHORT s[6]; expand_to_short(s);
        double invcpd = 1.0/cell.cpd;
        pos[0] = (double)(s[0]/2000.0+cell.i+0.5)*invcpd-0.5;
        pos[1] = (double)(s[1]/2000.0+cell.j+0.5)*invcpd-0.5;
        pos[2] = (double)(s[2]/2000.0+cell.k+0.5)*invcpd-0.5;
        vel[0] = (double)s[3]/2000.0*cell.vscale*invcpd;
        vel[1] = (double)s[4]/2000.0*cell.vscale*invcpd;
        vel[2] = (double)s[5]/2000.0*cell.vscale*invcpd;
        *id = unpack_id();
        return;
    }


    // This is for cell-centered positions, unit length per box.
    // Don't forget to scale the velocities to the same units!
    void pack(float pos[3], float vel[3], uint64_t id, cell_header cell) {
        SHORT s[6];
        pack_id(id);
        float norm = 2000.0*cell.cpd;
        s[0] = lround(pos[0]*norm);
        s[1] = lround(pos[1]*norm);
        s[2] = lround(pos[2]*norm);
        norm /= cell.vscale;
        s[3] = lround(vel[0]*norm);
        s[4] = lround(vel[1]*norm);
        s[5] = lround(vel[2]*norm);
        pack_from_short(s);
        return;
    }
    void pack(double pos[3], double vel[3], uint64_t id, cell_header cell) {
        SHORT s[6];
        pack_id(id);
        double norm = 2000.0*cell.cpd;
        s[0] = lround(pos[0]*norm);
        s[1] = lround(pos[1]*norm);
        s[2] = lround(pos[2]*norm);
        norm /= cell.vscale;
        s[3] = lround(vel[0]*norm);
        s[4] = lround(vel[1]*norm);
        s[5] = lround(vel[2]*norm);
        pack_from_short(s);
        return;
    }

    // Global positions must be in units of the unit box.
    // We assume that these are in the range -0.5..+0.5
    void pack_global(double pos[3], double vel[3], uint64_t id, cell_header cell) {
        // Get to local positions
        double p[3];
        double invcpd = 1.0/cell.cpd;
        p[0] = pos[0]+0.5 - (0.5 + cell.i)*invcpd;
        p[1] = pos[1]+0.5 - (0.5 + cell.j)*invcpd;
        p[2] = pos[2]+0.5 - (0.5 + cell.k)*invcpd;
        pack(p, vel, id, cell);
        return;
    }
    
    void pack_global(float pos[3], float vel[3], uint64_t id, cell_header cell) {
        // Get to local positions
        float p[3];
        float invcpd = 1.0/cell.cpd;
        p[0] = pos[0]+0.5 - (0.5 + cell.i)*invcpd;
        p[1] = pos[1]+0.5 - (0.5 + cell.j)*invcpd;
        p[2] = pos[2]+0.5 - (0.5 + cell.k)*invcpd;
        pack(p, vel, id, cell);
        return;
    }


#ifdef __THREEVECTOR_CC__
    void unpack(double3 &pos, double3 &vel, uint64_t &id, cell_header cell) {
        unpack((double*) &pos, (double*) &vel, &id, cell);
    }

    void unpack(float3 &pos, float3 &vel, uint64_t &id, cell_header cell) {
        unpack((float*) &pos, (float*) &vel, &id, cell);
    }

    void pack(double3 &pos, double3 &vel, uint64_t id, cell_header cell) {
        pack((double*) &pos, (double*) &vel, id, cell);
    }
    void pack(float3 &pos, float3 &vel, uint64_t id, cell_header cell) {
        pack((float*) &pos, (float*) &vel, id, cell);
    }

    void pack_global(double3 &pos, double3 &vel, uint64_t id, cell_header cell) {
        pack_global((double*) &pos, (double*) &vel, id, cell);
    }
    void pack_global(float3 &pos, float3 &vel, uint64_t id, cell_header cell) {
        pack_global((float*) &pos, (float*) &vel, id, cell);
    }
#endif
};


#ifdef TEST
// Here's a program to check the bit-packing
#include <stdlib.h>
#include <math.h>

#define WORST(a,b) do { double _w = fabs(b); if (_w>a) a = _w; } while(0)

void endian() {
    uint64_t id;
    id = 1;
    id = id*16+2;
    id = id*16+3;
    id = id*16+4;
    id = id*16+5;
    id = id*16+6;
    id = id*16+7;
    id = id*16+8;
    id = id*16+9;
    id = id*16+10;
    printf("Testing %llx\nTesting ",id);
    unsigned char *c = (unsigned char *)&id;
    for (int j=0; j<8;j++) printf("%02x ", c[j]);
    printf("\n\n");
    // We find that on Intel, a long long int has the MS byte at the end.

    uint32_t fid;
    fid = 1;
    fid = fid*16+2;
    fid = fid*16+3;
    fid = fid*16+4;
    fid = fid*16+5;
    fid = fid*16+6;
    fid = fid*16+7;
    fid = fid*16+8;
    printf("Testing %x\nTesting ",fid);
    c = (unsigned char *)&fid;
    for (int j=0; j<4;j++) printf("%02x ", c[j]);
    printf("\n\n");
    // Similarly for 32-bit ints.
}

int main() {
    int cpd = 2000, c[3], tc[3];
    double pos[3], vel[3];
    // double tpos[3], tvel[3];   // To test double precision, change these two lines
    float tpos[3], tvel[3];
    double worst_c = 0, worst_pos = 0, worst_vel = 0, worst_id = 0;
    uint64_t id, tid;
    pack14 cell, part;
    cell_header ch;
    srand48(1);
    endian();

    for (int n=0;n<10;n++) {
        // Make a particle
        c[0] = int(drand48()*cpd);
        c[1] = int(drand48()*cpd);
        c[2] = int(drand48()*cpd);
        // This uses cell-centered positions
        pos[0] = (drand48()-0.5)/cpd;
        pos[1] = (drand48()-0.5)/cpd;
        pos[2] = (drand48()-0.5)/cpd;
        vel[0] = (drand48()-0.5)*10/cpd;
        vel[1] = (drand48()-0.5)*10/cpd;
        vel[2] = (drand48()-0.5)*10/cpd;
        id = drand48()*1e11;

        // Pack it
        ch = cell.pack_cell(c,cpd,32);
        part.pack(pos,vel,id,ch);
        // Unpack it
        assert(cell.iscell());
        assert(part.isparticle());
        ch = cell.unpack_cell();
        assert(ch.islegal());
        part.unpack(tpos,tvel,&tid,ch);
        ch.cellid(tc);

        // Test it.  Converting to box-centered positions.
        pos[0] += (c[0]+0.5)/cpd-0.5;
        pos[1] += (c[1]+0.5)/cpd-0.5;
        pos[2] += (c[2]+0.5)/cpd-0.5;
        printf("Original %d %d %d  %f %f %f  %f %f %f %lld\n",
                c[0], c[1], c[2], 
                pos[0], pos[1], pos[2],
                vel[0], vel[1], vel[2], id);
        printf("Restored %d %d %d  %f %f %f  %f %f %f %lld\n",
                tc[0], tc[1], tc[2], tpos[0], tpos[1], tpos[2], tvel[0], tvel[1], tvel[2], id);
        printf("Residual %d %d %d  %f %f %f  %f %f %f %lld\n",
                tc[0]-c[0], tc[1]-c[1], tc[2]-c[2], 
                tpos[0]-pos[0], tpos[1]-pos[1], tpos[2]-pos[2], 
                tvel[0]-vel[0], tvel[1]-vel[1], tvel[2]-vel[2], id-tid);
        printf("\n");
        WORST(worst_c,tc[0]-c[0]);
        WORST(worst_c,tc[1]-c[1]);
        WORST(worst_c,tc[0]-c[0]);
        WORST(worst_pos,tpos[0]-pos[0]);
        WORST(worst_pos,tpos[1]-pos[1]);
        WORST(worst_pos,tpos[2]-pos[2]);
        WORST(worst_vel,tvel[0]-vel[0]);
        WORST(worst_vel,tvel[1]-vel[1]);
        WORST(worst_vel,tvel[2]-vel[2]);
        WORST(worst_id,tid-id);
    }
    printf("Largest position residual is %e (box), %f (cell)\n", worst_pos, worst_pos*cpd);
    printf("Largest velocity residual is %e (box), %f (cell)\n", worst_vel, worst_vel*cpd);
    printf("Largest id residual is %d\n", worst_id);
    return 0; 
}

#endif //  TEST

#endif //  PACKED_CPP


#ifdef DONOTCOMPILE
// Sample code to read a file:

pack14 particle;
cell_header cellhead;
double pos[3], vel[3];
uint64 id;

while (fread(&particle, sizeof(pack14), 1, fp)!=EOF) {
    if (particle.iscell()) {
        cellhead = particle.unpack_cell();
    } else {
        assert(cellhead.islegal());  
        particle.unpack(pos,vel,cellhead);
        // And now do something with pos & vel, 
        // such as multiply by Boxsize/CPD and store or output
    }
}


// Sample code to write one abacus cell to a file:

// Fetch the data and setup
Cell c = PP->GetCell(i,j,k);
int vscale = 32;  // Choose this based on the max velocity
float vel_to_zspace = 1.0;   // Apply the cosmological factor to get to zspace units

// Write the cell:
cell_header cell;
pack14 particle;
cell = particle.pack_cell(i,j,k,PP->cpd,vscale);
fwrite(&particle, sizeof(pack14), 1, fp);

// And now the particles:
for (int p=0; p<c.count(); p++) {
    FLOAT3 pos, vel;
    pos = c.pos[p];
    vel = c.vel[p];
    id = c.aux[p].pid();
    pos *= P.cpd;                // Must scale to cell being unit length 
    vel *= P.cpd*vel_to_zspace;
    particle.pack(pos,vel,id,cell);    // Or pack_global
    fwrite(&particle, sizeof(pack14), 1, fp);
}

#endif 
