/// Euler16.c

/* This code takes two orthogonal axes (notionally the major and minor)
and compresses them (very lossily) into a single uint16_t.

The result has about 4 degree pixels, so typical rms resolution is 
below 2 degrees.  That is more than enough for correlation studies.
The major axis has 1.6 deg of 2-d scatter.

What is happening is that the major (first) axis is being stored
in a variant of spherical coords, and then the second is having its
rotation around that axis stored.  In this sense, it is like storing
the Euler angles of the rotation over this coordinate system.

Note that this is treating these as unsigned axes, i.e., we
assume that parity-inverting an input doesn't matter.  That
saves us a factor of 4 in bins.

The good thing about this pixelization is that it uses very little
trigonometry (but 4 square roots and ~10 conditionals) and it
keeps the pixels relatively square.  A bad thing is that the area
of the pixels varies by about 1.6 peak-to-peak (and 15% rms).

*/

#define EULER_TEST
#ifdef EULER_TEST
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#endif

/// Given 2 vectors that must be orthogonal, return an 16-bit binning.
// The vectors need not be unit normalized.
// The orthogonality needs to be close, but not perfect.
uint16_t pack_euler16(float major[3], float minor[3]) {
    // First, we need to deal with the major axis.
    // We begin by deciding which coordinate is biggest
    int cap;
    // We have to use unit vectors
    float norm = 1.0/sqrt(major[0]*major[0]+major[1]*major[1]+major[2]*major[2]);
    major[0] *= norm; major[1] *= norm; major[2] *= norm;
    float xx, yy, zz;   // The cap-based coordinate system
    if (fabs(major[0])>fabs(major[1]) && fabs(major[0])>fabs(major[2])) {
        // x axis is biggest
        cap = 0; 
        if (major[0]<0) {    // Use the positive cap
            major[0] = -major[0]; major[1] = -major[1]; major[2] = -major[2]; 
        }
        zz = major[0];
        if (fabs(major[1])>fabs(major[2])) {    // y is bigger than z
            if (major[1]>0) { yy = major[1]; xx = major[2]; }
                       else { yy = -major[1]; xx = major[2]; cap += 1; }
        } else {
            if (major[2]>0) { yy = major[2]; xx = major[1]; cap += 2; }
                       else { yy = -major[2]; xx = major[1]; cap += 3; }
        }
    } else if (fabs(major[1])>fabs(major[2])) {
        // y axis is biggest
        cap = 4; 
        if (major[1]<0) {    // Use the positive cap
            major[0] = -major[0]; major[1] = -major[1]; major[2] = -major[2]; 
        }
        zz = major[1];
        if (fabs(major[2])>fabs(major[0])) {    // z is bigger than x
            if (major[2]>0) { yy = major[2]; xx = major[0]; }
                       else { yy = -major[2]; xx = major[0]; cap += 1; }
        } else {
            if (major[0]>0) { yy = major[0]; xx = major[2]; cap += 2; }
                       else { yy = -major[0]; xx = major[2]; cap += 3; }
        }
    } else {
        // z axis is biggest
        cap = 8; 
        if (major[2]<0) {    // Use the positive cap
            major[0] = -major[0]; major[1] = -major[1]; major[2] = -major[2]; 
        }
        zz = major[2];
        if (fabs(major[0])>fabs(major[1])) {    // x is bigger than y
            if (major[0]>0) { yy = major[0]; xx = major[1]; }
                       else { yy = -major[0]; xx = major[1]; cap += 1; }
        } else {
            if (major[1]>0) { yy = major[1]; xx = major[0]; cap += 2; }
                       else { yy = -major[1]; xx = major[0]; cap += 3; }
        }
    }


    // Now we have 12 "triangular" caps, with zz>=yy>=abs(xx).
    // The value s = yy/zz is [0,1], and r = xx/yy is [-1,1].
    // We transform to t = sqrt( 1 - 1/sqrt(1+s*s) ) / sqrt(1-1/sqrt(2))
    // t is also in the [0,1] range, but it is closer to uniform density 
    // for |r|<<1.
    #define EULER_NORM 1.8477590650225735122    // 1/sqrt(1-1/sqrt(2))
    #define EULER_TBIN 11
    #define EULER_ABIN 45

    float r = xx/(yy+1e-10);
    float t = yy/zz;
    t = EULER_NORM * sqrt(1-1/sqrt(1+t*t));

    // Now we grid t into 11 equal bins, numbered [0,11)
    // In each, we divide r into 2*i+1 bins.
    // This creates a fairly uniform division of the triangular area,
    // while avoiding elongated regions.
    // In total, there are 121 bins per cap.
    // Since we had about 20,600 sq deg of hemisphere, 
    // this is about 14.2 sq deg per location of the major axis.
    // That's equiv to a circle of 4.2 degree diameter.

    int it = floor(t*EULER_TBIN); if (it>=EULER_TBIN) it = EULER_TBIN-1;
    int ir = floor(0.5*(r+1.0)*(2*it+1));
    if (ir>=2*it+1) ir = 2*it;
    
    int bin = (EULER_TBIN*EULER_TBIN)*cap + it*it + ir;
    // There are 121*12 = 1452 such bins.  

    // printf("cap = %d, r = %f, t = %f, ir = %d, it = %d, bin = %d\n", cap, r, t, ir, it, bin);

    // Next, we need to track the rotation of the minor axis.
    // This will be divided into 45 even bins.
    // Since the rotation is only defined to 180 degrees, this is 4 deg.
    // There are 12*121*45 = 65340 possible bins, nearly all of 16 bits.

    norm = 1.0/sqrt(minor[0]*minor[0]+minor[1]*minor[1]+minor[2]*minor[2]);
    minor[0] *= norm; minor[1] *= norm; minor[2] *= norm;

    // We know that the minor axis is orthogonal to the major axis.
    // And we know where the major axis was pointing, from the cap.
    // We are just going to project the minor axis into the plane normal
    // to the coord axis closest to the major axis.  Then measure the 
    // angle in that plane.
    if (cap>=8) {
        // Major axis was closest to the z axis.
        xx = minor[0]; yy = minor[1];
    } else if (cap>=4) {
        // Major axis was closest to the y axis.
        xx = minor[2]; yy = minor[0];
    } else {
        // Major axis was closest to the x axis.
        xx = minor[1]; yy = minor[2];
    }
    if (yy<0) { xx = -xx; yy = -yy; }
    float az = atan2(yy,xx);   // This is [0,PI]
    int iaz = floor(EULER_ABIN*az/M_PI);
    if (iaz>=EULER_ABIN) iaz = 0;
    bin = EULER_ABIN*bin + iaz;
    // printf("iaz = %d, bin = %d\n", iaz, bin);
    return (uint16_t)bin;
}



/// Given a 16-bit bin, return two unit vectors.
/// The supplied pointers must point to float[3] memory each.
void unpack_euler16(uint16_t bin, float *major, float *minor) {
    int cap = bin/EULER_ABIN;
    int iaz = bin - cap*EULER_ABIN;   // This is the minor axis bin
    bin = cap;
    cap = bin/(EULER_TBIN*EULER_TBIN);   // This is the cap
    bin = bin - cap*(EULER_TBIN*EULER_TBIN);
    int it = floor(sqrt(bin));   
    int ir = bin - it*it;

    float t = (it+0.5)*(1.0/EULER_TBIN);   // [0,1]
    float r = (ir+0.5)/(it+0.5)-1.0;            // [-1,1]

    // We need to undo the transformation of t to get back to yy/zz
    t *= (1.0/EULER_NORM);
    t = t * sqrt(2.0-t*t)/(1.0-t*t);   // Now we have yy/zz

    // printf("cap = %d, ir = %d, it = %d, r = %f, t = %f\n", cap, ir, it, r, t);

    float yy = t;
    float xx = r*t;
    // and zz=1
    float norm = 1.0/sqrt(1.0+xx*xx+yy*yy);
    float zz = norm;
    yy *= norm; xx *= norm;  // These are now a unit vector

    switch(cap) {
       case  0: major[0]= zz; major[1]= yy; major[2]= xx; break;
       case  1: major[0]= zz; major[1]=-yy; major[2]= xx; break;
       case  2: major[0]= zz; major[1]= xx; major[2]= yy; break;
       case  3: major[0]= zz; major[1]= xx; major[2]=-yy; break;

       case  4: major[1]= zz; major[2]= yy; major[0]= xx; break;
       case  5: major[1]= zz; major[2]=-yy; major[0]= xx; break;
       case  6: major[1]= zz; major[2]= xx; major[0]= yy; break;
       case  7: major[1]= zz; major[2]= xx; major[0]=-yy; break;

       case  8: major[2]= zz; major[0]= yy; major[1]= xx; break;
       case  9: major[2]= zz; major[0]=-yy; major[1]= xx; break;
       case 10: major[2]= zz; major[0]= xx; major[1]= yy; break;
       case 11: major[2]= zz; major[0]= xx; major[1]=-yy; break;
    }

    // Next, we can get the minor axis
    float az = (iaz+0.5)*(1.0/EULER_ABIN)*M_PI;
    xx = cos(az);
    yy = sin(az);
    // printf("az = %f, %f, %f\n", az, xx, yy);
    // We have to derive the 3rd coord, using the fact that the two axes
    // are perpendicular. 
    switch (cap/4) {
        case 2: minor[0] = xx; minor[1] = yy; 
            minor[2] = (minor[0]*major[0]+minor[1]*major[1])/(-major[2]);
            break;
        case 0: minor[1] = xx; minor[2] = yy; 
            minor[0] = (minor[1]*major[1]+minor[2]*major[2])/(-major[0]);
            break;
        case 1: minor[2] = xx; minor[0] = yy; 
            minor[1] = (minor[2]*major[2]+minor[0]*major[0])/(-major[1]);
            break;
    }
    norm = 1.0/sqrt(minor[0]*minor[0]+minor[1]*minor[1]+minor[2]*minor[2]);
    minor[0] *= norm;
    minor[1] *= norm;
    minor[2] *= norm;
    return;
}

#ifdef EULER_TEST

float test(float ax, float ay, float az, float bx, float by, float bz) {
    float norm, ret;
    // printf("Dot = %f\n", ax*bx+ay*by+az*bz);
    norm = 1.0/sqrt(ax*ax+ay*ay+az*az);
    ax *= norm; ay *= norm; az *= norm;
    norm = 1.0/sqrt(bx*bx+by*by+bz*bz);
    bx *= norm; by *= norm; bz *= norm;
    float major[3], minor[3];
    major[0] = ax;
    major[1] = ay;
    major[2] = az;
    minor[0] = bx;
    minor[1] = by;
    minor[2] = bz;

    uint16_t bin = pack_euler16(major, minor);
    major[0] = 0.0; 
    major[1] = 0.0; 
    major[2] = 0.0; 
    minor[0] = 0.0; 
    minor[1] = 0.0; 
    minor[2] = 0.0; 
    unpack_euler16(bin, major, minor);
    // printf("%d\n", bin);

    float dx, dy, dz;
    dx = major[0]-ax;
    dy = major[1]-ay;
    dz = major[2]-az;
    norm = sqrt(dx*dx+dy*dy+dz*dz);
    if (norm>1) { 
        dx = major[0]+ax;
        dy = major[1]+ay;
        dz = major[2]+az;
        norm = sqrt(dx*dx+dy*dy+dz*dz);
    }
    // printf("%f %f %f -> %f %f %f = %f\n", ax,ay,az, major[0], major[1], major[2], norm);
    ret = norm;

    dx = minor[0]-bx;
    dy = minor[1]-by;
    dz = minor[2]-bz;
    norm = sqrt(dx*dx+dy*dy+dz*dz);
    if (norm>1) {
        dx = minor[0]+bx;
        dy = minor[1]+by;
        dz = minor[2]+bz;
        norm = sqrt(dx*dx+dy*dy+dz*dz);
    }
    // printf("%f %f %f -> %f %f %f = %f\n\n", bx,by,bz, minor[0], minor[1], minor[2], norm);
    // ret = norm;
    return ret;
}


float random_test() {
    float ax, ay, az, norm;
    float bx, by, bz;
    do {
        ax = drand48()*2-1;
        ay = drand48()*2-1;
        az = drand48()*2-1;
        norm = sqrt(ax*ax+ay*ay+az*az);
    } while (norm>1 || norm<0.1);
    ax /= norm;
    ay /= norm;
    az /= norm;
    // The major axis
    do {
        bx = drand48()*2-1;
        by = drand48()*2-1;
        bz = drand48()*2-1;
        norm = sqrt(bx*bx+by*by+bz*bz);
    } while (norm>1 || norm<0.1);
    norm = ax*bx+ay*by+az*bz;
    bx -= norm*ax;
    by -= norm*ay;
    bz -= norm*az;
    norm = sqrt(bx*bx+by*by+bz*bz);
    bx /= norm;
    by /= norm;
    bz /= norm;
    // The minor axis
    return test(ax,ay,az, bx,by,bz);
}
    

int main() {
    srand48(123);
    // test(1.0, 0.9, 0.8, -1.7, 1.0, 1.0);
    // test(0.9, -1.0, 0.8, 1.0, 0.1, -1.0);
    float sum = 0.0, max = 0.0, val;
    int iter = 50;
    // iter = 5e6;
    for (int j=0; j<iter; j++) { 
        val = random_test(); 
        sum += val; if (val>max) max=val; 
    }
    printf("RMS in degrees = %f, max %f\n", sum/iter*180/M_PI, max*180.0/M_PI);
    return 0;
}

#endif
