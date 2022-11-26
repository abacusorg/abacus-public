#ifndef __IC_BASE_H
#define __IC_BASE_H

/* This class contains the ICFile abstract base class; it defines the interface
 * that each IC format class must implement.  New IC formats should define a
 * class in loadIC.cpp that inherits ICFile and register itself in the 
 * ICFile::FromFormat() factory function at the bottom of that file.
 *
 * New IC classes must override the following function and define a constructor:
 * - unpack(): load the particles into the insert list
 * - Constructor: must record the number of particles in Npart
 *
 * New IC classes may override the following two functions:
 * - read_nonblocking(): start a nonblocking read of IC files
 * - check_read_done(): report whether the read is done
 *
 * The unpack functions will be provided the pos and vel unit conversions
 * that the user specified to get the particles into code units.
 */

#include <memory>  // for unique_ptr

class ICFile {         // This is an abstract base class.
public:
    // Used to generate sequential PIDs if we were not given them
    static uint64 next_pid; // 0 by default
    uint64 Npart;  // initialized by each subclass constructor
    uint64 foffset;  // where to start reading the file
    uint64 fbytes;  // size of the data to read

    int slab;
    int zsplit;

    FLOAT lpt_vel_scale;

    uint64 NsubsampleA, NsubsampleB;  // subsample A and B counts

    // This factory function is the primary entry point to construct an ICFile object
    static unique_ptr<ICFile> FromFormat(const char *format, int _slab);

    // Unpack the particles in the arena and push them to the Insert List
    uint64 unpack_to_IL(double convert_pos, double convert_vel);

    // Start a read via SB
    virtual void read_nonblocking();

    // If a read is in progress, return 0
    virtual int check_read_done();

    ICFile(int _slab, int _zsplit);

    virtual ~ICFile(void) = default;

    inline static void set_taggable_bits(auxstruct &aux, uint64 &sumA, uint64 &sumB);

private:
    virtual uint64 unpack(double convert_pos, double convert_vel) = 0;
};

void get_IC_unit_conversions(double &convert_pos, double &convert_vel);

#endif
