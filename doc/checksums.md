# File Checksums in Abacus

We checksum various output files in Abacus so that we can verify their integrity
as they are moved around for post-processing and analysis.  The checksums are
stored in a file called `checksums.crc32` in each directory that contains checksummed
files.  As the file name implies, these are CRC32 checksums (32-bit) and can be
verified with the ubiquitous GNU `cksum` util.  Or, for much faster alternative,
we distribute a util called `fast_cksum` with Abacus.

We are indebted to Stephan Brumme's Fast CRC32 code as the core of our CRC32
functionality: https://create.stephan-brumme.com/crc32/


## Motivation
While checksumming files is rarely a bad idea, we had not implemented it as an
Abacus feature until we were confronted with a rather dramatic example of file
corruption.  While analyzing one our first 7000^3 simulations, we observed 24 files
(out of 1875) in a 4.5 TB slice corrupted during a simple `cp` across filesystems!
A binary diff showed the errors occurring in 1 MB blocks in the middle of the files
(the file sizes were the same).

Why implement this on-the-fly and not as a post-processing step?  The first concern
is that reading data is expensive!  It's already passing through the IO thread
on its way to disk, which seems like the perfect time to checksum it.  And secondly,
once the data is on disk, a subsequent read might exhibit corruption.  So checksumming
in memory before it ever hits disk is preferred.

While MD5 is more common for user-facing file integrity checks, CRC32 is also quite
common and significantly cheaper computationally.  The GNU coreutils `md5sum` and
`cksum` utilities were actually the same speed in our tests, but a quick search
revealed that much more efficient implementations of CRC are possible.


## Abacus integration
The `include/crc32_fast.cpp` and `include/crc32_fast.h` files contain our CRC32
implementation.  The IO thread in `singlestep/Output/io_thread.cpp` runs the
checksums if the `iorequest` asked for it.  This is determined by `SlabBuffer::WantChecksum()`,
which has an enumeration of the `SlabType`s that opt in to checksumming.  In general,
we checksum outputs but not state.

There is a parameter file option called `NoChecksum` to globally disable checksumming.

The IO module support checksumming even for incremental writes (i.e. appending to a file).
This is useful for the lightcones.  We keep a running checksum and then "finalize" the
checksum (see "GNU cksum") before writing it out.

## How to verify checksums

Change directories so the filenames in the checksums file won't have a path pre-pended:
```bash
$ cd /path/to/files
```

Run the checksum (on \*.dat, for example):
```bash
$$ABACUS/util/fast_cksum *.dat > cksum.out
```
Or with the slower GNU utility:
```bash
$ cksum /path/to/files/\*.dat > cksum.out
```

Then compare the checksum files:
```bash
$ diff cksum.out checksums.crc32
```

The files look like:
```
1758675646 144185735 DESI_L760_N2660_prototype.z0.100.slab0000.field_pack9.dat
1906361843 123906288 DESI_L760_N2660_prototype.z0.100.slab0000.field_pack9_pids.dat
3087645703 100255241 DESI_L760_N2660_prototype.z0.100.slab0000.L0_pack9.dat
2775744001 88244488 DESI_L760_N2660_prototype.z0.100.slab0000.L0_pack9_pids.dat
```

The format is checksum, file size, and filename.

## Performance
High performance in CRC32 is obtained by precomputing a lookup table with 256 values,
one for each possible outcome of 1 byte convolved with the CRC32 generating polynomial.
And even higher performance is obtained by precomputing 16 such tables from the
composition of lookups, allowing for processing 16 bytes at a time.  This gets our
speed up to several GB/s on one CPU core, which is fast enough that it doesn't limit
our IO.

We also include an extension of the table scheme to 32 bytes, but it has virtually no
impact on performance (while going from 8 to 16 is about a factor of 2).

Checksumming does mean that the IO threads have more work to do than they did before,
so we'll have to keep a careful eye on their performance.  We may not be able to get
away with stacking all the IO threads on one core, as we sometimes did in the past.

## GNU cksum
The GNU coreutils' cksum uses an unfortuantely non-standard variant of CRC32, but cksum's
ubiquity on Linux system means we probably want to conform to that variant.  The generating
polynomial is the usual one, just bit-reversed (see: https://en.wikipedia.org/wiki/Cyclic_redundancy_check#Polynomial_representations_of_cyclic_redundancy_checks).  The byte processing order is also swapped.  We regenerated the lookup
tables in Stephan Brumme's CRC32 code in this convention, and modified the algorithms
similarly.

This page describes GNU's cksum variant: https://pubs.opengroup.org/onlinepubs/009695399/utilities/cksum.html

## Checksums in parallel runs
In multi-node parallel runs, each node computes checksums on its own output files.  Each node
writes its own checksum file, called `checksums.NNNN.crc32`.  The `merge_checksum_files()`
function in `abacus.py` concatenates these files (ordering the lines by filename).

We don't formally signal which directories might contain checksum files, so `merge_checksum_files()`
just looks in the most likely directories.  It's not robust, but it's also unlikely to change
very often.

We could imagine doing an MPI reduction of the checksums, but that would be rather verbose
and inelegant if implemented in C++.
