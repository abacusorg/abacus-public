# Abacus Ramdisk Interface
Author: Lehman Garrison

Abacus is ramdisk-aware, meaning some of its memory allocation and file IO behaviors
are different when files are on ramdisk instead of a normal file system.  This document
explains these behaviors and the design considerations that went into determining them.

## Background

Abacus is designed to operate on problems that do not fit in memory, usually by storing
particles in slab files on disk while they are not needed.  Sometimes the problem does
fit in memory, but it is convenient to maintain the illusion of files being written to
disk.

We achieve this with a ramdisk.  A ramdisk is a file system that exists in RAM.  It is accessible as a normal directory;
on Linux, it is automatically available under "/dev/shm/" <sup name="a1">[1](#f1)</sup>
("shm" stands for "shared memory").  By default, half of a system's RAM is available as
a ramdisk.

Files written to a ramdisk will stay there even after the writing process exits.  They can
then be read by a new process.  The ramdisk memory is managed by the kernel.

The ramdisk model offers many advantages over a pure in-memory version of Abacus (in which
the executable would be invoked once and run all timesteps before exiting):

1) The singlestep executable can be responsible for exactly one time step, meaning global
variables and objects can assume cosmology values and tuning parameters relevant to the current epoch.

2) Writing/restoring a backup is as simple as copying files to/from the ramdisk directory.

3) The code can remain almost identical to the non-ramdisk version, so we can continue
to support normal disk and ramdisk with almost no top-level logic changes.

## Overview

### Implementation

The basic steps of writing to a ramdisk are:

1) Open a file descriptor to the file on ramdisk with `open()` or `shm_open()`
2) Set the desired size of the file with `ftruncate()`
3) Map that file to the process's address space with `mmap()`.

This will return a pointer to a region of memory, just like `malloc()` would.  One
can write to this memory and immediately see the changes reflected in the `/dev/shm` file.

Don't forget to un-map the memory with `munmap()` when you're done!  Shared memory is only
freed back to the system when all mappings and the `/dev/shm` file handle are deleted.  Of
course, all mappings will be removed when a process exits.

### Ramdisk Overwrites


### Corner Cases
Multipole recovery
Group finding outputs




<b name="f1">1:</b> This path seems quite standard; see https://www.kernel.org/doc/Documentation/filesystems/tmpfs.txt. [↩](#a1)

## Resources
- https://github.com/abacusorg/abacus/issues/215
- https://github.com/abacusorg/abacus/pull/218
- http://man7.org/training/download/posix_shm_slides.pdf
- http://man7.org/linux/man-pages/man7/shm_overview.7.html
- http://man7.org/linux/man-pages/man3/shm_open.3p.html
