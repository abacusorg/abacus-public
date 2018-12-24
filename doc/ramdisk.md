# Abacus Ramdisk Interface
Author: Lehman Garrison

Abacus is ramdisk-aware, meaning some of its memory allocation and file IO behaviors
are different when files are on ramdisk instead of a normal file system.  This document
explains the design tradeoffs that went into determining these behaviors and interfaces.

## Background

Abacus is designed to operate on problems that do not fit in memory, usually by storing
particles in slab files on disk while they are not needed.  Sometimes the problem does
fit in memory, but it is convenient to maintain the illusion of files being written to
disk.

We achieve this with a ramdisk.  A ramdisk is a file system that exists in RAM.  It is accessible as a normal directory;
on Linux, it is automatically available under "/dev/shm/" <sup name="a1">[1](#f1)</sup>
("shm" stands for "shared memory").  By default, half of a system's RAM is available as
a ramdisk.

The ramdisk model offers many advantages over a pure in-memory version of Abacus (in which
the executable would be invoked once and run all timesteps before exiting):

1) The singlestep executable can be responsible for exactly one time step, meaning global
objects 

2) 


<b name="f1">1</b> This path seems quite standard; see https://www.kernel.org/doc/Documentation/filesystems/tmpfs.txt. [↩](#a1)
