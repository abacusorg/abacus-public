# Abacus Ramdisk Interface
Author: Lehman Garrison

[Pull Request #218](https://github.com/abacusorg/abacus/pull/218)

Abacus is ramdisk-aware, meaning some of its memory allocation and IO behaviors
are different when files are on ramdisk instead of a normal file system.  This document
explains these behaviors and the overall design of the Abacus ramdisk interface.


## Background

Abacus is designed to operate on problems that do not fit in memory, usually by storing
particles in slab files on disk while they are not needed.  Sometimes the problem *does*
fit in memory, but it is convenient to maintain the illusion of files being written to disk.

We achieve this with a ramdisk.  A ramdisk is a file system that exists in RAM.  It is a
Linux kernel service that is accessible as a normal directory;
by default, it is available under `/dev/shm/` <sup name="a1">[1](#f1)</sup>
("shm" stands for "shared memory").  By default, half of a system's RAM is available as
a ramdisk.

Files written to a ramdisk will stay there even after the writing process exits.  They can
then be read by a new process.  The ramdisk files stay in memory as long as the Linux kernel
is active.

The ramdisk model offers many advantages over a pure in-memory version of Abacus (in which
the executable would be invoked once and run all timesteps before exiting):

1) The singlestep executable can be responsible for exactly one time step, meaning global
variables and objects can assume cosmology values and tuning parameters relevant to the current epoch.

2) Writing/restoring a backup is as simple as copying files to/from the ramdisk directory.

3) The code can remain almost identical to the non-ramdisk version, so we can continue
to support normal disk and ramdisk with almost no top-level logic changes.


## Overview

In singlestep, arena allocations are automatically performed in shared memory if the
file path for that arena is on ramdisk (i.e. starts with `/dev/shm/`). "Reading" an
arena consists of attaching a pointer, so the "IO" is effectively instantaneous.
There are no `memcpy()`s involved, unlike old versions of Abacus.

In convolution, the multipoles are automatically mapped from shared memory if they
exist on the ramdisk.  One should use `Conv_IOMode = "overwrite"` for maximum efficiency;
that way, the Taylors can write directly into the multipoles memory.  The convolution
driver renames the multipoles files to Taylors at the end of a successful convolution.

**Important**: the `RamDisk` flag in the parameter file is now deprecated in favor of
direct detection of whether a file path is on the ramdisk.

**Caution**: one should probably use `StateIOMode = "overwrite"` on the ramdisk.
This is because the original slabs will get modified regardless, so one should
probably embrace it.


## singlestep

Some arenas we want to allocate in shared memory and some we don't.  When a new arena
allocation is requested, `SlabBuffer` checks whether the arena type is
a known "ramdisk type", like positions or multipoles.  In-memory slabs like `AccSlab`
are not considered ramdisk types.

A ramdisk slab can either be a read-type or a write-type (`RAMDISK_READ` or `RAMDISK_WRITE`).
A read type means that we should attach to existing shared memory; a write-type means
we should make a new allocation.  Again, these are explicitly itemized in `SlabBuffer::IsRamdiskSlab()`.

When `SlabBuffer` asks `ArenaAllocator` to make a new allocation, it may pass the ramdisk
read/write type flag to request a shared-memory allocation.  If it does so, it must also pass
the path at which to read/write the slab.  This path is usually just the canonical path from 
`Read`/`WriteSlabPath()`.

Some slabs may not be allocated in shared memory but end up with an output request to ramdisk,
like `TimeSlice` or group finding outputs.  Unfortunately, direct IO will break if one tries
to write to ramdisk!  Previously, one had to turn on the `RamDisk` parameter to tell the IO
module to never used direct IO and instead use the slower `fopen()`/`fwrite()` interface.
Now, when an IO request is created, the path is checked and flagged as ramdisk/non-ramdisk.
The IO thread then uses the appropriate IO mode.

This automatic detection and switching is only implemented for the IO thread module; it's
likely we want to deprecate the other IO modules (`io_fopen` and `io_dio`).  The IO thread
module has a blocking mode, so we're not losing flexibility that way.

The `ArenaAllocator` keeps a flag indicating whether an arena exists in shared memory.
Those arenas differ from normal arenas in a few ways: they do not have guard bytes, they
are never "oversized", and the allocations are never recycled.  The guard bytes would be written
to disk and complicate file size checks and direct examination of state files.  The oversizing/
recycling isn't applicable in the ramdisk model.

`SlabBuffer::IsRamdiskSlab(int type, int hint)` accepts a hint.  This is for cases where
we want to override normal defaults about a slab type being read or write, usually because
we are calling from a context like `ReadArena(...)` that gives us high confidence that a slab
is indeed a read type.  We use this in multipole recovery mode (see below).

### Ramdisk slab life cycle

To understand the ramdisk interface, consider how we would load a `PosSlab`:

1) `FetchSlabAction(...)` starts the read with `SlabBuffer::LoadArenaNonBlocking(PosSlab,slab)`.
2) `LoadArenaNonBlocking(...)` calls `SlabBuffer::AllocateArena(...)`, which triggers `SlabBuffer::AllocateSpecificSize(...)`.
3) `AllocateSpecificSize(...)` checks `SlabBuffer::IsRamdiskSlab(...)`.  This returns `RAMDISK_READSLAB` because the `PosSlab` path starts with `/dev/shm`.
4) `AllocateSpecificSize(...)`, seeing `RAMDISK_READSLAB`, calls `SlabBuffer::ReadSlabPath(...)` to get the slab file name.
5) `AllocateSpecificSize(...)` asks for the allocation with `ArenaAllocator::Allocate(...)`, passing the `RAMDISK_READSLAB` flag and the file name.
6) `Allocate(...)` receives the ramdisk flag and calls `open(...)` and `mmap(...)` to attach the existing allocation.
7) Having completed the allocation, `LoadArenaNonBlocking(...)` tries to call `SlabBuffer::ReadArenaNonBlocking(...)`.
This triggers a check of `ArenaAllocator::ArenaRamdiskType(...)`, which returns `RAMDISK_READSLAB`.  Thus, `ReadArena` knows to skip the read and instead calls `SetIOCompleted(...)`.
8) The arena is now ready to use.

When the `MergePosSlab` is allocated, it follows the same path following from a direct call
to `AllocateSpecificSize(...)` from `FillMergeSlab(...)`.  It is detected as a `RAMDISK_WRITESLAB`,
so `Allocate(...)` deletes the original `PosSlab` file on disk if overwriting (see "State file
overwrites" below).

Afterwards, we release the original `PosSlab` memory:

1) `FinishAction(...)` calls `SlabBuffer::DeAllocate(PosSlab, slab)`, which immediately triggers `ArenaAllocator::DeAllocateArena(...)`.
2) `DeAllocate(...)` sees that this is a shared memory arena and calls `ArenaAllocator::DiscardArena(...)`
3) `DiscardArena(...)` unmaps the memory with `munmap(...)`.

### Multipole recovery

If we have positions but not multipoles (e.g. they were lost), then we run "Multipole Recovery"
mode.  This just loads the positions into the merge slabs and outputs multipoles.  When on the
ramdisk, we would like to be able to map the positions directly into the merge slabs so we don't
have to make an expensive `memcpy`.  So we add the merge slabs to `ReadSlabPath`, but with
the paths from the *read* state, not the write.

But the merge slabs would still be considered "write slabs" by `IsRamdiskSlab()`; this is where
the "hint" system comes in.  `LoadArenaBlocking()` always hints `RAMDISK_AUTO_READSLAB`, which
indicates to `IsRamdiskSlab()` that if the path is a ramdisk path, it should be treated as a
read slab.  Thus, `AA` will read from an existing allocation (in this case, the `PosSlab`).

### State file overwrites

When overwriting a slab, we're writing a new "file" on the ramdisk with the same name.  But we
might not be done with that slab; for example, the orginal `CellInfo` slabs are still used after
we write the `MergeCellInfo` slab.  Fortunately, the memory remains available even if one deletes
the file on ramdisk, as long as there's an open `mmap`-ing.  This is POSIX guaranteed <sup name="a2">[2](#f2)</sup>.
So one just has to delete/unlink the file and create a new one with the same name; this will
produce a new file descriptor and thus a new shared memory allocation.

### Performance

`munmap()` seems to be a fairly expensive operation, probably because it holds a global lock
for the duration of the syscall.  `tcmalloc` probably can't help with this, since the call
needs to operate on kernel-held memory.  We don't strictly need to `munmap()` everything though,
only the slabs that are getting deleted.  So we could skip the `munmap()` for the merge
slabs (for example) and save about half the calls.  The rest will be unmapped at process exit.
We would have to test if that actually saves time or just defers it to the end.

### Design justifications

We explicitly itemize the slabs that go on the ramdisk.  One could imagine just determining
this based on whether the slab file path begins with `/dev/shm/`, but there ended up being
too many exceptions, like `NearAccSlab`, which we write in `ForceOutputDebug` mode but
never any other time.

Similarly, we explicitly mark slab types as being either read- or write-type.  One could
imagine determining this by presence in the `ReadSlabPath` or `WriteSlabPath`, but there
are again too many exceptions.  Plus, we'd probably prefer to hard-fail if we accidentally
try to get a read path for a write slab (and vice versa).

One could also imagine determining read-/write-type based on the presence of a file on
ramdisk.  One would read an existing slab or create it if it didn't exist.  That works
fine until we start overwriting.

Even if a slab is "ramdisk-able", we need to decide whether Abacus is using a ramdisk or not.
Historically, we've used the `RamDisk` flag in the `.par` file.  That's a global setting
and is too clumsy to support scenarios where some paths are on ramdisk and some aren't
(e.g. state is on ramdisk but the time slice output is going to RAID).  So now, we look at
the file path of the slab and check if it starts with `/dev/shm/`.  This seems incredibly
hacky, but there doesn't seem to be any other way to determine whether a given path is on a ramdisk.
Also, `/dev/shm/` seems to be a highly standard choice across Linux kernels.  In Abacus, we've made
this path configurable with `./configure --with-ramdisk=PATH`.

The "ramdisk aware" features are mostly coordinated from `ArenaAllocator` and `SlabBuffer`.
One could imagine coordinating it from the IO module instead; `ArenaAllocator` would set up
an arena with a NULL pointer, and "reading" the arena would instead consist of attaching
the shared memory.  That way, `ArenaAllocator` would never need to know about file paths
which would be nice.

However, when we allocate an arena for writing (e.g. `MergePosSlab`), we expect the memory
to exist right away without ever passing through the IO module. So that needs to be done
from `ArenaAllocator`.  So rather than split the ramdisk features between the `AA` and IO
modules, we opted to consolidate them into `AA`.  This also preserves the property that
the arena memory is immediately present after an allocation.

We chose to explicitly pass the ramdisk path to the `ArenaAllocator`.  Another option
would be to have `AA` query the `SlabBuffer` for the path.  That's a circular dependency
though, so `SlabBuffer` would have to be declared in a separate header file.  Do-able,
but it breaks with Abacus convention.

Furthermore, `AA` only knows about "IDs" and nothing about slab types and numbers.  Looking
up a slab's file name requires the type and number, though.  We could work around this,
but it would be ugly.

When we call "write" on a ramdisk slab, we query `AA` to decide if it's actually a ramdisk
slab.  If so, we don't pass the command to the IO module.  We could teach the IO module
to skip such commands, but the release of the slab might then get held up behind real read/write
commands. It's not that hard to teach `SlabBuffer` to just deallocate and mark the slab as IO complete.

Same goes for reading!  No sense in passing it through the IO module.

The double-listing of some slabs (like the merge slabs) in `ReadSlabPath` and `WriteSlabPath`
doesn't seem ideal.  But it seems to be the best way to support the multipole recovery pattern
of reading `PosSlab` directly in `MergePosSlab`.  A fully flexible solution would allow reading
of arbitrary paths (or slab types) into any slab type, but the plumbing for that would be quite
verbose.  
And any bugs arising from the double listing (like accidentally reading
instead of writing) should come to light quickly.


## convolution

The convolution ramdisk interface is fairly simple.  It is all handled through the block IO
manager, with one exception (see below).  A block is a number of rows from all files.  The
block IO manager checks whether the multipoles and Taylors directories are on ramdisk (they must
either both or neither be) and treats reads and writes appropriately.

When a block read is
requested, if we are overwriting Taylors onto multipoles, the memory is mapped directly.
If we are not overwriting, then the new Taylors shared memory is created (note: this means
we are creating the output file during the read!) and the multipoles are read into this new
shared memory.  This is effectively a `memcpy()`.  The rest of the code then operates on this
new shared memory block as if this was an overwrite.

In either case, the next step after the read is the "swizzle".  The swizzle reads a single plane
from the shared memory into a new buffer.  The convolution
operates on this buffer, and the result is "unswizzled" into the original location in the
shared memory.  Thus, we go plane-by-plane overwriting Taylors onto multipoles.

At the end of the convolution, if we were overwriting we thus have CPD files named `Multipole_XXXX` that are in
fact not filled with multipoles but Taylors!  The final step of convolution is thus to
rename these files.  This is always true if we are overwriting, regardless of whether it's a
ramdisk overwrite.

The exception to the rule that only the block IO manager knows about the ramdisk is when
we select the `zwidth`.  The `zwidth` is the number of planes to read at once.  If we
are using the ramdisk, then we know the problem fits in memory.  Thus, we want to set
`zwidth` to process the whole file at once.  Even more importantly, operating in blocks
requires reading at arbitrary offsets, which is tricky given the page-alignment requirements
of `mmap()` (and not currently supported by the convolution code).

It is highly convenient to select the zwidth in the top-level code, so for now we're stuck
with this exception to the rule.  The top-level code will override the user-selected `zwidth`
if a ramdisk is in use.

Since the ramdisk is handled through the block IO manager, one can either use threaded or blocking
IO.  One probably doesn't need an IO thread if doing a ramdisk overwrite, but it also probably
won't hurt—all operations will be essentially instantaneous, so one can safely put a compute
thread on the same core as the IO thread.

### Design justifications

One could imagine making all the Taylors files symlinks to the multipoles files (or vice versa),
so no renaming step would be necessary.  But then one can't check the success of the convolution
easily; it will always look like it was successful.  We rely on multipoles file sizes to trigger
multipole recovery, too, so we don't want to trick that mechanism by linking to Taylors.

We could also avoid the renaming by using some unified name like `MultipoleTaylors_XXXX`, but
that would have the same difficulties as the symlinking.


## General Notes on the Linux Ramdisk Interface

The basic steps of writing to a ramdisk from C/C++ are:

1) Open a file descriptor to the shared memory on ramdisk with `open()`
2) Set the desired size of the shared memory with `ftruncate()`
3) Get a pointer to the shared memory with `mmap()`.

This will return a pointer just like `malloc()` would.  One can write to this memory and
immediately see the changes reflected in the corresponding file under `/dev/shm/`.

Don't forget to un-map the memory with `munmap()` when you're done!  Shared memory is only
freed back to the system when all mappings and the `/dev/shm` file handle are deleted.  Of
course, all mappings will be removed when a process exits.

`open()` and `shm_open()` both open a file descriptor that one can pass to `mmap()`.
`shm_open()` is designed to be a bit more portable and takes a "name" instead of a full file
path.  This doesn't allow for subdirectories though, so in practice we just use `open()`.
Accordingly, we delete files with `unlink()` instead of `shm_unlink()`.

Ramdisk mappings must start at page boundaries.  This only matters for covnolution, where
we operate in chunks and thus start reading/mapping at progressive file offsets.  In practice,
we enforce mapping the whole file at once if we're on ramdisk which effectively enforces 0 offset.

### Errors

One can receive the unusual `SIGBUS` (bus error, signal 7) while using a ramdisk.  Usually this
just means the ramdisk has run out of space.

## Footnotes
[↩](#a1) <b name="f1">1:</b> This path seems quite standard.  From https://www.kernel.org/doc/Documentation/filesystems/tmpfs.txt: 
> glibc 2.2 and above expects tmpfs to be mounted at /dev/shm for
> POSIX shared memory

[↩](#a2) <b name="f2">2:</b> From http://man7.org/linux/man-pages/man7/shm_overview.7.html:
> POSIX shared memory objects have kernel persistence: a shared memory
> object will exist until the system is shut down, or until all
> processes have unmapped the object and it has been deleted with
> shm_unlink(3).

## Resources
- https://github.com/abacusorg/abacus/issues/215
- https://github.com/abacusorg/abacus/pull/218
- http://man7.org/training/download/posix_shm_slides.pdf
- http://man7.org/linux/man-pages/man7/shm_overview.7.html
- http://man7.org/linux/man-pages/man3/shm_open.3p.html
