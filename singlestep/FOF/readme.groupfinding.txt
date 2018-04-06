This is documentation for the group-finding.


Conceptual Strategy:

We are implementing the usual Friends-of-Friends (FOF) algorithm in a
fast multi-threaded method using the slab and cell decomposition of the 
Abacus code.  The algorithm consists of 5 steps: 

1) In each cell, do the usual FOF algorithm, resulting in CellGroups.
We permute the particles so that each CellGroup is a contiguous run
of particles.  Singlets are CellGroups too.  Mark the CellGroup as
not closed.

CellGroups have a unique ID number called a LinkID.  The LinkID can
be parsed to return the cell index and CellGroup number within that
cell.  Further, the LinkID are sortable, with the property that sorting
collects them into cell-contiguous order.  There is also a LinkID
value that a CellGroup can never generate, which will be used to mark
Links for deletion.

2) For each cell, search the 13 neighbors (choosing an orientation so
that all 26 neighbors are considered by one of the pair).  Find all
pairs of CellGroups that have at least one pair of particles within
the FOF linking length.  For each, we create a GroupLink, which is the
pair of LinkID's, and the recipricol GroupLink, just reversing the pair.  
We add both to a global list of GroupLinks.

3) Sort the GroupLink list by the first LinkID.  Scan the list to
find, for each cell, the starting index and the segment length for
the GroupLinks whose first LinkID belong to that cell.  In that way,
when given a cell, we can quickly locate all GroupLinks that involve
a CellGroup that is in that cell.  Recall that each GroupLink has its
recipricol in the list, so both ends of the link will be present in the
catalog.

4) Traverse the GroupLink list to associate CellGroups into
GlobalGroups.  This is just like FOF, but with the linked objects
being explicitly enumerated in the GroupLink list.  In particular,
we proceed sequentially through the CellGroups.  If a CellGroup is
not yet closed, then we use it as the nucleation point for a Global
Group.  We create a list of LinkIDs, starting with this first
CellGroup.  For the next element in the list, we unpack the cell
number from the LinkID and retrieve all of the GroupLinks involving
that cell.  We search that list for all links originating from the
current LinkID.  If the destination of the link is not already in
our list, then we add it to the end.  Either way, we mark the
GroupLink for deletion by changing the number of the first entry
to the deletion value.  When this search is complete, we mark the
CellGroup as closed.  We then continue to the next element in the
LinkID list, until there are no more elements.  The resulting list
of LinkID (and hence CellGroups) is stored as the next GlobalGroup.

When this process is complete, one does a partition sort on the 
GroupLink list to move all of the elements marked for deletion
to the end, and then shrinks the list to delete them.

Note that each CellGroup appears in only one GlobalGroup, and hence
is closed only once.  And each GroupLink is considered only once, 
after which it is marked for deletion.

As a simple two-group example, CellGroup A would discover its link to 
B and add B to the list.  Then when it is B's turn, we would discover
the link to A, but find that A was already on our list.  At that 
the process would end, with both A & B marked as closed and the
two links, A->B and B->A, marked for deletion.

5) Use the GlobalGroups.  To access the particles, one uses the
list of LinkIDs to fetch the CellGroup, which in turn points to the
contiguous segment of particles that belong to the group.  This 
gathering has a well-defined order and therefore can be easily 
reversed to scatter the particles back to the cells if one has
altered them.


Slab-based Processing

This processing can happen in a slab-based manner.  Step (2) can
happen for a slab when it and one neighbor have completed step (1).
In practice, we do (3) and (4) together.  If a GlobalGroup is
required to span no more than 2*R+1 slabs (where R is the GroupRadius),
then we have to have completed (2) over all slabs from -2*R to +2*R
around a chosen one.  We then process (4) using only CellGroups
from the central slab as nucleation points.  This will result in
all CellGroups in that slab being closed and all GroupLinks involving
that slab being deleted (as well as some CellGroups in other slabs).
One then uses the GlobalGroups (step 5) on that slab.  After that,
one can discard the CellGroups from that slab.  

At the end of the processing, the GroupLink list should be empty!
Also, it is an error if a GlobalGroup ever tries to close a CellGroup
that is already closed; that probably means that a GlobalGroup has
percolated to span more than 2*R+1 slabs.

Note that the GlobalGroup "belongs" to the first-closed slab that 
contains one of its CellGroups.  That's *not* the center of mass.


Implementation

There are a number of optimizations to try to make this fast
and understandable.


Finding CellGroups:

This is the usual FOF algorithm: one keeps a list of all friends
found, priming the pump with the first particle in the yet-unassigned
set.  For the next particle in the list, one searches all unassigned
particles to see if they're closer than distance b.  If yes, add
to the end of the friend list.  When one reaches the end of the
list, that FOF group is done and one uses the next unassigned
particle to start the next group.  The pending list refers to all
friends that remain in the list but that haven't yet been searched.

As is common, one does this in place: the current particle and the
end of the friend list are just indices in the list of particles.
When one finds a new friend, one increments the end pointer and
swaps the particle there with the friend.  Completed groups accumulate
at the front of the list and one is done when one gets to the end.

We accelerate the search by further partitioning the pending and
unassigned lists.  For a first particle, we compute the distance
to all pending and unassigned particles.  Pending particles closer
than R are swapped to the front of the pending list.  Unassigned
particles are separated into <R, R to R+b, and >R+b sets.  This is
a small extra time for this particle search.  But having done this,
we can search the <R pending list knowing that only unassigned list
out to R+b need to be considered.  New friends found from the
unassigned <R list are added to the pending <R list; those from R
to R+b are put on the pending >R list.  Only when we have exhausted
the pending R<b list do we do another partitioning, this time
starting from the first particle on the pending R>b list.

The choice of R is heuristic and is allowed to adapt so that that the
size of the nearby lists is kept small, but not too small.  One
wants to have enough nearby particles to get some re-use of the 
partitioning work, but not so many that the list of unassigned 
particles to search is getting big.  It's possible that a better 
heuristic would offer more speed.  For example, we're not considering
the ratio of the sizes of the nearby to distant unassigned particle sets.

In practice, we find that this is faster for cells with more than
about 70 particles.  For smaller groups, we revert to an algorithm
that doesn't do partitioning.  This has less overhead from bookkeeping,
and the AVX-enhanced distance calculation is fast enough that N^2 
doesn't kill us.

The primary advantage of this partitioning algorithm is that when one
has many particles in one group as well as many particles exterior 
(e.g., in a second group), one avoids having to do the full N^2 
comparison of all particles in the two groups.  Having done one 
partition, one will be able to resolve many other nearby cases
without touching the distant unassigned set.

When processing a cell, we first copy the positions into SSE-aligned
float4's called FOFparticle, placing the initial index position
into the 4th element.  We permute these FOFparticles as the CellGroups
are being built, rather than swapping the larger pos/vel/aux lists.
At the end, we use the indices to permute the input pos/vel/aux
arrays into the final order. 

In practice, we have chosen to adjust the final ordering so that 
all multiplets are first, followed by all singlets that are within
b of the edge of the cell, followed by all singlets that are 
in the interior of the cell.  This is not essential and might be
re-considered. 

We use AVX to compute the square distance to 4 particles at once,
returning a vector of 4 results.  This is useful in FOF, since for
a given particle, we need to search all pending particles to find
new friends.  We multiply the float3 positions by a large number
(1e12), so that we can compute float4 norms to get the distance
between two particles.

We also track the Cartesian XYZ min/max bounding boxes for each group.
This is used to decide which groups are close to the cell edges.

Each cell is independent, so we give each z-pencil of cells to a thread.
This algorithm runs in its own class, called FOFcell, for which there
is one instantiation per thread.

After run, we translate the results into a list of classes called
CellGroups, which contain the starting particle index in the cell
and the number of particles in the cell group.  We then overload
the upper 7 bits of the number for later processing: 6 bits are
used, being set if the CellGroup Bounding Box (or the singlet
particle position) is within b of each of the 6 faces.  The 7th bit
says whether the group is open or closed.  This means that a CellGroup
can not exceed a 25-bit number, which is multiplicity of 33 million.
This is sufficient for the anticipated use in Abacus.

The FOF implementation is in fof_sublist.cpp, while the CellGroup 
code is in cellgroup.cpp and in groupfinding.cpp.

There is a unit test for the fof_sublist.cpp code that sets up a
toy distribution in a cell and runs FOF, with timings.  One can
compare the results between different algorithms, including a
bare-bones reference C implementation.


Dynamical Cell-indexed Lists

Unlike most other cases in Abacus, in the group finding, we cannot
predict the size of the slab in advance.  This leads to us needing
a way to build multi-threaded lists while retaining the cell grid.
We do this in the templated SlabAccum<T> class, in the slab_accum.cpp
file.  Each instance will contain the objects of type T for that slab.

This code is very rigid about z-pencils.  Each z-pencil must be
handled by a single thread and must be processed in z=[0,cpd) order.
Further, it is built to first have a build phase, in which one
appends objects of type T to the current cell, and then a read-only
use phase.  The methods to access information during the appending
are limited -- really only the method to return the current size
of the pencil has been tested.

But the underlying class can accept arbitrary appends without needing to 
pre-specify the total memory.  The data for individual pencils is kept
contiguous; multiple pencils are not guaranteed contiguous (nor indeed
in any order).  The class builds a grid of start indices and sizes for
all of the cells in the slab.

For the read-only phase, the class provides a series of [] operators
so that one can access a pencil as s[j], a cell as s[j][k], and an
object in the cell as s[j][k][n].

The group finding code uses a fair number of these SlabAccum classes
to store the CellGroups and GlobalGroups, as well as the various
intermediate objects needed to search for GroupLinks.

There is a unit test for this code that loads a SlabAccum<int> with a
big set of objects and verifies that each object is found where it 
should be in the read-only phase.  This stops short of testing every
possible use of the code, however.


Finding Links between Cell Groups:

This is potentially expensive, because every cell has 26 neighbors
and therefore requires 13 directional searches.  If we had to do 
every particle, the N^2 would dominate the CellGroup FOF.

But of course we don't: we only need to consider particles near 
the boundaries.  The code implements this restriction only for
faces, not edges or corners.  We expect that the search even for
a face will be fast enough that the gain of further cutting down that 
search would overwhelm the overheads of making the new set and 
scanning through very short lists.  A single face has only ~10% 
of the particles on average.  

We precompute the applicable faces (separately for each of the 6
directions) for this slab and for the directionally adjacent slab.
Each x-based face will get used 9 times; each y-based face will get
used 3 times; each z-based face gets used once.  For the primary
slab, we're making 5 faces -- these get done for each cell,
so that we get to re-use the particles before they leave cache.

To make a face for a given cell, we scan all CellGroups that were
marked as having a bounding box overlapping this face.  The particles
in each CellGroup that are in fact within b of the cell boundary
are called faceParticles.  If a group has only one faceParticle,
then is promoted to a pseudoParticle with zero Radius and with an
index equal to the CellGroup number.  If a group has more than one
faceParticle, then we create a faceGroup to point to the start and
length of that chunk of faceParticles and to record the CellGroup
number.  We compute the Cartesian min/max bounding box of the 
associated faceParticles; we take the midpoint as a pseudoParticle,
with index equal to the faceGroup number, and the semi-diagonal 
length as the Radius.

Then searching between two faces is quick: we consider all pairs
of pseudoParticles and require the relative distance to be less
than b plus the sum of the radii.  If either or both is actually a
faceGroup, then a successful first pass means that we have to open
up the faceParticles, looking at pairs until we find a pair that
is closer than b.  Any successful pair triggers the addition of a
GroupLink and its recipricol to the global list.

Searching must account for the cell-centered position convention.

In detail, pseudoParticles and faceParticles are both just typedef'ed
to FOFparticles, so as to re-use the SSE alignment.  But this is not 
essential and might be simplified.

This code is in findgrouplink.cpp.


Managing the GroupLink list:

The LinkID class is a 64-bit integer that is a mashup of the cell (i,j,k) 
(12 bits each) and the CellGroup number (20 bits).  There are methods
to get info in and out, plus a function to return the corresponding
CellGroup from the global list of those.

The GroupLink class is a simple pair of those, plus a comparison operator
for sorting.

The GroupLinkList class is built on the template of a MultiAccumList,
which reuses the InsertList infrastructure.  This allows multiple 
threads to append to a list.  At present, there is no provision to
expand the underlying storage if we run out, so the initial allocation
has to be generous and wise.  We could alter that if needed.

This code is in grouplink.cpp and multiappendlist.cpp.
There is a unit test embedded in grouplink.cpp.


Forming Global Groups:

To close all of the groups on a given slab, we sort the GroupLinkList and
then identify the starting locations for cells in this slab and 2*R adjacent.
This indexing is done by doing a binary search to the start of each pencil,
and then having each thread do a linear sweep to set up the start and
extent of the GroupLinks for each cell.

Then we split the slab into pencils into sets that are at least
4*R+1 apart.  This requires some care around the periodic wrap.
Each pencil is assigned to one thread, which then works cell by
cell.

In a given cell, we prime a LinkID list with the first non-closed CellGroup, 
then repeatedly search the GroupLinkList for relevant GroupLinks from this
CellGroup.  Destination LinkIDs are added only if they are new to the list.
These searches over the LinkIDs in a Cell and over the current LinkID list
are just done as linear searches; we expect that the number of GroupLinks
in a cell will be small enough that it's faster just to scan through rather
than try for some bisection.

As CellGroups are considered in the list, we mark them as closed.
All GroupLinks considered are marked for deletion.

GlobalGroups end up represented as two SlabAccum's, assigned to the
cell of the first CellGroup in the list (i.e., the nucleation point
that closed the GlobalGroup).  The first contains the LinkID's for
the CellGroups, contiguous for each GlobalGroup.  The second contains
the start and length of that segment, offset from the start of the
pencil.  It also contains the number of particles in the group and
a cumulation of number of particles since the beginning of the
pencil (not including this group).  This latter is later updated to
the number from the beginning of the slab.

This code is in globalgroups.cpp.


Using Global Groups:

We provide a GatherGlobalGroups() function to copy all of the particles
from these groups into contiguous arrays.  This includes adjustments
to bring the positions into the cell-centered frame of the first CellGroup.
The start index and length of these segments of particles is available
in the GlobalGroup class.  This is in groupfinding.cpp.

We will provide a similar scatter routine to return the particles.

Eventually, there may be reason to copy these into multiple arrays,
separated by group multiplicity and the destination for microstepping.

We will likely seek to run sub-group finding on these groups as well;
that will require an alteration to the FOFcell class, as we probably
don't want to do the final permutation of particles, but rather just
chase index pointers to compile stats.  But the plan is to this directly
from the contiguous lists, not return to Cells and Links, as the 
partitioning FOF should be ok for performance.


Controlling the Group Finding:

There is a class GroupFindingControl that should have one global
instance.  This contains the basic external parameters needed by
the code, as well as the GroupLinkList and pointers to [0,cpd)
pointers to slabs of CellGroups.  This is in groupfinding.cpp.

There will need to be a set of Dependencies written.


Timing

The FOFcell class has timings embedded in it.  Nothing else has been
timed yet.


Test code:

There is a test_driver.cpp code that sets up a simple toy problem.
This does give a reproducible answer, and one can check that the
resulting particles are close together.  That said, it stops short
of a full test, either that all found groups are real or that 
all real groups are found.  For that, we likely want to put in a
particle distribution that has been run with an external periodic
FOF code and check that the results agree.

