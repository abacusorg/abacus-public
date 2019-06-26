'''
This is an interface to write particle data in some of our Abacus
binary formats.

TODO: add a pack14 writer

See ReadAbacus.py for equivalent interfaces to read these formats.
'''

import numpy as np
import numba

#import line_profiler

def write(outfn, particles=None, pos=None, vel=None, pid=None, ppd=None, write_zel=False, write_pid=False, append=False, dtype=np.float32):
    '''
    Write particle data (pos, vel, pid, etc).  The format is compatible with the
    Zeldovich/RV* family of formats.  See the Abacus readme for details.

    Parameters
    ----------
    particles: ndarray, optional
        The array containing the particles, ready to write.
        No copies will be made and the array will be written
        as-is to disk.  This overrides the other options.
        Default: None

    '''
    if particles is not None:
        a = particles
    else:
        # Input validation
        assert not (write_zel and write_pid)  # none of our formats use both
        if write_zel or write_pid:
            assert pid is not None
        if write_zel:
            assert ppd

        pos = np.atleast_2d(pos)
        assert pos.shape[-1] == 3

        if vel is not None:
            vel = np.atleast_2d(vel)
            assert vel.shape == pos.shape
        if pid is not None:
            pid = np.atleast_1d(pid)
            assert len(pid) == len(pos)

        # Construct the output array
        base_dt = dtype; del dtype
        dt = []
        if write_zel:
            dt += [('ijk',np.uint16,3)]
        dt += [('pos',base_dt,3)]
        if vel is not None:
            dt += [('vel',base_dt,3)]
        if write_pid:
            dt += [('pid',np.uint64)]
        dt = np.dtype(dt, align=True)
        a = np.empty(len(pos), dtype=dt)
        
        # Populate the output array
        a['pos'] = pos
        if vel is not None:
            a['vel'] = vel
        if write_pid:
            a['pid'] = pid
        if write_zel:
            a['ijk'][:,0] = pid // ppd**2
            a['ijk'][:,1] = (pid % ppd**2) // ppd
            a['ijk'][:,2] = (pid % ppd**2) % ppd

    # Write the output array
    mode = 'b'
    if append:
        mode = 'a' + mode
    with open(outfn, mode) as fp:
        a.tofile(fp)


# x-coord to slab
# warning: this function may give different answers for scalar and array arguments
# because numpy upcasts scalars, but not arrays
@numba.njit(parallel=True)
def x2s(x, box, cpd):
    N = len(x)
    slab = np.empty(N, dtype=np.int64)
    
    for i in numba.prange(N):
        slab[i] = (x[i]/box + .5)*cpd % cpd
        assert 0 <= slab[i] < cpd
    return slab


class SlabWriter:
    '''
    Sometimes we have particle data that is not organized according to Abacus slabs (or we want to change the
    CPD of an existing set of files).  This class takes these files one at a time, sorts the particles by
    slab, and writes each contiguous chunk to disk.  Thus only one file has to be in memory at once.

    TODO: we could wait until filling a buffer to write particles to disk.  That would mean larger chunks
    of IO, which could be more efficient, at the cost of using more memory.

    See `change_cpd.py` or `gadget_to_abacus.py` for example usage.

    '''
    # TODO: clean up format support with new write function
    writers = {'RVdouble': lambda *args,**kwargs: write(*args,Vel=True,Double=True,**kwargs),
               'RVdoubleZel': lambda *args,**kwargs: write(*args,Vel=True,Double=True,Zel=True,**kwargs),
               'RVZel': lambda *args,**kwargs: write(*args,Vel=True,Zel=True,**kwargs),
               'RVTag': lambda *args,**kwargs: write(*args, write_pid=True, **kwargs),
               'Zeldovich': lambda *args,**kwargs: write(*args,Zel=True,Double=True,**kwargs),
               'same':lambda *args,**kwargs: write(*args,**kwargs)
               }

    def x2s(self, x):
        return x2s(x, self.boxsize, self.cpd)

    # Assigns particles to slabs
    def ingest(self, particles=None, pos=None, vel=None, pid=None, slabinfo=None):
        if particles is not None:
            assert pos is None and vel is None and pid is None
            pos = particles['pos']
            vel = particles['vel']
            pid = particles['pid']
            # maybe we should forcibly pack together arrays passed separately?

        # Do we already know the pos->slab assignments?
        if slabinfo:
            slabs,s1,s2 = slabinfo
        else:
            slabs = self.x2s(pos[:,0])
            s1 = slabs.min()
            s2 = slabs.max()
            if self.verbose:
                print("\tHas {} particles spanning slabs {} to {}".format(len(pos),s1,s2))

        # Check for the recursion base case: all positions within one slab
        if s1 == s2:
            n = len(pos)
            outfn = self.output_fn_fmt.format(s1)
            if particles is not None:
                self.writer(outfn, particles=particles, ppd=self.ppd, append=True)
            else:
                self.writer(outfn, pos=pos, vel=vel, pid=pid, ppd=self.ppd, append=True)

            self.num_written[s1] += n
        else:
            # If partitioned, then split into slabs and recurse
            if (slabs[1:] >= slabs[:-1]).all():
                splits = np.where(slabs[1:] != slabs[:-1])[0] + 1  # indices of the first particle in each new slab
                splits = [0] + list(splits) + [len(pos)]
                for lo,hi in zip(splits[:-1],splits[1:]):
                    if particles is not None:
                        self.ingest(particles=particles[lo:hi],
                                    slabinfo=(slabs[lo:hi], slabs[lo], slabs[lo]))
                    else:
                        self.ingest(pos=pos[lo:hi],
                                    vel=vel[lo:hi],
                                    pid=pid[lo:hi],
                                    slabinfo=(slabs[lo:hi], slabs[lo], slabs[lo]))

            # Not sorted; sort and recurse
            else:
                inds = np.argsort(slabs)
                slabs = slabs[inds]

                if particles is not None:
                    particles = particles[inds]  # this is the performance bottleneck, could make parallel
                    self.ingest(particles=particles, slabinfo=(slabs,s1,s2))
                else:
                    pos = pos[inds]
                    vel = vel[inds]
                    pid = pid[inds]
                    self.ingest(pos=pos, vel=vel, pid=pid, slabinfo=(slabs,s1,s2))
                
    def __del__(self):
        if self.num_written.sum() != self.NP:
            raise RuntimeError('Only wrote {} particles; expected {}'.format(self.num_written.sum(), self.NP))
        elif self.verbose:
            print('\tFinished; wrote {} particles.'.format(self.num_written.sum()))
                
    def __init__(self, NP, cpd, boxsize, outdir, format, verbose=True):
        self.NP = NP
        self.ppd = int(round(NP**(1./3)))
        if self.ppd**3 != NP:
            self.ppd = None
        self.cpd = cpd
        self.boxsize = boxsize
        self.num_written = np.zeros(cpd, dtype=np.int)
        self.verbose = verbose
                
        self.output_fn_fmt = outdir+"/ic_{}"
        self.writer = SlabWriter.writers[format]
