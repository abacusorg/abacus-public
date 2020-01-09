# Debugging Abacus with Valgrind

Abacus behaves fairly well under valgrind, with a few caveats.
If one is using multiple threads (i.e. for IO and GPU directs),
then one should use valgrind's `--fair-sched=yes` option, which will
prevent the main thread from hogging the (single) valgrind
execution stream.

`tcmalloc` may also confuse valgrind (see https://github.com/gperftools/gperftools/issues/792).
One may turn this off by configuring with `./configure --disable-tcmalloc`, or
removing the `-ltcmalloc_minimal` option from
the build options in `common.mk`, or linking against `-ltcmalloc_minimal_debug`
instead.  Be warned that this debug library is very slow, even outside
of valgrind!

The ParseHeader library also does many operations
that make valgrind very upset (produce a lot of warnings).
These warnings are harmless as far as we've been able to tell.
The [`valgrind.supp`](https://github.com/abacusorg/abacus/blob/master/valgrind.supp)
suppressions file in the root of the Abacus repository should hide
most of these. (TODO: haven't checked the suppressions file
in a while, may need updating).

Usually one will want to invoke valgrind directly on the 
`singlestep` or `convolution` executable instead of the
Python wrapper, although it is possible to do either.  One should
use `--trace-children=yes` in this case to follow the exectuables
that the Python wrapper spawns.  Python will likely generate additional
warnings for which one will want a suppressions file; such files
are available online for different versions of Python
(e.g. https://svn.python.org/projects/python/trunk/Misc/valgrind-python.supp).
Similarly, CUDA generates a lot of warnings and may need a suppressions file.

Valgrind currently does not support AVX-512, so one may need to remove `-march=native`
from the `gcc` compilation invocation in `common.mk` (or `-xHost` for `icc`).
Issue tracker for AVX-512 support here: https://bugs.kde.org/show_bug.cgi?id=383010




A sample valgrind invocation for `singlestep` might be:
```
valgrind --leak-check=yes --fair-sched=yes --suppressions=$ABACUS/valgrind.supp $ABACUS/singlestep/singlestep /path/to/abacus.par 0
```
