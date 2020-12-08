The New Ewald Test
==================

The "Ewald" test has been our stock-standard Abacus force accuracy test since z~inf.
The New Ewald test is the same test, with a re-written top-level Python
driver applying modern Abacus best practices.

It may be worth developing extensions to run with more and/or arbitrary particle
distributions.  This would allow us to run Ewald in parallel with realistic ppc.
The homogeneous lattice test mostly fills this need, but has no intrinsic force
scale against which to compare the results because the right answer is zero!  So
it's a little harder to quote a result.

The `reference_ewald_<pos|acc>.double3` files are a renaming of the original `posinfile`
and `ewald` files.

Usage
-----
```console
$ ./run_newald.py --help
```
