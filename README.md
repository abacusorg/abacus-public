<h1 align="center">abacus</h1>

<div align="center">
<img src="doc/icon_red.png" width="100px" alt="Abacus Logo">
</div>

<div align="center">
A high-performance N-body code for cosmological simulations.

[![GPU Tests](https://jenkins.flatironinstitute.org/buildStatus/icon?job=abacus-public%2Fmain&subject=GPU%20Tests)](https://jenkins.flatironinstitute.org/job/abacus-public/job/main/)
</div>

## Obtaining the Code
Clone the code from GitHub, including submodules recursively:
```console
$ git clone --recursive git@github.com:abacusorg/abacus-public.git
```

If you forget the recursive clone, you can initialize the submodules later:
```console
$ git submodule update --init
```

The submodules are all located in `external/`.

## Dependencies
Abacus requires:
- A C++17 compiler
- CUDA Toolkit (for GPU directs)
- GSL
- tcmalloc (bundled in `external/gperftools/`)
- TBB or oneTBB (bundled in `external/oneTBB`)
- flex and bison, for `ParseHeader`

The Abacus Python interface requires a number of (mostly) standard Python packages. It's a good idea to install them in a `venv`:
```console
$ python -m venv myvenv
$ . myvenv/bin/activate
$ pip install -r requirements.txt
```

Conda is not recommended, as it may end up installing libraries that are incompatible with the compilers used to build Abacus.

A lot of these instructions will change if/when we port to a CMake/Meson structure.

### zeldovich-PLT
The IC code [zeldovich-PLT](https://github.com/abacusorg/zeldovich-PLT/) is bundled with Abacus. It needs to be compiled separately (optionally with buffering of state on disk):
```console
$ cd external/zeldovich-PLT
$ meson setup build [-DDISK=true]
$ meson compile -C build
```

### oneTBB
oneTBB is bundled with Abacus for convenience. One can use an external installation or module, or the bundled version. The bundled version can be built with:
```console
$ cd external/oneTBB
$ cmake -B build -G Ninja -DTBB_TEST=OFF
$ cmake --build build
```

The bundled oneTBB environment variables can be set by `vars.sh`, which lives in a directory similar to the following:
```console
$ . oneTBB/build/gnu_12.3_cxx11_64_relwithdebinfo/vars.sh
```

Be sure to replace `gnu_12.3_cxx11_64_relwithdebinfo` with your platform directory!

### tcmalloc
tcmalloc is bunded in `external/gperftools`. One can use an external installation or module, or the bundled version. The bundled version can be built with:
```console
$ cd external/gperftools
$ ./autogen.sh
$ ./configure --prefix=$(pwd)/install --enable-minimal --with-tcmalloc-pagesize=256
$ make && make install
```

The following environment variables should be set:
```console
$ export    LIBRARY_PATH=$ABACUS/external/gperftools/install/lib64:$LIBRARY_PATH
$ export LD_LIBRARY_PATH=$ABACUS/external/gperftools/install/lib64:$LD_LIBRARY_PATH
$ export CPATH=$ABACUS/external/gperftools/install/include:$CPATH
```

The `$ABACUS` variable needs to point to the repo root (see [Environment Setup](#environment-setup)).

## Environment Setup
### Prepare an environment script
The `env` directory contains example environment files from different clusters. Copying and modifying one of them will probably be helpful in setting up a working environment. We recommend copying it to the repo root for convenient activation. For example:

```console
$ cp env/perlmutter.sh env.sh
$ . env.sh
```

Remember to replace `perlmutter.sh` with the file relevant to your cluster.

You'll very likely need to modify the environment file, as cluster software environments change regularly.

> [!TIP]
> Keeping the Abacus environment in its own file and activating it as needed is highly recommended, rather than putting it in your shell startup files (e.g. `~/.bashrc`).

### Symlink the site file
We keep performance parameters and other options on a given system that are commonly shared between simulations in "site files". Users should create a symlink called `site.def` that links to the site file for the current system as follows (using perlmutter as an example):

```console
$ cd Production/site_files
$ ln -s perlmutter.def site.def
```

Be sure to replace `perlmutter.def` with the appropriate site file! If one doesn't exist, you can copy and modify a site file for a similar system.

Note that some site files have `Parallel = 1` or `Parallel = 0` to make running tests more convenient. But you will need to toggle this depending on whether Abacus was compiled with MPI or not.

## Building the Code
To build the non-MPI version of the code, suitable for running on a single node:
```console
$ ./configure
$ make
```

For the MPI code:
```console
$ ./configure --enable-parallel CXX=mpicxx
```

where `mpicxx` is the MPI compiler wrapper (if applicable). This may be different (or not needed) on your system, e.g. `CXX=CC` for Cray compiler wrappers. Alternatively, one can set `CXX` in the shell environment and that will be respected.

The configure script will print out a block of configuration options. Make sure everything looks as expected! More options can be seen by running `./configure --help`.

To start over at any point, use `make distclean`.

## Running a Simulation

The `Production` directory contains examples and instructions for running simulations. See [`Production/README.md`](Production/README.md).

## Tests
Tests are in the `Tests` directory. See [`Tests/README.md`](Tests/README.md).

## Source Layout

* **singlestep**: the primary abacus timestepping code

* **python**: drivers for the code and related utilies

* **Test**: code for small test cases

* **Convolution**: code to compute the multipoles, Taylors, and convolutions

* **Derivatives**: code to generate derivatives used by the far force

* **ParseHeader**: library to read/write headers on each output file

* **include**: small utility functions (timing, threevectors, files, cosmology)

* **zeldovich**: the Zel'dovich initial condition generator

* **Analysis**: code for generating power spectra, running halo finders, etc

* **clibs**: libraries that codes may want to link against for various analysis tasks

* **Configuration**: files related to the autoconf build system

* **Production**: various configurations used for major production runs

* **util**: utilities to deal with the pack14 output format

* **env**: example environment files

## Documentation Layout

* **doc**: Code documentation

* **notes**: Musings about the algorithms, computers, etc.

* **papers**: Journal & conference papers (make subdirectories for each one)

## Related Projects
- [abacusutils](https://github.com/abacusorg/abacusutils): Python code to interface with halo catalogs and other Abacus N-body data products
- [AbacusSummit](https://abacussummit.readthedocs.io/): information about the AbacusSummit simulation suite

## Citation
If you use or reference Abacus in a published work, please cite [Garrison, et al. (2021)](https://doi.org/10.1093/mnras/stab2482). A BibTeX entry in included below for convenience. Depending on which components or data are used, more citations may be appropriate. The [AbacusSummit Papers & Citation page](https://abacussummit.readthedocs.io/en/latest/citation.html) lists most of the core Abacus papers.

<details>
<summary>BibTeX entry</summary>

```bibtex
@ARTICLE{10.1093/mnras/stab2482,
    author = {Garrison, Lehman H and Eisenstein, Daniel J and Ferrer, Douglas and Maksimova, Nina A and Pinto, Philip A},
    title = "{The abacus cosmological N-body code}",
    journal = {Monthly Notices of the Royal Astronomical Society},
    volume = {508},
    number = {1},
    pages = {575-596},
    year = {2021},
    month = {09},
    abstract = "{We present abacus, a fast and accurate cosmological N-body code based on a new method for calculating the gravitational potential from a static multipole mesh. The method analytically separates the near- and far-field forces, reducing the former to direct 1/r2 summation and the latter to a discrete convolution over multipoles. The method achieves 70 million particle updates per second per node of the Summit supercomputer, while maintaining a median fractional force error of 10−5. We express the simulation time-step as an event-driven ‘pipeline’, incorporating asynchronous events such as completion of co-processor work, input/output, and network communication. abacus has been used to produce the largest suite of N-body simulations to date, the abacussummit suite of 60 trillion particles, incorporating on-the-fly halo finding. abacus enables the production of mock catalogues of the volume and resolution required by the coming generation of cosmological surveys.}",
    issn = {0035-8711},
    doi = {10.1093/mnras/stab2482},
    url = {https://doi.org/10.1093/mnras/stab2482},
    eprint = {https://academic.oup.com/mnras/article-pdf/508/1/575/40458823/stab2482.pdf},
}
```
</details>

## License
Abacus is licensed under the [GNU General Public License v3.0](./LICENSE) or any later version (SPDX-License-Identifier: GPL-3.0-or-later).

Some sources are included from external projects under the terms of those projects' licenses. These sources retain their original license notices.

### The Abacus Name

To avoid confusion in the scientific community, we kindly ask that you:

- Feel free to mention Abacus when describing work that uses this software
- Consider alternative naming schemes for derivative works that don't incorporate "Abacus" directly in the title (e.g., instead of "AbacusMocks," perhaps "Mocks with Abacus" or your own unique name)

This helps users distinguish between the core Abacus project and the community contributions built around it.
