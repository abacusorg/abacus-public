<h1 align="center">abacus</h1>

<div align="center">
<img src="doc/icon_red.png" width="100px" alt="Abacus Logo">
</div>

<div align="center">
A high-performance N-body code for cosmological simulations.

[![GPU Tests](https://jenkins.flatironinstitute.org/buildStatus/icon?job=abacus%2Fdevelopment&subject=GPU%20Tests)](https://jenkins.flatironinstitute.org/job/abacus/job/development/)
</div>

## Obtaining the Code
Clone the code from GitHub, including submodules recursively:
```console
$ git clone --recursive git@github.com:abacusorg/abacus.git
```

If you forget the recursive clone, you can initialize the submodules later:
```console
$ git submodule update --init
```

The submodules are all located in `external/`.

## Dependencies
Abacus requires:
- A C++17 compiler (e.g. GCC >= 8, Intel C++ Compiler Classic >= 19, or any Intel oneAPI C++ Compiler)
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
$ cmake -B build -G Ninja --install-prefix=install -Dgperftools_build_minimal=TRUE -Dgperftools_tcmalloc_pagesize=256
$ cmake --build build
$ cmake --install build
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
$ make
```

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

## Papers
The [AbacusSummit Papers & Citation page](https://abacussummit.readthedocs.io/en/latest/citation.html) lists all the relevant Abacus papers.

## Related Projects
- [abacusutils](https://github.com/abacusorg/abacusutils): Python code to interface with halo catalogs and other Abacus N-body data products
- [AbacusSummit](https://abacussummit.readthedocs.io/): information about the AbacusSummit simulation suite
