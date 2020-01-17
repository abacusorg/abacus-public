# This is the autoconf template file for ./configure
# Run ./autogen-config.sh to re-create ./configure from this file

# This approach to Abacus's configure system was inspired by
# Athena (https://trac.princeton.edu/Athena/)

# Initialization
AC_INIT([abacus], [2.0], [], [], [github.com/abacusorg/abacus])

if test -f ./config.site; then
    echo "configure: loading site script ./config.site"
    source ./config.site
fi

# TODO: any way to use mpicxx if --enable-parallel?
AC_LANG([C++])  # does this work?

AC_CONFIG_HEADERS([./include/config.h])  # this places DEFS in a global header file instead of Makefiles
AC_CONFIG_AUX_DIR([./Configuration])
#AC_PROG_CC
AC_PROG_CXX
AC_HEADER_STDC
AC_PROG_INSTALL

# figure out intel,ppc,etc
AC_CANONICAL_BUILD

AS_CASE([$build_cpu],
  [powerpc64*], [ARCH=powerpc],
  [x86_64], [ARCH=intel],
  [AC_MSG_ERROR("Unknown CPU architecture $build_cpu")]
)
AC_SUBST([ARCH])

# These files will be generated from their Makefile.in counterparts
# All ./configure information goes into common.mk
# Makefiles must inherit any ./configure information from there
AC_CONFIG_FILES([
 common.mk
])

# Detect the underlying compiler when we are using a wrapper like mpicxx
# TODO: better way to do this?
AC_MSG_CHECKING([underlying compiler])
AC_SUBST(REAL_CXX)
REAL_CXX=$($CXX --version) 
if \grep -qP '\bicc\b' <<< "$REAL_CXX"; then
    REAL_CXX=icc
elif \grep -qP '\bXL\b' <<< "$REAL_CXX"; then
    REAL_CXX=xlc++
elif \grep -qP '\bg++\b' <<< "$REAL_CXX"; then
    REAL_CXX=g++
else
    AC_MSG_WARN([Could not determine REAL_CXX from --version])
    REAL_CXX=$CXX  # fallback
fi
AC_MSG_RESULT([$REAL_CXX])
REAL_CC="$REAL_CXX -x c"

# git version
GITVERSION=\"$(git rev-parse HEAD)\"
AC_DEFINE_UNQUOTED(GITVERSION, $GITVERSION, "The current git revision")

# Imports
m4_include(Configuration/m4/ax_check_cuda.m4)
m4_include(Configuration/m4/ax_check_compile_flag.m4)
m4_include(Configuration/m4/ax_gcc_x86_avx_xgetbv.m4)
m4_include(Configuration/m4/ax_gcc_x86_cpuid.m4)
m4_include(Configuration/m4/ax_ext.m4)

# The following uses gsl-config to find GSL.  That sounds like
# a good thing, but Anaconda puts gsl-config on the PATH, but not
# the library paths.
# So we can't resolve the libraries at link/runtime if we do this.
# Putting the Anaconda lib dir on the search path is dangerous because
# it includes many versions of system libraries, and we might worry that
# their performance (and even behavior) would be different from
# the "real" system libraries.
#m4_include(Configuration/m4/gsl.m4)
#AX_PATH_GSL(,[], AC_MSG_ERROR(could not find required version of GSL))

# Calling make from configure... pretty sketchy!
# But we want tcmalloc to be available for configure to find
# This will only be run if libtcmalloc isn't already present
make --quiet tcmalloc

# These will automatically add the appropriate -L and -l flags to LIBS
AC_SEARCH_LIBS([cos],[m])
AC_SEARCH_LIBS(cblas_dgemm, [gslcblas], [], [AC_MSG_ERROR([unable to find gslcblas!])])
AC_SEARCH_LIBS(gsl_rng_uniform, [gsl], [], [AC_MSG_ERROR([unable to find gsl!])])
AC_CHECK_LIB(numa, move_pages, [], [AC_MSG_WARN([unable to find numa!])])


# With the Intel compiler, the builtin TBB is accessed via the -tbb flag (set in common.mk)
if test x$REAL_CXX != xicc; then
    AC_SEARCH_LIBS(TBB_runtime_interface_version, [tbb], [], [AC_MSG_ERROR([unable to find tbb!])])
fi

# Check for some utilities we use
AC_CHECK_PROGS(LEX, flex)
AC_CHECK_PROGS(YACC, bison)
AC_CHECK_PROGS(AR, ar)
AC_CHECK_PROGS(SWIG, swig)

# Check for AVX
AX_EXT

# We prefer yes/no, not yes/blank
if test x$ax_cv_support_avx_ext == xyes; then
    have_avx=yes
else
    have_avx=no
fi

if test x$ax_cv_support_avx512f_ext == xyes; then
    have_avx512f=yes
else
    have_avx512f=no
fi

if test x$ax_cv_have_vsx_ext == xyes; then
    have_vsx=yes
else
    have_vsx=no
fi

# Define the ./configure options
AC_ARG_ENABLE(double,
[AS_HELP_STRING([--enable-double], [Double precision [no]])],
	[], enable_double="no")
if test x$enable_double == xyes; then
	AC_DEFINE_UNQUOTED(DOUBLEPRECISION, [], "Do all computations in double precision instead of single")
fi

AC_SUBST(TCMALLOC)
TCMALLOC="TCMALLOC_DISABLED"
AC_ARG_ENABLE(tcmalloc,
[AS_HELP_STRING([--enable-tcmalloc], [TCMalloc heap allocator [auto]])],
    [], enable_tcmalloc="auto")
if test x$enable_tcmalloc == xyes || test x$enable_tcmalloc == xauto; then
    AC_CHECK_LIB(tcmalloc_minimal, malloc, [], [AC_MSG_WARN([unable to find tcmalloc!])])
fi
if test x$ac_cv_lib_tcmalloc_minimal_malloc == xyes; then
    TCMALLOC="TCMALLOC_ENABLED"
fi

AC_ARG_ENABLE(mkl,
[AS_HELP_STRING([--enable-mkl], [Intel Math Kernel Library FFT [yes if using icc]])],
	[], enable_mkl="check")
if test x$enable_mkl == xcheck; then
    if test x$REAL_CXX == xicc; then
        enable_mkl="yes"
    else
        enable_mkl="no"
    fi
fi

AC_ARG_WITH(radius,
	[AS_HELP_STRING([--with-radius=NFR], [Near field radius [2]])],
	[], with_radius="2")
AC_DEFINE_UNQUOTED(NFRADIUS,$with_radius, "Max Radius of cells to compute with near-field force")

AC_ARG_WITH(block-size,
[AS_HELP_STRING([--with-block-size=B], [Near field block size [64]])],
	[], with_block_size="64")
AC_DEFINE_UNQUOTED(NFBlockSize,$with_block_size, "Near field block size")

AC_ARG_ENABLE(avx,
[AS_HELP_STRING([--enable-avx], [All AVX (overridden by individual options) [yes if AVX present]])],
	[], enable_avx=$have_avx)

AC_ARG_ENABLE(avx512,
[AS_HELP_STRING([--enable-avx512], [All AVX-512 (overridden by individual options) [yes if AVX-512 present]])],
	[], enable_avx512=$have_avx512f)

AC_ARG_ENABLE(vsx,
[AS_HELP_STRING([--enable-vsx], [All VSX (overridden by individual options) [yes if VSX present]])],
    [], enable_vsx=$have_vsx)

AC_ARG_ENABLE(avx-directs,
[AS_HELP_STRING([--enable-avx-directs], [AVX directs [no]])],
	[], enable_avx_directs=no)
if test x$enable_avx_directs == xyes; then
	AC_DEFINE_UNQUOTED(AVXDIRECT, [], "Use AVX directs")
fi

AC_ARG_ENABLE(avx512-directs,
[AS_HELP_STRING([--enable-avx512-directs], [AVX-512 directs [no]])],
	[], enable_avx512_directs=no)
if test x$enable_avx512_directs == xyes; then
	AC_DEFINE_UNQUOTED(AVX512DIRECT, [], "Use AVX-512 directs")
fi


# Begin multipoles SIMD options

AC_ARG_WITH(multipoles-simd,
[AS_HELP_STRING([--with-multipoles-simd={avx512 | avx | vsx | unrolled | none}], [Multipoles/Taylors implementation [highest supported]])],
    [], with_multipoles_simd="auto")

if test "x$with_multipoles_simd" == xauto; then
    if test "x$enable_avx512" == xyes; then
        with_multipoles_simd="avx512"
    elif test "x$enable_avx" == xyes; then
        with_multipoles_simd="avx"
    elif test "x$enable_vsx" == xyes; then
        with_multipoles_simd="vsx"
    else
        with_multipoles_simd="unrolled"
    fi
fi

AC_SUBST(MULTIPOLES_SIMD)
MULTIPOLES_SIMD="MULTIPOLES_SIMD_NONE"
if test "x$with_multipoles_simd" == xavx512; then
    AC_DEFINE_UNQUOTED(AVX512MULTIPOLES, [], "AVX-512 multipoles/Taylors")
    MULTIPOLES_SIMD="MULTIPOLES_SIMD_AVX512"
elif test "x$with_multipoles_simd" == xavx; then
    AC_DEFINE_UNQUOTED(AVXMULTIPOLES, [], "AVX multipoles/Taylors")
    MULTIPOLES_SIMD="MULTIPOLES_SIMD_AVX"
elif test "x$with_multipoles_simd" == xvsx; then
    AC_DEFINE_UNQUOTED(VSXMULTIPOLES, [], "VSX (PowerPC) multipoles/Taylors")
    MULTIPOLES_SIMD="MULTIPOLES_SIMD_VSX"
elif test "x$with_multipoles_simd" == xunrolled; then
    AC_DEFINE_UNQUOTED(UNROLLEDMULTIPOLES, [], "Unrolled plain-C multipoles/Taylors")
    MULTIPOLES_SIMD="MULTIPOLES_SIMD_UNROLLED"
elif test "x$with_multipoles_simd" != xnone; then
    AC_MSG_ERROR([Unknown multipoles simd option "$with_multipoles_simd"!])
fi

# end multipoles SIMD

AC_ARG_ENABLE(avx-fof,
[AS_HELP_STRING([--enable-avx-fof], [AVX FOF [yes if AVX present]])],
    [], enable_avx_fof=$enable_avx)
if test x$enable_avx_fof == xyes; then
    AC_DEFINE_UNQUOTED(AVXFOF, [], "Use AVX for FOF and SO")
fi

AC_ARG_WITH(max-cpd,
[AS_HELP_STRING([--with-max-cpd=CPD], [Maximum cells-per-dimension [8192]])],
	[], with_max_cpd="8192")
AC_DEFINE_UNQUOTED(MAXCPD,$with_max_cpd, "Maximum cells-per-dimension")

AC_ARG_WITH(max-source-length,
[AS_HELP_STRING([--with-max-source-length=S], [Maximum source length [1048576]])],
	[], with_max_source_length="1048576")
AC_DEFINE_UNQUOTED(MAXSOURCELENGTH,$with_max_source_length, "Maximum source length")

AC_ARG_WITH(max-order,
[AS_HELP_STRING([--with-max-order=O], [Maximum multipole order [8]])],
    [], with_max_order="8")
MAXORDER=$with_max_order
AC_SUBST(MAXORDER)
AC_DEFINE_UNQUOTED(MAXORDER,$with_max_order, "Maximum multipole order")

AC_ARG_ENABLE(threaded-io,
[AS_HELP_STRING([--enable-threaded-io], [Threaded IO [yes]])],
	[], enable_threaded_io="yes")
if test x$enable_threaded_io == xyes; then
	AC_DEFINE_UNQUOTED(IOTHREADED, [], "Use a separate thread for IO")
fi

AC_ARG_ENABLE(one-sided-group-finding,
[AS_HELP_STRING([--enable-one-sided-group-finding], [Perform L0 finding one-sided [yes]])],
	[], enable_one_sided_group_finding="yes")
if test x$enable_one_sided_group_finding == xyes; then
	AC_DEFINE_UNQUOTED(ONE_SIDED_GROUP_FINDING, [], "Use one-sided algorithm for L0 group finding")
fi

AC_ARG_ENABLE(gpu-directs,
[AS_HELP_STRING([--enable-gpu-directs], [GPU directs [yes if CUDA detected]])],
	[], enable_gpu_directs="check")
if test x$enable_gpu_directs == xcheck -o x$enable_gpu_directs == xyes; then
	AX_CHECK_CUDA
	if test x$enable_gpu_directs == xcheck; then
		enable_gpu_directs=$VALID_CUDA
	fi
fi
AC_SUBST(CUDA_ENABLED)
CUDA_ENABLED="CUDA_DISABLED"
if test x$enable_gpu_directs == xyes; then
	AC_DEFINE_UNQUOTED(CUDADIRECT, [], "Do the near force on the GPU")
	CUDA_ENABLED="CUDA_ENABLED"
fi

AC_ARG_ENABLE(gpu-convolution,
[AS_HELP_STRING([--enable-gpu-convolution], [Convolution FFT on GPU [no]])],
	[], enable_gpu_convolution="no")
AC_SUBST(GPU_CONVOLUTION_ENABLED)
GPU_CONVOLUTION_ENABLED="GPU_CONVOLUTION_DISABLED"
if test x$enable_gpu_convolution == xyes; then
	AC_DEFINE_UNQUOTED(GPUFFT, [], "Do the convolution FFT on the GPU")
    GPU_CONVOLUTION_ENABLED="GPU_CONVOLUTION_ENABLED"
fi

AC_ARG_ENABLE(fof_density,
[AS_HELP_STRING([--enable-fof-density], [Compute FOF-scale density [yes]])],
	[], enable_fof_density="yes")
if test x$enable_fof_density == xyes; then
	AC_DEFINE_UNQUOTED(COMPUTE_FOF_DENSITY, [], "Compute FOF-scale density on the GPU, used to accelerate the L0 FOF computation and seed the L1/L2 halo finding.  This promotes FLOAT3 accelerations to FLOAT4.")
fi

AC_ARG_ENABLE(spherical-overdensity,
[AS_HELP_STRING([--enable-spherical-overdensity], [Use SO finder [yes if fof-density enabled]])],
    [], enable_spherical_overdensity=$enable_fof_density)
if test x$enable_spherical_overdensity == xyes; then
    AC_DEFINE_UNQUOTED(SPHERICAL_OVERDENSITY, [], "Use SO finder for L1/L2 halos, rather than FOF")
fi

AC_ARG_ENABLE(parallel,
[AS_HELP_STRING([--enable-parallel], [Generate MPI parallel code [no]])],
	[], enable_parallel="no")
if test x$enable_parallel == xyes; then
	AC_DEFINE_UNQUOTED(PARALLEL, [], "Generate MPI-parallel code.")
fi

# Query the GPUs to find the number and compute capability
if test x$enable_gpu_directs == xyes; then
    AC_MSG_CHECKING([compute capability of GPUs])

    make --quiet -C singlestep/Direct detect_gpus MAKE_CONTEXT=AUTOCONF NVCC=$NVCC CXX=$CXX
    DETECT_GPUS_OUTPUT=$(./singlestep/Direct/detect_gpus)

    NGPU=$(echo $DETECT_GPUS_OUTPUT | \grep -oP '(?<=Device count: )\d+')
    CUDA_ARCH=$(echo $DETECT_GPUS_OUTPUT | \grep -oP '(?<=Compute capability: )\d+')

    if test x$NGPU == x""; then
        NGPU=0
    fi
    if test x$CUDA_ARCH == x""; then
        CUDA_ARCH=0
    fi
    AC_MSG_RESULT([compute capability $CUDA_ARCH])
    AC_SUBST(CUDA_ARCH)
fi

AC_ARG_WITH(max-gpus,
[AS_HELP_STRING([--with-max-gpus=M], [Maximum number of GPUs to use [auto]])],
	[], with_max_gpus=32)
if test x$with_max_gpus == xauto; then
    with_max_gpus=$NGPU
fi
AC_DEFINE_UNQUOTED(MAX_GPUS,$with_max_gpus, "Limit on the number of GPUs")

AC_ARG_WITH(softening,
[AS_HELP_STRING([--with-softening={plummer | single_spline | cubic_spline | cubic_plummer}], [Near-force softening technique [single_spline]])],
	[], with_softening="single_spline")
if test "x$with_softening" == xcubic_plummer; then
	AC_DEFINE_UNQUOTED(DIRECTCUBICPLUMMER, [], "Cubic variant of Plummer softening")
elif test "x$with_softening" == xcubic_spline; then
	AC_DEFINE_UNQUOTED(DIRECTCUBICSPLINE, [], "Cubic Hernquist spline softening")
elif test "x$with_softening" == xsingle_spline; then
	AC_DEFINE_UNQUOTED(DIRECTSINGLESPLINE, [], "Single-branch spline softening")
elif test "x$with_softening" == xplummer; then
	AC_DEFINE_UNQUOTED(DIRECTPLUMMER, [], "Plummer softening")
else
    AC_MSG_ERROR([Unknown softening "$with_softening"!])
fi

# If not using MKL, look for FFTW
if test x$enable_mkl != xyes; then
    AC_SEARCH_LIBS(fftw_execute, [fftw3], [], [AC_MSG_ERROR([unable to find fftw3!])])
    AC_SEARCH_LIBS(fftw_init_threads, [fftw3_omp], [], [AC_MSG_ERROR([unable to find fftw3_omp!])])

    # Single-precision libs are only used for analysis routines right now
    AC_SEARCH_LIBS(fftwf_execute, [fftw3f], [], [AC_MSG_WARN([unable to find fftw3f!])])
    AC_SEARCH_LIBS(fftwf_init_threads, [fftw3f_omp], [], [AC_MSG_WARN([unable to find fftw3f_omp!])])
fi

AC_ARG_WITH(cache-line-size,
[AS_HELP_STRING([--with-cache-line=N], [Cache line size in bytes [auto]])],
    [], with_cache_line="auto")
if test "x$with_cache_line" == xauto; then
    # This is the only way to do this I've found that works on x86 and ppc
    with_cache_line=$(cat /sys/devices/system/cpu/cpu0/cache/index0/coherency_line_size)
fi
AC_DEFINE_UNQUOTED(CACHE_LINE_SIZE, $with_cache_line, "Cache line size in bytes")

AC_ARG_WITH(page-size,
[AS_HELP_STRING([--with-page-size=N], [Page size in bytes [auto]])],
    [], with_page_size="auto")
if test "x$with_page_size" == xauto; then
    with_page_size=$(getconf PAGESIZE)
fi
AC_DEFINE_UNQUOTED(PAGE_SIZE, $with_page_size, "Page size in bytes")

AC_ARG_ENABLE(numa-for,
[AS_HELP_STRING([--enable-numa-for], [Enable custom NUMA_FOR OpenMP scheduler [yes]])],
    [], enable_numa_for="yes")
if test x$enable_numa_for == xyes; then
    AC_DEFINE_UNQUOTED(ENABLE_NUMA_FOR, [], "Enable custom NUMA_FOR OpenMP scheduler")
fi

# If the user is using Anaconda Python, we have to be careful
# about how we build C/C++ libraries that will be loaded by Python.
# Anaconda Python has the Python executable rpath set to find
# Anaconda libs before system ones, so this may resolve shared
# libraries to different versions than the compiler expected.
# Often this is not too harmful, but we've observed OpenMP/-lgomp
# crashes as a result.  The solution is to always compile
# Python extensions with Anaconda's compilers if the user
# is using Anaconda Python.
# You may only want the impl packages; the activation scripts will mess with configure
PYTHON_EXEC=$(command -v python)
PYTHON_CC=$REAL_CC
PYTHON_CXX=$REAL_CXX
PYTHON_LDFLAGS=
PYTHON_CPPFLAGS=
PYTHON_COMPILER_WARNING=
AC_SUBST(PYTHON_CC)
AC_SUBST(PYTHON_CXX)
AC_SUBST(PYTHON_LDFLAGS)
AC_SUBST(PYTHON_CPPFLAGS)
# Check if the user is using Anaconda Python (better way to do this? sys.version doesn't work for me)
AC_MSG_CHECKING([for Anaconda Python])
if test "${PYTHON_EXEC#*conda}" != "$PYTHON_EXEC"; then
    AC_MSG_RESULT([yes])
    # Now look in the Anaconda/bin directory
    CONDA_BINDIR=$(dirname $PYTHON_EXEC)
    CONDA_BINDIRNAME=$(basename $CONDA_BINDIR)
    CONDA_LIBDIR=$(dirname $CONDA_BINDIR)/lib
    CONDA_INCLUDEDIR=$(dirname $CONDA_BINDIR)/include

    # sanity check
    if test xbin != x$CONDA_BINDIRNAME; then
        AC_MSG_WARN('bin' not found in directory containing 'python': "$CONDA_BINDIR" ("$CONDA_BINDIRNAME"). May fail to find conda compilers.)
    fi

    # There should be only one gcc and one g++ compiler in the bin directory.
    # All the ways to construct the exact compiler name seemed extremely complicated,
    # but we should keep an eye out for a better way to do this.
    # Directly invoking conda-build is a no-go: takes forever and won't compile an extension in-place
    GCC_FILES=$(ls -1q $CONDA_BINDIR/*gcc 2> /dev/null)
    # we count words here instead of lines; fixes a problem with empty strings and newlines but will break if there are spaces in the filenames
    NGCC_FILES=$(echo -n "$GCC_FILES" | wc -w)
    if test $NGCC_FILES -ge 1; then
        :
        # TODO: breaks systems with different binutils version (https://github.com/ContinuumIO/anaconda-issues/issues/11152)
        #PYTHON_CC=$(echo $GCC_FILES | awk '{print $1}')
    fi
    if test $NGCC_FILES -ne 1; then
        AC_MSG_WARN(Could not find a unique *gcc file in $CONDA_BINDIR. Python extensions may not use the correct compiler for this platform!  Install the Anaconda gcc_impl_linux-64 package (or the equivalent for your platform) to get rid of this message.)
        AC_MSG_WARN(Using: "$PYTHON_CC")
    fi

    GXX_FILES=$(ls -1q $CONDA_BINDIR/*g++ 2> /dev/null)
    NGXX_FILES=$(echo -n "$GXX_FILES" | wc -w)
    if test $NGXX_FILES -ge 1; then
        :
        # TODO: breaks systems with different binutils version (https://github.com/ContinuumIO/anaconda-issues/issues/11152)
        #PYTHON_CXX=$(echo $GXX_FILES | awk '{print $1}')
    fi
    if test $NGXX_FILES -ne 1; then
        AC_MSG_WARN(Could not find a unique *g++ file in $CONDA_BINDIR. Python extensions may not use the correct compiler for this platform!  Install the Anaconda gxx_impl_linux-64 package (or the equivalent for your platform) to get rid of this message.)
        AC_MSG_WARN(Using: "$PYTHON_CXX")
    fi

    # Note Anaconda PYTHON_CC is always a gcc version, so any FLAGS will have to respect that
    # Let's assume that Anaconda Python always does the rpath linkage (which I think is true),
    # so we might as well supply the lib at link time, no matter if we actually found a conda gcc
    PYTHON_LDFLAGS="-L $CONDA_LIBDIR"
    PYTHON_CPPFLAGS="-I $CONDA_INCLUDEDIR"

    # Most likely the user is using Python compiled with GCC.  So if we are building extensions
    # using icc (which shouldn't happen if the user installed the conda compilers), link some
    # Intel libs
    if test x$PYTHON_CXX == xicc; then
        PYTHON_LDFLAGS+=' -lirc'
    fi
else
    AC_MSG_RESULT([no])
fi

AC_ARG_WITH(ramdisk,
    [AS_HELP_STRING([--with-ramdisk=PATH], [Path to ramdisk [/dev/shm/]])],
    [], with_ramdisk="/dev/shm/")
AC_CHECK_FILE($with_ramdisk, [], [AC_MSG_WARN(Ramdisk path "$with_ramdisk" does not exist.)])
AC_DEFINE_UNQUOTED(RAMDISK_PATH, "$with_ramdisk", "Path to the ramdisk")

# This will be placed at the top of *.in files
AC_SUBST(WARNING)
WARNING="# WARNING! This file has been automatically generated by configure.
# Any changes to it will be overwritten the next time configure is run.
# You may want to modify the corresponding *.in file, or use configure
# flags to achieve your goal instead (see ./configure --help)."

echo
echo "------------------------------------------------------"
echo "Abacus has been configured with the following options:"
echo "Double precision:                $enable_double"
echo "Near field max radius:           $with_radius"
echo "Near field block size:           $with_block_size"
echo "AVX directs:                     $enable_avx_directs"
echo "AVX-512 directs:                 $enable_avx512_directs"
echo "SIMD multipoles:                 $with_multipoles_simd"
echo "AVX FOF:                         $enable_avx_fof"
echo "Spherical Overdensity:           $enable_spherical_overdensity"
echo "One-sided L0 Finding:            $enable_one_sided_group_finding"
echo "Threaded IO:                     $enable_threaded_io"
echo "GPU directs:                     $enable_gpu_directs"
echo "Compute FOF-scale density:       $enable_fof_density"
echo "MPI parallel code:               $enable_parallel"
echo "Near-force softening technique:  $with_softening"
echo "GPU Convolution FFT:             $enable_gpu_convolution"
echo "Cache line size:                 $with_cache_line"
echo "Page size:                       $with_page_size"
echo "CXX:                             $CXX"
echo "MKL FFT:                         $enable_mkl"
echo "Maximum order:                   $with_max_order"
echo "NUMA For OpenMP Scheduler:       $enable_numa_for"
echo "------------------------------------------------------"
echo

# These are configurable but uncommon to change
#echo "Maximum CPD:                     $with_max_cpd"
#echo "Maximum source length:           $with_max_source_length"
#echo "Maximum GPUs:                    $with_max_gpus"
#echo "Ramdisk path:                    $with_ramdisk"

AC_OUTPUT
