# This is the autoconf template file for ./configure
# Run ./autogen.sh to re-create ./configure from this file

# This approach to Abacus's configure system was inspired by
# Athena (https://trac.princeton.edu/Athena/)

# Initialization
AC_INIT([abacus], [0.9], [], [], [github.com/nbody/abacus])
AC_PROG_CC
AC_PROG_CXX
AC_HEADER_STDC
AC_PROG_INSTALL

# These files will be generated from their Makefile.in counterparts
AC_CONFIG_FILES([
 singlestep/Makefile
 singlestep/Direct/Makefile
 Analysis/PowerSpectrum/Makefile
 Convolution/Makefile
])

# Imports
m4_include(m4/ax_check_cuda.m4)
m4_include(m4/ax_check_compile_flag.m4)
m4_include(m4/ax_gcc_x86_avx_xgetbv.m4)
m4_include(m4/ax_gcc_x86_cpuid.m4)
m4_include(m4/ax_ext.m4)

# These will automatically add the appropriate -L and -l flags to LIBS
AC_SEARCH_LIBS(fftw_execute, [fftw3], [], [AC_MSG_ERROR([unable to find fftw3!])])
AC_SEARCH_LIBS(fftw_init_threads, [fftw3_omp], [], [AC_MSG_ERROR([unable to find fftw3_omp!])])
AC_SEARCH_LIBS([cos],[m])
AC_SEARCH_LIBS(cblas_dgemm, [gslcblas], [], [AC_MSG_ERROR([unable to find gslcblas!])])
AC_SEARCH_LIBS(gsl_rng_uniform, [gsl], [], [AC_MSG_ERROR([unable to find gsl!])])
AC_SEARCH_LIBS(ITT_DoOneTimeInitialization, [tbb], [], [AC_MSG_ERROR([unable to find tbb!])])

# Check for some utilities we use
AC_CHECK_PROGS(LEX, flex)
AC_CHECK_PROGS(YACC, bison)
AC_CHECK_PROGS(AR, ar)

# Check for AVX
AX_EXT([no-define])

# Define the ./configure options
AC_ARG_ENABLE(double,
[AS_HELP_STRING([--enable-double], [Double precision [no]])],
	[], enable_double="no")
if test x$enable_double == xyes; then
	AC_DEFINE_UNQUOTED(DOUBLEPRECISION)
fi

AC_ARG_WITH(radius,
	[AS_HELP_STRING([--with-radius=NFR], [Near field radius [2]])],
	[], with_radius="2")
AC_DEFINE_UNQUOTED(NFRADIUS,$with_radius)

AC_ARG_WITH(block-size,
[AS_HELP_STRING([--with-block-size=B], [Near field block size [64]])],
	[], with_block_size="64")
AC_DEFINE_UNQUOTED(NFBlockSize,$with_block_size)

AC_ARG_ENABLE(avx,
[AS_HELP_STRING([--enable-avx], [All AVX (overridden by individual options) [yes if AVX present]])],
	[], enable_avx=$ax_cv_have_avx_cpu_ext)

AC_ARG_ENABLE(avx-directs,
[AS_HELP_STRING([--enable-avx-directs], [AVX directs [yes if AVX present]])],
	[], enable_avx_directs=$enable_avx)
if test x$enable_avx_directs == xyes; then
	AC_DEFINE_UNQUOTED(AVXDIRECT)
fi

AC_ARG_ENABLE(avx-multipoles,
[AS_HELP_STRING([--enable-avx-multipoles], [AVX multipoles [yes if AVX present]])],
	[], enable_avx_multipoles=$enable_avx)
if test x$enable_avx_multipoles == xyes; then
	AC_DEFINE_UNQUOTED(AVXMULTIPOLES)
fi

AC_ARG_WITH(max-cpd,
[AS_HELP_STRING([--with-max-cpd=CPD], [Maximum cells-per-dimension [8192]])],
	[], with_max_cpd="8192")
AC_DEFINE_UNQUOTED(MAXCPD,$with_max_cpd)

AC_ARG_WITH(max-source-length,
[AS_HELP_STRING([--with-max-source-length=S], [Maximum source length [1048576]])],
	[], with_max_source_length="1048576")
AC_DEFINE_UNQUOTED(MAXSOURCELENGTH,$with_max_source_length)

AC_ARG_ENABLE(threaded-io,
[AS_HELP_STRING([--enable-thread-io], [Threaded IO [yes]])],
	[], enable_threaded_io="yes")
if test x$enable_threaded_io == xyes; then
	AC_DEFINE_UNQUOTED(IOTHREADED)
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
AC_SUBST(CUDA_MODE)
CUDA_MODE="CUDA_DISABLED"
if test x$enable_gpu_directs == xyes; then
	AC_DEFINE_UNQUOTED(CUDADIRECT)
	CUDA_MODE="CUDA_ENABLED"
fi

AC_ARG_ENABLE(gpu-convolution,
[AS_HELP_STRING([--enable-gpu-convolution], [Convolution FFT on GPU [no]])],
	[], enable_gpu_convolution="no")
if test x$enable_gpu_convolution == xyes; then
	AC_DEFINE_UNQUOTED(GPUFFT)
fi

# Get the number of GPUs
if test x$enable_gpu_directs == xyes; then
	AC_MSG_CHECKING([number of GPUs])
	NGPU=m4_esyscmd([nvidia-smi -q | grep 'Attached GPUs'| awk -F ':' '{print $NF}' | tr -d [[:space:]]])
    if test x$NGPU == x""; then
        NGPU=0
    fi
	AC_MSG_RESULT([$NGPU found])
fi

# If we don't have any GPUs, pretend there's one for BPD
if test x$NGPU == x -o x$NGPU == x0; then
    NGPU=1
fi
DEFAULT_BPD=$((2*$with_radius/$NGPU + 1))
AC_ARG_WITH(direct-bpd,
[AS_HELP_STRING([--with-direct-bpd=B],
[Direct buffers-per-device [ceil((2*NFR+1)/NGPU) = 3]])],
	[], with_direct_bpd=$DEFAULT_BPD)
AC_DEFINE_UNQUOTED(DirectBPD,$with_direct_bpd)


AC_ARG_WITH(softening,
[AS_HELP_STRING([--with-softening={plummer | spline | cubic}], [Near-force softening technique [plummer]])],
	[], with_softening="plummer")
if test "x$with_softening" == xcubic; then
	AC_DEFINE_UNQUOTED(DIRECTCUBIC)
elif test "x$with_softening" == xspline; then
	AC_DEFINE_UNQUOTED(DIRECTSPLINE)
fi

echo
echo "------------------------------------------------------"
echo "Abacus has been configured with the following options:"
echo "Double precision:                $enable_double"
echo "Near field radius:               $with_radius"
echo "Near field block size:           $with_block_size"
echo "AVX directs:                     $enable_avx_directs"
echo "AVX multipoles:                  $enable_avx_multipoles"
echo "Maximum CPD:                     $with_max_cpd"
echo "Maximum source length:           $with_max_source_length"
echo "Threaded IO:                     $enable_threaded_io"
echo "Direct buffers-per-device:       $with_direct_bpd"
echo "GPU directs:                     $enable_gpu_directs"
echo "Near-force softening technique:  $with_softening"
echo "GPU Convolution FFT:             $enable_gpu_convolution"
echo "------------------------------------------------------"
echo

AC_OUTPUT