##### 
#
# SYNOPSIS
#
# AX_CHECK_CUDA
#
# DESCRIPTION
#
# Figures out if CUDA Driver API/nvcc is available, i.e. existence of:
# 	cuda.h
#   libcuda.so
#   nvcc
#
# This version is hacked to use -lcudart intead of -lcuda
#
# If something isn't found, fails straight away.
#
# Locations of these are included in 
#   CUDA_CFLAGS and 
#   CUDA_LDFLAGS.
# Path to nvcc is included as
#   NVCC_PATH
# in config.h
# VALID_CUDA is set depending on whether CUDA was found.
# NVCC is exported to .in files.
#
# 
# The author is personally using CUDA such that the .cu code is generated
# at runtime, so don't expect any automake magic to exist for compile time
# compilation of .cu files.
#
# LICENCE
# Public domain
#
# AUTHOR
# wili
#
##### 

AC_DEFUN([AX_CHECK_CUDA], [

# There are a few semi-canonical variables that indicate the installed CUDA prefix
# Just take the first one that exists
# We could loop through these if we wanted
default_cuda_prefix=$CUDA_HOME
default_cuda_prefix=${default_cuda_prefix:-$CUDA_DIR}
default_cuda_prefix=${default_cuda_prefix:-$CUDAPATH}
default_cuda_prefix=${default_cuda_prefix:-$CUDA_BASE}
default_cuda_prefix=${default_cuda_prefix:-/usr/local/cuda}

# Provide your CUDA path with this

AC_ARG_WITH(cuda, [AS_HELP_STRING([--with-cuda=PREFIX], [Prefix of your CUDA installation [$CUDA_HOME or $CUDA_DIR or $CUDAPATH or /usr/local/cuda]])], [cuda_prefix=$withval], [cuda_prefix=$default_cuda_prefix])

# Checking for nvcc
AC_MSG_CHECKING([nvcc in $cuda_prefix/bin])
if test -x "$cuda_prefix/bin/nvcc"; then
	AC_MSG_RESULT([found])
	AC_DEFINE_UNQUOTED([NVCC_PATH], ["$cuda_prefix/bin/nvcc"], [Path to nvcc binary])
	# We need to add the CUDA search directories for header and lib searches

	AC_SUBST([CUDA_MAJOR_VER],[$(cat $cuda_prefix/version.txt |cut -d' ' -f3 | cut -d. -f1)])

	CUDA_CFLAGS=""

	# Saving the current flags
	ax_save_CFLAGS="${CFLAGS}"
	ax_save_LDFLAGS="${LDFLAGS}"

	# Announcing the new variables
	AC_SUBST([CUDA_CFLAGS])
	AC_SUBST([CUDA_LDFLAGS])
	AC_SUBST([NVCC],[$cuda_prefix/bin/nvcc])
	AC_CHECK_FILE([$cuda_prefix/lib64],[lib64_found=yes],[lib64_found=no])
	if test "x$lib64_found" = xno ; then
		AC_CHECK_FILE([$cuda_prefix/lib],[lib32_found=yes],[lib32_found=no])
		if test "x$lib32_found" = xyes ; then
			AC_SUBST([CUDA_LIBDIR],[$cuda_prefix/lib])
		else
			AC_MSG_WARN([Could not find cuda lib directory])
			VALID_CUDA=no
		fi
	else
		AC_CHECK_SIZEOF([long])
		if test "x$ac_cv_sizeof_long" = "x8" ; then
			AC_SUBST([CUDA_LIBDIR],[$cuda_prefix/lib64])
			CUDA_CFLAGS+=" -m64"
		elif test "x$ac_cv_sizeof_long" = "x4" ; then
			AC_CHECK_FILE([$cuda_prefix/lib32],[lib32_found=yes],[lib32_found=no])
			if test "x$lib32_found" = xyes ; then
				AC_SUBST([CUDA_LIBDIR],[$cuda_prefix/lib])
				CUDA_CFLAGS+=" -m32"
			else
				AC_MSG_WARN([Could not find cuda lib directory])
				VALID_CUDA=no
			fi
		else
			AC_MSG_ERROR([Could not determine size of long variable type])
		fi
	fi

	if test "x$VALID_CUDA" != xno ; then
		CUDA_CFLAGS+=" -I$cuda_prefix/include"
		CFLAGS="$CUDA_CFLAGS $CFLAGS"
		CUDA_LDFLAGS="-L$CUDA_LIBDIR"
		LDFLAGS="$CUDA_LDFLAGS $LDFLAGS"

		# And the header and the lib
		#AC_CHECK_HEADER([cuda.h], [],
		#	AC_MSG_WARN([Could not find cuda.h])
		#	VALID_CUDA=no
		#	,[#include <cuda.h>])
		if test "x$VALID_CUDA" != "xno" ; then
			AC_CHECK_LIB([cuda], [cuInit], [VALID_CUDA=yes], [AC_MSG_WARN([Could not find libcuda])
            VALID_CUDA=no])
		fi
	fi
	# Returning to the original flags
	CFLAGS=${ax_save_CFLAGS}
	LDFLAGS=${ax_save_LDFLAGS}
else
	AC_MSG_RESULT([not found!])
	AC_MSG_WARN([nvcc was not found in $cuda_prefix/bin])
	VALID_CUDA=no
fi

if test "x$enable_cuda" = xyes && test x$VALID_CUDA = xyes ; then 
	AC_MSG_NOTICE([Building with CUDA bindings])
elif test "x$enable_cuda" = xyes && test x$VALID_CUDA = xno ; then 
	AC_MSG_ERROR([Cannot build CUDA bindings. Check errors])
fi


])