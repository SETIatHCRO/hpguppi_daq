#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([FULL-PACKAGE-NAME],[VERSION],[BUG-REPORT-ADDRESS])
AM_INIT_AUTOMAKE([foreign subdir-objects])
LT_INIT
AM_SILENT_RULES([yes])
AC_CONFIG_SRCDIR([hpguppi_databuf.h])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_MACRO_DIR([m4])

# Set CFLAGS to nothing if it is not set by the user.  This prevents AC_PROG_CC
# from setting the (supposedly reserved-for-the-user!) variable CFLAGS in
# Makefile, which prevents AM_CFLAGS in Makefile.am from setting an
# optimization level.  For more details, see
# http://lists.gnu.org/archive/html/autoconf/2006-04/msg00007.html
AS_VAR_SET_IF(CFLAGS,[],[CFLAGS=])

# Checks for programs.
AC_PROG_CC
AC_PROG_CXX

# Check processor flags, set empty X86_FEATURE_CFLAGS if not x86 arch
AC_RUN_IFELSE(
    [AC_LANG_PROGRAM([#include <stdlib.h>],
        [__builtin_cpu_init ();]
    )],
    [AX_CHECK_X86_FEATURES],
    [X86_FEATURE_CFLAGS=]
)

AM_CONDITIONAL([HAVE_AVX512F],
               [test "x$ax_cv_gcc_x86_cpu_supports_avx512f" = xyes])

# Check for hashpipe
AX_CHECK_HASHPIPE

# Check for ibverbs
AX_CHECK_IBVERBS

# Check for sla
AX_CHECK_SLA

# Check for rawspec
# AX_CHECK_RAWSPEC

# Check for xgpu
AX_CHECK_XGPU

# Check for cuda
AX_CHECK_CUDA

# Check for blade
AX_CHECK_BLADE

# Check for uvh5
AX_CHECK_UVH5C99

# Check for radiointerferometry
AX_CHECK_RADIOINTERFEROMETRY99

# Check for filterbankh5
AX_CHECK_FILTERBANKC99

# Checks for libraries.
#AC_CHECK_LIB([pthread], [pthread_create])
AC_CHECK_LIB([m], [cos])
#AC_CHECK_LIB([rt], [clock_gettime])
#AC_CHECK_LIB([z], [crc32])

## Checks for header files.
#AC_CHECK_HEADERS([netdb.h stdint.h stdlib.h string.h sys/socket.h sys/time.h unistd.h zlib.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_INLINE
AC_TYPE_INT32_T
AC_TYPE_INT64_T
AC_TYPE_OFF_T
AC_TYPE_SIZE_T
AC_TYPE_UINT32_T
AC_TYPE_UINT64_T
AC_TYPE_UINT8_T

## Checks for library functions.
#AC_FUNC_MALLOC
#AC_CHECK_FUNCS([clock_gettime memset socket crc32])

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
